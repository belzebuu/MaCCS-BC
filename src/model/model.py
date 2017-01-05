from gurobipy import *
import common.constants as constants
import callback
import common.data as d
import common.preprocessing as p
import common.result as r
import time
import copy


class Solver:
    def __init__(self, mode, k, exclusive=False, verbosity=5, time_limit=None, preprocessing=True, 
                data=None, base_obj=None, sol_set=None):
        self.mode = mode
        self.k = k
        self.exclusive = exclusive
        self.verbose = verbosity
        self.preprocessing = preprocessing
        self.time_limit = time_limit
        self.data = data
        self.components = d.Data()
        self.base_obj = base_obj
        self.sol_set = sol_set
        self.result = r.Result(mode, exclusive, time_limit)

    def solve(self):

        if self.preprocessing:
            print "Preprocessing instance..."            
            #preprocess = p.Preprocess()
            components = p.Preprocess.preprocess(self.data, self.k)

            if len(components) == 1:
                print "single component"
                self.components = components[0]
            else:
                print "%d components" % (len(components))


        objective, solution_set, termination_indicator = 0, set(), 0
        main_model, main_variables = None, list()

        for component in components:
            self.data = component
            model, variables = self.define_model()
            obj, sol, ti = self.solve_model(model, variables)

            if obj > objective or (self.mode == constants.MODE_VERSUS and ti > termination_indicator):
                objective = obj
                solution_set = sol
                termination_indicator = ti
                main_model = model
                main_variables = variables

        #main_model, main_variables = self.define_model()
        #objective, solution_set, termination_indicator = self.solve_model(main_model, main_variables)


        if self.mode == constants.MODE_SOLVE:
            return objective, solution_set, termination_indicator, self.result

        elif self.mode == constants.MODE_VERSUS and self.base_obj:
            return objective, solution_set, termination_indicator, self.result

        elif self.mode == constants.MODE_RESOLVE:
            solution_list = list()
            while solution_set:
                previous_sol = solution_set
                main_model = self.clear_statistics(main_model)
                main_model = self.add_solution_cut(main_model, main_variables, solution_set, objective)
                main_model._resolve = True
                print "model resolve ", main_model._resolve
                objective, solution_set, termination_indicator = self.solve_model(main_model, main_variables)
                if solution_set:
                    print "#%d,%s" % (objective, solution_set)
                    if solution_set == previous_sol:
                        print "Found the same solution"
                        return objective, solution_set, termination_indicator, self.result
                    else:
                        solution_list.append(solution_set)

            return objective, solution_list, termination_indicator, self.result
        else:
            print "Wrong mode or no solution supplied."
            return None

    # HELPER FUNCTIONS
    def solve_model(self, m, variables):
        x, y, coverage, exclusive = variables
        solution_set = set()
        objective = 0.0
        termination_indicator = 0

        time_start = time.time()
        m.optimize(callback.node_sep)
        runtime = time.time() - time_start

        if m.status == GRB.OPTIMAL:
            # get statistics and return solution

            if (self.verbose==0):
                print "Optimal solution found"
                print "Explored %d nodes" % m.NodeCount
                print "Best objective "+str(m.ObjVal)+", best bound "+str(m.ObjBound)+", gap "+str(m.MIPGap)+"%"


            self.result.solve_measurement.append(r.SolveMeasurement(m.NodeCount, m._callback_count, m._n_lazy,
                                                                    m._root_callback_count, m._n_root_lazy,
                                                                    m._root_objbnd, m._root_runtime, m._t_callback,
                                                                    runtime, m.ObjVal))
            objective = m.ObjVal
            for i in self.data.nodes:
                if x[i].X >= 1-m.params.IntFeasTol:
                    solution_set.add(i)
            if self.mode == constants.MODE_VERSUS:
                if m.ObjVal >= self.base_obj:
                    termination_indicator = 1

        elif m.status == GRB.TIME_LIMIT:
            print 'exceeded time limit'
        elif GRB.INFEASIBLE:
            print 'model infeasible'
        elif m.status == GRB.INTERRUPTED:
            if self.mode == constants.MODE_VERSUS:
                termination_indicator = 1
            else:
                print "Unknown interruption"

        return objective, solution_set, termination_indicator

    def define_model(self):
        if self.weights:
            model, variables = self.node_sep_w()
        else:
            model, variables = self.node_sep()
        return model, variables

        # done

    def node_indices(self):
        indices = dict()
        # for i in range(len(self.data.nodes)):
        #     indices[self.data.nodes[i]] = i
        for idx, val in enumerate(self.data.nodes):
            indices[val] = idx
        return indices

    # MODELS
    def node_sep(self):
        m = Model("node_separation_st")

        # # # SOLVER PARAMETERS
        m.params.Threads = 1
        m.params.LazyConstraints = 1
        m.params.ScaleFlag = 0
        m.params.IntFeasTol = 1e-9
        m.params.Method = 3
        m.params.MIPFocus = 1
        if self.verbose>=1:
            m.params.OutputFlag = 1
            m.params.DisplayInterval = self.verbose
        else:
            m.params.OutputFlag = 0
            m.params.DisplayInterval = 1

        if self.mode == constants.MODE_SOLVE:
            if self.time_limit:
                m.params.TimeLimit = self.time_limit

        # # # VARIABLES
        x = {}
        for i in self.data.nodes:
            x[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="x" + str(i))
        y = {}
        for j in self.data.patients:
            y[j] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="y" + str(j))

        coverage = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="coverage")
        exclusive = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="alternative")

        m.update()

        # # # CONSTRAINTS
        print "Adding constraints"
        m.addConstr(coverage == quicksum(y[j] for j in self.data.patients), "coverage")
        m.addConstr(exclusive == 2 * coverage - quicksum(x[i] * len(self.data.node_mutations[i])
                                                         for i in self.data.node_mutations.keys()), "alternative")
        # k vertices
        m.addConstr(quicksum(x[i] for i in self.data.nodes) == self.k)

        print "Adding coverage constraints"
        # constraints either satisfied or y = 0
        for j in self.data.patients:
            m.addConstr(
                (1 - y[j] + quicksum(x[i] for i in self.data.mutated_genes[j])) >= 1
            )

        print "Constraints posted"
        if self.exclusive:
            m.setObjective(exclusive, GRB.MAXIMIZE)
        else:
            m.setObjective(coverage, GRB.MAXIMIZE)

        # variables accessible in callback
        m._callback_count = 0
        m._interactions = self.data.interactions
        m._g_count = len(self.data.nodes)
        m._p_count = len(self.data.patients)
        m._x_values = None
        m._vars = m.getVars()
        m._violations = list()
        m._k = self.k
        m._node_indices = self.node_indices()
        m._node_list = self.data.nodes
        # additional for this model
        m._node_mut = self.data.node_mutations
        m._node_neighbors = self.data.node_neighbors
        m._base_obj = self.base_obj
        m._termination_indicator = None
        m._t_callback = 0.0
        m._n_lazy = 0
        m._n_root_lazy = 0
        m._root_solcnt = 0
        m._root_obj = 0
        m._root_objbnd = 0.0
        m._root_callback_count = 0
        m._active_callback_count = 0
        m._active_root_callback_count = 0
        m._root_runtime = 0.0
        m._resolve = False
        m._add_all = True # whether to add node separators computed from both sides (True) or only one 
        return m, [x, y, coverage, exclusive]

    def node_sep_w(self):
        m = Model("node_separation_w")

        # # # SOLVER PARAMETERS
        m.params.Threads = 1
        m.params.LazyConstraints = 1
        m.params.ScaleFlag = 0
        m.params.IntFeasTol = 1e-9
        m.params.Method = 3
        m.params.MIPFocus = 1

        if self.verbose>=1:
            m.params.OutputFlag = 1
            m.params.DisplayInterval = self.verbose
        else:
            m.params.OutputFlag = 0
            m.params.DisplayInterval = 1

        if self.mode == constants.MODE_SOLVE:
            if self.time_limit:
                m.params.TimeLimit = self.time_limit

        # # # VARIABLES
        x = {}
        for i in self.data.nodes:
            x[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="x" + str(i))
        y = {}
        z = {}
        covering_indices = {}
        for j in self.data.patients:
            y[j] = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="p" + str(j))
            j_links = 0
            covering_indices[j] = list()
            z[j] = list()
            for i in self.data.mutated_genes[j]:
                covering_indices[j].append(i)
                j_links += 1
                z[j].append(m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="z" + str(j) + "," + str(i)))

        coverage = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="coverage")

        m.update()

        # # # CONSTRAINTS
        print "Adding constraints"
        for i in self.data.nodes:
            x[i].setAttr(GRB.Attr.BranchPriority, 10)

        # weighted coverage - objective constraint
        m.addConstr(coverage == quicksum(y[j] for j in self.data.patients), "coverage")

        # k vertices constraint
        m.addConstr(quicksum(x[i] for i in self.data.nodes) == self.k)

        # weighted coverage
        print "Adding coverage constraints"
        if self.data.weights:
            print "adding weights constraints"
            for j in self.data.patients:
                # constraint # 7
                m.addConstr(
                        y[j] <= quicksum(z[j][idx] for idx in range(len(covering_indices[j])))
                )
                for idx, i in enumerate(covering_indices[j]):
                    # constraint # 5
                    m.addConstr(
                        x[i] >= z[j][idx]
                    )
                    # constraint # 6
                    m.addConstr(
                        y[j] <= (1 - z[j][idx]) + self.data.weights[i]*z[j][idx]
                    )

        print "Constraints posted"
        m.setObjective(coverage, GRB.MAXIMIZE)

        # variables accessible in callback
        m._callback_count = 0
        m._interactions = self.data.interactions
        m._g_count = len(self.data.nodes)
        m._p_count = len(self.data.patients)
        m._x_values = None
        m._vars = m.getVars()
        m._violations = list()
        m._k = self.k
        m._node_indices = self.node_indices()
        m._node_list = self.data.nodes
        # additional for this model
        m._node_mut = self.data.node_mutations
        m._node_neighbors = self.data.node_neighbors
        m._base_obj = self.base_obj
        m._termination_indicator = None
        m._t_callback = 0.0
        m._n_lazy = 0
        m._n_root_lazy = 0
        m._root_solcnt = 0
        m._root_obj = 0
        m._root_objbnd = 0.0
        m._root_callback_count = 0
        m._active_callback_count = 0
        m._active_root_callback_count = 0
        m._root_runtime = 0.0
        m._resolve = False

        return m, [x, y, coverage, z]

    def add_solution_cut(self, m, variables, sol_set, objective=None):
        m.params.OutputFlag = 0  # TODO ?
        x, _, coverage, alternative = variables
        print "Adding solution cut and objective value bound"
        assert len(sol_set) == self.k
        m.addConstr(quicksum(x[i] for i in sol_set) <= self.k - 1)
        if objective:
            if self.exclusive:
                m.addConstr(exclusive == objective)
            else:
                m.addConstr(coverage >= objective)
        m.update()
        return m

    @staticmethod
    def clear_statistics(m):
        m._callback_count = 0
        m._x_values = None
        m._violations = list()
        m._termination_indicator = None
        m._t_callback = 0.0
        m._n_lazy = 0
        m._n_root_lazy = 0
        m._root_solcnt = 0
        m._root_obj = 0
        m._root_objbnd = 0.0
        m._root_callback_count = 0
        m._active_callback_count = 0
        m._active_root_callback_count = 0
        m._root_runtime = 0.0
        m._resolve = False
        return m
