import sys
import os.path
import json
import data as d
import networkx as nx
import common.constants as const


class Result:
    def __init__(self, mode, i, m, w=None, x=None, time_limit=None):
        """
        :param mode:
        :param i: (str)             interaction file path
        :param m: (str)             mutation file path
        :param k: (int)             k
        :param w: (str)             weight file path
        :param o: (str)             output file path
        :param x: (bool)            use exclusive objective
        :param time_limit: (int)    time limit in seconds, passed as a parameter to Gurobi
        :return:
        """
        self.mode = mode
        self.interaction_file = i
        self.mutation_file = m
        self.weights = w if w else ' '
        self.exclusive = x if x else ' '
        self.time_limit = time_limit if time_limit else sys.maxint
        self.solve_measurement = list()


    def __str__(self):
        s = Result.header() + '#%d,%s,%s,%s,%r,%d\n' % (self.mode, self.interaction_file, self.mutation_file,
                                                              self.weights, self.exclusive,
                                                              self.time_limit)
        m = ""
        if self.solve_measurement:
            m = SolveMeasurement.header()
            for i in self.solve_measurement:
                m += str(i)
        return s+m

    @staticmethod
    def header():
        return '#mode,interaction_file,mutation_file,k,weights,output,exclusive,time_limit,solve_measurement\n'


class SolveMeasurement:
    def __init__(self, node_cnt, callback_cnt, lazy_cnt, root_callback_cnt, root_lazy, root_objbnd,
                 root_runtime, callback_time, runtime, coverage):
        """
        Single run execution information, returned by solver (in model.model.py)
        :param node_cnt: (int)          number of nodes explored in branch and bound tree
        :param callback_cnt: (int)      number of callbacks from Gurobi
        :param lazy_cnt: (int)          number of lazy constraints added
        :param root_callback_cnt: (int) number of callbacks at the root node in B&B tree
        :param root_lazy: (int)         number of lazy constraints added at the root node in B&B tree
        :param root_objbnd: (float)     upper bound on the objective function at the root node in B&B tree
        :param root_runtime: (float)    time spent for solving the root node
        :param callback_time: (float)   time spent in callback functions
        :param runtime: (float)         running time of the solver
        :param coverage: (float)        found objective
        :return:
        """

        self.nodes = node_cnt
        self.callbacks = callback_cnt
        self.lazy_constraints = lazy_cnt
        self.callbacks_at_root = root_callback_cnt
        self.lazy_at_root = root_lazy
        self.ub_at_root = float(root_objbnd)
        self.time_at_root = root_runtime
        self.callback_time = callback_time
        self.total_time = runtime
        self.coverage = coverage

    @staticmethod
    def header():
        return '#nodes,callbacks,lazy constr.,callbacks at root,lazy at root,ub at root,time at root,callback time,' \
               'total time,coverage\n'

    def __str__(self):
        return "%d,%d,%d,%d,%d,%.2f,%.2f,%.2f,%.2f,%d" % (self.nodes, self.callbacks, self.lazy_constraints,
                                                          self.callbacks_at_root, self.lazy_at_root,
                                                          self.ub_at_root, self.time_at_root,
                                                          self.callback_time, self.total_time, self.coverage)


class Solution:
    def __init__(self, mode, output_file, interactions, mutations, weights=None, exclusive=None):
        """
        :param mode:                    solving mode: single solve, all solutions or bounded solver
        :param interactions:            interaction file path
        :param mutations:               mutation file path
        :param k:                       k
        :param output_file:             output file path
        :param weights:                 weight file path
        :param exclusive:               use exclusive objective
        :return:
        """
        self.data = d.Data.from_file(interactions, mutations, weights)
        self.output_file = output_file

        self.key = '_'.join([os.path.splitext(os.path.basename(interactions))[0],
                             os.path.splitext(os.path.basename(mutations))[0]])
        if weights:
            self.key += '_w'
        if exclusive:
            self.key += '_x'
        self.key = unicode(self.key, "utf-8")
        self.mode = mode

    def get_solutions(self):
        if not os.path.isfile(self.output_file):
            return dict()
        sol_file = open(self.output_file, 'r')
        solutions = json.load(sol_file)
        sol_file.close()
        return solutions

    def split_solution(self, solutions, key):
        """
        :return: set on nodes in the solution that has been saved for an instance of the problem
        """
        return solutions[key]['value'], solutions[key]['genes']

        #if self.key in solutions.keys():
        #    if self.mode == const.MODE_RESOLVE:
        #        return solutions[self.key]
        #    return solutions[self.key]
        #return None

    def write_solution(self, coverage, solution_set):
        solutions = dict() 

        sol_file = open(self.output_file, 'w')
        solutions[self.key] = {'value': coverage, 'genes': list(solution_set)}
        json.dump(solutions, sol_file, indent=4)
        sol_file.close()
        print "Solution saved"
        return 0

    def add_solution(self, coverage, solution_set):
        """
        :param coverage: (int)              value of the objective function
        :param solution_set: ({int})        set of nodes in the solution
        :return: 0 if solution was saved successfully, 1 otherwise
        """
        solutions = self.get_all_solutions()

        if self.key in solutions.keys():
            print "Solution already exists. To add, remove the old solution first. \nOld: %s \nNew:%s" \
                  % (set(solutions[self.key]), solution_set)
            return 1
        else:
            sol_file = open(self.output_file, 'w')
            solutions[self.key] = (coverage, list(solution_set))
            json.dump(solutions, sol_file)
            sol_file.close()
            print "Solution saved"
        return 0
    
    def add_solutions(self, solution_dict):
        solutions = self.get_all_solutions()

        if self.key in solutions.keys():
            print "Solution already exists. To add, remove the old solution first. \nOld: %s \nNew:%s" \
                  % (set(solutions[self.key]), solution_dict)
            return 1
        else:
            sol_file = open(self.output_file, 'w')
            solutions[self.key] = solution_dict
            json.dump(solutions, sol_file)
            sol_file.close()
            print "Solution saved"
        return 0

    def coverage(self, solution_set):
        coverage = set()

        if solution_set:
            for n in solution_set:
                if n in self.data.node_mutations:
                    coverage |= self.data.node_mutations[n]
        return coverage



    def weighted_coverage(self, solution_set):
        if not self.data.weights:
            return None

        w_coverage = dict()
        if solution_set:
            for n in solution_set:
                if n in data.node_mutations:
                    for p in data.node_mutations[n]:
                        if p in w_coverage:
                            if w_coverage[p] < data.weights[n]:
                                w_coverage[p] = data.weights[n]
                        else:
                            w_coverage[p] = data.weights[n]
        cov_sum = 0
        for p in w_coverage:
            cov_sum += w_coverage[p]
        return cov_sum

    def check_connectivity(self, solution_set):
        """
        :param solution_set: ({int})                set of nodes, that are believed to be a valid solution
        :return: True if the grap induced by solution set is connected, False otherwise
        """
        if solution_set:            
            G = self.data.make_graph()
            G_sub = G.subgraph(solution_set)
            con = nx.is_connected(G_sub)            
            return con
        else:
            return False

