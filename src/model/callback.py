from gurobipy import GRB, quicksum, GurobiError
from . import helpers as h
from . import sep as s
import time


def mipsol_info(model):
    nodecnt = int(model.cbGet(GRB.callback.MIPSOL_NODCNT))
    solcnt = model.cbGet(GRB.callback.MIPSOL_SOLCNT)
    obj = model.cbGet(GRB.callback.MIPSOL_OBJ)
    objbnd = model.cbGet(GRB.callback.MIPSOL_OBJBND)
    runtime = model.cbGet(GRB.Callback.RUNTIME)
    return nodecnt, solcnt, obj, objbnd, runtime


def mipnode_info(model):
    nodecnt = int(model.cbGet(GRB.callback.MIPNODE_NODCNT))
    solcnt = model.cbGet(GRB.callback.MIPNODE_SOLCNT)
    obj = model.cbGet(GRB.callback.MIPNODE_OBJBST)
    objbnd = model.cbGet(GRB.callback.MIPNODE_OBJBND)
    status = model.cbGet(GRB.Callback.MIPNODE_STATUS)
    runtime = model.cbGet(GRB.Callback.RUNTIME)

    return nodecnt, solcnt, obj, objbnd, status, runtime


def print_mipnode(model):
    nodecnt, solcnt, obj, objbnd, status, runtime = mipnode_info(model)
    if model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.Status.OPTIMAL:
        print('MIPNODE NODE : %d    SOLCNT:  %d   OBJBST:  %d   OBJBND:  %s      STATUS:  %d' \
                      % (nodecnt, solcnt, obj, objbnd, status))


def node_sep(model, where):
    time_start = time.time()
    if where == GRB.callback.MIPSOL:
        nodecnt, solcnt, obj, objbnd, runtime = mipsol_info(model)
        if solcnt > 0 or model._resolve:
            # UPPER BOUND LOWER THAN BASIC SOLUTION. TERMINATE.
            if model._base_obj:
                if objbnd < model._base_obj:
                    print("TERMINATED. UPPER BOUND < BASE SOL, %s < %s" % (objbnd, model._base_obj))
                    model._termination_indicator = 0
                    model.terminate()
            # Get statistics at root
            if nodecnt == 0:
                model._root_callback_count += 1
                model._root_solcnt = solcnt
                model._root_obj = obj
                model._root_objbnd = objbnd
                model._root_runtime = runtime

            model._callback_count += 1
            # Get values of current solution
            variables = model.cbGetSolution(model._vars)
            x_sol_keys = h.tolerant_int_solution_keys(variables[0:model._g_count],
                                                      model.params.IntFeasTol, model._node_list)
            # Find connected components and compute node separators
            violations = s.violated_separators(x_sol_keys, model._node_neighbors, model._add_all)

            if violations:
                for viol in violations:
                    try:
                        model.cbLazy(
                            model._vars[model._node_indices[viol.node_i]]
                            + model._vars[model._node_indices[viol.node_j]] - 1
                            <= quicksum(model._vars[model._node_indices[u]] for u in viol.separator)
                        )
                    except GurobiError as e:
                        print("Error %d" % e.errno)
                        print(e.message)
                model._active_callback_count += 1
                model._n_lazy += len(violations)
                if nodecnt == 0:
                    model._active_root_callback_count += 1
                    model._n_root_lazy += len(violations)
                print("Integer solution not feasible: found %d lazy constraints violated to be added" % len(violations))
                if model._base_obj:
                    if obj > model._base_obj:
                        print("TERMINATED. LOWER BOUND > BASE SOL, BND: %s" % obj)
                        model._termination_indicator = 1
                        model.terminate()
    elif where == GRB.callback.MIPNODE:
        print_mipnode(model)
    model._t_callback += time.time() - time_start
