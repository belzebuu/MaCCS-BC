#!/usr/bin/python
from gurobipy import *
import model.model as model
import argparse
import common.constants as const
from time import time
import common.result as cr


def main(args):
    if args.allsolutions:
        mode = const.MODE_RESOLVE
    else:
        mode = const.MODE_SOLVE

    try:
        start_time = time()
        solver = model.Solver(mode, args.interactions, args.mutations, args.genes, args.k, args.weights,
                              args.exclusive, args.verbose, args.time_limit, args.preprocessing)
        objective, solution_set, _, result = solver.solve()
        print result
        print 'Execution finished in %.3f seconds' % (time() - start_time)

        if solution_set and args.output_file is not None:
            solution_saver = cr.Solution(mode, args.output_file, args.interactions, args.mutations,
                                         args.weights, args.exclusive)
            if args.append:
                solution_saver.add_solution(objective, solution_set)
            else:
                solution_saver.write_solution(objective, solution_set)
                
    except GurobiError as e:
            print "Error %d" % e.errno
            print e.message


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Choose an instance to solve')
    parser.add_argument('-i', '--interactions', type=str,
                        help='a path to the file containing interactions', required=True)
    parser.add_argument('-m', '--mutations', type=str, help='a path to the file containing mutations', required=True)
    parser.add_argument('-g', '--genes', type=str, help='a path to the file containing genes', required=True)
    parser.add_argument('-k', type=int, help='number of nodes to choose from the gene nodes', required=True)
    parser.add_argument('-w', '--weights', type=str, help='a path to the file containing node (gene) weights')
    parser.add_argument('-o', '--output_file', type=str,
                        help='file for saving solutions. An extension (.json) will be added')
    parser.add_argument('-x', '--exclusive', dest='exclusive', action='store_true',
                        help='use exclusive objective')
    parser.add_argument('-v', '--verbose', dest='verbose', type=int,
                        help='set verbose level (0 = min verbose ... MAXINT = max verbose)')
    parser.add_argument('-t', '--time_limit', type=int, help='solver time limit (in seconds)')
    parser.add_argument('-p', '--preprocessing', dest="preprocessing", action="store_false",
                        help='preprocessing of the instance')
    parser.add_argument('-a', '--all', dest="allsolutions", action="store_true",help='find all solutions')
    parser.add_argument('-e', '--append', dest="append", action="store_true",help='append solution to json file')
    parser.set_defaults(exclusive=False)
    parser.set_defaults(verbose=False)
    parser.set_defaults(preprocessing=True)
    parser.set_defaults(allsolutions=False)
    parser.set_defaults(verbose=0)
    args = parser.parse_args()
    main(args)
