#!/usr/bin/python
import argparse
import common.result as r


def main(args):
    """
    Program for verification of solution
    :return: True if the solution induces a connected subgraph, False otherwise
    """
    mode = 0 # TODO add some default value
    solution_manager = r.Solution(mode, args.solution_file, args.interactions, args.mutations, args.weights,
                                args.exclusive)
    solutions = solution_manager.get_solutions()
    for s in solutions:
        value, subgraph = solution_manager.split_solution(solutions, s)
        con = solution_manager.check_connectivity(subgraph)
        print "Number of nodes: %d; connectivity check: %r" % (len(subgraph),con)
        cov = len(solution_manager.coverage(subgraph))
        print "Coverage declared: %d; coverage verified: %d" %(value,cov)
        
    #print "Coverage: ",len(solution_saver.solution_coverage(solution_set))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Choose an instance to solve')
    parser.add_argument('-i', '--interactions', type=str, help='a path to the file containing interactions', required=True)
    parser.add_argument('-m', '--mutations', type=str, help='a path to the file containing mutations', required=True)
    parser.add_argument('-w', '--weights', type=str, help='a path to the file containing node (gene) weights')
    parser.add_argument('-o', '--solution_file', type=str, help='file where the solution is saved. (json format)', required=True)
    parser.add_argument('-x', '--exclusive', dest='exclusive', action='store_true', help='use exclusive objective')
    parser.set_defaults(exclusive=False)
    parser.set_defaults(allsolutions=False)
    args = parser.parse_args()
    main(args)

