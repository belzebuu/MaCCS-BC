#!/usr/bin/python
import argparse
import logging
import time
from multiprocessing import Pool, cpu_count
import common.constants as const
import common.data as d
import model.model as model


def perm_test(solver):
    try:
        return solver.solve()
    except Exception:
        logging.exception("solve %s failed" % solver.instance_name)


def main(args):
    mode = const.MODE_SOLVE
    solver = model.Solver(mode, args.interactions, args.mutations, args.k, args.weights,
                          args.exclusive, args.verbose, args.time_limit, args.preprocessing)
    base_obj, sol_set, _, _ = solver.solve()

    if not base_obj:
        print("Solution not found")
        return

    mode = const.MODE_VERSUS
    threads = cpu_count()
    if args.threads:
        threads = args.threads

    pool = Pool(processes=threads)
    a_args = list(range(args.n))
    b_args = [model.Solver(mode, args.interactions, args.mutations, args.k, args.weights,
                           args.exclusive, args.verbose, args.time_limit, args.preprocessing,
                           data=d.Data.from_file(args.interactions, args.mutations, args.weights).permute_mutations(seed),
                           base_obj=base_obj,
                           sol_set=sol_set) for seed in a_args]
    start_time = time.time()
    result = pool.map(perm_test, b_args)
    pool.close()
    pool.join()

    print(result)
    count = 0
    for i in result:
        count += i
    del a_args
    del b_args
    print("# Parallel test completed in: %s" % (time.time()-start_time))

    print("count = %d" % count)
    print("p value %s" % ((1.0 + count) / (1.0 + args.n)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Choose an instance to solve')
    parser.add_argument('-i', '--interactions', type=str,
                        help='a path to the file containing interactions', required=True)
    parser.add_argument('-m', '--mutations', type=str, help='a path to the file containing mutations', required=True)
    parser.add_argument('-k', type=int, help='number of nodes to choose from the gene nodes', required=True)
    parser.add_argument('-w', '--weights', type=str, help='a path to the file containing node (gene) weights')
    parser.add_argument('-o', '--output_file', type=str,
                        help='file for saving solutions. An extension (.json) will be added')
    parser.add_argument('-n', type=int, help='number of instances in permutation test', required=True)
    parser.add_argument('-x', '--exclusive', dest='exclusive', action='store_true', help='use exclusive objective')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='show verbose output')
    parser.add_argument('-t', '--time_limit', type=int, help='solver time limit (in seconds)')
    parser.add_argument('-p', '--preprocessing', dest="preprocessing", action="store_false",
                        help='preprocessing of the instance')
    parser.add_argument('-a', '--all', dest="allsolutions", action="store_true",help='find all solutions')
    parser.add_argument('--threads', type=int,
                        help='number of processes to execute simultaneously (default number of cores)', required=False)
    parser.set_defaults(exclusive=False)
    parser.set_defaults(verbose=False)
    parser.set_defaults(preprocessing=True)
    parser.set_defaults(allsolutions=False)
    args = parser.parse_args()
    main(args)
