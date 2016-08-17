init:
	pip install -r requirements.txt

test:
	src/solve.py -i laml -k 15

testall:
	src/solve_all.py -i laml -k 15

testtest:
	src/permutation_test.py -i laml -k 15


.PHONY: init test testall test
