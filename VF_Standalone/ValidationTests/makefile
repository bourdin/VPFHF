include ../makefile.include	
CFLAGS='-I..'

all: test1 test2
runall: runtest1

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules


test1: test1.o ${VFOBJ} chkopts
	@${CLINKER} -o test1 test1.o ${VFOBJ} ${PETSC_LIB}

test2: test2.o ${VFOBJ} chkopts
	@${CLINKER} -o test2 test2.o ${VFOBJ} ${PETSC_LIB}

test3: test3.o ${VFOBJ} chkopts
	@${CLINKER} -o test3 test3.o ${VFOBJ} ${PETSC_LIB}

test4: test4.o ${VFOBJ} chkopts
	@${CLINKER} -o test4 test4.o ${VFOBJ} ${PETSC_LIB}

runtest2: test2
	@echo "Running " $@
	mpiexec -n 2 ./test2 -nx 41 -ny 41 -nz 6 -lz .25 -length .4 -orientation 1 -nu 0 -epsilon .1 -eta 1e-8 -p TEST2 | tee runtest2.out
	@diff runtest2.out results/test2-ref.out

runtest3: test3
	@echo "Running " $@
	@mpiexec -n 2 ./test3 -nx 41 -ny 41 -nz 6 -lz .25 -length .4 -orientation 1 -nu 0 -p TEST3 | tee runtest3.out
	@diff runtest3.out results/test3-ref.out

runtest4: test4
	@echo "Running " $@
	@mpiexec -n 2 ./test4 -nx 41 -ny 41 -nz 6 -lz .25 -lx .5 -lx 2 -ly 2 -length .4 -orientation 1 -nu 0 -epsilon .1 -eta 1e-8 -p TEST4 | tee runtest4.out
	@diff runtest4.out results/test4-ref.out

clean::
	-rm -f test1