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

test5: test5.o ${VFOBJ} chkopts
	@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

runtest2: test2
	@echo "Running " $@
	mpiexec -n 2 ./$< -nx 41 -ny 41 -nz 6 -lz .25 -length .4 -orientation 1 -nu 0 -epsilon .05 -eta 1e-8 -p $@ | tee $@.out
	@diff $@.out results/$@-ref.out

runtest3: test3
	@echo "Running " $@
	mpiexec -n 2 ./$< -nx 41 -ny 41 -nz 6 -lz .25 -length .4 -orientation 1 -nu 0 -p $@ | tee $@-out 
	@diff $@.out results/$@-ref.out

runtest4: test4
	@echo "Running " $@
	@mpiexec -n 2 ./$< -nx 41 -ny 41 -nz 6 -lz .25 -lx 2 -ly 2 -length .4 -orientation 1 -nu 0 -epsilon .1 -eta 1e-8 -p $@ | tee $@.out
	@diff $@.out results/$@-ref.out

runtest5: test5
	@echo "Running " $@
	@mpiexec -n 2 ./$< -nx 41 -ny 11 -nz 11 -lx 1 -ly .25 -lz .25 -orientation 1 -nu 0 -epsilon .05 -eta 1e-8 -p $@ | tee $@.out
	@tail -20 $@.out | diff - results/$@-ref.out

clean::
	-rm -f test1 test2 test3 test4 test5
	-rm -f runtest*
	-rm -f TEST*
