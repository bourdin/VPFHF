include ../makefile.include	
CFLAGS='-I..'

all: test1 test2 test3 test4 test5 test6 test7
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

test6: test6.o ${VFOBJ} chkopts
	@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test7: test7.o ${VFOBJ} chkopts
	@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

EPSILON =.05
LENGTH  =.4

runtest2: test2
	@echo "Running " $@
	mpiexec -n 2 ./$< -nx 41 -ny 41 -nz 6 -lx 1 -ly 1 -lz .25 -length ${LENGTH} -orientation 3 -nu 0 -epsilon ${EPSILON} -eta 1e-8 -p $@ | tee $@.out
	#@diff $@.out results/$@-ref.out

runtest3: test3
	@echo "Running " $@
	mpiexec -n 2 ./$< -nx 41 -ny 41 -nz 6 -lx 2 -ly 2 -lz .25 -length .8 -orientation 1 -nu 0 -p $@ | tee $@.out 
	@diff $@.out results/$@-ref.out

runtest4: test4
	@echo "Running " $@
	mpiexec -n 2 ./$< -nx 81 -ny 81 -nz 6 -lx 2 -ly 2 -lz .25 -length ${LENGTH} -orientation 3 -nu 0 -epsilon ${EPSILON} -eta 1e-8 -p $@ | tee $@.out
	@diff $@.out results/$@-ref.out

runtest5: test5
	@echo "Running " $@
	for or in 0 1 2 3 4 5; do \
	mpiexec -n 2 ./$< -nx 21 -ny 21 -nz 21 -lx 1 -ly 1 -lz 1 -nu 0 -epsilon ${EPSILON} -eta 1e-8 -orientation $$or -p $@-$$or  | tee $@-$$or.out ;\
	@diff $@-$$or.out results/$@-$$or-ref.out ;\
	done

runtest6: test6
	@echo "Running " $@
	for or in 0 1 2; do \
	mpiexec -n 2 ./$< -nx 21 -ny 21 -nz 21 -lx 1 -ly 1 -lz 1 -nu 0.3 -orientation $$or -p $@-$$or ;\
	@diff $@-$$or.out results/$@-$$or-ref.out ;\
	done

runtest7: test7
	@echo "Running " $@
	for or in 0 1 2; do \
	mpiexec -n 2 ./$< -nx 21 -ny 21 -nz 21 -lx 1 -ly 1 -lz 1 -nu 0.3 -orientation $$or -p $@-$$or ;\
	@diff $@-$$or.out results/$@-$$or-ref.out ;\
	done

clean::
	-rm -f test1 test2 test3 test4 test5
	-rm -f runtest*
	-rm -f TEST*
