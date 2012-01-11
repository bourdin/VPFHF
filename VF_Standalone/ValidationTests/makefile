include ../makefile.include	
CFLAGS='-I..'

all: test1 test2 test3 test4 test5 test6 test7
test: runtest2 runtest3 runtest4 ruuntest5 runtest6 runtest7

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules


test1: test1.o ${VFOBJ} chkopts
	-${CLINKER} -o test1 test1.o ${VFOBJ} ${PETSC_LIB}

test2: test2.o ${VFOBJ} chkopts
	-${CLINKER} -o test2 test2.o ${VFOBJ} ${PETSC_LIB}

test3: test3.o ${VFOBJ} chkopts
	-${CLINKER} -o test3 test3.o ${VFOBJ} ${PETSC_LIB}

test4: test4.o ${VFOBJ} chkopts
	-@${CLINKER} -o test4 test4.o ${VFOBJ} ${PETSC_LIB}

test5: test5.o ${VFOBJ} chkopts
	-${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test6: test6.o ${VFOBJ} chkopts
	-${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test7: test7.o ${VFOBJ} chkopts
	-${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

EPSILON =.05
LENGTH  =.4

runtest2: test2
	@echo "Running " $@
	${MPIEXEC} -n 2 ./$< -nx 41 -ny 41 -nz 6 -lx 1 -ly 1 -lz .25 -length ${LENGTH} -orientation 3 -nu 0 -epsilon ${EPSILON} -eta 1e-8 -p $@ > $@.out 2>&1;\
	if (${DIFF} results/$@-ref.out $@.out) then true; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest3: test3
	@echo "Running" $@
	-@ARGS='-nx 41 -ny 41 -nz 6 -lx 2 -ly 2 -lz .25 -length .8 -orientation 1 -nu 0' ;\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@> $@.out 2>&1;\
	if (${DIFF} results/$@-ref.out $@.out) then echo "\t$< succeeded\n========================================="; \
	else echo ${PWD} ; echo "\tPossible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest4: test4
	@echo "Running " $@
	-@ARGS='-nx 81 -ny 81 -nz 6 -lx 2 -ly 2 -lz .25 -length ${LENGTH} -orientation 3 -nu 0 -epsilon ${EPSILON} -eta 1e-8';\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@> $@.out 2>&1;\
	if (${DIFF} results/$@-ref.out $@.out) then echo "\t$< succeeded\n========================================="; \
	else echo ${PWD} ; echo "\tPossible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest5: test5
	@echo "Running " $@
	-@for or in 0 1 2 3 4 5; do \
	ARGS="-nx 21 -ny 21 -nz 21 -lx 1 -ly 1 -lz 1 -nu 0 -epsilon ${EPSILON} -eta 1e-8 -orientation $$or";\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@-$$or > $@-$$or.out 2>&1;\
	if (${DIFF} results/$@-$$or-ref.out $@-$$or.out) then echo "\t$<-$$or succeeded\n========================================="; \
	else echo ${PWD} ; echo "\tPossible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@-$$or.out;\
	done

runtest6: test6
	@echo "Running " $@
	-@for or in 0 1 2; do \
	ARGS="-nx 21 -ny 21 -nz 21 -lx 1 -ly 1 -lz 1 -nu 0.3 -orientation $$or";\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@-$$or > $@-$$or.out 2>&1;\
	if (${DIFF} results/$@-$$or-ref.out $@-$$or.out) then echo "\t$<-$$or succeeded\n========================================="; \
	else echo ${PWD} ; echo "\tPossible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@-$$or.out;\
	done
	
runtest7: test7
	@echo "Running " $@
	-@for or in 0 1 2; do \
	ARGS="-nx 21 -ny 21 -nz 21 -lx 1 -ly 1 -lz 1 -nu 0.3 -orientation $$or";\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@-$$or > $@-$$or.out 2>&1;\
	if (${DIFF} results/$@-$$or-ref.out $@-$$or.out) then echo "\t$<-$$or succeeded\n========================================="; \
	else echo ${PWD} ; echo "\tPossible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@-$$or.out;\
	done

clean::
	-rm -f test1 test2 test3 test4 test5
	-rm -f runtest*
	-rm -f TEST*
