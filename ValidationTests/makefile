include ../makefile.include	
CFLAGS='-I..'

all: test1 test2 test3 test4 test5 test6 test7 test8 test9 test10 test11 test12 test13 test14 test15 test16 test17 test18 test19 test20 test21
test: runtest1 runtest2 runtest3 runtest4 runtest5 runtest6 runtest7 runtest8 runtest9 runtest10 runtest11 runtest12 test13 runtest14 runtest15 runtest16 runtest17 runtest18 runtest19 runtest20 runtest21

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules


test1: test1.o ${VFOBJ} chkopts
	-@${CLINKER} -o test1 test1.o ${VFOBJ} ${PETSC_LIB}

test2: test2.o ${VFOBJ} chkopts
	-@${CLINKER} -o test2 test2.o ${VFOBJ} ${PETSC_LIB}

test3: test3.o ${VFOBJ} chkopts
	-@${CLINKER} -o test3 test3.o ${VFOBJ} ${PETSC_LIB}

test4: test4.o ${VFOBJ} chkopts
	-@${CLINKER} -o test4 test4.o ${VFOBJ} ${PETSC_LIB}

test5: test5.o ${VFOBJ} chkopts
	-@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test6: test6.o ${VFOBJ} chkopts
	-@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test7: test7.o ${VFOBJ} chkopts
	-@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test8: test8.o ${VFOBJ} chkopts
	-@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test9: test9.o ${VFOBJ} chkopts
	-@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test10: test10.o ${VFOBJ} chkopts
	-@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test11: test11.o ${VFOBJ} chkopts
	-@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test12: test12.o ${VFOBJ} chkopts
	-@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test13: test13.o ${VFOBJ} chkopts
	-@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test14: test14.o ${VFOBJ} chkopts
	-@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test15: test15.o ${VFOBJ} chkopts
	-@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test16: test16.o ${VFOBJ} chkopts
	-@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test17: test17.o ${VFOBJ} chkopts
	-@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test18: test18.o ${VFOBJ} chkopts
	-@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test19: test19.o ${VFOBJ} chkopts
	-@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test20: test20.o ${VFOBJ} chkopts
	-@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}

test21: test21.o ${VFOBJ} chkopts
	-@${CLINKER} -o $@ $@.o ${VFOBJ} ${PETSC_LIB}
	
runtest1: test1
	@echo "Running " $@
	-@ARGS='-n 41,41,81 -l 1,1,2 -radius .4 -orientation 3 -nu 0 -epsilon 0.05 -eta 1e-8';\
	echo ${MPIEXEC} -n 2 ./$< $${ARGS} -p $@;\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@ | grep -v -e '      \#' > $@.out 2>&1;\
	if (${DIFF} -B results/$@.out $@.out) then echo "\t$< succeeded\n========================================="; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest2: test2
	@echo "Running " $@
	-@ARGS='-n 41,81,6 -l 1,2,.25 -length .4 -orientation 3 -nu 0 -epsilon 0.05 -eta 1e-8';\
	echo ${MPIEXEC} -n 2 ./$< $${ARGS} -p $@;\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@ | grep -v -e '      \#' > $@.out 2>&1;\
	if (${DIFF} -B results/$@.out $@.out) then echo "\t$< succeeded\n========================================="; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest2-large: test2
	@echo "Running " $@
	-@ARGS='-n 161,321,6 -l 4,8,.25 -length .4 -orientation 3 -nu 0 -epsilon 0.05 -eta 1e-8';\
	echo ${MPIEXEC} -n 2 ./$< $${ARGS} -p $@;\
	${MPIEXEC} -n 4 ./$< $${ARGS} -p $@ | grep -v -e '      \#' > $@.out 2>&1;\
	if (${DIFF} -B results/$@.out $@.out) then echo "\t$< succeeded\n========================================="; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest3: test3
	@echo "Running" $@
	-@ARGS='-n 41,41,6 -l 2,2,.25 -length .4 -orientation 1 -nu 0 -eta 1e-8' ;\
	echo ${MPIEXEC} -n 2 ./$< $${ARGS} -p $@;
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@ | grep -v -e '      \#' > $@.out 2>&1;\
	if (${DIFF} -B results/$@.out $@.out) then echo "\t$< succeeded\n========================================="; \
	else echo ${PWD} ; echo "\tPossible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest4: test4
	@echo "Running " $@
	-@ARGS='-n 81,81,6 -l 2,2,.25 -length .4 -orientation 3 -nu 0 -epsilon ${EPSILON} -eta 1e-8';\
	echo ${MPIEXEC} -n 2 ./$< $${ARGS} -p $@;\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@ | grep -v -e '      \#' > $@.out 2>&1;\
	if (${DIFF} -B results/$@.out $@.out) then echo "\t$< succeeded\n========================================="; \
	else echo ${PWD} ; echo "\tPossible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest5: test5
	@echo "Running " $@
	-@for or in 0 1 2 3 4 5; do \
	ARGS="-n 21,21,21 -l 1,1,1 -nu 0 -epsilon .05 -eta 1e-8 -orientation $$or";\
	echo "Running " ${MPIEXEC} -n 2 ./$< $${ARGS} -p $@-$$or ;\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@-$$or | grep -v -e '      \#' > $@-$$or.out 2>&1;\
	if (${DIFF} -B results/$@-$$or.out $@-$$or.out) then echo "\t$<-$$or succeeded\n========================================="; \
	else echo ${PWD} ; echo "\tPossible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@-$$or.out;\
	done

runtest6: test6
	@echo "Running " $@
	-@for or in 0 1 2; do \
	ARGS="-n 21,21,21 -l 1,1,1 -nu 0.3 -orientation $$or";\
	echo ${MPIEXEC} -n 2 ./$< $${ARGS} -p $@-$$or;\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@-$$or | grep -v -e '      \#' > $@-$$or.out 2>&1;\
	if (${DIFF} -B results/$@-$$or.out $@-$$or.out) then echo "\t$<-$$or succeeded\n========================================="; \
	else echo ${PWD} ; echo "\tPossible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@-$$or.out;\
	done
	
runtest7: test7
	@echo "Running " $@
	-@for or in 0 1 2; do \
	ARGS="-n 21,21,21 -l 1,1,1 -nu 0.3 -orientation $$or";\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@-$$or | grep -v -e '      \#' > $@-$$or.out 2>&1;\
	if (${DIFF} -B  results/$@-$$or.out $@-$$or.out) then echo "\t$<-$$or succeeded\n========================================="; \
	else echo ${PWD} ; echo "\tPossible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@-$$or.out;\
	done

runtest8: test8
	@echo "Running " $@
	-@ARGS='-n 41,41,81 -l 1,1,2 -radius .75 -orientation 3 -nu 0 -epsilon 0.05 -eta 1e-8';\
	echo ${MPIEXEC} -n 2 ./$< $${ARGS} -p $@;\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@ | grep -v -e '      \#' > $@.out 2>&1;\
	if (${DIFF} -B results/$@.out $@.out) then true; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest8-large: test8
	@echo "Running " $@
	-@ARGS='-n 81,81,161 -l 1,1,2 -radius .75 -orientation 3 -nu 0.3 -epsilon 0.025 -eta 1e-8';\
	echo ${MPIEXEC} -n 8 ./$< $${ARGS} -p $@;\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@ | grep -v -e '      \#' > $@.out 2>&1;\
	if (${DIFF} -B results/$@.out $@.out) then true; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest9: test9
	@echo "Running " $@
	-@ARGS='-n 41,41,41 -epsilon .05 -nc 2 -c0_center .5,.5,.5 -c0_r .5 -c0_theta 45 -c0_phi 45 -c0_name crack1 -c1_r 0.5';\
	echo ${MPIEXEC} -n 2 ./$< $${ARGS} -p $@;\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@ | grep -v -e '      \#' > $@.out 2>&1;\
	if (${DIFF} -B results/$@.out $@.out) then true; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest10: test10
	@echo "Running " $@
	-@for or in 0; do \
	ARGS="-n 21,21,21 -l 1,1,1 -nu 0.3 -length .2 -orientation $$or";\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@-$$or | grep -v -e '      \#' > $@-$$or.out 2>&1;\
	if (${DIFF} -B  results/$@-$$or.out $@-$$or.out) then echo "\t$<-$$or succeeded\n========================================="; \
	else echo ${PWD} ; echo "\tPossible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@-$$or.out;\
	done

runtest11: test11
	@echo "Running " $@
	-@ARGS='-n 41,41,6 -l 2,2,.25 -length .4 -nu 0 -epsilon .05 -eta 1e-8 -pinc 1e-6 -maxtimestep 1000 -Gc 4.e-9';\
	echo ${MPIEXEC} -n 8 ./$< $${ARGS} -p $@;\
	${MPIEXEC} -n 8 ./$< $${ARGS} -p $@ | grep -v -e '      \#' > $@.out 2>&1;\
	if (${DIFF} -B results/$@.out $@.out) then true; \
	else echo ${PWD} ; echo "\tPossible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@-$$or.out;\

runtest12: test12
	@echo "Running " $@
	-@ARGS='-n 21,21,21 -l 1,1,1 -height .1 -length .1 -orientation 1 -nu 0 -epsilon 0.05 -eta 1e-8';\
	echo ${MPIEXEC} -n 2 ./$< $${ARGS} -p $@;\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@ | grep -v -e '      \#' > $@.out 2>&1;\
	if (${DIFF} -B results/$@.out $@.out) then true; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest13: test13
	@echo "Running " $@
	-@ARGS='-n 21,21,21 -l 1,1,1 -length .4 -orientation 1 -nu 0 -epsilon 0.05 -eta 1e-8';\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@ > $@.out 2>&1;\
	if (${DIFF} -I '^#' results/$@.out $@.out) then true; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest14: test14
	@echo "Running " $@
	-@ARGS='-n 21,21,21 -l 1,1,1 -height .1 -length .1 -orientation 1 -nu 0 -epsilon 0.05 -eta 1e-8';\
	echo ${MPIEXEC} -n 2 ./$< $${ARGS} -p $@;\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@ | grep -v -e '      \#' > $@.out 2>&1;\
	if (${DIFF} -B results/$@.out $@.out) then true; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest15: test15
	@echo "Running " $@
	-@ARGS='-n 21,21,21 -l 1,1,1 -height .1 -length .1 -orientation 1 -nu 0 -epsilon 0.05 -eta 1e-8';\
	echo ${MPIEXEC} -n 2 ./$< $${ARGS} -p $@;\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@ | grep -v -e '      \#' > $@.out 2>&1;\
	if (${DIFF} -B results/$@.out $@.out) then true; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest16: test16
	@echo "Running " $@
	-@ARGS='-n 21,21,21 -l 1,1,1 -radius 0.1 -orientation 1 -nu 0.2 -epsilon 0.05 -eta 1e-8 -maxtimestep 1000 -rate 1e-10 -Gc 1.22e-13 -pinit 1e-6';\
	echo ${MPIEXEC} -n 2 ./$< $${ARGS} -p $@;\
	${MPIEXEC} -n 8 ./$< $${ARGS} -p $@ | grep -v -e '      \#' > $@.out 2>&1;\
	if (${DIFF} -B results/$@.out $@.out) then true; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest17: test17
	@echo "Running " $@
	-@ARGS='-n 21,21,21 -l 1,1,1 -length 1. -length1 .5 -orientation 2 -nu 0 -epsilon 0.05 -eta 1e-8';\
	echo ${MPIEXEC} -n 2 ./$< $${ARGS} -p $@;\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@ | grep -v -e '      \#' > $@.out 2>&1;\
	if (${DIFF} -B results/$@.out $@.out) then true; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest18: test18
	@echo "Running " $@
	-@ARGS='-n 21,21,21 -l 1,1,1 -length 1. -length1 .5 -orientation 2 -nu 0 -epsilon 0.05 -eta 1e-8';\
	echo ${MPIEXEC} -n 2 ./$< $${ARGS} -p $@;\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@ | grep -v -e '      \#' > $@.out 2>&1;\
	if (${DIFF} -B results/$@.out $@.out) then true; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runtest19: test19
	@echo "Running " $@
	-@for or in 0; do \
	ARGS="-n 21,21,21 -l 1,1,1 -nu 0.3 -length .2 -orientation $$or";\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@-$$or | grep -v -e '      \#' > $@-$$or.out 2>&1;\
	if (${DIFF} -B  results/$@-$$or.out $@-$$or.out) then echo "\t$<-$$or succeeded\n========================================="; \
	else echo ${PWD} ; echo "\tPossible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@-$$or.out;\
	done

runtest20: test20
	@echo "Running " $@
	-@for or in 0; do \
	ARGS="-n 21,21,21 -l 1,1,1 -nu 0.3 -length .2 -orientation $$or";\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@-$$or | grep -v -e '      \#' > $@-$$or.out 2>&1;\
	if (${DIFF} -B  results/$@-$$or.out $@-$$or.out) then echo "\t$<-$$or succeeded\n========================================="; \
	else echo ${PWD} ; echo "\tPossible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@-$$or.out;\
	done

runtest21: test21
	@echo "Running " $@
	-@for or in 0; do \
	ARGS="-n 21,21,21 -l 1,1,1 -nu 0.3 -length .2 -orientation $$or";\
	${MPIEXEC} -n 2 ./$< $${ARGS} -p $@-$$or | grep -v -e '      \#' > $@-$$or.out 2>&1;\
	if (${DIFF} -B  results/$@-$$or.out $@-$$or.out) then echo "\t$<-$$or succeeded\n========================================="; \
	else echo ${PWD} ; echo "\tPossible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@-$$or.out;\
	done

clean::
	-rm -f test1 test2 test3 test4 test5 test6 test7 test8 test9 test10 test11 test12 test13 test14 test15 test16 test17 test18 test19 test20 test21 pressure.txt
	-rm -f runtest*
	-rm -f TEST*
