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

clean::
	-rm -f test1