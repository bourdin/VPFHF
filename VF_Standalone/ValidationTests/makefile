include ../makefile.include	
CFLAGS='-I..'

all: test1
runall: runtest1

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules


test1: test1.o ${VFOBJ} chkopts
	@${CLINKER} -o test1 test1.o ${VFOBJ} ${PETSC_LIB}

clean::
	-rm -f test1