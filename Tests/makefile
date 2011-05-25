include ../Utils/makefile.include

all: testwell

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

testwell: testwell.o ../VF/VFWell.o chkopts
	-${CLINKER} -o testwell testwell.o ../VF/VFWell.o ../Utils/xdmf.o ${PETSC_LIB}
