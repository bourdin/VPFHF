include makefile.include	
#include ../Utils/makefile.include

all: VF_Chevron

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

VF_Chevron: VF_Chevron.o ${VFOBJ} ${VFUTILSOBJ} chkopts
	-${CLINKER} -o VF_Chevron VF_Chevron.o ${VFOBJ} ${VFUTILSOBJ} ${PETSC_LIB}

