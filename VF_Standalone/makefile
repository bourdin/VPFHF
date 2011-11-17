include makefile.include	

all: VF_Chevron

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

VF_Chevron: VF_Chevron.o ${VFOBJ} chkopts
	@${CLINKER} -o VF_Chevron VF_Chevron.o ${VFOBJ} ${PETSC_LIB}

