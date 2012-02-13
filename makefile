include makefile.include	

all: VF_Chevron #Utils

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

binarch: bin
	@mkdir -p ${VFDIR}/bin/${PETSC_ARCH}
bin:
	@mkdir ${VFDIR}/bin

VF_Chevron: binarch VF_Chevron.o ${VFOBJ} chkopts
	@${CLINKER} -o bin/${PETSC_ARCH}/VF_Chevron VF_Chevron.o ${VFOBJ} ${PETSC_LIB}

Utils: binarch
	@echo Making all in ${VFDIR}/Utils
	@make -C ${VFDIR}/Utils

test:
	@echo Running all tests in ${VFDIR}/ValidationTests
	@make -C ${VFDIR}/ValidationTests test
