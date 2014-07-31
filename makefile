all: VF_Chevron

include makefile.include	


include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

.PHONY: Utils test dirs

dirs:
	@mkdir -p bin/${PETSC_ARCH}

VF_Chevron: dirs ${VFOBJ} VF_Chevron.o chkopts
	${CLINKER} -o bin/${PETSC_ARCH}/VF_Chevron VF_Chevron.o ${VFOBJ} ${PETSC_TAO_LIB}

Utils: dirs
	@echo Making all in ${VFDIR}/Utils
	@make -C ${VFDIR}/Utils

test:
	@echo Running all tests in ${VFDIR}/ValidationTests
	@make -C ${VFDIR}/ValidationTests test

clean::
	@make -C ${VFDIR}/ValidationTests clean