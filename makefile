include ${VFDIR}/Utils/makefile.include
include ${VFDIR}/VF/makefile.include
include ${VFDIR}/FakeGMRS/makefile.include
include ${VFDIR}/FakeVF/makefile.include

all: Utils FakeGMRS_FakeVF FakeGMRS_VF

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules


Utils: ${VFUTILSOBJ}
	@echo Making all in ${VFDIR}/Utils
	@make -C Utils

FakeVF: ${FAKEVFOBJ}
	@echo Making all in ${VFDIR}/FakeVF
	@make -C FakeVF
	
VF: ${VFOBJ}
	@echo Making all in ${VFDIR}/VF
	@make -C VF

FakeGMRS: ${FAKEGMRSOBJ}
	@echo Making all in ${VFDIR}/FakeGMRS
	@make -C FakeGMRS

FakeGMRS_FakeVF: Utils FakeVF FakeGMRS chkopts
	@${FLINKER} -o bin/FakeGMRS_FakeVF ${VFUTILSOBJ} ${FAKEVFOBJ} ${FAKEGMRSOBJ} ${PETSC_LIB}
	
FakeGMRS_VF: Utils VF FakeGMRS chkopts
	@${FLINKER} -o bin/FakeGMRS_VF ${VFUTILSOBJ} ${VFOBJ} ${FAKEGMRSOBJ} ${PETSC_LIB}
	
clean::
	@make -C Utils clean
	@make -C FakeVF clean
	@make -C FakeGMRS clean
	@make -C VF clean
	@${RM} bin/FakeGMRS_FakeVF bin/FakeGMRS_VF bin/h5export

