include ${VFDIR}/Utils/makefile.include
include ${VFDIR}/VF/makefile.include
include ${VFDIR}/FakeGMRS/makefile.include
include ${VFDIR}/FakeGMRS2/makefile.include
include ${VFDIR}/FakeVF/makefile.include

all: Utils FakeGMRS_FakeVF FakeGMRS_VF FakeGMRS2_FakeVF FakeGMRS2_VF bin

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

bin:
	@mkdir -p ${VFDIR}/bin/${PETSC_ARCH}
	
Utils: bin ${VFUTILSOBJ}
	@echo Making all in ${VFDIR}/Utils
	@make -C Utils

FakeVF: bin ${FAKEVFOBJ}
	@echo Making all in ${VFDIR}/FakeVF
	@make -C FakeVF
	
VF: bin ${VFOBJ}
	@echo Making all in ${VFDIR}/VF
	@make -C VF

FakeGMRS: bin ${FAKEGMRSOBJ}
	@echo Making all in ${VFDIR}/FakeGMRS
	@make -C FakeGMRS

FakeGMRS_FakeVF: bin Utils FakeVF FakeGMRS chkopts
	@${FLINKER} -o bin/${PETSC_ARCH}/FakeGMRS_FakeVF ${VFUTILSOBJ} ${FAKEVFOBJ} ${FAKEGMRSOBJ} ${PETSC_LIB}
	
FakeGMRS_VF: bin Utils VF FakeGMRS chkopts
	@${FLINKER} -o bin/${PETSC_ARCH}/FakeGMRS_VF ${VFUTILSOBJ} ${VFOBJ} ${FAKEGMRSOBJ} ${PETSC_LIB}
	
FakeGMRS2: bin ${FAKEGMRS2OBJ}
	@echo Making all in ${VFDIR}/FakeGMRS2
	@make -C FakeGMRS2

FakeGMRS2_FakeVF: bin Utils FakeVF FakeGMRS2 chkopts
	@${FLINKER} -o bin/${PETSC_ARCH}/FakeGMRS2_FakeVF ${VFUTILSOBJ} ${FAKEVFOBJ} ${FAKEGMRS2OBJ} ${PETSC_LIB}
	
FakeGMRS2_VF: bin Utils VF FakeGMRS2 chkopts
	@${FLINKER} -o bin/${PETSC_ARCH}/FakeGMRS2_VF ${VFUTILSOBJ} ${VFOBJ} ${FAKEGMRS2OBJ} ${PETSC_LIB}
	
clean::
	@make -C Utils clean
	@make -C FakeVF clean
	@make -C FakeGMRS clean
	@make -C FakeGMRS2 clean
	@make -C VF clean
	@${RM} -f bin/${PETSC_ARCH}/FakeGMRS_FakeVF bin/${PETSC_ARCH}/FakeGMRS_VF
	@${RM} -f bin/${PETSC_ARCH}/FakeGMRS2_FakeVF bin/${PETSC_ARCH}/FakeGMRS2_VF
	@${RM} -f bin/${PETSC_ARCH}/h5export bin/${PETSC_ARCH}/vtkexport

