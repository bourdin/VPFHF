include makefile.include

all: h5export vtkexport 

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

h5export: ${VFUTILSOBJ} h5export.o chkopts
	@${CLINKER} -o ${VFDIR}/bin/${PETSC_ARCH}/h5export ${VFUTILSOBJ} h5export.o ${PETSC_LIB}
	@${RM} h5export.o

vtkexport: ${VFUTILSOBJ} vtkexport.o chkopts
	@${CLINKER} -o ${VFDIR}/bin/${PETSC_ARCH}/vtkexport vtkexport.o ${PETSC_LIB}
	@${RM} vtkexport.o

clean::
	@${RM} h5export
