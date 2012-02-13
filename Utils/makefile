all: h5export vtkexport 

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

h5export: h5export.o chkopts
	@${CLINKER} -o ${VFDIR}/bin/${PETSC_ARCH}/h5export h5export.o ${VFDIR}/xdmf.o ${PETSC_LIB}
	@${RM} h5export.o

vtkexport: vtkexport.o chkopts
	@${CLINKER} -o ${VFDIR}/bin/${PETSC_ARCH}/vtkexport vtkexport.o ${PETSC_LIB}
	@${RM} vtkexport.o

clean::
	@${RM} h5export vtkexport
