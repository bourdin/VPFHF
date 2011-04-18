include makefile.include

all: h5export vtkexport 

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
	
h5export: ${VFUTILSOBJ} h5export.o chkopts
	@${CLINKER} -o ${VFDIR}/bin/h5export ${VFUTILSOBJ} h5export.o ${PETSC_LIB}
	
vtkexport: ${VFUTILSOBJ} vtkexport.o chkopts
	@${CLINKER} -o ${VFDIR}/bin/vtkexport vtkexport.o ${PETSC_LIB}

clean::
	@${RM} h5export
