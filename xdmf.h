#ifndef XDMF_H
#define XDMF_H

extern PetscErrorCode XDMFuniformgridInitialize(PetscViewer viewer,PetscReal time,const char gridname[]);
extern PetscErrorCode XDMFtopologyAdd(PetscViewer viewer,PetscInt nx,PetscInt ny,PetscInt nz,const char h5filename[],const char coordname[]);
extern PetscErrorCode XDMFattributeAdd(PetscViewer viewer,PetscInt nx,PetscInt ny,PetscInt nz,PetscInt nfields,const char fieldtype [],const char location [],const char h5filename[],const char fieldname[]);
extern PetscErrorCode XDMFuniformgridFinalize(PetscViewer viewer);

extern PetscErrorCode XDMFmultistepInitialize(PetscViewer viewer);
extern PetscErrorCode XDMFmultistepAddstep(PetscViewer viewer,const char filename[]);
extern PetscErrorCode XDMFmultistepFinalize(PetscViewer viewer);
#endif /* XDMF_H */
