#ifndef VFCINTERFACE_H
#define VFCINTERFACE_H

typedef enum {
  FILEFORMAT_BIN,
  FILEFORMAT_HDF5
  } VFFileFormatType;
static const char *VFFileFormatName[] = {"bin","hdf5","VFFileFormatName","",0};

typedef struct {
  PetscInt          nx,ny,nz;
  DA                daVect,daScal;
  Vec               coordinates,pres,tempr,pmult;
  Vec               U,V;
  PetscReal         E,nu;
  char              prefix[PETSC_MAX_PATH_LEN];
  PetscViewer       logviewer,XDMFviewer;
  VFFileFormatType  fileformat;
} VFCtx;

extern VFCtx user;

extern PetscErrorCode VFInitialize(PetscInt nx,PetscInt ny,PetscInt nz,PetscReal *dx,PetscReal *dy,PetscReal *dz);
extern PetscErrorCode vendit(PetscInt rank,PetscReal tim,PetscInt nstep,PetscInt nfout,PetscInt nfbug);
extern PetscErrorCode vtdata(PetscInt rank,PetscReal tim,PetscInt nstep,PetscInt nfout,PetscInt nfbug);
extern PetscErrorCode vperm(PetscInt rank,PetscReal *pres,PetscReal *tempr,PetscReal *pmult,PetscReal tim,PetscInt nstep,PetscInt newt,PetscInt nfout,PetscInt nfbug,PetscInt n);
extern PetscErrorCode vstdout(PetscInt rank,PetscReal tim,PetscInt nstep,PetscInt nfout,PetscInt nfbug);

#endif /* VFCINTERFACE_H */
