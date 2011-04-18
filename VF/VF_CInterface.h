#ifndef VFCINTERFACE_H
#define VFCINTERFACE_H

/*
extern VFCtx user;
*/
extern PetscErrorCode VFInitialize(PetscInt nx,PetscInt ny,PetscInt nz,PetscReal *dx,PetscReal *dy,PetscReal *dz);
extern PetscErrorCode vendit(PetscInt rank,PetscReal tim,PetscInt nstep,PetscInt nfout,PetscInt nfbug);
extern PetscErrorCode vtdata(PetscInt rank,PetscReal tim,PetscInt nstep,PetscInt nfout,PetscInt nfbug);
extern PetscErrorCode vperm(PetscInt rank,PetscReal *pres,PetscReal *tempr,PetscReal *pmult,PetscReal tim,PetscInt nstep,PetscInt newt,PetscInt nfout,PetscInt nfbug,PetscInt n);
extern PetscErrorCode vstdout(PetscInt rank,PetscReal tim,PetscInt nstep,PetscInt nfout,PetscInt nfbug);

#endif /* VFCINTERFACE_H */
