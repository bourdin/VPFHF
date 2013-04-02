/*
  VFV.h
  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/

#ifndef VFV_H
#define VFV_H

extern PetscErrorCode BCVInit(BC *BC,VFPreset preset);
extern PetscErrorCode BCVUpdate(BC *BC,VFPreset preset);
extern PetscErrorCode VF_VEnergy3D(PetscReal *SurfaceEnergy,VFFields *fields,VFCtx *ctx);
extern PetscErrorCode VF_StepV(VFFields *fields,VFCtx *ctx);
/*
  These functions are not meant to be called outside of VFV, 
  but since snesV is initialized outside of VFV, I have no other choice
*/
extern PetscErrorCode VF_VIJacobian(SNES snes,Vec V,Mat *Jac,Mat *Jac1,MatStructure *str,void *user);
extern PetscErrorCode VF_VResidual(SNES snes,Vec V,Vec Func,void *user);
extern PetscErrorCode VF_VSNESMonitor(SNES snes,PetscInt its,PetscReal fnorm,void* ptr);


#endif /* VFV_H */
