/*
  VFV.h
  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/

#ifndef VFV_H
#define VFV_H

extern PetscErrorCode BCVInit(BC *BC,VFPreset preset);
extern PetscErrorCode BCVUpdate(BC *BC,VFPreset preset);

extern PetscErrorCode VF_IrrevEQ(Mat K,Vec RHS,Vec V,VFProp *vfprop,VFCtx *ctx);

extern PetscErrorCode VF_VAssembly3D(Mat K,Vec RHS,MatProp *matprop,VFProp *vfprop,VFCtx *ctx);
extern PetscErrorCode VF_VEnergy3D(PetscReal *SurfaceEnergy,VFFields *fields,VFCtx *ctx);

extern PetscErrorCode VF_StepV(VFFields *fields,VFCtx *ctx);
#endif /* VFV_H */
