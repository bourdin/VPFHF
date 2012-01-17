/*
  VFV.h
  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/

#ifndef VFV_H
#define VFV_H
#include "PetscFixes.h"

extern PetscErrorCode BCVInit(BC *BC,VFPreset preset);
extern PetscErrorCode BCVUpdate(BC *BC,VFPreset preset);

//extern PetscErrorCode VF_MatVAT2Surface3D_local(PetscReal *Mat_local,MatProp *matprop,VFProp *vfprop,CartFE_Element3D *e);
// extern PetscErrorCode VF_MatVCoupling3D_local(PetscReal *Mat_local,PetscReal ****U_array,PetscReal ***theta_array,MatProp *matprop,VFProp *vfprop,CartFE_Element3D *e);

//extern PetscErrorCode VF_RHSVAT2Surface3D_local(PetscReal *RHS_local,MatProp *matprop,VFProp *vfprop,CartFE_Element3D *e);

extern PetscErrorCode VF_IrrevEQ(Mat K,Vec RHS,Vec V,VFProp *vfprop,VFCtx *ctx);

extern PetscErrorCode VF_VAssembly3D(Mat K,Vec RHS,MatProp *matprop,VFProp *vfprop,VFCtx *ctx);
extern PetscErrorCode VF_VEnergy3D(PetscReal *SurfaceEnergy,VFFields *fields,VFCtx *ctx);

extern PetscErrorCode VF_StepV(VFFields *fields,VFCtx *ctx);

extern PetscReal DistanceToDisk(PetscReal *dist,PetscReal *x,PetscReal *xc,PetscReal phi,PetscReal theta,PetscReal r);
#endif /* VFV_H */
