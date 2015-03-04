/*
  VFWell.h
  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
#include "VFCartFE.h"
#include "VFCommon.h"

#ifndef VFWELL_H
#define VFWELL_H

extern PetscErrorCode VFWellGet(const char prefix[],VFWell *well);
extern PetscErrorCode VFWellCreate(VFWell *well);
extern PetscErrorCode VFWellView(VFWell *well,PetscViewer viewer);
extern PetscErrorCode VFWellSetName(VFWell *well,const char name[]);
extern PetscErrorCode VFDistanceToWell(PetscReal *d,PetscReal *x,VFWell *well);
extern PetscErrorCode VFWellBuildVAT2(Vec V,VFWell *well,VFCtx *ctx);
extern PetscErrorCode VFWellBuildVAT1(Vec V,VFWell *well,VFCtx *ctx);
extern PetscErrorCode VFRegDiracDeltaFunction(Vec RegV,VFWell *well,VFPennyCrack *crack,VFRectangularCrack *rcrack,VFCtx *ctx, Vec V);
extern PetscErrorCode VFRegRateScalingFactor(PetscReal *InjectedVolume, Vec Rate, Vec V, VFCtx *ctx);
#endif /* VFWELL_H */