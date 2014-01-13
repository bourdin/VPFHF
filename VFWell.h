/*
  VFWell.h
  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
#include "CartFE.h"
#include "VFCommon.h"

#ifndef VFWELL_H
#define VFWELL_H

extern PetscErrorCode VFWellGet(const char prefix[],VFWell *well);
extern PetscErrorCode VFWellCreate(VFWell *well);
extern PetscErrorCode VFWellView(VFWell *well,PetscViewer viewer);
extern PetscErrorCode VFWellSetName(VFWell *well,const char name[]);
extern PetscErrorCode VFDistanceToWell(PetscReal *d,PetscReal *x,VFWell *well);
extern PetscErrorCode VFWellBuildVAT2(Vec V,VFWell *well,VFCtx *ctx);
extern PetscErrorCode VFRegDiracDeltaFunction1(Vec V,VFWell *well,VFPennyCrack *crack,VFCtx *ctx);
extern PetscErrorCode VFRegDiracDeltaFunction2(Vec RegV,VFWell *well,VFPennyCrack *crack,VFCtx *ctx,Vec V);
#endif /* VFWELL_H */