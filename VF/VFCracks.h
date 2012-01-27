/*
  VFCracks.h:
    Initializers for various shaped cracks
  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
#include "VFCommon.h"

#ifndef VFCRACKS_H
#define VFCRACKS_H

extern PetscErrorCode VFPennyCrackGet(const char prefix[],VFPennyCrack *PennyCrack);
extern PetscErrorCode VFPennyCrackCreate(VFPennyCrack *PennyCrack);
extern PetscErrorCode VFPennyCrackView(VFPennyCrack *PennyCrack,PetscViewer viewer);
extern PetscErrorCode VFPennyCrackSetName(VFPennyCrack *PennyCrack,const char name[]);
extern PetscErrorCode VFDistanceToPennyCrack(PetscReal *d,PetscReal *x,VFPennyCrack *PennyCrack);

extern PetscErrorCode VFPennyCrackBuildVAT2(Vec V,VFPennyCrack *crack,VFCtx *ctx);
#endif /* VFCRACKS_H */