/*
  VFCracks.h:
    Initializers for various shaped cracks
  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
#include "VFCommon.h"

#if !defined(VFCRACKS_H)
#define VFCRACKS_H

extern PetscErrorCode VFPennyCrackGet(const char prefix[],VFPennyCrack *PennyCrack);
extern PetscErrorCode VFPennyCrackCreate(VFPennyCrack *PennyCrack);
extern PetscErrorCode VFPennyCrackView(VFPennyCrack *PennyCrack,PetscViewer viewer);
extern PetscErrorCode VFPennyCrackSetName(VFPennyCrack *PennyCrack,const char name[]);
extern PetscErrorCode VFDistanceToPennyCrack(PetscReal *d,PetscReal *x,VFPennyCrack *PennyCrack);
extern PetscErrorCode VFPennyCrackBuildVAT2(Vec V,VFPennyCrack *crack,VFCtx *ctx);


extern PetscErrorCode VFRectangularCrackGet(const char prefix[],VFRectangularCrack *RectangularCrack);
extern PetscErrorCode VFRectangularCrackCreate(VFRectangularCrack *RectangularCrack);
extern PetscErrorCode VFRectangularCrackView(VFRectangularCrack *RectangularCrack,PetscViewer viewer);
extern PetscErrorCode VFRectangularCrackSetName(VFRectangularCrack *RectangularCrack,const char name[]);
extern PetscErrorCode VFDistanceToRectangularCrack(PetscReal *d,PetscReal *x,VFRectangularCrack *RectangularCrack);
extern PetscErrorCode VFRectangularCrackBuildVAT2(Vec V,VFRectangularCrack *crack,VFCtx *ctx);

#endif /* VFCRACKS_H */