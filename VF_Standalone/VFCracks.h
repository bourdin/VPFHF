/*
  VFCracks.h:
    Initializers for various shaped cracks
  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/

#ifndef VFCRACKS_H
#define VFCRACKS_H

typedef struct {
  char          name[256];
  PetscReal    *center;
  PetscReal     r,phi,theta;
} VFPennyCrack;

extern PetscErrorCode VFPennyCrackGet(const char prefix[],VFPennyCrack *PennyCrack);
extern PetscErrorCode VFPennyCrackCreate(VFPennyCrack *PennyCrack);
extern PetscErrorCode VFPennyCrackDestroy(VFPennyCrack *PennyCrack);
extern PetscErrorCode VFPennyCrackView(VFPennyCrack *PennyCrack,PetscViewer viewer);
extern PetscErrorCode VFPennyCrackSetName(VFPennyCrack *PennyCrack,const char name[]);
extern PetscErrorCode VFDistanceToPennyCrack(PetscReal *d,PetscReal *x,VFPennyCrack *PennyCrack);
#endif /* VFCRACKS_H */