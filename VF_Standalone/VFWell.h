/*
  VFWell.h
  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
#include "CartFE.h"

#ifndef VFWELL_H
#define VFWELL_H

typedef struct {
  char         name[256];
  PetscReal    *top;
  PetscReal    *bottom;
  PetscReal    rate;
  BCTYPE       BCV;
} VFWell;

extern PetscErrorCode VFWellGet(const char prefix[],VFWell *well);
extern PetscErrorCode VFWellCreate(VFWell *well);
extern PetscErrorCode VFWellDestroy(VFWell *well);
extern PetscErrorCode VFWellView(VFWell *well,PetscViewer viewer);
extern PetscErrorCode VFWellSetName(VFWell *well,const char name[]);
extern PetscErrorCode VFDistanceToWell(PetscReal *d,PetscReal *x,VFWell *well);
#endif /* VFWELL_H */