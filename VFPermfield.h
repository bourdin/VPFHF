/*
   VFPermfield.h
   (c) 2012	C. Chukwudozie, LSU
*/

#ifndef VFPERMFIELD_H
#define VFPERMFIELD_H


/* 
  Rename and check if all these need to be public
*/
extern PetscErrorCode CellToNodeInterpolation(Vec cell_vec, VFCtx *ctx);
extern PetscErrorCode VolumetricCrackOpening3D_local(PetscReal *CrackVolume_local, PetscReal ***volcrackopening_array, PetscReal ****displ_array, PetscReal ***vfield_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFE_Element3D *e);
extern PetscErrorCode VolumetricCrackOpening(PetscReal *CrackVolume, VFCtx *ctx, VFFields *fields);
extern PetscErrorCode PostProcessing(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode Permeabilityfield(PetscReal *COD_local, PetscReal ***volcrackopening_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFE_Element2D *e, FACE face);

#endif 
