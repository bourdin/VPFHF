/*
   VFPermfield.h
   (c) 2012	C. Chukwudozie, LSU
*/

#ifndef VFPERMFIELD_H
#define VFPERMFIELD_H


/* 
  Rename and check if all these need to be public
*/
extern PetscErrorCode TrialFunctionCompute(PetscReal *FunctionValue, VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VolumetricFunction_local(PetscReal *Function_local, PetscReal ***pmult_array, PetscReal ****u_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFE_Element3D *e);
extern PetscErrorCode CellToNodeInterpolation(Vec node_vec,Vec cell_vec,VFCtx *ctx);
extern PetscErrorCode VolumetricCrackOpening(PetscReal *CrackVolume, VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VolumetricCrackOpening3D_local(PetscReal *CrackVolume_local, PetscReal ***volcrackopening_array, PetscReal ****displ_array, PetscReal ***vfield_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFE_Element3D *e);
extern PetscErrorCode PermeabilityUpDate(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode Permeabilityfield(PetscReal *COD_local, PetscReal ***volcrackopening_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFE_Element2D *e, FACE face);
extern PetscErrorCode VolumetricLeakOffRate_local(PetscReal *LeakoffRate_local, PetscReal ***volleakoffrate_array, PetscReal ****q_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFE_Element3D *e);
extern PetscErrorCode VolumetricLeakOffRate(PetscReal *LeakOffRate, VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VF_PermeabilityUpDate(VFCtx *ctx, VFFields *fields);
#endif 
