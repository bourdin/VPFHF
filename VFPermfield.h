/*
   VFPermfield.h
   (c) 2012	C. Chukwudozie, LSU
*/

#ifndef VFPERMFIELD_H
#define VFPERMFIELD_H


/* 
  Rename and check if all these need to be public
*/

extern PetscErrorCode VFCheckVolumeBalance(PetscReal *ModulusVolume, PetscReal *DivVolume, PetscReal *SurfVolume, PetscReal *SumWellRate,PetscReal *SumSourceRate,PetscReal *VolStrainVolume,VFCtx *ctx, VFFields *fields);
extern PetscErrorCode TrialFunctionCompute(PetscReal *FunctionValue, VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VolumetricFunction_local(PetscReal *Function_local, PetscReal ***pmult_array, PetscReal ****u_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e);
extern PetscErrorCode CellToNodeInterpolation(Vec node_vec,Vec cell_vec,VFCtx *ctx);
extern PetscErrorCode VolumetricCrackOpening(PetscReal *CrackVolume, VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VolumetricCrackOpening3D_local(PetscReal *CrackVolume_local, PetscReal ***volcrackopening_array, PetscReal ****displ_array, PetscReal ***vfield_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e);
extern PetscErrorCode VolumetricCrackOpeningCC(PetscReal *CrackVolume, VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VolumetricCrackOpening3D_localCC(PetscReal *CrackVolume_local, PetscReal ***volcrackopening_array, PetscReal ***udotn_array, PetscReal ****u_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e);
extern PetscErrorCode PermeabilityUpDate(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode Permeabilityfield(PetscReal *COD_local, PetscReal ***volcrackopening_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement2D *e, FACE face);
extern PetscErrorCode VolumetricLeakOffRate_local(PetscReal *LeakoffRate_local, PetscReal ***volleakoffrate_array, PetscReal ****q_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e);
extern PetscErrorCode VolumetricLeakOffRate(PetscReal *LeakOffRate, VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VF_PermeabilityUpDate(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode SourceVolume_local(PetscReal *SrcVolume_local,PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e,PetscReal ***src_array,PetscReal ***v_array);
extern PetscErrorCode DivergenceVolume_local(PetscReal *DivVolume_local,PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e,PetscReal ****vel_array,PetscReal ***v_array);
extern PetscErrorCode SurfaceFluxVolume_local(PetscReal *mysurfVolumeLocal,PetscInt ek, PetscInt ej, PetscInt ei, FACE face, VFCartFEElement2D *e, PetscReal ****vel_array, PetscReal ***v_array);
extern PetscErrorCode ModulusVolume_local(PetscReal *ModVolume_local,PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e, VFFlowProp *flowpropty, PetscReal ***press_diff_array, PetscReal ***v_array);
extern PetscErrorCode VF_ComputeAverageGradient(PetscReal *grad_array, PetscInt dof, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e);
extern PetscErrorCode VF_ComputeAverageVField_local(PetscReal ***v_c_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e);
extern PetscErrorCode VolumetricStrainVolume_local(PetscReal *VolStrainVolume_local,PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e, VFFlowProp *flowpropty, PetscReal ****u_diff_array, PetscReal ***v_array);
extern PetscErrorCode VolumetricFractureWellRate_local(PetscReal *InjVolume_local, PetscReal ***regrate_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e);
extern PetscErrorCode VF_IntegrateOnBoundary(PetscReal *SumnIntegral,Vec vec, FACE face, VFCtx *ctx);

#endif 
