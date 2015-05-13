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
extern PetscErrorCode VolumetricFunction_local(PetscReal *Function_local, PetscReal ***pmult_array, PetscReal ****u_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e);
extern PetscErrorCode CellToNodeInterpolation(Vec node_vec,Vec cell_vec,VFCtx *ctx);
extern PetscErrorCode VolumetricCrackOpening(PetscReal *CrackVolume, VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VolumetricCrackOpening3D_local(PetscReal *CrackVolume_local, PetscReal ***volcrackopening_array, PetscReal ****displ_array, PetscReal ***vfield_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e);
extern PetscErrorCode VolumetricCrackOpening3D_localCC(PetscReal *CrackVolume_local, PetscReal ***volcrackopening_array, PetscReal ***udotn_array, PetscReal ****u_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e);
extern PetscErrorCode VolumetricLeakOffRate_local(PetscReal *LeakoffRate_local, PetscReal ***volleakoffrate_array, PetscReal ****q_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e);
extern PetscErrorCode VolumetricLeakOffRate(PetscReal *LeakOffRate, VFCtx *ctx, VFFields *fields);
extern PetscErrorCode SourceVolume_local(PetscReal *SrcVolume_local,PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e,PetscReal ***src_array,PetscReal ***v_array);
extern PetscErrorCode DivergenceVolume_local(PetscReal *DivVolume_local,PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e,PetscReal ****vel_array,PetscReal ***v_array);
extern PetscErrorCode SurfaceFluxVolume_local(PetscReal *mysurfVolumeLocal,PetscInt ek, PetscInt ej, PetscInt ei, FACE face, VFCartFEElement2D *e, PetscReal ****vel_array, PetscReal ***v_array);
extern PetscErrorCode ModulusVolume_local(PetscReal *ModVolume_local,PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e, PetscReal m_inv, PetscReal ***press_diff_array, PetscReal ***v_array);
extern PetscErrorCode VF_ComputeCellCenterVGradient_local(PetscReal *ave_v, PetscReal *grad_cc, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFEElement3D *s);
extern PetscErrorCode VF_ComputeAverageVField_local(PetscReal ***v_c_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e);
extern PetscErrorCode VolumetricStrainVolume_local(PetscReal *VolStrainVolume_local,PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e, VFMatProp *matprop, PetscReal ****u_diff_array, PetscReal ***v_array);
extern PetscErrorCode VolumetricFractureWellRate_local(PetscReal *InjVolume_local, PetscReal ***regrate_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e);
extern PetscErrorCode VF_IntegrateOnBoundary(PetscReal *SumnIntegral,Vec vec, FACE face, VFCtx *ctx);
extern PetscErrorCode VF_FastFourierTransforms(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VolumetricFractureWellRate(PetscReal *InjectedVolume, VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VF_ComputeRegularizedFracturePressure_local(PetscReal ***press_c_array,  PetscReal ***press_array, PetscReal ****u_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e);
extern PetscErrorCode VF_ComputeRegularizedFracturePressure(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode IntegrateUcdotGradVlocal(PetscReal *w_ave, PetscReal *cod, VFCartFEElement1D *e);
extern PetscErrorCode UpdateFractureWidth(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode ComputeUcdotGradVlocal1(PetscReal *cod, PetscReal *v_elem, PetscReal *n_elem, PetscReal ****u_array, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFEElement3D *s);
extern PetscErrorCode VolumeFromWidth(PetscReal *CrackVolume, VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VolumeFromWidth_local(PetscReal *CrackVolume_local, PetscReal w, PetscReal ***v_array, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e);

#endif
