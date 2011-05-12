/*
  VFU.h
  (c) 2010 Blaise Bourdin bourdin@lsu.edu
*/

#ifndef VFU_H
#define VFU_H
#include "PetscFixes.h"

extern PetscErrorCode BCUInit(BC *BC,VFPreset preset);

extern PetscErrorCode ElasticEnergyDensity3D_local(PetscReal *ElasticEnergyDensity_local,PetscReal ****u_array,PetscReal ***theta_array,PetscReal ***thetaRef_array,MatProp *matprop,PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e);
extern PetscErrorCode ElasticEnergyDensitySphericalDeviatoric3D_local(PetscReal *ElasticEnergyDensityS_local,PetscReal *ElasticEnergyDensityD_local,PetscReal ****u_array,PetscReal ***theta_array,PetscReal ***thetaRef_array,MatProp *matprop,PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e);
extern PetscErrorCode ModifiedElasticEnergyDensity3D_local(PetscReal *ElasticEnergyDensity_local,PetscReal ****u_array,PetscReal ***theta_array,PetscReal ***thetaRef_array,MatProp *matprop,PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e);

extern PetscErrorCode VF_UAssembly3D(Mat K,Vec RHS,VFFields *fields,VFCtx *ctx);
extern PetscErrorCode VF_UEnergy3D(PetscReal *ElasticEnergy,PetscReal *OverbdnWork,VFFields *fields,VFCtx *ctx);

extern PetscErrorCode VF_StepU(VFFields *fields,VFCtx *ctx);
#endif /* VFU_H */
