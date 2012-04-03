/*
  VFU.h
  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/

#ifndef VFU_H
#define VFU_H

extern PetscErrorCode BCUInit(BC *BC,VFPreset preset);
extern PetscErrorCode BCUUpdate(BC *BC,VFPreset preset);
extern PetscErrorCode ElasticEnergyDensity3D_local(PetscReal *ElasticEnergyDensity_local,
                                                   PetscReal ****u_array,
                                                   PetscReal ***theta_array,PetscReal ***thetaRef_array,
                                                   PetscReal ***pressure_array,PetscReal ***pressureRef_array,
                                                   MatProp *matprop,PetscInt ek,PetscInt ej,PetscInt ei,
                                                   CartFE_Element3D *e);
extern PetscErrorCode ElasticEnergyDensitySphericalDeviatoric3D_local(PetscReal *ElasticEnergyDensityS_local,
                                                                      PetscReal *ElasticEnergyDensityD_local,
                                                                      PetscReal ****u_array,
                                                                      PetscReal ***theta_array,PetscReal ***thetaRef_array,
                                                                      PetscReal ***pressure_array,PetscReal ***pressureRef_array,
                                                                      MatProp *matprop,PetscInt ek,PetscInt ej,PetscInt ei,
                                                                      CartFE_Element3D *e);

extern PetscErrorCode VF_UAssembly3D(Mat K,Vec RHS,VFFields *fields,VFCtx *ctx);
extern PetscErrorCode VF_UEnergy3D(PetscReal *ElasticEnergy,PetscReal *OverbdnWork,PetscReal *PressureWork,VFFields *fields,VFCtx *ctx);

extern PetscErrorCode VF_StepU(VFFields *fields,VFCtx *ctx);
extern PetscErrorCode VF_ComputeBCU(VFFields *fields,VFCtx *ctx);
#endif /* VFU_H */
