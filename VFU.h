/*
  VFU.h
  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
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

extern PetscErrorCode VF_ComputeBCU(VFFields *fields,VFCtx *ctx);
extern PetscErrorCode VF_UIJacobian(SNES snes,Vec U,Mat *Jac,Mat *Jac1,MatStructure *str,void *user);
extern PetscErrorCode VF_UResidual(SNES snes,Vec U,Vec Func,void *user);
extern PetscErrorCode VF_USNESMonitor(SNES snes,PetscInt U_its,PetscReal fnorm,void* ptr);
extern PetscErrorCode VF_StepU(VFFields *fields,VFCtx *ctx);



#endif /* VFU_H */
