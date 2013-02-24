/*
   VFHeat.h
   Generic interface to flow solvers

     (c) 2010-2012 Blaise Bourdin, LSU. bourdin@lsu.edu
                    Keita Yoshioka, Chevron ETC
*/

#ifndef VFHEAT_H
#define VFHEAT_H

extern PetscErrorCode BCTInit(BC *BCT,VFCtx *ctx);
extern PetscErrorCode HeatSolverFinalize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode HeatSolverInitialize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode BCQQInit(BC *BCq,VFCtx *ctx);
extern PetscErrorCode VFHeatTimeStep(VFCtx *ctx,VFFields *fields);



#endif /* VFFLOW_H */

