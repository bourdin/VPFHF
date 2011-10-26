/*
   VFFlow.c
   Generic interface to flow solvers

     (c) 2010-2011 Blaise Bourdin, LSU. bourdin@lsu.edu
                    Keita Yoshioka, Chevron ETC
*/

#ifndef VFFLOW_H
#define VFFLOW_H
#include "PetscFixes.h"

extern PetscErrorCode VFFlowTimeStep(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode BCPInit(BC *BC,VFCtx *ctx);


#endif /* VFFLOW_H */
