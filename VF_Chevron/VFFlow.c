#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"

#undef __FUNCT__
#define __FUNCT__ "VFFlowTimeStep"


extern PetscErrorCode VFFlowTimeStep(VFCtx *ctx,VFFields *fields)
{
  PetscErrorCode     ierr;
//  KSPConvergedReason reason;  

  PetscFunctionBegin;
//  ierr = VF_FFAssembly3D(...);CHKERRQ(ierr);
//  ierr = KSPSolve(...);CHKERRQ(ierr);
  
  ierr = VF_FakeFlow(ctx,fields,&ctx->resprop);  

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_FakeFlow"

/* 
   Fake flow solver for VF_Chevron.c test
*/

extern PetscErrorCode VF_FakeFlow(VFCtx *ctx, VFFields *fields,ResProp *resprop)
{
  PetscErrorCode ierr;
  PetscReal      time;
  PetscReal      pres,temp;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       i,j,k;
  Vec            theta_localVec;
  Vec            pressure_localVec;
  PetscReal      ***theta_array;
  PetscReal      ***pressure_array;
  PetscReal      ****coords_array;
  PetscReal      x,y,z,r;
  PetscReal      Tinit,Pinit,P1,P2,T1,T2;

  PetscFunctionBegin;

  ierr = DAGetInfo(ctx->daVect,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DAGetCorners(ctx->daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  /*
     Get coordinate, temperature array, and pressure array
  */
  ierr = DAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DAVecGetArray(ctx->daScal,fields->theta,&theta_array);CHKERRQ(ierr);
  ierr = DAVecGetArray(ctx->daScal,fields->pressure,&pressure_array);CHKERRQ(ierr);

  time = ctx->timevalue;
  Tinit = resprop->Tinit;
  Pinit = resprop->Pinit;
  P1 = 30.;
  P2 = 1.e6;
  T1 = 2.;
  T2 = 1.e6;

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        x = coords_array[k][j][i][0]+0.5;
        y = coords_array[k][j][i][1]+0.5;
        z = coords_array[k][j][i][2]+0.5;
        r = sqrt(x*x+y*y);
        theta_array[k][j][i] = Tinit - T1*log(T2*time/(r*r));
        pressure_array[k][j][i] = Pinit + P1*log(P2*time/(r*r));
      }
    }
  } 
  
  /*
     Cleanup
  */
  ierr = DAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(ctx->daScal,fields->theta,&theta_array);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(ctx->daScal,fields->pressure,&pressure_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


  
