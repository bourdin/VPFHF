/*
 test1.c:
 1D tension experiment
 
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */

#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFMech.h"
#include "VFFlow.h"
#include "VFCracks.h"
#include "VFPermfield.h"

#undef __FUNCT__
#define __FUNCT__ "boundaryConditionsInitialize"
PetscErrorCode boundaryConditionsInitialize(VFCtx *ctx,PetscInt orientation){
  PetscInt       i,j,c;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  for (i = 0; i < 6; i++) {
    ctx->bcV[0].face[i] = NONE;
    for (j = 0; j < 3; j++) {
      ctx->bcU[j].face[i] = NONE;
    }
  }
  for (i = 0; i < 12; i++) {
    ctx->bcV[0].edge[i] = NONE;
    for (j = 0; j < 3; j++) {
      ctx->bcU[j].edge[i] = NONE;
    }
  }
  for (i = 0; i < 8; i++) {
    ctx->bcV[0].vertex[i] = NONE;
    for (j = 0; j < 3; j++) {
      ctx->bcU[j].vertex[i] = NONE;
    }
  }
  switch (orientation) {
    case 0:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Applying traction Dirichlet conditions on faces X0 X1\n");CHKERRQ(ierr);
      for (c = 0; c < 3; c++) {
        ctx->bcU[c].face[X0] = FIXED;
        ctx->bcU[c].face[X1] = FIXED;
      }
      ctx->bcV[0].face[X0] = ONE;
      ctx->bcV[0].face[X1] = ONE;
      break;
    case 1:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Applying traction Dirichlet conditions on faces Y0 Y1\n");CHKERRQ(ierr);
      for (c = 0; c < 3; c++) {
        ctx->bcU[c].face[Y0] = FIXED;
        ctx->bcU[c].face[Y1] = FIXED;
      }
      ctx->bcV[0].face[Z0] = ONE;
      ctx->bcV[0].face[Z1] = ONE;
      break;
    case 2:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Applying traction Dirichlet conditions on faces Z0 Z1\n");CHKERRQ(ierr);
      for (c = 0; c < 3; c++) {
        ctx->bcU[c].face[Z0] = FIXED;
        ctx->bcU[c].face[Z1] = FIXED;
      }
      ctx->bcV[0].face[Z0] = ONE;
      ctx->bcV[0].face[Z1] = ONE;
      break;
    default:
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Orientation should be between 0 and 2, got %i\n",orientation);
      break;
  }
      PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "boundaryDisplacementInitialize"
PetscErrorCode boundaryDisplacementInitialize(VFCtx *ctx,PetscInt orientation,PetscReal t,PetscReal *bc,Vec BCU) {
  PetscErrorCode ierr;
  PetscReal      ****bcu_array;
  PetscInt       i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daVect,NULL,&nx,&ny,&nz,NULL,NULL,NULL,
                     NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daVect,BCU,&bcu_array);CHKERRQ(ierr);
  switch (orientation) {
    case 0:
      for (k = zs; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
          for (i = xs; i < xs+xm; i++) {
            if (i == 0) {
              bcu_array[k][j][i][0] = -t * bc[0];
              bcu_array[k][j][i][1] = -t * bc[1];
              bcu_array[k][j][i][2] = -t * bc[2];
            }
            if (i == nx-1) {
              bcu_array[k][j][i][0] = t * bc[0];
              bcu_array[k][j][i][1] = t * bc[1];
              bcu_array[k][j][i][2] = t * bc[2];
            }
          }
        }
      }
      
      break;
      
    case 1:
      for (k = zs; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
          for (i = xs; i < xs+xm; i++) {
            if (j == 0) {
              bcu_array[k][j][i][0] = -t * bc[0];
              bcu_array[k][j][i][1] = -t * bc[1];
              bcu_array[k][j][i][2] = -t * bc[2];
            }
            if (j == ny-1) {
              bcu_array[k][j][i][0] = t * bc[0];
              bcu_array[k][j][i][1] = t * bc[1];
              bcu_array[k][j][i][2] = t * bc[2];
            }
          }
        }
      }
      break;
      
    case 2:
      for (k = zs; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
          for (i = xs; i < xs+xm; i++) {
            if (k == 0) {
              bcu_array[k][j][i][0] = -t * bc[0];
              bcu_array[k][j][i][1] = -t * bc[1];
              bcu_array[k][j][i][2] = -t * bc[2];
            }
            if (k == nz-1) {
              bcu_array[k][j][i][2] = t * bc[0];
              bcu_array[k][j][i][2] = t * bc[1];
              bcu_array[k][j][i][2] = t * bc[2];
            }
          }
        }
      }
      break;
      
    default:
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Orientation should be between 0 and 2, got %i\n",orientation);
      break;
  }
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,BCU,&bcu_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  VFCtx          ctx;
  VFFields       fields;
  PetscErrorCode ierr;
  Vec            Vold;
  PetscReal      errV;
  PetscInt       altminit;
  PetscInt       nCycle = 0;
  PetscReal      boundaryDisplacement[3] = {1,0,0};
  PetscInt       nopt=3;
  PetscBool      flg;  
  PetscReal      p = 0.,crackVolume=0;

  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  
  PetscInt  orientation = 2;
  ierr = PetscOptionsGetInt(NULL,"-orientation",&orientation,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,"-pressure",&p,NULL);CHKERRQ(ierr);
  ctx.hasCrackPressure = PETSC_TRUE;
  ierr = PetscOptionsGetRealArray(NULL,"-bcu",boundaryDisplacement,&nopt,&flg);CHKERRQ(ierr);
  ctx.matprop[0].beta  = 0.;
  ctx.matprop[0].alpha = 0.;
  
  ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
  ierr = VecCopy(fields.VIrrev,fields.V);CHKERRQ(ierr);
  ierr = VecSet(fields.U,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.BCU,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
  ierr = VecSet(fields.BCU,0.0);CHKERRQ(ierr);

  ierr = boundaryConditionsInitialize(&ctx,orientation);
  
  for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
    if (ctx.maxtimestep > 1) {
      ctx.timevalue = ctx.mintimevalue + (ctx.maxtimevalue - ctx.mintimevalue) / (PetscReal) (ctx.maxtimestep-1) * (PetscReal) (ctx.timestep);
    } else {
      ctx.timevalue = ctx.maxtimevalue;
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"==== time step %i: t = %e\n",ctx.timestep,ctx.timevalue);CHKERRQ(ierr);
    ierr = boundaryDisplacementInitialize(&ctx,orientation,ctx.timevalue,boundaryDisplacement,ctx.fields->BCU);
    ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
    
    altminit = 0;
    errV = 1.e+10;
    do {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Step %i / %i\n",ctx.timestep,altminit);CHKERRQ(ierr);
      ierr = VecCopy(fields.V,Vold);
      ierr = VF_StepU(&fields,&ctx);
      ierr = VF_StepV(&fields,&ctx);

      ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
      ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"      Max. change on V: %e\n",errV);CHKERRQ(ierr);
      
      if (altminit%10 == 0) {
        ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);
      }
      altminit++;
    } while (errV >= ctx.altmintol && altminit <= ctx.altminmaxit);
    
    ierr = VolumetricCrackOpening(&crackVolume,&ctx,&fields);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Total crack opening: %f\n",crackVolume);CHKERRQ(ierr);
    ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
    
    ctx.SurfaceEnergy = 0.;
    ctx.ElasticEnergy = 0;
    ctx.InsituWork    = 0;
    ctx.PressureWork  = 0.;
    
    ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,&fields,&ctx);CHKERRQ(ierr);
    ierr = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);
    
    ctx.TotalEnergy = ctx.ElasticEnergy-ctx.InsituWork-ctx.PressureWork;
    
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic Energy:            %e\n",ctx.ElasticEnergy);CHKERRQ(ierr);
    if (ctx.hasCrackPressure) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of pressure forces:   %e\n",ctx.PressureWork);CHKERRQ(ierr);
    }
    if (ctx.hasInsitu) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of insitu stresses:   %e\n",ctx.InsituWork);CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Total Mechanical energy:  %e\n",ctx.ElasticEnergy-ctx.InsituWork-ctx.PressureWork);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface  energy:          %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);
    ctx.TotalEnergy = ctx.ElasticEnergy-ctx.InsituWork-ctx.PressureWork + ctx.SurfaceEnergy;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Total  energy:            %e\n",ctx.TotalEnergy);CHKERRQ(ierr);
		ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%d \t\t%e \t%e \t%e \t%e \t%e \t%e\n",ctx.timestep,ctx.timevalue,ctx.SurfaceEnergy,ctx.ElasticEnergy,ctx.PressureWork,ctx.InsituWork,ctx.TotalEnergy);
    
    ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);
  }
  ierr = VecDestroy(&Vold);CHKERRQ(ierr);
  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}

