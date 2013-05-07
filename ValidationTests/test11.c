/*
 test11.c:
 Validate crack opening computation: Solves for displacement field for sliding and pulling experiment in the presence of crack that divided material into two equal halves.

 (c) 2010-2012 chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
 ./test11 -n 101,101,2 -l 1,1,0.01 -epsilon 0.04

 
 ./test11 -n 101,101,2 -l 1,1,0.01 -npc 1 -pc0_r 0.6 -pc0_center 0.5,0.5,0.005 -pc0_thickness 0.04 -epsilon 0.02 -pc0_theta 0 -pc0_phi 90 -orientation 0
 
 thickness determines result more than epsilon. A thickness of atleast 3*resolution
 
 
 ./test11 -n 101,101,2 -l 1,1,0.01 -npc 1 -pc0_r 0.8 -pc0_center 0.5,0.5,0.005 -pc0_thickness 0.04 -epsilon 0.02 -pc0_theta 45 -pc0_phi 90 -orientation 0
 ./test11 -n 201,201,2 -l 1,1,0.01 -npc 1 -pc0_r 0.8 -pc0_center 0.5,0.5,0.005 -pc0_thickness 0.02 -epsilon 0.01 -pc0_theta 45 -pc0_phi 90 -orientation 3
 ./test11 -n 101,101,2 -l 1,1,0.01 -npc 1 -pc0_r 0.8 -pc0_center 0.5,0.5,0.005 -pc0_thickness 0.04 -epsilon 0.02 -pc0_theta 0 -pc0_phi 90 -orientation 0
 
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFV.h"
#include "VFU.h"
#include "VFFlow.h"

VFCtx    ctx;
VFFields fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  VFCtx                 ctx;
  VFFields              fields;
  PetscErrorCode ierr;

  PetscReal             length      = .2;
  PetscInt              orientation = 0;
  PetscInt              nopts       = 3;
  PetscInt              i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal             ****coords_array;
  PetscReal             ****bcu_array;
  PetscReal             BBmin[3],BBmax[3];
  PetscReal             ElasticEnergy = 0;
  PetscReal             InsituWork    = 0;
  PetscReal             SurfaceEnergy = 0;
  char                  filename[FILENAME_MAX];
  PetscReal             bc = .2;
  PetscReal             ***v_array;
  PetscReal             lx,ly,lz;
  PetscReal             ***pmult_array;
  PetscReal             ****vfperm_array;
	PetscInt              altminit=1;
	Vec                   Vold;
	PetscReal             errV=1e+10;

  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);

  ierr = PetscOptionsGetInt(PETSC_NULL,"-orientation",&orientation,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);

  lz = BBmax[2];
  ly = BBmax[1];
  lx = BBmax[0];

  printf("\nBounding box: lx = %f\tly = %f\tlz = %f\n",lx,ly,lz);

  ctx.matprop[0].beta  = 0.;
  ctx.matprop[0].alpha = 0.;

  ctx.timestep  = 1;
  ctx.timevalue = 1.;

  ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
//  ierr = VecSet(fields.V,1.0);CHKERRQ(ierr);
  ierr = VecSet(fields.U,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(ctx.daScal,fields.V,&v_array);CHKERRQ(ierr);

  ierr = VecSet(fields.BCU,0.0);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx.daVect,fields.BCU,&bcu_array);CHKERRQ(ierr);

  /*
   Reset all BC for U and V
   */
  for (i = 0; i < 6; i++) {
    ctx.bcV[0].face[i] = NONE;
    for (j = 0; j < 3; j++) {
      ctx.bcU[j].face[i] = NONE;
    }
  }
  for (i = 0; i < 12; i++) {
    ctx.bcV[0].edge[i] = NONE;
    for (j = 0; j < 3; j++) {
      ctx.bcU[j].edge[i] = NONE;
    }
  }
  for (i = 0; i < 8; i++) {
    ctx.bcV[0].vertex[i] = NONE;
    for (j = 0; j < 3; j++) {
      ctx.bcU[j].vertex[i] = NONE;
    }
  }
  switch (orientation) {
  case 0:
    ierr                = PetscPrintf(PETSC_COMM_WORLD,"Applying traction Dirichlet conditions on faces X0 X1\n");CHKERRQ(ierr);
    ctx.bcU[0].face[X0] = FIXED;ctx.bcU[0].face[X1] = FIXED;
    ctx.bcU[1].face[X0] = ZERO;ctx.bcU[1].face[X1] = ZERO;
    ctx.bcU[2].face[X0] = ZERO;ctx.bcU[2].face[X1] = ZERO;

    ctx.bcU[2].face[Z0] = ZERO;ctx.bcU[2].face[Z1] = ZERO;
    ctx.bcU[1].face[Z0] = ZERO;ctx.bcU[1].face[Z1] = ZERO;
      
    ctx.bcV[0].face[X0] = ONE;
		ctx.bcV[0].face[X1] = ONE;
    for (k = zs; k < zs+zm; k++) {
      for (j = ys; j < ys+ym; j++) {
        for (i = xs; i < xs+xm; i++) {
          if (i == 0) {
            bcu_array[k][j][i][0] = -bc;
          }
          if (i == nx-1) {
            bcu_array[k][j][i][0] = bc;
          }
        }
      }
    }
    break;
	  case 1:
		  ierr                = PetscPrintf(PETSC_COMM_WORLD,"Applying traction Dirichlet conditions on faces X0 X1\n");CHKERRQ(ierr);
		  ctx.bcU[0].face[X0] = ZERO;ctx.bcU[0].face[X1] = ZERO;
		  ctx.bcU[1].face[X0] = FIXED;ctx.bcU[1].face[X1] = FIXED;
		  ctx.bcU[2].face[X0] = ZERO;ctx.bcU[2].face[X1] = ZERO;
		  
		  for (k = zs; k < zs+zm; k++) {
			  for (j = ys; j < ys+ym; j++) {
				  for (i = xs; i < xs+xm; i++) {
					  if (((i == nx/2)) || (i == nx/2-1)) {
						  v_array[k][j][i] = 0.;
					  }
					  if (i == 0) {
						  bcu_array[k][j][i][1] = -bc;
					  }
					  if (i == nx-1) {
						  bcu_array[k][j][i][1] = bc;
					  }
				  }
			  }
		  }
		  break;
	  case 2:
		  ierr                = PetscPrintf(PETSC_COMM_WORLD,"Applying traction Dirichlet conditions on faces X0 X1\n");CHKERRQ(ierr);
		  ctx.bcU[0].face[X0] = ZERO;ctx.bcU[0].face[X1] = ZERO;
		  ctx.bcU[1].face[X0] = ZERO;ctx.bcU[1].face[X1] = ZERO;
		  ctx.bcU[2].face[X0] = FIXED;ctx.bcU[2].face[X1] = FIXED;
		  
		  for (k = zs; k < zs+zm; k++) {
			  for (j = ys; j < ys+ym; j++) {
				  for (i = xs; i < xs+xm; i++) {
					  if (((i == nx/2)) || (i == nx/2-1)) {
						  v_array[k][j][i] = 0.;
					  }
					  if (i == 0) {
						  bcu_array[k][j][i][2] = -bc;
					  }
					  if (i == nx-1) {
						  bcu_array[k][j][i][2] = bc;
					  }
				  }
			  }
		  }
		  break;
	  case 3:
		  ierr                = PetscPrintf(PETSC_COMM_WORLD,"Applying traction Dirichlet conditions on faces X0 X1\n");CHKERRQ(ierr);
		  ctx.bcU[0].face[X0] = FIXED;ctx.bcU[0].face[X1] = FIXED;
		  ctx.bcU[1].face[X0] = FIXED;ctx.bcU[1].face[X1] = FIXED;
		  ctx.bcU[2].face[X0] = ZERO;ctx.bcU[2].face[X1] = ZERO;
		  
		  ctx.bcU[2].face[Z0] = ZERO;ctx.bcU[2].face[Z1] = ZERO;

		  for (k = zs; k < zs+zm; k++) {
			  for (j = ys; j < ys+ym; j++) {
				  for (i = xs; i < xs+xm; i++) {
					  if (i == 0) {
						  bcu_array[k][j][i][0] = -bc;
						  bcu_array[k][j][i][1] = -bc;
					  }
					  if (i == nx-1) {
						  bcu_array[k][j][i][0] = bc;
						  bcu_array[k][j][i][1] = bc;
					  }
				  }
			  }
		  }
		  break;
  default:
    SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Orientation should be between 0 and 2, got %i\n",orientation);
    break;
  }
  ierr = DMDAVecRestoreArray(ctx.daScal,fields.V,&v_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,fields.BCU,&bcu_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ctx.hasCrackPressure = PETSC_FALSE;

  ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
	do {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"alt min step %i\n",altminit);CHKERRQ(ierr);
		ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
		ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
		ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);
		ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
		ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
		ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"   Max. change on V: %e\n",errV);CHKERRQ(ierr);
		altminit++;
	} while (errV > ctx.altmintol && altminit <= ctx.altminmaxit);
  ctx.ElasticEnergy = 0;
  ctx.InsituWork    = 0;
  ctx.PressureWork  = 0.;
  ierr              = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,&fields,&ctx);CHKERRQ(ierr);
  ierr              = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);

  ctx.TotalEnergy = ctx.ElasticEnergy-ctx.InsituWork-ctx.PressureWork;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface energy:            %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic Energy:            %e\n",ctx.ElasticEnergy);CHKERRQ(ierr);
  if (ctx.hasCrackPressure) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of pressure forces:   %e\n",ctx.PressureWork);CHKERRQ(ierr);
  }
  if (ctx.hasInsitu) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of surface forces:    %e\n",ctx.InsituWork);CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Total Mechanical energy:  %e\n",ctx.ElasticEnergy-InsituWork-ctx.PressureWork);CHKERRQ(ierr);

	ierr = VolumetricCrackOpening(&ctx.CrackVolume, &ctx, &fields);CHKERRQ(ierr);
  PetscReal   functionvalue;
	ierr = TrialFunctionCompute(&functionvalue,&ctx,&fields);CHKERRQ(ierr);

    ierr = FieldsH5Write(&ctx,&fields);
    ierr = FieldsH5Write(&ctx,&fields);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"###################################################################\n\n\n");CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#        Actual crack volume change = %f\t      \n\n\n\n", (lz*ly*bc*2));CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#        VF crack volume change = %e\t      \n\n", ctx.CrackVolume);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#        FunctionValue = %e\t      \n", functionvalue);CHKERRQ(ierr);

  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}

