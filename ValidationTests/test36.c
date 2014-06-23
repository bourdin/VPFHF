/*
 test36.c: 1D. Flow problem with source term = 1 and Homogeneous pressure boundary conditions on all sides. Analytical solution is p = x(x-1)/2
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
./test36 -n 101,2,101 -l 1,0.01,1 -m_inv 0 -E 17 -0 -nu 0.2 -npc 1 -pc0_r 0.1 -pc0_center 0.3,0.005,0.5 -epsilon 0.03eta 0 -pc0_phi 40 -pc0_thickness 0.01 -atnum 1 -flowsolver FLOWSOLVER_snesmixeDFEM -perm 1e-14 -miu 1e-13 -num 2 -timestepsize 0.1  -Gc 5e6
 
 
 ./test36 -n 51,2,51 -l 1,0.01,1 -m_inv 0 -E 17 -0 -nu 0.2 -npc 1 -pc0_r 0.1 -pc0_center 0.5,0.005,0.5 -epsilon 0.02 -pc0_theta 0 -pc0_phi 30 -pc0_thickness 0.04 -atnum 2 -flowsolver FLOWSOLVER_snesmixeDFEM -perm 1e-14 -miu 1e-13 -num 2 -timestepsize 0.1  -Gc 5e6 -permmax 1e-12
 
 ./test36 -n 101,2,101 -l 1,0.01,1 -m_inv 0 -E 17 -0 -nu 0.2 -npc 1 -pc0_r 0.1 -pc0_center 0.5,0.005,0.5 -epsilon 0.04 -pc0_theta 0 -pc0_phi 0 -pc0_thickness 0.02 -atnum 2 -flowsolver FLOWSOLVER_snesmixeDFEM -perm 1e-14 -miu 1e-13 -num 2 -timestepsize 0.1  -Gc 5e6 -permmax 1e-12
 
 ./test36 -n 101,2,101 -l 1,0.01,1 -m_inv 0 -E 17 -0 -nu 0.2 -npc 2 -pc0_r 0.1 -pc0_center 0.45,0.005,0.5 -epsilon 0.04 -pc0_theta 0 -pc0_phi 90 -pc0_thickness 0.02 -pc1_r 0.1 -pc1_center 0.55,0.005,0.5 -pc1_theta 0 -pc1_phi 90 -pc1_thickness 0.02 -atnum 2 -flowsolver FLOWSOLVER_snesmixeDFEM -perm 1e-14 -miu 1e-13 -num 2 -timestepsize 0.1  -Gc 5e6 -permmax 1e-12
 
 
 ./test36 -n 101,2,101 -l 1,0.01,1  -E 17 -0 -nu 0.2 -npc 1 -pc0_r 0.2 -pc0_center 0.5,0.005,0.5 -epsilon 0.01 -pc0_theta 0 -pc0_phi 45 -pc0_thickness 0.02 -atnum 2 -flowsolver FLOWSOLVER_snesmixeDFEM  -Gc 5e6
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFMech.h"
#include "VFFlow.h"
#include "VFPermfield.h"


VFCtx               ctx;
VFFields            fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{	
	PetscErrorCode  ierr;
	PetscViewer		viewer;
	PetscViewer     logviewer;
	char			filename[FILENAME_MAX];
	PetscInt		i,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm,ite;
	PetscReal		BBmin[3],BBmax[3];
	PetscReal		****coords_array;
	PetscReal		hx,hy,hz;
	PetscReal		lx,ly,lz;
  PetscInt    xs1,xm1,ys1,ym1,zs1,zm1;
  PetscReal      ****bcu_array;


	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx.daVFperm,&xs1,&ys1,&zs1,&xm1,&ym1,&zm1);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	ierr = VecSet(fields.U,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.U_old,0.);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  ierr = VecSet(fields.BCU,0.0);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx.daVect,fields.BCU,&bcu_array);CHKERRQ(ierr);
	lz = BBmax[2]-BBmin[2];
	ly = BBmax[1]-BBmin[1];
	lx = BBmax[0]-BBmin[0];
	hx = lx/(nx-1);
	hy = ly/(nx-1);
	hz = lz/(nz-1);
  
  for (j = ys; j < ys+ym; j++) {
    for (i = xs; i < xs+xm; i++) {
      bcu_array[0][j][i][2] = -0.1;
      bcu_array[0][j][i][0] = -0.;
      
      bcu_array[nz-1][j][i][2] = 0.5;
      bcu_array[nz-1][j][i][0] = 0.;
    }
  }
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,fields.BCU,&bcu_array);CHKERRQ(ierr);
/*      Mechanical model settings       */
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
  ctx.bcU[0].face[X0]= ZERO;
  ctx.bcU[1].face[X0]= ZERO;
  
  ctx.bcU[0].face[X1]= ZERO;
  ctx.bcU[1].face[X1]= ZERO;
  
  
  ctx.bcU[0].face[Z0]= ZERO;
  ctx.bcU[1].face[Z0]= ZERO;
  ctx.bcU[2].face[Z0]= FIXED;
  
  ctx.bcU[0].face[Z1]= ZERO;
  ctx.bcU[1].face[Z1]= ZERO;
  ctx.bcU[2].face[Z1]= FIXED;

  ctx.bcU[1].face[Y0]= ZERO;
  ctx.bcU[1].face[Y1]= ZERO;
  
  ctx.bcV[0].face[Z0] = ONE;
  ctx.bcV[0].face[Z1] = ONE;
  
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(ctx.U_old,0.);CHKERRQ(ierr);

  
  
  
  ctx.timestep = 0;
  ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
  ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
  ierr = FieldsH5Write(&ctx,&fields);

  
  
  
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

