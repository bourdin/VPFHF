/*
 test37.c: 2D SNES. Single well problem. Source term implemented as dirac function and analytical solution is 2D Green's function. All pressure boundary condition.
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
 ./test37  -n 201,201,2 -l 1,1,1 -maxtimestep 1 timestepsize 1 -theta 1 -nc 2 -w0_coords 0.50001,0.50001,0.00 -w0_Qw 0.5 -w0_constraint Rate -w0_rw 0.1 -w0_type INJECTOR -w1_coords 0.50001,0.50001,1 -w1_Qw 0.5 -w1_constraint Rate -w1_rw 0.1 -w1_type INJECTOR -m_inv 0
 
 ./test37  -n 101,101,2 -l 1,1,1 -maxtimestep 1 timestepsize 1 -theta 1 -nc 2 -w0_coords 0.50001,0.50001,0.00 -w0_Qw 0.5 -w0_constraint Rate -w0_rw 0.1 -w0_type INJECTOR -w1_coords 0.50001,0.50001,1 -w1_Qw 0.5 -w1_constraint Rate -w1_rw 0.1 -w1_type INJECTOR -m_inv 0
 
 
 ./test37  -n 101,101,2 -l 1,1,1 -maxtimestep 1 timestepsize 1 -theta 1 -nc 1 -w0_coords 0.50001,0.50001,0.5 -w0_Qw 1 -w0_constraint Rate -w0_rw 0.1 -w0_type producer -m_inv 0
 
 ./test37  -n 101,101,2 -l 1,1,1 -maxtimestep 1 timestepsize 10 -theta 1 -nc 1 -w0_coords 0.5,0.5,0.5 -w0_Qw 1 -w0_constraint Rate -w0_rw 0.1 -w0_type producer -m_inv 0
 
 ./test37  -n 101,101,2 -l 1000,1000,10 -maxtimestep 3 -timestepsize 8640 -theta 1 -nc 1 -w0_coords 500,500,5 -w0_Qw 10 -w0_constraint Rate -w0_rw 0.1 -w0_type injector -m_inv 1
 
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFMech.h"
#include "VFFlow.h"
#include "VFWell.h"
#include "VFCracks.h"


VFCtx               ctx;
VFFields            fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode  ierr;
	PetscInt		i,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm,xs1,xm1,ys1,ym1,zs1,zm1;
	PetscReal		BBmin[3],BBmax[3];
	PetscReal		***presbc_array;
	PetscReal		****velbc_array;
	PetscReal		****coords_array;
	PetscReal		gx,gy,gz;
	PetscReal		gamma, beta, rho, mu;
	PetscReal		****perm_array;
	PetscReal		****velnpre_array;
	PetscReal		***pre_array;
	char			prefix[PETSC_MAX_PATH_LEN+1];
	PetscReal		pi,dist;
	
	
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ctx.flowsolver = FLOWSOLVER_SNESMIXEDFEM;
  ctx.fractureflowsolver = FRACTUREFLOWSOLVER_NONE;
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	ierr = VecSet(ctx.PresBCArray,10.);CHKERRQ(ierr);
	ierr = VecSet(ctx.VelBCArray,0.);CHKERRQ(ierr);
	ierr = VecSet(fields.VelnPress,10.);CHKERRQ(ierr);
	ierr = VecSet(fields.vfperm,0.);CHKERRQ(ierr);
	ctx.hasFlowWells = PETSC_TRUE;
	ctx.hasFluidSources = PETSC_FALSE;
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daFlow,fields.VelnPress,&velnpre_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daVFperm,&xs1,&ys1,&zs1,&xm1,&ym1,&zm1);CHKERRQ(ierr);
	ctx.numWells = 1;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-nc",&ctx.numWells,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscMalloc(ctx.numWells*sizeof(VFWell),&ctx.well);CHKERRQ(ierr);
	if(ctx.hasFlowWells == PETSC_TRUE){
		for (i = 0; i < ctx.numWells; i++) {
			ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN,"w%d_",i);CHKERRQ(ierr);
			ierr = VFWellCreate(&ctx.well[i]);CHKERRQ(ierr);
			ierr = VFWellGet(prefix,&ctx.well[i]);CHKERRQ(ierr);
			ierr = VFWellView(&ctx.well[i],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
		}
		for (k = zs1; k < zs1+zm1; k++) {
			for (j = ys1; j < ys1+ym1; j++) {
				for (i = xs1; i < xs1+xm1; i++) {
					perm_array[k][j][i][0] = 1;
					perm_array[k][j][i][1] = 1;
					perm_array[k][j][i][2] = 1;
          
          if( (j == ny/2-1 || j == ny/2) && (i > nx/8 && i < 6*nx/8)){
              //          if( (j == ny/2-1 || j == ny/2) /*&& (i > nx/8 && i < 6*nx/8)*/){
              //            if( (i == nx/2-1 || i == nx/2) /*&& (i > nx/8 && i < 6*nx/8)*/){
            perm_array[k][j][i][0] = 10;
              //            perm_array[k][j][i][1] = 10;
              //            perm_array[k][j][i][2] = 10;
          }
					perm_array[k][j][i][3] = 0;
					perm_array[k][j][i][4] = 0;
					perm_array[k][j][i][5] = 0;
          
				}
			}
		}
	}
	pi = 6.*asin(0.5);
	ctx.flowprop.rho = 1000.;
	ctx.flowprop.mu = 1;
	beta = ctx.flowprop.beta;
	gamma = ctx.flowprop.gamma;
  gx = ctx.flowprop.g[0];
  gy = ctx.flowprop.g[1];
  gz = ctx.flowprop.g[2];
	/*
	 Reset all Flow BC for velocity and P
	 */
	for (i = 0; i < 6; i++) {
		ctx.bcP[0].face[i] = NONE;
		for (c = 0; c < 3; c++) {
			ctx.bcQ[c].face[i] = NONE;
		}
	}
	for (i = 0; i < 12; i++) {
		ctx.bcP[0].edge[i] = NONE;
		for (c = 0; c < 3; c++) {
			ctx.bcQ[c].edge[i] = NONE;
		}
	}
	for (i = 0; i < 8; i++) {
		ctx.bcP[0].vertex[i] = NONE;
		for (c = 0; c < 3; c++) {
			ctx.bcQ[c].vertex[i] = NONE;
		}
	}
  /*
   ctx.bcQ[0].face[X0] = VALUE;
   ctx.bcQ[0].face[X1] = VALUE;
   ctx.bcQ[1].face[Y0] = VALUE;
   ctx.bcQ[1].face[Y1] = VALUE;
   ctx.bcQ[2].face[Z0] = VALUE;
   ctx.bcQ[2].face[Z1] = VALUE;
   */
  
  ctx.bcP[0].face[X0] = VALUE;
  ctx.bcP[0].face[X1] = VALUE;
  ctx.bcP[0].face[Y0] = VALUE;
  ctx.bcP[0].face[Y1] = VALUE;
  ctx.bcQ[2].face[Z0] = VALUE;
  ctx.bcQ[2].face[Z1] = VALUE;
  for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				velbc_array[k][j][i][0] = 0.;
				velbc_array[k][j][i][1] = 0.;
				velbc_array[k][j][i][2] = 0.;
			}
		}
	}
	for(k = zs; k < zs+zm; k++){
		for(j = ys; j < ys+ym; j++){
			presbc_array[k][j][0] = 10;
      presbc_array[k][j][nx-1] = 10;
			
		}
	}
	for(k = zs; k < zs+zm; k++){
		for(i = xs; i < xs+xm; i++){
			presbc_array[k][0][i] = 10;
      presbc_array[k][ny-1][i] = 10;
			
		}
	}
	ierr = DMDAVecRestoreArrayDOF(ctx.daFlow,fields.VelnPress,&velnpre_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-theta",&ctx.flowprop.theta,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-timestepsize",&ctx.flowprop.timestepsize,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
	for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
		ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
		ierr = FieldsH5Write(&ctx,&fields);
	}
	ierr = DMDAVecGetArray(ctx.daScal,fields.pressure,&pre_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx.daScal,fields.pressure,&pre_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	++ctx.timestep;
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

