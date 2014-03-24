/*
 test41.c: 3D SNES. Single well problem. Source term implemented as dirac function and analytical solution is 3D Green's function. All pressure boundary condition.
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
 ./test41  -n 41,41,41 -l 1,1,1 -maxtimestep 1 timestepsize 1 -theta 1 -nc 1 -w0_coords 0.50001,0.50001,0.50001 -w0_Qw 1 -w0_constraint Rate -w0_rw 0.1 -w0_type INJECTOR -m_inv 0.
 */

#include "petsc.h"
#include "VFCartFE.h"
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
	PetscViewer		viewer;
	PetscViewer     logviewer;
	char			filename[FILENAME_MAX];
	PetscInt		i,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm,xs1,xm1,ys1,ym1,zs1,zm1;
	PetscReal		BBmin[3],BBmax[3];
	PetscReal		***presbc_array;
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
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	ierr = VecSet(ctx.PresBCArray,0.);CHKERRQ(ierr);
	ierr = VecSet(fields.VelnPress,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.Source,0.);CHKERRQ(ierr);
	ctx.hasFlowWells = PETSC_TRUE;
	ctx.hasFluidSources = PETSC_FALSE;
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);	
	ierr = DMDAVecGetArrayDOF(ctx.daFlow,fields.VelnPress,&velnpre_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr); 
	ierr = DMDAGetCorners(ctx.daVFperm,&xs1,&ys1,&zs1,&xm1,&ym1,&zm1);CHKERRQ(ierr);
	ctx.numWells = 0;	
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
					for(c = 3; c < 6; c++){
						perm_array[k][j][i][c] = 0.;
					}
				}
			}
		}
	}
	pi = 6.*asin(0.5);
	rho = ctx.flowprop.rho;									 
	mu = ctx.flowprop.mu;     
	beta = ctx.flowprop.beta;		
	gamma = ctx.flowprop.gamma;									
    gx = ctx.flowprop.g[0];
    gy = ctx.flowprop.g[1];
    gz = ctx.flowprop.g[2];
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
	ctx.bcP[0].face[X0] = FIXED;
	ctx.bcP[0].face[X1] = FIXED;
	ctx.bcP[0].face[Y0] = FIXED;
	ctx.bcP[0].face[Y1] = FIXED;
	ctx.bcP[0].face[Z0] = FIXED;
	ctx.bcP[0].face[Z1] = FIXED;
	for(k = zs; k < zs+zm; k++){
		for(j = ys; j < ys+ym; j++){
			dist = sqrt(pow((ctx.well[0].coords[0]-coords_array[k][j][0][0]),2)+pow((ctx.well[0].coords[1]-coords_array[k][j][0][1]),2)+pow((ctx.well[0].coords[2]-coords_array[k][j][0][2]),2));
			presbc_array[k][j][0] = -1./(4.*pi*dist);
			
			dist = sqrt(pow((ctx.well[0].coords[0]-coords_array[k][j][nx-1][0]),2)+pow((ctx.well[0].coords[1]-coords_array[k][j][nx-1][1]),2)+pow((ctx.well[0].coords[2]-coords_array[k][j][nx-1][2]),2));
			presbc_array[k][j][nx-1] = -1./(4.*pi*dist);
			
		}
	}
	for(k = zs; k < zs+zm; k++){
		for(i = xs; i < xs+xm; i++){
			dist = sqrt(pow((ctx.well[0].coords[0]-coords_array[k][0][i][0]),2)+pow((ctx.well[0].coords[1]-coords_array[k][0][i][1]),2)+pow((ctx.well[0].coords[2]-coords_array[k][0][i][2]),2));
			presbc_array[k][0][i] = -1./(4.*pi*dist);
			
			dist = sqrt(pow((ctx.well[0].coords[0]-coords_array[k][ny-1][i][0]),2)+pow((ctx.well[0].coords[1]-coords_array[k][ny-1][i][1]),2)+pow((ctx.well[0].coords[2]-coords_array[k][ny-1][i][2]),2));
			presbc_array[k][ny-1][i] = -1./(4.*pi*dist);
			
		}
	}
	for(j = ys; j < ys+ym; j++){
		for(i = xs; i < xs+xm; i++){
			dist = sqrt(pow((ctx.well[0].coords[0]-coords_array[0][j][i][0]),2)+pow((ctx.well[0].coords[1]-coords_array[0][j][i][1]),2)+pow((ctx.well[0].coords[2]-coords_array[0][j][i][2]),2));
			presbc_array[0][j][i] = -1./(4.*pi*dist);
			
			dist = sqrt(pow((ctx.well[0].coords[0]-coords_array[nz-1][j][i][0]),2)+pow((ctx.well[0].coords[1]-coords_array[nz-1][j][i][1]),2)+pow((ctx.well[0].coords[2]-coords_array[nz-1][j][i][2]),2));
			presbc_array[nz-1][j][i] = -1./(4.*pi*dist);
			
		}
	}
	ierr = PetscOptionsGetReal(PETSC_NULL,"-theta",&ctx.flowprop.theta,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-timestepsize",&ctx.flowprop.timestepsize,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
	for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
		ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
		ierr = FieldsH5Write(&ctx,&fields);
	}
	ierr = VecSet(fields.VelnPress,0.);CHKERRQ(ierr);
	ierr = VecSet(fields.velocity,0.);CHKERRQ(ierr);
	ierr = VecSet(fields.pressure,0.);CHKERRQ(ierr);
	ierr = VecSet(fields.U,0.);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx.daScal,fields.pressure,&pre_array);CHKERRQ(ierr);
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				dist = sqrt(pow((ctx.well[0].coords[0]-coords_array[k][j][i][0]),2)+pow((ctx.well[0].coords[1]-coords_array[k][j][i][1]),2)+pow((ctx.well[0].coords[2]-coords_array[k][j][i][2]),2));
				pre_array[k][j][i] = -1./(4.*pi*dist);
			}
		}
	}
	ierr = DMDAVecRestoreArray(ctx.daScal,fields.pressure,&pre_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daFlow,fields.VelnPress,&velnpre_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);	
	++ctx.timestep;
	ierr = FieldsH5Write(&ctx,&fields);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

