/*
 test33.c: 3D Flow problem with source term [pressure = 1/(3*pi^2)*cos(pi*x)*cos(pi*y)*cos(pi*z)]. All velocity boundary condition
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFV.h"
#include "VFU.h"
#include "VFFlow.h"

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
	PetscInt		i,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm;
	PetscReal		BBmin[3],BBmax[3];
	PetscReal		****flowbc_array;
	PetscReal		***src_array;
	PetscReal		****coords_array;
	PetscReal		hx,hy,hz;
	PetscReal		gx,gy,gz;
	PetscReal		gamma, beta, rho, mu;
	PetscReal		pi;
		
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = FlowSolverInitialize(&ctx,&fields);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	ierr = VecSet(fields.FlowBCArray,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.Source,0.);CHKERRQ(ierr);
	
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);	
	ierr = DMDAVecGetArrayDOF(ctx.daFlow,fields.FlowBCArray,&flowbc_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx.daScal,ctx.Source,&src_array);CHKERRQ(ierr);

	pi = 6.*asin(0.5);
	hx = 1./(nx-1);
	hy = 1./(nx-1);
	hz = 1./(nz-1);	
	rho = ctx.flowprop.rho;									 
	mu = ctx.flowprop.mu;     
	beta = ctx.flowprop.beta;		
	gamma = ctx.flowprop.gamma;									
    gx = ctx.flowprop.g[0];
    gy = ctx.flowprop.g[1];
    gz = ctx.flowprop.g[2];
	/*
	 Reset all Flow BC for velocity and P
	 */
	for (i = 0; i < 6; i++) {
		for (c = 0; c < 4; c++) {
			ctx.bcFlow[c].face[i] = NOBC;
		}
	}
	for (i = 0; i < 12; i++) {
		for (c = 0; c < 4; c++) {
			ctx.bcFlow[c].edge[i] = NOBC;
		}
	}
	for (i = 0; i < 8; i++) {
		for (c = 0; c < 4; c++) {
			ctx.bcFlow[c].vertex[i] = NOBC;
		}
	}
	ctx.bcFlow[0].face[X0] = VELOCITY;
	ctx.bcFlow[0].face[X1] = VELOCITY;
	ctx.bcFlow[1].face[Y0] = VELOCITY;
	ctx.bcFlow[1].face[Y1] = VELOCITY;
	ctx.bcFlow[2].face[Z0] = VELOCITY;
	ctx.bcFlow[2].face[Z1] = VELOCITY;		
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				flowbc_array[k][j][i][0] = beta/mu*(sin(pi*i*hx)*cos(pi*j*hy)*cos(pi*k*hz)-gamma*rho*gx)/(3.*pi);
				flowbc_array[k][j][i][1] = beta/mu*(cos(pi*i*hx)*sin(pi*j*hy)*cos(pi*k*hz)-gamma*rho*gy)/(3.*pi);
				flowbc_array[k][j][i][2] = beta/mu*(cos(pi*i*hx)*cos(pi*j*hy)*sin(pi*k*hz)-gamma*rho*gz)/(3.*pi);
			}
		}
	}	
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				src_array[k][j][i] = beta/mu*cos(pi*k*hz)*cos(pi*j*hy)*cos(pi*i*hx); 
			}
		}
	}
	ierr = DMDAVecRestoreArray(ctx.daScal,ctx.Source,&src_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daFlow,fields.FlowBCArray,&flowbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	/* 
	 Now done with all initializations
	 */
	ctx.maxtimestep = 1;
	for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nProcessing step %i.\n",ctx.timestep);CHKERRQ(ierr);
		ctx.timevalue = ctx.timestep * ctx.maxtimevalue / (ctx.maxtimestep-1.);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\ntime value %f \n",ctx.timevalue);CHKERRQ(ierr);
		/*
		 Do flow solver step 
		 */
		ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
		/*
		 Save fields and write statistics about current run
		 */    
		switch (ctx.fileformat) {
			case FILEFORMAT_HDF5:       
				ierr = FieldsH5Write(&ctx,&fields);
				break;
			case FILEFORMAT_BIN:
				ierr = FieldsBinaryWrite(&ctx,&fields);
				break; 
		} 
		ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.log",ctx.prefix);CHKERRQ(ierr);
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&logviewer);CHKERRQ(ierr);
		ierr = PetscLogView(logviewer);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&logviewer);
	}
	ierr = FlowSolverFinalize(&ctx,&fields);CHKERRQ(ierr);
//	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

