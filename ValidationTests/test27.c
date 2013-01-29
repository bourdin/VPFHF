/*
 test27.c: 3D KSP. Flow problem with source term [pressure = 1/(3*pi^2)*cos(pi*x)*cos(pi*y)*cos(pi*z)]. All pressure boundary condition
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
 ./test27 -n 11,11,11 -l 1,1,1 -m_inv 0 -theta 1 -velp_pc_type jacobi -maxtimestep 1
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
	PetscReal		***presbc_array;
	PetscReal		***src_array;
	PetscReal		****coords_array;
	PetscReal		hx,hy,hz;
	PetscReal		gx,gy,gz;
	PetscReal		lx,ly,lz;
	PetscReal		gamma, beta, rho, mu;
	PetscReal		pi;
		
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ctx.flowsolver = FLOWSOLVER_KSPMIXEDFEM;
	ierr = FlowSolverInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	ierr = VecSet(ctx.PresBCArray,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.Source,0.);CHKERRQ(ierr);
	ctx.hasFluidSources = PETSC_TRUE;
	ctx.hasFlowWells = PETSC_FALSE;
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);	
	ierr = DMDAVecGetArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx.daScal,ctx.Source,&src_array);CHKERRQ(ierr);
	pi = 6.*asin(0.5);
	rho = ctx.flowprop.rho;									 
	mu = ctx.flowprop.mu;     
	beta = ctx.flowprop.beta;		
	gamma = ctx.flowprop.gamma;									
    gx = ctx.flowprop.g[0];
    gy = ctx.flowprop.g[1];
    gz = ctx.flowprop.g[2];
	lz = BBmax[2]-BBmin[2];
	ly = BBmax[1]-BBmin[1];
	lx = BBmax[0]-BBmin[0];
	hx = lx/(nx-1);
	hy = ly/(nx-1);
	hz = lz/(nz-1);	
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
	for (i = 0; i < 6; i++) {
		ctx.bcP[0].face[i] = VALUE;
	}	
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				presbc_array[k][j][i] = (cos(pi*i*hx/lx)*cos(pi*j*hy/ly)*cos(pi*k*hz/lz)-gamma*rho*gz)/(3.*pi*pi);
			}
		}
	}	
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				src_array[k][j][i] = beta/mu*cos(pi*k*hz/lz)*cos(pi*j*hy/ly)*cos(pi*i*hx/lx)/3.*((1./(lx*lx))+(1./(ly*ly))+(1./(lz*lz))); 
			}
		}
	}
	ierr = DMDAVecRestoreArray(ctx.daScal,ctx.Source,&src_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-theta",&ctx.flowprop.theta,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-timestepsize",&ctx.flowprop.timestepsize,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
	ctx.maxtimestep = 1;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-maxtimestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr);
	for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nProcessing step %i.\n",ctx.timestep);CHKERRQ(ierr);
		ctx.timevalue = ctx.timestep * ctx.maxtimevalue / (ctx.maxtimestep-1.);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\ntime value %f \n",ctx.timevalue);CHKERRQ(ierr);
		ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
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
	Vec error;
	PetscReal norm_1,norm_2,norm_inf;
	ierr = VecDuplicate(fields.VelnPress,&error);
	ierr = VecWAXPY(error,-1.0,fields.VelnPress,fields.FlowBCArray);
	ierr = VecNorm(error,NORM_1,&norm_1);
	ierr = VecNorm(error,NORM_2,&norm_2);
	ierr = VecNorm(error,NORM_INFINITY,&norm_inf);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n1_NORM = %f \n 2_norm = %f \n inf_norm = %f \n",norm_1, norm_2,norm_inf);CHKERRQ(ierr);	
	ierr = FlowSolverFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

