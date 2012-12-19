/*
 test30.c: 3D TS. Flow problem with source term [pressure = sin(2*pi*x)*sin(2*pi*y)*sin(2(pi*z)]. All velocity boundary condition
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
	/*	Set flow solver type	*/
	ctx.flowsolver = FLOWSOLVER_TSMIXEDFEM;
	ierr = PetscOptionsEnum("-flowsolver","\n\tFlow solver","",VFFlowSolverName,(PetscEnum)ctx.flowsolver,(PetscEnum*)&ctx.flowsolver,PETSC_NULL);CHKERRQ(ierr);
	
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
	ctx.bcQ[0].face[X0] = VALUE;
	ctx.bcQ[0].face[X1] = VALUE;
	ctx.bcQ[1].face[Y0] = VALUE;
	ctx.bcQ[1].face[Y1] = VALUE;
	ctx.bcQ[2].face[Z0] = VALUE;
	ctx.bcQ[2].face[Z1] = VALUE;		
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				flowbc_array[k][j][i][0] = -beta/mu*(2.*pi*cos(2.*pi*i*hx)*sin(2.*pi*j*hy)*sin(2.*pi*k*hz)-gamma*rho*gx);
				flowbc_array[k][j][i][1] = -beta/mu*(2.*pi*sin(2.*pi*i*hx)*cos(2.*pi*j*hy)*sin(2.*pi*k*hz)-gamma*rho*gy);
				flowbc_array[k][j][i][2] = -beta/mu*(2.*pi*sin(2.*pi*i*hx)*sin(2.*pi*j*hy)*cos(2.*pi*k*hz)-gamma*rho*gz);
				flowbc_array[k][j][i][3] = (sin(2.*pi*i*hx)*sin(2.*pi*j*hy)*sin(2.*pi*k*hz));
			}
		}
	}	
	
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				src_array[k][j][i] = 4.*pi*pi*3.*beta/mu*sin(2.*pi*k*hz)*sin(2.*pi*j*hy)*sin(2.*pi*i*hx);
			}
		}
	}
	ierr = DMDAVecRestoreArray(ctx.daScal,ctx.Source,&src_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daFlow,fields.FlowBCArray,&flowbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	
	PetscViewer     viewer;
	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"BCValues.txt",&viewer);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
	ierr = VecView(fields.FlowBCArray,viewer);CHKERRQ(ierr);

	
	/* Setting time parameters	*/
	ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
	ctx.maxtimestep = 2;
	ctx.maxtimevalue = 60.;
	ctx.timevalue = 1;
	/*	Do flow solver step	*/
	ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
	/*	Save fields and write statistics about current run	*/    
	ierr = FieldsH5Write(&ctx,&fields);
	ierr = FlowSolverFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

