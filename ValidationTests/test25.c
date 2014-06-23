/*
 test25.c: 3D KSP. Flow problem with source term [pressure = sin(2*pi*x)*sin(2*pi*y)*sin(2(pi*z)]. All pressure boundary condition.
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
 mpiexec -np 2 ./test25 -n 11,11,11 -l 1,1,1 -flowsolver FLOWSOLVER_SNESMIXEDFEM -M_inv 0
 ./test25 -n 11,11,11 -l 1,1,1 -m_inv 0 -ts_type beuler -ts_dt 1 -ts_max_steps 2 -flowsolver FLOWSOLVER_TSMIXEDFEM

 */

#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFMech.h"
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
  PetscReal   vol,vol1,vol2,vol3,vol4,vol5;
	PetscReal   ****perm_array;
	PetscInt    xs1,xm1,ys1,ym1,zs1,zm1;
		
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ctx.flowsolver = FLOWSOLVER_KSPMIXEDFEM;
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
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
	
	ierr = DMDAGetCorners(ctx.daVFperm,&xs1,&ys1,&zs1,&xm1,&ym1,&zm1);CHKERRQ(ierr);	
	ierr = DMDAVecGetArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr); 
	for (k = zs1; k < zs1+zm1; k++) {
		for (j = ys1; j < ys1+ym1; j++) {
				for (i = xs1; i < xs1+xm1; i++) {
				perm_array[k][j][i][3] = 0.;
				perm_array[k][j][i][4] = 0.;			
				perm_array[k][j][i][5] = 0.;				
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
		ctx.bcP[0].face[i] = FIXED;
	}	
	
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				presbc_array[k][j][i] = sin(2.*pi*i*hx/lx)*sin(2.*pi*j*hy/ly)*sin(2.*pi*k*hz/lz);
			}
		}
	} 
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				src_array[k][j][i] = 4.*pi*pi*beta/mu*sin(2.*pi*k*hz/lz)*sin(2.*pi*j*hy/ly)*sin(2.*pi*i*hx/lx)*( (1/(lx*lx))+(1/(ly*ly))+(1/(lz*lz)) ); 
			}
		}
	} 
	ierr = DMDAVecRestoreArray(ctx.daScal,ctx.Source,&src_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);		
	ierr = PetscOptionsGetReal(PETSC_NULL,"-theta",&ctx.flowprop.theta,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-timestepsize",&ctx.flowprop.timestepsize,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
	ctx.maxtimestep = 1;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-maxtimestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr);
  if(ctx.flowsolver == FLOWSOLVER_TSMIXEDFEM || ctx.flowsolver == FLOWSOLVER_TS){
  	ctx.maxtimestep = 1;
    ctx.maxtimevalue = 60.;
    ctx.timevalue = 1.;
    ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
    ierr = FieldsH5Write(&ctx,&fields);
  }
  else{
    /*Initialization Set initial flow field values. This case is zero. This will have to be called an initialization function*/
    ierr = VecSet(ctx.PreFlowFields,0.);CHKERRQ(ierr);
    ierr = VecSet(ctx.RHSVelPpre,0.);CHKERRQ(ierr);
    ierr = VecSet(ctx.pressure_old,0.);CHKERRQ(ierr);
    ierr = VecSet(ctx.RHSPpre,0.);CHKERRQ(ierr);
    for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nProcessing step %i.\n",ctx.timestep);CHKERRQ(ierr);
      ctx.timevalue = ctx.timestep * ctx.maxtimevalue / (ctx.maxtimestep-1.);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\ntime value %f \n",ctx.timevalue);CHKERRQ(ierr);
      ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
      ierr = FieldsH5Write(&ctx,&fields);
      /*This will have to be called "an update function"*/
      ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
      ierr = VecCopy(ctx.RHSVelP,ctx.RHSVelPpre);CHKERRQ(ierr);
      ierr = VecCopy(fields.pressure,ctx.pressure_old);CHKERRQ(ierr);
      ierr = VecCopy(ctx.RHSP,ctx.RHSPpre);CHKERRQ(ierr);
    }
    Vec error;
    PetscReal norm_1,norm_2,norm_inf;
    ierr = VecDuplicate(fields.VelnPress,&error);
    ierr = VecWAXPY(error,-1.0,fields.VelnPress,fields.FlowBCArray);
    ierr = VecNorm(error,NORM_1,&norm_1);
    ierr = VecNorm(error,NORM_2,&norm_2);
    ierr = VecNorm(error,NORM_INFINITY,&norm_inf);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n1_NORM = %f \n 2_norm = %f \n inf_norm = %f \n",norm_1, norm_2,norm_inf);CHKERRQ(ierr);
  }
  ierr = VFCheckVolumeBalance(&vol,&vol1,&vol2,&vol3,&vol4,&vol5,&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n modulus_volume = %g\n",vol);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," divergence_volume = %g\n",vol1);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," surface_flux_volume = %g\n",vol2);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," well_volume = %g\n",vol3);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," source_volume = %g\n",vol4);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," vol.strain_volume = %g\n",vol5);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Volume Balance ::: RHS = %g \t LHS = %g \n",vol+vol1,vol3+vol4+vol5);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

