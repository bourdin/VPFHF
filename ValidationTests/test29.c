/*
 test29.c: 1D. Test for transient solver
 http://tutorial.math.lamar.edu/Classes/DE/SolvingHeatEquation.aspx
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
 ./test29 -n 51,51,2 -l 1,1,0.01 -m_inv 1 -ts_dt 0.005 -ts_max_steps 20
 ./test29 -n 51,51,2 -l 1,1,0.01 -m_inv 1 -maxtimestep 20 -flowsolver FLOWSOLVER_snesMIXEDFEM -timestepsize 0.005
 
 */

#include "petsc.h"
#include "CartFE.h"
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
	PetscInt		i,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm,xs1,xm1,ys1,ym1,zs1,zm1;
	PetscReal		BBmin[3],BBmax[3];
	PetscReal		***presbc_array;
	PetscReal		****velbc_array;
	PetscReal		***pres_ini_array;
	PetscReal		****vel_ini_array;
	PetscReal		***src_array;
	PetscReal		****coords_array;
	PetscReal		hx,hy,hz;
	PetscReal		lx,ly,lz;
	PetscReal		gx,gy,gz;
	PetscReal		gamma, beta, rho, mu;
	PetscReal		pi;
	PetscReal		****perm_array;
	PetscReal		****velnpre_array;
		
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ctx.fractureflowsolver = FRACTUREFLOWSOLVER_NONE;
	ctx.flowsolver = FLOWSOLVER_TSMIXEDFEM;
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	ierr = VecSet(fields.velocity,0.);CHKERRQ(ierr);
	ierr = VecSet(fields.pressure,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.PresBCArray,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.VelBCArray,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.Source,0.);CHKERRQ(ierr);
	ctx.hasFluidSources = PETSC_TRUE;
	ctx.hasFlowWells = PETSC_FALSE;	
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);	
	ierr = DMDAVecGetArray(ctx.daScal,fields.pressure,&pres_ini_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVect,fields.velocity,&vel_ini_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daFlow,fields.VelnPress,&velnpre_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx.daScal,ctx.Source,&src_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr); 
	ierr = DMDAGetCorners(ctx.daVFperm,&xs1,&ys1,&zs1,&xm1,&ym1,&zm1);CHKERRQ(ierr);
	for (k = zs1; k < zs1+zm1; k++) {
		for (j = ys1; j < ys1+ym1; j++) {
				for (i = xs1; i < xs1+xm1; i++) {
				perm_array[k][j][i][2] = 0.;
				perm_array[k][j][i][3] = 0.;
				perm_array[k][j][i][4] = 0.;			
				perm_array[k][j][i][5] = 0.;					
			}
		}
	}
	pi = 6.*asin(0.5);
	lz = BBmax[2]-BBmin[2];
	ly = BBmax[1]-BBmin[1];
	lx = BBmax[0]-BBmin[0];	
	hx = lx/(nx-1);
	hy = ly/(nx-1);
	hz = lz/(nz-1);	
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
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				vel_ini_array[k][j][i][0] = 0.;
				vel_ini_array[k][j][i][1] = 0.;
				vel_ini_array[k][j][i][2] = 0.;
				pres_ini_array[k][j][i] = 6.*sin(pi*hx*i/lx);
			}
		}
	}		

	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				velnpre_array[k][j][i][0] = vel_ini_array[k][j][i][0];
				velnpre_array[k][j][i][1] = vel_ini_array[k][j][i][1];
				velnpre_array[k][j][i][2] = vel_ini_array[k][j][i][2];
				velnpre_array[k][j][i][3] = pres_ini_array[k][j][i];
			}
		}
	}			
	ctx.bcP[0].face[X0] = FIXED;
	ctx.bcP[0].face[X1] = FIXED;
	ctx.bcQ[1].face[Y0] = FIXED;
	ctx.bcQ[1].face[Y1] = FIXED;
	ctx.bcQ[2].face[Z0] = FIXED;
	ctx.bcQ[2].face[Z1] = FIXED;
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				velbc_array[k][j][i][0] = 0.;
				velbc_array[k][j][i][1] = 0.;
				velbc_array[k][j][i][2] = 0.;
				presbc_array[k][j][i] = 0.;
			}
		}
	}	
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				src_array[k][j][i] = 0.;
			}
		}
	}
	ierr = DMDAVecRestoreArrayDOF(ctx.daFlow,fields.VelnPress,&velnpre_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx.daScal,ctx.Source,&src_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx.daScal,fields.pressure,&pres_ini_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,fields.velocity,&vel_ini_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);	
	/* Setting time parameters	*/
  ctx.maxtimestep = 20;
	ctx.flowprop.timestepsize = 1;
	ierr = PetscOptionsGetReal(PETSC_NULL,"-theta",&ctx.flowprop.theta,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-timestepsize",&ctx.flowprop.timestepsize,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
	ctx.maxtimestep = 1;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-maxtimestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr);
  if(ctx.flowsolver == FLOWSOLVER_TSMIXEDFEM || ctx.flowsolver == FLOWSOLVER_TS){
  	ctx.maxtimestep = 20;
    ctx.maxtimevalue = 100.;
    ctx.timevalue = 0.1;
    ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
    ierr = FieldsH5Write(&ctx,&fields);
  }
  else{
    /*Initialization Set initial flow field values. This case is zero. This will have to be called an initialization function*/
    ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
    ierr = VecSet(ctx.RHSVelPpre,0.);CHKERRQ(ierr);
    for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nProcessing step %i.\n",ctx.timestep);CHKERRQ(ierr);
      ctx.timevalue = ctx.timestep * ctx.maxtimevalue / (ctx.maxtimestep-1.);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\ntime value %f \n",ctx.timevalue);CHKERRQ(ierr);
      ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
      ierr = FieldsH5Write(&ctx,&fields);
      /*This will have to be called "an update function"*/
      ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
      ierr = VecCopy(ctx.RHSVelP,ctx.RHSVelPpre);CHKERRQ(ierr);
      
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
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

