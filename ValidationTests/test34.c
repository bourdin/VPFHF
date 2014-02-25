/*
 test34.c: 1D Mandel's problem. Coupling of flow and displacmement solver for geomechanics application.
 See http://scholar.lib.vt.edu/theses/available/etd-07012008-115136/unrestricted/Dissertation_ImsooLee.pdf

 
 This example also shows the instability of the fixed-strain algorithm
 (c) 2010-2013 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
./test34 -l 5,0.01,1 -n 21,2,11 -flowsolver FLOWSOLVER_snesstandarDFEM -E 1. -nu 0.2 -maxtimestep 5 -timestepsize 0.02 -resflowmechcoupling fixedstress
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
	PetscInt		i,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm;
	PetscReal		BBmin[3],BBmax[3];
	PetscReal		***presbc_array;
	PetscReal		****coords_array;
	PetscReal		hx,hy,hz;
	PetscReal		gx,gy,gz;
	PetscReal		lx,ly,lz;
	PetscReal		gamma, beta, rho, mu;
	PetscReal		pi,dist;
  PetscReal		***pre_array;
	PetscReal   ****perm_array;
	PetscInt    xs1,xm1,ys1,ym1,zs1,zm1;
  Vec         V_hold;
  PetscReal   ****bcu_array;
	PetscReal		****velbc_array;
  PetscReal   vol,vol1,vol2,vol3,vol4,vol5;


		
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ctx.fractureflowsolver = FRACTUREFLOWSOLVER_NONE;
	ctx.flowsolver = FLOWSOLVER_KSPMIXEDFEM;
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx.daVFperm,&xs1,&ys1,&zs1,&xm1,&ym1,&zm1);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx.daScal,&V_hold);CHKERRQ(ierr);
	ierr = VecSet(ctx.PresBCArray,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.VelBCArray,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.Source,0.);CHKERRQ(ierr);
	ctx.hasFlowWells = PETSC_TRUE;
	ctx.hasFluidSources = PETSC_FALSE;
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  ierr = VecSet(fields.BCU,0.0);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx.daVect,fields.BCU,&bcu_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);
  for (k = zs1; k < zs1+zm1; k++) {
    for (j = ys1; j < ys1+ym1; j++) {
      for (i = xs1; i < xs1+xm1; i++) {
        perm_array[k][j][i][0] = 1;
        perm_array[k][j][i][1] = 0;
        perm_array[k][j][i][2] = 1;
        perm_array[k][j][i][3] = 0.;
        perm_array[k][j][i][4] = 0.;
        perm_array[k][j][i][5] = 0.;
      }
    }
  }
  ierr = DMDAVecRestoreArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);

								
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
  ctx.bcQ[0].face[X0] = FIXED;
  ctx.bcQ[0].face[X1] = FIXED;
  ctx.bcQ[1].face[Y0] = FIXED;
  ctx.bcQ[1].face[Y1] = FIXED;
  ctx.bcQ[2].face[Z0] = FIXED;
	ctx.bcQ[2].face[Z1] = FIXED;
  ierr = VecSet(ctx.VelBCArray,0);CHKERRQ(ierr);
  ierr = VecSet(ctx.PresBCArray,0);CHKERRQ(ierr);
 ctx.flowprop.alphabiot = 	ctx.matprop[0].beta = 1;									//biot's constant
 ctx.flowprop.theta = 1.;
  ctx.flowprop.M_inv = 0.;
 ctx.matprop[0].nu = 0.2;										//Poisson's ratio's modulus
 ctx.flowprop.rho = 1;									 //density in lb/ft^3
  ctx.flowprop.mu = 0.2;                    //viscosity in cp
 ctx.flowprop.cf = 1.;										  //compressibility in psi^{-1}
 ctx.flowprop.beta = 1.;										  //flow rate conversion constant
 ctx.flowprop.gamma = 1.;										  //pressue conversion constant
 ctx.flowprop.alpha = 1.;										  //volume conversion constatnt
 ctx.flowprop.g[0] = 0.;										  //x-component of gravity. unit is ft/s^2
 ctx.flowprop.g[1] = 0.;										  //y-component of gravity. unit is ft/s^2
 ctx.flowprop.g[2] = 0.;								  //z-component of gravity. unit is ft/s^2
  ctx.flowprop.K_dr = ctx.matprop[0].E/(2*(1+ctx.matprop[0].nu)*(1-2*ctx.matprop[0].nu));   // For 2D


  /*
   Beginning of mechanical part of code*/
    /*
	 Reset all BC for U and V
	 */
	for (i = 0; i < 6; i++) {
		ctx.bcV[0].face[i]=NONE;
		for (j = 0; j < 3; j++) {
			ctx.bcU[j].face[i] = NONE;
		}
	}
	for (i = 0; i < 12; i++) {
		ctx.bcV[0].edge[i]=NONE;
		for (j = 0; j < 3; j++) {
			ctx.bcU[j].edge[i] = NONE;
		}
	}
	for (i = 0; i < 8; i++) {
		ctx.bcV[0].vertex[i]=NONE;
		for (j = 0; j < 3; j++) {
			ctx.bcU[j].vertex[i] = NONE;
		}
	}
    //We will use Z and X direction
  /*	face X0	*/
  ctx.bcU[0].face[X0]= ZERO;
  ctx.bcU[1].face[X0]= ZERO;
  /*	face X1	*/
    //  ctx.bcU[0].face[X1]= ZERO;
  ctx.bcU[1].face[X1]= ZERO;
  /*	face Y0	*/
  ctx.bcU[1].face[Y0]= ZERO;
  /*	face Y1	*/
  ctx.bcU[1].face[Y1]= ZERO;
  /*	face Z0	*/
  ctx.bcU[1].face[Z0]= ZERO;
  ctx.bcU[2].face[Z0]= ZERO;
  /*	face Z1	*/
  ctx.bcU[1].face[Z1]= ZERO;
  ctx.bcU[2].face[Z1]= NONE;
  for (i = 0; i < 6; i++) {
    ctx.insitumax[i] = 0.;
    ctx.insitumin[i] = 0.;
  }
  ctx.insitumax[2] = -1.0;
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,fields.BCU,&bcu_array);CHKERRQ(ierr);
  ctx.hasInsitu        = PETSC_TRUE;
  ctx.FlowDisplCoupling = PETSC_TRUE;
	ctx.hasCrackPressure = PETSC_TRUE;
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);

	ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.pressure,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.V,1.0);CHKERRQ(ierr);
  ierr = VecSet(ctx.U_old,0.);CHKERRQ(ierr);
 /*End of mechanical part of code*/

  ierr = PetscOptionsGetReal(PETSC_NULL,"-theta",&ctx.flowprop.theta,PETSC_NULL);CHKERRQ(ierr);
  ctx.flowprop.timestepsize = 0.;

  PetscReal displ_p_tol = 1e-6;
  Vec error;
  Vec PreIteSol;
	PetscReal norm_1 = 1e+3,norm_2 = 1e+3,norm_inf = 1e+3;
  PetscInt displ_iter = 0;
  ierr = VecDuplicate(fields.pressure,&error);
  ierr = VecDuplicate(fields.pressure,&PreIteSol);
  
  /*Initialization Set initial flow field values. This case is zero. This will have to be called an initialization function*/
  
  ierr = VecSet(ctx.PreFlowFields,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.RHSVelPpre,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.pressure_old,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.RHSPpre,0.);CHKERRQ(ierr);
  ctx.timestep = 0;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  Computing initial time step solution\n");CHKERRQ(ierr);
    ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
    ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
    ierr = VecCopy(fields.pressure,PreIteSol);CHKERRQ(ierr);
    while (norm_inf > displ_p_tol) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"   Iteration Step: %d\n",displ_iter);CHKERRQ(ierr);
      ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
      ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
      displ_iter++;
      ierr = VecWAXPY(error,-1.0,fields.pressure,PreIteSol);
      ierr = VecCopy(fields.pressure,PreIteSol);CHKERRQ(ierr);
      ierr = VecNorm(error,NORM_INFINITY,&norm_inf);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n inf_norm = %f \n",norm_inf);CHKERRQ(ierr);
    }
    ierr = FieldsH5Write(&ctx,&fields);
  
  

  ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
  ierr = VecCopy(ctx.RHSVelP,ctx.RHSVelPpre);CHKERRQ(ierr);
  ierr = VecCopy(fields.pressure,ctx.pressure_old);CHKERRQ(ierr);
  ierr = VecCopy(ctx.RHSP,ctx.RHSPpre);CHKERRQ(ierr);
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
  ctx.bcQ[0].face[X0] = FIXED;
  ctx.bcP[0].face[X1] = FIXED;
  ctx.bcQ[1].face[Y0] = FIXED;
  ctx.bcQ[1].face[Y1] = FIXED;
  ctx.bcQ[2].face[Z0] = FIXED;
  ctx.bcQ[2].face[Z1] = FIXED;
  ierr = VecSet(ctx.PresBCArray,0);CHKERRQ(ierr);
  ierr = VecSet(ctx.VelBCArray,0);CHKERRQ(ierr);

  ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
  ctx.flowprop.M_inv = 1.0;

  
  ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
  ierr = VecSet(ctx.RHSVelPpre,0.);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-theta",&ctx.flowprop.theta,PETSC_NULL);CHKERRQ(ierr);
  ctx.flowprop.timestepsize = 1.;
	ierr = PetscOptionsGetReal(PETSC_NULL,"-timestepsize",&ctx.flowprop.timestepsize,PETSC_NULL);CHKERRQ(ierr);
  ctx.maxtimestep = 10;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-maxtimestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr);
  for (ctx.timestep = 1; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"      ##########################################################\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"      #                                                        #\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"      #                          STAGE %d!!!                    #\n",ctx.timestep );CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"      #                                                        #\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"      #        Start of drained consolidation steps            #\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"      #                                                        #\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"      ##########################################################\n\n\n");CHKERRQ(ierr);

    
    ierr = VecCopy(fields.U,ctx.U_old);CHKERRQ(ierr);
    norm_inf = 1e+3;
    displ_iter = 0;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  Computing solution at %d time step solution \n", ctx.timestep);CHKERRQ(ierr);
    ierr = VecCopy(fields.pressure,PreIteSol);CHKERRQ(ierr);
  while (norm_inf > displ_p_tol) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  Step %d, Iteration Step: %d\n",ctx.timestep, displ_iter);CHKERRQ(ierr);
    ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
    ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
    displ_iter++;
    ierr = VecWAXPY(error,-1.0,fields.pressure,PreIteSol);
    ierr = VecCopy(fields.pressure,PreIteSol);CHKERRQ(ierr);
    ierr = VecNorm(error,NORM_INFINITY,&norm_inf);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n inf_norm = %f \n",norm_inf);CHKERRQ(ierr);
    }
    
    ierr = VFCheckVolumeBalance(&vol,&vol1,&vol2,&vol3,&vol4,&vol5,&ctx,&fields);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n modulus_volume = %g\n",vol);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," divergence_volume = %g\n",vol1);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," surface_flux_volume = %g\n",vol2);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," well_volume = %g\n",vol3);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," source_volume = %g\n",vol4);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," vol.strain_volume = %g\n",vol5);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Volume Balance ::: RHS = %g \t LHS = %g \n",vol+vol1,vol3+vol4+vol5);CHKERRQ(ierr);

    ierr = FieldsH5Write(&ctx,&fields);
    ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSVelP,ctx.RHSVelPpre);CHKERRQ(ierr);
    ierr = VecCopy(fields.pressure,ctx.pressure_old);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSP,ctx.RHSPpre);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&PreIteSol);CHKERRQ(ierr);
  ierr = VecDestroy(&error);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

