/*
 test30.c: 2D. Terzaghi problem. Coupling of flow and displacmement solver for geomechanics application.
 
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
./test30 -l 2,1,3 -n 11,2,11 -flowsolver FLOWSOLVER_snesMIXEDFEM
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

		
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ctx.flowsolver = FLOWSOLVER_KSPMIXEDFEM;
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx.daScal,&V_hold);CHKERRQ(ierr);
	ierr = VecSet(ctx.PresBCArray,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.VelBCArray,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.Source,0.);CHKERRQ(ierr);
	ctx.hasFlowWells = PETSC_TRUE;
	ctx.hasFluidSources = PETSC_FALSE;
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daVFperm,&xs1,&ys1,&zs1,&xm1,&ym1,&zm1);CHKERRQ(ierr);	
	ierr = DMDAVecGetArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
  for (k = zs1; k < zs1+zm1; k++) {
    for (j = ys1; j < ys1+ym1; j++) {
      for (i = xs1; i < xs1+xm1; i++) {
//        perm_array[k][j][i][0] = 2.e-13;
//        perm_array[k][j][i][1] = 2.e-13;
//        perm_array[k][j][i][2] = 2.e-13;
        perm_array[k][j][i][0] = 1;
        perm_array[k][j][i][1] = 1;
        perm_array[k][j][i][2] = 0.;
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
  ctx.bcQ[2].face[Z0] = FIXED;
	ctx.bcQ[2].face[Z1] = FIXED;
    for (k = zs; k < zs+zm; k++) {
      for (j = ys; j < ys+ym; j++) {
        for (i = xs; i < xs+xm; i++) {
          if(i == 0){
            velbc_array[k][j][i][0] = 0.;
            presbc_array[k][j][i] = 0.;
        }
          if(i == nx-1){
            velbc_array[k][j][i][0] = 0.;
            presbc_array[k][j][i] = 0.;
          }
          if(k == 0){
            velbc_array[k][j][i][2] = 0.;
            presbc_array[k][j][i] = 0.;
          }
          if(k == nz-1){
            velbc_array[k][j][i][2] = 0.;
            presbc_array[k][j][i] = 0.;
          }
      }
    }  
  }
  
  
  
  
  
/*
  ctx.flowprop.alphabiot = 	ctx.matprop[0].beta = .9;									//biot's constant
  ctx.matprop[0].rho = ctx.flowprop.rho = 2.2;
  ctx.flowprop.theta = 1.;
	ctx.matprop[0].E = 100;										//Young's modulus
	ctx.flowprop.M_inv = 1./100.0;
	ctx.matprop[0].nu = 0.35;										//Poisson's ratio's modulus
    //	ctx.flowprop.Kf = 0.1;
    //	ctx.flowprop.Ks = 0.001;
    //	ctx.flowprop.por = 0.2;
    //	ctx.flowprop.rho = 940e-6;									 //density in lb/ft^3
    //	ctx.flowprop.rhor = 2000e-6;									 //rock density in lb/ft^3
	ctx.flowprop.mu = 10;                    //viscosity in cp
	ctx.flowprop.cf = 1.;										  //compressibility in psi^{-1}
	ctx.flowprop.beta = 1.;										  //flow rate conversion constant
	ctx.flowprop.gamma = 1.;										  //pressue conversion constant
	ctx.flowprop.alpha = 1.;										  //volume conversion constatnt
  ctx.flowprop.g[0] = 0.;										  //x-component of gravity. unit is ft/s^2
  ctx.flowprop.g[1] = 0.;										  //y-component of gravity. unit is ft/s^2
  ctx.flowprop.g[2] = 0.;								  //z-component of gravity. unit is ft/s^2
*/

  
  
 
 ctx.flowprop.alphabiot = 	ctx.matprop[0].beta = .79;									//biot's constant
 ctx.flowprop.theta = 1.;
 ctx.matprop[0].E = 1.44e4;										//Young's modulus
 ctx.flowprop.M_inv = 1./12300.0;
 ctx.matprop[0].nu = 0.2;										//Poisson's ratio's modulus
 //	ctx.flowprop.Kf = 0.1;
 //	ctx.flowprop.Ks = 0.001;
 //	ctx.flowprop.por = 0.2;
 //	ctx.flowprop.rho = 940e-6;									 //density in lb/ft^3
 ctx.flowprop.rho = 1788e-6;									 //density in lb/ft^3
 //	ctx.flowprop.rhor = 2000e-6;									 //rock density in lb/ft^3
 ctx.flowprop.mu = 1.222e-07;                    //viscosity in cp
 ctx.flowprop.cf = 1.;										  //compressibility in psi^{-1}
 ctx.flowprop.beta = 1.;										  //flow rate conversion constant
 ctx.flowprop.gamma = 1.;										  //pressue conversion constant
 ctx.flowprop.alpha = 1.;										  //volume conversion constatnt
 ctx.flowprop.g[0] = 0.;										  //x-component of gravity. unit is ft/s^2
 ctx.flowprop.g[1] = 0.;										  //y-component of gravity. unit is ft/s^2
 ctx.flowprop.g[2] = 0.;								  //z-component of gravity. unit is ft/s^2

 

  
  
  
  /*Beginning of mechanical part of code*/
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
  ierr = VecSet(fields.BCU,0.0);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx.daVect,fields.BCU,&bcu_array);CHKERRQ(ierr);

//We will use Z and X direction
  /*	face X0	*/
  ctx.bcU[0].face[X0]= ZERO;
  ctx.bcU[1].face[X0]= ZERO;
  /*	face X1	*/
  ctx.bcU[0].face[X1]= ZERO;
  ctx.bcU[1].face[X1]= ZERO;
  /*	face Y0	*/
  ctx.bcU[1].face[Y0]= ZERO;
  /*	face Y1	*/
  ctx.bcU[1].face[Y1]= ZERO;
  /*	face Z0	*/
  ctx.bcU[0].face[Z0]= ZERO;
  ctx.bcU[1].face[Z0]= ZERO;
  ctx.bcU[2].face[Z0]= ZERO;
  /*	face Z1	*/
  ctx.bcU[0].face[Z1]= ZERO;
  ctx.bcU[1].face[Z1]= ZERO;
  ctx.bcU[2].face[Z1]= NONE;
  for (k = zs; k < zs + zm; k++){
		for (j = ys; j < ys+ym; j++){
			for (i = xs; i < xs+xm; i++) {
        if(k == nz-1){
//          bcu_array[k][j][i][2] = -5.1e-4;
        }
			}
		}
	}
  for (i = 0; i < 6; i++) {
    ctx.insitumax[i] = 0.;
    ctx.insitumin[i] = 0.;
  }
  ctx.insitumax[2] = ctx.insitumin[2] = -4.;

  ctx.hasInsitu        = PETSC_TRUE;
  ctx.FlowDisplCoupling = PETSC_TRUE;
	ctx.hasCrackPressure = PETSC_FALSE;
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);

	ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.pressure,0);CHKERRQ(ierr);
	ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.V,1.0);CHKERRQ(ierr);
  ierr = VecSet(ctx.U_old,0.);CHKERRQ(ierr);
 /*End of mechanical part of code*/

	ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-theta",&ctx.flowprop.theta,PETSC_NULL);CHKERRQ(ierr);
  ctx.flowprop.timestepsize = 0.;
	ierr = PetscOptionsGetReal(PETSC_NULL,"-timestepsize",&ctx.flowprop.timestepsize,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
	ctx.maxtimestep = 1;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-maxtimestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr);
  
  
  PetscReal displ_p_tol = 1e1;
  Vec error;
  Vec PreIteSol;
	PetscReal norm_1 = 1e+3,norm_2 = 1e+3,norm_inf = 1e+3;
  PetscInt displ_iter = 0;
  ierr = VecDuplicate(fields.VelnPress,&error);
  ierr = VecDuplicate(fields.VelnPress,&PreIteSol);

  
  
  /*Initialization Set initial flow field values. This case is zero. This will have to be called an initialization function*/
  ierr = VecSet(ctx.PreFlowFields,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.RHSVelPpre,0.);CHKERRQ(ierr);
  
  ctx.timestep = 0;
  ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
//  ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
  ierr = FieldsH5Write(&ctx,&fields);
  
  ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
  ierr = VecCopy(ctx.RHSVelP,ctx.RHSVelPpre);CHKERRQ(ierr);
  /*
	for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
    ierr = PetscPrintf(PETSC_COMM_WORLD,"   Step: %d\n",ctx.timestep);CHKERRQ(ierr);
    
    
    ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
    ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
    ierr = VecCopy(fields.VelnPress,PreIteSol);CHKERRQ(ierr);
    
    while (norm_inf > displ_p_tol) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"   Iteration Step: %d\n",displ_iter);CHKERRQ(ierr);
      ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
      ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
      displ_iter++;
      ierr = VecWAXPY(error,-1.0,fields.VelnPress,PreIteSol);
      ierr = VecCopy(fields.VelnPress,PreIteSol);CHKERRQ(ierr);
      ierr = VecNorm(error,NORM_1,&norm_1);
      ierr = VecNorm(error,NORM_2,&norm_2);
      ierr = VecNorm(error,NORM_INFINITY,&norm_inf);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n1_NORM = %f \n 2_norm = %f \n inf_norm = %f \n",norm_1, norm_2,norm_inf);CHKERRQ(ierr);
    }
    
    
    ierr = FieldsH5Write(&ctx,&fields);
//    This will have to be called "an update function"
    ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSVelP,ctx.RHSVelPpre);CHKERRQ(ierr);
	}
*/
  
  ierr = VecDestroy(&PreIteSol);CHKERRQ(ierr);
  ierr = VecDestroy(&error);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,fields.BCU,&bcu_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

