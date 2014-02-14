/*
 test50.c: 2D SNES. Single well problem. Source term implemented as dirac function and analytical solution is 2D Green's function. All pressure boundary condition.
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
 
 ./test50  -n 101,2,101 -l 4,0.04,4 -maxtimestep 1 timestepsize 1 -theta 1 -m_inv 0 -flowsolver FLOWSOLVER_snesstandardFEM -npc 1 -pc0_r 0.5 -pc0_center 2.,0.02,2 -epsilon 0.28 -pc0_theta 0 -pc0_phi 50 -nfw 1 -fracw0_coords 2,0.02,2.  -fracw0_constraint Rate -fracw0_type injector -fracw0_Qw 0.001 -pc0_thickness 0.1 -pre 0.0001 -perm 5e-7 -FlowStSnes_pc_type lu -atnum 2 -num 3
 
 
 ./test50  -n 101,2,101 -l 4,0.04,4 -maxtimestep 1 timestepsize 1 -theta 1 -m_inv 0 -flowsolver FLOWSOLVER_snesstandardFEM -npc 1 -pc0_r 2.5 -pc0_center 2.,0.02,2 -epsilon 0.28 -pc0_theta 0 -pc0_phi 30 -nw 1 -w0_coords 1,0.02,2.  -w0_constraint Rate -w0_type injector -w0_Qw 0.00 -pc0_thickness 0.08 -pre 0.2 -perm 1e-3 -atnum 2 -num 3 -FlowSnes_pc_type lu -permmax 0.34
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
  PetscViewer    viewer;
  PetscViewer     logviewer;
  char      filename[FILENAME_MAX];
  PetscInt    i,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal    BBmin[3],BBmax[3];
  PetscReal    ***presbc_array;
  PetscReal    ***src_array;
  PetscReal    ****coords_array;
  PetscReal    hx,hy,hz;
  PetscReal    gx,gy,gz;
  PetscReal    lx,ly,lz;
  PetscReal    gamma, beta, rho, mu;
  PetscReal    pi,dist;
  PetscReal    ***pre_array;
  PetscReal   ****perm_array;
  PetscReal   ****velbc_array;
  PetscInt    xs1,xm1,ys1,ym1,zs1,zm1;
  Vec                 Vold;
  PetscInt            altminit=1;
  PetscReal    perm = 1e-2;
  PetscReal    pre = 1e-2;
  
  
  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ctx.flowsolver = FLOWSOLVER_SNESMIXEDFEM;
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
  ierr = VecSet(ctx.PresBCArray,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.VelBCArray,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.Source,0.);CHKERRQ(ierr);
  ctx.hasFlowWells = PETSC_TRUE;
  ctx.hasFluidSources = PETSC_FALSE;
  ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx.daScal,ctx.Source,&src_array);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx.daVFperm,&xs1,&ys1,&zs1,&xm1,&ym1,&zm1);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);

  ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-perm",&perm,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-pre",&pre,PETSC_NULL);CHKERRQ(ierr);
  
  for (k = zs1; k < zs1+zm1; k++) {
    for (j = ys1; j < ys1+ym1; j++) {
      for (i = xs1; i < xs1+xm1; i++) {
        perm_array[k][j][i][0] = perm;
        perm_array[k][j][i][1] = perm;
        perm_array[k][j][i][2] = perm;
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
  
//  ctx.bcQ[0].face[X0] = FIXED;
  ctx.bcP[0].face[X0] = FIXED;
  ctx.bcP[0].face[X1] = FIXED;
  ctx.bcP[0].face[Z0] = FIXED;
  ctx.bcP[0].face[Z1] = FIXED;
  ctx.bcQ[1].face[Y0] = FIXED;
  ctx.bcQ[1].face[Y1] = FIXED;
  for (k = zs; k < zs+zm; k++) {
    for (j = ys; j < ys+ym; j++) {
      for (i = xs; i < xs+xm; i++) {
        presbc_array[k][j][i] = 0.;
        velbc_array[k][j][i][0] = 0.1;
      }
    }
  }
  
  ierr = DMDAVecRestoreArray(ctx.daScal,ctx.Source,&src_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-theta",&ctx.flowprop.theta,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-timestepsize",&ctx.flowprop.timestepsize,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
  
  
  
  
  
  
  
  
  
  
  
  
  
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
//  ctx.bcU[2].face[X0]= ZERO;
  
  ctx.bcU[0].face[X1]= ZERO;
  ctx.bcU[1].face[X1]= ZERO;
//  ctx.bcU[2].face[X1]= ZERO;
  
  ctx.bcU[0].face[Z0]= ZERO;
  ctx.bcU[1].face[Z0]= ZERO;
  ctx.bcU[2].face[Z0]= ZERO;
  
  ctx.bcU[0].face[Z1]= ZERO;
  ctx.bcU[1].face[Z1]= ZERO;
  ctx.bcU[2].face[Z1]= ZERO;
  
  
  
  ctx.bcU[1].face[Y0]= ZERO;
  ctx.bcU[1].face[Y1]= ZERO;
  
//  ctx.bcV[0].face[X0] = ONE;
//  ctx.bcV[0].face[X1] = ONE;
//  ctx.bcV[0].face[Z0] = ONE;
//  ctx.bcV[0].face[Z1] = ONE;
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,pre);CHKERRQ(ierr);
  ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(ctx.U_old,0.);CHKERRQ(ierr);

  ctx.flowprop.alphabiot = 	ctx.matprop[0].beta = 0;									//biot's constant
  ctx.hasCrackPressure = PETSC_TRUE;
  ctx.FlowDisplCoupling = PETSC_FALSE;
  ctx.ResFlowMechCoupling == FIXEDSTRAIN;
  ctx.FractureFlowCoupling = PETSC_TRUE;
  
  ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
  altminit = 0.;
  
  
  PetscReal  Q_inj = 1.42e-7;
  PetscReal  tolV = 1e-4;
  PetscReal  tolP = 1e-4;
  PetscReal  timestepsize_o = 0;
  PetscReal  time = 0;
  PetscReal  time_o = 0;
  PetscReal  vol_o = 0;
  PetscReal      errV,errP;
  PetscInt    num = 10;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-num",&num,PETSC_NULL);CHKERRQ(ierr);
  ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
  ctx.timestep = 0;

  
  ctx.flowprop.alphabiot = 	ctx.matprop[0].beta = 0;									//biot's constant
  ierr = VF_StepV(&fields,&ctx);
  ierr = VF_StepU(&fields,&ctx);
  ierr = VecCopy(fields.U,ctx.U_old);
  ierr = VecCopy(fields.V,ctx.V_old);
  ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
  ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
  ctx.timestep++;
  ierr = FieldsH5Write(&ctx,&fields);


  Vec V_o;
  ierr = VecDuplicate(fields.V,&V_o);CHKERRQ(ierr);
  ierr = VecCopy(fields.V,V_o);CHKERRQ(ierr);
  ierr = VecCopy(fields.V,V_o);CHKERRQ(ierr);
  ierr = VF_PermeabilityUpDate(&ctx,&fields);CHKERRQ(ierr);
  ierr = VecSet(fields.V,1.0);CHKERRQ(ierr);
//  ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
  ierr = VecCopy(V_o,fields.V);CHKERRQ(ierr);
  ctx.timestep++;
  ierr = FieldsH5Write(&ctx,&fields);

  
 
    
  
  
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}
