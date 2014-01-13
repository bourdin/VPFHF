/*
 test38.c: 2D SNES. Single well problem. Source term implemented as dirac function and analytical solution is 2D Green's function. All pressure boundary condition.
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
 
 
 
 matrix location = i +j*nx+nx*ny*k
 
 
 ./test38  -n 101,2,101 -l 50,1,50  -theta 1 -m_inv 0 -flowsolver FLOWSOLVER_snesstandardFEM -npc 1 -pc0_r 2.5 -pc0_center 25.,0.5,25 -epsilon 3.0 -pc0_theta 0 -pc0_phi 90 -nfw 1 -fracw0_coords 25,0.5,25  -fracw0_constraint Rate -fracw0_type injector -fracw0_Qw 5e-4 -pc0_thickness 1.6 -pre 0.0001 -perm 1e-14 -permmax 1e-8 -atnum 2 -num 0 -num1 7  -timestepsize 0.01 -E 17 -nu 0.2 -m_inv 0.1 -miu 1e-10 -FlowStSnes_pc_type lu -nw 1
 
 ./test38  -n 101,2,101 -l 50,1,50  -theta 1 -m_inv 0 -flowsolver FLOWSOLVER_snesstandardFEM -npc 1 -pc0_r 2.5 -pc0_center 25.,0.5,25 -epsilon 3.0 -pc0_theta 0 -pc0_phi 90 -nfw 1 -fracw0_coords 25,0.5,25  -fracw0_constraint Rate -fracw0_type injector -fracw0_Qw 5e-4 -pc0_thickness 1.6 -pre 0.0001 -perm 1e-14 -permmax 1e-8 -atnum 2 -num 0 -num1 6  -timestepsize 0.01 -E 17 -nu 0.2 -m_inv 0.2 -miu 1e-10 -FlowStSnes_pc_type lu -nw 1
 
 
 ./test38  -n 101,2,101 -l 50,1,50  -theta 1 -m_inv 0 -flowsolver FLOWSOLVER_snesstandardFEM -npc 1 -pc0_r 3.0 -pc0_center 25.,0.5,25 -epsilon 4.0 -pc0_theta 0 -pc0_phi 90 -nfw 1 -fracw0_coords 25,0.5,25  -fracw0_constraint Rate -fracw0_type injector -fracw0_Qw 5e-1 -pc0_thickness 1.6 -pre 0.0001 -perm 1e-16 -permmax 1e-8 -atnum 2 -num 3 -num1 0  -timestepsize 0.1 -E 17 -nu 0.2 -m_inv 10 -miu 1e-13  -nw 1 -FlowStSnes_pc_type lu -Gc 0.005
 
 
 ./test38  -n 101,2,101 -l 50,1,50  -theta 1 -m_inv 0 -flowsolver FLOWSOLVER_snesmixedFEM -npc 1 -pc0_r 3.0 -pc0_center 25.,0.5,25 -epsilon 4.0 -pc0_theta 0 -pc0_phi 90 -nfw 1 -fracw0_coords 25,0.5,25  -fracw0_constraint Rate -fracw0_type injector -fracw0_Qw 5e-1 -pc0_thicknes
 
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
  PetscInt    xs1,xm1,ys1,ym1,zs1,zm1;
  Vec                 Vold;
  PetscInt            altminit=1;
  PetscReal    perm = 1e-2;
  PetscReal    p = 1e-5;
  
  
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
  
  ierr = PetscOptionsGetReal(PETSC_NULL,"-perm",&perm,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-pre",&p,PETSC_NULL);CHKERRQ(ierr);
  
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
      }
    }
  }
  
  ierr = DMDAVecRestoreArray(ctx.daScal,ctx.Source,&src_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-theta",&ctx.flowprop.theta,PETSC_NULL);CHKERRQ(ierr);
  ctx.flowprop.timestepsize = 1;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-timestepsize",&ctx.flowprop.timestepsize,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
    //  ctx.flowprop.K_dr = ctx.flowprop.M_inv*0.4;
  
  
  
  
  
  
  
  
  
  
  
  
  
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
  ctx.bcU[2].face[X0]= ZERO;
  
  ctx.bcU[0].face[X1]= ZERO;
  ctx.bcU[1].face[X1]= ZERO;
  ctx.bcU[2].face[X1]= ZERO;
  
  ctx.bcU[0].face[Z0]= ZERO;
  ctx.bcU[1].face[Z0]= ZERO;
  ctx.bcU[2].face[Z0]= ZERO;
  
  ctx.bcU[0].face[Z1]= ZERO;
  ctx.bcU[1].face[Z1]= ZERO;
  ctx.bcU[2].face[Z1]= ZERO;
  
  
  
  ctx.bcU[1].face[Y0]= ZERO;
  ctx.bcU[1].face[Y1]= ZERO;
  
  ctx.bcV[0].face[X0] = ONE;
  ctx.bcV[0].face[X1] = ONE;
  ctx.bcV[0].face[Z0] = ONE;
  ctx.bcV[0].face[Z1] = ONE;
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
  ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(ctx.U_old,0.);CHKERRQ(ierr);
  
  ctx.flowprop.alphabiot = 	ctx.matprop[0].beta = 0;									//biot's constant
  ctx.hasCrackPressure = PETSC_TRUE;
  ctx.FlowDisplCoupling = PETSC_TRUE;
  ctx.ResFlowMechCoupling = FIXEDSTRESS;
  ctx.ResFlowMechCoupling = FIXEDSTRAIN;
  ctx.FractureFlowCoupling = PETSC_TRUE;
  
  ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
  altminit = 0.;
  
  
//  PetscReal  Q_inj = ctx.fracwell[0].Qw;
  PetscReal  Q_inj;
  PetscReal  tolV = 7e-2;
  PetscReal  tolP = 7e-2;
  PetscReal  timestepsize_o = 0;
  PetscReal  time = 0;
  PetscInt    num = 10;
  PetscInt    num1 = 10;
  Vec         Pold;
  PetscReal  pmax;
  PetscReal vmax;

  PetscReal  TotalLeakOff_o = 0;
  PetscReal  timevalue_o = 0;
  PetscReal InjVolrate;

  
  ierr = VecDuplicate(fields.pressure,&Pold);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-num",&num,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-num1",&num1,PETSC_NULL);CHKERRQ(ierr);
  ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
  
  PetscReal p_old;
  PetscReal      errP=1e+10,errV=1e+10;
  PetscReal      prestol = 5.e-4;
  
  ctx.flowprop.alphabiot = 	ctx.matprop[0].beta = 0.75;									//biot's constant
  ctx.flowprop.mu = 1;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-miu",&ctx.flowprop.mu,PETSC_NULL);CHKERRQ(ierr);
  /*
   while (errP >= prestol){
   p_old = p;
   ierr = PetscPrintf(PETSC_COMM_WORLD," Time = %g \t alt min step %i with pressure %g\n",ctx.flowprop.timestepsize,altminit, p);CHKERRQ(ierr);
   ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
   ierr = VecScale(fields.U,1./p);CHKERRQ(ierr);
   ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
   p = Q_inj*ctx.flowprop.timestepsize/ctx.CrackVolume;
   ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
   ierr = VecScale(fields.U,p);CHKERRQ(ierr);
   ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
   ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);
   
   ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
   ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"     Updated crack pressure: %e\n",p);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"     Rel. change on p: %e\n",errP);CHKERRQ(ierr);
   
   errP = PetscAbs((p-p_old)/p);
   altminit++;
   
   }
   */
  
  ierr = VolumetricFractureWellRate(&InjVolrate,&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," \n\n\n    Injected Rate: %e\n\n\n\n",InjVolrate);
  Q_inj = InjVolrate;
//  Q_inj = ctx.fracwell[0].Qw;
  
  p = 1e-9;
  ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
  ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
  ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);
  ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
  ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
  ctx.timestep = 0;
  ierr = FieldsH5Write(&ctx,&fields);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"     Initial crack pressure: %e\n",p);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Initial crack volume =  %e\n",ctx.CrackVolume);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n............End of Inititalization stage\n\n\n");CHKERRQ(ierr);
  
  ctx.flowprop.alphabiot = 	ctx.matprop[0].beta = 0.75;									//biot's constant
  TotalLeakOff_o = 0;
  timevalue_o = 0;
  for(i = 1; i < num; i++){
    ctx.timestep++;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n ########### Time step =  %i \n",i);CHKERRQ(ierr);
    ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
    ierr = VF_StepU(&fields,&ctx);
    ierr = VolumetricLeakOffRate(&ctx.LeakOffRate,&ctx,&fields);CHKERRQ(ierr);
    ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
    ctx.timestep++;
    ierr = FieldsH5Write(&ctx,&fields);
    ctx.timevalue = (ctx.CrackVolume+TotalLeakOff_o-ctx.LeakOffRate*timevalue_o)/(Q_inj-ctx.LeakOffRate);
//    ctx.flowprop.timestepsize = ctx.timevalue-timevalue_o;
    ierr = PetscPrintf(PETSC_COMM_WORLD,".........Time = %i.......prior to iteration \t timevalue = %e \t errP = %e \t errV = %e \t InjVol = %e\n\n",ctx.timestep, ctx.flowprop.timestepsize, errP, errV,Q_inj);CHKERRQ(ierr);
    errP  = 1e+10;
    errV  = 1e+10;
    ierr = PetscPrintf(PETSC_COMM_WORLD," crack volume =  %e \n vol. leak-off rate =  %e \n errP = %e \t InjVol = %e\n",ctx.CrackVolume,ctx.LeakOffRate, errP,Q_inj);CHKERRQ(ierr);
    while (errV >= tolV){
      errP  = 1e+10;
      ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
      while (errP >= tolP){
        ierr = VecCopy(fields.pressure,Pold);CHKERRQ(ierr);

        ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
        ierr = VF_StepU(&fields,&ctx);
        ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
        ierr = VolumetricLeakOffRate(&ctx.LeakOffRate,&ctx,&fields);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD," crack volume =  %e \n vol. leak-off rate =  %e \n errP = %e\n",ctx.CrackVolume,ctx.LeakOffRate, errP);CHKERRQ(ierr);
        ierr = VecAXPY(Pold,-1.,fields.pressure);CHKERRQ(ierr);
        ierr = VecNorm(Pold,NORM_1,&errP);CHKERRQ(ierr);
        ierr = VecMax(fields.pressure,PETSC_NULL,&pmax);CHKERRQ(ierr);
        
        errP = errP/pmax;
        ctx.timestep++;
        ctx.timevalue = (ctx.CrackVolume+TotalLeakOff_o-ctx.LeakOffRate*timevalue_o)/(Q_inj-ctx.LeakOffRate);
//        ctx.flowprop.timestepsize = ctx.timevalue-timevalue_o;
        ierr = PetscPrintf(PETSC_COMM_WORLD,".........Time = %i.......Iteration step =  %i \t timevalue = %e \t errP = %e \t errV = %e \t InjVol = %e\n\n",i,ctx.timestep, ctx.flowprop.timestepsize, errP, errV,Q_inj);CHKERRQ(ierr);
        ierr = FieldsH5Write(&ctx,&fields);
      }
      ctx.timestep++;
      ierr = VF_StepV(&fields,&ctx);
      ierr = FieldsH5Write(&ctx,&fields);

      ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
      ierr = VecNorm(Vold,NORM_1,&errV);CHKERRQ(ierr);
      ierr = VecMax(fields.V,PETSC_NULL,&vmax);CHKERRQ(ierr);
      
      errV = errV/vmax;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n.........V STEP ERROR = %e\n\n", errV);CHKERRQ(ierr);

    }
    ierr      = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
    ierr = VecCopy(fields.pressure,ctx.pressure_old);CHKERRQ(ierr);
    ierr = VecCopy(fields.U,ctx.U_old);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSP,ctx.RHSPpre);CHKERRQ(ierr);
    timevalue_o = ctx.timevalue;
    TotalLeakOff_o = TotalLeakOff_o + ctx.LeakOffRate*ctx.flowprop.timestepsize;
  }
  /*
   ctx.timestep++;
   ierr = VecScale(fields.U,0.15);CHKERRQ(ierr);
   ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
   ierr = FieldsH5Write(&ctx,&fields);
   */
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"     \n\n\n### End of Fracture Part of Flow Example ### \n\n\n");
  ierr = VF_PermeabilityUpDate(&ctx,&fields);CHKERRQ(ierr);
  ctx.FractureFlowCoupling = PETSC_FALSE;
  ctx.hasCrackPressure = PETSC_FALSE;
  
  ctx.FlowDisplCoupling = PETSC_TRUE;
  ctx.hasFlowWells = PETSC_FALSE;
  ctx.hasFluidSources = PETSC_TRUE;
  ierr = VecCopy(ctx.RegFracWellFlowRate,ctx.Source);
  ierr = VecSet(fields.V,1.0);CHKERRQ(ierr);
  ierr = VecSet(fields.U,0.0);CHKERRQ(ierr);
  ctx.well[0].type = INJECTOR;
  ctx.well[0].Qw = 1.86e-4;
  ctx.well[0].coords[0] = 25;
  ctx.well[0].coords[1] = 0.5;
  ctx.well[0].coords[2] = 25;
  
  ierr = VecSet(ctx.pressure_old,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.RHSP,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.RHSPpre,0.);CHKERRQ(ierr);
  ierr = VecSet(Pold,0.);CHKERRQ(ierr);
  
  for(i = 1; i < num1; i++){
    errP  = 1e+10;
    while (errP >= prestol){
      ierr = VecCopy(fields.pressure,Pold);CHKERRQ(ierr);
      ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
      ierr = VF_StepU(&fields,&ctx);
      ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Initial crack volume =  %e\n\n",ctx.CrackVolume);CHKERRQ(ierr);
      ierr = VecAXPY(Pold,-1.,fields.pressure);CHKERRQ(ierr);
      ierr = VecNorm(Pold,NORM_INFINITY,&errP);CHKERRQ(ierr);
      ierr = VecMax(fields.pressure,PETSC_NULL,&pmax);CHKERRQ(ierr);
      errP = errP/pmax;
      ctx.timestep++;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n.........Time = %i.......Iteration step =  %i \t timevalue = %e \t errP = %e\n",i,ctx.timestep, ctx.flowprop.timestepsize, errP);CHKERRQ(ierr);    ierr = FieldsH5Write(&ctx,&fields);
    }
    ierr = VecCopy(fields.pressure,ctx.pressure_old);CHKERRQ(ierr);
    ierr = VecCopy(fields.U,ctx.U_old);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSP,ctx.RHSPpre);CHKERRQ(ierr);
  }
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}

