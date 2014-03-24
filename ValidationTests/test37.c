/*
 test28.c: 1D. Flow problem with source term = 1 and Homogeneous pressure boundary conditions on all sides. Analytical solution is p = x(x-1)/2
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
./test38 -n 51,2,11 -l 1,0.01,1 -m_inv 10 -flowsolver FLOWSOLVER_snesstandarDFEM -perm 1e-16 -miu 1e-16 -num 10 -timestepsize 1
 
 ./test38 -n 101,2,101 -l 50,1,50 -m_inv 10 -E 17 -nu 0.2 -npc 1 -pc0_r 3.0 -pc0_center 25.,0.5,25 -epsilon 4.0 -pc0_theta 0 -pc0_phi 90  -atnum 2 -pc0_thickness 0.7 -nfw 1 -fracw0_coords 25,0.5,25 -fracw0_constraint Rate -fracw0_type injector -fracw0_Qw 5e-1  -flowsolver FLOWSOLVER_snesstandarDFEM -perm 1e-16 -miu 1e-13 -num 10 -timestepsize 0.1 -FlowStSnes_pc_type lu -Gc 0.005
 
 
 ./test38 -n 101,2,101 -l 50,1,50 -m_inv 10 -E 17 -nu 0.2 -npc 1 -pc0_r 3.0 -pc0_center 25.,0.5,25 -epsilon 4.0 -pc0_theta 0 -pc0_phi 90  -atnum 2 -pc0_thickness 1.6 -nfw 1 -fracw0_coords 25,0.5,25 -fracw0_constraint Rate -fracw0_type injector -fracw0_Qw 5e-3  -flowsolver FLOWSOLVER_snesstandarDFEM -perm 1e-13 -miu 1e-13 -num 10 -timestepsize 0.10  -Gc 0.005 -FlowStSnes_ksp_type gmres
 
 
 ./test38 -n 101,2,101 -l 50,1,50 -m_inv 20 -E 17 -nu 0.2 -npc 1 -pc0_r 3.0 -pc0_center 25.,0.5,25 -epsilon 2.0 -pc0_theta 0 -pc0_phi 90  -atnum 1 -pc0_thickness 1.6 -nfw 1 -fracw0_coords 25,0.5,25 -fracw0_constraint Rate -fracw0_type injector -fracw0_Qw 5e-4  -flowsolver FLOWSOLVER_snesstandarDFEM -perm 1e-16 -miu 1e-14 -num 4 -timestepsize .5  -Gc 5e-3 -FlowstSnes_ksp_type gmres
 
 
 
 
 
 ./test38 -n 201,2,201 -l 50,1,50 -m_inv 20 -E 17 -nu 0.2 -npc 1 -pc0_r 3.0 -pc0_center 25.,0.5,25 -epsilon 2.0 -pc0_theta 0 -pc0_phi 90  -atnum 1 -pc0_thickness 1.1 -nfw 1 -fracw0_coords 25,0.5,25 -fracw0_constraint Rate -fracw0_type injector -fracw0_Qw 5e-2  -flowsolver FLOWSOLVER_snesstandarDFEM -perm 1e-16 -miu 1e-13 -num 1000 -timestepsize .1  -Gc 5e-3
 
 
 
 ./test38 -n 200,2,201 -l 50,1,50 -m_inv 1 -E 17 -nu 0.2 -npc 1 -pc0_r 3.0 -pc0_center 25.,0.5,25 -epsilon 1.6 -pc0_theta 0 -pc0_phi 90  -atnum 2 -pc0_thickness 0.3 -nfw 1 -fracw0_coords 25,0.5,25 -fracw0_constraint Rate -fracw0_type injector -fracw0_Qw 5e-2  -flowsolver FLOWSOLVER_snesstandarDFEM -perm 1e-16 -miu 1e-13 -num 1000 -timestepsize 0.1  -Gc 5e-3

 scp cchukw1@head.schur.math.lsu.edu:/archive/cchukw1/new5/TEST.xmf TEST.xmf
 
 scp cchukw1@stampede.tacc.utexas.edu:/home1/01968/cchukw1/vf-chevron/new/TEST.xmf TEST.xmf
 
 ./test38 -n 201,2,201 -l 100,1,100 -m_inv 1 -E 17 -nu 0.2 -npc 2 -pc0_r 5.0 -pc0_center 45.,0.5,50 -epsilon 3.2 -pc0_theta 0 -pc0_phi 90  -atnum 2 -pc0_thickness 1.2 -pc1_r 5.0 -pc1_center 55.,0.5,50 -pc1_theta 0 -pc1_phi 90 -pc1_thickness 1.2 -nfw 2 -fracw0_coords 45,0.5,50 -fracw0_constraint Rate -fracw0_type injector -fracw0_Qw 5e-2 -fracw1_coords 55,0.5,50 -fracw1_constraint Rate -fracw1_type injector -fracw1_Qw 5e-2 -flowsolver FLOWSOLVER_snesstandarDFEM -perm 1e-16 -miu 1e-13 -num 1000 -timestepsize 0.1  -Gc 5e-6
 
 */

#include "petsc.h"
#include "VFCartFE.h"
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
	PetscInt		i,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm,ite;
	PetscReal		BBmin[3],BBmax[3];
	PetscReal		****velbc_array;
	PetscReal		***src_array;
	PetscReal		****coords_array;
	PetscReal		hx,hy,hz;
	PetscReal		lx,ly,lz;
  PetscReal		***presbc_array;
  PetscReal    perm = 1;
  PetscReal   ****perm_array;
  PetscInt    xs1,xm1,ys1,ym1,zs1,zm1;
  PetscReal  tolP = 7e-2,tolV = 7e-2;
  PetscInt    num = 10;
  PetscReal      errP=1e+10,errV=1e+10;
  Vec         Pold, Vold;
  PetscReal  pmax,vmax;
  PetscReal   InjVolrate, Q_inj;
  PetscReal  TotalLeakOff_o = 0;
  PetscReal  timevalue_o = 0;
  PetscReal  crackvolume_o = 0;
  PetscReal   vol,vol1,vol2,vol3,vol4,vol5;
  PetscReal   p = 1e-6,p_old;
  PetscInt    altminit = 0;

		
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ctx.fractureflowsolver = FRACTUREFLOWSOLVER_NONE;
	ctx.flowsolver = FLOWSOLVER_KSPMIXEDFEM;
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx.daVFperm,&xs1,&ys1,&zs1,&xm1,&ym1,&zm1);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	ierr = VecSet(ctx.VelBCArray,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.Source,1.);CHKERRQ(ierr);
	ierr = VecSet(fields.U,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.U_old,0.);CHKERRQ(ierr);
  ierr = VecDuplicate(fields.pressure,&Pold);CHKERRQ(ierr);
  ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);

	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);	
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);
  
  ierr = PetscOptionsGetReal(PETSC_NULL,"-perm",&perm,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-miu",&ctx.flowprop.mu,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-theta",&ctx.flowprop.theta,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-timestepsize",&ctx.flowprop.timestepsize,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-maxtimestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-num",&num,PETSC_NULL);CHKERRQ(ierr);
  
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
  /*
   ctx.bcP[0].face[X0] = FIXED;
   ctx.bcP[0].face[X1] = FIXED;
   ctx.bcQ[1].face[Y0] = FIXED;
   ctx.bcQ[1].face[Y1] = FIXED;
   ctx.bcQ[2].face[Z0] = FIXED;
   ctx.bcQ[2].face[Z1] = FIXED;
  */
  
   ctx.bcP[0].face[X0] = FIXED;
   ctx.bcP[0].face[X1] = FIXED;
   ctx.bcQ[1].face[Y0] = FIXED;
   ctx.bcQ[1].face[Y1] = FIXED;
   ctx.bcP[0].face[Z0] = FIXED;
   ctx.bcP[0].face[Z1] = FIXED;
   
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

//  Mechanical part
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
  
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(ctx.U_old,0.);CHKERRQ(ierr);
  
  ctx.flowprop.alphabiot = 	ctx.matprop[0].beta = 0.75;									//biot's constant
  ctx.flowprop.K_dr = ctx.matprop[0].E/(3*(1-2*ctx.matprop[0].nu));   //For 3D

  ctx.hasCrackPressure = PETSC_TRUE;
  ctx.FlowDisplCoupling = PETSC_TRUE;
  ctx.ResFlowMechCoupling = FIXEDSTRESS;
  
  ctx.FractureFlowCoupling = PETSC_TRUE;
  ctx.hasFluidSources = PETSC_FALSE;
  ctx.hasFlowWells = PETSC_TRUE;

  ierr = VecSet(ctx.PreFlowFields,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.RHSVelPpre,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.pressure_old,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.RHSPpre,0.);CHKERRQ(ierr);
  
  ierr = VolumetricFractureWellRate(&InjVolrate,&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," \n\n\n    Injected Rate: %e\n\n\n\n",InjVolrate);
  Q_inj = InjVolrate;
  
  
  ierr = PetscPrintf(PETSC_COMM_WORLD," INITIALIZATION TO COMPUTE INITIAL PRESSURE TO CREATE FRACTURE WITH INITIAL VOLUME ......\n");CHKERRQ(ierr);
  p = 1e-6;
  ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
  ctx.timestep = 0;
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
  ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
  ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);
  ierr = FieldsH5Write(&ctx,&fields);
  
  ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);

  ierr = VolumetricCrackOpening(&ctx.CrackVolume, &ctx, &fields);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Initial fracture pressure =  %e  Initial fracture volume = %e \n ",p,ctx.CrackVolume);CHKERRQ(ierr);
  ierr = FieldsH5Write(&ctx,&fields);


  ctx.timestep = 0;
  for(i = 1; i < num; i++){
    ite = 0;
    ctx.timestep++;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nPROCESSING STEP %i. \t\t Iteration = %i \t vtkfiletime = %i\n",i, ite,ctx.timestep);CHKERRQ(ierr);
    ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
    ierr = VF_StepU(&fields,&ctx);
    ierr = VolumetricLeakOffRate(&ctx.LeakOffRate,&ctx,&fields);CHKERRQ(ierr);
    ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
    ierr = FieldsH5Write(&ctx,&fields);
    ctx.flowprop.timestepsize = (ctx.CrackVolume-crackvolume_o)/(Q_inj-ctx.LeakOffRate);
    errP  = 1e+10;
    errV  = 1e+10;
    while (errV >= tolV){
      ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
    errP  = 1e+10;
    ierr = PetscPrintf(PETSC_COMM_WORLD," crack volume =  %e \n vol. leak-off rate =  %e \n errP = %e \n InjVol = %e\n",ctx.CrackVolume,ctx.LeakOffRate, errP,Q_inj);CHKERRQ(ierr);
    while (errP >= tolP){
      ierr = VecCopy(fields.pressure,Pold);CHKERRQ(ierr);
      ite++;
      ctx.timestep++;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\nTime step = %d .......Iteration step =  %i \t timevalue = %e \t timestepsize = %e \t errP = %e \t vtkfiletime = %i\n",i,ite, ctx.flowprop.timestepsize,ctx.flowprop.timestepsize, errP, ctx.timestep);CHKERRQ(ierr);
   		ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
      ierr = VF_StepU(&fields,&ctx);
      ierr = VolumetricLeakOffRate(&ctx.LeakOffRate,&ctx,&fields);CHKERRQ(ierr);
      ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
      ierr = FieldsH5Write(&ctx,&fields);
      ctx.flowprop.timestepsize = (ctx.CrackVolume-crackvolume_o)/(Q_inj-ctx.LeakOffRate);
      ierr = VecAXPY(Pold,-1.,fields.pressure);CHKERRQ(ierr);
      ierr = VecNorm(Pold,NORM_1,&errP);CHKERRQ(ierr);
      ierr = VecMax(fields.pressure,PETSC_NULL,&pmax);CHKERRQ(ierr);
      errP = errP/pmax;
      ierr = PetscPrintf(PETSC_COMM_WORLD," crack volume =  %e \n vol. leak-off rate =  %e \n errP = %e \n InjVol = %e\n",ctx.CrackVolume,ctx.LeakOffRate, errP,Q_inj);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," Pressure errP = %e\n", errP);CHKERRQ(ierr);

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
    crackvolume_o = ctx.CrackVolume;
    ierr = VFCheckVolumeBalance(&vol,&vol1,&vol2,&vol3,&vol4,&vol5,&ctx,&fields);CHKERRQ(ierr);
    
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n modulus_volume = %g\n",vol);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," divergence_volume = %g\n",vol1);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," surface_flux_volume = %g\n",vol2);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," well_volume = %g\n",vol3);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," source_volume = %g\n",vol4);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," vol.strain_volume = %g\n",vol5);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Volume Balance: RHS = %g \t LHS = %g \n",vol+vol1,vol3+vol4+vol5+Q_inj);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Reservoir Volume Balance: RHS = %g \t LHS = %g \n",vol+vol2,vol3+vol4+vol5+ctx.LeakOffRate);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Fracture Volume Balance: RHS = %g \t LHS = %g \n",ctx.LeakOffRate*ctx.flowprop.timestepsize+ctx.CrackVolume,Q_inj*ctx.flowprop.timestepsize);CHKERRQ(ierr);
    /*This will have to be called "an update function"*/
    ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSVelP,ctx.RHSVelPpre);CHKERRQ(ierr);
    ierr = VecCopy(fields.pressure,ctx.pressure_old);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSP,ctx.RHSPpre);CHKERRQ(ierr);
    ierr = VecCopy(fields.U,ctx.U_old);CHKERRQ(ierr);
    ierr = VecCopy(fields.V,ctx.V_old);CHKERRQ(ierr);
	}
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

