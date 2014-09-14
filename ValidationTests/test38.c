/*
 test28.c: 1D. Flow problem with source term = 1 and Homogeneous pressure boundary conditions on all sides. Analytical solution is p = x(x-1)/2
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 

 
 
 
 ./test38 -n 100,2,100 -l 50,1,50 -m_inv 1 -E 17 -nu 0.2 -npc 1 -pc0_r 3.0 -pc0_center 25,0.5,25 -epsilon 2 -pc0_theta 0 -pc0_phi 90  -atnum 2 -pc0_thickness 0.6 -nfw 1 -fracw0_coords 25,0.5,25 -fracw0_constraint Rate -fracw0_type injector -fracw0_Qw 5e-2  -flowsolver FLOWSOLVER_snesstandarDFEM -perm 1e-16 -miu 1e-13 -num 1000 -timestepsize 0.1  -eta 1e-8 -Gc 5e-6 -V_X0_BC ONE -V_X1_BC ONE -V_Z0_BC ONE -V_Z1_BC ONE -U_Y0_BC_1 ZERO -U_Y1_BC_1 ZERO -U_X0Y0Z0_BC_0 ZERO -U_X0Y0Z0_BC_2 ZERO -U_X1Y0Z0_BC_2 ZERO

 
 
 ./test38 -n 100,2,100 -l 50,1,50 -m_inv 1 -E 17 -nu 0.2 -npc 1 -pc0_r 3.0 -pc0_center 25,0.5,25 -epsilon 1 -pc0_theta 0 -pc0_phi 90  -atnum 2 -pc0_thickness 0.6 -nfw 1 -fracw0_coords 25,0.5,25 -fracw0_constraint Rate -fracw0_type injector -fracw0_Qw 5e-2  -flowsolver FLOWSOLVER_snesstandarDFEM -perm 1e-16 -miu 1e-13 -num 2 -timestepsize 0.5  -eta 1e-8 -Gc 5e-6  -V_X0_BC ONE -V_X1_BC ONE -V_Z0_BC ONE -V_Z1_BC ONE -U_Y0_BC_1 ZERO -U_Y1_BC_1 ZERO -U_X0Y0Z0_BC_0 ZERO -U_X0Y0Z0_BC_2 ZERO -U_X1Y0Z0_BC_2 ZERO
 
  ./test38 -n 100,2,100 -l 50,1,50 -m_inv 1 -E 17 -nu 0.2 -npc 1 -pc0_r 3.0 -pc0_center 25,0.5,25 -epsilon 1 -pc0_theta 0 -pc0_phi 90  -atnum 2 -pc0_thickness 0.6 -nfw 1 -fracw0_coords 25,0.5,25 -fracw0_constraint Rate -fracw0_type injector -fracw0_Qw 5e-2  -flowsolver FLOWSOLVER_snesstandarDFEM -perm 1e-16 -miu 1e-13 -num 2 -timestepsize 0.5  -eta 1e-8 -Gc 5e-6  -V_X0_BC ONE -V_X1_BC ONE -V_Z0_BC ONE -V_Z1_BC ONE -U_Y0_BC_1 ZERO -U_Y1_BC_1 ZERO -U_X0_BC_0 ZERO -U_X0_BC_1 ZERO -U_X0_BC_2 ZERO -U_X1_BC_0 ZERO -U_X1_BC_1 ZERO -U_X1_BC_2 ZERO -U_Z0_BC_0 ZERO -U_Z0_BC_1 ZERO -U_Z0_BC_2 ZERO -U_Z1_BC_0 ZERO -U_Z1_BC_1 ZERO -U_Z1_BC_2 ZERO

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
	PetscViewer     viewer,volviewer;
	PetscViewer     logviewer;
	char            filename[FILENAME_MAX];
  PetscInt        mode=0;
	PetscInt        i,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm,ite;
	PetscReal       BBmin[3],BBmax[3];
	PetscReal       ****velbc_array;
	PetscReal       ***src_array;
	PetscReal       ****coords_array;
	PetscReal       hx,hy,hz;
	PetscReal       lx,ly,lz;
  PetscReal       ***presbc_array;
  PetscReal       perm = 1;
  PetscReal       ****perm_array;
  PetscInt        xs1,xm1,ys1,ym1,zs1,zm1;
  PetscReal       tolP = 1e-5,tolV = 1e-5,tolTime = 1e-5;
  PetscInt        num = 10;
  PetscReal       errP=1e+10,errV=1e+10,errVP=1e+10,errTime=1e+10;
  Vec             Pold,Pold1,Vold;
  PetscReal       pmax,vmax;
  PetscReal       InjVolrate, Q_inj;
  PetscReal       TotalLeakOff = 0;
  PetscReal       timevalue_o = 0;
  PetscReal       crackvolume_o = 0;
  PetscReal       vol,vol1,vol2,vol3,vol4,vol5;
  PetscReal       p = 1e-6,p_old;
  PetscInt        altminit = 0;
  PetscReal       timesize_o_ite = 0;
  PetscReal       TotalVolInj = 0;
  PetscReal       PreInt = 0;
  PetscReal       pm1 = 0,pm2 = 0, pmax1;

  
  
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
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
  ierr = VecDuplicate(fields.pressure,&Pold1);CHKERRQ(ierr);
  ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
  
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-mode",&mode,PETSC_NULL);CHKERRQ(ierr);
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
  
  switch (mode) {
    case 0:
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
      break;
    case 1:
      ctx.bcU[1].face[Y0]= ZERO;
      ctx.bcU[1].face[Y1]= ZERO;
      break;
    default:
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: mode should be one of {0,1}, got %i\n",mode);
      break;
  }
  
  
  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.pres",ctx.prefix);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,filename);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"Time step \t Timestepsize \t TotalTime \t Pressure \t VolumeInj  \t FracVolume \t LeakOffVolume \t SurfaceEnergy \t ElasticEnergy \t PressureWork \t TotalEnergy\n");CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&volviewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(volviewer,PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.vol",ctx.prefix);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(volviewer,filename);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(volviewer,"Time step \t Timestepsize \t TotalTime \t Pressure \t Pressure1 \t PressureInt \t VolumeInj  \t FracVolume \t ModVolume \t StrainVolume \t SurfFlux \t VolBalance \t LeakOffVolume \n");CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,0);CHKERRQ(ierr);
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
  ierr = VolumetricCrackOpening(&ctx.CrackVolume, &ctx, &fields);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Initial fracture pressure =  %e  Initial fracture volume = %e \n ",p,ctx.CrackVolume);CHKERRQ(ierr);
  ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
  
  ctx.timestep = 0;
  for(i = 1; i < num; i++){
    ite = 0;
    ctx.timestep = i;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nPROCESSING STEP %i. \t\t Iteration = %i \t vtkfiletime = %i\n",i, ite,ctx.timestep);CHKERRQ(ierr);
    ierr = VolumetricLeakOffRate(&ctx.LeakOffRate,&ctx,&fields);CHKERRQ(ierr);
    ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
    errV  = 1e+10;
    errVP = 1e+10;
    altminit = 0;
    while ((errV >= tolV && errVP >= tolV && altminit < 300)){
      altminit++;
      ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
      ierr = VecCopy(fields.pressure,Pold1);CHKERRQ(ierr);
      errP  = 1e+10;
      ierr = PetscPrintf(PETSC_COMM_WORLD," crack volume =  %e \n vol. leak-off rate =  %e \n errP = %e \n InjVol = %e\n",ctx.CrackVolume,ctx.LeakOffRate, errP,Q_inj);CHKERRQ(ierr);
      errTime  = 1e+10;
      errP  = 1e+10;
      while (errP >= tolP){
        ierr = VecCopy(fields.pressure,Pold);CHKERRQ(ierr);
        ite++;
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\nTime step = %d .......Iteration step =  %i \t timevalue = %e \t timestepsize = %e \t errP = %e \t vtkfiletime = %i\n",i,ite, ctx.timevalue,ctx.flowprop.timestepsize, errP, ctx.timestep);CHKERRQ(ierr);
        ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
        ierr = VF_StepU(&fields,&ctx);
        ierr = VolumetricLeakOffRate(&ctx.LeakOffRate,&ctx,&fields);CHKERRQ(ierr);
        ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
        ierr = VecAXPY(Pold,-1.,fields.pressure);CHKERRQ(ierr);
        ierr = VecNorm(Pold,NORM_INFINITY,&errP);CHKERRQ(ierr);
        ierr = VecMax(fields.pressure,PETSC_NULL,&pmax);CHKERRQ(ierr);
        errP = errP/pmax;
        ierr = PetscPrintf(PETSC_COMM_WORLD," Pressure norm = %e, max pressure = %e\n", errP,pmax);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD," crack volume =  %e \n vol. leak-off vol =  %e \n errP = %e \n InjVol = %e\n",ctx.CrackVolume,TotalLeakOff+ctx.flowprop.timestepsize*ctx.LeakOffRate, errP,TotalVolInj+ctx.flowprop.timestepsize*Q_inj);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD," Pressure errP = %e\n", errP);CHKERRQ(ierr);
        ctx.timestep = ite;
      }
      ierr = VF_StepV(&fields,&ctx);
      ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
      ierr = VecAXPY(Pold1,-1.,fields.pressure);CHKERRQ(ierr);
      ierr = VecNorm(Pold1,NORM_INFINITY,&errVP);CHKERRQ(ierr);
      ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
      ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n.........V STEP ERROR errV = %e \t errVP = %e\n\n", errV,errVP);CHKERRQ(ierr);
    }
    TotalLeakOff = TotalLeakOff + ctx.flowprop.timestepsize*ctx.LeakOffRate;
    TotalVolInj = TotalVolInj + ctx.flowprop.timestepsize*Q_inj;
    ierr = VFCheckVolumeBalance(&vol,&vol1,&vol2,&vol3,&vol4,&vol5,&ctx,&fields);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n modulus_volume = %g\n",vol);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," divergence_volume = %g\n",vol1);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," surface_flux_volume = %g\n",vol2);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," well_volume = %g\n",vol3);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," source_volume = %g\n",vol4);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," vol.strain_volume = %g\n",vol5);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Volume Balance: RHS = %g \t LHS = %g \n",vol+vol1,vol3+vol4+vol5);CHKERRQ(ierr);
    
    
    ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSVelP,ctx.RHSVelPpre);CHKERRQ(ierr);
    ierr = VecCopy(fields.pressure,ctx.pressure_old);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSP,ctx.RHSPpre);CHKERRQ(ierr);
    ierr = VecCopy(fields.U,ctx.U_old);CHKERRQ(ierr);
    ierr = VecCopy(fields.V,ctx.V_old);CHKERRQ(ierr);
    ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
    ierr = VecMax(fields.pressure,PETSC_NULL,&pmax);CHKERRQ(ierr);
    ctx.timestep = i;
    ctx.ElasticEnergy = 0;
    ctx.InsituWork    = 0;
    ctx.PressureWork  = 0.;
    ctx.SurfaceEnergy = 0.;
    ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,fields.U,&ctx);CHKERRQ(ierr);
    ierr = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);
    ctx.TotalEnergy   = ctx.ElasticEnergy - ctx.InsituWork - ctx.PressureWork + ctx.SurfaceEnergy;
    ierr = PetscViewerASCIIPrintf(viewer,"%d \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e\n",i,ctx.flowprop.timestepsize,i*ctx.flowprop.timestepsize,pmax,TotalVolInj,ctx.CrackVolume,TotalLeakOff,ctx.SurfaceEnergy,ctx.ElasticEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);
    ierr = VF_ComputeRegularizedFracturePressure(&PreInt,&ctx,&fields);
    ierr = VecMax(fields.theta,PETSC_NULL,&pmax1);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(volviewer,"%d \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e\n",i,ctx.flowprop.timestepsize,i*ctx.flowprop.timestepsize,pmax,pmax1,PreInt,ctx.flowprop.timestepsize*Q_inj,ctx.CrackVolume,vol,vol5,vol2,vol+vol5+vol2+ctx.CrackVolume-crackvolume_o,TotalLeakOff);CHKERRQ(ierr);
    crackvolume_o = ctx.CrackVolume;
    ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);

	}
  ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = VecDestroy(&Vold);CHKERRQ(ierr);
	ierr = VecDestroy(&Pold);CHKERRQ(ierr);
	ierr = VecDestroy(&Pold1);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

