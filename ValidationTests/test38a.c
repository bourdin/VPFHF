/*
 test38.c: 2D/3D. Coupled reservoir and fracture flow example. 
 (c) 2012-2014 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu

 ./test38 -options_file test38.opts
 ./test38 -options_file test38_1.opts  -U_snes_monitor   -u_ksp_monitor -alpha 0 -beta 0 -insitumin 0,0,0,0,0,0 -insitumax 0,0,0,0,0,0 -eta 1.e-3 -unilateral nocompression 

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
	char            filename[FILENAME_MAX];
  PetscReal       errVP=1e+10;
  Vec             Pold,Vold;
  PetscReal       pmax;
  PetscReal       InjVolrate, Q_inj;
  PetscReal       crackvolume_old = 0;
  PetscReal       vol,vol1,vol2,vol3,vol4,vol5;
  PetscReal       p = 1e-6;
  PetscInt        altminit = 0;
  PetscReal       pw = 0,pw_old = 0;
  PetscReal       ini_pressure;
  PetscReal       volume;

	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  ctx.hasInsitu        = PETSC_TRUE;
  ierr = VecDuplicate(fields.pressure,&Pold);CHKERRQ(ierr);
  ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.pres",ctx.prefix);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,filename);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"Time step \t Timestepsize \t TotalTime \t Pressure \t MaxPressure \t TotVolumeInj  \t FracVolume \t SurfaceEnergy \t ElasticEnergy \t PressureWork \t TotalEnergy\n");CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&volviewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(volviewer,PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.vol",ctx.prefix);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(volviewer,filename);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(volviewer,"Time step \t Timestepsize \t TotalTime \t Pressure \t MaxPressure \t VolumeInj  \t FracVolume \t ModVolume \t StrainVolume \t SurfFlux \t VolBalance \t LeakOffVolume \t Volume4rmW \n");CHKERRQ(ierr);
  ierr = VolumetricFractureWellRate(&InjVolrate,&ctx,&fields);CHKERRQ(ierr);
  Q_inj = InjVolrate;
  ierr = PetscPrintf(PETSC_COMM_WORLD," \n\n INITIALIZATION TO COMPUTE INITIAL PRESSURE TO CREATE FRACTURE WITH INITIAL VOLUME AT INJECTION RATE = %e ......\n",InjVolrate);CHKERRQ(ierr);
  p = 1e-6;
  ctx.timestep = 0;
  ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
  ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
  ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);
  ierr = UpdateFractureWidth(&ctx,&fields);CHKERRQ(ierr);
  ierr = VolumetricCrackOpening(&ctx.CrackVolume, &ctx, &fields);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Initial fracture pressure =  %e  Initial fracture volume = %e \n ",p,ctx.CrackVolume);CHKERRQ(ierr);
  ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
  ini_pressure = 0.0;
	ierr = PetscOptionsGetReal(NULL,"-ini_pressure",&ini_pressure,NULL);CHKERRQ(ierr);
  ierr = VecSet(ctx.PresBCArray,ini_pressure);CHKERRQ(ierr);
  ierr = VecSet(ctx.pressure_old,ini_pressure);CHKERRQ(ierr);
   for (ctx.timestep = 1; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nPROCESSING STEP %i ........................................................... \n",ctx.timestep);CHKERRQ(ierr);
      altminit = 1;
      errVP  = 1e+10;
      do
      {
        altminit++;
        pw_old = pw;
        ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
        ierr = VecCopy(fields.pressure,Pold);CHKERRQ(ierr);
        ierr = VF_StepP(&fields,&ctx);
        ierr = VF_StepU(&fields,&ctx);
        ierr = VF_StepV(&fields,&ctx);
        ierr = UpdateFractureWidth(&ctx,&fields);CHKERRQ(ierr);
        ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
        ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,fields.U,&ctx);CHKERRQ(ierr);
        pw = ctx.PressureWork/ctx.CrackVolume;
        errVP = PetscAbs((pw-pw_old)/pw);
        ierr = PetscPrintf(PETSC_COMM_WORLD," Time step %i, U-P-V-loop alt min step %i, U_P_V_ERROR = %e \t pw = %e; vol(udotv) = %e\n",ctx.timestep,altminit,errVP, pw,ctx.CrackVolume);CHKERRQ(ierr);
      }
      while (errVP >= ctx.altmintol  && altminit <= ctx.altminmaxit);
    ierr = VolumetricLeakOffRate(&ctx.LeakOffRate,&ctx,&fields);CHKERRQ(ierr);
    ierr = VolumeFromWidth(&volume, &ctx, &fields);CHKERRQ(ierr);
    ierr = VFCheckVolumeBalance(&vol,&vol1,&vol2,&vol3,&vol4,&vol5,&ctx,&fields);CHKERRQ(ierr);
    ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSVelP,ctx.RHSVelPpre);CHKERRQ(ierr);
    ierr = VecCopy(fields.pressure,ctx.pressure_old);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSP,ctx.RHSPpre);CHKERRQ(ierr);
    ierr = VecCopy(fields.U,ctx.U_old);CHKERRQ(ierr);
    ierr = VecCopy(fields.V,ctx.V_old);CHKERRQ(ierr);
    ierr = VecCopy(fields.widthc,ctx.widthc_old);CHKERRQ(ierr);
    ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
    ierr = VecMax(fields.pressure,NULL,&pmax);CHKERRQ(ierr);
    ctx.ElasticEnergy = 0;
    ctx.InsituWork    = 0;
    ctx.PressureWork  = 0.;
    ctx.SurfaceEnergy = 0.;
    ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,fields.U,&ctx);CHKERRQ(ierr);
    ierr = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);
    ctx.TotalEnergy   = ctx.ElasticEnergy - ctx.InsituWork - ctx.PressureWork + ctx.SurfaceEnergy;
    ierr = PetscViewerASCIIPrintf(viewer,"%d \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e\n",ctx.timestep,ctx.timevalue,ctx.timestep*ctx.timevalue,pw,pmax,ctx.timevalue*Q_inj,ctx.CrackVolume,ctx.SurfaceEnergy,ctx.ElasticEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);
    ierr = VF_ComputeRegularizedFracturePressure(&ctx,&fields);
    ierr = PetscViewerASCIIPrintf(volviewer,"%d \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e\n",ctx.timestep,ctx.timevalue,ctx.timestep*ctx.timevalue,pw,pmax,ctx.timevalue*Q_inj,ctx.CrackVolume,vol,vol5,vol2,vol+vol5+vol2+ctx.CrackVolume-crackvolume_old,ctx.timevalue*ctx.LeakOffRate,volume);CHKERRQ(ierr);
    crackvolume_old = ctx.CrackVolume;
    ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
	}
  ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = VecDestroy(&Vold);CHKERRQ(ierr);
	ierr = VecDestroy(&Pold);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

