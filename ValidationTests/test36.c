/*
 *  test38.c: 2D/3D. Coupled reservoir and fracture flow example.
 *   (c) 2012-2014 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 *
 *    ./test38 -options_file test38.opts
 *     ./test38 -options_file test38_1.opts  -U_snes_monitor   -u_ksp_monitor -alpha 0 -beta 0 -insitumin 0,0,0,0,0,0 -insitumax 0,0,0,0,0,0 -eta 1.e-3 -unilateral nocompression
 *
 *      */

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
  PetscReal       errV=1e+10;
  Vec             U_1,P_1;
  
  
  
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  ctx.FlowDisplCoupling = PETSC_TRUE;
  ctx.ResFlowMechCoupling = FIXEDSTRESS;
  ctx.hasInsitu        = PETSC_FALSE;
  ctx.FractureFlowCoupling = PETSC_TRUE;
  ctx.hasFluidSources = PETSC_FALSE;
  ctx.hasFlowWells = PETSC_TRUE;
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
  p = 0;
  ctx.timestep = 0;
  ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
  
  ierr = VecDuplicate(fields.U,&U_1);CHKERRQ(ierr);
  ierr = VecDuplicate(fields.pressure,&P_1);CHKERRQ(ierr);
  ctx.hasInsitu        = PETSC_TRUE;

  ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
  ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
  ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);
  ierr = UpdateFractureWidth(&ctx,&fields);CHKERRQ(ierr);
  ierr = VolumetricCrackOpening(&ctx.CrackVolume, &ctx, &fields);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Initial fracture pressure =  %e  Initial fracture volume = %e \n ",p,ctx.CrackVolume);CHKERRQ(ierr);
  
  
//  
  ctx.hasInsitu        = PETSC_TRUE;
//  ierr  = VecSet(fields.pressure,0.);CHKERRQ(ierr);
//  ierr = PetscPrintf(PETSC_COMM_WORLD,"     Solving for U_s");CHKERRQ(ierr);
//  ierr = VF_StepU(&fields,&ctx);
//  ierr = VF_StepV(&fields,&ctx);
//  ierr = VF_StepU(&fields,&ctx);
//  
//  
//  ierr = UpdateFractureWidth(&ctx,&fields);CHKERRQ(ierr);
//  ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
  ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD," CrackVolume =  %e \n ",ctx.CrackVolume);CHKERRQ(ierr);

  

  ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = VecDestroy(&Vold);CHKERRQ(ierr);
	ierr = VecDestroy(&Pold);CHKERRQ(ierr);
  ierr = VecDestroy(&U_1);CHKERRQ(ierr);
  ierr = VecDestroy(&P_1);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}


