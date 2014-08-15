/*
 test1.c:
 1D tension experiment
 
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */

#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFMech.h"
#include "VFFlow.h"
#include "VFCracks.h"
#include "VFPermfield.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  VFCtx          ctx;
  VFFields       fields;
  PetscErrorCode ierr;
  Vec            Vold;
  PetscReal      errV;
  PetscInt       altminit;
  PetscInt       nopt=3;
  PetscBool      flg;  
  PetscReal      p = 0.,crackVolume=0;

  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  
  PetscInt  orientation = 2;
  ierr = PetscOptionsGetReal(NULL,"-pressure",&p,NULL);CHKERRQ(ierr);
  ctx.hasCrackPressure = PETSC_TRUE;
  ctx.matprop[0].beta  = 0.;
  ctx.matprop[0].alpha = 0.;
  
  ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
  ierr = VecCopy(fields.VIrrev,fields.V);CHKERRQ(ierr);
  ierr = VecSet(fields.U,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.BCU,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
  ierr = VecSet(fields.BCU,0.0);CHKERRQ(ierr);

  for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
    if (ctx.maxtimestep > 1) {
      ctx.timevalue = ctx.mintimevalue + (ctx.maxtimevalue - ctx.mintimevalue) / (PetscReal) (ctx.maxtimestep-1) * (PetscReal) (ctx.timestep);
    } else {
      ctx.timevalue = ctx.maxtimevalue;
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"==== time step %i: t = %e\n",ctx.timestep,ctx.timevalue);CHKERRQ(ierr);
    ierr = VecSetFromBC(fields.BCU,&ctx.bcU[0]);CHKERRQ(ierr);
    ierr = VecScale(fields.BCU,ctx.timevalue);CHKERRQ(ierr);
    ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
    
    altminit = 0;
    errV = 1.e+10;
    do {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Step %i / %i\n",ctx.timestep,altminit);CHKERRQ(ierr);
      ierr = VecCopy(fields.V,Vold);
      ierr = VF_StepU(&fields,&ctx);
      ierr = VF_StepV(&fields,&ctx);

      ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
      ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"      Max. change on V: %e\n",errV);CHKERRQ(ierr);
      
      if (altminit%10 == 0) {
        ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);
      }
      altminit++;
    } while (errV >= ctx.altmintol && altminit <= ctx.altminmaxit);
    
    ierr = VolumetricCrackOpening(&crackVolume,&ctx,&fields);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Total crack opening: %f\n",crackVolume);CHKERRQ(ierr);
    ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
    
    ctx.SurfaceEnergy = 0.;
    ctx.ElasticEnergy = 0;
    ctx.InsituWork    = 0;
    ctx.PressureWork  = 0.;
    
    ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,fields.U,&ctx);CHKERRQ(ierr);
    ierr = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);
    
    ctx.TotalEnergy = ctx.ElasticEnergy-ctx.InsituWork-ctx.PressureWork;
    
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic Energy:            %e\n",ctx.ElasticEnergy);CHKERRQ(ierr);
    if (ctx.hasCrackPressure) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of pressure forces:   %e\n",ctx.PressureWork);CHKERRQ(ierr);
    }
    if (ctx.hasInsitu) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of insitu stresses:   %e\n",ctx.InsituWork);CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Total Mechanical energy:  %e\n",ctx.ElasticEnergy-ctx.InsituWork-ctx.PressureWork);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface  energy:          %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);
    ctx.TotalEnergy = ctx.ElasticEnergy-ctx.InsituWork-ctx.PressureWork + ctx.SurfaceEnergy;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Total  energy:            %e\n",ctx.TotalEnergy);CHKERRQ(ierr);
		ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%d \t\t%e \t%e \t%e \t%e \t%e \t%e\n",ctx.timestep,ctx.timevalue,ctx.SurfaceEnergy,ctx.ElasticEnergy,ctx.PressureWork,ctx.InsituWork,ctx.TotalEnergy);
    
    ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);
  }
  ierr = VecDestroy(&Vold);CHKERRQ(ierr);
  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}

