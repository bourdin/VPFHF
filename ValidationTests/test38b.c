/*
 *  test38b.c: 2D/3D. Coupled reservoir and fracture flow example.
 *   (c) 2012-2015 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 *
 *    ./test38b -options_file test38.opts
 *     
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
	char            filename[FILENAME_MAX];
  PetscReal       errVP=1e+10;
  Vec             Pold,Vold;
  PetscReal       pmax;
  PetscInt        altminit = 0;
  PetscReal       pw = 0,pw_old = 0;
  PetscReal       ini_pressure;
  PetscReal       errV=1e+10;



	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  ctx.FlowDisplCoupling = PETSC_TRUE;
  ctx.ResFlowMechCoupling = FIXEDSTRESS;
  ctx.hasInsitu        = PETSC_TRUE;
  ctx.FractureFlowCoupling = PETSC_FALSE;
  ctx.hasFluidSources = PETSC_FALSE;
  ctx.hasFlowWells = PETSC_TRUE;

  ctx.timestep = 0;
  ierr = VecDuplicate(fields.pressure,&Pold);CHKERRQ(ierr);
  ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
  
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
  ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
  ini_pressure = 0.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,"-ini_pressure",&ini_pressure,PETSC_NULL);CHKERRQ(ierr);
  ierr = VecSet(ctx.PresBCArray,ini_pressure);CHKERRQ(ierr);
  ierr = VecSet(ctx.pressure_old,ini_pressure);CHKERRQ(ierr);
  
   for (ctx.timestep = 1; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nPROCESSING STEP %i ........................................................... \n",ctx.timestep);CHKERRQ(ierr);
      altminit = 1;
      errVP  = 1e+10;
     do{
       ierr = UpdatePermeablitysingMultipliers(&ctx,&fields);CHKERRQ(ierr);
      do
      {
        altminit++;
        pw_old = pw;

        ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
        ierr = VecSet(fields.V,1.0);CHKERRQ(ierr);
        ierr = VecCopy(fields.pressure,Pold);CHKERRQ(ierr);
        ierr = VF_StepP(&fields,&ctx);
        ierr = VecCopy(Vold,fields.V);CHKERRQ(ierr);
        ierr = VF_StepU(&fields,&ctx);
        ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
        ierr = UpdatePermeablitysingMultipliers(&ctx,&fields);CHKERRQ(ierr);
        ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
        ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,fields.U,&ctx);CHKERRQ(ierr);
        pw = ctx.PressureWork/ctx.CrackVolume;
        errVP = PetscAbs((pw-pw_old)/pw);
        ierr = PetscPrintf(PETSC_COMM_WORLD," Time step %i, U-P-loop alt min step %i, U_P_V_ERROR = %e \t pw = %e; vol(udotv) = %e\n",ctx.timestep,altminit,errVP, pw,ctx.CrackVolume);CHKERRQ(ierr);
        }
      while (errVP >= ctx.altmintol  && altminit <= ctx.altminmaxit);
       ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
       ierr = VF_StepV(&fields,&ctx);
       ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
       ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
       ierr = PetscPrintf(PETSC_COMM_WORLD," Time step %i \t V_ERROR = %e\n",ctx.timestep,errV);CHKERRQ(ierr);
     }
     while(errV >= ctx.altmintol  && altminit <= ctx.altminmaxit);

    ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSVelP,ctx.RHSVelPpre);CHKERRQ(ierr);
    ierr = VecCopy(fields.pressure,ctx.pressure_old);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSP,ctx.RHSPpre);CHKERRQ(ierr);
    ierr = VecCopy(fields.U,ctx.U_old);CHKERRQ(ierr);
    ierr = VecCopy(fields.V,ctx.V_old);CHKERRQ(ierr);
    ierr = VecCopy(fields.widthc,ctx.widthc_old);CHKERRQ(ierr);
    ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
    ierr = VecMax(fields.pressure,PETSC_NULL,&pmax);CHKERRQ(ierr);
    ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
	}
  ierr = VecDestroy(&Vold);CHKERRQ(ierr);
	ierr = VecDestroy(&Pold);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}


