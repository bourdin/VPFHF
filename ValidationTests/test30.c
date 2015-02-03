/*
 test30.c: 1D Terzaghi problem from Zheng et al. Coupling of flow and displacmement solver for geomechanics application.
 
 Reservoir simultion with the finite element method using Biot poroelastic approach
 Zheng, Y., Burridge R. and Burns D.
 
 (c) 2012-2014 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
./test30  -options_file test30.opts
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
  Vec         V_hold;
  PetscReal   vol,vol1,vol2,vol3,vol4,vol5;
  PetscReal   displ_p_tol = 1e-6;
  Vec         error;
  Vec         PreIteSol;
	PetscReal   norm_inf = 1e+3;
  PetscInt    displ_iter = 0;
  
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx.daScal,&V_hold);CHKERRQ(ierr);
  ierr = VecDuplicate(fields.pressure,&error);
  ierr = VecDuplicate(fields.pressure,&PreIteSol);
  ctx.hasFlowWells = PETSC_FALSE;
	ctx.hasFluidSources = PETSC_FALSE;
  ctx.hasInsitu        = PETSC_TRUE;
  ctx.FlowDisplCoupling = PETSC_TRUE;
	ctx.hasCrackPressure = PETSC_FALSE;
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
  ctx.timestep = 0;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  Computing initial time step solution\n");CHKERRQ(ierr);
  
  ierr = VecSet(fields.pressure,0.);CHKERRQ(ierr);
  ierr = VecSet(fields.theta,0.);CHKERRQ(ierr);
  ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
  

  ierr = VF_StepP(&fields,&ctx);
  
  ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);

  ierr = VecCopy(fields.pressure,PreIteSol);CHKERRQ(ierr);
  while (norm_inf > displ_p_tol) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"   Iteration Step: %d\n",displ_iter);CHKERRQ(ierr);
    ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
    ierr = VF_StepP(&fields,&ctx);
    displ_iter++;
    ierr = VecWAXPY(error,-1.0,fields.pressure,PreIteSol);
    ierr = VecCopy(fields.pressure,PreIteSol);CHKERRQ(ierr);
    ierr = VecNorm(error,NORM_INFINITY,&norm_inf);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n inf_norm = %f \n",norm_inf);CHKERRQ(ierr);
    ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);

  }
  ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
  ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
  ierr = VecCopy(ctx.RHSVelP,ctx.RHSVelPpre);CHKERRQ(ierr);
  ierr = VecCopy(fields.pressure,ctx.pressure_old);CHKERRQ(ierr);
  ierr = VecCopy(ctx.RHSP,ctx.RHSPpre);CHKERRQ(ierr);
  
  
  
  
  ctx.bcQ[2].face[Z0] = NONE;
  ctx.bcP[0].face[Z1] = FIXED;
  ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
  ierr = VecSet(ctx.RHSVelPpre,0.);CHKERRQ(ierr);
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
      ierr = VF_StepP(&fields,&ctx);
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
    ierr = VF_FastFourierTransforms(&ctx,&fields);
    ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
    ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSVelP,ctx.RHSVelPpre);CHKERRQ(ierr);
    ierr = VecCopy(fields.pressure,ctx.pressure_old);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSP,ctx.RHSPpre);CHKERRQ(ierr);
  }
  
  
  
  
  ierr = VecDestroy(&PreIteSol);CHKERRQ(ierr);
  ierr = VecDestroy(&error);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}