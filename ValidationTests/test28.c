/*
 test28.c: 1D. Flow problem with source term = 1 and Homogeneous pressure boundary conditions on all sides. Analytical solution is p = x(x-1)/2
 (c) 2012-2014 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
 ./test28 -n 51,11,2 -l 1,1,0.01 -m_inv 10 -maxtimestep 20 -flowsolver FLOWSOLVER_snesMIXEDFEM -P_X0_BC FIXED -P_X1_BC FIXED -Q_Y0_BC_1 FIXED -Q_Y1_BC_1 FIXED -Q_Z0_BC_2 FIXED -Q_Z1_BC_2 FIXED -theta 1 -timevalue 1 -maxtimestep 20

 ./test28 -n 51,11,2 -l 1,1,0.01 -m_inv 10 -maxtimestep 20 -flowsolver FLOWSOLVER_snesstandarDFEM -P_X0_BC FIXED -P_X1_BC FIXED -Q_Y0_BC_1 FIXED -Q_Y1_BC_1 FIXED -Q_Z0_BC_2 FIXED -Q_Z1_BC_2 FIXED -theta 1 -timevalue 1 -maxtimestep 20 -kx 1. -ky 1. -kz 1 -g 0,0,0 -rhof 0. -mu 1 -theta 1.
 
 ./test28 -n 51,51,2 -l 1,1,0.01 -m_inv 10 -ts_dt 1 -ts_max_steps 20 -flowsolver FLOWSOLVER_tsMIXEDFEM -P_X0_BC FIXED -P_X1_BC FIXED -Q_Y0_BC_1 FIXED -Q_Y1_BC_1 FIXED -Q_Z0_BC_2 FIXED -Q_Z1_BC_2 FIXED -theta 1 -timevalue 1 -maxtimestep 20
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
  PetscReal   vol,vol1,vol2,vol3,vol4,vol5;
  Vec         error;
	PetscReal   norm_1,norm_2,norm_inf;

	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  ctx.FlowDisplCoupling = PETSC_FALSE;
	ctx.hasFluidSources = PETSC_TRUE;
	ctx.hasFlowWells = PETSC_FALSE;
	ierr = VecSet(ctx.Source,1.);CHKERRQ(ierr);

  if(ctx.flowsolver == FLOWSOLVER_TSMIXEDFEM || ctx.flowsolver == FLOWSOLVER_TS){
/*  	ctx.maxtimestep = 20;
    ctx.maxtimevalue = 100.;
    ctx.timevalue = 0.1;*/
    ierr = VF_StepP(&fields,&ctx);
    ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
}
  else{
	for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nProcessing step %i.\n",ctx.timestep);CHKERRQ(ierr);
		ctx.timevalue = ctx.timestep * ctx.maxtimevalue / (ctx.maxtimestep-1.);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\ntime value %f \n",ctx.timevalue);CHKERRQ(ierr);
    ierr = VF_StepP(&fields,&ctx);
    ierr = VF_FastFourierTransforms(&ctx,&fields);
    ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
    /*This will have to be called "an update function"*/
    ierr = VFCheckVolumeBalance(&vol,&vol1,&vol2,&vol3,&vol4,&vol5,&ctx,&fields);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n modulus_volume = %g\n",vol);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," divergence_volume = %g\n",vol1);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," surface_flux_volume = %g\n",vol2);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," well_volume = %g\n",vol3);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," source_volume = %g\n",vol4);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," vol.strain_volume = %g\n",vol5);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Volume Balance ::: RHS = %g \t LHS = %g \n",vol+vol1,vol3+vol4+vol5);CHKERRQ(ierr);
    ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSVelP,ctx.RHSVelPpre);CHKERRQ(ierr);
    ierr = VecCopy(fields.pressure,ctx.pressure_old);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSP,ctx.RHSPpre);CHKERRQ(ierr);
	}
	ierr = VecDuplicate(fields.VelnPress,&error);
	ierr = VecWAXPY(error,-1.0,fields.VelnPress,fields.FlowBCArray);
	ierr = VecNorm(error,NORM_1,&norm_1);
	ierr = VecNorm(error,NORM_2,&norm_2);
	ierr = VecNorm(error,NORM_INFINITY,&norm_inf);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n1_NORM = %f \n 2_norm = %f \n inf_norm = %f \n",norm_1, norm_2,norm_inf);CHKERRQ(ierr);
  }
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

