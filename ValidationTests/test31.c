/*
 test31.c: 2D Mandel's problem. Coupling of flow and displacmement solver for geomechanics application.
 See http://scholar.lib.vt.edu/theses/available/etd-07012008-115136/unrestricted/Dissertation_ImsooLee.pdf
  
 ./test31 -options_file test31.opts
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
	char			filename[FILENAME_MAX];
	PetscReal		BBmin[3],BBmax[3];
	PetscReal		lx,ly,lz;
  PetscReal   vol,vol1,vol2,vol3,vol4,vol5;
  PetscReal   pressure;
  PetscReal   stress;
  PetscReal   displ;
  PetscReal   SumnIntegral = 0;
  PetscReal   displ_p_tol = 1e-6;
  Vec         error;
  Vec         PreIteSol;
	PetscReal   norm_inf = 1e+3;
  PetscInt    displ_iter = 0;

	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  ctx.hasFlowWells = PETSC_FALSE;
	ctx.hasFluidSources = PETSC_FALSE;
  ctx.hasInsitu        = PETSC_FALSE;
  ctx.FlowDisplCoupling = PETSC_TRUE;
	ctx.hasCrackPressure = PETSC_FALSE;
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
  lz = BBmax[2]-BBmin[2];
	ly = BBmax[1]-BBmin[1];
	lx = BBmax[0]-BBmin[0];
  stress = -1.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,"-stress",&stress,PETSC_NULL);CHKERRQ(ierr);
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%smandel.txt",ctx.prefix);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,filename);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"#Time \t Displs \t Pressure \n");CHKERRQ(ierr);
  ierr = VecDuplicate(fields.pressure,&error);
  ierr = VecDuplicate(fields.pressure,&PreIteSol);
  ctx.timestep = 0;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  Computing initial time step solution..................\n");CHKERRQ(ierr);
  displ = lz*(1.+ctx.matprop[0].nu)*((1.-ctx.matprop[0].nu)*stress*lx*ly+ctx.matprop[0].beta*SumnIntegral*(1.-2*ctx.matprop[0].nu))/(lx*ly*ctx.matprop[0].E);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"   Initial guess displacement: %g\n",displ);CHKERRQ(ierr);
  ierr = VecSet(fields.BCU,displ);CHKERRQ(ierr);  
  ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
  ierr = VF_StepP(&fields,&ctx);
  ierr = VF_IntegrateOnBoundary(&SumnIntegral,fields.pressure,Z1,&ctx);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n Integrand = %f \n",SumnIntegral);CHKERRQ(ierr);
  ierr = VecCopy(fields.pressure,PreIteSol);CHKERRQ(ierr);
  while (norm_inf > displ_p_tol) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"   Iteration Step: %d\n",displ_iter);CHKERRQ(ierr);
    displ = lz*(1.+ctx.matprop[0].nu)*((1.-ctx.matprop[0].nu)*stress*lx*ly+ctx.matprop[0].beta*SumnIntegral*(1.-2*ctx.matprop[0].nu))/(lx*ly*ctx.matprop[0].E);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"   Update displacement: %g \t Displ_error = %g\n",displ,norm_inf);CHKERRQ(ierr);
    ierr = VecSet(fields.BCU,displ);CHKERRQ(ierr);
    ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
    ierr = VF_StepP(&fields,&ctx);
    ierr = VF_IntegrateOnBoundary(&SumnIntegral,fields.pressure,Z1,&ctx);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n Integrand = %f \n",SumnIntegral);CHKERRQ(ierr);
    displ_iter++;
    ierr = VecWAXPY(error,-1.0,fields.pressure,PreIteSol);
    ierr = VecCopy(fields.pressure,PreIteSol);CHKERRQ(ierr);
    ierr = VecNorm(error,NORM_INFINITY,&norm_inf);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n inf_norm = %f \n",norm_inf);CHKERRQ(ierr);
  }
  ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
  ierr = VecMax(fields.pressure,PETSC_NULL,&pressure);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"%e \t %e \t %e \n",ctx.timestep*ctx.timevalue,displ,pressure);CHKERRQ(ierr);
  ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
  ierr = VecCopy(ctx.RHSVelP,ctx.RHSVelPpre);CHKERRQ(ierr);
  ierr = VecCopy(fields.pressure,ctx.pressure_old);CHKERRQ(ierr);
  ierr = VecCopy(ctx.RHSP,ctx.RHSPpre);CHKERRQ(ierr);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  COMPUTING TIME EVOLUTION OF SUBSEQUENT DEFORMATION STAGES..................\n");CHKERRQ(ierr);
  ctx.bcQ[0].face[X0] = NONE;
  ctx.bcQ[0].face[X1] = NONE;
  ctx.bcP[0].face[X0] = FIXED;
  ctx.bcP[0].face[X1] = FIXED;
  ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
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
      ierr = VF_IntegrateOnBoundary(&SumnIntegral,fields.pressure,Z1,&ctx);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n Integrand = %f \n",SumnIntegral);CHKERRQ(ierr);
      displ = lz*(1.+ctx.matprop[0].nu)*((1.-ctx.matprop[0].nu)*stress*lx*ly+ctx.matprop[0].beta*SumnIntegral*(1.-2*ctx.matprop[0].nu))/(lx*ly*ctx.matprop[0].E);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"   Initial guess displacement: %g\n",displ);CHKERRQ(ierr);
      ierr = VecSet(fields.BCU,displ);CHKERRQ(ierr);
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
    
    ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
    ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSVelP,ctx.RHSVelPpre);CHKERRQ(ierr);
    ierr = VecCopy(fields.pressure,ctx.pressure_old);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSP,ctx.RHSPpre);CHKERRQ(ierr);
    if(ctx.timestep % 50 == 0){
      ierr = VecMax(fields.pressure,PETSC_NULL,&pressure);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"%e \t %e \t %e \n",ctx.timestep*ctx.timevalue,displ,pressure);CHKERRQ(ierr);
    }
  }
  ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = VecDestroy(&PreIteSol);CHKERRQ(ierr);
  ierr = VecDestroy(&error);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

