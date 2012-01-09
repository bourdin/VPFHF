/*
  VFracture.c
  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFV.h"
#include "VFU.h"
#include "VFFlow.h"
/*
  Remove this as soon as the functions are moved in the proper location
*/
#include "VFFlow_MixedFEM.h"

VFCtx               ctx;
VFFields            fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  VFCtx               ctx;
  VFFields            fields;
  PetscErrorCode      ierr;
  char                filename[FILENAME_MAX];

  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  ierr = FlowSolverInitialize(&ctx,&fields);CHKERRQ(ierr);

  for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nProcessing step %i.\n",ctx.timestep);CHKERRQ(ierr);
    ctx.timevalue = ctx.timestep * ctx.maxtimevalue / (ctx.maxtimestep-1.);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\ntime value %f \n",ctx.timevalue);CHKERRQ(ierr);

    /*
      Do flow solver step 
    */
    ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);

    /*
      Do mechanics step
    */
    ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
    switch (ctx.mechsolver) {
      case ELASTICITY:
        ierr = VFElasticityTimeStep(&ctx,&fields);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic Energy:            %e\n",ctx.ElasticEnergy);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of surface forces:    %e\n",ctx.InsituWork);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Total energy:              %e\n",ctx.TotalEnergy);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%i   \t%e   \t%e   \t%e\n",ctx.timestep,ctx.ElasticEnergy,ctx.TotalEnergy);CHKERRQ(ierr); 
      break;
      case FRACTURE:
        ierr = VFFractureTimeStep(&ctx,&fields);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic Energy:            %e\n",ctx.ElasticEnergy);CHKERRQ(ierr);
        if (ctx.hasInsitu) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of surface forces:    %e\n",ctx.InsituWork);CHKERRQ(ierr);
        }
        if (ctx.hasCrackPressure) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of pressure forces:   %e\n",ctx.PressureWork);CHKERRQ(ierr);
        }
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface energy:            %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Total energy:              %e\n",ctx.TotalEnergy);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%i   \t%e   \t%e   \t%e   \t%e   \t%e\n",ctx.timestep,ctx.ElasticEnergy,
          ctx.InsituWork,ctx.SurfaceEnergy,ctx.TotalEnergy);CHKERRQ(ierr); 
      break;
      case NOMECH:
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Skipping mechanics step\n");CHKERRQ(ierr);
      break;
    }


    /*
      Save fields and write statistics about current run
    */    
    switch (ctx.fileformat) {
      case FILEFORMAT_HDF5:       
        ierr = FieldsH5Write(&ctx,&fields);
      break;
      case FILEFORMAT_BIN:
        ierr = FieldsBinaryWrite(&ctx,&fields);
      break; 
    } 
    ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.log",ctx.prefix);CHKERRQ(ierr);
    ierr = PetscLogPrintSummary(PETSC_COMM_WORLD,filename);CHKERRQ(ierr);

  }

	ierr = FlowSolverFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
  return(0);
}
  
