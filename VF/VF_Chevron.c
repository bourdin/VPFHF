/*
  VFracture.c
  (c) 2010 Blaise Bourdin bourdin@lsu.edu
*/

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFV.h"
#include "VFU.h"

VFCtx               ctx;
VFFields            fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode      ierr;
  
  /*
    Iteration
  */
  PetscReal           ElasticEnergy = 0;
  PetscReal           InsituWork = 0;
  PetscReal           SurfaceEnergy = 0;
  FILE                *file;
  char                H5filename[FILENAME_MAX],filename[FILENAME_MAX];

  /*
    This is essentially VIADAT
  */
  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  /* end VIADAT */

  ctx.timestep = 1 + ctx.skip;
  ierr = PetscSNPrintf(H5filename,FILENAME_MAX,"%s.%.5i.h5",ctx.prefix,ctx.timestep);CHKERRQ(ierr);
  while ( file = fopen(H5filename,"r") ) {
    fclose(file);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nProcessing step %i (%s).\n",ctx.timestep,H5filename);CHKERRQ(ierr);
    ctx.timevalue = (PetscReal) ctx.timestep;
    /*
      This is called in VTDATA
    */
    ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
    /*
      End VTDATA
    */
    
    /*
      This is essentially VPERM (minus pmult update, and transfer from GMRS)
    */
    switch (ctx.mode) {
      case ELASTICITY:
        SurfaceEnergy = 0.;
        ierr = VFElasticityTimeStep(&ctx,&fields,&ElasticEnergy,&InsituWork);CHKERRQ(ierr);
      break;
      case FRACTURE:
        ElasticEnergy = 0.;
        InsituWork = 0.;
        SurfaceEnergy = 0.;
        ierr = VFFractureTimeStep(&ctx,&fields,&ElasticEnergy,&InsituWork,&SurfaceEnergy);CHKERRQ(ierr);
      break;
    }

    /*
      This is essentially VSTOUT
    */
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Elastic Energy:            %e\n",ElasticEnergy);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Work of surface forces:    %e\n",InsituWork);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Surface energy:            %e\n",SurfaceEnergy);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Total energy:              %e\n",ElasticEnergy+SurfaceEnergy-InsituWork);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%i   \t%e   \t%e   \t%e   \t%e\n", ctx.timestep, ElasticEnergy, InsituWork, SurfaceEnergy, 
                                  ElasticEnergy - InsituWork + SurfaceEnergy);CHKERRQ(ierr); 

     switch (ctx.fileformat) {
      case EXPORT_HDF5:       
        ierr = FieldsH5Write(H5filename,&fields,&ctx);
      break;
      case EXPORT_BIN:
        ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.%.5i.h5",ctx.prefix,ctx.timestep);CHKERRQ(ierr);
        ierr = FieldsBinaryWrite(filename,&fields,&ctx);
      break; 
    } 
    ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.log",ctx.prefix);CHKERRQ(ierr);
    ierr = PetscLogPrintSummary(PETSC_COMM_WORLD, filename);CHKERRQ(ierr);

    /* 
      flushes out statistics about the current run
    */
    
    /*
      build filename of next step and loops
    */
    ctx.timestep++;
    ierr = PetscSNPrintf(H5filename,FILENAME_MAX,"%s.%.5i.h5",ctx.prefix,ctx.timestep);CHKERRQ(ierr);
  }

  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}
  