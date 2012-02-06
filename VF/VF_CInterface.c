#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFCracks.h"
#include "VFWell.h"
#include "VFU.h"
#include "../Utils/xdmf.h"
#include "../Utils/VFGMRS.h"

VFCtx               ctx;
VFFields            fields;

#undef __FUNCT__
#define __FUNCT__ "VFInitialize"
/*
  VFInitialize: Initialize the VF code. Called by the fortran implementation of VIADAT

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFInitialize(PetscInt nx,PetscInt ny,PetscInt nz,PetscReal *dx,PetscReal *dy,PetscReal *dz)
{
  PetscErrorCode ierr;

  PetscTruth      printhelp;
  FILE           *file;
  char            filename[FILENAME_MAX];
  Vec             V;
  PetscInt        c;
  VFPennyCrack   *crack;
  VFWell         *well;
  char            prefix[PETSC_MAX_PATH_LEN+1];
  
  PetscFunctionBegin;
  ierr = PetscPrintf(PETSC_COMM_WORLD,banner);CHKERRQ(ierr);
#if defined(PETSC_USE_DEBUG)
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      ##########################################################\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      #                                                        #\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      #                          WARNING!!!                    #\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      #                                                        #\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      #   This code was compiled with a debugging option,      #\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      #   For production runs, use a petsc compiled with       #\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      #   optimization, the performance will be generally      #\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      #   two or three times faster.                           #\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      #                                                        #\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      ##########################################################\n\n\n");CHKERRQ(ierr);
#endif
  

  ierr = PetscOptionsHasName(PETSC_NULL,"-help",&printhelp);CHKERRQ(ierr);

  ierr = VFLogInitialize(&ctx.vflog);CHKERRQ(ierr);

  ctx.ncellx = nx-1;
  ctx.ncelly = ny-1;
  ctx.ncellz = nz-1;
  ierr = VFCtxGet(&ctx);CHKERRQ(ierr);
  ierr = VFPropGet(&ctx.vfprop);CHKERRQ(ierr);

  ierr = PetscMalloc(ctx.nlayer * sizeof(MatProp),&ctx.matprop);CHKERRQ(ierr);
  ierr = VFMatPropGet(ctx.matprop,ctx.nlayer);CHKERRQ(ierr);

  if (printhelp) {
    ierr = PetscFinalize();
    return(-1);
  }

  ierr = VFGeometryInitialize(&ctx,dx,dy,dz);CHKERRQ(ierr);
  ierr = VFFieldsInitialize(&ctx,&fields);CHKERRQ(ierr);
  ierr = VFBCInitialize(&ctx);CHKERRQ(ierr);
  ierr = VFSolversInitialize(&ctx);CHKERRQ(ierr);

  /*
    Save command line options to a file
  */
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.txt",ctx.prefix,0);CHKERRQ(ierr);
  file = fopen(filename,"w");
  ierr = PetscOptionsPrint(file);CHKERRQ(ierr);
  fclose(file);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Option table:\n");CHKERRQ(ierr);
  ierr = PetscOptionsPrint(stdout);CHKERRQ(ierr);

  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.ener",ctx.prefix);CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&ctx.energyviewer);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"#i,Elastic Energy,InsituWork,Surface Energy,Total Energy\n");CHKERRQ(ierr);
  ierr = PetscViewerFlush(ctx.energyviewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "vendit"
/*
  vendit:

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode vendit(PetscInt rank,PetscReal tim,PetscInt nstep,PetscInt nfout,PetscInt nfbug)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "vtdata"
/*
  vtdata: 

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode vtdata(PetscInt rank,PetscReal tim,PetscInt nstep,PetscInt nfout,PetscInt nfbug)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "vperm"
/*
  vperm: This is where the permeability is updated from the pressure + temperature, i.e. the main hook for the fracture computation

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode vperm(PetscInt rank,PetscReal *pres,PetscReal *tempr,PetscReal *pmult,PetscReal tim,PetscInt nstep,PetscInt newt,PetscInt nfout,PetscInt nfbug,PetscInt n)
{
  PetscReal      fieldmin,fieldmax;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nThis is the PETSc / C implementation of %s\n",__FUNCT__);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"   TIME  = %15.5E\n",tim);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"   NSTEP = %10i\n",nstep);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"   NEWT  = %10i\n",newt);CHKERRQ(ierr);
  
  ctx.timestep  = nstep;
  ctx.timevalue = tim;
  
  /*
    Compute boundary displacements if necessary
  */
  if (nstep == 1 && newt == 1 && ctx.hasInsitu) {
    ierr = VF_ComputeBCU(&fields,&ctx);CHKERRQ(ierr);
  }
  
  if (ctx.coupling != COUPLING_NONE) {
    /*
      Get Pressure and temperature from GMRS
    */
    ierr = VecGMRSToPetsc(pres,fields.pressure);CHKERRQ(ierr);  
    ierr = VecGMRSToPetsc(tempr,fields.theta);CHKERRQ(ierr);  
    ierr = VecGMRSToPetsc(pmult,fields.pmult);CHKERRQ(ierr);
    
    /*
	  Scaling pressure from GMRS - may need adjustment
	*/
    ierr = VecScale(fields.pressure, 1.e-6);CHKERRQ(ierr);	
	
    if (ctx.verbose > 0) {
      ierr = VecMin(fields.theta,PETSC_NULL,&fieldmin);CHKERRQ(ierr);
      ierr = VecMax(fields.theta,PETSC_NULL,&fieldmax);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Temperature range: %g - %g\n",fieldmin,fieldmax);CHKERRQ(ierr);
      ierr = VecMin(fields.pressure,PETSC_NULL,&fieldmin);CHKERRQ(ierr);
      ierr = VecMax(fields.pressure,PETSC_NULL,&fieldmax);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Pressure range:    %g - %g\n",fieldmin,fieldmax);CHKERRQ(ierr);
      ierr = VecMin(fields.pmult,PETSC_NULL,&fieldmin);CHKERRQ(ierr);
      ierr = VecMax(fields.pmult,PETSC_NULL,&fieldmax);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Pmult range:       %g - %g\n",fieldmin,fieldmax);CHKERRQ(ierr);
    }
    /*
      If this is the first Newton iteration of the first time step, set the reference temperature:
    */
    if (nstep == 1 && newt == 1) {
      ierr = VecCopy(fields.theta,fields.thetaRef);CHKERRQ(ierr);
      ierr = VecCopy(fields.pressure,fields.pressureRef);CHKERRQ(ierr);
    }
    /* 
      Minimize the fracture energy in order to get Displacement and Fracture
    */
    /*
      Compute boundary displacements if necessary
    */
    if (nstep == 1 && newt == 1 && ctx.hasInsitu) {
      ierr = VF_ComputeBCU(&fields,&ctx);CHKERRQ(ierr);
    }
    if (newt == 1) {
      /*
        We are starting a new time step
      */
      ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
    }

    switch (ctx.mode) {
      case ELASTICITY:
        ierr = VFElasticityTimeStep(&ctx,&fields);CHKERRQ(ierr);
      break;
      case FRACTURE:
        ierr = VFFractureTimeStep(&ctx,&fields);CHKERRQ(ierr);
      break;
    }
    
    if (ctx.coupling == COUPLING_FULL) {
      /*
        Send permeability multiplier back to GMRS
      */
      if (ctx.verbose > 0) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Updating permeability multipliers and sending them back to GMRS\n");CHKERRQ(ierr)
      }
      ierr = PermUpdate(fields.V,fields.pmult,&ctx.vfprop,&ctx);
      ierr = VecPetscToGMRS(fields.pmult,pmult);CHKERRQ(ierr);
      //ierr = VecPetscToGMRS(fields.V,pmult);CHKERRQ(ierr);
    }
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic Energy:            %e\n",ctx.ElasticEnergy);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of insitu forces:     %e\n",ctx.InsituWork);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface energy:            %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Total energy:              %e\n",ctx.TotalEnergy);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%i \t%e \t%e \t%e \t%e \t%e\n",nstep,tim,ctx.ElasticEnergy,
                                ctx.InsituWork,ctx.SurfaceEnergy,ctx.TotalEnergy);CHKERRQ(ierr); 
  ierr = PetscViewerFlush(ctx.energyviewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "vstdout"
/*
  vstdout: 

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode vstdout(PetscInt rank,PetscReal tim,PetscInt nstep,PetscInt nfout,PetscInt nfbug)
{
  PetscErrorCode ierr;
  char           filename[FILENAME_MAX];
   
  PetscFunctionBegin;
  /*
    Since we already save the energies at the end of each time step, there is no need to do it again
  */
  /*
  ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%i \t%e \t%e \t%e \t%e \t%e\n",nstep,tim,ctx.ElasticEnergy,ctx.InsituWork,ctx.SurfaceEnergy,
                                ctx.TotalEnergy);CHKERRQ(ierr); 
  ierr = PetscViewerFlush(ctx.energyviewer);CHKERRQ(ierr);
  */
  ierr = PetscLogStagePush(ctx.vflog.VF_IOStage);CHKERRQ(ierr);
  switch (ctx.fileformat) {
    case FILEFORMAT_HDF5:       
      ierr = FieldsH5Write(&ctx,&fields);
    break;
    case FILEFORMAT_BIN:
      ierr = FieldsBinaryWrite(&ctx,&fields);
    break; 
  } 
  ierr = PetscLogStagePop();CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.log",ctx.prefix);CHKERRQ(ierr);
  ierr = PetscLogPrintSummary(PETSC_COMM_WORLD,filename);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
