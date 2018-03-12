/*
 test16a.c: Solves for the displacement and v-field in a volume loaded line crack in 2d (Sneddon 2D)
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 test16a -options_file test16a.opts is a small but relevant example
 */

#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFMech.h"
#include "VFPermfield.h"
#include "VFCracks.h"
#include "VFWell.h"
VFCtx    ctx;
VFFields fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscViewer    viewer;
  PetscViewer    logviewer;
  char           filename[FILENAME_MAX],filenameC[FILENAME_MAX];
  PetscReal      p,po;
  PetscReal      prestol = 1.e-2;
  PetscInt       altminit  =1;
  Vec            Vold;
  PetscReal      errV=1e+10;
  PetscReal      targetVol;
  PetscBool      debug=PETSC_FALSE;
  PetscBool      saveall=PETSC_FALSE;
  int            i;
  
  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,"-prestol",&prestol,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,"-debug",&debug,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,"-saveall",&saveall,NULL);CHKERRQ(ierr);
  /*
   Overwrite ctx.maxtimestep with something more reasonable
   */
  ctx.maxtimestep = 150;
  ierr            = PetscOptionsGetInt(NULL,"-maxtimestep",&ctx.maxtimestep,NULL);CHKERRQ(ierr);
  po = 0.0001;
  ierr = PetscOptionsGetReal(NULL,"-ini_pressure",&po,NULL);CHKERRQ(ierr);  
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
  
  ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
  
  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.pres",ctx.prefix);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,filename);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"#Time step \t Volume \t Pressure \t SurfaceEnergy \t ElasticEnergy \t PressureForces \t TotalMechEnergy \n");CHKERRQ(ierr);
  
  ierr  = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr  = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr  = VecSet(fields.pressure,p);CHKERRQ(ierr);
  ierr  = VecSet(fields.U,0.0);CHKERRQ(ierr);
  
  for (i=0;i< ctx.nlayer;i++){
    ctx.matprop[i].beta  = 0.;
    ctx.matprop[i].alpha = 0.;
  }
  ctx.timevalue        = 0;
  
  for (ctx.timestep = 1; ctx.timestep < ctx.maxtimestep; ctx.timestep++) {    
    ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Time step %i. Targeting injected volume of %g\n",ctx.timestep,targetVol);CHKERRQ(ierr);
    ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
    altminit = 1;
    do {
      p = po*ctx.timestep;
      ierr  = PetscPrintf(PETSC_COMM_WORLD,"  Time step %i, alt min step %i with pressure %g\n",ctx.timestep,altminit,p);CHKERRQ(ierr);
      
      ctx.hasCrackPressure = PETSC_TRUE;
      ctx.hasInsitu        = PETSC_FALSE;
      
      ierr = PetscPrintf(PETSC_COMM_WORLD,"     Solving for U");CHKERRQ(ierr);
      ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
      ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
      
      ierr = PetscPrintf(PETSC_COMM_WORLD,"     Solving for V  ");CHKERRQ(ierr);
      ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
      ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
      ierr = VF_StepV(&fields,&ctx);
      
      ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);

      
      ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
      ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"      Max. change on V: %e\n",errV);CHKERRQ(ierr);
      
      if (debug || saveall) {
        switch (ctx.fileformat) {
          case FILEFORMAT_BIN:
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing fields in binary file\n");CHKERRQ(ierr);
            ierr = FieldsBinaryWrite(&ctx,&fields);CHKERRQ(ierr);
            break;
          case FILEFORMAT_VTK:
            ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s_nodal.%.5i-%.5i.vts",ctx.prefix,ctx.timestep,altminit);CHKERRQ(ierr);
            ierr = PetscSNPrintf(filenameC,FILENAME_MAX,"%s_cell.%.5i-%.5i.vts",ctx.prefix,ctx.timestep,altminit);CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing fields in VTK file %s %s\n",filename,filenameC);CHKERRQ(ierr);
            ierr = FieldsVTKWrite(&ctx,&fields,filename,filenameC);CHKERRQ(ierr);
            break;
        }
        ctx.ElasticEnergy = 0;
        ctx.InsituWork    = 0;
        ctx.PressureWork  = 0.;
        ctx.SurfaceEnergy = 0.;
        
        ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,fields.U,&ctx);CHKERRQ(ierr);
        ierr = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);
        
        ctx.TotalEnergy   = ctx.ElasticEnergy - ctx.InsituWork - ctx.PressureWork + ctx.SurfaceEnergy;
        
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface energy: %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer,"%d \t\t %e \t %e \t %e \t %e \t %e \t %e\n",ctx.timestep,ctx.CrackVolume,p,ctx.SurfaceEnergy,ctx.ElasticEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%i   \t%e   \t%e   \t%e   \t%e   \t%e\n",ctx.timestep,ctx.ElasticEnergy,
                                      ctx.InsituWork,ctx.SurfaceEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);
      }
      altminit++;
    } while ((errV >= ctx.altmintol) && altminit <= ctx.altminmaxit);
    switch (ctx.fileformat) {
      case FILEFORMAT_BIN:
        ierr = FieldsBinaryWrite(&ctx,&fields);CHKERRQ(ierr);
        break;
      case FILEFORMAT_VTK:
        ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
        break;
    }
    ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.log",ctx.prefix);CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&logviewer);CHKERRQ(ierr);
    ierr = PetscLogView(logviewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&logviewer);CHKERRQ(ierr);
    
    ctx.ElasticEnergy = 0;
    ctx.InsituWork    = 0;
    ctx.PressureWork  = 0.;
    ctx.SurfaceEnergy = 0.;
    ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
    ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,fields.U,&ctx);CHKERRQ(ierr);
    ierr = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);
    
    ctx.TotalEnergy   = ctx.ElasticEnergy - ctx.InsituWork - ctx.PressureWork + ctx.SurfaceEnergy;
    
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface energy: %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"%d \t\t %e \t %e \t %e \t %e \t %e \t %e\n",ctx.timestep,ctx.CrackVolume,p,ctx.SurfaceEnergy,ctx.ElasticEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%i   \t%e   \t%e   \t%e   \t%e   \t%e\n",ctx.timestep,ctx.ElasticEnergy,
                                  ctx.InsituWork,ctx.SurfaceEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);
  }
  ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  
  ierr = VecDestroy(&Vold);CHKERRQ(ierr);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");
  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}