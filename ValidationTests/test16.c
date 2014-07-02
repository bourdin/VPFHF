/*
 test16.c: Solves for the displacement and v-field in a volume loaded line crack in 2d (Sneddon 2D)
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 test16 -options_file test16.opts is a small but relevant example

mpiexec -n 8 ./test16 -n 100,100,2 -l 1,1,.1 -E 1 -nu 0 -U_snes_monitor -p runtest16                  \
             -options_file ../BCfiles/2DXYRigidMotion.txt                                             \
             -V_X0_BC ONE -V_X1_BC ONE -V_Y0_BC ONE -V_Y1_BC ONE                                      \
             -insitumin -1,0,0,0,0,0  -insitumax -1,0,0,0,0,0 -npc 1                                    \
             -minvol 0 -maxvol .1 -maxtimestep 11                                                     \
             -pc0_r .2 -pc0_thickness .015 -pc0_center 0.5,0.5,0.01 -pc0_phi 90  -pc0_theta 00        \
             -U_pc_type hypre  -u_pc_hypre_type boomeramg -u_pc_hypre_boomeramg_strong_threshold 0.7  \
             -epsilon .03 -verbose 0 -eta 1e-8 -atnum 1 -unilateral none

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
  PetscReal      p;
  PetscReal      p_old;
  PetscReal      prestol = 1.e-2;
  PetscInt       altminit  =1;
  Vec            Vold,U_s,U_1;
  PetscReal      vol_s,vol_1;
  PetscReal      errV=1e+10,errP;
  PetscReal      flowrate,maxvol = .03,minvol = 0.;
  PetscReal      targetVol;
  PetscBool      debug=PETSC_FALSE;
  PetscBool      saveall=PETSC_FALSE;

  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-maxvol",&maxvol,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-minvol",&minvol,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-prestol",&prestol,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-debug",&debug,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-saveall",&saveall,PETSC_NULL);CHKERRQ(ierr);
  /*
    Overwrite ctx.maxtimestep with something more reasonable
  */
  ctx.maxtimestep = 150;
  ierr            = PetscOptionsGetInt(PETSC_NULL,"-maxtimestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr);
  flowrate        = (maxvol - minvol) / (ctx.maxtimestep-1);

  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);

  ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
  ierr = VecDuplicate(fields.U,&U_s);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) U_s,"U_s");CHKERRQ(ierr);
  ierr = VecDuplicate(fields.U,&U_1);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) U_1,"U_1");CHKERRQ(ierr);

  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.pres",ctx.prefix);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,filename);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"#Time step \t Volume \t Pressure \t SurfaceEnergy \t ElasticEnergy \t PressureForces \t TotalMechEnergy \n");CHKERRQ(ierr);

  p     = 1.e-5;
  p_old = 1.e-5;
  ierr  = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr  = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr  = VecSet(fields.pressure,p);CHKERRQ(ierr);
  ierr  = VecSet(fields.U,0.0);CHKERRQ(ierr);
  ierr  = VecSet(U_s,0.0);CHKERRQ(ierr);
  ierr  = VecSet(U_1,0.0);CHKERRQ(ierr);

  ctx.matprop[0].beta  = 0.;
  ctx.matprop[0].alpha = 0.;
  ctx.timevalue        = 0;

  for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++) {
    targetVol = minvol + flowrate * ctx.timestep;
    altminit  = 0.;

    ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Time step %i. Targeting injected volume of %g\n",ctx.timestep,targetVol);CHKERRQ(ierr);
    ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
    do {
      p_old = p;
      ierr  = PetscPrintf(PETSC_COMM_WORLD,"  Time step %i, alt min step %i with pressure %g\n",ctx.timestep,altminit,p);CHKERRQ(ierr);

      /*
        Update the pressure based on the relation
            V = p vol_1 + vol_s
        with
            vol_1 = \int U_1 \cdot \nabla V
        where U_1 is the displacement field associated with null in-situ stress and unit pressure, and
            vol_s = \int U_s \cdot \nabla V
        where U_s is the displacement field associated with null pressure and in-situ stress
      */
      ctx.hasCrackPressure = PETSC_TRUE;
      ctx.hasInsitu        = PETSC_FALSE;

      ierr = PetscPrintf(PETSC_COMM_WORLD,"     Solving for U_1");CHKERRQ(ierr);
      ierr = VecSet(fields.pressure,1.0);CHKERRQ(ierr);
      ierr = VecCopy(U_1,fields.U);CHKERRQ(ierr);
      ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
      ierr = VolumetricCrackOpening(&vol_1,&ctx,&fields);CHKERRQ(ierr);
      ierr = VecCopy(fields.U,U_1);CHKERRQ(ierr);

      ctx.hasCrackPressure = PETSC_FALSE;
      ctx.hasInsitu        = PETSC_TRUE;

      ierr = PetscPrintf(PETSC_COMM_WORLD,"     Solving for U_s");CHKERRQ(ierr);
      ierr = VecCopy(U_s,fields.U);CHKERRQ(ierr);
      ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
      ierr = VolumetricCrackOpening(&vol_s,&ctx,&fields);CHKERRQ(ierr);
      ierr = VecCopy(fields.U,U_s);CHKERRQ(ierr);

      /*
         This will fail if vol_1 = 0, which should only happen when there are no cracks
      */
      p    = (targetVol - vol_s) / vol_1;
      errP = PetscAbs((p-p_old)/p);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"     Updated crack pressure: %e\n",p);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"     Rel. change on p: %e\n",errP);CHKERRQ(ierr);
      ierr = VecAXPY(fields.U,p,U_1);CHKERRQ(ierr);

      ctx.CrackVolume      = targetVol;
      ctx.hasCrackPressure = PETSC_TRUE;

      ierr = PetscPrintf(PETSC_COMM_WORLD,"     Solving for V  ");CHKERRQ(ierr);
      ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
      ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
      ierr = VF_StepV(&fields,&ctx);

      ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
      ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"      Max. change on V: %e\n",errV);CHKERRQ(ierr);

      if (debug || saveall) {
        switch (ctx.fileformat) {
        case FILEFORMAT_BIN:
          ierr = FieldsBinaryWrite(&ctx,&fields);CHKERRQ(ierr);
          break;
        case FILEFORMAT_VTK:
          ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s_nodal.%.5i-%.5i.vts",ctx.prefix,ctx.timestep,altminit);CHKERRQ(ierr);
          ierr = PetscSNPrintf(filenameC,FILENAME_MAX,"%s_cell.%.5i-%.5i.vts",ctx.prefix,ctx.timestep,altminit);CHKERRQ(ierr);
          ierr = FieldsVTKWrite(&ctx,&fields,filename,filenameC);CHKERRQ(ierr);
          break;
        }
        ctx.ElasticEnergy = 0;
        ctx.InsituWork    = 0;
        ctx.PressureWork  = 0.;
        ctx.SurfaceEnergy = 0.;

        ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,&fields,&ctx);CHKERRQ(ierr);
        ierr = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);

        ctx.TotalEnergy   = ctx.ElasticEnergy - ctx.InsituWork - ctx.PressureWork + ctx.SurfaceEnergy;

        ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface energy: %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer,"%d \t\t %e \t %e \t %e \t %e \t %e \t %e\n",ctx.timestep,ctx.CrackVolume,p,ctx.SurfaceEnergy,ctx.ElasticEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%i   \t%e   \t%e   \t%e   \t%e   \t%e\n",ctx.timestep,ctx.ElasticEnergy,
                                      ctx.InsituWork,ctx.SurfaceEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);
      }
      altminit++;
    } while ((errP >= prestol || errV >= ctx.altmintol) && altminit <= ctx.altminmaxit);
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
    ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,&fields,&ctx);CHKERRQ(ierr);
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
  ierr = VecDestroy(&U_s);CHKERRQ(ierr);
  ierr = VecDestroy(&U_1);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");
  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}

