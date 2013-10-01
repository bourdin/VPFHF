/*
 test16.c: Solves for the displacement and v-field in a volume loaded line crack in 2d (Sneddon 2D)
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 test16 -options_file test16.opts is a small but relevant example
 */

#include "petsc.h"
#include "CartFE.h"
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
  PetscErrorCode ierr,ierrStepV;
  PetscViewer    viewer;
  PetscViewer    logviewer;
  PetscInt       mode=1;
  PetscInt       i,j;
  PetscInt       c;
  char           filename[FILENAME_MAX];
  PetscReal      p;
  PetscReal      p_old;
  PetscReal      prestol = 1.e-2;
  PetscInt       altminit  =1;
  Vec            Vold,U_s,U_1;
  PetscReal      vol_s,vol_1;
  PetscReal      errV=1e+10,errP;
  PetscReal      flowrate,maxvol = .03,minvol = 0.;
  PetscReal      targetVol;
  PetscBool      debug=PETSC_FALSE,StepVFailed=PETSC_FALSE;
  PetscInt       currentstep,savestep=0;

  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-mode",&mode,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-maxvol",&maxvol,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-minvol",&minvol,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-prestol",&prestol,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-debug",&debug,PETSC_NULL);CHKERRQ(ierr);
  /*
    Overwrite ctx.maxtimestep with something more reasonable
  */
  ctx.maxtimestep = 150;
  ierr            = PetscOptionsGetInt(PETSC_NULL,"-maxtimestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr);
  flowrate        = (maxvol - minvol) / (ctx.maxtimestep-1);

  /*
   Reset all BC for U and V
   */
  for (i = 0; i < 6; i++) {
    ctx.bcV[0].face[i] =NONE;
    for (j = 0; j < 3; j++) ctx.bcU[j].face[i] = NONE;
  }
  for (i = 0; i < 12; i++) {
    ctx.bcV[0].edge[i] =NONE;
    for (j = 0; j < 3; j++) ctx.bcU[j].edge[i] = NONE;
  }
  for (i = 0; i < 8; i++) {
    ctx.bcV[0].vertex[i] =NONE;
    for (j = 0; j < 3; j++) ctx.bcU[j].vertex[i] = NONE;
  }
  switch (mode) {
  case 0:
    ierr                = PetscPrintf(PETSC_COMM_WORLD,"Doing a 3D computation with null-displacement on all faces\n");CHKERRQ(ierr);
    ctx.bcU[0].face[X0] = ZERO; ctx.bcU[1].face[X0] = ZERO; ctx.bcU[2].face[X0] = ZERO; ctx.bcV[0].face[X0] = ONE;
    ctx.bcU[0].face[X1] = ZERO; ctx.bcU[1].face[X1] = ZERO; ctx.bcU[2].face[X1] = ZERO; ctx.bcV[0].face[X1] = ONE;
    ctx.bcU[0].face[Y0] = ZERO; ctx.bcU[1].face[Y0] = ZERO; ctx.bcU[2].face[Y0] = ZERO; ctx.bcV[0].face[Y0] = ONE;
    ctx.bcU[0].face[Y1] = ZERO; ctx.bcU[1].face[Y1] = ZERO; ctx.bcU[2].face[Y1] = ZERO; ctx.bcV[0].face[Y1] = ONE;
    ctx.bcU[0].face[Z0] = ZERO; ctx.bcU[1].face[Z0] = ZERO; ctx.bcU[2].face[Z0] = ZERO; ctx.bcV[0].face[Z0] = ONE;
    ctx.bcU[0].face[Z1] = ZERO; ctx.bcU[1].face[Z1] = ZERO; ctx.bcU[2].face[Z1] = ZERO; ctx.bcV[0].face[Z1] = ONE;
    break;
  case 1:
    ierr                = PetscPrintf(PETSC_COMM_WORLD,"Doing a 2D computation with null-displacement on all faces\n");CHKERRQ(ierr);
    ctx.bcU[0].face[X0] = ZERO; ctx.bcU[1].face[X0] = ZERO; ctx.bcU[2].face[X0] = ZERO; ctx.bcV[0].face[X0] = ONE;
    ctx.bcU[0].face[X1] = ZERO; ctx.bcU[1].face[X1] = ZERO; ctx.bcU[2].face[X1] = ZERO; ctx.bcV[0].face[X1] = ONE;
    ctx.bcU[1].face[Y0] = ZERO;
    ctx.bcU[1].face[Y1] = ZERO;
    ctx.bcU[0].face[Z0] = ZERO; ctx.bcU[1].face[Z0] = ZERO; ctx.bcU[2].face[Z0] = ZERO; ctx.bcV[0].face[Z0] = ONE;
    ctx.bcU[0].face[Z1] = ZERO; ctx.bcU[1].face[Z1] = ZERO; ctx.bcU[2].face[Z1] = ZERO; ctx.bcV[0].face[Z1] = ONE;
    break;
  case 2:
    ierr                      = PetscPrintf(PETSC_COMM_WORLD,"Doing a 3D computation blocking only rigid motions (may be instable)\n");CHKERRQ(ierr);
    ctx.bcU[0].vertex[X0Y0Z0] = ZERO; ctx.bcU[1].vertex[X0Y0Z0] = ZERO; ctx.bcU[2].vertex[X0Y0Z0] = ZERO;
    ctx.bcU[2].vertex[X1Y0Z0] = ZERO;
    ctx.bcU[0].vertex[X0Y1Z0] = ZERO;
    ctx.bcU[1].vertex[X0Y0Z1] = ZERO;
    ctx.bcV[0].face[X0]       = ONE;
    ctx.bcV[0].face[X1]       = ONE;
    ctx.bcV[0].face[Y0]       = ONE;
    ctx.bcV[0].face[Y1]       = ONE;
    ctx.bcV[0].face[Z0]       = ONE;
    ctx.bcV[0].face[Z1]       = ONE;
    break;
  case 3:
    ierr                      = PetscPrintf(PETSC_COMM_WORLD,"Doing a 2D computation blocking only rigid motions (may be instable)\n");CHKERRQ(ierr);
    ctx.bcU[0].vertex[X0Y0Z0] = ZERO; ctx.bcU[2].vertex[X0Y0Z0] = ZERO;
    ctx.bcU[2].vertex[X1Y0Z0] = ZERO;
    ctx.bcU[1].face[Y0]       = ZERO;
    ctx.bcU[1].face[Y1]       = ZERO;
    ctx.bcV[0].face[X0]       = ONE;
    ctx.bcV[0].face[X1]       = ONE;
    ctx.bcV[0].face[Z0]       = ONE;
    ctx.bcV[0].face[Z1]       = ONE;
    break;
  case 4:
    ierr                = PetscPrintf(PETSC_COMM_WORLD,"Doing a 3D computation blocking normal displacement on 3 planes (may lead to asymmetric solutions)\n");CHKERRQ(ierr);
    ctx.bcU[0].face[X0] = ZERO; ctx.bcU[1].face[Y0] = ZERO; ctx.bcU[2].face[Z0] = ZERO;
    ctx.bcV[0].face[X0] = ONE;
    ctx.bcV[0].face[X1] = ONE;
    ctx.bcV[0].face[Y0] = ONE;
    ctx.bcV[0].face[Y1] = ONE;
    ctx.bcV[0].face[Z0] = ONE;
    ctx.bcV[0].face[Z1] = ONE;
    break;
  case 5:
    ierr                = PetscPrintf(PETSC_COMM_WORLD,"Doing a 2D computation blocking normal displacement on 2 planes (may lead to asymmetric solutions)\n");CHKERRQ(ierr);
    ctx.bcU[0].face[X0] = ZERO;
    ctx.bcU[2].face[Z0] = ZERO;
    ctx.bcU[1].face[Y0] = ZERO;
    ctx.bcU[1].face[Y1] = ZERO;
    ctx.bcV[0].face[X0] = ONE;
    ctx.bcV[0].face[X1] = ONE;
    ctx.bcV[0].face[Z0] = ONE;
    ctx.bcV[0].face[Z1] = ONE;
    break;
  case 6:
    ierr                = PetscPrintf(PETSC_COMM_WORLD,"Doing a 3D computation blocking normal displacement on all planes but Z1\n");CHKERRQ(ierr);
    ctx.bcU[0].face[X0] = ZERO; ctx.bcV[0].face[X0] = ONE;
    ctx.bcU[0].face[X1] = ZERO; ctx.bcV[0].face[X1] = ONE;
    ctx.bcU[1].face[Y0] = ZERO; ctx.bcV[0].face[Y0] = ONE;
    ctx.bcU[1].face[Y1] = ZERO; ctx.bcV[0].face[Y1] = ONE;
    ctx.bcU[2].face[Z0] = ZERO; ctx.bcV[0].face[Z0] = ONE;
    ctx.bcV[0].face[Z1] = ONE;
    break;
  case 7:
    ierr                = PetscPrintf(PETSC_COMM_WORLD,"Doing a 2D computation blocking normal displacement on all planes but Z1\n");CHKERRQ(ierr);
    ctx.bcU[0].face[X0] = ZERO; ctx.bcV[0].face[X0] = ONE;
    ctx.bcU[0].face[X1] = ZERO; ctx.bcV[0].face[X1] = ONE;
    ctx.bcU[1].face[Y0] = ZERO;
    ctx.bcU[1].face[Y1] = ZERO;
    ctx.bcU[2].face[Z0] = ZERO; ctx.bcV[0].face[Z0] = ONE;
    ctx.bcV[0].face[Z1] = ONE;
    break;
  case 8:
    ierr                = PetscPrintf(PETSC_COMM_WORLD,"Doing a 3D computation with symmetry on x,y,z \n");CHKERRQ(ierr);
    ctx.bcU[0].face[X0] = ZERO; ctx.bcU[1].face[Y0] = ZERO; ctx.bcU[2].face[Z0] = ZERO;
    /*
      ctx.bcV[0].face[X0] = ONE;
    */
    ctx.bcV[0].face[X1] = ONE;
    /*
      ctx.bcV[0].face[Y0] = ONE;
    */
    ctx.bcV[0].face[Y1] = ONE;
    /*
      ctx.bcV[0].face[Z0] = ONE;
    */
    ctx.bcV[0].face[Z1] = ONE;
    break;
  case 9:
    ierr                = PetscPrintf(PETSC_COMM_WORLD,"Doing a 3D computation with symmetry on z \n");CHKERRQ(ierr);
    ctx.bcU[0].face[X0] = ZERO; ctx.bcU[1].face[Y0] = ZERO; ctx.bcU[2].face[Z0] = ZERO;
    ctx.bcU[0].face[X1] = ZERO; ctx.bcU[1].face[Y1] = ZERO; ctx.bcU[2].face[Z1] = ZERO;
    ctx.bcV[0].face[X0] = ONE;
    ctx.bcV[0].face[X1] = ONE;
    ctx.bcV[0].face[Y0] = ONE;
    ctx.bcV[0].face[Y1] = ONE;
    /*
      ctx.bcV[0].face[Z0] = ONE;
    */
    ctx.bcV[0].face[Z1] = ONE;
    break;
  default:
    SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: mode should be one of {0,1,2,3,4,5,6,7},got %i\n",mode);
    break;
  }

  /*
  for (c = 0; c < ctx.numPennyCracks; c++) {
    ierr = VecSet(fields.VIrrev,1.0);CHKERRQ(ierr);
    ierr = VFPennyCrackBuildVAT2(fields.VIrrev,&ctx.pennycrack[c],&ctx);CHKERRQ(ierr);
    ierr = VecPointwiseMin(fields.V,fields.VIrrev,fields.V);CHKERRQ(ierr);
  }
  for (c = 0; c < ctx.numRectangularCracks; c++) {
    ierr = VecSet(fields.VIrrev,1.0);CHKERRQ(ierr);
    ierr = VFRectangularCrackBuildVAT2(fields.VIrrev,&ctx.rectangularcrack[c],&ctx);CHKERRQ(ierr);
    ierr = VecPointwiseMin(fields.V,fields.VIrrev,fields.V);CHKERRQ(ierr);
  }
  for (c = 0; c < ctx.numWells; c++) {
    ierr = VecSet(fields.VIrrev,1.0);CHKERRQ(ierr);
    ierr = VFWellBuildVAT2(fields.VIrrev,&ctx.well[c],&ctx);CHKERRQ(ierr);
    ierr = VecPointwiseMin(fields.V,fields.VIrrev,fields.V);CHKERRQ(ierr);
  }
  ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
  */
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);

  ierr = VecCopy(fields.VIrrev,fields.V);CHKERRQ(ierr);CHKERRQ(ierr);  
  ctx.hasCrackPressure = PETSC_TRUE;
  ierr                 = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
  ierr                 = VecDuplicate(fields.U,&U_s);CHKERRQ(ierr);
  ierr                 = PetscObjectSetName((PetscObject) U_s,"U_s");CHKERRQ(ierr);
  ierr                 = VecDuplicate(fields.U,&U_1);CHKERRQ(ierr);
  ierr                 = PetscObjectSetName((PetscObject) U_1,"U_1");CHKERRQ(ierr);

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
  ierr  = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
  ierr  = VecSet(fields.U,0.0);CHKERRQ(ierr);

  ierr                = VecSet(U_s,0.0);CHKERRQ(ierr);
  ierr                = VecSet(U_1,0.0);CHKERRQ(ierr);
  ctx.matprop[0].beta = 0.;
  ctx.timevalue       = 0;

  ctx.hasCrackPressure = PETSC_FALSE;
  ctx.hasInsitu        = PETSC_TRUE;

  ierr = VecCopy(U_s,fields.U);CHKERRQ(ierr);
  ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
  ierr = VolumetricCrackOpening(&vol_s,&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  minvol %g\n",vol_s);CHKERRQ(ierr);
  
  switch (ctx.fileformat) {
  case FILEFORMAT_HDF5:
    ierr = FieldsH5Write(&ctx,&fields);CHKERRQ(ierr);
    break;
  case FILEFORMAT_BIN:
    ierr = FieldsBinaryWrite(&ctx,&fields);CHKERRQ(ierr);
    break;
  case FILEFORMAT_VTK:
    ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
    break;
  }


  for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++) {
    ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
    targetVol = minvol + flowrate * ctx.timestep;
    ierr      = PetscPrintf(PETSC_COMM_WORLD,"Time step %i. Targeting injected volume of %g\n",ctx.timestep,targetVol);CHKERRQ(ierr);
    ierr      = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
    altminit  = 0.;
    if (StepVFailed) break;
    do {
      p_old = p;
      ierr  = PetscPrintf(PETSC_COMM_WORLD,"  Time step %i, alt min step %i with pressure %g\n",ctx.timestep,altminit,p);CHKERRQ(ierr);

      switch (ctx.fileformat) {
      case FILEFORMAT_HDF5:
        ierr = FieldsH5Write(&ctx,&fields);CHKERRQ(ierr);
        break;
      case FILEFORMAT_BIN:
        ierr = FieldsBinaryWrite(&ctx,&fields);CHKERRQ(ierr);
        break;
      case FILEFORMAT_VTK:
        ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
        break;
      }

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
      ierrStepV = VF_StepV(&fields,&ctx);

      ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
      ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"      Max. change on V: %e\n",errV);CHKERRQ(ierr);

      if (debug) {
        currentstep = ctx.timestep;
        savestep++;
        ctx.timestep = savestep;
        switch (ctx.fileformat) {
        case FILEFORMAT_HDF5:
          ierr = FieldsH5Write(&ctx,&fields);CHKERRQ(ierr);
          break;
        case FILEFORMAT_BIN:
          ierr = FieldsBinaryWrite(&ctx,&fields);CHKERRQ(ierr);
          break;
        case FILEFORMAT_VTK:
          ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
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

        ctx.timestep = currentstep;
      }
      if (ierrStepV  < 0) {
        StepVFailed = PETSC_TRUE;
        break;
      }
      altminit++;
    } while ((errP >= prestol || errV >= ctx.altmintol) && altminit <= ctx.altminmaxit);
    if (!StepVFailed) {
      switch (ctx.fileformat) {
      case FILEFORMAT_HDF5:
        ierr = FieldsH5Write(&ctx,&fields);CHKERRQ(ierr);
        break;
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

