/*
 test16.c: Solves for the displacement and v-field in a volume loaded line crack in 2d (Sneddon 2D)
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu

mpiexec -n 2 ../test16 -n 26,2,50 -l .5,.02,1 -epsilon .04 -mode 1 eta 1e-8 \
        -Gc 1e-1 -nu 0 -maxtimestep 2 -maxvol .001 -mechsolver ELASTICITY   \
        -nc 1 -c0_r .1 -c0_center .25,.01,.5 

Try:  
mpiexec -n 2 ./test16  -n 26,51,2 -l 0.5,1,0.02 -epsilon 0.02 -length 2. \
      -center 0,0.5,0 -mode 2 -nu 0. -eta 1e-8 -maxtimestep 10 -maxvol .03 \
      -Gc 1.e-1  -insitumin 0.,-1.e-4,0.,0.,0.,0. -insitumax 0.,-1.e-4,0.,0.,0.,0.

mpiexec -n 2 ./test16  -n 26,51,2 -l 0.5,1,0.02 -epsilon 0.02 -length .1 \
      -center 0,0.5,0 -mode 2 -nu 0. -eta 1e-8 -maxtimestep 10 -maxvol .03 
      -Gc 1.e+1  -insitumin 0.,-1.e-4,0.,0.,0.,0. -insitumax 0,-1e-4,0.,0.,0.,0.

mpiexec -n 2 ./test16  -n 26,51,2 -l 0.5,1,0.02 -epsilon 0.02 -length .1 \
      -center 0,0.5,0 -mode 2 -nu 0. -eta 1e-8 -maxtimestep 10 -maxvol .03 \
      -Gc 1.e+1  -insitumin 0.,-2.e+0,0.,0.,0.,0. -insitumax 0,-2e+0,0.,0.,0.,0. 
      
An example that doesn't work:
mpiexec -n 4 ./test16  -n 26,51,2 -l 0.5,1,0.02 -epsilon 0.02 -length .1 \
      -center 0,0.5,0 -mode 2 -nu 0. -eta 1e-8 -maxtimestep 21 -maxvol .03 \
      -Gc 1.e+1  -insitumin 0.,-4.e+0,0.,0.,0.,0. -insitumax 0,-4e+0,0.,0.,0.,0. \
      -U_ksp_atol 1e-7 -U_ksp_rtol 1e-7 -U_ksp_type cg 
      
Unsurprisingly, setting eta to 1e-10 helps a lot
Using a finer mesh solves the issue:
mpiexec -n 4 ./test16  -n 52,101,2 -l 0.5,1,0.02 -epsilon 0.01 -length .1 \
        -center 0,0.5,0 -mode 2 -nu 0. -eta 1e-8 -maxtimestep 21 -maxvol .03 \
        -Gc 1.e+1  -insitumin 0.,-4.e+0,0.,0.,0.,0. -insitumax 0,-4e+0,0.,0.,0.,0. \
        -U_ksp_atol 1e-7 -U_ksp_rtol 1e-7 -U_ksp_type cg -eta 1e-10

Same using gamg:
mpiexec -n 4 ../test16  -n 52,101,2 -l 0.5,1,0.02 -epsilon 0.01 -length .1 \
        -center 0,0.5,0 -mode 2 -nu 0. -eta 1e-8 -maxtimestep 21 -maxvol .03 \
        -Gc 1.e+1  -insitumin 0.,-4.e+0,0.,0.,0.,0. -insitumax 0,-4e+0,0.,0.,0.,0. \
        -U_ksp_atol 1e-7 -U_ksp_rtol 1e-7 -U_ksp_type cg -eta 1e-10 \
        -U_ksp_monitor_true_residual  -U_mg_levels_ksp_chebyshev_estimate_eigenvalues 0,0.1,0,1.1 \
        -U_mg_levels_ksp_type chebyshev -U_mg_levels_pc_type sor -U_pc_gamg_agg_nsmooths 1 \
        -U_pc_gamg_threshold 0 -U_pc_gamg_type agg -U_pc_type gamg        

An even thiner mesh is required for higher in-situ stresses
mpiexec -n 6 ../test16  -n 101,201,2 -l 0.5,1,0.02 -epsilon 0.005 -length .1 \
        -center 0,0.5,0 -mode 2 -nu 0. -eta 1e-8 -maxtimestep 21 -maxvol .03 \
        -Gc 1.e+1  -insitumin 0.,-8.e+0,0.,0.,0.,0. -insitumax 0,-8e+0,0.,0.,0.,0. \
        -U_ksp_atol 1e-7 -U_ksp_rtol 1e-7 -eta 1e-10
 
Going back to default KSP is OK
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFV.h"
#include "VFU.h"
#include "VFPermfield.h"
#include "VFCracks.h"

VFCtx               ctx;
VFFields            fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode  ierr;
  PetscViewer     viewer;
  PetscViewer     logviewer;
  PetscInt        mode=1;
  PetscInt        i,j;
  PetscInt        c;
  char            filename[FILENAME_MAX];
  PetscReal       p;
  PetscReal       p_old;
  PetscReal       p_epsilon = 1.e-2;
  PetscInt        altminit=1;
  Vec             Vold,U_s,U_1;
  PetscReal       vol_s,vol_1;
  PetscReal       errV=1e+10;
  PetscReal       flowrate,maxvol = .03,minvol = 0.;
  
  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);

  ierr = PetscOptionsGetInt(PETSC_NULL,"-mode",&mode,PETSC_NULL);CHKERRQ(ierr);
  
  ierr = PetscOptionsGetReal(PETSC_NULL,"-maxvol",&maxvol,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-minvol",&minvol,PETSC_NULL);CHKERRQ(ierr);
  /*
   Overwrite ctx.maxtimestep with something more reasonable
   */
  ctx.maxtimestep = 151;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-maxtimestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr);
  flowrate = (maxvol - minvol) / (ctx.maxtimestep-1);
  
  /*
   Reset all BC for U and V
   */
  for (i = 0; i < 6; i++) {
    ctx.bcV[0].face[i] =NONE;
    for (j = 0; j < 3; j++) {
      ctx.bcU[j].face[i] = NONE;
    }
  }
  for (i = 0; i < 12; i++) {
    ctx.bcV[0].edge[i] =NONE;
    for (j = 0; j < 3; j++) {
      ctx.bcU[j].edge[i] = NONE;
    }
  }
  for (i = 0; i < 8; i++) {
    ctx.bcV[0].vertex[i] =NONE;
    for (j = 0; j < 3; j++) {
      ctx.bcU[j].vertex[i] = NONE;
    }
  }
  switch (mode) {
    case 0:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Doing a 3D experiment with null-displacement on all faces\n");CHKERRQ(ierr);     
      /*  face X0 */
      ctx.bcU[0].face[X0] = ZERO;
      ctx.bcU[1].face[X0] = ZERO;
      ctx.bcU[2].face[X0] = ZERO;
      /*  face X1 */
      ctx.bcU[0].face[X1] = ZERO;
      ctx.bcU[1].face[X1] = ZERO;
      ctx.bcU[2].face[X1] = ZERO;
      /*  face Y0 */
      ctx.bcU[0].face[Y0] = ZERO;
      ctx.bcU[1].face[Y0] = ZERO;
      ctx.bcU[2].face[Y0] = ZERO;
      /*  face Y1 */
      ctx.bcU[0].face[Y1] = ZERO;     
      ctx.bcU[1].face[Y1] = ZERO;     
      ctx.bcU[2].face[Y1] = ZERO;     
      /*  face Z0 */
      ctx.bcU[0].face[Z0] = ZERO;
      ctx.bcU[1].face[Z0] = ZERO;
      ctx.bcU[2].face[Z0] = ZERO;
      /*  face Z1 */
      ctx.bcU[0].face[Z0] = ZERO;
      ctx.bcU[1].face[Z0] = ZERO;
      ctx.bcU[2].face[Z1] = ZERO;
      break;
    case 1:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Doing a 2D experiment with null-displacement on all faces\n");CHKERRQ(ierr);     
      /*  face X0 */
      ctx.bcU[0].face[X0] = ZERO;
      ctx.bcU[1].face[X0] = ZERO;
      ctx.bcU[2].face[X0] = ZERO;
      /*  face X1 */
      ctx.bcU[0].face[X1] = ZERO;
      ctx.bcU[1].face[X1] = ZERO;
      ctx.bcU[2].face[X1] = ZERO;
      /*  face Y0 */
      ctx.bcU[1].face[Y0] = ZERO;
      /*  face Y1 */
      ctx.bcU[1].face[Y1] = ZERO;     
      /*  face Z0 */
      ctx.bcU[0].face[Z0] = ZERO;
      ctx.bcU[1].face[Z0] = ZERO;
      ctx.bcU[2].face[Z0] = ZERO;
      /*  face Z1 */
      ctx.bcU[0].face[Z1] = ZERO;
      ctx.bcU[1].face[Z1] = ZERO;
      ctx.bcU[2].face[Z1] = ZERO;
      break;
    case 2:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Doing a 3D experiment blocking only rigid motions (may be instable)\n");CHKERRQ(ierr);     
      ctx.bcU[0].vertex[X0Y0Z0] = ZERO;
      ctx.bcU[1].vertex[X0Y0Z0] = ZERO;
      ctx.bcU[2].vertex[X0Y0Z0] = ZERO;

      ctx.bcU[2].vertex[X1Y0Z0] = ZERO;
      ctx.bcU[0].vertex[X0Y1Z0] = ZERO;
      ctx.bcU[1].vertex[X0Y0Z1] = ZERO;
      break;
    case 3:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Doing a 2D experiment blocking only rigid motions (may be instable)\n");CHKERRQ(ierr);     
      ctx.bcU[0].vertex[X0Y0Z0] = ZERO;
      ctx.bcU[2].vertex[X0Y0Z0] = ZERO;

      ctx.bcU[2].vertex[X1Y0Z0] = ZERO;

      ctx.bcU[1].face[Y0] = ZERO;
      ctx.bcU[1].face[Y1] = ZERO;
      break;
    case 4:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Doing a 3D blocking normal displacement on 3 planes (may lead to asymmetric solutions)\n");CHKERRQ(ierr);      
      ctx.bcU[0].face[X0] = ZERO;
      ctx.bcU[1].face[Y0] = ZERO;
      ctx.bcU[2].face[Z0] = ZERO;
      break;
    case 5:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Doing a 2D blocking normal displacement on 2 planes (may lead to asymmetric solutions)\n");CHKERRQ(ierr);      
      ctx.bcU[0].face[X0] = ZERO;
      ctx.bcU[2].face[Z0] = ZERO;
      ctx.bcU[1].face[Y0] = ZERO;
      ctx.bcU[1].face[Y1] = ZERO;
      break;
    default:
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: mode should be one of {0,1,2,3,4,5,6,7},got %i\n",mode);
      break;
  }  

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
  ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
  
  //  ctx.hasCrackPressure = PETSC_TRUE;
  ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
  ierr = VecDuplicate(fields.U,&U_s);CHKERRQ(ierr);
  ierr = VecDuplicate(fields.U,&U_1);CHKERRQ(ierr);
  
  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
  //ierr = PetscViewerFileSetMode(viewer,FILE_MODE_APPEND);CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.pres",ctx.prefix);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,filename);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"#Time step \t Volume \t Pressure \t SurfaceEnergy \t ElasticEnergy \t PressureForces \t TotalMechEnergy \n");CHKERRQ(ierr);
  
  p = 1.e-5;
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
  ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
  
  ierr = VecSet(U_s,0.0);CHKERRQ(ierr);
  ierr = VecSet(U_1,0.0);CHKERRQ(ierr);
  ctx.matprop[0].beta = 0.;
  ctx.timevalue = 0;
  
  for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Time step %i, injected volume %g\n",ctx.timestep,flowrate*ctx.timestep);CHKERRQ(ierr);
    ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr); 
    altminit = 0.;
    do {
      p_old = p;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  Time step %i, alt min step %i with pressure %g, crackvolume %g\n",ctx.timestep,altminit, p, ctx.CrackVolume);CHKERRQ(ierr);
      /* 
       Update the pressure based on the relation
       V = p vol_1 + vol_s 
       with
       vol_1 = \int U_1 \cdot \nabla V
       where U_1 is the displacement field associated with null in-situ stress and unit pressure, and
       vol_s = \int U_s \cdot \nabla V
       where U_s is the displacement field associated with null pressure and in-situ stress
       */
      ctx.hasCrackPressure = PETSC_FALSE;
      ctx.hasInsitu = PETSC_TRUE;
      ierr = VecCopy(U_s,fields.U);CHKERRQ(ierr);
      ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
      ierr = VolumetricCrackOpening(&vol_s,&ctx,&fields);CHKERRQ(ierr);
      ierr = VecCopy(fields.U,U_s);CHKERRQ(ierr);
      
      ctx.hasCrackPressure = PETSC_TRUE;
      ctx.hasInsitu = PETSC_FALSE;
      //ierr = VecCopy(U_1,fields.U);CHKERRQ(ierr);
      ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
      ierr = VecScale(fields.U,1./p);CHKERRQ(ierr);
      ierr = VolumetricCrackOpening(&vol_1,&ctx,&fields);CHKERRQ(ierr);         
      ierr = VecCopy(fields.U,U_1);CHKERRQ(ierr);

      ctx.hasCrackPressure = PETSC_FALSE;
      ctx.hasInsitu = PETSC_TRUE;
      ierr = VecCopy(U_s,fields.U);CHKERRQ(ierr);
      ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
      ierr = VolumetricCrackOpening(&vol_s,&ctx,&fields);CHKERRQ(ierr);   
      ierr = VecCopy(fields.U,U_s);CHKERRQ(ierr);
      
      // This will fail if vol_1 = 0, which should only happen when there are no cracks
      p = (minvol + flowrate*ctx.timestep - vol_s) / vol_1;
      ierr = VecAXPY(fields.U,p,U_1);CHKERRQ(ierr);

      ctx.CrackVolume = vol_s + p * vol_1;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"      vol_s: %e vol_1 %e volume %e (target vol: %e)\n",vol_s,vol_1,ctx.CrackVolume,minvol + flowrate*ctx.timestep);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"      Updated crack pressure: %e (was %e)\n",p,p_old);
      
      ctx.hasCrackPressure = PETSC_TRUE;
      ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
      ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
      ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);
      
      ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
      ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"      Max. change on V: %e\n",errV);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"      Max. change on p: %e\n",PetscAbs(p-p_old));CHKERRQ(ierr);
      altminit++;
    } while (PetscAbs(p-p_old) >= p_epsilon && altminit <= ctx.altminmaxit);
    ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);   
    switch (ctx.fileformat) {
      case FILEFORMAT_HDF5:       
        ierr = FieldsH5Write(&ctx,&fields);CHKERRQ(ierr);
        break;
      case FILEFORMAT_BIN:
        ierr = FieldsBinaryWrite(&ctx,&fields);CHKERRQ(ierr);
        break; 
    }
    ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.log",ctx.prefix);CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&logviewer);CHKERRQ(ierr);
    ierr = PetscLogView(logviewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&logviewer);CHKERRQ(ierr);
    
    
    ctx.ElasticEnergy=0;
    ctx.InsituWork=0;
    ctx.PressureWork = 0.;
    ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,&fields,&ctx);CHKERRQ(ierr);
    ierr = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);
    ctx.TotalEnergy = ctx.ElasticEnergy - ctx.InsituWork - ctx.PressureWork + ctx.SurfaceEnergy;
    ierr = PetscViewerASCIIPrintf(viewer,"%d \t\t %e \t %e \t %e \t %e \t %e \t %e\n",ctx.timestep ,ctx.CrackVolume,p,ctx.SurfaceEnergy,ctx.ElasticEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%i   \t%e   \t%e   \t%e   \t%e   \t%e\n",ctx.timestep,ctx.ElasticEnergy,
                    ctx.InsituWork,ctx.SurfaceEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);
  }
  ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&logviewer);CHKERRQ(ierr);
  ierr = VecDestroy(&Vold);CHKERRQ(ierr);
  ierr = VecDestroy(&U_s);CHKERRQ(ierr);
  ierr = VecDestroy(&U_1);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");
  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}

