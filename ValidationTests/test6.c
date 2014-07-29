/*
  test6.c:
  Validate elasticity solver 

  (c) 2010-2014 Blaise Bourdin bourdin@lsu.edu
  
  Try: 
mpiexec -n 4 ./test6 -n 101,2,101 -l 1,.1,1 -E 1 -nu 0 -U_snes_monitor -orientation 2 \
             -insitumin -1e5,0,0,0,0,0  -insitumax -1e5,0,0,0,0,0 -npc 1              \
             -pc0_r .4 -pc0_thickness .02 -pc0_center 0.5,0.01,0.5 -pc0_phi 90        \
             -epsilon .03 -verbose 0 -eta 0 -bcu 0,0,0  -eta .001 -atnum 1            \
             -unilateral nocompression -pressure 5e5

mpiexec -n 4 ./test6 -n 101,2,101 -l 1,.1,1 -E 1 -nu 0 -U_snes_monitor -orientation 2 \
             -insitumin -1,0,0,0,0,0  -insitumax -1,0,0,0,0,0 -npc 3                  \
             -pc0_r .4 -pc0_thickness .02 -pc0_center 0.5,0.01,0.5 -pc0_phi 45        \
             -pc1_r .4 -pc1_thickness .02 -pc1_center 0.4,0.01,0.5 -pc1_phi 90        \
             -pc2_r .4 -pc2_thickness .02 -pc2_center 0.6,0.01,0.5 -pc2_phi 90        \
             -epsilon .03 -verbose 0 -eta 0 -bcu 0,0,-.2  -eta .001 -atnum 1          \
             -unilateral nocompression -U_pc_type hypre
             
mpiexec -n 4 ./test6 -n 101,2,101 -l 1,.1,1 -E 1 -nu 0 -U_snes_monitor -orientation 2 \
             -insitumin -1,0,0,0,0,0  -insitumax -1,0,0,0,0,0 -npc 1                  \
             -pc0_r .4 -pc0_thickness .02 -pc0_center 0.5,0.01,0.5 -pc0_phi 90        \
             -epsilon .03 -verbose 0 -eta 0 -bcu 0,0,0  -eta .001 -atnum 1            \
             -unilateral nocompression -U_pc_type hypre

mpiexec -n 4 ./test6 -n 101,2,101 -l 1,.1,1 -E 1 -nu 0 -U_snes_monitor                \
             -options_file ../BCfiles/2DXZRigidMotion.txt                             \
             -insitumin -1,0,0,0,0,0  -insitumax -1,0,0,0,0,0 -npc 1                  \
             -pc0_r .4 -pc0_thickness .02 -pc0_center 0.5,0.01,0.5 -pc0_phi 90        \
             -epsilon .03 -verbose 0 -eta .001 -atnum 1                               \
             -unilateral nocompression -U_pc_type hypre

mpiexec -n 4 ./test6 -n 101,2,101 -l 1,.1,1 -E 1 -nu 0 -U_snes_monitor                \
             -options_file ../BCfiles/2DXZRigidMotion.txt                             \
             -insitumin -1,0,0,0,0,0  -insitumax -1,0,0,0,0,0 -npc 1                  \
             -pc0_r .25 -pc0_thickness .02 -pc0_center 0.5,0.01,0.5 -pc0_phi 45       \
             -epsilon .03 -verbose 0 -eta 1e-8 -atnum 1                               \
             -unilateral nocompression -U_pc_type hypre -pressure .8

mpiexec -n 4 ./test6 -n 101,2,101 -l 1,.1,1 -E 1 -nu 0 -U_snes_monitor                \
             -options_file ../BCfiles/2DXZRigidMotion.txt                             \
             -insitumin -1,0,-1,0,0,0  -insitumax -1,0,-1,0,0,0 -npc 3                \
             -pc0_r .3 -pc0_thickness .02 -pc0_center 0.5,0.01,0.5 -pc0_phi 00        \
             -pc1_r .2 -pc1_thickness .02 -pc1_center 0.65,0.01,0.5 -pc1_phi 90       \
             -pc2_r .2 -pc2_thickness .02 -pc2_center 0.35,0.01,0.5 -pc2_phi 90       \
             -epsilon .03 -verbose 0 -eta 1e-8 -atnum 1                               \
             -unilateral nocompression -U_pc_type hypre -pressure .8

mpiexec -n 4 ./test6 -n 101,2,101 -l 1,.1,1 -E 1 -nu 0 -U_snes_monitor                \
             -options_file ../BCfiles/2DXZRigidMotion.txt                             \
             -insitumin -1,0,-1.5,0,0,0  -insitumax -1,0,-1.5,0,0,0 -npc 2            \
             -pc0_r .3 -pc0_thickness .02 -pc0_center 0.5,0.01,0.5 -pc0_phi 00        \
             -pc1_r .3 -pc1_thickness .02 -pc1_center 0.5,0.01,0.5 -pc1_phi 90        \
             -epsilon .03 -verbose 0 -eta 1e-4 -atnum 1                               \
             -unilateral nocompression -U_pc_type hypre -pressure 1.25

mpiexec -n 4 ./test6 -n 100,100,2 -l 1,1,.1 -E 1 -nu 0 -U_snes_monitor                                \
             -options_file ../BCfiles/2DXYRigidMotion.txt                                             \
             -insitumin -1,-2,0,0,0,0  -insitumax -1,-2,0,0,0,0 -npc 2                                \
             -pc0_r .3 -pc0_thickness .015 -pc0_center 0.5,0.5,0.01 -pc0_phi 90  -pc0_theta 00        \
             -pc1_r .3 -pc1_thickness .015 -pc1_center 0.5,0.5,0.01 -pc1_phi 90  -pc1_theta 90       \
             -U_pc_type hypre  -u_pc_hypre_type boomeramg -u_pc_hypre_boomeramg_strong_threshold 0.7  \
             -epsilon .03 -verbose 0 -eta 1e-8 -atnum 1 -unilateral nocompression  -pressure 1.5

mpiexec -n 4 ./test6 -n 100,100,2 -l 1,1,.1 -E 1 -nu 0 -U_snes_monitor                                \
             -options_file ../BCfiles/2DXYRigidMotion.txt                                             \
             -insitumin -1,0,0,0,0,0  -insitumax -1,0,0,0,0,0 -npc 1                                  \
             -pc0_r .3 -pc0_thickness .015 -pc0_center 0.5,0.5,0.01 -pc0_phi 90  -pc0_theta 00        \
             -U_pc_type hypre  -u_pc_hypre_type boomeramg -u_pc_hypre_boomeramg_strong_threshold 0.7  \
             -epsilon .03 -verbose 0 -eta 1e-8 -atnum 1 -unilateral nocompression  -pressure 1.5      \
             -maxtimestep 2 
*/

#include "petsc.h"
#include "petsctao.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFMech.h"
#include "VFFlow.h"
#include "VFPermfield.h"


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  VFCtx          ctx;
  VFFields       fields;
  PetscErrorCode ierr;

  PetscReal p = 0.;
  PetscBool flg;
  PetscReal crackVolume;

  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);

  ierr = PetscOptionsGetReal(NULL,"-pressure",&p,&flg);CHKERRQ(ierr);
  if (flg) {
    ctx.hasCrackPressure = PETSC_TRUE;
  }

  ctx.matprop[0].alpha = 0.;
  ierr = VecSet(fields.U,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);

  /*
    VF BC for U and V are now setup in VFCommon.c
  */

  ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"#p,volume,Elastic Energy,InsituWork,Pressure Work,Total Energy\n");CHKERRQ(ierr);
  for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++) {
    if (ctx.maxtimestep > 1) {
      ctx.timevalue = (PetscReal) ctx.timestep / (ctx.maxtimestep - 1.) * p;
    } else {
      ctx.timevalue = p;
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Solving for p = %g\n",ctx.timevalue);CHKERRQ(ierr);
    ierr = VecSet(fields.pressure,ctx.timevalue);CHKERRQ(ierr);
    ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
    ierr = VF_StepU(&fields,&ctx);

    ctx.ElasticEnergy = 0;
    ctx.InsituWork    = 0;
    ctx.PressureWork  = 0.;
    ierr              = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,&fields,&ctx);CHKERRQ(ierr);
    ctx.TotalEnergy   = ctx.ElasticEnergy-ctx.InsituWork-ctx.PressureWork;

    ierr = PetscPrintf(PETSC_COMM_WORLD,  "Elastic Energy:          \t%e\n",ctx.ElasticEnergy);CHKERRQ(ierr);
    if (ctx.hasCrackPressure) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of pressure forces: \t%e\n",ctx.PressureWork);CHKERRQ(ierr);
    }
    if (ctx.hasInsitu) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of surface forces:  \t%e\n",ctx.InsituWork);CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,  "Total Mechanical energy: \t%e\n",ctx.ElasticEnergy-ctx.InsituWork-ctx.PressureWork);CHKERRQ(ierr);

    /*
      Compute and display total crack opening
    */
    ierr = VolumetricCrackOpening(&crackVolume,&ctx,&fields);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,   "Total crack opening:    \t%e\n",crackVolume);CHKERRQ(ierr);
  
    /*
      Save fields and write statistics about current run
    */
    ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%e \t%e \t%e \t%e \t%e \t%e\n",ctx.timevalue,crackVolume,ctx.ElasticEnergy,ctx.InsituWork,ctx.PressureWork,ctx.ElasticEnergy-ctx.InsituWork-ctx.PressureWork);CHKERRQ(ierr);
  }
  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}

