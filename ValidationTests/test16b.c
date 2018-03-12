/* 

Example

 
   y
   |
   |_____x
   /
  /
 z                             
       ________________________
    ->|                       |<-      
    ->|                       |<-      
    ->|         \             |<-      
    ->|          \            |<-
    ->|           \           |<-
    ->|                       |<-
    ->|_______________________|<-

--------------------------------------------------------------------------------------------------------
srun  test16b -n 100,100,2 -l  -l 1,1,.1 -E 1 -nu 0 -U_snes_monitor -p test16b         \
             -U_X0Y0Z0_BC_0 ZERO  -U_X0Y0Z0_BC_1 ZERO  -U_X1Y0Z0_BC_1 ZERO  -U_Z0_BC_2 ZERO  -U_Z1_BC_2 ZERO                               \
             -U_tao_type ntr -U_tao_monitor -U_tao_fatol 1e-9 -U_tao_frtol 1e-9                           \
             -type lmvm -U_tao_monitor -U_tao_fatol 1e-9 -U_tao_frtol 1e-9                          \
             -V_X0_BC ONE -V_X1_BC ONE -V_Y0_BC ONE -V_Y1_BC ONE                                      \
             -insitumin -1,0,0,0,0,0  -insitumax -1,0,0,0,0,0 -npc 1                                    \
             -minvol .002 -maxvol .03  -maxtimestep   16                                                  \
             -pc0_r .2 -pc0_thickness .015 -pc0_center .5,.5,0.05 -pc0_phi 90  -pc0_theta 30        \
             -epsilon .03 -verbose 0 -eta 1e-8 -atnum 1 -unilateral nocompression -alpha 0 -beta 0 -U_snes_monitor

---------------------------------------------------------------------------------------------

Modified by Erwan Tanne (erwan.tanne@gmail.com) , include a secant method to determine the pressure of the penny crack subject to external pressure

Remarks, 
-Do not works well for vol min equal to zero
-Should not have volume increment to small
-Job running took a long time

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
  PetscViewer    secant_viewer;
  PetscViewer    logviewer;
  char           filename[FILENAME_MAX],filenameC[FILENAME_MAX];
  PetscReal      p;
  PetscReal      Tol=0.01;
  /*
  new variables, 
  pa , vol_a  ,   U_a
  pb , vol_b  ,   U_b
  P  , vol    ,   U_s
  */
  PetscReal      pa,pb,pc;
  PetscInt       k_secant=0;
  PetscReal      p_old;
  PetscReal      prestol = 1.e-2;
  PetscInt       altminit  =1;
  Vec            Vold,U_b,U_a,U_s;
  PetscReal      vol_b,vol_a,vol;
  PetscReal      errV=1e+10,errP;
  PetscReal      flowrate,maxvol = .03,minvol = 0.;
  PetscReal      targetVol;
  PetscBool      debug=PETSC_FALSE;
  PetscBool      saveall=PETSC_FALSE;

  /* PetscBool      secant=PETSC_FALSE;       definir une boucle qui fait la secante    */

  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,"-maxvol",&maxvol,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,"-minvol",&minvol,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,"-prestol",&prestol,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,"-debug",&debug,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,"-saveall",&saveall,NULL);CHKERRQ(ierr);
  /*
    Overwrite ctx.maxtimestep with something more reasonable
  */
  ctx.maxtimestep = 150;
  ierr            = PetscOptionsGetInt(NULL,"-maxtimestep",&ctx.maxtimestep,NULL);CHKERRQ(ierr);
  flowrate        = (maxvol - minvol) / (ctx.maxtimestep-1);

  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);

  ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
  ierr = VecDuplicate(fields.U,&U_b);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) U_b,"U_b");CHKERRQ(ierr);
  ierr = VecDuplicate(fields.U,&U_a);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) U_a,"U_a");CHKERRQ(ierr);

  ierr = VecDuplicate(fields.U,&U_s);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) U_s,"U_s");CHKERRQ(ierr);

  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.pres",ctx.prefix);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,filename);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"#Time step \t Target_volume \t Volume \t Pressure \t SurfaceEnergy \t ElasticEnergy \t PressureForces \t TotalMechEnergy \n");CHKERRQ(ierr);

  p     = 1.e-5;
  p_old = 1.e-5;
  ierr  = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr  = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr  = VecSet(fields.pressure,p);CHKERRQ(ierr);
  ierr  = VecSet(fields.U,0.0);CHKERRQ(ierr);
  ierr  = VecSet(U_b,0.0);CHKERRQ(ierr);
  ierr  = VecSet(U_a,0.0);CHKERRQ(ierr);
  ierr  = VecSet(U_s,0.0);CHKERRQ(ierr);

  ctx.matprop[0].beta  = 0.;
  ctx.matprop[0].alpha = 0.;
  ctx.timevalue        = 0;

  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&secant_viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(secant_viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.sec",ctx.prefix);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(secant_viewer,filename);CHKERRQ(ierr);




  for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++) {
    targetVol = minvol + flowrate * ctx.timestep;
    altminit  = 0.;
    ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Time step %i. Targeting injected volume of %g\n",ctx.timestep,targetVol);CHKERRQ(ierr);
    ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
    
    /*
    Initialisation pressure for the secant method
    */

    pc=p;


    if (ctx.timestep != 0)
    {
    pa=.9*pc;
    pb=1.1*pc;
    }
    else
    pa=0.0;
    pb=2.0;



    do {
      p_old = p;
      ierr  = PetscPrintf(PETSC_COMM_WORLD,"  Time step %i, alt min step %i with pressure %g\n",ctx.timestep,altminit,p);CHKERRQ(ierr);

      /*----------------------
          Start Secant method
       ------------------------     
            ->Solving V_a for pa

            ->Solving V_b for pb

            ->estimation of p

            ->Solving V for p

            -> While ( Tol>0.01 ) 
                      -> if V_b=V_a  then multipy pb by two. (because p is divided by zero)
                      -> do it again with the new pressures values 
      */



      ctx.hasCrackPressure = PETSC_TRUE;
      ctx.hasInsitu        = PETSC_TRUE;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"     Solving for Vol_a \n");CHKERRQ(ierr);
      ierr = VecSet(fields.pressure,pa);CHKERRQ(ierr);
      ierr = VecCopy(U_a,fields.U);CHKERRQ(ierr);
      ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
      ierr = VolumetricCrackOpening(&vol_a,&ctx,&fields);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," pa \t %e vol_a \t %e \n",pa,vol_a);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(secant_viewer,"timestep \t %d \t altminit \t %d \t  pa \t %e \t vol_a \t %e \n ",ctx.timestep,altminit,pa,vol_a);CHKERRQ(ierr);
      ierr = VecCopy(fields.U,U_a);CHKERRQ(ierr);




      ctx.hasCrackPressure = PETSC_TRUE;
      ctx.hasInsitu        = PETSC_TRUE;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"     Solving for Vol_b");CHKERRQ(ierr);
      ierr = VecSet(fields.pressure,pb);CHKERRQ(ierr);
      ierr = VecCopy(U_b,fields.U);CHKERRQ(ierr);
      ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr); 
      ierr = VolumetricCrackOpening(&vol_b,&ctx,&fields);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," pb \t %e vol_b \t %e \n",pb,vol_b);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(secant_viewer,"timestep \t %d \t altminit \t %d \t  pb \t %e \t vol_b \t %e \n ",ctx.timestep,altminit,pb,vol_b);CHKERRQ(ierr);
      ierr = VecCopy(fields.U,U_b);CHKERRQ(ierr);
  
      p = (targetVol-vol_a)/(vol_b-vol_a) * (pb-pa) + pa;

      ctx.hasCrackPressure = PETSC_TRUE;
      ctx.hasInsitu        = PETSC_TRUE;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"     Solving for Vol");CHKERRQ(ierr);
      ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
      ierr = VecCopy(U_s,fields.U);CHKERRQ(ierr);
      ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
      ierr = VolumetricCrackOpening(&vol,&ctx,&fields);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," p \t %e vol \t %e \n",p,vol);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(secant_viewer,"timestep \t %d \t altminit \t %d \t  p \t %e \t vol \t %e \n ",ctx.timestep,altminit,p,vol);CHKERRQ(ierr);
      ierr = VecCopy(fields.U,U_s);CHKERRQ(ierr);

      pa=pb;
      vol_a=vol_b;
      pb=p;
      vol_b=vol;
  
      ierr = PetscPrintf(PETSC_COMM_WORLD,"     Vol  \t %e  targetVol  \t %e",vol,targetVol);

      
      Tol=PetscAbs((vol - targetVol)/targetVol);

      ierr = PetscPrintf(PETSC_COMM_WORLD,"     Tol \t %e ",Tol);
  
      while ((k_secant<100) && Tol > .01)
        { 
        k_secant++;

        while (vol_b==vol_a)
          {
          pb=2*pb;
          ierr = PetscPrintf(PETSC_COMM_WORLD,"  vol_a==vol_b   Solving for Vol_b");CHKERRQ(ierr);
          ierr = VecSet(fields.pressure,pb);CHKERRQ(ierr);
          ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
          ierr = VolumetricCrackOpening(&vol_b,&ctx,&fields);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD," pb \t %e vol_b \t %e \n",pb,vol_b);CHKERRQ(ierr);
          }

        p = (targetVol-vol_b)/(vol_b-vol_a) * (pb-pa) + pb;

        ctx.hasCrackPressure = PETSC_TRUE;
        ctx.hasInsitu        = PETSC_TRUE;
        ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
        ierr = VecCopy(U_s,fields.U);CHKERRQ(ierr);
        ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
        ierr = VolumetricCrackOpening(&vol,&ctx,&fields);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Dans la boucle p \t %e vol \t %e \n",p,vol);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(secant_viewer,"timestep \t %d \t altminit \t %d \t  p \t %e \t vol \t %e \n ",ctx.timestep,altminit,p,vol);CHKERRQ(ierr);
        ierr = VecCopy(fields.U,U_s);CHKERRQ(ierr);

        pa=pb;
        vol_a=vol_b;
        pb=p;
        vol_b=vol;
        Tol=PetscAbs((vol - targetVol)/targetVol);
        }
       
      /*----------------------
          End Secant method
      ----------------------*/

      /*
         This will fail if vol_a = 0, which should only happen when there are no cracks
      */
      /* p    = (targetVol - vol_b) / vol_a;*/
      errP = PetscAbs((p-p_old)/p);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"     Updated crack pressure: %e\n",p);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"     Rel. change on p: %e\n",errP);CHKERRQ(ierr);

      /*ierr = VecAXPY(fields.U,p,U_a);CHKERRQ(ierr);*/

      ctx.CrackVolume      = targetVol;
      ctx.hasCrackPressure = PETSC_TRUE;

      ierr = PetscPrintf(PETSC_COMM_WORLD,"     Solving for V  ",ctx.CrackVolume);CHKERRQ(ierr);
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

        ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,fields.U,&ctx);CHKERRQ(ierr);
        ierr = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);

        ctx.TotalEnergy   = ctx.ElasticEnergy - ctx.InsituWork - ctx.PressureWork + ctx.SurfaceEnergy;

        ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface energy: %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer,"%d \t\t %e \t %e \t %e \t %e \t %e \t %e \t %e\n",ctx.timestep,targetVol,ctx.CrackVolume,p,ctx.SurfaceEnergy,ctx.ElasticEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);
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
    ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,fields.U,&ctx);CHKERRQ(ierr);
    ierr = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);

    ctx.TotalEnergy   = ctx.ElasticEnergy - ctx.InsituWork - ctx.PressureWork + ctx.SurfaceEnergy;

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface energy: %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"%d \t\t %e \t  %e \t %e \t %e \t %e \t %e \t %e\n",ctx.timestep,targetVol ,ctx.CrackVolume,p,ctx.SurfaceEnergy,ctx.ElasticEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%i   \t%e   \t%e   \t%e   \t%e   \t%e\n",ctx.timestep,ctx.ElasticEnergy,
                                  ctx.InsituWork,ctx.SurfaceEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);
  }
  ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ierr = VecDestroy(&Vold);CHKERRQ(ierr);
  ierr = VecDestroy(&U_b);CHKERRQ(ierr);
  ierr = VecDestroy(&U_a);CHKERRQ(ierr);
  ierr = VecDestroy(&U_s);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");
  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}
