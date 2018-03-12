/*
 test17.c: Solves for displacement and v-field in volume driven propagation of arbitrary number of line cracks in 2d using VFRectangularCrack
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */

#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFMech.h"
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
  PetscInt        orientation=1;
  PetscInt        i,j,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal   ****coords_array;
  PetscReal       BBmin[3],BBmax[3];
  char            filename[FILENAME_MAX];
  PetscReal       p;
  PetscReal       p_old;
  PetscReal       p_epsilon = 1.e-5;
  PetscInt        altminit=1;
  Vec             Vold;
  PetscReal       errV=1e+10;
  PetscReal       lx,ly,lz;
  PetscReal       q,maxvol = .06;
  PetscInt        nc = 0;

  
  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  
  ierr = PetscOptionsGetReal(NULL,NULL,"-q",&q,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-orientation",&orientation,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-maxvol",&maxvol,NULL);CHKERRQ(ierr);
  ierr = DMDAGetInfo(ctx.daScal,NULL,&nx,&ny,&nz,NULL,NULL,NULL,
                     NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DMDAGetBoundingBox(ctx.daScal,BBmin,BBmax);CHKERRQ(ierr);  
  /*
    Overwrite ctx.maxtimestep with something more reasonable
  */
  ctx.maxtimestep = 150;
  ierr = PetscOptionsGetInt(NULL,NULL,"-maxtimestep",&ctx.maxtimestep,NULL);CHKERRQ(ierr);
  q = maxvol / ctx.maxtimestep;
  
  ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  lz = BBmax[2];
  ly = BBmax[1];
  lx = BBmax[0];  
  /*
   Reset all BC for U and V
   */
  for (i = 0; i < 6; i++) {
    ctx.bcV[0].face[i]=NONE;
    for (j = 0; j < 3; j++) {
      ctx.bcU[j].face[i] = NONE;
    }
  }
  for (i = 0; i < 12; i++) {
    ctx.bcV[0].edge[i]=NONE;
    for (j = 0; j < 3; j++) {
      ctx.bcU[j].edge[i] = NONE;
    }
  }
  for (i = 0; i < 8; i++) {
    ctx.bcV[0].vertex[i]=NONE;
    for (j = 0; j < 3; j++) {
      ctx.bcU[j].vertex[i] = NONE;
    }
  }
  /*Initializing the v-field*/
  switch (orientation) {
    case 1:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Building %g number of cracks line of arbitrary orientation\n",
                 nc);CHKERRQ(ierr);     
      /*  face X0 */
      ctx.bcU[0].face[X0]= ZERO;
      ctx.bcU[1].face[X0]= ZERO;
      ctx.bcU[2].face[X0]= ZERO;
      /*  face X1 */
      ctx.bcU[0].face[X1]= ZERO;
      ctx.bcU[1].face[X1]= ZERO;
      ctx.bcU[2].face[X1]= ZERO;
      /*  face Y0 */
      ctx.bcU[0].face[Y0]= ZERO;
      ctx.bcU[1].face[Y0]= ZERO;
      ctx.bcU[2].face[Y0]= ZERO;
      /*  face Y1 */
      ctx.bcU[0].face[Y1]= ZERO;      
      ctx.bcU[1].face[Y1]= ZERO;      
      ctx.bcU[2].face[Y1]= ZERO;      
      /*  face Z0 */
      ctx.bcU[2].face[Z0]= ZERO;
      /*  face Z1 */
      ctx.bcU[2].face[Z1]= ZERO;    
      break;
    case 2:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Building %g number of cracks line of arbitrary orientation\n",
                 nc);CHKERRQ(ierr);     
      /*  face X0 */
      ctx.bcU[0].face[X0]= ZERO;
      ctx.bcU[2].face[X0]= ZERO;
      /*  face X1 */
      ctx.bcU[0].face[X1]= ZERO;
      ctx.bcU[2].face[X1]= ZERO;
      /*  face Y0 */
      ctx.bcU[1].face[Y0]= ZERO;
      ctx.bcU[2].face[Y0]= ZERO;
      /*  face Y1 */
      ctx.bcU[1].face[Y1]= ZERO;
      ctx.bcU[2].face[Y1]= ZERO;      
      /*  face Z0 */
      ctx.bcU[2].face[Z0]= ZERO;
      /*  face Z1 */
      ctx.bcU[2].face[Z1]= ZERO;
      break;
    case 3:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Building %g number of cracks line of arbitrary orientation\n",
                 nc);CHKERRQ(ierr);     
      /*  face X0 */
      ctx.bcU[0].face[X0]= ZERO;
      ctx.bcU[1].face[X0]= ZERO;
      ctx.bcU[2].face[X0]= ZERO;
      /*  face X1 */
      ctx.bcU[0].face[X1]= ZERO;
      ctx.bcU[1].face[X1]= ZERO;
      ctx.bcU[2].face[X1]= ZERO;
      /*  face Y0 */
      ctx.bcU[1].face[Y0]= ZERO;
      ctx.bcU[2].face[Y0]= ZERO;
      /*  face Y1 */
      ctx.bcU[1].face[Y1]= ZERO;      
      ctx.bcU[2].face[Y1]= ZERO;      
      /*  face Z0 */
      ctx.bcU[2].face[Z0]= ZERO;
      /*  face Z1 */
      ctx.bcU[2].face[Z1]= ZERO;
      break;
    case 4:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Building %g number of cracks line of arbitrary orientation\n",
                 nc);CHKERRQ(ierr);     
      /*  face X0 */
      ctx.bcU[0].face[X0]= ZERO;
      ctx.bcU[2].face[X0]= ZERO;
      /*  face X1 */
      ctx.bcU[0].face[X1]= ZERO;
      ctx.bcU[2].face[X1]= ZERO;
      /*  face Y0 */
      ctx.bcU[0].face[Y0]= ZERO;
      ctx.bcU[1].face[Y0]= ZERO;
      ctx.bcU[2].face[Y0]= ZERO;
      /*  face Y1 */
      ctx.bcU[0].face[Y1]= ZERO;      
      ctx.bcU[1].face[Y1]= ZERO;      
      ctx.bcU[2].face[Y1]= ZERO;      
      /*  face Z0 */
      ctx.bcU[2].face[Z0]= ZERO;
      /*  face Z1 */
      ctx.bcU[2].face[Z1]= ZERO;
      break;
    case 5:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Building %g number of cracks line of arbitrary orientation\n",
                 nc);CHKERRQ(ierr);     
      /*  face X0 */
      ctx.bcU[0].face[X0]= ZERO;
      ctx.bcU[2].face[X0]= ZERO;
      /*  face X1 */
      ctx.bcU[0].face[X1]= ZERO;
      ctx.bcU[2].face[X1]= ZERO;
      /*  face Y0 */
        /*.....FREE.......*/
      /*  face Y1 */
        /*.....FREE.......*/
      /*  face Z0 */
      ctx.bcU[2].face[Z0]= ZERO;
      /*  face Z1 */
      ctx.bcU[2].face[Z1]= ZERO;
      break;
    case 6:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Building %g number of cracks line of arbitrary orientation\n",
                 nc);CHKERRQ(ierr);     
      /*  face X0 */
      ctx.bcU[0].face[X0]= ZERO;
      ctx.bcU[1].face[X0]= ZERO;
      ctx.bcU[2].face[X0]= ZERO;
      /*  face X1 */
      ctx.bcU[0].face[X1]= ZERO;
      ctx.bcU[1].face[X1]= ZERO;
      ctx.bcU[2].face[X1]= ZERO;
      /*  face Y0 */
      /*.....FREE.......*/
      /*  face Y1 */
      /*.....FREE.......*/
      /*  face Z0 */
      ctx.bcU[2].face[Z0]= ZERO;
      /*  face Z1 */
      ctx.bcU[2].face[Z1]= ZERO;
      break;
    case 7:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Building %g number of cracks line of arbitrary orientation\n",
                 nc);CHKERRQ(ierr);     
      /*  face X0 */
      ctx.bcU[0].face[X0]= ZERO;
      ctx.bcV[0].face[X0]= ONE;
      //ctx.bcU[1].face[X0]= ZERO;
      //ctx.bcU[2].face[X0]= ZERO;
      /*  face X1 */
      ctx.bcU[0].face[X1]= ZERO;
      ctx.bcV[0].face[X1]= ONE;
      //ctx.bcU[1].face[X1]= ZERO;
      //ctx.bcU[2].face[X1]= ZERO;
      /*  face Y0 */
      ctx.bcU[1].face[Y0]= ZERO;
      //ctx.bcV[0].face[Y0]= FIXED;
      /*  face Y1 */
      ctx.bcU[1].face[Y1]= ZERO;
      //ctx.bcV[0].face[Y1]= FIXED;
      /*  face Z0 */
      ctx.bcU[2].face[Z0]= ZERO;
      ctx.bcV[0].face[Z0]= ONE;
      /*  face Z1 */
      //ctx.bcU[2].face[Z1]= ZERO;
      ctx.bcV[0].face[Z1]= ONE;
      break;
    default:
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Orientation should be one of {1,2,3,4,5,6,7}, got %i\n",orientation);
      break;
  }  
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr); 
  ctx.hasCrackPressure = PETSC_TRUE;
  ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_APPEND);CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.pres",ctx.prefix);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,filename);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer, "#Time step \t Volume \t Pressure \t SurfaceEnergy \t ElasticEnergy \t PressureForces \t TotalMechEnergy \n");CHKERRQ(ierr);
  
  p = 1.e-5;
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
  ctx.matprop[0].beta = 0.;
  ctx.timevalue = 0;
  ctx.maxtimestep = 2;
  for (ctx.timestep = 1; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Time step %i, injected volume %g\n",ctx.timestep,q*ctx.timestep);CHKERRQ(ierr);
    ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr); 
    do {
      p_old = p;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  Time step %i, alt min step %i with pressure %g\n",ctx.timestep,altminit, p);CHKERRQ(ierr);
      ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
      ierr = VecScale(fields.U,1./p);CHKERRQ(ierr);
      ierr = VolumetricCrackOpening(&ctx.CrackVolume, &ctx, &fields);CHKERRQ(ierr);   
      p    = q*ctx.timestep / ctx.CrackVolume;
      ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
      ierr = VecScale(fields.U,p);CHKERRQ(ierr);
      ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
      ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);

      ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
      ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"    Max. change on V: %e\n",errV);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"    Max. change on p: %e\n", PetscAbs(p-p_old));CHKERRQ(ierr);
      altminit++;
    } while (PetscAbs(p-p_old) >= p_epsilon && altminit <= ctx.altminmaxit);
    ierr = VolumetricCrackOpening(&ctx.CrackVolume, &ctx, &fields);CHKERRQ(ierr);   
    switch (ctx.fileformat) {
      case FILEFORMAT_VTK:       
        ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
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
    ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,fields.U,&ctx);CHKERRQ(ierr);
    ierr = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);
    ctx.TotalEnergy = ctx.ElasticEnergy - ctx.InsituWork - ctx.PressureWork + ctx.SurfaceEnergy;
    ierr = PetscViewerASCIIPrintf(viewer, "%d \t %e \t %e \t %e \t %e \t %e \t %e\n", ctx.timestep , ctx.CrackVolume, p, ctx.SurfaceEnergy, ctx.ElasticEnergy, ctx.PressureWork, ctx.TotalEnergy);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%i   \t%e   \t%e   \t%e   \t%e   \t%e\n",ctx.timestep,ctx.ElasticEnergy,
                                  ctx.InsituWork,ctx.SurfaceEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);
    altminit = 0.;
  }
  ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&logviewer);CHKERRQ(ierr);
  ierr = VecDestroy(&Vold);CHKERRQ(ierr);

  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}

