/*
  test14.c:
  SPE Case2: Crack initiation test

  (c) 2010-2012 Keita Yoshioka yoshk@chevron.com
*/

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFV.h"
#include "VFU.h"
#include "VFFlow.h"

VFCtx    ctx;
VFFields fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  VFCtx          ctx;
  VFFields       fields;
  PetscErrorCode ierr;

  PetscReal length      = .2;
  PetscInt  orientation = 0;
  PetscInt  nopts       = 3;
  PetscInt  i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal p,pinc;
  PetscReal ****coords_array;
  PetscReal ****bcu_array;
  PetscReal ***v_array;
  PetscReal BBmin[3],BBmax[3];
  PetscReal ElasticEnergy = 0;
  PetscReal InsituWork    = 0;
  PetscReal SurfaceEnergy = 0;
  PetscReal pinit;
  char      filename[FILENAME_MAX];



  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);

  ierr            = PetscOptionsGetReal(PETSC_NULL,"-length",&length,PETSC_NULL);CHKERRQ(ierr);
  pinc            = 1e-6;
  ierr            = PetscOptionsGetReal(PETSC_NULL,"-pinc",&pinc,PETSC_NULL);CHKERRQ(ierr);
  pinit            = 0.;
  ierr            = PetscOptionsGetReal(PETSC_NULL,"-pinit",&pinit,PETSC_NULL);CHKERRQ(ierr);
  ctx.maxtimestep = 1;
  ierr            = PetscOptionsGetInt(PETSC_NULL,"-maxtimestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr);
  ierr            = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                                PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);

  p                    = pinit;
  ctx.matprop[0].beta  = 0.;
  ctx.matprop[0].alpha = 0.;

  ctx.timestep  = 0;
  ctx.timevalue = 1.;

  ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);

  ierr = VecSet(fields.V,1.0);CHKERRQ(ierr);
  ierr = VecSet(fields.VIrrev,1.0);CHKERRQ(ierr);
  ierr = VecSet(fields.U,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
  ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(ctx.daScal,fields.VIrrev,&v_array);CHKERRQ(ierr);


  ctx.bcU[0].face[X0] = ZERO;ctx.bcU[1].face[X0] = NONE;ctx.bcU[2].face[X0] = NONE;ctx.bcV[0].face[X0] = NONE;
  ctx.bcU[0].face[X1] = ZERO;ctx.bcU[1].face[X1] = NONE;ctx.bcU[2].face[X1] = NONE;ctx.bcV[0].face[X1] = NONE;
  ctx.bcU[0].face[Y0] = NONE;ctx.bcU[1].face[Y0] = ZERO;ctx.bcU[2].face[Y0] = NONE;ctx.bcV[0].face[Y0] = NONE;
  ctx.bcU[0].face[Y1] = NONE;ctx.bcU[1].face[Y1] = ZERO;ctx.bcU[2].face[Y1] = NONE;ctx.bcV[0].face[Y1] = NONE;
  ctx.bcU[0].face[Z0] = NONE;ctx.bcU[1].face[Z0] = NONE;ctx.bcU[2].face[Z0] = ZERO;ctx.bcV[0].face[Z0] = NONE;
  ctx.bcU[0].face[Z1] = NONE;ctx.bcU[1].face[Z1] = NONE;ctx.bcU[2].face[Z1] = ZERO;ctx.bcV[0].face[Z1] = NONE;

  switch (orientation) {
  case 0:
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a transverse rectangular crack of length %g parallel to faces Y0/Y1  along <0,0,1>\n",
                       length);CHKERRQ(ierr);
    for (k = zs; k < zs+zm; k++) {
      for (j = ys; j < ys+ym; j++) {
        for (i = xs; i < xs+xm; i++) {
          if (((j == ny/2) || (j == ny/2-1)) && (PetscAbs(coords_array[k][j][i][0]-(BBmin[0]+BBmax[0])/2.) < length)) {
            v_array[k][j][i] = 0.;
          }
        }
      }
    }
    break;

  default:
    SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Orientation should be between 0 and 2, got %i\n",orientation);
    break;
  }

  ierr = DMDAVecRestoreArray(ctx.daScal,fields.VIrrev,&v_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);


  /*
    Initial Computation
  */

  ierr = VecCopy(fields.VIrrev,fields.V);CHKERRQ(ierr);
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
  ierr = VF_StepV(&fields,&ctx);

  ctx.hasCrackPressure = PETSC_TRUE;
  ierr                 = VF_StepU(&fields,&ctx);

  ctx.ElasticEnergy = 0;
  ctx.InsituWork    = 0;
  ctx.PressureWork  = 0.;
  ierr              = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,&fields,&ctx);CHKERRQ(ierr);
  ierr              = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);
  ctx.TotalEnergy   = ctx.ElasticEnergy+ctx.SurfaceEnergy-ctx.InsituWork-ctx.PressureWork;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic Energy:            %e\n",ctx.ElasticEnergy);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface energy:            %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);
  if (ctx.hasCrackPressure) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of pressure forces:   %e\n",ctx.PressureWork);CHKERRQ(ierr);
  }
  if (ctx.hasInsitu) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of surface forces:    %e\n",ctx.InsituWork);CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Total energy:              %e\n",ctx.TotalEnergy);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%i   \t%e   \t%e   \t%e   \t%e   \t%e   \t%e\n",ctx.timestep,p,ctx.ElasticEnergy,
                                ctx.InsituWork,ctx.SurfaceEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);

  /*
    Save fields and write statistics
  */
  switch (ctx.fileformat) {
  case FILEFORMAT_HDF5:
    ierr = FieldsH5Write(&ctx,&fields);
    break;

  case FILEFORMAT_BIN:
    ierr = FieldsBinaryWrite(&ctx,&fields);
    break;
  }


  /*
    Start of pressure increase by pinc each step
  */

  for (ctx.timestep = 1; ctx.timestep < ctx.maxtimestep; ctx.timestep++) {
    p    = p+pinc;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n Processing time step %i. The current pressure is %e\n",ctx.timestep,p);CHKERRQ(ierr);
    ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
    ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
    VFFractureTimeStep(&ctx,&fields);CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic Energy:            %e\n",ctx.ElasticEnergy);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface energy:            %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);
    if (ctx.hasInsitu) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of surface forces:    %e\n",ctx.InsituWork);CHKERRQ(ierr);
    }
    if (ctx.hasCrackPressure) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of pressure forces:   %e\n",ctx.PressureWork);CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Total energy:              %e\n",ctx.TotalEnergy);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%i   \t%e   \t%e   \t%e   \t%e   \t%e   \t%e\n",ctx.timestep,p,ctx.ElasticEnergy,
                                  ctx.InsituWork,ctx.SurfaceEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);

  }
/* end of time step */

  switch (ctx.fileformat) {
  case FILEFORMAT_HDF5:
    ierr = FieldsH5Write(&ctx,&fields);
    break;

  case FILEFORMAT_BIN:
    ierr = FieldsBinaryWrite(&ctx,&fields);
    break;
  }

  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}
 

