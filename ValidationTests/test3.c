/*
test3.c:  Computes the elastic energy associated with boundary displacement
          given by Sneddon for a pressurized crack in a 2d domain

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/

#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFMech.h"
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
  PetscInt  orientation = 2;
  PetscInt  nopts       = 3;
  PetscInt  i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal x,y,z,u;
  PetscReal E,nu,p;
  PetscReal ****coords_array;
  PetscReal ****bcu_array;
  PetscReal BBmin[3],BBmax[3];
  PetscReal ElasticEnergy = 0;
  PetscReal InsituWork    = 0;
  PetscReal SurfaceEnergy = 0;
  char      filename[FILENAME_MAX];

  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);

  if (ctx.nlayer > 1) {
    SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: This example only makes sense for 1 layer, got %i\n",ctx.nlayer);
  }
  E  = ctx.matprop[0].E;
  nu = ctx.matprop[0].nu;
  p  = 1.e-3;

  ierr = PetscOptionsGetReal(PETSC_NULL,"-length",&length,PETSC_NULL);CHKERRQ(ierr);

  ierr = PetscOptionsGetInt(PETSC_NULL,"-orientation",&orientation,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);

  ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  ierr = VecSet(fields.BCU,0.0);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx.daVect,fields.BCU,&bcu_array);CHKERRQ(ierr);

  /*
    Reset all BC for U and V
  */
  for (i = 0; i < 6; i++) {
    ctx.bcV[0].face[i] = NONE;
    for (j = 0; j < 3; j++) {
      ctx.bcU[j].face[i] = NONE;
    }
  }
  for (i = 0; i < 12; i++) {
    ctx.bcV[0].edge[i] = NONE;
    for (j = 0; j < 3; j++) {
      ctx.bcU[j].edge[i] = NONE;
    }
  }
  for (i = 0; i < 8; i++) {
    ctx.bcV[0].vertex[i] = NONE;
    for (j = 0; j < 3; j++) {
      ctx.bcU[j].vertex[i] = NONE;
    }
  }
  switch (orientation) {
  case 0:
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a transverse rectangular crack of length %g with normal vector <1,0,0> along <0,1,0>\n",
                       length);CHKERRQ(ierr);
    ctx.bcU[0].face[X0] = FIXED;ctx.bcU[1].face[X0] = NONE;ctx.bcU[2].face[X0] = NONE;
    ctx.bcU[0].face[X1] = ZERO;ctx.bcU[1].face[X1] = ZERO;ctx.bcU[2].face[X1] = ZERO;
    ctx.bcU[0].face[Y0] = NONE;ctx.bcU[1].face[Y0] = ZERO;ctx.bcU[2].face[Y0] = NONE;
    ctx.bcU[0].face[Y1] = NONE;ctx.bcU[1].face[Y1] = ZERO;ctx.bcU[2].face[Y1] = NONE;
    ctx.bcU[0].face[Z0] = NONE;ctx.bcU[1].face[Z0] = NONE;ctx.bcU[2].face[Z0] = ZERO;
    ctx.bcU[0].face[Z1] = ZERO;ctx.bcU[1].face[Z1] = ZERO;ctx.bcU[2].face[Z1] = ZERO;
    if (xs == 0) {
      i = 0;
      for (k = zs; k < zs+zm; k++) {
        z = coords_array[k][ys][xs][2];
        for (j = ys; j < ys+ym; j++) {
          if (z < length) {
            bcu_array[k][j][i][0] = 2.*(1-nu*nu)/E*p*sqrt(length*length-z*z);
          }
        }
      }
    }
    break;

  case 1:
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a transverse rectangular crack of length %g with normal vector <1,0,0> along <0,0,1>\n",
                       length);CHKERRQ(ierr);
    ctx.bcU[0].face[X0] = FIXED;ctx.bcU[1].face[X0] = NONE;ctx.bcU[2].face[X0] = NONE;
    ctx.bcU[0].face[X1] = ZERO;ctx.bcU[1].face[X1] = ZERO;ctx.bcU[2].face[X1] = ZERO;
    ctx.bcU[0].face[Y0] = NONE;ctx.bcU[1].face[Y0] = ZERO;ctx.bcU[2].face[Y0] = NONE;
    ctx.bcU[0].face[Y1] = NONE;ctx.bcU[1].face[Y1] = NONE;ctx.bcU[2].face[Y1] = NONE;
    ctx.bcU[0].face[Z0] = NONE;ctx.bcU[1].face[Z0] = NONE;ctx.bcU[2].face[Z0] = ZERO;
    ctx.bcU[0].face[Z1] = NONE;ctx.bcU[1].face[Z1] = NONE;ctx.bcU[2].face[Z1] = ZERO;
    if (xs == 0) {
      i = 0;
      for (k = zs; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
          y = coords_array[zs][j][xs][1];
          if (y < length) {
            bcu_array[k][j][i][0] =  2.*(1-nu*nu)/E*p*sqrt(length*length-y*y);
          }
        }
      }
    }
    break;

  case 2:
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a transverse rectangular crack of length %g with normal vector <0,1,0> along <1,0,0>\n",
                       length);CHKERRQ(ierr);
    ctx.bcU[0].face[X0] = ZERO;ctx.bcU[1].face[X0] = NONE;ctx.bcU[2].face[X0] = NONE;
    ctx.bcU[0].face[X1] = ZERO;ctx.bcU[1].face[X1] = NONE;ctx.bcU[2].face[X1] = NONE;
    ctx.bcU[0].face[Y0] = NONE;ctx.bcU[1].face[Y0] = FIXED;ctx.bcU[2].face[Y0] = NONE;
    ctx.bcU[0].face[Y1] = ZERO;ctx.bcU[1].face[Y1] = ZERO;ctx.bcU[2].face[Y1] = ZERO;
    ctx.bcU[0].face[Z0] = NONE;ctx.bcU[1].face[Z0] = NONE;ctx.bcU[2].face[Z0] = ZERO;
    ctx.bcU[0].face[Z1] = ZERO;ctx.bcU[1].face[Z1] = ZERO;ctx.bcU[2].face[Z1] = ZERO;
    if (ys == 0) {
      j = 0;
      for (i = xs; i < xs+xm; i++) {
        for (k = zs; k < zs+zm; k++) {
          if (coords_array[k][j][i][2] < length) {
            /* bcu_array[k][j][i] = 0.; */
          }
        }
      }
    }
    break;

  case 3:
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a transverse rectangular crack of length %g with normal vector <0,1,0> along <0,0,1>\n",
                       length);CHKERRQ(ierr);
    ctx.bcU[0].face[X0] = ZERO;ctx.bcU[1].face[X0] = NONE;ctx.bcU[2].face[X0] = NONE;
    ctx.bcU[0].face[X1] = ZERO;ctx.bcU[1].face[X1] = ZERO;ctx.bcU[2].face[X1] = ZERO;
    ctx.bcU[0].face[Y0] = NONE;ctx.bcU[1].face[Y0] = FIXED;ctx.bcU[2].face[Y0] = NONE;
    ctx.bcU[0].face[Y1] = ZERO;ctx.bcU[1].face[Y1] = ZERO;ctx.bcU[2].face[Y1] = ZERO;
    ctx.bcU[0].face[Z0] = NONE;ctx.bcU[1].face[Z0] = NONE;ctx.bcU[2].face[Z0] = ZERO;
    ctx.bcU[0].face[Z1] = NONE;ctx.bcU[1].face[Z1] = NONE;ctx.bcU[2].face[Z1] = ZERO;
    if (ys == 0) {
      j = 0;
      for (i = xs; i < xs+xm; i++) {
        for (k = zs; k < zs+zm; k++) {
          if (coords_array[k][j][i][0] < length) {
            /* bcu_array[k][j][i] = 0.; */
          }
        }
      }
    }
    break;

  case 4:
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a transverse rectangular crack of length %g with normal vector <0,0,1> along <1,0,0>\n",
                       length);CHKERRQ(ierr);
    ctx.bcU[0].face[X0] = ZERO;ctx.bcU[1].face[X0] = NONE;ctx.bcU[2].face[X0] = NONE;
    ctx.bcU[0].face[X1] = ZERO;ctx.bcU[1].face[X1] = NONE;ctx.bcU[2].face[X1] = NONE;
    ctx.bcU[0].face[Y0] = NONE;ctx.bcU[1].face[Y0] = ZERO;ctx.bcU[2].face[Y0] = NONE;
    ctx.bcU[0].face[Y1] = ZERO;ctx.bcU[1].face[Y1] = ZERO;ctx.bcU[2].face[Y1] = ZERO;
    ctx.bcU[0].face[Z0] = NONE;ctx.bcU[1].face[Z0] = NONE;ctx.bcU[2].face[Z0] = FIXED;
    ctx.bcU[0].face[Z1] = ZERO;ctx.bcU[1].face[Z1] = ZERO;ctx.bcU[2].face[Z1] = ZERO;
    if (zs == 0) {
      k = 0;
      for (j = ys; j < ys+ym; j++) {
        for (i = xs; i < xs+xm; i++) {
          if (coords_array[k][j][i][1] < length) {
            /* bcu_array[k][j][i] = 0.; */
          }
        }
      }
    }
    break;

  case 5:
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a transverse rectangular crack of length %g with normal vector <0,0,1> along <0,1,0>\n",
                       length);CHKERRQ(ierr);
    ctx.bcU[0].face[X0] = ZERO;ctx.bcU[1].face[X0] = NONE;ctx.bcU[2].face[X0] = NONE;
    ctx.bcU[0].face[X1] = ZERO;ctx.bcU[1].face[X1] = ZERO;ctx.bcU[2].face[X1] = ZERO;
    ctx.bcU[0].face[Y0] = NONE;ctx.bcU[1].face[Y0] = ZERO;ctx.bcU[2].face[Y0] = NONE;
    ctx.bcU[0].face[Y1] = NONE;ctx.bcU[1].face[Y1] = ZERO;ctx.bcU[2].face[Y1] = NONE;
    ctx.bcU[0].face[Z0] = NONE;ctx.bcU[1].face[Z0] = NONE;ctx.bcU[2].face[Z0] = FIXED;
    ctx.bcU[0].face[Z1] = ZERO;ctx.bcU[1].face[Z1] = ZERO;ctx.bcU[2].face[Z1] = ZERO;
    if (zs == 0) {
      k = 0;
      for (j = ys; j < ys+ym; j++) {
        for (i = xs; i < xs+xm; i++) {
          if (coords_array[k][j][i][0] < length) {
            /* bcu_array[k][j][i] = 0.; */
          }
        }
      }
    }
    break;

  default:
    SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Orientation should be between 0 and 5, got %i\n",orientation);
    break;
  }


  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,fields.BCU,&bcu_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);

  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  /* ierr = VecSet(fields.pressure,1);CHKERRQ(ierr); */
  ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);

  /* ierr = BCUUpdate(&ctx.bcU[0],ctx.preset);CHKERRQ(ierr); */
  ctx.hasCrackPressure = PETSC_FALSE;
  ierr                 = VF_StepU(&fields,&ctx);
  ctx.ElasticEnergy    = 0;
  ctx.InsituWork       = 0;
  ctx.PressureWork     = 0.;
  ierr                 = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,&fields,&ctx);CHKERRQ(ierr);
  ctx.TotalEnergy      = ctx.ElasticEnergy-ctx.InsituWork-ctx.PressureWork;



  ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic Energy:            %e\n",ctx.ElasticEnergy);CHKERRQ(ierr);
  if (ctx.hasCrackPressure) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of pressure forces:   %e\n",ctx.PressureWork);CHKERRQ(ierr);
  }
  if (ctx.hasInsitu) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of surface forces:    %e\n",ctx.InsituWork);CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Total energy:              %e\n",ctx.ElasticEnergy-InsituWork-ctx.PressureWork);CHKERRQ(ierr);



  /*
    Save fields and write statistics about current run
  */
  switch (ctx.fileformat) {
  case FILEFORMAT_HDF5:
    ierr = FieldsH5Write(&ctx,&fields);
    ctx.timestep++;
    ctx.timevalue += 1.;
    ierr           = FieldsH5Write(&ctx,&fields);
    break;

  case FILEFORMAT_BIN:
    ierr = FieldsBinaryWrite(&ctx,&fields);
    break;
  }

  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}

