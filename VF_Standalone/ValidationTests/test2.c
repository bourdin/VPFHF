/*
  test2.c: 
    solve for the displacement in a pressurized rectangular crack in 3d
    uses symmetry to solve on 1/4 domain

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFV.h"
#include "VFU.h"
#include "VFFlow.h"

VFCtx               ctx;
VFFields            fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  VFCtx               ctx;
  VFFields            fields;
  PetscErrorCode      ierr;
  
  PetscReal           length = .4;
  PetscInt            orientation=2;
  PetscInt            nopts=3;
  PetscInt            i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal       ****coords_array;
  PetscReal        ***v_array;  
  PetscReal           BBmin[3],BBmax[3];
  PetscReal           ElasticEnergy = 0;
  PetscReal           InsituWork = 0;
  PetscReal           SurfaceEnergy = 0;
  char                filename[FILENAME_MAX];
  PetscReal           p = 1.;

  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-length",&length,PETSC_NULL);CHKERRQ(ierr);

  ierr = PetscOptionsGetInt(PETSC_NULL,"-orientation",&orientation,PETSC_NULL);CHKERRQ(ierr);
  ierr = DAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	
  ctx.matprop[0].E     = 1.;
  ctx.matprop[0].nu    = 0.;
  ctx.matprop[0].beta  = 0.;
  ctx.matprop[0].alpha = 0.;
  
  ctx.timestep = 0;
  ctx.timevalue = 1.;

  ierr = DAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  ierr = VecSet(fields.V,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.VIrrev,1.0);CHKERRQ(ierr);
  ierr = VecSet(fields.U,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
  ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
  
  ierr = DAVecGetArray(ctx.daScal,fields.VIrrev,&v_array);CHKERRQ(ierr);    

  switch (orientation) {
    case 0:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a transverse rectangular crack of length %g on face X0 along <0,1,0>\n",
                         length);CHKERRQ(ierr);     
      ctx.bcU[0].face[X0]=ZERO; ctx.bcU[1].face[X0]=NONE; ctx.bcU[2].face[X0]=NONE; ctx.bcV[0].face[X0]=NONE;
      ctx.bcU[0].face[X1]=NONE; ctx.bcU[1].face[X1]=NONE; ctx.bcU[2].face[X1]=NONE; ctx.bcV[0].face[X1]=NONE;
      ctx.bcU[0].face[Y0]=NONE; ctx.bcU[1].face[Y0]=ZERO; ctx.bcU[2].face[Y0]=NONE; ctx.bcV[0].face[Y0]=NONE;
      ctx.bcU[0].face[Y1]=NONE; ctx.bcU[1].face[Y1]=ZERO; ctx.bcU[2].face[Y1]=NONE; ctx.bcV[0].face[Y1]=NONE;
      ctx.bcU[0].face[Z0]=NONE; ctx.bcU[1].face[Z0]=NONE; ctx.bcU[2].face[Z0]=ZERO; ctx.bcV[0].face[Z0]=NONE;
      ctx.bcU[0].face[Z1]=NONE; ctx.bcU[1].face[Z1]=NONE; ctx.bcU[2].face[Z1]=NONE; ctx.bcV[0].face[Z1]=NONE;
      if (xs == 0) { 
        i = 0;
        for (k = zs; k < zs+zm; k++) {
          for (j = ys; j < ys+ym; j++) {
            if ( coords_array[k][j][i][2] < length) {
              v_array[k][j][i] = 0.;
            }
          }
        }
      }      
      break;
    case 1:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a transverse rectangular crack of length %g on face X0 along <0,0,1>\n",
                         length);CHKERRQ(ierr);      
      ctx.bcU[0].face[X0]=ZERO; ctx.bcU[1].face[X0]=NONE; ctx.bcU[2].face[X0]=NONE; ctx.bcV[0].face[X0]=NONE;
      ctx.bcU[0].face[X1]=NONE; ctx.bcU[1].face[X1]=NONE; ctx.bcU[2].face[X1]=NONE; ctx.bcV[0].face[X1]=NONE;
      ctx.bcU[0].face[Y0]=NONE; ctx.bcU[1].face[Y0]=ZERO; ctx.bcU[2].face[Y0]=NONE; ctx.bcV[0].face[Y0]=NONE;
      ctx.bcU[0].face[Y1]=NONE; ctx.bcU[1].face[Y1]=NONE; ctx.bcU[2].face[Y1]=NONE; ctx.bcV[0].face[Y1]=NONE;
      ctx.bcU[0].face[Z0]=NONE; ctx.bcU[1].face[Z0]=NONE; ctx.bcU[2].face[Z0]=ZERO; ctx.bcV[0].face[Z0]=NONE;
      ctx.bcU[0].face[Z1]=NONE; ctx.bcU[1].face[Z1]=NONE; ctx.bcU[2].face[Z1]=ZERO; ctx.bcV[0].face[Z1]=NONE;
      if (xs == 0) { 
        i = 0;
        for (k = zs; k < zs+zm; k++) {
          for (j = ys; j < ys+ym; j++) {
            if ( coords_array[k][j][i][1] < length) {
              v_array[k][j][i] = 0.;
            }
          }
        }
      }      
      break;
    case 2:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a transverse rectangular crack of length %g on face Y0 along <1,0,0>\n",
                         length);CHKERRQ(ierr);      
      ctx.bcU[0].face[X0]=ZERO; ctx.bcU[1].face[X0]=NONE; ctx.bcU[2].face[X0]=NONE; ctx.bcV[0].face[X0]=NONE;
      ctx.bcU[0].face[X1]=ZERO; ctx.bcU[1].face[X1]=NONE; ctx.bcU[2].face[X1]=NONE; ctx.bcV[0].face[X1]=NONE;
      ctx.bcU[0].face[Y0]=NONE; ctx.bcU[1].face[Y0]=ZERO; ctx.bcU[2].face[Y0]=NONE; ctx.bcV[0].face[Y0]=NONE;
      ctx.bcU[0].face[Y1]=NONE; ctx.bcU[1].face[Y1]=NONE; ctx.bcU[2].face[Y1]=NONE; ctx.bcV[0].face[Y1]=NONE;
      ctx.bcU[0].face[Z0]=NONE; ctx.bcU[1].face[Z0]=NONE; ctx.bcU[2].face[Z0]=ZERO; ctx.bcV[0].face[Z0]=NONE;
      ctx.bcU[0].face[Z1]=NONE; ctx.bcU[1].face[Z1]=NONE; ctx.bcU[2].face[Z1]=NONE; ctx.bcV[0].face[Z1]=NONE;
      if (ys == 0) { 
        j = 0;
        for (i = xs; i < xs+xm; i++) {
          for (k = zs; k < zs+zm; k++) {
            if ( coords_array[k][j][i][2] < length) {
              v_array[k][j][i] = 0.;
            }
          }
        }
      }      
      break;
    case 3:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a transverse rectangular crack of length %g on face Y1 along <0,0,1>\n",
                         length);CHKERRQ(ierr);      
      ctx.bcU[0].face[X0]=ZERO; ctx.bcU[1].face[X0]=NONE; ctx.bcU[2].face[X0]=NONE; ctx.bcV[0].face[X0]=NONE;
      ctx.bcU[0].face[X1]=ZERO; ctx.bcU[1].face[X1]=NONE; ctx.bcU[2].face[X1]=NONE; ctx.bcV[0].face[X1]=NONE;
      ctx.bcU[0].face[Y0]=NONE; ctx.bcU[1].face[Y0]=ZERO; ctx.bcU[2].face[Y0]=NONE; ctx.bcV[0].face[Y0]=NONE;
      ctx.bcU[0].face[Y1]=NONE; ctx.bcU[1].face[Y1]=ZERO; ctx.bcU[2].face[Y1]=NONE; ctx.bcV[0].face[Y1]=NONE;
      ctx.bcU[0].face[Z0]=NONE; ctx.bcU[1].face[Z0]=NONE; ctx.bcU[2].face[Z0]=ZERO; ctx.bcV[0].face[Z0]=NONE;
      ctx.bcU[0].face[Z1]=NONE; ctx.bcU[1].face[Z1]=NONE; ctx.bcU[2].face[Z1]=ZERO; ctx.bcV[0].face[Z1]=NONE;
      if (ys ==0) { 
        j = 0;
        for (i = xs; i < xs+xm; i++) {
          for (k = zs; k < zs+zm; k++) {
            if ( coords_array[k][j][i][0] < length) {
              v_array[k][j][i] = 0.;
            }
          }
        }
      }      
      break;
    case 4:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a transverse rectangular crack of length %g on face Z0 along <1,0,0>\n",
                         length);CHKERRQ(ierr);      
      ctx.bcU[0].face[X0]=ZERO; ctx.bcU[1].face[X0]=NONE; ctx.bcU[2].face[X0]=NONE; ctx.bcV[0].face[X0]=NONE;
      ctx.bcU[0].face[X1]=NONE; ctx.bcU[1].face[X1]=NONE; ctx.bcU[2].face[X1]=NONE; ctx.bcV[0].face[X1]=NONE;
      ctx.bcU[0].face[Y0]=NONE; ctx.bcU[1].face[Y0]=ZERO; ctx.bcU[2].face[Y0]=NONE; ctx.bcV[0].face[Y0]=NONE;
      ctx.bcU[0].face[Y1]=NONE; ctx.bcU[1].face[Y1]=ZERO; ctx.bcU[2].face[Y1]=NONE; ctx.bcV[0].face[Y1]=NONE;
      ctx.bcU[0].face[Z0]=NONE; ctx.bcU[1].face[Z0]=NONE; ctx.bcU[2].face[Z0]=ZERO; ctx.bcV[0].face[Z0]=NONE;
      ctx.bcU[0].face[Z1]=NONE; ctx.bcU[1].face[Z1]=NONE; ctx.bcU[2].face[Z1]=NONE; ctx.bcV[0].face[Z1]=NONE;
      if (zs == 0) { 
        k = 0;
        for (j = ys; j < ys+ym; j++) {
          for (i = xs; i < xs+xm; i++) {
            if ( coords_array[k][j][i][1] < length) {
              v_array[k][j][i] = 0.;
            }
          }
        }
      }      
      break;
    case 5:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a transverse rectangular crack of length %g on face Z1 along <0,1,0>\n",
                         length);CHKERRQ(ierr);      
      ctx.bcU[0].face[X0]=ZERO; ctx.bcU[1].face[X0]=NONE; ctx.bcU[2].face[X0]=NONE; ctx.bcV[0].face[X0]=NONE;
      ctx.bcU[0].face[X1]=ZERO; ctx.bcU[1].face[X1]=NONE; ctx.bcU[2].face[X1]=NONE; ctx.bcV[0].face[X1]=NONE;
      ctx.bcU[0].face[Y0]=NONE; ctx.bcU[1].face[Y0]=ZERO; ctx.bcU[2].face[Y0]=NONE; ctx.bcV[0].face[Y0]=NONE;
      ctx.bcU[0].face[Y1]=NONE; ctx.bcU[1].face[Y1]=NONE; ctx.bcU[2].face[Y1]=NONE; ctx.bcV[0].face[Y1]=NONE;
      ctx.bcU[0].face[Z0]=NONE; ctx.bcU[1].face[Z0]=NONE; ctx.bcU[2].face[Z0]=ZERO; ctx.bcV[0].face[Z0]=NONE;
      ctx.bcU[0].face[Z1]=NONE; ctx.bcU[1].face[Z1]=NONE; ctx.bcU[2].face[Z1]=NONE; ctx.bcV[0].face[Z1]=NONE;
      if (zs == 0) { 
        k = 0;
        for (j = ys; j < ys+ym; j++) {
          for (i = xs; i < xs+xm; i++) {
            if ( coords_array[k][j][i][0] < length) {
              v_array[k][j][i] = 0.;
            }
          }
        }
      }    
      break;  
    default:
      SETERRQ1(PETSC_ERR_USER,"ERROR: Orientation should be between 0 and 5, got %i\n",orientation);
      break;
  }  
  ierr = DAVecRestoreArray(ctx.daScal,fields.VIrrev,&v_array);CHKERRQ(ierr);
  ierr = DAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);

  ierr = VecCopy(fields.VIrrev,fields.V);CHKERRQ(ierr);

  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
  ierr = VF_StepV(&fields,&ctx);
  
  ctx.hasCrackPressure = PETSC_TRUE;
  ierr = VF_StepU(&fields,&ctx);
  ctx.ElasticEnergy=0;
  ctx.InsituWork=0;
  ctx.PressureWork = 0.;
  ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,&fields,&ctx);CHKERRQ(ierr);
  ierr = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);
  ctx.TotalEnergy = ctx.ElasticEnergy - ctx.InsituWork - ctx.PressureWork;
    
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface energy:            %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic Energy:            %e\n",ctx.ElasticEnergy);CHKERRQ(ierr);
  if (ctx.hasCrackPressure) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of pressure forces:   %e\n",ctx.PressureWork);CHKERRQ(ierr);
  }
  if (ctx.hasInsitu) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of surface forces:    %e\n",ctx.InsituWork);CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Total Mechanical energy:  %e\n",ctx.ElasticEnergy-InsituWork-ctx.PressureWork);CHKERRQ(ierr);

  /*
    Save fields and write statistics about current run
  */    
  switch (ctx.fileformat) {
    case FILEFORMAT_HDF5:       
      ierr = FieldsH5Write(&ctx,&fields);
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
  
