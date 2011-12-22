/*
  test6.c: 
  Validate elasticity solver by applying pure tension boundary condition

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
  
  PetscReal           length = .2;
  PetscInt            orientation=2;
  PetscInt            nopts=3;
  PetscInt            i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal       ****coords_array;
  PetscReal       ****bcu_array;  
  PetscReal           BBmin[3],BBmax[3];
  PetscReal           ElasticEnergy = 0;
  PetscReal           InsituWork = 0;
  PetscReal           SurfaceEnergy = 0;
  char                filename[FILENAME_MAX];
  PetscReal           bc = .5;

  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  
  ierr = PetscOptionsGetInt(PETSC_NULL,"-orientation",&orientation,PETSC_NULL);CHKERRQ(ierr);
  ierr = DAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	

  ctx.matprop[0].beta  = 0.;
  ctx.matprop[0].alpha = 0.;
  
  ctx.timestep  = 1;
  ctx.timevalue = 1.;

  ierr = DAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  ierr = VecSet(fields.V,1.0);CHKERRQ(ierr);
  ierr = VecSet(fields.VIrrev,1.0);CHKERRQ(ierr);
  ierr = VecSet(fields.U,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
  
	ierr = VecSet(fields.BCU,0.0);CHKERRQ(ierr);
	ierr = DAVecGetArrayDOF(ctx.daVect,fields.BCU,&bcu_array);CHKERRQ(ierr);    

  ctx.bcU[0].face[X0]=NONE; ctx.bcU[1].face[X0]=NONE; ctx.bcU[2].face[X0]=NONE; ctx.bcV[0].face[X0]=NONE;
  ctx.bcU[0].face[X1]=NONE; ctx.bcU[1].face[X1]=NONE; ctx.bcU[2].face[X1]=NONE; ctx.bcV[0].face[X1]=NONE;
  ctx.bcU[0].face[Y0]=NONE; ctx.bcU[1].face[Y0]=NONE; ctx.bcU[2].face[Y0]=NONE; ctx.bcV[0].face[Y0]=NONE;
  ctx.bcU[0].face[Y1]=NONE; ctx.bcU[1].face[Y1]=NONE; ctx.bcU[2].face[Y1]=NONE; ctx.bcV[0].face[Y1]=NONE;
  ctx.bcU[0].face[Z0]=NONE; ctx.bcU[1].face[Z0]=NONE; ctx.bcU[2].face[Z0]=NONE; ctx.bcV[0].face[Z0]=NONE;
  ctx.bcU[0].face[Z1]=NONE; ctx.bcU[1].face[Z1]=NONE; ctx.bcU[2].face[Z1]=NONE; ctx.bcV[0].face[Z1]=NONE;
  switch (orientation) {
    case 0:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Applying traction Dirichlet conditions on faces X0 X1\n");CHKERRQ(ierr);     
      ctx.bcU[0].face[X0]=FIXED;        ctx.bcU[0].face[X1]=FIXED;
      for (k = zs; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
          for (i = xs; i < xs+xm; i++) { 
            if (i == 0) {
              bcu_array[k][j][i][0] = -bc;
            }
            if (i == nx-1) {
              bcu_array[k][j][i][0] = bc;
            }
          }
        }
      }      

      break;
    case 1:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Applying traction Dirichlet conditions on faces Y0 Y1\n");CHKERRQ(ierr);     
      ctx.bcU[1].face[Y0]=FIXED;        ctx.bcU[1].face[Y1]=FIXED;
      for (k = zs; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
          for (i = xs; i < xs+xm; i++) { 
            if (j == 0) {
              bcu_array[k][j][i][1] = -bc;
            }
            if (j == ny-1) {
              bcu_array[k][j][i][1] = bc;
            }
          }
        }
      }      
      break;
    case 2:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Applying traction Dirichlet conditions on faces Z0 Z1\n");CHKERRQ(ierr);     
      ctx.bcU[2].face[Z0]=FIXED;        ctx.bcU[2].face[Z1]=FIXED;
      for (k = zs; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
          for (i = xs; i < xs+xm; i++) { 
            if (k == 0) {
              bcu_array[k][j][i][2] = -bc;
            }
            if (k == nz-1) {
              bcu_array[k][j][i][2] = bc;
            }
          }
        }
      }      
      break;
    default:
      SETERRQ1(PETSC_ERR_USER,"ERROR: Orientation should be between 0 and 2, got %i\n",orientation);
      break;
  }  
	ierr = DAVecRestoreArrayDOF(ctx.daVect,fields.BCU,&bcu_array);CHKERRQ(ierr);
  ierr = DAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);

  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);

  ierr = VF_StepU(&fields,&ctx);
  ctx.ElasticEnergy=0;
  ctx.InsituWork=0;
  ctx.PressureWork = 0.;
  ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,&fields,&ctx);CHKERRQ(ierr);
  ctx.TotalEnergy = ctx.ElasticEnergy - ctx.InsituWork - ctx.PressureWork;
    
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
  
