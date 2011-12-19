/*
  test1.c: solves for the displacement in a pressurized penny crack in 3d
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
  
  PetscReal           radius = .2;
  PetscReal           center[3]={0.,0.,.5};
  PetscInt            orientation=2;
  PetscInt            nopts=3;
  PetscInt            i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal       ****coords_array;
  PetscReal        ***v_array;  
  PetscReal           BBmin[3],BBmax[3];
  PetscReal           x,y,z;  
  PetscReal           ElasticEnergy = 0;
  PetscReal           InsituWork = 0;
  PetscReal           SurfaceEnergy = 0;
  char                filename[FILENAME_MAX];

  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  
  ierr = PetscOptionsGetReal(PETSC_NULL,"-radius",&radius,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetRealArray(PETSC_NULL,"-center",&center[0],&nopts,PETSC_NULL);CHKERRQ(ierr);

  ierr = PetscOptionsGetInt(PETSC_NULL,"-orientation",&orientation,PETSC_NULL);CHKERRQ(ierr);
  ierr = DAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	
	ierr = DAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VecSet(fields.VIrrev,1.0);CHKERRQ(ierr);
	ierr = DAVecGetArray(ctx.daScal,fields.VIrrev,&v_array);CHKERRQ(ierr);    

  switch (orientation) {
    case 1:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a penny-shaped crack of radius %g at (%g,%g,%g) with normal vector <1,0,0>\n",
                         radius,center[0],center[1],center[2]);CHKERRQ(ierr);
      i = (center[0] - BBmin[0]) / (BBmax[0] - BBmin[0]) * nx;      
      if (i >= xs && i < (xs+xm)) { 
        for (k = zs; k < zs+zm; k++) {
          for (j = ys; j < ys+ym; j++) {
            z = (coords_array[k][j][i][2] - center[2]) / radius;
            y = (coords_array[k][j][i][1] - center[1]) / radius;
            if ( z*z + y*y < 1.) {
              v_array[k][j][i] = 0.;
            }
          }
        }
      }      
      break;
    case 2:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a penny-shaped crack of radius %g at (%g,%g,%g) with normal vector <0,1,0>\n",
                         radius,center[0],center[1],center[2]);CHKERRQ(ierr);
      j = (center[1] - BBmin[1]) / (BBmax[1] - BBmin[1]) * ny;      
      if (j >= ys && j < (ys+ym)) { 
        for (k = zs; k < zs+zm; k++) {
          for (i = xs; i < xs+xm; i++) {
            z = (coords_array[k][j][i][2] - center[2]) / radius;
            x = (coords_array[k][j][i][0] - center[0]) / radius;
            if ( x*x + z*z < 1.) {
              v_array[k][j][i] = 0.;
            }
          }
        }
      }      
      break;
    case 3:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a penny-shaped crack of radius %g at (%g,%g,%g) with normal vector <0,0,1>\n",
                         radius,center[0],center[1],center[2]);CHKERRQ(ierr);
      k = (center[2] - BBmin[2]) / (BBmax[2] - BBmin[2]) * nz;      
      if (k >= zs && k < (zs+zm)) { 
        for (j = ys; j < ys+ym; j++) {
          for (i = xs; i < xs+xm; i++) {
            y = (coords_array[k][j][i][1] - center[1]) / radius;
            x = (coords_array[k][j][i][0] - center[0]) / radius;
            if ( x*x + y*y < 1.) {
              v_array[k][j][i] = 0.;
            }
          }
        }
      }      
      break;
    default:
      SETERRQ1(PETSC_ERR_USER,"ERROR: Orientation should be one of {1,2,3}, got %i\n",orientation);
      break;
  }  


	ierr = DAVecRestoreArray(ctx.daScal,fields.VIrrev,&v_array);CHKERRQ(ierr);
	ierr = DAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);

  ierr = VecCopy(fields.VIrrev,fields.V);CHKERRQ(ierr);
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
  ierr = VF_StepV(&fields,&ctx);
  
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,ctx.resprop.Pinit);CHKERRQ(ierr);
  ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);

  ierr = BCUUpdate(&ctx.bcU[0],ctx.preset);CHKERRQ(ierr);
  ctx.hasCrackPressure = PETSC_TRUE;
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
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Total energy:              %e\n",ctx.ElasticEnergy-InsituWork-ctx.PressureWork);CHKERRQ(ierr);



    /*
      Save fields and write statistics about current run
    */    
    switch (ctx.fileformat) {
      case FILEFORMAT_HDF5:       
        ierr = FieldsH5Write(&ctx,&fields);
        ctx.timestep++;
        ctx.timevalue += 1.;
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
  
