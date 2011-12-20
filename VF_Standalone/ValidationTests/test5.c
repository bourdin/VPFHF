/*
  test5.c: 
  solves for the displacement in brick subject to a pressure force in one of its faces

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
  PetscReal        ***v_array;  
  PetscReal           BBmin[3],BBmax[3];
  PetscReal           ElasticEnergy = 0;
  PetscReal           InsituWork = 0;
  PetscReal           SurfaceEnergy = 0;
  char                filename[FILENAME_MAX];

  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  
  ierr = PetscOptionsGetInt(PETSC_NULL,"-orientation",&orientation,PETSC_NULL);CHKERRQ(ierr);
  ierr = DAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	
	ierr = DAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VecSet(fields.VIrrev,1.0);CHKERRQ(ierr);
	ierr = DAVecGetArray(ctx.daScal,fields.VIrrev,&v_array);CHKERRQ(ierr);    

  switch (orientation) {
    case 0:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Applying pressure force on face X0\n",
                         length);CHKERRQ(ierr);     
      ctx.bcU[0].face[X0]=NONE; ctx.bcU[1].face[X0]=NONE; ctx.bcU[2].face[X0]=NONE; ctx.bcV[0].face[X0]=ZERO;
      ctx.bcU[0].face[X1]=ZERO; ctx.bcU[1].face[X1]=ZERO; ctx.bcU[2].face[X1]=ZERO; ctx.bcV[0].face[X1]=NONE;
      ctx.bcU[0].face[Y0]=NONE; ctx.bcU[1].face[Y0]=NONE; ctx.bcU[2].face[Y0]=NONE; ctx.bcV[0].face[Y0]=NONE;
      ctx.bcU[0].face[Y1]=NONE; ctx.bcU[1].face[Y1]=NONE; ctx.bcU[2].face[Y1]=NONE; ctx.bcV[0].face[Y1]=NONE;
      ctx.bcU[0].face[Z0]=NONE; ctx.bcU[1].face[Z0]=NONE; ctx.bcU[2].face[Z0]=NONE; ctx.bcV[0].face[Z0]=NONE;
      ctx.bcU[0].face[Z1]=NONE; ctx.bcU[1].face[Z1]=NONE; ctx.bcU[2].face[Z1]=NONE; ctx.bcV[0].face[Z1]=NONE;
      break;
    case 1:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Applying pressure force on face X1\n",
                         length);CHKERRQ(ierr);      
      ctx.bcU[0].face[X0]=ZERO; ctx.bcU[1].face[X0]=ZERO; ctx.bcU[2].face[X0]=ZERO; ctx.bcV[0].face[X0]=NONE;
      ctx.bcU[0].face[X1]=NONE; ctx.bcU[1].face[X1]=NONE; ctx.bcU[2].face[X1]=NONE; ctx.bcV[0].face[X1]=ZERO;
      ctx.bcU[0].face[Y0]=NONE; ctx.bcU[1].face[Y0]=NONE; ctx.bcU[2].face[Y0]=NONE; ctx.bcV[0].face[Y0]=NONE;
      ctx.bcU[0].face[Y1]=NONE; ctx.bcU[1].face[Y1]=NONE; ctx.bcU[2].face[Y1]=NONE; ctx.bcV[0].face[Y1]=NONE;
      ctx.bcU[0].face[Z0]=NONE; ctx.bcU[1].face[Z0]=NONE; ctx.bcU[2].face[Z0]=NONE; ctx.bcV[0].face[Z0]=NONE;
      ctx.bcU[0].face[Z1]=NONE; ctx.bcU[1].face[Z1]=NONE; ctx.bcU[2].face[Z1]=NONE; ctx.bcV[0].face[Z1]=NONE;
      break;
    case 2:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Applying pressure force on face Y0\n",
                         length);CHKERRQ(ierr);      
      ctx.bcU[0].face[X0]=NONE; ctx.bcU[1].face[X0]=NONE; ctx.bcU[2].face[X0]=NONE; ctx.bcV[0].face[X0]=NONE; 
      ctx.bcU[0].face[X1]=NONE; ctx.bcU[1].face[X1]=NONE; ctx.bcU[2].face[X1]=NONE; ctx.bcV[0].face[X1]=NONE;
      ctx.bcU[0].face[Y0]=NONE; ctx.bcU[1].face[Y0]=NONE; ctx.bcU[2].face[Y0]=NONE; ctx.bcV[0].face[Y0]=ZERO;
      ctx.bcU[0].face[Y1]=ZERO; ctx.bcU[1].face[Y1]=ZERO; ctx.bcU[2].face[Y1]=ZERO; ctx.bcV[0].face[Y1]=NONE;
      ctx.bcU[0].face[Z0]=NONE; ctx.bcU[1].face[Z0]=NONE; ctx.bcU[2].face[Z0]=NONE; ctx.bcV[0].face[Z0]=NONE;
      ctx.bcU[0].face[Z1]=NONE; ctx.bcU[1].face[Z1]=NONE; ctx.bcU[2].face[Z1]=NONE; ctx.bcV[0].face[Z1]=NONE;
      break;
    case 3:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Applying pressure force on face Y1\n",
                         length);CHKERRQ(ierr);      
      ctx.bcU[0].face[X0]=NONE; ctx.bcU[1].face[X0]=NONE; ctx.bcU[2].face[X0]=NONE; ctx.bcV[0].face[X0]=NONE; 
      ctx.bcU[0].face[X1]=NONE; ctx.bcU[1].face[X1]=NONE; ctx.bcU[2].face[X1]=NONE; ctx.bcV[0].face[X1]=NONE;
      ctx.bcU[0].face[Y0]=ZERO; ctx.bcU[1].face[Y0]=ZERO; ctx.bcU[2].face[Y0]=ZERO; ctx.bcV[0].face[Y0]=NONE;
      ctx.bcU[0].face[Y1]=NONE; ctx.bcU[1].face[Y1]=NONE; ctx.bcU[2].face[Y1]=NONE; ctx.bcV[0].face[Y1]=ZERO;
      ctx.bcU[0].face[Z0]=NONE; ctx.bcU[1].face[Z0]=NONE; ctx.bcU[2].face[Z0]=NONE; ctx.bcV[0].face[Z0]=NONE;
      ctx.bcU[0].face[Z1]=NONE; ctx.bcU[1].face[Z1]=NONE; ctx.bcU[2].face[Z1]=NONE; ctx.bcV[0].face[Z1]=NONE;
      break;
    case 4:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Applying pressure force on face Z0\n",
                         length);CHKERRQ(ierr);      
      ctx.bcU[0].face[X0]=NONE; ctx.bcU[1].face[X0]=NONE; ctx.bcU[2].face[X0]=NONE; ctx.bcV[0].face[X0]=NONE; 
      ctx.bcU[0].face[X1]=NONE; ctx.bcU[1].face[X1]=NONE; ctx.bcU[2].face[X1]=NONE; ctx.bcV[0].face[X1]=NONE;
      ctx.bcU[0].face[Y0]=NONE; ctx.bcU[1].face[Y0]=NONE; ctx.bcU[2].face[Y0]=NONE; ctx.bcV[0].face[Y0]=NONE;
      ctx.bcU[0].face[Y1]=NONE; ctx.bcU[1].face[Y1]=NONE; ctx.bcU[2].face[Y1]=NONE; ctx.bcV[0].face[Y1]=NONE;
      ctx.bcU[0].face[Z0]=NONE; ctx.bcU[1].face[Z0]=NONE; ctx.bcU[2].face[Z0]=NONE; ctx.bcV[0].face[Z0]=ZERO;
      ctx.bcU[0].face[Z1]=ZERO; ctx.bcU[1].face[Z1]=ZERO; ctx.bcU[2].face[Z1]=ZERO; ctx.bcV[0].face[Z1]=NONE;
      break;
    case 5:
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Applying pressure force on face Z1\n",
                         length);CHKERRQ(ierr);      
      ctx.bcU[0].face[X0]=NONE; ctx.bcU[1].face[X0]=NONE; ctx.bcU[2].face[X0]=NONE; ctx.bcV[0].face[X0]=NONE; 
      ctx.bcU[0].face[X1]=NONE; ctx.bcU[1].face[X1]=NONE; ctx.bcU[2].face[X1]=NONE; ctx.bcV[0].face[X1]=NONE;
      ctx.bcU[0].face[Y0]=NONE; ctx.bcU[1].face[Y0]=NONE; ctx.bcU[2].face[Y0]=NONE; ctx.bcV[0].face[Y0]=NONE;
      ctx.bcU[0].face[Y1]=NONE; ctx.bcU[1].face[Y1]=NONE; ctx.bcU[2].face[Y1]=NONE; ctx.bcV[0].face[Y1]=NONE;
      ctx.bcU[0].face[Z0]=ZERO; ctx.bcU[1].face[Z0]=ZERO; ctx.bcU[2].face[Z0]=ZERO; ctx.bcV[0].face[Z0]=NONE;
      ctx.bcU[0].face[Z1]=NONE; ctx.bcU[1].face[Z1]=NONE; ctx.bcU[2].face[Z1]=NONE; ctx.bcV[0].face[Z1]=ZERO;
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
  
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,1.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);

  //ierr = BCUUpdate(&ctx.bcU[0],ctx.preset);CHKERRQ(ierr);
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
  
