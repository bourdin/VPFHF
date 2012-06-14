/*
  test27.c: Check penny shape cracks surface energy computation
  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFV.h"
#include "VFCracks.h"

VFCtx    ctx;
VFFields fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  VFCtx          ctx;
  VFFields       fields;
  PetscErrorCode ierr;

  PetscReal radius      = .2;
  PetscReal center[3]   = {0.,0.,.5};
  PetscInt  nopts       = 3;
  PetscReal ****coords_array;
  PetscReal ***v_array;
  PetscReal BBmin[3],BBmax[3];
  PetscReal InsituWork    = 0;
  PetscReal p = 1e-3;
  char      prefix[PETSC_MAX_PATH_LEN];
  PetscInt  i,j,c;
  
  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);


  ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  ierr = VecSet(fields.VIrrev,1.0);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx.daScal,fields.VIrrev,&v_array);CHKERRQ(ierr);

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

  ctx.numCracks = 0;
  ierr = PetscOptionsInt("-nc","\n\tNumber of penny-shaped cracks to insert","",ctx.numCracks,&ctx.numCracks,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscMalloc(ctx.numCracks*sizeof(VFPennyCrack),&ctx.crack);CHKERRQ(ierr);
  for (i = 0; i < ctx.numCracks; i++) {
    ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN,"c%d_",i);CHKERRQ(ierr);
    ierr = VFPennyCrackCreate(&(ctx.crack[i]));CHKERRQ(ierr);
    ierr = VFPennyCrackGet(prefix, &(ctx.crack[i]));CHKERRQ(ierr);
    if (ctx.verbose > 0) {
      ierr = VFPennyCrackView(&(ctx.crack[i]),PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    }
  }

	ierr = VecSet(fields.V,1.0);CHKERRQ(ierr);
	ierr = VecSet(fields.VIrrev,1.0);CHKERRQ(ierr);
	for (c = 0; c < ctx.numCracks; c++) {
		ierr = VFPennyCrackBuildVAT2(fields.V,&(ctx.crack[c]),&ctx);CHKERRQ(ierr);
		ierr = VecPointwiseMin(fields.VIrrev,fields.V,fields.VIrrev);CHKERRQ(ierr);
	}

  ctx.SurfaceEnergy    = 0.;
  ierr                 = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface energy:              %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);

  ierr = VecCopy(fields.VIrrev,fields.V);CHKERRQ(ierr);
  switch (ctx.fileformat) {
  case FILEFORMAT_HDF5:
    ctx.timestep++;
    ctx.timevalue += 1.;
    ierr = FieldsH5Write(&ctx,&fields);
    break;

  case FILEFORMAT_BIN:
    ierr = FieldsBinaryWrite(&ctx,&fields);
    break;
  }

  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
	ctx.hasCrackPressure = PETSC_FALSE;
  ierr = VF_StepV(&fields,&ctx);

  ctx.SurfaceEnergy    = 0.;
  ierr                 = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface energy:              %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);



  /*
    Save fields and write statistics about current run
  */
  switch (ctx.fileformat) {
  case FILEFORMAT_HDF5:
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

