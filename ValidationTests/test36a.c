/*
 *  test36a.c: 2D/3D. Coupled reservoir and fracture flow example.
 *   (c) 2012-2014 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 *
 *    Verifies application of tip removal for a pressurized line fracture
 *
 *
 *      */

#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFMech.h"
#include "VFFlow.h"
#include "VFPermfield.h"


VFCtx               ctx;
VFFields            fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode  ierr;

	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);

  ctx.timestep = 0;
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
  ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,0.);CHKERRQ(ierr);

  ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);

//  Normal width computation without tip effect for translated domain
  ctx.removeTipEffect = PETSC_FALSE;
  ierr = PetscPrintf(PETSC_COMM_WORLD," STEP 0 \n ");CHKERRQ(ierr);
  ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
  ierr = VolumetricCrackOpening(&ctx.CrackVolume, &ctx, &fields);CHKERRQ(ierr);
  ierr = UpdateFractureWidth(&ctx,&fields);CHKERRQ(ierr);
  ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
  

	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}


