/*
   VFFlow_Fake.c
   A fake fluid flow returning the solution of transient flow between 2 walls
   
     (c) 2010-2011 Blaise Bourdin, C. Chukwudozie, LSU. bourdin@lsu.edu
*/
#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"

#undef __FUNCT__
#define __FUNCT__ "VFFlow_Fake"
/* 
   Fake flow solver for VF_Chevron.c test
*/

extern PetscErrorCode VFFlow_Fake(VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode ierr;
	PetscReal     x,y,z;
	PetscReal			radius_x,radius_y,radius_z;
	PetscReal			dist;
	PetscReal			incremnt = 3.;
	PetscInt			xs,xm,nx;
	PetscInt			ys,ym,ny;
	PetscInt			zs,zm,nz;
	PetscInt			ei,ej,ek;
	PetscReal     p;
	PetscReal			***pressure_array;
	PetscReal			****coords_array;
	
	
	/*
	  Fake flow where the pressure in the reservoir is 
	  constant in an ellipsoid, growing with time/
	*/
	
	
	PetscFunctionBegin;
	if (ctx->verbose >0) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "Entering flow solver %s implemented in %s\n",__FUNCT__,__FILE__);CHKERRQ(ierr);
	}
	
	ierr = DAGetInfo(ctx->daVect,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DAGetCorners(ctx->daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	
	
	ierr = DAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DAVecGetArray(ctx->daScal,fields->pressure,&pressure_array);CHKERRQ(ierr);    
	
	
  radius_x = (ctx->BoundingBox[1] - ctx->BoundingBox[0]) / 2. * (ctx->timevalue / ctx->maxtimevalue + .1);
  radius_y = (ctx->BoundingBox[3] - ctx->BoundingBox[2]) / 20.;
  radius_z = (ctx->BoundingBox[5] - ctx->BoundingBox[4]) / 4.;
	for (ek = zs; ek < zs+zm; ek++) {
    for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
    		z = (coords_array[ek][ej][ei][2] - (ctx->BoundingBox[4] + ctx->BoundingBox[5]) * .5) / radius_z;
        y = (coords_array[ek][ej][ei][1] - ctx->BoundingBox[2]) / radius_y;
				x = (coords_array[ek][ej][ei][0] - ctx->BoundingBox[0]) / radius_x;
				//ierr = PetscPrintf(PETSC_COMM_WORLD,"(%i,%i,%i): [%g,%g,%g], r=%g\n",ei,ej,ek,x,y,z,1. - x*x + y*y + z*z);
				if ( x*x + y*y + z*z < 1.) {
				  p = ctx->resprop.Pinit;
				} else {
				  p = ctx->resprop.Pinit * PetscExpScalar ( (1. - x*x - y*y - z*z) / 2. / ctx->vfprop.epsilon);
				}
				pressure_array[ek][ej][ei] = p;
			}
		}
	}
	
	ierr = DAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DAVecRestoreArray(ctx->daScal,fields->pressure,&pressure_array);CHKERRQ(ierr);

	PetscFunctionReturn(0);	
}
