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
	PetscReal			p;
	PetscReal			***pressure_array;
	PetscReal			****coords_array;
	PetscReal     BBmin[3],BBmax[3];
	
	
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
	  /*
    Get bounding box from petsc DA
  */
  ierr = DAGetBoundingBox(ctx->daVect,BBmin,BBmax);CHKERRQ(ierr);

	ierr = DAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DAVecGetArray(ctx->daScal,fields->pressure,&pressure_array);CHKERRQ(ierr);    
	
  radius_x = (BBmax[0]-BBmin[0]) / 2. * (ctx->timevalue / ctx->maxtimevalue + .1);
  radius_y = (BBmax[1]-BBmin[1]) / 20.;
  radius_z = (BBmax[2]-BBmin[2]) / 4.;
	for (ek = zs; ek < zs+zm; ek++) {
    for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
    		z = (coords_array[ek][ej][ei][2] - (BBmin[2] + BBmax[2]) * .5) / radius_z;
        y = (coords_array[ek][ej][ei][1] - BBmin[1]) / radius_y;
				x = (coords_array[ek][ej][ei][0] - BBmin[0]) / radius_x;
				if ( x*x + y*y + z*z < 1.) {
				  p = 2. * ctx->resprop.Pinit;
				} else {
				  p = ctx->resprop.Pinit * (1. + PetscExpScalar ( (1. - x*x - y*y - z*z) / 2. / ctx->vfprop.epsilon));
				}
				pressure_array[ek][ej][ei] = p;
			}
		}
	}
	ierr = DAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DAVecRestoreArray(ctx->daScal,fields->pressure,&pressure_array);CHKERRQ(ierr);

	PetscFunctionReturn(0);	
}
