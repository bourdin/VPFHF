/*
   VFFlow_Fake.c
   A fake fluid flow returning the solution of transient flow between 2 walls
   
     (c) 2010-2011 Blaise Bourdin, LSU. bourdin@lsu.edu
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
	PetscReal			time;
	PetscReal			radius_x, radius_y, alphax = 0.5;
	PetscReal			dist;
	PetscReal			incremnt = 3.;
	PetscInt			xs,xm,nx;
	PetscInt			ys,ym,ny;
	PetscInt			zs,zm,nz;
	PetscInt			ei,ej,ek;
	PetscReal			***pressure_array;
	PetscReal			****coords_array;
	PetscReal			Pinit;
	
	

	
	/*Fake flow where the pressure in the reservoir has an elliptic profile (profile is independent of z axis) 
	 so that the principal radii (min_radius and max_radius) of the profile are a functions of time*/
	
	
	PetscFunctionBegin;
	if (ctx->verbose >0) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "Entering flow solver %s implemented in %s\n",__FUNCT__,__FILE__);CHKERRQ(ierr);
	}
	
	//printf("/n................This is fake flow................../n");
	
	ierr = DAGetInfo(ctx->daVect,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DAGetCorners(ctx->daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	
	
	ierr = DAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DAVecGetArray(ctx->daScal,fields->pressure,&pressure_array);CHKERRQ(ierr);    
	
	
	time = ctx->timevalue;
	Pinit = ctx->resprop.Pinit;
	radius_x = 0.2*(1.0-exp(-time*alphax));
	radius_y = 0.4*(1.0-exp(-time*alphax));
	
	for (ek = zs; ek < zs + zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				dist = pow( ((coords_array[ek][ej][ei][0]-0.5)/radius_x), 2) + pow(((coords_array[ek][ej][ei][1]-0.5)/radius_y), 2);
				if( dist <= 1 )
					pressure_array[ek][ej][ei] = incremnt * Pinit;
				else
					pressure_array[ek][ej][ei] = Pinit;
			}
		}
	}
	
	ierr = DAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DAVecRestoreArray(ctx->daScal,fields->pressure,&pressure_array);CHKERRQ(ierr);

	PetscFunctionReturn(0);	
}
