/*
   VFFlow_Fake.c
   A fake fluid flow returning the solution of transient flow between 2 walls
   
     (c) 2010-2012 Blaise Bourdin, C. Chukwudozie, LSU. bourdin@lsu.edu
*/
#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"

#undef __FUNCT__
#define __FUNCT__ "VFFlow_Fake"
/* 
   Fake flow solver for VF_Chevron.c test
*/

extern PetscErrorCode VFFlow_Fake(VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode ierr;
	PetscInt			xs,xm,nx;
	PetscInt			ys,ym,ny;
	PetscInt			zs,zm,nz;
	PetscInt			ei,ej,ek;
	PetscReal			***pressure_array;
	PetscReal			****vel_array;
	PetscReal			pi,hx,hy,hz;

	
	/*
	  Fake flow where the pressure in the reservoir is 
	  constant in an ellipsoid, growing with time/
	*/
	
	
	PetscFunctionBegin;
	if (ctx->verbose >0) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "Entering flow solver %s implemented in %s\n",__FUNCT__,__FILE__);CHKERRQ(ierr);
	}
	
  ierr = DMDAGetInfo(ctx->daVect,NULL,&nx,&ny,&nz,NULL,NULL,NULL,
                    NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	  /*
    Get bounding box from petsc DA
  */

	ierr = DMDAVecGetArrayDOF(ctx->daVect,fields->velocity,&vel_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,fields->pressure,&pressure_array);CHKERRQ(ierr);    
	pi = 6.*asin(0.5);
	hx = 1./(nx-1.);
	hy = 1./(ny-1.);
	hz = 1./(nz-1.);
	
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				pressure_array[ek][ej][ei] = cos(pi*ei*hx)*cos(pi*ej*hy)*cos(pi*ek*hz)/(3.*pi*pi);
				vel_array[ek][ej][ei][0]=sin(pi*ei*hx)*cos(pi*ej*hy)*cos(pi*ek*hz)/(3.*pi);
				vel_array[ek][ej][ei][1]=cos(pi*ei*hx)*sin(pi*ej*hy)*cos(pi*ek*hz)/(3.*pi);
				vel_array[ek][ej][ei][2]=cos(pi*ei*hx)*cos(pi*ej*hy)*sin(pi*ek*hz)/(3.*pi);
			}
		}
	}
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,fields->velocity,&vel_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,fields->pressure,&pressure_array);CHKERRQ(ierr);

	PetscFunctionReturn(0);	
}
