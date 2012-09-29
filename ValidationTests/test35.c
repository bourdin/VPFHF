/*
 test35.c: Testing for different DA's for cell and node by creating DAComposite
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFV.h"
#include "VFU.h"
#include "VFPermfield.h"

VFCtx               ctx;
VFFields            fields;

typedef struct {
	DM				da_node,da_cell;
	DM				packer;
	PetscViewer		u_viewer,lambda_viewer;
	PetscViewer		fu_viewer,flambda_viewer;
} UserCtx;



#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode  ierr;
	PetscViewer		viewer;
	Vec				U;
	UserCtx			user;
	PetscInt        nx,ny,nz,xs,xm,ys,ym,zs,zm;
	PetscInt		xs1,xm1,ys1,ym1,zs1,zm1;
	PetscInt		xs2,xm2,ys2,ym2,zs2,zm2;	
	
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);

	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	 /* Create a global vector that includes two da arrays */
	ierr = DMCompositeCreate(PETSC_COMM_WORLD,&user.packer);CHKERRQ(ierr);
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,
						DMDA_STENCIL_BOX,nx,ny,nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,3,1,
						PETSC_NULL,PETSC_NULL,PETSC_NULL,&user.da_node);CHKERRQ(ierr);
	ierr = DMCompositeAddDM(user.packer,user.da_node);CHKERRQ(ierr);
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,
						DMDA_STENCIL_BOX,nx-1,ny-1,nz-1,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,3,1,
						PETSC_NULL,PETSC_NULL,PETSC_NULL,&user.da_cell);CHKERRQ(ierr);
	ierr = DMCompositeAddDM(user.packer,user.da_cell);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user.packer,&U);CHKERRQ(ierr);
	
	ierr = DMDAGetCorners(user.da_node,&xs1,&ys1,&zs1,&xm1,&ym1,&zm1);CHKERRQ(ierr);
	ierr = DMDAGetCorners(user.da_cell,&xs2,&ys2,&zs2,&xm2,&ym2,&zm2);CHKERRQ(ierr);
	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Node based DA info are:\n\t xs = %d \t ys = %d \t zs = %d \t xm = %d \t ym = %d \t zm = %d\n", xs1, ys1, zs1, xm1, ym1, zm1);CHKERRQ(ierr);
	ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Cell based DA info are:\n\t xs = %d \t ys = %d \t zs = %d \t xm = %d \t ym = %d \t zm = %d\n", xs2, ys2, zs2, xm2, ym2, zm2);CHKERRQ(ierr);
	ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = DMDestroy(&user.da_node);CHKERRQ(ierr);
	ierr = DMDestroy(&user.da_cell);CHKERRQ(ierr);
	ierr = DMDestroy(&user.packer);CHKERRQ(ierr);
	ierr = VecDestroy(&U);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

