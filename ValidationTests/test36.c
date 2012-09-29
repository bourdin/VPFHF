/*
 test36.c: Testing for different DA's for cell and node by creating cell DA from node DA information
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
	UserCtx			user;
	PetscInt        nx,ny,nz,xs,xm,ys,ym,zs,zm;
	PetscInt		xs2,xm2,ys2,ym2,zs2,zm2;	
	const PetscInt	*lx,*ly,*lz;
	PetscInt		x_nprocs,y_nprocs,z_nprocs,*olx,*oly,*olz;
	
	
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,&x_nprocs,&y_nprocs,&z_nprocs,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	
	ierr = DMDAGetOwnershipRanges(ctx.daScal,&lx,&ly,&lz);CHKERRQ(ierr);
	ierr = PetscMalloc(x_nprocs*sizeof(*olx),&olx);CHKERRQ(ierr);
	ierr = PetscMalloc(y_nprocs*sizeof(*oly),&oly);CHKERRQ(ierr);
	ierr = PetscMalloc(z_nprocs*sizeof(*olz),&olz);CHKERRQ(ierr);

	
	ierr = PetscMemcpy(olx,lx,x_nprocs*sizeof(*olx));CHKERRQ(ierr);
	ierr = PetscMemcpy(oly,ly,y_nprocs*sizeof(*oly));CHKERRQ(ierr);
	ierr = PetscMemcpy(olz,lz,z_nprocs*sizeof(*olz));CHKERRQ(ierr);

	olx[x_nprocs-1]--;
	oly[y_nprocs-1]--;
	olz[z_nprocs-1]--;
	
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,
						DMDA_STENCIL_BOX,nx-1,ny-1,nz-1,x_nprocs,y_nprocs,z_nprocs,3,1,
						olx,oly,olz,&user.da_cell);CHKERRQ(ierr);

	ierr = DMDAGetCorners(user.da_cell,&xs2,&ys2,&zs2,&xm2,&ym2,&zm2);CHKERRQ(ierr);
	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "node based DA info are:\n\t xs = %d \t ys = %d \t zs = %d \t xm = %d \t ym = %d \t zm = %d\n", xs, ys, zs, xm, ym, zm);CHKERRQ(ierr);
	ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);

	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Cell based DA info are:\n\t xs = %d \t ys = %d \t zs = %d \t xm = %d \t ym = %d \t zm = %d\n", xs2, ys2, zs2, xm2, ym2, zm2);CHKERRQ(ierr);
	ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
	
	
	ierr = DMDestroy(&user.da_cell);CHKERRQ(ierr);
	
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

