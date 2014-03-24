/*
 test47.c: 2D Heat problem, advection and diffusion with no heat source term. All temperature boundary condition implemented by specifying the analytical pressure values on the boundaries
 (c) 2010-2013 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 Taken from
 http://www.cs.uky.edu/~jzhang/pub/PAPER/adi4th.pdf   

 ./test47 -n 51,51,2 -l 2,2,0.01  -theta 1  -timestepsize 0.01 -maxtimestep 50 -condx 0.01 -condy 0.01 -velx 0.8 -vely 0.8
 
 ./test47 -n 101,101,2 -l 2,2,0.01  -theta 1  -timestepsize 0.001 -maxtimestep 200 -condx 0.01 -condy 0.01 -velx 01 -vely 01

 */

#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFMech.h"
#include "VFFlow.h"
#include "VFHeat.h"

VFCtx               ctx;
VFFields            fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{	
	PetscErrorCode  ierr;
	PetscViewer		viewer;
	PetscViewer     logviewer;
	char			filename[FILENAME_MAX];
	PetscInt		i,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm,xs1,xm1,ys1,ym1,zs1,zm1;
	PetscReal		BBmin[3],BBmax[3];
	PetscReal		****coords_array;
	PetscReal		hx,hy,hz;
	PetscReal		lx,ly,lz;
	PetscReal		***heatsrc_array;
	PetscReal		****heatfluxbc_array;
	PetscReal		****cond_array;
	PetscReal		****vel_array;
  PetscReal		***Tbc_array;
  PetscReal		***T_array;
	PetscReal		ux = 0, uy = 0,kx = 0, ky = 0.0, src = 0;
	
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ctx.flowsolver = FLOWSOLVER_SNESMIXEDFEM;
	ctx.flowsolver = HEATSOLVER_SNESFEM;
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	ierr = VecSet(fields.velocity,0.);CHKERRQ(ierr);
	ierr = VecSet(fields.theta,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.Cond,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.HeatSource,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.TBCArray,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.prevT,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.PresBCArray,0.);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);	
	ierr = DMDAVecGetArray(ctx.daScal,ctx.HeatSource,&heatsrc_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.HeatFluxBCArray,&heatfluxbc_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVFperm,ctx.Cond,&cond_array);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daVFperm,&xs1,&ys1,&zs1,&xm1,&ym1,&zm1);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVect,fields.velocity,&vel_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx.daScal,fields.theta,&T_array);CHKERRQ(ierr);

	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-condx",&kx,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-condy",&ky,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-velx",&ux,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-vely",&uy,PETSC_NULL);CHKERRQ(ierr);
	ctx.hasHeatSources = PETSC_FALSE;
	ctx.flowprop.theta = 1.;
	ctx.flowprop.timestepsize = 1.;
	lz = BBmax[2]-BBmin[2];
	ly = BBmax[1]-BBmin[1];
	lx = BBmax[0]-BBmin[0];
	hx = lx/(nx-1);
	hy = ly/(nx-1);
	hz = lz/(nz-1);
  ctx.matprop[0].rho = 0;
	ctx.matprop[0].Cp = 0;
	for (k = zs1; k < zs1+zm1; k++) {
		for (j = ys1; j < ys1+ym1; j++) {
			for (i = xs1; i < xs1+xm1; i++) {
				cond_array[k][j][i][0] = kx;
				cond_array[k][j][i][1] = ky;
			}
		}
	}
	for (i = 0; i < 6; i++) {
		ctx.bcT[0].face[i] = NONE;
    for(j = 0; j < 3; j++){
      ctx.bcQT[j].face[i] = NONE;
    }
	}
	for (i = 0; i < 12; i++) {
		ctx.bcT[0].edge[i] = NONE;
    for(j = 0; j < 3; j++){
      ctx.bcQT[j].edge[i] = NONE;
    }
	}
	for (i = 0; i < 8; i++) {
		ctx.bcT[0].vertex[i] = NONE;
    for(j = 0; j < 3; j++){
      ctx.bcQT[j].vertex[i] = NONE;
    }
	}
	ctx.bcT[0].face[X0] = FIXED;
	ctx.bcT[0].face[X1] = FIXED;
	ctx.bcT[0].face[Y0] = FIXED;
	ctx.bcT[0].face[Y1] = FIXED;

	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
        T_array[k][j][i] = exp(-(pow((hx*i-0.5),2))/kx-(pow((hy*j-0.5),2))/ky);
				heatsrc_array[k][j][i] = 0;
				vel_array[k][j][i][0] = ux;
				vel_array[k][j][i][1] = uy;
				vel_array[k][j][i][2] = 0.00;
			}
		}
	}
  ierr = DMDAVecRestoreArray(ctx.daScal,fields.theta,&T_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,fields.velocity,&vel_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVFperm,ctx.Cond,&cond_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.HeatFluxBCArray,&heatfluxbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx.daScal,ctx.HeatSource,&heatsrc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-theta",&ctx.flowprop.theta,PETSC_NULL);CHKERRQ(ierr);

	ierr = PetscOptionsGetReal(PETSC_NULL,"-timestepsize",&ctx.flowprop.timestepsize,PETSC_NULL);CHKERRQ(ierr);
	ctx.maxtimestep = 1;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-maxtimestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr);
  
  ctx.timestep = 0;
  ierr = FieldsH5Write(&ctx,&fields);

	for (ctx.timestep = 1; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nProcessing step %i.\n",ctx.timestep);CHKERRQ(ierr);
    
    ierr = DMDAVecGetArray(ctx.daScal,ctx.TBCArray,&Tbc_array);CHKERRQ(ierr);
    for (k = zs; k < zs+zm; k++) {
      for (j = ys; j < ys+ym; j++) {
        for (i = xs; i < xs+xm; i++) {
          Tbc_array[k][j][i] = 1/(4*ctx.timestep*ctx.flowprop.timestepsize+1)*exp(-(pow((hx*i-ux*ctx.timestep*ctx.flowprop.timestepsize-0.5),2))/(kx*(4*ctx.timestep*ctx.flowprop.timestepsize+1))-(pow((hy*j-uy*ctx.timestep*ctx.flowprop.timestepsize-0.5),2))/(ky*(4*ctx.timestep*ctx.flowprop.timestepsize+1)));
        }
      }
    }
    ierr = DMDAVecRestoreArray(ctx.daScal,ctx.TBCArray,&Tbc_array);CHKERRQ(ierr);
    ierr = VecCopy(fields.theta,ctx.prevT);CHKERRQ(ierr);
		ierr = VF_HeatTimeStep(&ctx,&fields);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHST,ctx.RHSTpre);CHKERRQ(ierr);
		ierr = FieldsH5Write(&ctx,&fields);
		ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.log",ctx.prefix);CHKERRQ(ierr);
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&logviewer);CHKERRQ(ierr);
		ierr = PetscLogView(logviewer);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&logviewer);
	}
  ctx.timestep++;
  ierr = VecCopy(ctx.TBCArray,fields.theta);CHKERRQ(ierr);
  ierr = FieldsH5Write(&ctx,&fields);

	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

