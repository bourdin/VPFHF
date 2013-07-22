/*
 test36.c: 2D. Single well problem. Source term implemented as dirac function and analytical solution is 3D Green's function. All pressure boundary condition.
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
./test36  -n 101,101,2 -l 1,1,1 -maxtimestep 1 timestepsize 10 -theta 1 -nw 1 -w0_coords 0.5,0.5,0.5 -w0_Qw 1 -w0_constraint Rate -w0_rw 0.1 -w0_type injector -m_inv 0 -flowsolver FLOWSOLVER_SNESMIXEDFEM
 
 ./test36  -n 101,101,2 -l 1,1,0.01 -nw 1 -w0_coords 0.5,0.5,0.005 -w0_Qw 1 -w0_constraint Rate -w0_rw 0.1 -w0_type injector -m_inv 0 -flowsolver FLOWSOLVER_SNESMIXEDFEM -Gc 0.6667 -npc 1 -pc0_r 0.2 -pc0_center 0.5,0.5,0.005 -pc0_thickness 0.05 -pc0_theta 0 -pc0_phi 90 -epsilon 0.02
 
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFMech.h"
#include "VFFlow.h"

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
	PetscInt		i,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm;
	PetscReal		BBmin[3],BBmax[3];
	PetscReal		***presbc_array;
	PetscReal		***src_array;
	PetscReal		****coords_array;
	PetscReal		hx,hy,hz;
	PetscReal		gx,gy,gz;
	PetscReal		lx,ly,lz;
	PetscReal		gamma, beta, rho, mu;
	PetscReal		pi,dist;
  PetscReal		***pre_array;
	PetscReal   ****perm_array;
	PetscInt    xs1,xm1,ys1,ym1,zs1,zm1;
  Vec         V_hold;

		
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ctx.flowsolver = FLOWSOLVER_KSPMIXEDFEM;
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx.daScal,&V_hold);CHKERRQ(ierr);
	ierr = VecSet(ctx.PresBCArray,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.VelBCArray,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.Source,0.);CHKERRQ(ierr);
	ctx.hasFlowWells = PETSC_TRUE;
	ctx.hasFluidSources = PETSC_FALSE;
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx.daScal,ctx.Source,&src_array);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daVFperm,&xs1,&ys1,&zs1,&xm1,&ym1,&zm1);CHKERRQ(ierr);	
	ierr = DMDAVecGetArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr); 
  for (k = zs1; k < zs1+zm1; k++) {
    for (j = ys1; j < ys1+ym1; j++) {
      for (i = xs1; i < xs1+xm1; i++) {
        perm_array[k][j][i][0] = 1;
        perm_array[k][j][i][1] = 1;
        perm_array[k][j][i][2] = 0.;
        perm_array[k][j][i][3] = 0.;
        perm_array[k][j][i][4] = 0.;
        perm_array[k][j][i][5] = 0.;
      }
    }
  }
  for (k = zs1; k < zs1+zm1; k++) {
    for (j = ys1; j < ys1+ym1; j++) {
      for (i = xs1; i < xs1+xm1; i++) {
        if((j == ny/2 || j == ny/2-1) && ( i > nx/8 && i < 8*nx/10)){
          perm_array[k][j][i][0] = 1.;
          perm_array[k][j][i][1] = 1.;
        }
      }
    }
  }
	pi = 6.*asin(0.5);
	rho = ctx.flowprop.rho;									 
	mu = ctx.flowprop.mu;     
	beta = ctx.flowprop.beta;		
	gamma = ctx.flowprop.gamma;									
	for (i = 0; i < 6; i++) {
		ctx.bcP[0].face[i] = NONE;
		for (c = 0; c < 3; c++) {
			ctx.bcQ[c].face[i] = NONE;
		}
	}
	for (i = 0; i < 12; i++) {
		ctx.bcP[0].edge[i] = NONE;
		for (c = 0; c < 3; c++) {
			ctx.bcQ[c].edge[i] = NONE;
		}
	}
	for (i = 0; i < 8; i++) {
		ctx.bcP[0].vertex[i] = NONE;
		for (c = 0; c < 3; c++) {
			ctx.bcQ[c].vertex[i] = NONE;
		}
	}
	for (i = 0; i < 4; i++) {
		ctx.bcP[0].face[i] = FIXED;	
	}	
  ctx.bcQ[2].face[Z0] = FIXED;
	ctx.bcQ[2].face[Z1] = FIXED;
    for (k = zs; k < zs+zm; k++) {
      for (j = ys; j < ys+ym; j++) {
        for (i = xs; i < xs+xm; i++) {
          if(i == 0){
            dist = sqrt(pow((ctx.well[0].coords[0]-coords_array[k][j][i][0]),2)+pow((ctx.well[0].coords[1]-coords_array[k][j][i][1]),2));
            presbc_array[k][j][i] = 1./(2.*pi)*log(dist);
            presbc_array[k][j][i] = 10.;
        }
          if(i == nx-1){
            dist = sqrt(pow((ctx.well[0].coords[0]-coords_array[k][j][i][0]),2)+pow((ctx.well[0].coords[1]-coords_array[k][j][i][1]),2));
            presbc_array[k][j][i] = 1./(2.*pi)*log(dist);
            presbc_array[k][j][i] = 10.;
          }
          if(j == 0){
            dist = sqrt(pow((ctx.well[0].coords[0]-coords_array[k][j][i][0]),2)+pow((ctx.well[0].coords[1]-coords_array[k][j][i][1]),2));
            presbc_array[k][j][i] = 1./(2.*pi)*log(dist);
            presbc_array[k][j][i] = 10.;
          }
          if(j == ny-1){
            dist = sqrt(pow((ctx.well[0].coords[0]-coords_array[k][j][i][0]),2)+pow((ctx.well[0].coords[1]-coords_array[k][j][i][1]),2));
            presbc_array[k][j][i] = 1./(2.*pi)*log(dist);
            presbc_array[k][j][i] = 10.;
          }
      }
    }  
  }
  
  
  
  
  
  /*Beginning of mechanical part of code*/
  
  PetscReal       p;
  /*
	 Reset all BC for U and V
	 */
	for (i = 0; i < 6; i++) {
		ctx.bcV[0].face[i]=NONE;
		for (j = 0; j < 3; j++) {
			ctx.bcU[j].face[i] = NONE;
		}
	}
	for (i = 0; i < 12; i++) {
		ctx.bcV[0].edge[i]=NONE;
		for (j = 0; j < 3; j++) {
			ctx.bcU[j].edge[i] = NONE;
		}
	}
	for (i = 0; i < 8; i++) {
		ctx.bcV[0].vertex[i]=NONE;
		for (j = 0; j < 3; j++) {
			ctx.bcU[j].vertex[i] = NONE;
		}
	}
  /*	face X0	*/
  ctx.bcU[0].face[X0]= ZERO;
  ctx.bcU[1].face[X0]= ZERO;
  ctx.bcU[2].face[X0]= ZERO;
  /*	face X1	*/
  ctx.bcU[0].face[X1]= ZERO;
  ctx.bcU[1].face[X1]= ZERO;
  ctx.bcU[2].face[X1]= ZERO;
  /*	face Y0	*/
  ctx.bcU[0].face[Y0]= ZERO;
  ctx.bcU[1].face[Y0]= ZERO;
  ctx.bcU[2].face[Y0]= ZERO;
  /*	face Y1	*/
  ctx.bcU[0].face[Y1]= ZERO;
  ctx.bcU[1].face[Y1]= ZERO;
  ctx.bcU[2].face[Y1]= ZERO;
  /*	face Z0	*/
  ctx.bcU[2].face[Z0]= ZERO;
  /*	face Z1	*/
  ctx.bcU[2].face[Z1]= ZERO;
	ctx.hasCrackPressure = PETSC_TRUE;  
  ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
	p = 1.e-5;
	ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
	ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
 /*End of mechanical part of code*/

  
  
  
  
  
  
		
	ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-theta",&ctx.flowprop.theta,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-timestepsize",&ctx.flowprop.timestepsize,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
	ctx.maxtimestep = 1;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-maxtimestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr);
	for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
    ierr = PetscPrintf(PETSC_COMM_WORLD,"   Step: %d\n",ctx.timestep);CHKERRQ(ierr);
    ierr = VF_PermeabilityUpDate(&ctx,&fields);CHKERRQ(ierr);
    ierr = VecCopy(fields.V,V_hold);CHKERRQ(ierr);
    ierr = VecSet(fields.V,1.0);CHKERRQ(ierr);
    ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
    ierr = VecCopy(V_hold,fields.V);CHKERRQ(ierr);
    ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);
    ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
    ierr = FieldsH5Write(&ctx,&fields);
	}
  
  
  
	ierr = DMDAVecRestoreArray(ctx.daScal,ctx.Source,&src_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

