/*
 test36.c: 1D. Flow problem with source term = 1 and Homogeneous pressure boundary conditions on all sides. Analytical solution is p = x(x-1)/2
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
./test36 -n 101,2,101 -l 1,0.01,1 -m_inv 0 -E 17 -0 -nu 0.2 -npc 1 -pc0_r 0.1 -pc0_center 0.3,0.005,0.5 -epsilon 0.03eta 0 -pc0_phi 40 -pc0_thickness 0.01 -atnum 1 -flowsolver FLOWSOLVER_snesmixeDFEM -perm 1e-14 -miu 1e-13 -num 2 -timestepsize 0.1  -Gc 5e6
 
 
 ./test36 -n 51,2,51 -l 1,0.01,1 -m_inv 0 -E 17 -0 -nu 0.2 -npc 1 -pc0_r 0.1 -pc0_center 0.5,0.005,0.5 -epsilon 0.02 -pc0_theta 0 -pc0_phi 30 -pc0_thickness 0.04 -atnum 2 -flowsolver FLOWSOLVER_snesmixeDFEM -perm 1e-14 -miu 1e-13 -num 2 -timestepsize 0.1  -Gc 5e6 -permmax 1e-12
 */

#include "petsc.h"
#include "CartFE.h"
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
	PetscViewer		viewer;
	PetscViewer     logviewer;
	char			filename[FILENAME_MAX];
	PetscInt		i,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm,ite;
	PetscReal		BBmin[3],BBmax[3];
	PetscReal		****velbc_array;
	PetscReal		***src_array;
	PetscReal		****coords_array;
	PetscReal		hx,hy,hz;
	PetscReal		lx,ly,lz;
  PetscReal		***presbc_array;
  PetscReal    perm = 1;
  PetscReal   ****perm_array;
  PetscInt    xs1,xm1,ys1,ym1,zs1,zm1;
  PetscReal  tolP = 7e-2,tolV = 7e-2;
  PetscInt    num = 10;
  PetscReal      errP=1e+10,errV=1e+10;
  Vec         Pold, Vold;
  PetscReal  pmax,vmax;
  PetscReal   InjVolrate, Q_inj;
  PetscReal  TotalLeakOff_o = 0;
  PetscReal  timevalue_o = 0;
  PetscReal  crackvolume_o = 0;
  PetscReal   vol,vol1,vol2,vol3,vol4,vol5;
  PetscReal   p = 1e-6,p_old;
  PetscInt    altminit = 0;

		
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ctx.fractureflowsolver = FRACTUREFLOWSOLVER_NONE;
	ctx.flowsolver = FLOWSOLVER_KSPMIXEDFEM;
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx.daVFperm,&xs1,&ys1,&zs1,&xm1,&ym1,&zm1);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	ierr = VecSet(ctx.VelBCArray,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.Source,1.);CHKERRQ(ierr);
	ierr = VecSet(fields.U,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.U_old,0.);CHKERRQ(ierr);
  ierr = VecDuplicate(fields.pressure,&Pold);CHKERRQ(ierr);
  ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);

	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);	
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);
  
  ierr = PetscOptionsGetReal(PETSC_NULL,"-perm",&perm,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-miu",&ctx.flowprop.mu,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-theta",&ctx.flowprop.theta,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-timestepsize",&ctx.flowprop.timestepsize,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-maxtimestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-num",&num,PETSC_NULL);CHKERRQ(ierr);
  
  for (k = zs1; k < zs1+zm1; k++) {
    for (j = ys1; j < ys1+ym1; j++) {
      for (i = xs1; i < xs1+xm1; i++) {
        perm_array[k][j][i][0] = perm;
        perm_array[k][j][i][1] = perm;
        perm_array[k][j][i][2] = perm;
        perm_array[k][j][i][3] = 0.;
        perm_array[k][j][i][4] = 0.;
        perm_array[k][j][i][5] = 0.;
      }
    }
  }
	lz = BBmax[2]-BBmin[2];
	ly = BBmax[1]-BBmin[1];
	lx = BBmax[0]-BBmin[0];
	hx = lx/(nx-1);
	hy = ly/(nx-1);
	hz = lz/(nz-1);	
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
  /*
   ctx.bcP[0].face[X0] = FIXED;
   ctx.bcP[0].face[X1] = FIXED;
   ctx.bcQ[1].face[Y0] = FIXED;
   ctx.bcQ[1].face[Y1] = FIXED;
   ctx.bcQ[2].face[Z0] = FIXED;
   ctx.bcQ[2].face[Z1] = FIXED;
  */
  
   ctx.bcQ[0].face[X0] = FIXED;
   ctx.bcP[0].face[X1] = FIXED;
   ctx.bcQ[1].face[Y0] = FIXED;
   ctx.bcQ[1].face[Y1] = FIXED;
   ctx.bcP[0].face[Z0] = FIXED;
   ctx.bcP[0].face[Z1] = FIXED;
   
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
        velbc_array[k][j][i][0] = 0.;
        velbc_array[k][j][i][1] = 0.;
        velbc_array[k][j][i][2] = 0.;
        if(i == 0){
          velbc_array[k][j][i][0] = 0.1;
          velbc_array[k][j][i][1] = 0.;
          velbc_array[k][j][i][2] = 0.;
        }
        if(i == nx-1){
          velbc_array[k][j][i][0] = 0.1;
          velbc_array[k][j][i][1] = 0.;
          velbc_array[k][j][i][2] = 0.;
        }
				presbc_array[k][j][i] = 0.;
			}
		}
	}	
/*      Mechanical model settings       */
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
  ctx.bcU[0].face[X0]= ZERO;
  ctx.bcU[1].face[X0]= ZERO;
  ctx.bcU[2].face[X0]= ZERO;
  
  ctx.bcU[0].face[X1]= ZERO;
  ctx.bcU[1].face[X1]= ZERO;
  ctx.bcU[2].face[X1]= ZERO;
  
  ctx.bcU[0].face[Z0]= ZERO;
  ctx.bcU[1].face[Z0]= ZERO;
  ctx.bcU[2].face[Z0]= ZERO;
  
  ctx.bcU[0].face[Z1]= ZERO;
  ctx.bcU[1].face[Z1]= ZERO;
  ctx.bcU[2].face[Z1]= ZERO;
  
  ctx.bcU[1].face[Y0]= ZERO;
  ctx.bcU[1].face[Y1]= ZERO;
  
  ctx.bcV[0].face[X0] = ONE;
  ctx.bcV[0].face[X1] = ONE;
  ctx.bcV[0].face[Z0] = ONE;
  ctx.bcV[0].face[Z1] = ONE;
  
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);
  ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,0);CHKERRQ(ierr);
  ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
  ierr = VecSet(ctx.U_old,0.);CHKERRQ(ierr);
  
  ctx.flowprop.alphabiot = 	ctx.matprop[0].beta = 0.75;									//biot's constant
  ctx.flowprop.K_dr = ctx.matprop[0].E/(3*(1-2*ctx.matprop[0].nu));   //For 3D

  ctx.hasCrackPressure = PETSC_TRUE;
  ctx.FlowDisplCoupling = PETSC_FALSE;
  ctx.ResFlowMechCoupling = FIXEDSTRESS;
  
  ctx.FractureFlowCoupling = PETSC_TRUE;
  ctx.hasFluidSources = PETSC_FALSE;
  ctx.hasFlowWells = PETSC_FALSE;

  ierr = VecSet(ctx.PreFlowFields,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.RHSVelPpre,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.pressure_old,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.RHSPpre,0.);CHKERRQ(ierr);
  

  ctx.timestep = 0;
  
  ierr = VecSet(fields.pressure,5e-2);CHKERRQ(ierr);
  ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
  ierr = FieldsH5Write(&ctx,&fields);


  ctx.timestep++;
//  ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
  ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
  ierr = FieldsH5Write(&ctx,&fields);

  
  
  
  ierr = VF_PermeabilityUpDate(&ctx,&fields);CHKERRQ(ierr);
  ctx.FractureFlowCoupling = PETSC_FALSE;
  ierr = VecSet(fields.pressure,0);CHKERRQ(ierr);

  ierr = VecSet(fields.V,1.0);CHKERRQ(ierr);
  ierr = VecSet(ctx.pressure_old,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.RHSP,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.RHSPpre,0.);CHKERRQ(ierr);
  ierr = VecSet(Pold,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.RHSVelP,0.);CHKERRQ(ierr);
  ierr = VecSet(ctx.RHSVelPpre,0.);CHKERRQ(ierr);

  
  ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
  ctx.timestep++;
  ierr = FieldsH5Write(&ctx,&fields);

  
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

