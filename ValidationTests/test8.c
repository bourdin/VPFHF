/*
 test8.c: Coupled fracture propagation and fracture fluid flow. 
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
 ./test8 -n 201,2,201 -l 4,0.01,4 -npc 1 -pc0_r 0.2 -pc0_center 2.,0.005,2 -pc0_thickness 0.05 -epsilon 0.04 -pc0_theta 0 -pc0_phi 150 -orientation 1 -nw 1
 ./test8 -n 101,2,101 -l 4,0.01,4 -npc 1 -pc0_r 0.2 -pc0_center 2.,0.005,2 -pc0_thickness 0.08 -epsilon 0.08 -pc0_theta 0 -pc0_phi 150 -orientation 1 -nw 1
 ./test8 -n 51,2,51 -l 4,0.01,4 -npc 1 -pc0_r 0.2 -pc0_center 2.,0.005,2 -pc0_thickness 0.16 -epsilon 0.16 -pc0_theta 0 -pc0_phi 150 -orientation 1 -nw 1
 
 
 ./test8 -n 21,2,21 -l 4,0.01,4 -npc 1 -pc0_r 0.2 -pc0_center 2.,0.005,2 -pc0_thickness 0.4 -epsilon 0.4 -pc0_theta 0 -pc0_phi 90 -orientation 1 -nw 1
 
./test8 -n 51,2,51 -l 4,0.01,4 -npc 1 -pc0_r 0.2 -pc0_center 2.,0.005,2 -pc0_thickness 0.16 -epsilon 0.16 -pc0_theta 0 -pc0_phi 150 -orientation 1 -nw 1 -Fracsnes_snes_max_linear_solve_fail  -fracsnes_snes_type qn -fracsnes_snes_max_linear_solve_fail
 
 
 
 ./test8 -n 101,2,101 -l 4,0.01,4 -npc 1 -pc0_r 0.2 -pc0_center 2.,0.005,2 -pc0_thickness 0.1 -epsilon 0.1 -pc0_theta 0 -pc0_phi 90 -orientation 1 -nw 1 -Fracsnes_snes_max_linear_solve_fail  -fracsnes_snes_type tr -fracsnes_snes_max_linear_solve_fail -fracsnes_pc_type lu -fracsnes_ksp_type preonly 
 
 
 ./test8 -n 51,2,51 -l 4,0.01,4 -npc 1 -pc0_r 0.2 -pc0_center 2.,0.005,2 -pc0_thickness 0.25 -epsilon 0.16 -pc0_theta 0 -pc0_phi 90 -orientation 1 -nw 1 -Fracsnes_snes_max_linear_solve_fail  -Fracsnes_snes_type tr -Fracsnes_snes_view  -Fracsnes_snes_monitor -Fracsnes_snes_max_it 100 
 
 ./test8 -n 101,2,101 -l 4,0.01,4 -npc 1 -pc0_r 0.2 -pc0_center 2.,0.005,2 -pc0_thickness 0.1 -epsilon 0.1 -pc0_theta 0 -pc0_phi 90 -orientation 1 -nw 1 -Fracsnes_snes_max_linear_solve_fail  -fracsnes_snes_type tr 
 
 
 ./test8 -n 101,2,101 -l 4,0.01,4 -npc 1 -pc0_r 0.2 -pc0_center 2.,0.005,2 -pc0_thickness 0.1 -epsilon 0.1 -pc0_theta 0 -pc0_phi 90 -orientation 1 -nw 1 -Fracsnes_snes_max_linear_solve_fail  -Fracsnes_snes_type tr -Fracsnes_snes_view  -Fracsnes_snes_monitor -Fracsnes_snes_max_it 100 
 
 ./test8 -n 201,2,201 -l 4,0.01,4 -npc 1 -pc0_r 0.2 -pc0_center 2.,0.005,2 -pc0_thickness 0.08 -epsilon 0.05 -pc0_theta 0 -pc0_phi 90 -orientation 1 -nw 1 -Fracsnes_snes_type qn
 
 ./test8 -n 101,2,101 -l 4,0.04,4 -npc 1 -pc0_r 0.2 -pc0_center 2.,0.005,2 -pc0_thickness 0.1 -epsilon 0.1 -pc0_theta 0 -pc0_phi 90 -orientation 1 -nw 1 -Fracsnes_snes_max_linear_solve_fail  -fracsnes_snes_type tr
 
 ./test8 -n 151,2,151 -l 4,0.01,4 -npc 1 -pc0_r 0.2 -pc0_center 2.,0.005,2 -pc0_thickness 0.1 -epsilon 0.1 -pc0_theta 0 -pc0_phi 90 -orientation 1 -nw 1 -Fracsnes_snes_max_linear_solve_fail  -fracsnes_snes_type tr
 
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFV.h"
#include "VFU.h"
#include "VFFlow.h"
#include "VFCracks.h"


VFCtx               ctx;
VFFields            fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	VFCtx               ctx;
	VFFields            fields;
	PetscErrorCode      ierr;
	PetscInt            orientation=2;
	PetscInt            nopts=3;
	PetscInt            i,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm,xs1,xm1,ys1,ym1,zs1,zm1;
	PetscReal           ****coords_array;
	PetscReal           BBmin[3],BBmax[3];
	PetscReal           x,y,z;  
	PetscReal           ElasticEnergy = 0;
	PetscReal           InsituWork = 0;
	PetscReal           SurfaceEnergy = 0;
	char                filename[FILENAME_MAX];
	PetscReal           lx,ly,lz;
	PetscReal           p = 1e-3;
	char                prefix[PETSC_MAX_PATH_LEN+1];
	PetscReal           errV=1e+10,errP;
	Vec                 Vold;
	PetscReal           p_epsilon = 1.e-4;
	PetscInt            altminit=1;

	PetscReal           ***presbc_array;
	PetscReal           ****velbc_array;
	PetscReal           ****fracvelbc_array;


	
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	
//	ctx.fractureflowsolver = FRACTUREFLOWSOLVER_NONE;
	ctx.fractureflowsolver = FRACTUREFLOWSOLVER_SNESMIXEDFEM;
//	ctx.flowsolver = FLOWSOLVER_SNESMIXEDFEM;
	ctx.flowsolver = FLOWSOLVER_NONE;

	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-orientation",&orientation,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daScal,BBmin,BBmax);CHKERRQ(ierr);	
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VecSet(fields.VIrrev,1.0);CHKERRQ(ierr);
	lz = BBmax[2]-BBmin[2];
	ly = BBmax[1]-BBmin[1];
	lx = BBmax[0]-BBmin[0];
	
// flow part
	ctx.hasFluidSources = PETSC_FALSE;
	ctx.hasFlowWells = PETSC_TRUE;
	ierr = DMDAVecGetArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.FracVelBCArray,&fracvelbc_array);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daVFperm,&xs1,&ys1,&zs1,&xm1,&ym1,&zm1);CHKERRQ(ierr);

//	ctx.numWells = 1;
	ctx.well[0].Qw = 0.01;
	ctx.well[0].coords[0] = lx/2.;
	ctx.well[0].coords[1] = ly/2.;
	ctx.well[0].coords[2] = lz/2.;
	ctx.well[0].condition = RATE;
	ctx.well[0].type = INJECTOR;
	ierr = PetscOptionsGetReal(PETSC_NULL,"-Qw",&ctx.well[0].Qw,PETSC_NULL);CHKERRQ(ierr);

	
  /*	Flow model settings	*/
  
  
   for (i = 0; i < 6; i++) {
     for (c = 0; c < 3; c++) {
       ctx.bcFracQ[c].face[i] = NONE;
     }
   }
   for (i = 0; i < 12; i++) {
     for (c = 0; c < 3; c++) {
       ctx.bcFracQ[c].edge[i] = NONE;
     }
   }
   for (i = 0; i < 8; i++) {
     for (c = 0; c < 3; c++) {
       ctx.bcFracQ[c].vertex[i] = NONE;
     }
   }
   ctx.bcFracQ[1].face[Y0] = VALUE;
   ctx.bcFracQ[1].face[Y1] = VALUE;
   for (k = zs; k < zs+zm; k++) {
     for (j = ys; j < ys+ym; j++) {
       for (i = xs; i < xs+xm; i++) {
         fracvelbc_array[k][j][i][1] = 0.;
       }
     }
   }
   

  
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
	ctx.bcP[0].face[X0] = VALUE;
	ctx.bcP[0].face[X1] = VALUE;
	ctx.bcP[0].face[Y0] = VALUE;
	ctx.bcP[0].face[Y1] = VALUE;
	ctx.bcQ[2].face[Z0] = VALUE;
	ctx.bcQ[2].face[Z1] = VALUE;
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				velbc_array[k][j][i][0] = 0.;
				velbc_array[k][j][i][1] = 0.;
				velbc_array[k][j][i][2] = 0.;
				presbc_array[k][j][i] = 0.;
			}
		}
	}
	
	ierr = DMDAVecRestoreArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.FracVelBCArray,&fracvelbc_array);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-theta",&ctx.flowprop.theta,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-timestepsize",&ctx.flowprop.timestepsize,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
	
//	ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);

	/*	Mechanical model settings	*/
	
	
	
	
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
	ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
	switch (orientation) {
		case 1:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a line crack\n");CHKERRQ(ierr);		  
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
			/*	face Y1	*/
			ctx.bcU[1].face[Y1]= ZERO; 
			/* BCV*/
			ctx.bcV[0].face[X0] = ONE;
			ctx.bcV[0].face[X1] = ONE;
			ctx.bcV[0].face[Z0] = ONE;
			ctx.bcV[0].face[Z1] = ONE;
			break;
		case 2:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a line crack\n");CHKERRQ(ierr);		  
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcU[1].face[X0]= ZERO;
			/*	face X1	*/
			ctx.bcU[0].face[X1]= ZERO;
			ctx.bcU[1].face[X1]= ZERO;
			/*	face Z0	*/
			ctx.bcU[1].face[Z0]= ZERO;
			ctx.bcU[2].face[Z0]= ZERO;
			/*	face Z1	*/
			ctx.bcU[1].face[Z1]= ZERO;
			ctx.bcU[2].face[Z1]= ZERO;		  
			/*	face Y0	*/
			ctx.bcU[1].face[Y0]= ZERO;
			/*	face Z1	*/
			ctx.bcU[1].face[Y1]= ZERO;  
			/* BCV*/
			ctx.bcV[0].face[X0] = ONE;
			ctx.bcV[0].face[X1] = ONE;
			ctx.bcV[0].face[Z0] = ONE;
			ctx.bcV[0].face[Z1] = ONE;
			break;
		case 3:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a line crack\n");CHKERRQ(ierr);		  
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcU[1].face[X0]= ZERO;
			ctx.bcU[2].face[X0]= ZERO;
			/*	face X1	*/
			ctx.bcU[0].face[X1]= ZERO;
			ctx.bcU[1].face[X1]= ZERO;
			ctx.bcU[2].face[X1]= ZERO;
			/*	face Z0	*/
			ctx.bcU[1].face[Z0]= ZERO;
			ctx.bcU[2].face[Z0]= ZERO;
			/*	face Z1	*/
			ctx.bcU[1].face[Z1]= ZERO;		  
			ctx.bcU[2].face[Z1]= ZERO;		  
			/*	face Y0	*/
			ctx.bcU[1].face[Y0]= ZERO;
			/*	face Y1	*/
			ctx.bcU[1].face[Y1]= ZERO; 
			/* BCV*/
			ctx.bcV[0].face[X0] = ONE;
			ctx.bcV[0].face[X1] = ONE;
			ctx.bcV[0].face[Z0] = ONE;
			ctx.bcV[0].face[Z1] = ONE;
			break;
		case 4:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a line crack\n");CHKERRQ(ierr);		  
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcU[1].face[X0]= ZERO;
			ctx.bcU[0].face[X1]= ZERO;
			ctx.bcU[1].face[X1]= ZERO;
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
			break;
		case 5:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a line crack\n");CHKERRQ(ierr);		  
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcU[1].face[X0]= ZERO;
			ctx.bcU[0].face[X1]= ZERO;
			ctx.bcU[1].face[X1]= ZERO;
			ctx.bcU[1].face[Y0]= ZERO;
			ctx.bcU[1].face[Y1]= ZERO;
			ctx.bcV[0].face[X0] = ONE;
			ctx.bcV[0].face[X1] = ONE;
			ctx.bcV[0].face[Z0] = ONE;
			ctx.bcV[0].face[Z1] = ONE;
			break;
		case 6:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a line crack\n");CHKERRQ(ierr);		  
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcU[1].face[X0]= ZERO;
			ctx.bcU[2].face[X0]= ZERO;
			ctx.bcU[0].face[X1]= ZERO;
			ctx.bcU[1].face[X1]= ZERO;
			ctx.bcU[2].face[X1]= ZERO;
			ctx.bcU[1].face[Y0]= ZERO;
			ctx.bcU[1].face[Y1]= ZERO;  
			ctx.bcV[0].face[X0] = ONE;
			ctx.bcV[0].face[X1] = ONE;
			ctx.bcV[0].face[Z0] = ONE;
			ctx.bcV[0].face[Z1] = ONE;
			break;
		default:
			SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Orientation should be one of {1,2,3,4,5,6}, got %i\n",orientation);
			break;
	} 
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
	ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
	ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-theta",&ctx.flowprop.theta,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-timestepsize",&ctx.flowprop.timestepsize,PETSC_NULL);CHKERRQ(ierr);
	ctx.matprop[0].beta = 0.;
	ctx.hasCrackPressure = PETSC_TRUE;
	ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
	altminit = 0.;
//	do {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"  alt min step %i with errorV %g\n",altminit,errV);CHKERRQ(ierr);

		ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
		ierr = VF_StepU(&fields,&ctx);
		ierr = VF_StepV(&fields,&ctx);
		
		ierr = PermeabilityUpDate(&ctx,&fields);CHKERRQ(ierr);
		ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
		ierr = VolumetricLeakOffRate(&ctx.LeakOffRate,&ctx,&fields);CHKERRQ(ierr);

		
		
		ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);		
		
		ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
		ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
		
		ierr = FieldsH5Write(&ctx,&fields);

//		altminit++;
//	} while (errV >= ctx.altmintol && altminit <= ctx.altminmaxit);
	
	ctx.ElasticEnergy=0;
	ctx.InsituWork=0;
	ctx.PressureWork = 0.;
	ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,&fields,&ctx);CHKERRQ(ierr);
	ctx.TotalEnergy = ctx.ElasticEnergy - ctx.InsituWork - ctx.PressureWork;	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic Energy:            %e\n",ctx.ElasticEnergy);CHKERRQ(ierr);
	if (ctx.hasCrackPressure) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of pressure forces:   %e\n",ctx.PressureWork);CHKERRQ(ierr);
	}
	if (ctx.hasInsitu) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of surface forces:    %e\n",ctx.InsituWork);CHKERRQ(ierr);
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Total energy:              %e\n",ctx.ElasticEnergy-InsituWork-ctx.PressureWork);CHKERRQ(ierr);
//	ierr = VolumetricCrackOpening(&ctx.CrackVolume, &ctx, &fields);CHKERRQ(ierr);   
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n Final Crack volume\t = %g, Pressure\t= %g\n\n", ctx.CrackVolume, p);CHKERRQ(ierr);

	
//	ierr = PermeabilityUpDate(&ctx,&fields);CHKERRQ(ierr);

	ierr = FieldsH5Write(&ctx,&fields);
	ierr = VecDestroy(&Vold);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

