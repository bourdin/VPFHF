/*
 test39.c: 2D SNES. Coupled flow and fracture problem using VFRectangularCrackCreate.
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
 ./test39 -n 101,101,2 -epsilon 0.01 -Gc 0.6667  -l 1,1,0.01 -m_inv 0. -Qw 0.1 -maxtimestep 5 -nrc 1 -rc0_corners 0.45,0.5,0,0.55,0.5,0.0,0.44,0.5,0.01 -rc0_thickness 0.005 -nw 1
 
 ./test39 -n 101,101,2 -epsilon 0.01 -Gc 0.6667  -l 1,1,0.01 -m_inv 5. -maxtimestep 5 -nrc 2 -rc0_corners 0.25,0.2,0,0.7,0.5,0.0,0.25,0.2,0.01 -rc0_thickness 0.008 -rc1_corners 0.5,0.6,0,0.7,0.3,0.,0.5,0.6,0.01 -rc1_thickness 0.008 -w0_coords 0.5,0.6,0.005 -w0_Qw 0.1 -maxtimestep 20 -nw 1
 
 ./test39 -n 101,101,2 -epsilon 0.01 -Gc 0.6667  -l 1,1,0.01 -m_inv 5. -maxtimestep 5 -nrc 2 -rc0_corners 0.25,0.2,0,0.7,0.5,0.0,0.25,0.2,0.01 -rc0_thickness 0.008 -rc1_corners 0.5,0.6,0,0.7,0.3,0.,0.5,0.6,0.01 -rc1_thickness 0.008 -w0_coords 0.5,0.6,0.005 -w0_Qw 0.1 -maxtimestep 20 -nw 1
 
 ./test39 -n 101,101,2 -epsilon 0.02 -Gc 0.6667  -l 1,1,0.01 -m_inv 5. -maxtimestep 10 -npc 1 -pc0_r 0.2 -pc0_center 0.5,0.5,0.005 -pc0_thickness 0.05 -epsilon 0.04 -pc0_theta 0 -pc0_phi 90 -nw 1 -npc 1
 
mpirun -np 4 $VFDIR/ValidationTests/test39 -n 101,101,2 -l 1000,1000,10 -m_inv 0 -maxtimestep 3 -nw 1 -timestepsize 8640 -w0_coords 500,500,5 -w0_Qw 1. -w0_constraint Rate -w0_rw 1 -w0_type injector\
 -w0_top 500,500,0 -w0_bottom 500,500,10 -flowsolver FLOWSOLVER_SNES -npc 1 -pc0_r 200 -pc0_center 500,500,5 -pc0_thickness 20 -epsilon 15 -pc0_theta 90 -pc0_phi 90

For TS
mpirun -np 1 $VFDIR/ValidationTests/test39 -n 101,101,2 -l 1000,1000,10 -m_inv 3.e-4 -nw 1 -w0_coords 500,500,5 -w0_Qw 1. -w0_constraint Rate -w0_rw 1 -w0_type injector -w0_top 500,500,0 -w0_bottom 500,500,10 \
-flowsolver FLOWSOLVER_TS -npc 1 -pc0_r 200 -pc0_center 500,500,5 -pc0_thickness 20 -epsilon 15 -pc0_theta 90 -pc0_phi 90 -timestepsize 864

 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFMech.h"
#include "VFFlow.h"
#include "VFWell.h"
#include "VFCracks.h"

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
	PetscReal		***presbc_array;
	PetscReal		****velbc_array;
	PetscReal		****coords_array;
	PetscReal		gx,gy,gz;
	PetscReal		gamma, beta, rho, mu;
	PetscReal		****perm_array;
	PetscReal       ****velnpre_array;
	char			prefix[PETSC_MAX_PATH_LEN+1];
	PetscReal       p;
	PetscReal       p_old;
	PetscReal       p_epsilon = 1.e-5;
	PetscInt        altminit=1;
	Vec	            Uold,Vold,Pold;
	PetscReal       errU,errV=1e+10;
	PetscReal       lx,ly,lz;

		
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	
    ctx.flowsolver = FLOWSOLVER_SNES;
    ctx.fractureflowsolver = FRACTUREFLOWSOLVER_NONE;	
	
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	ierr = VecSet(ctx.PresBCArray,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.VelBCArray,0.);CHKERRQ(ierr);
	ierr = VecSet(fields.VelnPress,0.);CHKERRQ(ierr);
	ierr = VecSet(fields.vfperm,0.);CHKERRQ(ierr);	
	ierr = VecSet(fields.pressure,0.);CHKERRQ(ierr);	
    ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
    ierr = VecDuplicate(fields.U,&Uold);CHKERRQ(ierr);	
    ierr = VecDuplicate(fields.pressure,&Pold);CHKERRQ(ierr);
	
	ctx.hasFluidSources = PETSC_FALSE;
	ctx.hasFlowWells = PETSC_TRUE;
	ierr = DMDAVecGetArrayDOF(ctx.daFlow,fields.VelnPress,&velnpre_array);CHKERRQ(ierr);	
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);	
	ierr = DMDAVecGetArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr); 
	ierr = DMDAGetCorners(ctx.daVFperm,&xs1,&ys1,&zs1,&xm1,&ym1,&zm1);CHKERRQ(ierr);
	lz = BBmax[2]-BBmin[2];
	ly = BBmax[1]-BBmin[1];
	lx = BBmax[0]-BBmin[0];	
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
	ctx.bcV[0].face[X0] = ONE;
	/*	face X1	*/
	ctx.bcU[0].face[X1]= ZERO;
	ctx.bcU[1].face[X1]= ZERO;
	ctx.bcU[2].face[X1]= ZERO;
	ctx.bcV[0].face[X1] = ONE;
	/*	face Y0	*/
	ctx.bcU[0].face[Y0]= ZERO;
	ctx.bcU[1].face[Y0]= ZERO;
	ctx.bcU[2].face[Y0]= ZERO;
	ctx.bcV[0].face[Y0] = ONE;	
	/*	face Y1	*/
	ctx.bcU[0].face[Y1]= ZERO;		  
	ctx.bcU[1].face[Y1]= ZERO;		  
	ctx.bcU[2].face[Y1]= ZERO;		  
	ctx.bcV[0].face[Y1] = ONE;
	/*	face Z0	*/
	ctx.bcU[2].face[Z0]= ZERO;
	/*	face Z1	*/
	ctx.bcU[2].face[Z1]= ZERO;
	
	
	ierr = PetscOptionsGetInt(PETSC_NULL,"-nc",&ctx.numWells,PETSC_NULL);CHKERRQ(ierr);	
	ierr = PetscMalloc(ctx.numWells*sizeof(VFWell),&ctx.well);CHKERRQ(ierr);
	if(ctx.hasFlowWells == PETSC_TRUE){
		for (i = 0; i < ctx.numWells; i++) {
			ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN,"w%d_",i);CHKERRQ(ierr);
			ierr = VFWellCreate(&ctx.well[i]);CHKERRQ(ierr);
			ierr = VFWellGet(prefix,&ctx.well[i]);CHKERRQ(ierr);
			ierr = VFWellView(&ctx.well[i],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
		}	
	} 
							
	/*
	 Reset all Flow BC for velocity and P
	 */
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
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);	
	ctx.maxtimestep = 1;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-maxtimestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-theta",&ctx.flowprop.theta,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-timestepsize",&ctx.flowprop.timestepsize,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-m_inv",&ctx.flowprop.M_inv,PETSC_NULL);CHKERRQ(ierr);
	
	ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
	ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
	ctx.hasCrackPressure = PETSC_TRUE;
	ctx.matprop[0].beta = 0.;	
	ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
	ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);

	ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
    
	ctx.current_time = 0.;
	ctx.timevalue = 0.;
//	PetscInt local_timestep =0;
    ierr = VF_PermeabilityUpDate(&ctx,&fields);CHKERRQ(ierr);	
	
	for (ctx.timestep = 0; ctx.timestep < 20; ctx.timestep++){	
//		ctx.timestep = 0;
        ctx.current_time = ctx.timestep*ctx.flowprop.timestepsize;
        ctx.timevalue = (ctx.timestep+1)*ctx.flowprop.timestepsize;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nProcessing step %i. Current time %g Timevalue %g timestepsize %g\n",ctx.timestep,ctx.current_time,ctx.timevalue,ctx.flowprop.timestepsize);CHKERRQ(ierr);
 
//        ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);	
		ierr = VecCopy(fields.pressure,Pold);CHKERRQ(ierr);	
      do{		
		do {
		    ierr = VecCopy(Pold,fields.pressure);CHKERRQ(ierr);	
		    ierr = VF_PermeabilityUpDate(&ctx,&fields);CHKERRQ(ierr);		
            ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);		
            ierr = VecCopy(fields.U,Uold);CHKERRQ(ierr);	
		    ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
            ierr = VecAXPY(Uold,-1.,fields.U);CHKERRQ(ierr);
            ierr = VecNorm(Uold,NORM_INFINITY,&errU);CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"   Max change on U: %e \n",errU);CHKERRQ(ierr);		
	        altminit++;
		} while ((errU >= ctx.altmintol) && altminit <= ctx.altminmaxit);
		   		
		
		ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
//		ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);
        ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
        ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"   Max change on V: %e\n",errV);CHKERRQ(ierr);
		
		} while ((errV >= ctx.altmintol) && altminit <= ctx.altminmaxit);
		
/*		ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr); 
		ierr = VF_PermeabilityUpDate(&ctx,&fields);CHKERRQ(ierr);
		ierr = VolumetricLeakOffRate(&ctx.LeakOffRate,&ctx,&fields);CHKERRQ(ierr);   
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\nInjectedVol = %e,\tCrackVolume =%e\tLeak-off = %e\n",0.01*ctx.flowprop.timestepsize,ctx.CrackVolume,ctx.LeakOffRate);CHKERRQ(ierr);*/
		ierr = FieldsH5Write(&ctx,&fields);	
	}
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

