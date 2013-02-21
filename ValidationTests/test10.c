/*
 test10.c: Solves for the displacement and v-field in a pressurized penny crack in 2d (Sneddon 2D) using VFPennyCracks and VFRectangularcracks
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
 ./test10 -n 201,2,201 -l 4,0.01,4 -npc 1 -pc0_r 0.2 -pc0_center 2.,0.005,2 -pc0_thickness 0.05 -epsilon 0.04 -pc0_theta 0 -pc0_phi 150 -orientation 1
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
	PetscInt            i,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm;
	PetscReal			****coords_array;
	PetscReal           BBmin[3],BBmax[3];
	PetscReal           x,y,z;  
	PetscReal           ElasticEnergy = 0;
	PetscReal           InsituWork = 0;
	PetscReal           SurfaceEnergy = 0;
	char                filename[FILENAME_MAX];
	PetscReal			lx,ly,lz;
	PetscReal           p = 1e-3;
	char				prefix[PETSC_MAX_PATH_LEN+1];
	PetscReal			errV=1e+10,errP;
	Vec					Vold;
	PetscReal			p_epsilon = 1.e-4;
	PetscInt			altminit=1;


	
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-orientation",&orientation,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daScal,BBmin,BBmax);CHKERRQ(ierr);	
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VecSet(fields.VIrrev,1.0);CHKERRQ(ierr);
	lz = BBmax[2];
	ly = BBmax[1];
	lx = BBmax[0];	
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
	/*Initializing the v-field*/	
	/*
	ierr = PetscOptionsGetInt(PETSC_NULL,"-nrc",&ctx.numRectangularCracks,PETSC_NULL);CHKERRQ(ierr);	
	ierr = PetscMalloc(ctx.numRectangularCracks*sizeof(VFRectangularCrack),&ctx.rectangularcrack);CHKERRQ(ierr);
	for (i = 0; i < ctx.numRectangularCracks; i++) {
		ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN,"c%d_",i);CHKERRQ(ierr);
		ierr = VFRectangularCrackCreate(&ctx.rectangularcrack[i]);CHKERRQ(ierr);
		ierr = VFRectangularCrackGet(prefix,&ctx.rectangularcrack[i]);CHKERRQ(ierr);
		ierr = VFRectangularCrackView(&ctx.rectangularcrack[i],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	}
	for (c = 0; c < ctx.numRectangularCracks; c++) {
		ierr = VFRectangularCrackBuildVAT2(fields.VIrrev,&ctx.rectangularcrack[c],&ctx);CHKERRQ(ierr);
		ierr = VecPointwiseMin(fields.V,fields.VIrrev,fields.V);CHKERRQ(ierr);
		
	}
	
	
	ierr = PetscOptionsGetInt(PETSC_NULL,"-npc",&ctx.numPennyCracks,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscMalloc(ctx.numPennyCracks*sizeof(VFPennyCrack),&ctx.pennycrack);CHKERRQ(ierr);
	for (i = 0; i < ctx.numPennyCracks; i++) {
		ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN,"pc%d_",i);CHKERRQ(ierr);
		ierr = VFPennyCrackCreate(&ctx.pennycrack[i]);CHKERRQ(ierr);
		ierr = VFPennyCrackGet(prefix,&ctx.pennycrack[i]);CHKERRQ(ierr);
		ierr = VFPennyCrackView(&ctx.pennycrack[i],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	}	
	for (c = 0; c < ctx.numPennyCracks; c++) {
		ierr = VecSet(fields.VIrrev,1.0);CHKERRQ(ierr);
		ierr = VFPennyCrackBuildVAT2(fields.VIrrev,&ctx.pennycrack[c],&ctx);CHKERRQ(ierr);
		ierr = VecPointwiseMin(fields.V,fields.VIrrev,fields.V);CHKERRQ(ierr);
		
	}
	
	*/
	
	
	ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
	switch (orientation) {
		case 1:
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
			ctx.bcU[0].face[Z0]= ZERO;
			ctx.bcU[1].face[Z0]= ZERO;
			ctx.bcU[2].face[Z0]= ZERO;
			/*	face Z1	*/
			ctx.bcU[0].face[Z1]= ZERO;		  
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
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcU[1].face[X0]= ZERO;
			/*	face X1	*/
			ctx.bcU[0].face[X1]= ZERO;
			ctx.bcU[1].face[X1]= ZERO;
			/*	face Z0	*/
			ctx.bcU[0].face[Z0]= ZERO;
			ctx.bcU[1].face[Z0]= ZERO;
			ctx.bcU[2].face[Z0]= ZERO;
			/*	face Z1	*/
			ctx.bcU[0].face[Z1]= ZERO;		  
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
		case 5:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a line crack\n");CHKERRQ(ierr);		  
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcU[1].face[X0]= ZERO;
			/*	face X1	*/
			ctx.bcU[0].face[X1]= ZERO;
			ctx.bcU[1].face[X1]= ZERO;
			/*	face Z0	*/
			/*.....FREE.......*/
			/*	face Z1	*/
			/*.....FREE.......*/
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
		case 6:
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
			/*.....FREE.......*/
			/*	face Z1	*/
			/*.....FREE.......*/
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
	ctx.matprop[0].beta = 0.;
	ctx.hasCrackPressure = PETSC_TRUE;
	ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
	altminit = 0.;
	do {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"  alt min step %i with errorV %g\n",altminit,errV);CHKERRQ(ierr);

		ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
		ierr = VF_StepU(&fields,&ctx);
		ierr = VF_StepV(&fields,&ctx);
		ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
		ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
		altminit++;
	} while (errV >= ctx.altmintol && altminit <= ctx.altminmaxit);
	
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
	ierr = VolumetricCrackOpening(&ctx.CrackVolume, &ctx, &fields);CHKERRQ(ierr);   
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n Final Crack volume\t = %g, Pressure\t= %g\n\n", ctx.CrackVolume, p);CHKERRQ(ierr);	
	/*
	 Save fields and write statistics about current run
	 */    
	
    switch (ctx.fileformat) {
		case FILEFORMAT_HDF5:       
			ierr = FieldsH5Write(&ctx,&fields);
			break;
		case FILEFORMAT_BIN:
			ierr = FieldsBinaryWrite(&ctx,&fields);
			break; 
    } 
	ierr = VecDestroy(&Vold);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

