/*
 test10.c:
 Validate crack opening computation
 
 (c) 2010-2012 chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFCracks.h"
#include "VFV.h"
#include "VFU.h"
#include "VFFlow.h"

VFCtx    ctx;
VFFields fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	VFCtx          ctx;
	VFFields       fields;
	PetscErrorCode ierr;
	
	PetscInt     i,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm;
	PetscReal    ****coords_array;
	PetscInt            orientation=1;
	PetscReal    ElasticEnergy = 0;
	PetscReal    InsituWork    = 0;
	PetscReal    SurfaceEnergy = 0;
	char         filename[FILENAME_MAX];
	PetscReal    p = 1e-3;
	PetscReal    x[3],xc[3];
	PetscReal    dist,d;
	VFPennyCrack *crack;
	PetscInt     nc = 0;
	char         prefix[PETSC_MAX_PATH_LEN+1];
	Vec          V;
	PetscReal ****bcu_array;
	PetscReal bcx = .5;
	PetscReal bcy = .05;
	PetscReal bcz = .05;

	PetscInt			altminit=1;
	Vec					Vold;
	PetscReal			errV=1e+10;
	
	
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	
	ierr = PetscOptionsGetInt(PETSC_NULL,"-nc",&nc,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	
	ierr = PetscMalloc(nc*sizeof(VFPennyCrack),&crack);CHKERRQ(ierr);
	for (i = 0; i < nc; i++) {
		ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN,"c%d_",i);CHKERRQ(ierr);
		ierr = VFPennyCrackCreate(&crack[i]);CHKERRQ(ierr);
		ierr = VFPennyCrackGet(prefix,&crack[i]);CHKERRQ(ierr);
		ierr = VFPennyCrackView(&crack[i],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	}
	
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(ctx.daScal,&V);CHKERRQ(ierr);
	ierr = VecSet(V,1.0);CHKERRQ(ierr);
	
	ierr = VecSet(fields.V,1.0);CHKERRQ(ierr);
	
	/*
	 Reset all BC for U and V
	 */
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

	for (c = 0; c < nc; c++) {
		ierr = VFPennyCrackBuildVAT2(V,&crack[c],&ctx);CHKERRQ(ierr);
		ierr = VecPointwiseMin(fields.V,V,fields.V);CHKERRQ(ierr);
	}
	ierr = VecDestroy(&V);CHKERRQ(ierr);
	ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
	
	ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
	ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
	
	switch (orientation) {
		case 1:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a penny crack of radius");CHKERRQ(ierr);		  
			/*	face X0						face X1	*/
			ctx.bcU[0].face[X0]= FIXED;ctx.bcU[0].face[X1]= FIXED;
			/*	face Y0						face Y1	*/
			ctx.bcU[1].face[Y0]= FIXED;ctx.bcU[1].face[Y1]= FIXED;
			/*	face Z0						face Z1	*/
			ctx.bcU[2].face[Z0]= FIXED;ctx.bcU[2].face[Z1]= FIXED;
			break;
		default:
			SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Orientation should be one of {1,2,3,4,5,6}, got %i\n",orientation);
			break;
	} 
	
	ctx.maxtimestep = 2;
	ctx.hasCrackPressure = PETSC_FALSE;
	ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
	for(ctx.timestep = 1; ctx.timestep < ctx.maxtimestep; ctx.timestep++){	
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Time step %i\n",ctx.timestep);CHKERRQ(ierr);
		ierr = VecSet(fields.BCU,0.0);CHKERRQ(ierr);
		ierr = DMDAVecGetArrayDOF(ctx.daVect,fields.BCU,&bcu_array);CHKERRQ(ierr);
		switch (orientation) {
			case 1:
				ierr                = PetscPrintf(PETSC_COMM_WORLD,"Applying normal displacement boundary condition on all faces\n");CHKERRQ(ierr);
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						for (i = xs; i < xs+xm; i++) {
							if (i == 0) {
								bcu_array[k][j][i][0] = -ctx.timestep*bcx;
							}
							if (i == nx-1) {
								bcu_array[k][j][i][0] = ctx.timestep*bcx;
							}
							if (j == 0) {
								bcu_array[k][j][i][1] = -ctx.timestep*bcy;
							}
							if (j == ny-1) {
								bcu_array[k][j][i][1] = ctx.timestep*bcy;
							}
							if (k == 0) {
								bcu_array[k][j][i][2] = -ctx.timestep*bcz;
							}
							if (k == nz-1) {
								bcu_array[k][j][i][2] = ctx.timestep*bcz;
							}
						}
					}
				}
				break;
			default:
				SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Orientation should be one of {1,2,3,4,5,6}, got %i\n",orientation);
				break;
		}
		ierr = DMDAVecRestoreArrayDOF(ctx.daVect,fields.BCU,&bcu_array);CHKERRQ(ierr);
		ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
		do {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"alt min step %i\n",altminit);CHKERRQ(ierr);
			ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
			ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
			ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);
			ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
			
			ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
			ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,"   Max. change on V: %e\n",errV);CHKERRQ(ierr);
			altminit++;
		} while (errV > ctx.altmintol && altminit <= ctx.altminmaxit);
		ctx.ElasticEnergy = 0;
		ctx.InsituWork    = 0;
		ctx.PressureWork  = 0.;
		ierr              = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,&fields,&ctx);CHKERRQ(ierr);
		ierr              = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);
		ctx.TotalEnergy = ctx.ElasticEnergy-ctx.InsituWork-ctx.PressureWork;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface energy:            %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic Energy:            %e\n",ctx.ElasticEnergy);CHKERRQ(ierr);
		if (ctx.hasCrackPressure) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of pressure forces:   %e\n",ctx.PressureWork);CHKERRQ(ierr);
		}
		if (ctx.hasInsitu) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of surface forces:    %e\n",ctx.InsituWork);CHKERRQ(ierr);
		}
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Total Mechanical energy:  %e\n",ctx.ElasticEnergy-InsituWork-ctx.PressureWork);CHKERRQ(ierr);
		ierr = VolumetricCrackOpening(&ctx.CrackVolume, &ctx, &fields);CHKERRQ(ierr);   
		switch (ctx.fileformat) {
			case FILEFORMAT_HDF5:
				ierr = FieldsH5Write(&ctx,&fields);
				break;
			case FILEFORMAT_BIN:
				ierr = FieldsBinaryWrite(&ctx,&fields);
				break;
		}
		altminit = 1;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"###################################################################\n\n\n");CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#        VF crack volume change = %e\t      \n", ctx.CrackVolume);CHKERRQ(ierr);
	}
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

