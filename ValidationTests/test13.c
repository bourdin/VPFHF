/*
 test13.c:
 Chevron internship test: Crack initiation and propagation under monotonically increasing tension, sliding and shear forces respectively, in the presence of weak (seed) points and pre-existing cracks.
 
 (c) 2010-2011 chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 */

#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFMech.h"
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
	PetscViewer			viewer;	
	PetscReal length      = .0;
	PetscInt  mode = 1;
	PetscInt  i,ii,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm;
	PetscReal ****coords_array;
	PetscReal ****bcu_array;
	PetscReal BBmin[3],BBmax[3];
	PetscReal ElasticEnergy = 0;
	PetscReal InsituWork    = 0;
	PetscReal SurfaceEnergy = 0;
	char      filename[FILENAME_MAX];
	PetscReal bc1 = 0.05;
	PetscReal bc2 = 0.0;
	PetscReal ***v_array;
	PetscReal lx,ly,lz;
	PetscInt			altminit=1;
	Vec					Vold;
	PetscReal			errV=1e+10;
	PetscInt	no_seed = 10;
	PetscInt	*seed_array_z;
	PetscInt	*seed_array_x;
	PetscInt	seedz_start;
	PetscInt	seedx_start;
	VFRectangularCrack *crack;
	PetscInt     nc = 0;
	char         prefix[PETSC_MAX_PATH_LEN+1];
	Vec          V;


	
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-no_seed",&no_seed,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-max_timestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-length",&length,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-mode",&mode,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-bcUx",&bc1,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-bcUz",&bc2,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-nc",&nc,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
	lz = BBmax[2];
	ly = BBmax[1];
	lx = BBmax[0];	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Bounding box: %e, %e, %e\n",lx,ly,lz);CHKERRQ(ierr);
	ctx.matprop[0].beta  = 0.;
	ctx.matprop[0].alpha = 0.;
	ctx.timestep  = 1;
	ctx.maxtimestep = 100;
	
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx.daScal,&V);CHKERRQ(ierr);
	ierr = VecSet(V,1.0);CHKERRQ(ierr);
	ierr = VecSet(fields.V,1.0);CHKERRQ(ierr);
	ierr = VecSet(fields.VIrrev,1.0);CHKERRQ(ierr);
	ierr = VecSet(fields.U,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.pressure,0.0);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx.daScal,fields.V,&v_array);CHKERRQ(ierr);
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
	/*Initializing the v-field*/	
	ierr = PetscMalloc(nc*sizeof(VFRectangularCrack),&crack);CHKERRQ(ierr);
	for (i = 0; i < nc; i++) {
		ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN,"c%d_",i);CHKERRQ(ierr);
		ierr = VFRectangularCrackCreate(&crack[i]);CHKERRQ(ierr);
		ierr = VFRectangularCrackGet(prefix,&crack[i]);CHKERRQ(ierr);
		ierr = VFRectangularCrackView(&crack[i],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	}
	for (c = 0; c < nc; c++) {
		ierr = VFRectangularCrackBuildVAT2(V,&crack[c],&ctx);CHKERRQ(ierr);
		ierr = VecPointwiseMin(fields.V,V,fields.V);CHKERRQ(ierr);
	}
	ierr = VecDestroy(&V);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx.daScal,fields.V,&v_array);CHKERRQ(ierr);
	/*Initializing seed points*/
	ierr = PetscMalloc2(no_seed,&seed_array_z,no_seed,&seed_array_x);CHKERRQ(ierr);
	seedz_start = seedx_start = 2;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-seedz_start",&seedz_start,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-seedx_start",&seedx_start,PETSC_NULL);CHKERRQ(ierr);
	
	srand((time(NULL)+1));
	for(i = 0; i < no_seed; i++){
		seed_array_z[i] = rand() % (nz-2*seedz_start)+seedz_start;
	}
	for(i = 0; i < no_seed; i++){
		seed_array_x[i] = rand() % (nx-2*seedx_start)+seedx_start;
	}
	
	PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
	PetscViewerSetType(viewer, PETSCVIEWERASCII);
	PetscViewerFileSetMode(viewer, FILE_MODE_APPEND);
	PetscViewerFileSetName(viewer, "seed_nodes.txt");
	PetscViewerASCIIPrintf(viewer, "%d \n", no_seed);
	PetscViewerASCIIPrintf(viewer, "seed \t z_nodes \t x_nodes \n");
	for(ii = 0; ii < no_seed; ii++){
		PetscViewerASCIIPrintf(viewer, "%d \t %d \t %d \n",ii, seed_array_z[ii],seed_array_x[ii] );
	}
	
	/*Initializing the v-field*/
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) { 
				if ( ((i == 3*nx/5) || (i == 3*nx/5-1)) && (coords_array[k][j][i][2] < 2*length) ) {
					v_array[k][j][i] = 0.;
				}
				if( ((i == 2*nx/5) || (i == 2*nx/5-1)) && coords_array[k][j][i][2] > (lz-2*length) ){
					v_array[k][j][i] = 0.;
				}
			}
		}
	}

	
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) { 
				for(ii = 0; ii < no_seed; ii++)
					if(k == seed_array_z[ii] && i == seed_array_x[ii]){
						v_array[k][j][i] = 0.;
					}
			}
		}
	}
	ierr = PetscFree2(seed_array_z,seed_array_x);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx.daScal,fields.VIrrev,&v_array);CHKERRQ(ierr);
	ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ctx.hasCrackPressure = PETSC_FALSE;
	for(ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++){	
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Time step %i\n",ctx.timestep);CHKERRQ(ierr);
		ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
		ierr = VecSet(fields.BCU,0.0);CHKERRQ(ierr);
		ierr = DMDAVecGetArrayDOF(ctx.daVect,fields.BCU,&bcu_array);CHKERRQ(ierr);
		switch (mode) {
			case 1:
				ierr                = PetscPrintf(PETSC_COMM_WORLD,"Applying normal displacement boundary condition on faces X0 X1 to simulate crack opening: Mode I \n");CHKERRQ(ierr);
				ctx.bcU[0].face[X0] = FIXED;ctx.bcU[0].face[X1] = FIXED;
				ctx.bcU[1].face[X0] = ZERO;ctx.bcU[1].face[X1] = ZERO;
				ctx.bcU[2].face[X0] = ZERO;ctx.bcU[2].face[X1] = ZERO;
				
				ctx.bcU[1].face[Y0] = ZERO;ctx.bcU[1].face[Y1] = ZERO;
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						for (i = xs; i < xs+xm; i++) {
							if (i == 0) {
								bcu_array[k][j][i][0] = -ctx.timestep*bc1;
							}
							if (i == nx-1) {
								bcu_array[k][j][i][0] = ctx.timestep*bc1;
							}
						}
					}
				}
				break;
			case 2:
				ierr                = PetscPrintf(PETSC_COMM_WORLD,"Applying normal displacement boundary condition on faces X0 X1 to simulate crack opening: Mode I \n");CHKERRQ(ierr);
				ctx.bcU[0].face[X0] = FIXED;ctx.bcU[0].face[X1] = FIXED;
				ctx.bcU[1].face[X0] = ZERO;ctx.bcU[1].face[X1] = ZERO;
				ctx.bcU[2].face[X0] = ZERO;ctx.bcU[2].face[X1] = ZERO;

				ctx.bcU[2].face[Z0] = FIXED;ctx.bcU[2].face[Z1] = FIXED;
				
				ctx.bcU[1].face[Y0] = ZERO;ctx.bcU[1].face[Y1] = ZERO;
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						for (i = xs; i < xs+xm; i++) {
							if (i == 0) {
								bcu_array[k][j][i][0] = -ctx.timestep*bc1;
							}
							if (i == nx-1) {
								bcu_array[k][j][i][0] = ctx.timestep*bc1;
							}
							if (k == 0) {
								bcu_array[k][j][i][2] = -ctx.timestep*bc2;
							}
							if (k == nz-1) {
								bcu_array[k][j][i][2] = ctx.timestep*bc2;
							}
						}
					}
				}
				break;
			case 3:
				ierr                = PetscPrintf(PETSC_COMM_WORLD,"Applying tangential displacement boundary condition on faces X0 X1 to simulate in-plane shear: Mode II \n");CHKERRQ(ierr);
				ctx.bcU[0].face[X0] = ZERO;ctx.bcU[0].face[X1] = ZERO;
				ctx.bcU[1].face[X0] = ZERO;ctx.bcU[1].face[X1] = ZERO;
				ctx.bcU[2].face[X0] = FIXED;ctx.bcU[2].face[X1] = FIXED;
				
				ctx.bcU[1].face[Y0] = ZERO;ctx.bcU[1].face[Y1] = ZERO;
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						for (i = xs; i < xs+xm; i++) {
							if (i == 0) {
								bcu_array[k][j][i][2] = -ctx.timestep*bc2;
							}
							if (i == nx-1) {
								bcu_array[k][j][i][2] = ctx.timestep*bc2;
							}
						}
					}
				}
				break;
			case 4:
				ierr                = PetscPrintf(PETSC_COMM_WORLD,"Applying displacement boundary condition on faces X0 X1 to simulate mixed mode: Mode I & II\n");CHKERRQ(ierr);
				ctx.bcU[0].face[X0] = FIXED;ctx.bcU[0].face[X1] = FIXED;
				ctx.bcU[1].face[X0] = ZERO;ctx.bcU[1].face[X1] = ZERO;
				ctx.bcU[2].face[X0] = FIXED;ctx.bcU[2].face[X1] = FIXED;
				
				ctx.bcU[1].face[Y0] = ZERO;ctx.bcU[1].face[Y1] = ZERO;
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						for (i = xs; i < xs+xm; i++) {
							if (i == 0) {
								bcu_array[k][j][i][0] = -ctx.timestep*bc1;
								bcu_array[k][j][i][2] = -ctx.timestep*bc2;
							}
							if (i == nx-1) {
								bcu_array[k][j][i][0] = ctx.timestep*bc1;
								bcu_array[k][j][i][2] = ctx.timestep*bc2;
							}
						}
					}
				}
				break;
			default:
				SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: mode should be between 1 & 3 got %i\n",mode);
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
		ierr              = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,fields.U,&ctx);CHKERRQ(ierr);
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
    ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);

		altminit = 1;
//		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n###################################################################\n");CHKERRQ(ierr);
//		ierr = PetscPrintf(PETSC_COMM_WORLD,"#        VF crack volume change = %f\t      \n",ctx.CrackVolume);CHKERRQ(ierr);
//		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n###################################################################\n");CHKERRQ(ierr);
	}
	PetscViewerFlush(viewer);CHKERRQ(ierr);
	PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

