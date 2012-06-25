/*
 test26.c: 3D:Solves for displacement and fracture given random distribution of seed points (i.e v = 0. at random points)
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFCracks.h"
#include "VFV.h"
#include "VFU.h"
#include "VFFlow.h"
#include "time.h"

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
	PetscInt     i,ii,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm;
	PetscReal    ****coords_array;
	PetscReal    ElasticEnergy = 0;
	PetscReal    InsituWork    = 0;
	PetscReal    SurfaceEnergy = 0;
	char         filename[FILENAME_MAX];
	PetscReal    p = 1e-3;
	VFPennyCrack *crack;
	VFRectangularCrack *rcrack;
	PetscInt     npc = 0;
	PetscInt     nrc = 0;
	char         prefix[PETSC_MAX_PATH_LEN+1];
	Vec          V;
	PetscReal bc0 = 0.05;
	PetscReal bc1 = 0.01;
	PetscReal bc2 = 0.02;
	PetscInt			altminit=1;
	Vec					Vold;
	PetscReal			errV=1e+10;
	PetscReal ****bcu_array;
	PetscInt  mode = 1;
	PetscReal BBmin[3],BBmax[3];
	PetscInt	no_seed = 10;
	PetscInt *seed_array_z;
	PetscInt *seed_array_y;
	PetscInt *seed_array_x;
	PetscInt	seedz_start;
	PetscInt	seedy_start;
	PetscInt	seedx_start;
	PetscReal ***v_array;
	
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-no_seed",&no_seed,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-mode",&mode,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-bcUx",&bc0,PETSC_NULL);CHKERRQ(ierr); 
	ierr = PetscOptionsGetReal(PETSC_NULL,"-bcUy",&bc1,PETSC_NULL);CHKERRQ(ierr); 
	ierr = PetscOptionsGetReal(PETSC_NULL,"-bcUz",&bc2,PETSC_NULL);CHKERRQ(ierr); 
	ierr = PetscOptionsGetInt(PETSC_NULL,"-npc",&npc,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-nrc",&nrc,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);	
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx.daScal,&V);CHKERRQ(ierr);
	ierr = VecSet(V,1.0);CHKERRQ(ierr);
	ierr = VecSet(fields.V,1.0);CHKERRQ(ierr);
	ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
	ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
	ierr = PetscMalloc(npc*sizeof(VFPennyCrack),&crack);CHKERRQ(ierr);
	ierr = PetscMalloc(nrc*sizeof(VFRectangularCrack),&rcrack);CHKERRQ(ierr);
/*nitialize and build penny shaped crack*/
	for (i = 0; i < npc; i++) {
		ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN,"cp%d_",i);CHKERRQ(ierr);
		ierr = VFPennyCrackCreate(&crack[i]);CHKERRQ(ierr);
		ierr = VFPennyCrackGet(prefix,&crack[i]);CHKERRQ(ierr);
		ierr = VFPennyCrackView(&crack[i],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	}
	for (c = 0; c < npc; c++) {
		ierr = VFPennyCrackBuildVAT2(V,&crack[c],&ctx);CHKERRQ(ierr);
		ierr = VecPointwiseMin(fields.V,V,fields.V);CHKERRQ(ierr);
	}
/*nitialize and build rectanular shaped crack*/
	for (i = 0; i < nrc; i++) {
		ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN,"cr%d_",i);CHKERRQ(ierr);
		ierr = VFRectangularCrackCreate(&rcrack[i]);CHKERRQ(ierr);
		ierr = VFRectangularCrackGet(prefix,&rcrack[i]);CHKERRQ(ierr);
		ierr = VFRectangularCrackView(&rcrack[i],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	}
	for (c = 0; c < nrc; c++) {
		ierr = VFRectangularCrackBuildVAT2(V,&rcrack[c],&ctx);CHKERRQ(ierr);
		ierr = VecPointwiseMin(fields.V,V,fields.V);CHKERRQ(ierr);
	}
	ierr = VecDestroy(&V);CHKERRQ(ierr);
	ierr = PetscMalloc3(no_seed,PetscInt,&seed_array_z,no_seed,PetscInt,&seed_array_y,no_seed,PetscInt,&seed_array_x);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx.daScal,fields.V,&v_array);CHKERRQ(ierr);
	
	/*Initializing the v-field*/	
	seedz_start = seedy_start = seedx_start = 1;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-seedz_start",&seedz_start,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-seedy_start",&seedy_start,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-seedx_start",&seedx_start,PETSC_NULL);CHKERRQ(ierr);

	srand((time(NULL)+1));
	for(i = 0; i < no_seed; i++){
		seed_array_z[i] = rand() % (nz-seedz_start)+seedz_start;
	}
	for(i = 0; i < no_seed; i++){
		seed_array_y[i] = rand() % (ny-seedy_start)+seedy_start;
	}
	for(i = 0; i < no_seed; i++){
		seed_array_x[i] = rand() % (nx-seedx_start)+seedx_start;
	}

	PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
	PetscViewerSetType(viewer, PETSCVIEWERASCII);
	PetscViewerFileSetMode(viewer, FILE_MODE_APPEND);
	PetscViewerFileSetName(viewer, "seed_nodes.txt");
	PetscViewerASCIIPrintf(viewer, "seed \t z_nodes \t y_nodes \t x_nodes \n");
	for(ii = 0; ii < no_seed; ii++){
		PetscViewerASCIIPrintf(viewer, "%d \t %d \t %d \t %d \n",ii, seed_array_z[ii],seed_array_y[ii],seed_array_x[ii] );
	}
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) { 
				for(ii = 0; ii < no_seed; ii++)
					if(k == seed_array_z[ii] && j == seed_array_y[ii] && i == seed_array_x[ii]){
						v_array[k][j][i] = 0.;
					}
			}
		}
	}
	ierr = PetscFree3(seed_array_z,seed_array_y,seed_array_x);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx.daScal,fields.V,&v_array);CHKERRQ(ierr);
	ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
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
	ctx.timestep  = 1;
	ctx.maxtimestep = 50;
	ctx.hasCrackPressure = PETSC_FALSE;
	for(ctx.timestep = 1; ctx.timestep < ctx.maxtimestep; ctx.timestep++){	
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Time step %i\n",ctx.timestep);CHKERRQ(ierr);
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
								bcu_array[k][j][i][0] = -ctx.timestep*bc0;
							}
							if (i == nx-1) {
								bcu_array[k][j][i][0] = ctx.timestep*bc0;
							}
						}
					}
				}
				break;
			case 2:
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
			case 3:
				ierr                = PetscPrintf(PETSC_COMM_WORLD,"Applying displacement boundary condition on faces X0 X1 to simulate mixed mode: Mode I & II\n");CHKERRQ(ierr);
				ctx.bcU[0].face[X0] = FIXED;ctx.bcU[0].face[X1] = FIXED;
				ctx.bcU[1].face[X0] = ZERO;ctx.bcU[1].face[X1] = ZERO;
				ctx.bcU[2].face[X0] = FIXED;ctx.bcU[2].face[X1] = FIXED;
				
				ctx.bcU[1].face[Y0] = ZERO;ctx.bcU[1].face[Y1] = ZERO;
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						for (i = xs; i < xs+xm; i++) {
							if (i == 0) {
								bcu_array[k][j][i][0] = -ctx.timestep*bc0;
								bcu_array[k][j][i][2] = -ctx.timestep*bc2;
							}
							if (i == nx-1) {
								bcu_array[k][j][i][0] = ctx.timestep*bc0;
								bcu_array[k][j][i][2] = ctx.timestep*bc2;
							}
						}
					}
				}
				break;
                        case 4:
                                ierr                = PetscPrintf(PETSC_COMM_WORLD,"Applying displacement boundary condition on faces X0 X1 to simulate mixed mode: Mode I & II\n");CHKERRQ(ierr);
                                ctx.bcU[0].face[X0] = FIXED;ctx.bcU[0].face[X1] = FIXED;
                                ctx.bcU[1].face[Y0] = FIXED;ctx.bcU[1].face[Y1] = FIXED;
                                ctx.bcU[2].face[Z0] = FIXED;ctx.bcU[2].face[Z1] = FIXED;

                                for (k = zs; k < zs+zm; k++) {
                                        for (j = ys; j < ys+ym; j++) {
                                                for (i = xs; i < xs+xm; i++) {
                                                        if (i == 0) {
                                                                bcu_array[k][j][i][0] = -ctx.timestep*bc0;
                                                        }
                                                        if (i == nx-1) {
                                                                bcu_array[k][j][i][0] = ctx.timestep*bc0;
                                                        }
                                                        if (j == 0) {
                                                                bcu_array[k][j][i][1] = -ctx.timestep*bc1;
                                                        }
                                                        if (j == ny-1) {
                                                                bcu_array[k][j][i][1] = ctx.timestep*bc1;
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
			default:
				SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Orientation should be between 1 & 4 got %i\n",mode);
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
	PetscViewerFlush(viewer);CHKERRQ(ierr);
	PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

