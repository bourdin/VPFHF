/*
 test13.c: solves for the displacement and v-field in a pressurized line crack in 2d (Sneddon 2D)
 (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFV.h"
#include "VFU.h"
#include "VFFlow.h"

VFCtx               ctx;
VFFields            fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode      ierr;
	PetscViewer			viewer;
	PetscViewer         logviewer;
	PetscReal           length = .3;
	PetscReal           length1 = .1;
	PetscReal           center[3]={0.,0.,.5};
	PetscInt            orientation=2;
	PetscInt            nopts=3;
	PetscInt            i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
	PetscReal			****coords_array;
	PetscReal			***v_array;  
	PetscReal           BBmin[3],BBmax[3];
	PetscReal           x,y,z;  
	PetscReal           ElasticEnergy = 0;
	PetscReal           InsituWork = 0;
	PetscReal           SurfaceEnergy = 0;
	char                filename[FILENAME_MAX];
	PetscReal           p;
	PetscReal           p_old;
	PetscReal           p_epsilon = 1.e-5;
	PetscInt			altminit=1;
	Vec					Vold;
	PetscReal			errV=1e+10;
	PetscReal			q=2.e-4;
	
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-length",&length,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-q",&q,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetRealArray(PETSC_NULL,"-center",&center[0],&nopts,PETSC_NULL);CHKERRQ(ierr);
	
	ierr = PetscOptionsGetInt(PETSC_NULL,"-orientation",&orientation,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daScal,BBmin,BBmax);CHKERRQ(ierr);
	
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VecSet(fields.VIrrev,1.0);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx.daScal,fields.VIrrev,&v_array);CHKERRQ(ierr);    
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
	switch (orientation) {
		case 1:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a penny-shaped crack of length %g at (%g,%g,%g) with normal vector <0,1,0>\n",
							   length,center[0],center[1],center[2]);CHKERRQ(ierr);		  
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
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) { 
						if ( ((j == 3*ny/4) || (j == 3*ny/4-1)) && PetscAbs(coords_array[k][j][i][0]-(BBmin[0]+BBmax[0])/2.) <= length ) {
							v_array[k][j][i] = 0.;
						}
						if ( ((i == nx/2) || (i == nx/2-1)) && PetscAbs(coords_array[k][j][i][1]-(BBmin[1]+BBmax[1])/2.) <= length1 ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		case 2:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a penny-shaped crack of length %g at (%g,%g,%g) with normal vector <0,1,0>\n",
							   length,center[0],center[1],center[2]);CHKERRQ(ierr);		  
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcU[2].face[X0]= ZERO;
			/*	face X1	*/
			ctx.bcU[0].face[X1]= ZERO;
			ctx.bcU[2].face[X1]= ZERO;
			/*	face Y0	*/
			ctx.bcU[1].face[Y0]= ZERO;
			ctx.bcU[2].face[Y0]= ZERO;
			/*	face Y1	*/
			ctx.bcU[1].face[Y1]= ZERO;
			ctx.bcU[2].face[Y1]= ZERO;		  
			/*	face Z0	*/
			ctx.bcU[2].face[Z0]= ZERO;
			/*	face Z1	*/
			ctx.bcU[2].face[Z1]= ZERO;
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) { 
						if ( ((j == 3*ny/4) || (j == 3*ny/4-1)) && PetscAbs(coords_array[k][j][i][0]-(BBmin[0]+BBmax[0])/2.) <= length ) {
							v_array[k][j][i] = 0.;
						}
						if ( ((i == nx/2) || (i == nx/2-1)) && PetscAbs(coords_array[k][j][i][1]-(BBmin[1]+BBmax[1])/2.) <= length1 ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		case 3:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a penny-shaped crack of length %g at (%g,%g,%g) with normal vector <0,1,0>\n",
							   length,center[0],center[1],center[2]);CHKERRQ(ierr);		  
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcU[1].face[X0]= ZERO;
			ctx.bcU[2].face[X0]= ZERO;
			/*	face X1	*/
			ctx.bcU[0].face[X1]= ZERO;
			ctx.bcU[1].face[X1]= ZERO;
			ctx.bcU[2].face[X1]= ZERO;
			/*	face Y0	*/
			ctx.bcU[1].face[Y0]= ZERO;
			ctx.bcU[2].face[Y0]= ZERO;
			/*	face Y1	*/
			ctx.bcU[1].face[Y1]= ZERO;		  
			ctx.bcU[2].face[Y1]= ZERO;		  
			/*	face Z0	*/
			ctx.bcU[2].face[Z0]= ZERO;
			/*	face Z1	*/
			ctx.bcU[2].face[Z1]= ZERO;
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) { 
						if ( ((j == 3*ny/4) || (j == 3*ny/4-1)) && PetscAbs(coords_array[k][j][i][0]-(BBmin[0]+BBmax[0])/2.) <= length ) {
							v_array[k][j][i] = 0.;
						}
						if ( ((i == nx/2) || (i == nx/2-1)) && PetscAbs(coords_array[k][j][i][1]-(BBmin[1]+BBmax[1])/2.) <= length1 ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		case 4:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a penny-shaped crack of length %g at (%g,%g,%g) with normal vector <0,1,0>\n",
							   length,center[0],center[1],center[2]);CHKERRQ(ierr);		  
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcU[2].face[X0]= ZERO;
			/*	face X1	*/
			ctx.bcU[0].face[X1]= ZERO;
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
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) { 
						if ( ((j == 3*ny/4) || (j == 3*ny/4-1)) && PetscAbs(coords_array[k][j][i][0]-(BBmin[0]+BBmax[0])/2.) <= length ) {
							v_array[k][j][i] = 0.;
						}
						if ( ((i == nx/2) || (i == nx/2-1)) && PetscAbs(coords_array[k][j][i][1]-(BBmin[1]+BBmax[1])/2.) <= length1 ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		case 5:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a penny-shaped crack of length %g at (%g,%g,%g) with normal vector <0,1,0>\n",
							   length,center[0],center[1],center[2]);CHKERRQ(ierr);		  
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcU[2].face[X0]= ZERO;
			/*	face X1	*/
			ctx.bcU[0].face[X1]= ZERO;
			ctx.bcU[2].face[X1]= ZERO;
			/*	face Y0	*/
				/*.....FREE.......*/
			/*	face Y1	*/
				/*.....FREE.......*/
			/*	face Z0	*/
			ctx.bcU[2].face[Z0]= ZERO;
			/*	face Z1	*/
			ctx.bcU[2].face[Z1]= ZERO;
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) { 
						if ( ((j == 3*ny/4) || (j == 3*ny/4-1)) && PetscAbs(coords_array[k][j][i][0]-(BBmin[0]+BBmax[0])/2.) <= length ) {
							v_array[k][j][i] = 0.;
						}
						if ( ((i == nx/2) || (i == nx/2-1)) && PetscAbs(coords_array[k][j][i][1]-(BBmin[1]+BBmax[1])/2.) <= length1 ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		case 6:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a penny-shaped crack of length %g at (%g,%g,%g) with normal vector <0,1,0>\n",
							   length,center[0],center[1],center[2]);CHKERRQ(ierr);		  
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcU[1].face[X0]= ZERO;
			ctx.bcU[2].face[X0]= ZERO;
			/*	face X1	*/
			ctx.bcU[0].face[X1]= ZERO;
			ctx.bcU[1].face[X1]= ZERO;
			ctx.bcU[2].face[X1]= ZERO;
			/*	face Y0	*/
			/*.....FREE.......*/
			/*	face Y1	*/
			/*.....FREE.......*/
			/*	face Z0	*/
			ctx.bcU[2].face[Z0]= ZERO;
			/*	face Z1	*/
			ctx.bcU[2].face[Z1]= ZERO;
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) { 
						if ( ((j == 3*ny/4) || (j == 3*ny/4-1)) && PetscAbs(coords_array[k][j][i][0]-(BBmin[0]+BBmax[0])/2.) <= length ) {
							v_array[k][j][i] = 0.;
						}
						if ( ((i == nx/2) || (i == nx/2-1)) && PetscAbs(coords_array[k][j][i][1]-(BBmin[1]+BBmax[1])/2.) <= length1 ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		default:
			SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Orientation should be one of {1,2,3}, got %i\n",orientation);
			break;
	}  
	ierr = DMDAVecRestoreArray(ctx.daScal,fields.VIrrev,&v_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VecCopy(fields.VIrrev,fields.V);CHKERRQ(ierr);
	ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
	
	ctx.hasCrackPressure = PETSC_TRUE;
	ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
	

	
	PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
	PetscViewerSetType(viewer, PETSCVIEWERASCII);
	PetscViewerFileSetMode(viewer, FILE_MODE_APPEND);
	PetscViewerFileSetName(viewer, "pressure.txt");
	PetscViewerASCIIPrintf(viewer, "Time step \t Volume \t Pressure \t SurfaceEnergy \t ElasticEnergy \t PressureForces \t TotalMechEnergy \n");
	
	p = 1.;
	ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
	ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
	ctx.timevalue = 0;
	ctx.maxtimestep = 150;
	for (ctx.timestep = 1; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
		do {
			p_old = p;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Time step %i, alt min step %i with pressure %g\n",ctx.timestep,altminit, p);CHKERRQ(ierr);
			ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
			ierr = VecScale(fields.U,1./p);CHKERRQ(ierr);
			ierr = VolumetricCrackOpening(&ctx.CrackVolume, &ctx, &fields);CHKERRQ(ierr);   
			p = q*ctx.timestep/ctx.CrackVolume;
			ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
			ierr = VecScale(fields.U,p);CHKERRQ(ierr);
			ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
			ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);

			ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
			ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,"   Max. change on V: %e\n",errV);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,"   Max. change on p: %e\n", PetscAbs(p-p_old));CHKERRQ(ierr);
			altminit++;
		} while (PetscAbs(p-p_old) >= p_epsilon && altminit <= ctx.altminmaxit);
		ierr = VolumetricCrackOpening(&ctx.CrackVolume, &ctx, &fields);CHKERRQ(ierr);   
		switch (ctx.fileformat) {
			case FILEFORMAT_HDF5:       
				ierr = FieldsH5Write(&ctx,&fields);
				break;
			case FILEFORMAT_BIN:
				ierr = FieldsBinaryWrite(&ctx,&fields);
				break; 
		}
		ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.log",ctx.prefix);CHKERRQ(ierr);
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&logviewer);CHKERRQ(ierr);
		ierr = PetscLogView(logviewer);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&logviewer);
		ctx.ElasticEnergy=0;
		ctx.InsituWork=0;
		ctx.PressureWork = 0.;
		ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,&fields,&ctx);CHKERRQ(ierr);
		ierr = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);
		ctx.TotalEnergy = ctx.ElasticEnergy - ctx.InsituWork - ctx.PressureWork;
		PetscViewerASCIIPrintf(viewer, "%d \t %g \t %g \t %g \t %g \t %g \t %g\n", ctx.timestep , ctx.CrackVolume, p, ctx.SurfaceEnergy, ctx.ElasticEnergy, ctx.PressureWork, ctx.TotalEnergy);
		altminit = 0.;
	}
	PetscViewerFlush(viewer);CHKERRQ(ierr);
	PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	PetscViewerDestroy(&logviewer);CHKERRQ(ierr);
	ierr = VecDestroy(&Vold);CHKERRQ(ierr);

	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}
