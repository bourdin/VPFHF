/*
 test31.c: 3D Flow problem with source term [pressure = sin(2*pi*x)*sin(2*pi*y)*sin(2(pi*z)]. All velocity boundary condition
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
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
	PetscErrorCode  ierr;
	PetscViewer		viewer;
	PetscViewer     logviewer;
	char			filename[FILENAME_MAX];
	PetscInt		i,j,k,c,nx,ny,nz,xs,xm,ys,ym,zs,zm;
	PetscReal		BBmin[3],BBmax[3];
	PetscReal		****flowbc_array;
	PetscReal		***src_array;
	PetscReal		****coords_array;
	PetscReal		hx,hy,hz;
	PetscReal		gx,gy,gz;
	PetscReal		gamma, beta, rho, mu;
	/*	Fracture mechanics part */
	PetscInt        orientation=1;
	PetscReal       length = .1;
	PetscReal		***v_array;  
	Vec				Vold;
	PetscReal       p;
	PetscReal		****perm_array;
	PetscReal			  lx,ly,lz;

	
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = FlowSolverInitialize(&ctx,&fields);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	ierr = VecSet(fields.FlowBCArray,0.);CHKERRQ(ierr);
	ierr = VecSet(ctx.Source,0.);CHKERRQ(ierr);

	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);	

	/*	Fracture mechanics part */
	lz = BBmax[2]-BBmin[2];
	ly = BBmax[1]-BBmin[1];
	lx = BBmax[0]-BBmin[0];	
	
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
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
				//			ctx.bcU[1].face[X0]= ZERO;
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
						if ( ((j == ny/2) || (j == ny/2-1) ) && (coords_array[k][j][i][0] > lx/2.-length) && (coords_array[k][j][i][0] < lx/2.+length ) ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		default:
			SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Orientation should be one of {1,2,3,4,5,6}, got %i\n",orientation);
			break;
	}  
	ierr = DMDAVecRestoreArray(ctx.daScal,fields.VIrrev,&v_array);CHKERRQ(ierr);
	ierr = VecCopy(fields.VIrrev,fields.V);CHKERRQ(ierr);

	
	
	
	
	hx = 1./(nx-1);
	hy = 1./(nx-1);
	hz = 1./(nz-1);	
	rho = ctx.flowprop.rho;									 
	mu = ctx.flowprop.mu;     
	beta = ctx.flowprop.beta;		
	gamma = ctx.flowprop.gamma;									
    gx = ctx.flowprop.g[0];
    gy = ctx.flowprop.g[1];
    gz = ctx.flowprop.g[2];
	/*
	 Reset all Flow BC for velocity and P
	 */
	for (i = 0; i < 6; i++) {
		for (c = 0; c < 4; c++) {
			ctx.bcFlow[c].face[i] = NOBC;
		}
	}
	for (i = 0; i < 12; i++) {
		for (c = 0; c < 4; c++) {
			ctx.bcFlow[c].edge[i] = NOBC;
		}
	}
	for (i = 0; i < 8; i++) {
		for (c = 0; c < 4; c++) {
			ctx.bcFlow[c].vertex[i] = NOBC;
		}
	}
	ctx.bcFlow[3].face[X0] = PRESSURE;
	ctx.bcFlow[3].face[X1] = PRESSURE;
	ctx.bcFlow[3].face[Y0] = PRESSURE;
	ctx.bcFlow[3].face[Y1] = PRESSURE;
	ctx.bcFlow[2].face[Z0] = VELOCITY;
	ctx.bcFlow[2].face[Z1] = VELOCITY;	
	/* 
	 Now done with all initializations
	 */
	ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
	ctx.hasCrackPressure = PETSC_TRUE;
	ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
	p = 1;
	ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
	ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
	ctx.timevalue = 0;
	ctx.maxtimestep = 1;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-maxtimestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-length",&length,PETSC_NULL);CHKERRQ(ierr);

	for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nProcessing step %i.\n",ctx.timestep);CHKERRQ(ierr);
		ctx.timevalue = ctx.timestep * ctx.maxtimevalue / (ctx.maxtimestep-1.);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\ntime value %f \n",ctx.timevalue);CHKERRQ(ierr);

		ierr = DMDAVecGetArrayDOF(ctx.daFlow,fields.FlowBCArray,&flowbc_array);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(ctx.daScal,ctx.Source,&src_array);CHKERRQ(ierr);
		for (k = zs; k < zs+zm; k++) {
			for (j = ys; j < ys+ym; j++) {
				for (i = xs; i < xs+xm; i++) {
					if( (i == nx/2) && ((j == ny/2) || (j == ny/2-1) )){
//					if( (i == nx/2) && ((j == ny/2) )){
						src_array[k][j][i] = 0.01 * (ctx.timestep + 1);
					}
				}
			}
		}
		ierr = DMDAVecRestoreArray(ctx.daScal,ctx.Source,&src_array);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArrayDOF(ctx.daFlow,fields.FlowBCArray,&flowbc_array);CHKERRQ(ierr);

		/*
		 Do mechanics step 
		 */
		ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
		ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
		ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr);
		ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);
		ierr = VolumetricCrackOpening(&ctx.CrackVolume, &ctx, &fields);CHKERRQ(ierr);   

		
		ierr = VecSet(fields.vfperm,0.);CHKERRQ(ierr);	
		ierr = DMDAVecGetArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);
		for (k = zs; k < zs + zm; k++) {
			for (j = ys; j < ys+ym; j++) {
				for (i = xs; i < xs+xm; i++) {
					for(c = 0; c < 2; c++){
						perm_array[k][j][i][c] = 2.e-5;
					}				 
				}
			}
		}
		ierr = DMDAVecRestoreArrayDOF(ctx.daVFperm,fields.vfperm,&perm_array);CHKERRQ(ierr);
		
		ierr = PostProcessing(&ctx,&fields);CHKERRQ(ierr);

		/*
		 Do flow solver step 
		 */
		ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
		
		
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
		ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.log",ctx.prefix);CHKERRQ(ierr);
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&logviewer);CHKERRQ(ierr);
		ierr = PetscLogView(logviewer);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&logviewer);
	}
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = FlowSolverFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

