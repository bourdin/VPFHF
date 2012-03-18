/*
 test10.c: 
 Validate crack opening computation
 
 (c) 2010-2011 chukwudi Chukwudozie cchukw1@tigers.lsu.edu
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
	VFCtx               ctx;
	VFFields            fields;
	PetscErrorCode      ierr;
	
	PetscReal           length = .2;
	PetscInt            orientation=0;
	PetscInt            nopts=3;
	PetscInt            i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
	PetscReal			****coords_array;
	PetscReal			****bcu_array;  
	PetscReal           BBmin[3],BBmax[3];
	PetscReal           ElasticEnergy = 0;
	PetscReal           InsituWork = 0;
	PetscReal           SurfaceEnergy = 0;
	char                filename[FILENAME_MAX];
	PetscReal           bc = .2;
	PetscReal           ***v_array;  
	PetscReal           lx,ly,lz;
	PetscReal           ***pmult_array;  
	PetscReal           ****vfperm_array;  
	
	
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	
	ierr = PetscOptionsGetInt(PETSC_NULL,"-orientation",&orientation,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
	
	lz = BBmax[2];
	ly = BBmax[1];
	lx = BBmax[0];
    
    printf("\nBounding box: lx = %f\tly = %f\tlz = %f\n", lx, ly, lz);
    
	ctx.matprop[0].beta  = 0.;
	ctx.matprop[0].alpha = 0.;
	
	ctx.timestep  = 1;
	ctx.timevalue = 1.;
	
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VecSet(fields.V,1.0);CHKERRQ(ierr);
	ierr = VecSet(fields.VIrrev,1.0);CHKERRQ(ierr);
	ierr = VecSet(fields.U,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.pressure,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
	
	ierr = DMDAVecGetArray(ctx.daScal,fields.VIrrev,&v_array);CHKERRQ(ierr);    
	
	ierr = VecSet(fields.BCU,0.0);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVect,fields.BCU,&bcu_array);CHKERRQ(ierr);    
	
	/*
	 Reset all BC for U and V
	 */
	for (i = 0; i < 6; i++) {
		ctx.bcV[0].face[i]=NONE;
		for (j = 1; j < 3; j++) {
			ctx.bcU[j].face[i] = ZERO;
		}
	}
	for (i = 0; i < 12; i++) {
		ctx.bcV[0].edge[i]=NONE;
		for (j = 1; j < 3; j++) {
			ctx.bcU[j].edge[i] = ZERO;
		}
	}
	for (i = 0; i < 8; i++) {
		ctx.bcV[0].vertex[i]=NONE;
		for (j = 1; j < 3; j++) {
			ctx.bcU[j].vertex[i] = ZERO;
		}
	}
	switch (orientation) {
		case 0:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Applying traction Dirichlet conditions on faces X0 X1\n");CHKERRQ(ierr);     
			ctx.bcU[0].face[X0]=FIXED;        ctx.bcU[0].face[X1]=FIXED;
			ctx.bcU[1].face[X0]=ZERO;        ctx.bcU[1].face[X1]=ZERO;
			ctx.bcU[2].face[X0]=ZERO;        ctx.bcU[2].face[X1]=ZERO;
			
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) { 
						if (((i == nx/2)) || (i == nx/2-1)) {
							v_array[k][j][i] = 0.;
						}
						if (i == 0) {
							bcu_array[k][j][i][0] = -bc;
						}
						if (i == nx-1) {
							bcu_array[k][j][i][0] = bc;
						}
					}
				}
			}      
			break;
		default:
			SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Orientation should be between 0 and 2, got %i\n",orientation);
			break;
	}  
	ierr = DMDAVecRestoreArray(ctx.daScal,fields.VIrrev,&v_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,fields.BCU,&bcu_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	
	ierr = VecCopy(fields.VIrrev,fields.V);CHKERRQ(ierr);
	
	ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
	
	ierr = VF_StepU(&fields,&ctx);
	ctx.hasCrackPressure = PETSC_FALSE;
	ierr = VF_StepV(&fields,&ctx);
	
	ctx.ElasticEnergy=0;
	ctx.InsituWork=0;
	ctx.PressureWork = 0.;
	ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,&fields,&ctx);CHKERRQ(ierr);
    ierr = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);
	
	ctx.TotalEnergy = ctx.ElasticEnergy - ctx.InsituWork - ctx.PressureWork;
	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface energy:            %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic Energy:            %e\n",ctx.ElasticEnergy);CHKERRQ(ierr);
	if (ctx.hasCrackPressure) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of pressure forces:   %e\n",ctx.PressureWork);CHKERRQ(ierr);
	}
	if (ctx.hasInsitu) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of surface forces:    %e\n",ctx.InsituWork);CHKERRQ(ierr);
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Total Mechanical energy:  %e\n",ctx.ElasticEnergy-InsituWork-ctx.PressureWork);CHKERRQ(ierr);
	
	ierr = CrackOpeningDisplacement(&ctx, &fields);CHKERRQ(ierr);
	
	ierr = DMDAVecGetArray(ctx.daScal,fields.pmult,&pmult_array);CHKERRQ(ierr);    
	ierr = DMDAVecGetArrayDOF(ctx.daVFperm,fields.vfperm,&vfperm_array);CHKERRQ(ierr);    
	
	
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) { 
				pmult_array[k][j][i] = vfperm_array[k][j][i][0];
			}
		}
	} 
	
	ierr = DMDAVecRestoreArray(ctx.daScal,fields.pmult,&pmult_array);CHKERRQ(ierr);    
	ierr = DMDAVecRestoreArrayDOF(ctx.daVFperm,fields.vfperm,&vfperm_array);CHKERRQ(ierr);

	ierr = CellToNodeInterpolation(ctx.daScal, fields.FVCellndof, fields.pmult, &ctx); CHKERRQ(ierr);

	
	/*
	 Save fields and write statistics about current run
	 */    
	switch (ctx.fileformat) {
		case FILEFORMAT_HDF5:       
			ierr = FieldsH5Write(&ctx,&fields);
			ierr = FieldsH5Write(&ctx,&fields);
			break;
		case FILEFORMAT_BIN:
			ierr = FieldsBinaryWrite(&ctx,&fields);
			break; 
	} 
    printf("#        Actual crack volume change = %f\t      \n\n", (lz*ly*0.4) );
    printf("\n###################################################################\n\n\n");
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

