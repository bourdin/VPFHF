/*
 test1.c: solves for the displacement and v-field in a pressurized penny crack in 3d (Sneddon 3D)
 (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFV.h"
#include "VFU.h"
#include "VFFlow.h"

VFCtx    ctx;
VFFields fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	VFCtx               ctx;
	VFFields            fields;
	PetscErrorCode      ierr;
	PetscReal           length = .5;
	PetscReal           height = .1;
	PetscReal           center[3]={0.,0.,.5};
	PetscInt            orientation=3;
	PetscInt            nopts=3;
	PetscInt			ek,ej,ei,c;
	PetscInt            i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
	PetscReal			****coords_array;
	PetscReal			***v_array;  
	PetscReal           BBmin[3],BBmax[3];
	PetscReal           x,y,z;  
	PetscReal           ElasticEnergy = 0;
	PetscReal           InsituWork = 0;
	PetscReal           SurfaceEnergy = 0;
	char                filename[FILENAME_MAX];
	PetscReal           p = 1.;
	PetscReal           ***pmult_array;  
	PetscReal           ****vfperm_array;
	PetscInt			altminit=1;
	Vec					Vold;
	PetscReal			errV=1e+10;
	PetscReal			q = 2.67e-4;
	PetscReal			****U_array;
	PetscReal			***FVnodof_array;
	PetscReal           ***pressure_array;  
	
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-length",&length,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-height",&height,PETSC_NULL);CHKERRQ(ierr);
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
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a penny-shaped crack of length %g and height %g at (%g,%g,%g) with normal vector <0,1,0>\n",
							   length,height,center[0],center[1],center[2]);CHKERRQ(ierr);	  
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcU[1].face[X0]= ZERO;
			/*	face X1	*/
			ctx.bcU[0].face[X1]= ZERO;
			ctx.bcU[1].face[X1]= ZERO;
			ctx.bcU[2].face[X1]= ZERO;
			/*	face Y0	*/
			ctx.bcU[1].face[Y0]= ZERO;
			/*	face Y1	*/
			ctx.bcU[1].face[Y1]= ZERO;		  
			/*	face Z0	*/
			ctx.bcU[0].face[Z0]= ZERO;
			ctx.bcU[1].face[Z0]= ZERO;
			ctx.bcU[2].face[Z0]= ZERO;
			/*	face Z1	*/
			ctx.bcU[0].face[Z1]= ZERO;
			ctx.bcU[1].face[Z1]= ZERO;
			ctx.bcU[2].face[Z1]= ZERO;
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) { 
						x = coords_array[k][j][i][0];
						y = coords_array[k][j][i][1];					  
						if ( ((k == nz/2) || (k == nz/2-1)) && (x*x+(y-0.5)*(y-0.5)) <= height*height ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		case 2:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a penny-shaped crack of length %g and height %g at (%g,%g,%g) with normal vector <0,1,0>\n",
							   length,height,center[0],center[1],center[2]);CHKERRQ(ierr);	  
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcU[1].face[X0]= ZERO;
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
						x = coords_array[k][j][i][0];
						y = coords_array[k][j][i][1];					  
						if ( ((k == nz/2) || (k == nz/2-1)) && (x*x+(y-0.5)*(y-0.5)) <= height*height ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		case 3:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a penny-shaped crack of length %g and height %g at (%g,%g,%g) with normal vector <0,1,0>\n",
							   length,height,center[0],center[1],center[2]);CHKERRQ(ierr);	  
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcU[1].face[X0]= ZERO;
			/*	face X1	*/
			ctx.bcU[0].face[X1]= ZERO;
			ctx.bcU[1].face[X1]= ZERO;
			ctx.bcU[2].face[X1]= ZERO;
			/*	face Y0	*/
			ctx.bcU[1].face[Y0]= ZERO;
			/*	face Y1	*/
			ctx.bcU[1].face[Y1]= ZERO;		  
			/*	face Z0	*/
			ctx.bcU[2].face[Z0]= ZERO;
			/*	face Z1	*/
			ctx.bcU[2].face[Z1]= ZERO;
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) { 
						x = coords_array[k][j][i][0];
						y = coords_array[k][j][i][1];					  
						if ( ((k == nz/2) || (k == nz/2-1)) && (x*x+(y-0.5)*(y-0.5)) <= height*height ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		case 4:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a penny-shaped crack of length %g and height %g at (%g,%g,%g) with normal vector <0,1,0>\n",
							   length,height,center[0],center[1],center[2]);CHKERRQ(ierr);	  
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcU[1].face[X0]= ZERO;
			/*	face X1	*/
			ctx.bcU[0].face[X1]= ZERO;
			/*	face Y0	*/
			ctx.bcU[1].face[Y0]= ZERO;
			/*	face Y1	*/
			ctx.bcU[1].face[Y1]= ZERO;		  
			/*	face Z0	*/
			ctx.bcU[2].face[Z0]= ZERO;
			/*	face Z1	*/
			ctx.bcU[2].face[Z1]= ZERO;
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) { 
						x = coords_array[k][j][i][0];
						y = coords_array[k][j][i][1];					  
						if ( ((k == nz/2) || (k == nz/2-1)) && (x*x+(y-0.5)*(y-0.5)) <= height*height ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		case 5:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a penny-shaped crack of length %g and height %g at (%g,%g,%g) with normal vector <0,1,0>\n",
							   length,height,center[0],center[1],center[2]);CHKERRQ(ierr);	  
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcU[1].face[X0]= ZERO;
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
			ctx.bcU[0].face[Z0]= ZERO;
			ctx.bcU[1].face[Z0]= ZERO;
			ctx.bcU[2].face[Z0]= ZERO;
			/*	face Z1	*/
			ctx.bcU[0].face[Z1]= ZERO;
			ctx.bcU[1].face[Z1]= ZERO;
			ctx.bcU[2].face[Z1]= ZERO;
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) { 
						x = coords_array[k][j][i][0];
						y = coords_array[k][j][i][1];					  
						if ( ((k == nz/2) || (k == nz/2-1)) && (x*x+(y-0.5)*(y-0.5)) <= height*height ) {
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

	ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
	ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);

	ctx.timevalue = 0;
	q = 2.e-3;
	ctx.maxtimestep = 20;
	for (ctx.timestep = 1; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
		p = 1.;
	do {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Time step %i, alt min step %i with pressure %g\n",ctx.timestep,altminit, p);CHKERRQ(ierr);

		ierr = VecSet(fields.pressure,1.);CHKERRQ(ierr);
		ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);

		ierr = VolumetricCrackOpening(&ctx.CrackVolume, &ctx, &fields);CHKERRQ(ierr);   
		p = q*ctx.timestep/ctx.CrackVolume;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Print pressure ................... = %g\n", p);CHKERRQ(ierr);
		ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
		
		
		ierr = VecScale(fields.U,p);CHKERRQ(ierr);
		ierr = VecSet(fields.pressure,0);CHKERRQ(ierr);
		/*
		ierr = DMDAVecGetArrayDOF(ctx.daVect,fields.U,&U_array);CHKERRQ(ierr);
		for (ek = zs; ek < zs+zm; ek++) {
			for (ej = ys; ej < ys+ym; ej++) {
				for (ei = xs; ei < xs+xm; ei++) {
					for (c = 0; c < 3; c++) {
						U_array[ek][ej][ei][c] = U_array[ek][ej][ei][c] * p;
					}
				}
			}
		}
		ierr = DMDAVecRestoreArrayDOF(ctx.daVect,fields.U,&U_array);CHKERRQ(ierr);	
		ierr = DMDAVecGetArray(ctx.daScal,fields.pressure,&pressure_array);CHKERRQ(ierr);
		for (ek = zs; ek < zs+zm; ek++) {
			for (ej = ys; ej < ys+ym; ej++) {
				for (ei = xs; ei < xs+xm; ei++) {
						pressure_array[ek][ej][ei] = p;
				}
			}
		}
		ierr = DMDAVecRestoreArray(ctx.daScal,fields.pressure,&pressure_array);CHKERRQ(ierr);
		*/
		ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);
		ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
		ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"   Max. change on V: %e\n",errV);CHKERRQ(ierr);
		altminit++;
	} while (errV > ctx.altmintol && altminit <= ctx.altminmaxit);
		switch (ctx.fileformat) {
			case FILEFORMAT_HDF5:       
				ierr = FieldsH5Write(&ctx,&fields);
				break;
			case FILEFORMAT_BIN:
				ierr = FieldsBinaryWrite(&ctx,&fields);
				break; 
		}
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n Final Crack volume\t = %g, Pressure\t= %g\n\n", p*ctx.CrackVolume, p);CHKERRQ(ierr);	
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
		altminit = 0.;
	}
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}
