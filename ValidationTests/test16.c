/*
 test16.c: test for hydraulic fracturing in 3D with no leakoff
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
	PetscViewer			viewer;
	VFCtx               ctx;
	VFFields            fields;
	PetscErrorCode      ierr;
	PetscReal           radius = .1;
	PetscReal           center[3]={0.5,0.5,0.5};
	PetscInt            orientation=1;
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
	PetscReal           p_old = 1.;
	PetscReal           p_epsilon = 1.e-6;
	PetscInt			altminit=1;
	Vec					Vold;
	PetscReal			errV=1e+10;
	PetscReal			q = 2.67e-4;
	PetscReal           p_read, q_read;
	PetscReal           vol_inj;
	PetscReal           max_it, p_conv;

	ctx.maxtimestep = 10;
	max_it = 200;
	p_conv = 1e-6;
	
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-radius",&radius,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-max_it",&max_it,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-p_conv",&p_conv,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-pinit",&p_read,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-rate",&q_read,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetRealArray(PETSC_NULL,"-center",&center[0],&nopts,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(PETSC_NULL,"-maxtimestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr); 	
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
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a penny-shaped crack of radius %g at (%g,%g,%g) with normal vector <0,1,0>\n",
							   radius,center[0],center[1],center[2]);CHKERRQ(ierr);	  
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
			/*	face X1	*/
			ctx.bcU[0].face[X1]= ZERO;
			/*	face Y0	*/

			/*	face Y1	*/

			/*	face Z0	*/
			ctx.bcU[2].face[Z0]= ZERO;
			/*	face Z1	*/
			ctx.bcU[2].face[Z1]= ZERO;

            ctx.bcU[0].vertex[X0Y0Z1] = ZERO;		
            ctx.bcU[1].vertex[X0Y0Z1] = ZERO;	
            ctx.bcU[2].vertex[X0Y0Z1] = ZERO;	
			
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) { 
						x = coords_array[k][j][i][0];
						y = coords_array[k][j][i][1];
                        z = coords_array[k][j][i][2];						
						if ( ((j == ny/2) || (j == ny/2+1)) && ((x-center[0])*(x-center[0])+(z-center[2])*(z-center[2])) <= radius*radius ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		case 2:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a penny-shaped crack of radius %g at (%g,%g,%g) with normal vector <0,1,0>\n",
							   radius,center[0],center[1],center[2]);CHKERRQ(ierr);	  
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
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

            ctx.bcU[0].vertex[X0Y0Z1] = ZERO;		
            ctx.bcU[1].vertex[X0Y0Z1] = ZERO;	
            ctx.bcU[2].vertex[X0Y0Z1] = ZERO;	
			
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) { 
						x = coords_array[k][j][i][0];
						y = coords_array[k][j][i][1];
                        z = coords_array[k][j][i][2];						
						if ( ((j == ny/2) || (j == ny/2+1)) && ((x-center[0])*(x-center[0])+(z-center[2])*(z-center[2])) <= radius*radius ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		default:
			SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Orientation specified is not defined yet, got %i\n",orientation);
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

	PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
	PetscViewerSetType(viewer, PETSCVIEWERASCII);
	PetscViewerFileSetName(viewer, "pressure.txt");
	PetscViewerASCIIPrintf(viewer, "Time step \t Volume \t Pressure\n");

	
	
	ctx.timevalue = 0;
	q = q_read;
	vol_inj = 0;
	p = p_read;
    ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
	for (ctx.timestep = 1; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
	  vol_inj += q;
	  do {
	    beginning:
		p_old = p;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Time step %i, alt min step %i with pressure %g rate %g\n",ctx.timestep,altminit, p,q);CHKERRQ(ierr);
        ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
		ierr = VecScale(fields.U,1./p);CHKERRQ(ierr);		
		ierr = VolumetricCrackOpening(&ctx.CrackVolume, &ctx, &fields);CHKERRQ(ierr);   
		p = vol_inj/ctx.CrackVolume;
		ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
		ierr = VecScale(fields.U,p);CHKERRQ(ierr);
		ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
		ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);
		ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
		ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"   Max. change on V: %e\n",errV);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"   Max. change on p: %e\n", PetscAbs(p-p_old));CHKERRQ(ierr);
		altminit++;
		if (altminit >= 200){
		  vol_inj -= q;
		  q = 0.5 * q;
		  altminit = 0;
		  goto beginning;
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Flow rate has been cut back... restarting the iteration \n");CHKERRQ(ierr);
		}
	  } while (PetscAbs((p-p_old)/p) >= 1e-5 && altminit <= ctx.altminmaxit);
			/*	} while (errV > ctx.altmintol && altminit <= ctx.altminmaxit && PetscAbs(p-p_old) >= 1e-7);	*/
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n Final Crack volume\t = %g, Pressure\t= %g\n\n", p*ctx.CrackVolume, p);CHKERRQ(ierr);	
		ierr = VolumetricCrackOpening(&ctx.CrackVolume, &ctx, &fields);CHKERRQ(ierr);   
//		PetscViewerASCIIPrintf(viewer, "%d \t %g \t %g\n", ctx.timestep , ctx.CrackVolume, p);		
		switch (ctx.fileformat) {
			case FILEFORMAT_HDF5:       
				ierr = FieldsH5Write(&ctx,&fields);
				break;
			case FILEFORMAT_BIN:
				ierr = FieldsBinaryWrite(&ctx,&fields);
				break; 
		}
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
		
		PetscViewerASCIIPrintf(viewer, "%d \t %g \t %g \t %g \t %g \t %g \t %g\n", ctx.timestep , ctx.CrackVolume, p, ctx.InsituWork, ctx.ElasticEnergy, ctx.PressureWork, ctx.TotalEnergy);

		altminit = 0.;
	}
	PetscViewerFlush(viewer);CHKERRQ(ierr);
	PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}
