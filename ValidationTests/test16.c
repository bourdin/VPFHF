/*
 test16.c: Solves for the displacement and v-field in a volume loaded line crack in 2d (Sneddon 2D)
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu

Try:  
mpiexec -n 2 ./test16  -n 26,51,2 -l 0.5,1,0.02 -epsilon 0.02 -length 2. \
      -center 0,0.5,0 -orientation 2 -nu 0. -eta 1e-8 -maxtimestep 10 -maxvol .03 \
      -Gc 1.e-1  -insitumin 0.,-1.e-4,0.,0.,0.,0. -insitumax 0.,-1.e-4,0.,0.,0.,0.

 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFV.h"
#include "VFU.h"
#include "VFPermfield.h"

VFCtx               ctx;
VFFields            fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode  ierr;
	PetscViewer			viewer;
	PetscViewer     logviewer;
	PetscReal       length = .2;
	PetscReal       center[3]={0.,0.,.5};
	PetscInt        orientation=1;
	PetscInt        nopts=3;
	PetscInt        i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
	PetscReal		****coords_array;
	PetscReal		 ***v_array;  
	PetscReal       BBmin[3],BBmax[3];
	char            filename[FILENAME_MAX];
	PetscReal       p;
	PetscReal       p_old;
	PetscReal       p_epsilon = 1.e-5;
	PetscInt			  altminit=1;
	Vec					    Vold,U_s,U_1;
	PetscReal       vol_s,vol_1;
	PetscReal			  errV=1e+10;
	PetscReal			  lx,ly,lz;
	PetscReal			  q,maxvol = .03;
	
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
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-maxvol",&maxvol,PETSC_NULL);CHKERRQ(ierr);
  /*
    Overwrite ctx.maxtimestep with something more reasonable
  */
  ctx.maxtimestep = 150;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-maxtimestep",&ctx.maxtimestep,PETSC_NULL);CHKERRQ(ierr);
	q = maxvol / ctx.maxtimestep;
	
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VecSet(fields.VIrrev,1.0);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx.daScal,fields.VIrrev,&v_array);CHKERRQ(ierr);    
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
	switch (orientation) {
		case 1:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a line crack of length %g at (%g,%g,%g) with normal vector <0,1,0>\n",
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
						if ( ((j == ny/2) || (j == ny/2-1)) && (coords_array[k][j][i][0] > lx/2.-length) && (coords_array[k][j][i][0] < lx/2.+length ) ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		case 2:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a line crack of length %g at (%g,%g,%g) with normal vector <0,1,0>\n",
							   length,center[0],center[1],center[2]);CHKERRQ(ierr);		  
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
			//ctx.bcU[2].face[X0]= ZERO;
			/*	face X1	*/
			ctx.bcU[0].face[X1]= ZERO;
			//ctx.bcU[2].face[X1]= ZERO;
			/*	face Y0	*/
			//ctx.bcU[1].face[Y0]= ZERO;
			//ctx.bcU[2].face[Y0]= ZERO;
			/*	face Y1	*/
			//ctx.bcU[1].face[Y1]= ZERO;
			//ctx.bcU[2].face[Y1]= ZERO;		  
			/*	face Z0	*/
			ctx.bcU[2].face[Z0]= ZERO;
			/*	face Z1	*/
			ctx.bcU[2].face[Z1]= ZERO;
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) { 
						if ( ((j == ny/2) || (j == ny/2-1)) && (coords_array[k][j][i][0] > lx/2.-length) && (coords_array[k][j][i][0] < lx/2.+length ) ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		case 3:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a line crack of length %g at (%g,%g,%g) with normal vector <0,1,0>\n",
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
						if ( ((j == ny/2) || (j == ny/2-1)) && (coords_array[k][j][i][0] > lx/2.-length) && (coords_array[k][j][i][0] < lx/2.+length ) ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		case 4:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a line crack of length %g at (%g,%g,%g) with normal vector <0,1,0>\n",
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
						if ( ((j == ny/2) || (j == ny/2-1)) && (coords_array[k][j][i][0] > lx/2.-length) && (coords_array[k][j][i][0] < lx/2.+length ) ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		case 5:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a line crack of length %g at (%g,%g,%g) with normal vector <0,1,0>\n",
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
						if ( ((j == ny/2) || (j == ny/2-1)) && (coords_array[k][j][i][0] > lx/2.-length) && (coords_array[k][j][i][0] < lx/2.+length ) ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		case 6:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a line crack of length %g at (%g,%g,%g) with normal vector <0,1,0>\n",
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
						if ( ((j == ny/2) || (j == ny/2-1)) && (coords_array[k][j][i][0] > lx/2.-length) && (coords_array[k][j][i][0] < lx/2.+length ) ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		case 7:
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Building a line crack of length %g at (%g,%g,%g) with normal vector <0,1,0>\n",
							   length,center[0],center[1],center[2]);CHKERRQ(ierr);		  
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting up BC for in situ stress in the <0,0,1> direction only\n");CHKERRQ(ierr);
			/*	face X0	*/
			ctx.bcU[0].face[X0]= ZERO;
			ctx.bcV[0].face[X0]= ONE;
			//ctx.bcU[1].face[X0]= ZERO;
			//ctx.bcU[2].face[X0]= ZERO;
			/*	face X1	*/
			ctx.bcU[0].face[X1]= ZERO;
			ctx.bcV[0].face[X1]= ONE;
			//ctx.bcU[1].face[X1]= ZERO;
			//ctx.bcU[2].face[X1]= ZERO;
			/*	face Y0	*/
			ctx.bcU[1].face[Y0]= ZERO;
			//ctx.bcV[0].face[Y0]= FIXED;
			/*	face Y1	*/
			ctx.bcU[1].face[Y1]= ZERO;
			//ctx.bcV[0].face[Y1]= FIXED;
			/*	face Z0	*/
			ctx.bcU[2].face[Z0]= ZERO;
			ctx.bcV[0].face[Z0]= ONE;
			/*	face Z1	*/
			//ctx.bcU[2].face[Z1]= ZERO;
			ctx.bcV[0].face[Z1]= ONE;
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) { 
						if ( ((k == nz/2) || (k == nz/2-1)) && (coords_array[k][j][i][0] > lx/2.-length) && (coords_array[k][j][i][0] < lx/2.+length ) ) {
							v_array[k][j][i] = 0.;
						}
					}
				}
			}      
			break;
		default:
			SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Orientation should be one of {1,2,3,4,5,6,7},got %i\n",orientation);
			break;
	}  
	ierr = DMDAVecRestoreArray(ctx.daScal,fields.VIrrev,&v_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VecCopy(fields.VIrrev,fields.V);CHKERRQ(ierr);
	ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
	
	ctx.hasCrackPressure = PETSC_TRUE;
	ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);
	ierr = VecDuplicate(fields.U,&U_s);CHKERRQ(ierr);
	ierr = VecDuplicate(fields.U,&U_1);CHKERRQ(ierr);

	ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
	ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
	ierr = PetscViewerFileSetMode(viewer,FILE_MODE_APPEND);CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.pres",ctx.prefix);CHKERRQ(ierr);
	ierr = PetscViewerFileSetName(viewer,filename);CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(viewer,"#Time step \t Volume \t Pressure \t SurfaceEnergy \t ElasticEnergy \t PressureForces \t TotalMechEnergy \n");CHKERRQ(ierr);
	
	p = 1.e-5;
	ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
	ierr = VecSet(fields.pressureRef,0.0);CHKERRQ(ierr);
	
	ierr = VecSet(U_s,0.0);CHKERRQ(ierr);
	ierr = VecSet(U_1,0.0);CHKERRQ(ierr);
	ctx.matprop[0].beta = 0.;
	ctx.timevalue = 0;
	//ctx.maxtimestep = 150;
	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"ctx.hasInsitu: %d\n",ctx.hasInsitu);CHKERRQ(ierr);
	
	for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Time step %i, injected volume %g\n",ctx.timestep,q*ctx.timestep);CHKERRQ(ierr);
		ierr = VecCopy(fields.V,fields.VIrrev);CHKERRQ(ierr); 
		do {
			p_old = p;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"  Time step %i,alt min step %i with pressure %g\n",ctx.timestep,altminit,p);CHKERRQ(ierr);
			
			/* 
			  Update the pressure based on the relation
			      V = p vol_1 + vol_s 
			  with
			      vol_1 = \int U_1 \cdot \nabla V
			  where U_1 is the displacement field associated with null in-situ stress and unit pressure, and
			      vol_s = \int U_s \cdot \nabla V
			  where U_s is the displacement field associated with null pressure and in-situ stress
			*/
    	ctx.hasCrackPressure = PETSC_TRUE;
      ctx.hasInsitu = PETSC_FALSE;
      ierr = VecSet(fields.pressure,1.0);CHKERRQ(ierr);
      ierr = VecCopy(U_1,fields.U);CHKERRQ(ierr);
			ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
			ierr = VolumetricCrackOpening(&vol_1,&ctx,&fields);CHKERRQ(ierr);         
      ierr = VecCopy(fields.U,U_1);CHKERRQ(ierr);

    	ctx.hasCrackPressure = PETSC_FALSE;
      ctx.hasInsitu = PETSC_TRUE;
      ierr = VecCopy(U_s,fields.U);CHKERRQ(ierr);
			ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
			ierr = VolumetricCrackOpening(&vol_s,&ctx,&fields);CHKERRQ(ierr);   
      ierr = VecCopy(fields.U,U_s);CHKERRQ(ierr);
      
      // This will fail if vol_1 = 0, which should only happen when there are no cracks
      p = (q*ctx.timestep - vol_s) / vol_1;
      ierr = VecAXPY(fields.U,p,U_1);CHKERRQ(ierr);

      ctx.CrackVolume = vol_s + p * vol_1;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"      vol_s: %e vol_1 %e volume\n",vol_s,vol_1,ctx.CrackVolume);
			ierr = PetscPrintf(PETSC_COMM_WORLD,"      Updated crack pressure: %e (was %e)\n",p,p_old);

			ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
			ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
			ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);

			ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
			ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,"      Max. change on V: %e\n",errV);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,"      Max. change on p: %e\n",PetscAbs(p-p_old));CHKERRQ(ierr);
			altminit++;
		} while (PetscAbs(p-p_old) >= p_epsilon && altminit <= ctx.altminmaxit);
		ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);   
		switch (ctx.fileformat) {
			case FILEFORMAT_HDF5:       
				ierr = FieldsH5Write(&ctx,&fields);CHKERRQ(ierr);
				break;
			case FILEFORMAT_BIN:
				ierr = FieldsBinaryWrite(&ctx,&fields);CHKERRQ(ierr);
				break; 
		}
		ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.log",ctx.prefix);CHKERRQ(ierr);
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&logviewer);CHKERRQ(ierr);
		ierr = PetscLogView(logviewer);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&logviewer);CHKERRQ(ierr);

		ctx.ElasticEnergy=0;
		ctx.InsituWork=0;
		ctx.PressureWork = 0.;
		ierr = VF_UEnergy3D(&ctx.ElasticEnergy,&ctx.InsituWork,&ctx.PressureWork,&fields,&ctx);CHKERRQ(ierr);
		ierr = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);
		ctx.TotalEnergy = ctx.ElasticEnergy - ctx.InsituWork - ctx.PressureWork + ctx.SurfaceEnergy;
		ierr = PetscViewerASCIIPrintf(viewer,"%d \t %e \t %e \t %e \t %e \t %e \t %e\n",ctx.timestep ,ctx.CrackVolume,p,ctx.SurfaceEnergy,ctx.ElasticEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%i   \t%e   \t%e   \t%e   \t%e   \t%e\n",ctx.timestep,ctx.ElasticEnergy,
                                  ctx.InsituWork,ctx.SurfaceEnergy,ctx.PressureWork,ctx.TotalEnergy);CHKERRQ(ierr);
		altminit = 0.;
	}
	ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&logviewer);CHKERRQ(ierr);
	ierr = VecDestroy(&Vold);CHKERRQ(ierr);
	ierr = VecDestroy(&U_s);CHKERRQ(ierr);
	ierr = VecDestroy(&U_1);CHKERRQ(ierr);

	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

