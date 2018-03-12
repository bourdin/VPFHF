/*
 test20.c: test for hydraulic fracturing in 3D with no leakoff
*/

#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFMech.h"
#include "VFFlow.h"
#include "VFPermfield.h"

VFCtx    ctx;
VFFields fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscViewer			    viewer;
	VFCtx               ctx;
	VFFields            fields;
	PetscErrorCode      ierr;
	PetscReal           radius = .1;
	PetscReal           center[3] = {0.5,0.5,0.5};
	PetscInt            orientation=1;
	PetscInt            nopts=3;
	/*PetscInt			      ek,ej,ei,c;*/
	PetscInt            i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
	PetscReal			****coords_array;
	PetscReal			***v_array;  
	PetscReal           BBmin[3],BBmax[3];
	PetscReal           x,y,z;  
	/*char                filename[FILENAME_MAX];*/
	PetscReal           p = 5.e-3;
	PetscReal           p_old = 0.;
 	/*PetscReal           p_epsilon = 1.e-6;*/
	PetscInt			altminit = 1;
	Vec					Vold;
	PetscReal			errV = 1e+10;
	PetscReal			q = 2.67e-4;
	//PetscReal           p_read,q_read=1e-4;
	PetscReal           vol_inj,volinit=0.;
	PetscReal           p_conv;
	PetscReal           CrackVolInit=0.;

  /*
    Do not do anything before PetscInitialize
  */	
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);

	ctx.maxtimestep = 10;
	p_conv = 1e-6;

	ierr = PetscOptionsGetReal(NULL,"-radius",&radius,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-p_conv",&p_conv,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-pinit",&p,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-volinit",&volinit,NULL);CHKERRQ(ierr);		
	ierr = PetscOptionsGetReal(NULL,"-rate",&q,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetRealArray(NULL,"-center",&center[0],&nopts,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,"-maxtimestep",&ctx.maxtimestep,NULL);CHKERRQ(ierr); 	
	ierr = PetscOptionsGetInt(NULL,"-orientation",&orientation,NULL);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(ctx.daScal,NULL,&nx,&ny,&nz,NULL,NULL,NULL,
					           NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
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
			ctx.bcU[1].face[X0]= ZERO;
			ctx.bcU[2].face[X0]= ZERO;			
			/*	face X1	*/
			ctx.bcU[0].face[X1]= ZERO;
			ctx.bcU[1].face[X0]= ZERO;
			ctx.bcU[2].face[X0]= ZERO;
			/*	face Y0	*/

			/*	face Y1	*/

			/*	face Z0	*/
			ctx.bcU[0].face[Z0]= ZERO;
			ctx.bcU[1].face[Z0]= ZERO;
			ctx.bcU[2].face[Z0]= ZERO;
			/*	face Z1	*/
			ctx.bcU[0].face[Z0]= ZERO;
			ctx.bcU[1].face[Z0]= ZERO;
			ctx.bcU[2].face[Z1]= ZERO;

			
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) { 
						x = coords_array[k][j][i][0];
						y = coords_array[k][j][i][1];
                        z = coords_array[k][j][i][2];						
						if ( ((j == ny/2) || (j == ny/2-1)) && ((x-center[0])*(x-center[0])+(z-center[2])*(z-center[2])) <= radius*radius ) {
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
			ctx.bcU[0].face[X0] = ZERO;
			//ctx.bcU[1].face[X0] = ZERO;
			/*	face X1	*/
			ctx.bcU[0].face[X1]= ZERO;
			ctx.bcU[1].face[X1]= ZERO;
			/*	face Y0	*/

			/*	face Y1	*/

			/*	face Z0	*/
			//ctx.bcU[1].face[Z0]= ZERO;
			ctx.bcU[2].face[Z0]= ZERO;
			/*	face Z1	*/
			ctx.bcU[1].face[Z1]= ZERO;
			ctx.bcU[2].face[Z1]= ZERO;

			
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
			SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: Orientation specified is not defined yet,got %i\n",orientation);
			break;
	}  	
	ierr = DMDAVecRestoreArray(ctx.daScal,fields.VIrrev,&v_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = VecCopy(fields.VIrrev,fields.V);CHKERRQ(ierr);
	ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);


	ierr = VecDuplicate(fields.V,&Vold);CHKERRQ(ierr);

	ierr = VecSet(fields.theta,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.thetaRef,0.0);CHKERRQ(ierr);
	ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
	for ( i=0; i < ctx.nlayer; i++) {
      ctx.matprop[i].alpha = 0.;
	  ctx.matprop[i].beta = 0.;
	}

	if(ctx.hasInsitu) {
	  ctx.hasCrackPressure = PETSC_FALSE;		
	  ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
	  ierr = VolumetricCrackOpening(&CrackVolInit,&ctx,&fields);CHKERRQ(ierr);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Crack Volume without pressure %g \n",CrackVolInit);CHKERRQ(ierr);
	}
	
	PetscViewerCreate(PETSC_COMM_WORLD,&viewer);
	PetscViewerSetType(viewer,PETSCVIEWERASCII);
	PetscViewerFileSetName(viewer,"pressure.txt");
	PetscViewerASCIIPrintf(viewer,"#step Volume        Pressure      Surf. energy  Elast. Energy Press. Work   Total Energy\n");

	
	ctx.hasCrackPressure = PETSC_TRUE;	
	ctx.timevalue = 0;
	vol_inj = volinit;
    ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
	for (ctx.timestep = 1; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
	  vol_inj += q;
//	  printf("Injected volume: %f\n",vol_inj);
//	  printf("Well pressure:   %f\n",p);
	  ierr = VecCopy(fields.VIrrev,fields.V);CHKERRQ(ierr);
	  do {
	    beginning:
        p_old = p;
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Time step %i,alt min step %i with pressure %g rate %g crack-vol %g\n",ctx.timestep,altminit,p,q,ctx.CrackVolume);CHKERRQ(ierr);
        ierr = VF_StepU(&fields,&ctx);CHKERRQ(ierr);
        ierr = VecScale(fields.U,1./p);CHKERRQ(ierr);		
        ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);
        if(ctx.hasInsitu){
          p = (vol_inj - CrackVolInit) / ctx.CrackVolume;
	      ierr = PetscPrintf(PETSC_COMM_WORLD,"Inj Vol %g Crack Vol wo pres %g Crack vol %g pres %g \n",vol_inj, CrackVolInit,ctx.CrackVolume,p);CHKERRQ(ierr);			  
        } else {	 
          p = vol_inj / ctx.CrackVolume;
        }
        ierr = VecCopy(fields.V,Vold);CHKERRQ(ierr);
        ierr = VecScale(fields.U,p);CHKERRQ(ierr);
        ierr = VecSet(fields.pressure,p);CHKERRQ(ierr);
        ierr = VF_StepV(&fields,&ctx);CHKERRQ(ierr);
        
		ierr = VecAXPY(Vold,-1.,fields.V);CHKERRQ(ierr);
        ierr = VecNorm(Vold,NORM_INFINITY,&errV);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"   Max. change on V: %e\n",errV);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"   Max. change on p: %e\n",PetscAbs(p-p_old));CHKERRQ(ierr);
        altminit++;
        if (altminit >= 200){
          vol_inj -= q;
          q = 0.5 * q;
          altminit = 0;
          vol_inj += q;
          goto beginning;
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Flow rate has been cut back... restarting the iteration \n");CHKERRQ(ierr);
        /*
          The model is rate independent, so changing the injection rate should have no effect on stability
        */
        }
	  } while (PetscAbs((p-p_old)/p) >= 1e-5 && altminit <= ctx.altminmaxit);

		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n Final Crack volume\t = %g,Pressure\t= %g\n\n",p*ctx.CrackVolume,p);CHKERRQ(ierr);	
		ierr = VolumetricCrackOpening(&ctx.CrackVolume,&ctx,&fields);CHKERRQ(ierr);   

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
		ierr = VF_VEnergy3D(&ctx.SurfaceEnergy,&fields,&ctx);CHKERRQ(ierr);
		ctx.TotalEnergy = ctx.ElasticEnergy - ctx.InsituWork - ctx.PressureWork;
		
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic Energy:            %e\n",ctx.ElasticEnergy);CHKERRQ(ierr);
		if (ctx.hasCrackPressure) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of pressure forces:   %e\n",ctx.PressureWork);CHKERRQ(ierr);
		}
		if (ctx.hasInsitu) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of surface forces:    %e\n",ctx.InsituWork);CHKERRQ(ierr);
		}
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Total energy:              %e\n",ctx.ElasticEnergy-ctx.InsituWork-ctx.PressureWork);CHKERRQ(ierr);
		
		PetscViewerASCIIPrintf(viewer,"%d \t\t%e \t%e \t%e \t%e \t%e \t%e\n",ctx.timestep ,ctx.CrackVolume,p,ctx.SurfaceEnergy,ctx.ElasticEnergy,ctx.PressureWork,ctx.TotalEnergy);

		altminit = 0.;
	}
	PetscViewerFlush(viewer);CHKERRQ(ierr);
	PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}
