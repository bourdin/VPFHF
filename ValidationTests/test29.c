/*
 test28.c: 1D. Transient flow problem with no source term and homogeneous pressure boundary conditions on all sides taken from http://tutorial.math.lamar.edu/Classes/DE/SolvingHeatEquation.aspx. Initial pressure is p (x,0)= 6sin(pi*x/Lx)
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
./test29 -n 51,51,2 -l 1,1,0.01 -m_inv 1 -ts_dt 0.005 -ts_max_steps 20 -flowsolver FLOWSOLVER_tsMIXEDFEM
./test29 -n 51,51,2 -l 1,1,0.01 -m_inv 1 -maxtimestep 20 -flowsolver FLOWSOLVER_snesMIXEDFEM -timestepsize 0.005
 
./test29 -n 51,51,2 -l 1,1,0.01 -maxtimestep 20 -flowsolver FLOWSOLVER_snesstandardFEM -m_inv 1 -maxtimestep 20 -flowsolver FLOWSOLVER_snesstandarDFEM -P_X0_BC FIXED -P_X1_BC FIXED -Q_Y0_BC_1 FIXED -Q_Y1_BC_1 FIXED -Q_Z0_BC_2 FIXED -Q_Z1_BC_2 FIXED -theta 1 -deltat 0.005 -maxtimestep 20 -theta 1 -k 1,1,1 -g 0,0,0 -rhof 0. -mu 1 

 */

#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFMech.h"
#include "VFFlow.h"
#include "VFPermfield.h"

VFCtx               ctx;
VFFields            fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{	
	PetscErrorCode  ierr;
	PetscInt		i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
	PetscReal		BBmin[3],BBmax[3];
	PetscReal		****coords_array;
	PetscReal		lx;
  PetscReal		****velnpre_array;
	PetscReal		***pres_ini_array;
  PetscReal   vol,vol1,vol2,vol3,vol4,vol5;

	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ctx.flowsolver = FLOWSOLVER_SNESSTANDARDFEM;
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
  ctx.FlowDisplCoupling = PETSC_FALSE;
	ctx.hasFluidSources = PETSC_FALSE;
	ctx.hasFlowWells = PETSC_FALSE;
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);	
	lx = BBmax[0]-BBmin[0];
	ierr = DMDAVecGetArrayDOF(ctx.daFlow,fields.VelnPress,&velnpre_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx.daScal,fields.pressure,&pres_ini_array);CHKERRQ(ierr);
  for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				pres_ini_array[k][j][i] = 6.*sin(PETSC_PI*coords_array[k][j][i][0]/lx);
			}
		}
	}
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				velnpre_array[k][j][i][3] = pres_ini_array[k][j][i];
			}
		}
	}
  
	ierr = DMDAVecRestoreArray(ctx.daScal,fields.pressure,&pres_ini_array);CHKERRQ(ierr);  
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);    
  if(ctx.flowsolver == FLOWSOLVER_TSMIXEDFEM || ctx.flowsolver == FLOWSOLVER_TS){
  	ctx.maxtimestep = 20;
    ctx.maxtimevalue = 100.;
    ctx.timevalue = 0.1;
    ierr = VF_StepP(&fields,&ctx);
    ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
}
  else{
  ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
  ierr = VecCopy(fields.pressure,ctx.pressure_old);CHKERRQ(ierr);
	for (ctx.timestep = 0; ctx.timestep < ctx.maxtimestep; ctx.timestep++){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nProcessing step %i.\n",ctx.timestep);CHKERRQ(ierr);
    ierr = VF_StepP(&fields,&ctx);
    ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
    ierr = VFCheckVolumeBalance(&vol,&vol1,&vol2,&vol3,&vol4,&vol5,&ctx,&fields);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n modulus_volume = %g\n",vol);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," divergence_volume = %g\n",vol1);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," surface_flux_volume = %g\n",vol2);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," well_volume = %g\n",vol3);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," source_volume = %g\n",vol4);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," vol.strain_volume = %g\n",vol5);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Volume Balance ::: RHS = %g \t LHS = %g \n",vol+vol1,vol3+vol4+vol5);CHKERRQ(ierr);

    /*This will have to be called "an update function"*/
    ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSVelP,ctx.RHSVelPpre);CHKERRQ(ierr);
    ierr = VecCopy(fields.pressure,ctx.pressure_old);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSP,ctx.RHSPpre);CHKERRQ(ierr);
	}
	Vec error;
	PetscReal norm_1,norm_2,norm_inf;
	ierr = VecDuplicate(fields.VelnPress,&error);
	ierr = VecWAXPY(error,-1.0,fields.VelnPress,fields.FlowBCArray);
	ierr = VecNorm(error,NORM_1,&norm_1);
	ierr = VecNorm(error,NORM_2,&norm_2);
	ierr = VecNorm(error,NORM_INFINITY,&norm_inf);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n1_NORM = %f \n 2_norm = %f \n inf_norm = %f \n",norm_1, norm_2,norm_inf);CHKERRQ(ierr);
  }
	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

