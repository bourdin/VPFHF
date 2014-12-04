/*
 test24.c: 3D KSP. Flow problem with source term [pressure = sin(2*pi*x)*sin(2*pi*y)*sin(2(pi*z)]. All normal velocity boundary condition.
 (c) 2012-2014 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
 ./test24 -n 11,11,11 -l 1,1,1 -flowsolver FLOWSOLVER_SNESMIXEDFEM
 ./test24 -n 11,11,11 -l 1,1,1 -m_inv 0 -ts_type beuler -ts_dt 1 -ts_max_steps 2 -flowsolver FLOWSOLVER_tSMIXEDFEM
 ./test24 -n 11,11,11 -l 1,1,1 -m_inv 0 -flowsolver FLOWSOLVER_SNESstandarDFEM -Q_X0_BC_0 FIXED -Q_X1_BC_0 FIXED -Q_Y0_BC_1 FIXED -Q_Y1_BC_1 FIXED -Q_Z0_BC_2 FIXED -Q_Z1_BC_2 FIXED -k 1,1,1 -g 0,0,0 -rhof 0. -mu 1
 
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
	PetscReal		****velbc_array;
	PetscReal		***src_array;
	PetscReal		****coords_array;
	PetscReal		lx,ly,lz;
  PetscReal   vol,vol1,vol2,vol3,vol4,vol5;
  
	ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
	ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
	ierr = DMDAGetInfo(ctx.daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
  ctx.FlowDisplCoupling = PETSC_FALSE;
	ctx.hasFluidSources = PETSC_TRUE;
	ctx.hasFlowWells = PETSC_FALSE;
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx.daScal,ctx.Source,&src_array);CHKERRQ(ierr);
	lz = BBmax[2]-BBmin[2];
	ly = BBmax[1]-BBmin[1];
	lx = BBmax[0]-BBmin[0];
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				velbc_array[k][j][i][0] = -1./ctx.flowprop.mu*(2.*PETSC_PI/lx*cos(2.*PETSC_PI*coords_array[k][j][i][0]/lx)*sin(2.*PETSC_PI*coords_array[k][j][i][1]/ly)*sin(2.*PETSC_PI*coords_array[k][j][i][2]/lz)-ctx.flowprop.rho*ctx.flowprop.g[0]);
				velbc_array[k][j][i][1] = -1./ctx.flowprop.mu*(2.*PETSC_PI/ly*sin(2.*PETSC_PI*coords_array[k][j][i][0]/lx)*cos(2.*PETSC_PI*coords_array[k][j][i][1]/ly)*sin(2.*PETSC_PI*coords_array[k][j][i][2]/lz)-ctx.flowprop.rho*ctx.flowprop.g[1]);
				velbc_array[k][j][i][2] = -1./ctx.flowprop.mu*(2.*PETSC_PI/lz*sin(2.*PETSC_PI*coords_array[k][j][i][0]/lx)*sin(2.*PETSC_PI*coords_array[k][j][i][1]/ly)*cos(2.*PETSC_PI*coords_array[k][j][i][2]/lz)-ctx.flowprop.rho*ctx.flowprop.g[2]);
			}
		}
	}
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
         src_array[k][j][i] = 4.*pow(PETSC_PI,2)*1./ctx.flowprop.mu*sin(2.*PETSC_PI*coords_array[k][j][i][2]/lz)*sin(2.*PETSC_PI*coords_array[k][j][i][1]/ly)*sin(2.*PETSC_PI*coords_array[k][j][i][0]/lx)*( (1/(lx*lx))+(1/(ly*ly))+(1/(lz*lz)) );
			}
		}
	}  
	ierr = DMDAVecRestoreArray(ctx.daScal,ctx.Source,&src_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.VelBCArray,&velbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);  
  if(ctx.flowsolver == FLOWSOLVER_TSMIXEDFEM || ctx.flowsolver == FLOWSOLVER_TS){
  	ctx.maxtimestep = 1;
    ctx.maxtimevalue = 60.;
    ctx.timevalue = 1.;
    ierr = VF_StepP(&fields,&ctx);
    ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
  }
  else{
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
      ierr = PetscPrintf(PETSC_COMM_WORLD," Volume Balance: RHS = %g \t LHS = %g \n",vol+vol1,vol3+vol4+vol5);CHKERRQ(ierr);
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
  
  ierr = VF_FastFourierTransforms(&ctx,&fields);
  ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);

	ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
	ierr = PetscFinalize();
	return(0);
}

