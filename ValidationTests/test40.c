/*
 test40.c: 2D SNES. Single well problem. Source term implemented as dirac function and analytical solution is 2D Green's function. All pressure boundary condition.
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 
 ./test40  -n 201,201,2 -l 1,1,1 -maxtimestep 1 timestepsize 1 -theta 1 -nw 2 -w0_coords 0.50001,0.50001,0.00 -w0_Qw 0.5 -w0_constraint Rate -w0_rw 0.1 -w0_type INJECTOR -w1_coords 0.50001,0.50001,1 -w1_Qw 0.5 -w1_constraint Rate -w1_rw 0.1 -w1_type INJECTOR -m_inv 0
 
 ./test40  -n 101,101,2 -l 1,1,1 -maxtimestep 1 timestepsize 1 -theta 1 -nw 2 -w0_coords 0.50001,0.50001,0.00 -w0_Qw 0.5 -w0_constraint Rate -w0_rw 0.1 -w0_type INJECTOR -w1_coords 0.50001,0.50001,1 -w1_Qw 0.5 -w1_constraint Rate -w1_rw 0.1 -w1_type INJECTOR -m_inv 0
 
 
 ./test40  -n 101,101,2 -l 1,1,1 -maxtimestep 1 timestepsize 1 -theta 1 -nw 1 -w0_coords 0.50001,0.50001,0.5 -w0_Qw 1 -w0_constraint Rate -w0_rw 0.1 -w0_type producer -m_inv 0
 
 ./test40  -n 101,101,2 -l 1,1,1 -maxtimestep 1 timestepsize 10 -theta 1 -nw 1 -w0_coords 0.5,0.5,0.5 -w0_Qw 1 -w0_constraint Rate -w0_rw 0.1 -w0_type producer -m_inv 0 -flowsolver FLOWSOLVER_SNESMIXEDFEM -flowsnes_snes_type tr -flowsnes_snes_view
 
 ./test40  -n 41,41,2 -l 1,1,0.01 -maxtimestep 1 timestepsize 10 -theta 1 -nw 1 -w0_coords 0.5,0.5,0.005 -w0_Qw 1 -w0_constraint Rate -w0_rw 0.1 -w0_type injector -m_inv 0 -flowsolver FLOWSOLVER_kspMIXEDFEM -velp_ksp_view
 
 
 
 ./test40  -n 101,101,2 -l 1,1,1 -maxtimestep 1 timestepsize 1 -theta 1 -nw 1 -w0_coords 0.50001,0.50001,0.5 -w0_Qw 1 -w0_constraint Rate -w0_rw 0.1 -w0_type pRoducer -m_inv 0. -flowsolver FLOWSOLVER_snesstandarDFEM
 
 
 
 
 
 
 
 
 
 
 ./test40  -n 101,101,2 -l 1,1,1 -maxtimestep 1 timestepsize 1 -theta 1 -nw 1 -w0_coords 0.50001,0.50001,0.5 -w0_Qw 1 -w0_constraint Rate -w0_rw 0.1 -w0_type pRoducer -m_inv 0 -flowsolver FLOWSOLVER_tsMIXEDFEM -ts_dt 1 -ts_max_steps 1
 
./test40  -n 101,101,2 -l 1,1,1 -maxtimestep 1 timevalue 1 -theta 1 -nw 2 -w0_coords 0.5,0.5,0.5 -w0_Qw 1 -w0_constraint Rate -w0_rw 0.1 -w0_type producer -flowsolver FLOWSOLVER_snesStandarDFEM -ks 1e+30 -kf 1e+30 -flowwells -P_X0_BC FIXED -P_X1_BC FIXED -P_Y0_BC FIXED -P_Y1_BC FIXED -Q_Z0_BC_2 FIXED -Q_Z1_BC_2 FIXED -w1_coords 0.5,0.75,0.5 -w1_Qw 1 -w1_constraint Rate -w1_rw 0.1 -w1_type producer 
 
 ./test40  -n 101,101,2 -l 1,1,1 -maxtimestep 1 timevalue 1 -theta 1 -nw 1 -w0_coords 0.5,0.5,0.5 -w0_Qw 1 -w0_constraint Rate -w0_rw 0.1 -w0_type producer -flowsolver FLOWSOLVER_snesStandarDFEM -ks 1e+30 -kf 1e+30 -flowwells -P_X0_BC FIXED -P_X1_BC FIXED -P_Y0_BC FIXED -P_Y1_BC FIXED -Q_Z0_BC_2 FIXED -Q_Z1_BC_2 FIXED
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
  PetscInt     i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal    BBmin[3],BBmax[3];
  PetscReal    ***presbc_array;
  PetscReal    ****coords_array;
  PetscReal    pi,dist;
  PetscReal    ***pre_array;
  PetscInt     xs1,xm1,ys1,ym1,zs1,zm1;
  PetscReal    vol,vol1,vol2,vol3,vol4,vol5;

  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = VFInitialize(&ctx,&fields);CHKERRQ(ierr);
  ierr = DMDAGetInfo(ctx.daScal,NULL,&nx,&ny,&nz,NULL,NULL,NULL,
                     NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx.daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DMDAGetBoundingBox(ctx.daVect,BBmin,BBmax);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx.daVFperm,&xs1,&ys1,&zs1,&xm1,&ym1,&zm1);CHKERRQ(ierr);
  pi = 6.*asin(0.5);
  for (k = zs; k < zs+zm; k++) {
    for (j = ys; j < ys+ym; j++) {
      for (i = xs; i < xs+xm; i++) {
        if(i == 0){
          dist = sqrt(pow((ctx.well[0].coords[0]-coords_array[k][j][i][0]),2)+pow((ctx.well[0].coords[1]-coords_array[k][j][i][1]),2));
          presbc_array[k][j][i] = 1./(2.*pi)*log(dist);
            //            presbc_array[k][j][i] = 10.;
        }
        if(i == nx-1){
          dist = sqrt(pow((ctx.well[0].coords[0]-coords_array[k][j][i][0]),2)+pow((ctx.well[0].coords[1]-coords_array[k][j][i][1]),2));
          presbc_array[k][j][i] = 1./(2.*pi)*log(dist);
            //            presbc_array[k][j][i] = 10.;
        }
        if(j == 0){
          dist = sqrt(pow((ctx.well[0].coords[0]-coords_array[k][j][i][0]),2)+pow((ctx.well[0].coords[1]-coords_array[k][j][i][1]),2));
          presbc_array[k][j][i] = 1./(2.*pi)*log(dist);
            //            presbc_array[k][j][i] = 10.;
        }
        if(j == ny-1){
          dist = sqrt(pow((ctx.well[0].coords[0]-coords_array[k][j][i][0]),2)+pow((ctx.well[0].coords[1]-coords_array[k][j][i][1]),2));
          presbc_array[k][j][i] = 1./(2.*pi)*log(dist);
            //            presbc_array[k][j][i] = 10.;
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(ctx.daScal,ctx.PresBCArray,&presbc_array);CHKERRQ(ierr);
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
  
    
    ierr = VecCopy(fields.VelnPress,ctx.PreFlowFields);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSVelP,ctx.RHSVelPpre);CHKERRQ(ierr);
    ierr = VecCopy(fields.pressure,ctx.pressure_old);CHKERRQ(ierr);
    ierr = VecCopy(ctx.RHSP,ctx.RHSPpre);CHKERRQ(ierr);
  }
  ierr = VecSet(fields.VelnPress,0.);CHKERRQ(ierr);
  ierr = VecSet(fields.velocity,0.);CHKERRQ(ierr);
  ierr = VecSet(fields.pressure,0.);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx.daScal,fields.pressure,&pre_array);CHKERRQ(ierr);
  for (k = zs; k < zs+zm; k++) {
    for (j = ys; j < ys+ym; j++) {
      for (i = xs; i < xs+xm; i++) {
        dist = sqrt(pow((ctx.well[0].coords[0]-coords_array[k][j][i][0]),2)+pow((ctx.well[0].coords[1]-coords_array[k][j][i][1]),2));
        pre_array[k][j][i] = 1./(2.*pi)*log(dist);
        if(i == 0 && j == 0 && k == 0){
//          ierr = PetscPrintf(PETSC_COMM_WORLD,"dist = %g; pressure = %g \n",dist,pre_array[k][j][i]);CHKERRQ(ierr);
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(ctx.daScal,fields.pressure,&pre_array);CHKERRQ(ierr);
  ++ctx.timestep;
  ierr = FieldsVTKWrite(&ctx,&fields,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx.daVect,ctx.coordinates,&coords_array);CHKERRQ(ierr);

  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}

