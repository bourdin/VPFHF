/*
 VFFlow_FractureFlow.c
 A mixed finite elements Darcy solver for fracture flow based on the method in
 Masud, A. and Hughes, T. J. (2002). A stabilized mixed finite element method for
 Darcy flow. Computer Methods in Applied Mechanics and Engineering, 191(3940):43414370.
 
 (c) 2011-2012 C. Chukwudozie, LSU
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFFlow.h"
#include "VFFlow_SNESMixedFEM.h"
#include "VFFlow_KSPMixedFEM.h"
#include "VFWell.h"
#include "VFMech.h"
#include "VFFlow_FractureFlow.h"

/*
 ################################################################################################################
 SNES ROUTINE
 ################################################################################################################
 */
#undef __FUNCT__
#define __FUNCT__ "MixedFractureFlowSolverInitialize"
extern PetscErrorCode MixedFractureFlowSolverInitialize(VFCtx *ctx, VFFields *fields)
{
  PetscMPIInt    comm_size;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  /*
  Moving this into VFCtxGet since it is an attribute of VFCtx
  If we split VFCtx into VFFlowCtx, VFMechCtx , and VFHeatCtx, it will have to go into VFFlow.c
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"","");CHKERRQ(ierr);
  {
    ctx->units    = UnitaryUnits;
    ierr          = PetscOptionsEnum("-flowunits","\n\tFlow solver","",VFUnitName,(PetscEnum)ctx->units,(PetscEnum*)&ctx->units,PETSC_NULL);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  */
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&comm_size);CHKERRQ(ierr);
  ierr = DMCreateMatrix(ctx->daFlow,MATAIJ,&ctx->KFracVelP);CHKERRQ(ierr);
  ierr = DMCreateMatrix(ctx->daFlow,MATAIJ,&ctx->KFracVelPlhs);CHKERRQ(ierr);
  ierr = DMCreateMatrix(ctx->daFlow,MATAIJ,&ctx->JacFracVelP);CHKERRQ(ierr);
  ierr = MatZeroEntries(ctx->JacFracVelP);CHKERRQ(ierr);
  ierr = MatZeroEntries(ctx->KFracVelPlhs);CHKERRQ(ierr);
  ierr = MatZeroEntries(ctx->KFracVelP);CHKERRQ(ierr);
  ierr = MatSetOption(ctx->KFracVelP,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatSetOption(ctx->KFracVelPlhs,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatSetOption(ctx->JacFracVelP,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(ctx->daFlow,&ctx->RHSFracVelP);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)ctx->RHSFracVelP,"RHS vector of fracture flow equation");CHKERRQ(ierr);
  ierr = VecSet(ctx->RHSFracVelP,0.0);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(ctx->daFlow,&ctx->RHSFracVelPpre);CHKERRQ(ierr);
  ierr = VecSet(ctx->RHSFracVelPpre,0.);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(ctx->daFlow,&ctx->PreFracFlowFields);CHKERRQ(ierr);
  ierr = VecSet(ctx->PreFracFlowFields,0.0);CHKERRQ(ierr);
//  
//  ierr = PetscPrintf(PETSC_COMM_WORLD," Step:I am initializing \n\n");CHKERRQ(ierr);

  
  ierr = DMCreateGlobalVector(ctx->daFlow,&ctx->FracResidual);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)ctx->FracResidual,"RHS of fracture SNES flow solver");CHKERRQ(ierr);
  ierr = VecSet(ctx->FracResidual,0.0);CHKERRQ(ierr);
  
  ierr = SNESCreate(PETSC_COMM_WORLD,&ctx->snesFracVelP);CHKERRQ(ierr);
  ierr = SNESAppendOptionsPrefix(ctx->snesFracVelP,"FracSnes_");CHKERRQ(ierr);
  ierr = SNESSetFromOptions(ctx->snesFracVelP);CHKERRQ(ierr);
  ierr = BCFracQInit(&ctx->bcFracQ[0],ctx);
  ierr = GetFlowProp(&ctx->flowprop,ctx->units,ctx->resprop);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MixedFractureFlowSolverFinalize"
extern PetscErrorCode MixedFractureFlowSolverFinalize(VFCtx *ctx,VFFields *fields)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = MatDestroy(&ctx->KFracVelP);CHKERRQ(ierr);
  ierr = MatDestroy(&ctx->KFracVelPlhs);CHKERRQ(ierr);
  ierr = MatDestroy(&ctx->JacFracVelP);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->RHSFracVelP);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->RHSFracVelPpre);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->FracResidual);CHKERRQ(ierr);
  ierr = SNESDestroy(&ctx->snesFracVelP);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->PreFracFlowFields);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MixedFracFlowSNESSolve"
extern PetscErrorCode MixedFracFlowSNESSolve(VFCtx *ctx,VFFields *fields)
{
  PetscErrorCode     ierr;
  SNESConvergedReason reason;
  PetscReal       ****VelnPress_array;
  PetscReal        ***Press_array;
  PetscReal       ****vel_array;
  PetscInt            i,j,k,c,veldof = 3;
  PetscInt            xs,xm,ys,ym,zs,zm;
  PetscInt            its;
  PetscReal           VelPmin,VelPmax;
  
  
  PetscFunctionBegin;
  
//  ierr = VecCopy(fields->fracVelnPress,ctx->RHSFracVelP);CHKERRQ(ierr);
  ierr = VecCopy(fields->fracVelnPress,ctx->PreFracFlowFields);CHKERRQ(ierr);
  ierr = FormFracMatricesnVector(ctx->KFracVelP,ctx->KFracVelPlhs,ctx->RHSFracVelP,ctx,fields);CHKERRQ(ierr);
  
  /*
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"RHSvec.txt",&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
  ierr = VecView(ctx->RHSFracVelP,viewer);CHKERRQ(ierr);
  
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"Matrix.txt",&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
  ierr = MatView(ctx->KFracVelP,viewer);CHKERRQ(ierr);
  */
  ierr = SNESSetFunction(ctx->snesFracVelP,ctx->FracResidual,FormFracSNESIFunction,ctx);CHKERRQ(ierr);
  ierr = SNESSetJacobian(ctx->snesFracVelP,ctx->JacFracVelP,ctx->JacFracVelP,FormFracSNESIJacobian,ctx);CHKERRQ(ierr);
  if (ctx->verbose > 1) {
    ierr = SNESMonitorSet(ctx->snesFracVelP,FEMSNESMonitor,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  }
  ierr = SNESSolve(ctx->snesFracVelP,PETSC_NULL,fields->fracVelnPress);CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(ctx->snesFracVelP,&reason);CHKERRQ(ierr);
  if (reason < 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"[ERROR] Fracture snesVelP diverged with reason %d\n",(int)reason);CHKERRQ(ierr);
  } else {
    ierr = SNESGetIterationNumber(ctx->snesFracVelP,&its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Fracture snesVelP converged in %d iterations %d.\n",(int)its,(int)reason);CHKERRQ(ierr);
  }
  ierr = VecMin(fields->fracVelnPress,PETSC_NULL,&VelPmin);CHKERRQ(ierr);
  ierr = VecMax(fields->fracVelnPress,PETSC_NULL,&VelPmax);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"     fracture VelP min / max:     %e %e\n",VelPmin,VelPmax);CHKERRQ(ierr);
  ierr = VecCopy(ctx->RHSFracVelP,ctx->RHSFracVelPpre);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,fields->fracpressure,&Press_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daVect,fields->fracvelocity,&vel_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daFlow,fields->fracVelnPress,&VelnPress_array);CHKERRQ(ierr);
  for (k = zs; k < zs+zm; k++) {
    for (j = ys; j < ys+ym; j++) {
      for (i = xs; i < xs+xm; i++) {
        Press_array[k][j][i] = VelnPress_array[k][j][i][3];
        for (c = 0; c < veldof; c++) {
          vel_array[k][j][i][c] =  VelnPress_array[k][j][i][c];
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(ctx->daScal,fields->fracpressure,&Press_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,fields->fracvelocity,&vel_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daFlow,fields->fracVelnPress,&VelnPress_array);CHKERRQ(ierr);
  /*
   ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"RHSvec.txt",&viewer);CHKERRQ(ierr);
   ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
   ierr = VecView(ctx->RHSFracVelP,viewer);CHKERRQ(ierr);
   
   ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"Matrix.txt",&viewer);CHKERRQ(ierr);
   ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
   ierr = MatView(ctx->JacFracVelP,viewer);CHKERRQ(ierr);
   
   */
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormFracSNESIJacobian"
extern PetscErrorCode FormFracSNESIJacobian(SNES snes,Vec VelnPress,Mat *Jac,Mat *Jacpre,MatStructure *str,void *user)
{
  PetscErrorCode ierr;
  VFCtx       *ctx=(VFCtx*)user;
  
  PetscFunctionBegin;
  *str = DIFFERENT_NONZERO_PATTERN;
  ierr = MatZeroEntries(*Jac);CHKERRQ(ierr);
  ierr = MatCopy(ctx->KFracVelP,*Jac,*str);
  ierr = MatAssemblyBegin(*Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (*Jac != *Jacpre) {
    ierr = MatAssemblyBegin(*Jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*Jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormFracSNESIFunction"
extern PetscErrorCode FormFracSNESIFunction(SNES snes,Vec VelnPress,Vec Residual,void *user)
{
  PetscErrorCode ierr;
  VFCtx                   *ctx=(VFCtx*)user;
  PetscReal               theta,timestepsize;
  PetscReal               dt_dot_theta,dt_dot_one_minus_theta;
  Vec             VecRHS;
  
  PetscFunctionBegin;
  
  timestepsize = ctx->flowprop.timestepsize;
  theta = ctx->flowprop.theta;
  dt_dot_theta = timestepsize*theta;
  dt_dot_one_minus_theta = timestepsize*(1.-theta);
  
  ierr = VecSet(Residual,0.0);CHKERRQ(ierr);
  ierr = VecDuplicate(ctx->RHSFracVelP,&VecRHS);CHKERRQ(ierr);
  ierr = VecCopy(ctx->RHSFracVelP,VecRHS);CHKERRQ(ierr);
  ierr = VecAXPBY(VecRHS,dt_dot_one_minus_theta,dt_dot_theta,ctx->RHSFracVelPpre);CHKERRQ(ierr);
  ierr = MatMultAdd(ctx->KFracVelPlhs,ctx->PreFracFlowFields,VecRHS,VecRHS);CHKERRQ(ierr);
  ierr = VecApplySNESVelocityBC(VecRHS,ctx->FracVelBCArray,&ctx->bcFracQ[0],ctx);CHKERRQ(ierr);
  ierr = MatMult(ctx->KFracVelP,VelnPress,Residual);CHKERRQ(ierr);
  ierr = VecAXPY(Residual,-1.0,VecRHS);CHKERRQ(ierr);
  ierr = VecDestroy(&VecRHS);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormFracMatricesnVector"
extern PetscErrorCode FormFracMatricesnVector(Mat K,Mat Klhs,Vec RHS,VFCtx *ctx, VFFields *fields)
{
  PetscErrorCode ierr;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       ek,ej,ei;
  PetscInt       i,j,k,l,ii;
  PetscInt       veldof = 3;
  PetscInt       c;
  PetscReal      ****coords_array;
  PetscReal      ****RHS_array;
  PetscReal      *RHS_local;
  Vec            RHS_localVec;
  PetscReal      hx,hy,hz;
  PetscReal      *KA_local,*KB_local,*KD_local,*KBTrans_local, *KS_local;
  PetscReal      beta_c,alpha_c,mu,gx,gy,gz;
  PetscReal      theta,timestepsize;
  PetscInt       nrow = ctx->e3D.nphix*ctx->e3D.nphiy*ctx->e3D.nphiz;
  MatStencil     *row,*row1;
  PetscReal      time_theta,time_one_minus_theta;
  MatStructure    flg;
  PetscReal      hwx, hwy, hwz;
  PetscReal      ***V_array;
  Vec            V_local;
  PetscReal      ***one_minus_v_array;
  Vec            One_minus_v_local;
  Vec            One_minus_V;
  PetscReal      ****U_array;
  Vec            U_local;
  Mat            K1,K2;
  PetscReal      Kw;
  
  
  PetscFunctionBegin;
  flg = SAME_NONZERO_PATTERN;
  beta_c = ctx->flowprop.beta;
  alpha_c = ctx->flowprop.alpha;
  theta = ctx->flowprop.theta;
  timestepsize = ctx->flowprop.timestepsize;
  time_theta = theta * timestepsize;
  time_one_minus_theta = -1.*(1.-theta)*timestepsize;
  mu     = ctx->flowprop.mu;
  gx     = ctx->flowprop.g[0];
  gy     = ctx->flowprop.g[1];
  gz     = ctx->flowprop.g[2];
  Kw     = ctx->flowprop.Kw;
  
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = MatZeroEntries(K);CHKERRQ(ierr);
  ierr = MatZeroEntries(Klhs);CHKERRQ(ierr);
  ierr = MatDuplicate(K,MAT_COPY_VALUES,&K1);CHKERRQ(ierr);
  ierr = MatDuplicate(K,MAT_COPY_VALUES,&K2);CHKERRQ(ierr);
  ierr = VecSet(RHS,0.);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daFlow,&RHS_localVec);CHKERRQ(ierr);
  ierr = VecSet(RHS_localVec,0.);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daFlow,RHS_localVec,&RHS_array);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(ctx->daScal,&One_minus_V);CHKERRQ(ierr);
  ierr = VecSet(One_minus_V,1.0);CHKERRQ(ierr);
  ierr = VecAXPY(One_minus_V,-1.0,fields->V);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ctx->daScal,&One_minus_v_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,One_minus_V,INSERT_VALUES,One_minus_v_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,One_minus_V,INSERT_VALUES,One_minus_v_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,One_minus_v_local,&V_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,One_minus_v_local,&one_minus_v_array);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daScal,&V_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,V_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,V_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,V_local,&V_array);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ctx->daVect,&U_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daVect,fields->U,INSERT_VALUES,U_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daVect,fields->U,INSERT_VALUES,U_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daVect,U_local,&U_array);CHKERRQ(ierr);
  
  ierr = PetscMalloc5(nrow*nrow,PetscReal,&KA_local,
                      nrow*nrow,PetscReal,&KB_local,
                      nrow*nrow,PetscReal,&KD_local,
                      nrow*nrow,PetscReal,&KS_local,
                      nrow*nrow,PetscReal,&KBTrans_local);CHKERRQ(ierr);
  ierr = PetscMalloc3(nrow,PetscReal,&RHS_local,
                      nrow,MatStencil,&row,
                      nrow,MatStencil,&row1);CHKERRQ(ierr);
  for (ek = zs; ek < zs+zm; ek++) {
    for (ej = ys; ej < ys+ym; ej++) {
      for (ei = xs; ei < xs+xm; ei++) {
        hx   = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
        hy   = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
        hz   = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
        ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
          //This computes the local contribution of the global A matrix
        ierr = FracFlow_MatA(KA_local,&ctx->e3D,ek,ej,ei,U_array,V_array);CHKERRQ(ierr);
        for (c = 0; c < veldof; c++) {
          ierr = FracFlow_MatB(KB_local,&ctx->e3D,ek,ej,ei,c,U_array,V_array);CHKERRQ(ierr);
          ierr = FracFlow_MatBTranspose(KBTrans_local,&ctx->e3D,ek,ej,ei,c,ctx->flowprop,U_array,V_array);CHKERRQ(ierr);
          for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
            for (j = 0; j < ctx->e3D.nphiy; j++) {
              for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                row[l].i  = ei+i;row[l].j = ej+j;row[l].k = ek+k;row[l].c = c;
                row1[l].i = ei+i;row1[l].j = ej+j;row1[l].k = ek+k;row1[l].c = 3;
              }
            }
          }
          ierr = MatSetValuesStencil(K1,nrow,row,nrow,row,KA_local,ADD_VALUES);CHKERRQ(ierr);
          ierr = MatSetValuesStencil(K1,nrow,row1,nrow,row,KB_local,ADD_VALUES);CHKERRQ(ierr);
          ierr = MatSetValuesStencil(K1,nrow,row,nrow,row1,KBTrans_local,ADD_VALUES);CHKERRQ(ierr);
        }
        ierr = FracFlow_MatS(KS_local,&ctx->e3D,ek,ej,ei,U_array,V_array);CHKERRQ(ierr);
        ierr = FracFlow_MatD(KD_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,U_array,V_array);CHKERRQ(ierr);
        for (l = 0; l < nrow*nrow; l++) {
          KS_local[l] = KS_local[l]/Kw;
        }
        for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
          for (j = 0; j < ctx->e3D.nphiy; j++) {
            for (i = 0; i < ctx->e3D.nphix; i++,l++) {
              row[l].i = ei+i;row[l].j = ej+j;row[l].k = ek+k;row[l].c = 3;
            }
          }
        }
        ierr = MatSetValuesStencil(K2,nrow,row,nrow,row,KS_local,ADD_VALUES);CHKERRQ(ierr);
        ierr = MatSetValuesStencil(K1,nrow,row,nrow,row,KD_local,ADD_VALUES);CHKERRQ(ierr);
        
        
        
        /*
         if(ei == 0 && ej == 0 && ek == 0){
         for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
         for (j = 0; j < ctx->e3D.nphiy; j++) {
         for (i = 0; i < ctx->e3D.nphix; i++,l++) {
         ierr = PetscPrintf(PETSC_COMM_WORLD,"ek = %d \t ej = %d \t ei = %d \t KA[%d]= %g \n",ek,ej,ei,l,KA_local[l]);CHKERRQ(ierr);
         
         }
         }
         }
         }
         */
        
        /*
         if(ei == nx/2 && ej == ny/2 && ek == nz/2){
         for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
         for (j = 0; j < ctx->e3D.nphiy; j++) {
         for (i = 0; i < ctx->e3D.nphix; i++,l++) {
         ierr = PetscPrintf(PETSC_COMM_WORLD,"ek = %d \t ej = %d \t ei = %d \t KA[%d]= %g \t cod[%d][%d][%d] = %g  \n",ek,ej,ei,l,KA_local[l],ek+k,ej+j,ei+i, crackopening_array[ek+k][ej+j][ei+i]);CHKERRQ(ierr);
         
         }
         }
         }
         
         }
         */
        
        
        /*Assembling the righthand side vector g*/
        ierr = FracFlow_Vecg(RHS_local,&ctx->e3D,ek,ej,ei,U_array,V_array);CHKERRQ(ierr);
        for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
          for (j = 0; j < ctx->e3D.nphiy; j++) {
            for (i = 0; i < ctx->e3D.nphix; i++,l++) {
              RHS_array[ek+k][ej+j][ei+i][3] += RHS_local[l]/timestepsize;
              
              /*
               if(ei == 0 && ej == 0 && ek == 0){
               ierr = PetscPrintf(PETSC_COMM_WORLD,"RHS[%d]= %g  \n",l,RHS_local[l]);CHKERRQ(ierr);}
               
               
               if(ei == nx/2 && ej == ny/2 && ek == nz/2){
               ierr = PetscPrintf(PETSC_COMM_WORLD,"RHS[%d]= %g  \n",l,RHS_local[l]);CHKERRQ(ierr);}
               */
            }
          }
        }
      }
    }
  }
  if(ctx->hasFlowWells){
    for(ii = 0; ii < ctx->numWells; ii++){
      for (ek = zs; ek < zs+zm; ek++) {
        for (ej = ys; ej < ys+ym; ej++) {
          for (ei = xs; ei < xs+xm; ei++) {
            if(
               ((coords_array[ek][ej][ei+1][0] >= ctx->well[ii].coords[0]) && (coords_array[ek][ej][ei][0] <= ctx->well[ii].coords[0] ))    &&
               ((coords_array[ek][ej+1][ei][1] >= ctx->well[ii].coords[1]) && (coords_array[ek][ej][ei][1] <= ctx->well[ii].coords[1] ))    &&
               ((coords_array[ek+1][ej][ei][2] >= ctx->well[ii].coords[2]) && (coords_array[ek][ej][ei][2] <= ctx->well[ii].coords[2] ))
               )
            {
              hx   = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
              hy   = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
              hz   = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
              ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
              hwx = (ctx->well[ii].coords[0]-coords_array[ek][ej][ei][0])/hx;
              hwy = (ctx->well[ii].coords[1]-coords_array[ek][ej][ei][1])/hy;
              hwz = (ctx->well[ii].coords[2]-coords_array[ek][ej][ei][2])/hz;
              if(ctx->well[ii].condition == RATE){
                ierr = VecApplyWellFlowRate(RHS_local,&ctx->e3D,ctx->well[ii].Qw,hwx,hwy,hwz,ek,ej,ei,one_minus_v_array);CHKERRQ(ierr);
                for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
                  for (j = 0; j < ctx->e3D.nphiy; j++) {
                    for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                      if(ctx->well[ii].type == INJECTOR){
                        RHS_array[ek+k][ej+j][ei+i][3] += RHS_local[l];
                      }
                      else if(ctx->well[ii].type == PRODUCER)
                      {
                        RHS_array[ek+k][ej+j][ei+i][3] += -1*RHS_local[l];
                      }
                    }
                  }
                }
              }
              else if(ctx->well[ii].condition == PRESSURE){
                
              }
              goto outer;
            }
          }
        }
      }
    outer:;
    }
  }
  ierr = MatAssemblyBegin(K2,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K2,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(K1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAXPY(K,time_theta,K1,flg);CHKERRQ(ierr);
  ierr = MatAXPY(K,1.0,K2,flg);CHKERRQ(ierr);
  ierr = MatAXPY(Klhs,time_one_minus_theta,K1,flg);CHKERRQ(ierr);
  ierr = MatAXPY(Klhs,1.0,K2,flg);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatApplyDirichletBC(K,ctx->daVect,&ctx->bcFracQ[0]);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,U_local,&U_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daVect,&U_local);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,V_local,&V_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&V_local);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,V_local,&one_minus_v_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&One_minus_v_local);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daFlow,RHS_localVec,&RHS_array);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(ctx->daFlow,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ctx->daFlow,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daFlow,&RHS_localVec);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = PetscFree5(KA_local,KB_local,KD_local,KS_local,KBTrans_local);CHKERRQ(ierr);
  ierr = PetscFree3(RHS_local,row,row1);CHKERRQ(ierr);
  ierr = VecDestroy(&One_minus_V);CHKERRQ(ierr);
  ierr = MatDestroy(&K1);CHKERRQ(ierr);
  ierr = MatDestroy(&K2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FracFlow_MatB"
extern PetscErrorCode FracFlow_MatB(PetscReal *KB_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscInt dof,PetscReal ****u_array,PetscReal ***v_array)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,l,c;
  PetscInt       ii,jj,kk;
  PetscReal      *dv_elem[3],*u_elem[3],*n_elem[3],*n_mag_elem;
  PetscInt       eg;  
  
  PetscFunctionBegin;
  ierr = PetscMalloc6(e->ng,PetscReal,&dv_elem[0],e->ng,PetscReal,&dv_elem[1],e->ng,PetscReal,&dv_elem[2],e->ng,PetscReal,&u_elem[0],e->ng,PetscReal,&u_elem[1],e->ng,PetscReal,&u_elem[2]);CHKERRQ(ierr);
  ierr = PetscMalloc4(e->ng,PetscReal,&n_elem[0],e->ng,PetscReal,&n_elem[1],e->ng,PetscReal,&n_elem[2],e->ng,PetscReal,&n_mag_elem);CHKERRQ(ierr);
  
  for (eg = 0; eg < e->ng; eg++){
    for(c = 0; c < 3; c++){
      dv_elem[c][eg] = 0;
      u_elem[c][eg] = 0;
      n_elem[c][eg] = 0;
    }
    n_mag_elem[eg] = 0.;
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          for (c = 0; c < 3; c++){
            dv_elem[c][eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
            u_elem[c][eg] += u_array[ek+k][ej+j][ei+i][c] * e->phi[k][j][i][eg];
          }
        }
      }
    }
  }
  for (eg = 0; eg < e->ng; eg++){
    n_mag_elem[eg] = sqrt((pow(dv_elem[0][eg],2))+(pow(dv_elem[1][eg],2))+(pow(dv_elem[2][eg],2)));
    for(c = 0; c < 3; c++){
      n_elem[0][eg] = dv_elem[0][eg]/n_mag_elem[eg];
      n_elem[1][eg] = dv_elem[1][eg]/n_mag_elem[eg];
      n_elem[2][eg] = dv_elem[2][eg]/n_mag_elem[eg];
    }
    {
      n_elem[0][eg] = n_elem[1][eg] = n_elem[2][eg] = n_mag_elem[eg] = 0;
    }
  }
  for (l = 0,k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (kk = 0; kk < e->nphiz; kk++) {
          for (jj = 0; jj < e->nphiy; jj++) {
            for (ii = 0; ii < e->nphix; ii++,l++) {
              KB_ele[l] = 0.;
              if (dof == 0){
                for (eg = 0; eg < e->ng; eg++) {
                  for(c = 0; c < 3; c++){
                    KB_ele[l] += (e->phi[kk][jj][ii][eg]*(-e->dphi[k][j][i][0][eg]*(1.-pow(n_elem[0][eg],2))
                                                          +e->dphi[k][j][i][1][eg]*n_elem[0][eg]*n_elem[1][eg]
                                                          +e->dphi[k][j][i][2][eg]*n_elem[0][eg]*n_elem[2][eg])*pow((dv_elem[c][eg]*u_elem[c][eg]),1))*e->weight[eg]
                    +0.5*e->dphi[k][j][i][0][eg]*pow((dv_elem[c][eg]*u_elem[c][eg]),1)*e->weight[eg];
                  }
                }
              }
              if (dof == 1){
                for (eg = 0; eg < e->ng; eg++) {
                  for(c = 0; c < 3; c++){
                    KB_ele[l] += (e->phi[kk][jj][ii][eg]*(+e->dphi[k][j][i][0][eg]*n_elem[0][eg]*n_elem[1][eg]
                                                          -e->dphi[k][j][i][1][eg]*(1.-pow(n_elem[1][eg],2))
                                                          +e->dphi[k][j][i][2][eg]*n_elem[1][eg]*n_elem[2][eg])*pow((dv_elem[c][eg]*u_elem[c][eg]),1))*e->weight[eg]
                    +0.5*e->dphi[k][j][i][1][eg]*pow((dv_elem[c][eg]*u_elem[c][eg]),1)*e->weight[eg];
                  }
                }
              }
              if (dof == 2){
                for (eg = 0; eg < e->ng; eg++) {
                  for(c = 0; c < 3; c++){
                    KB_ele[l] += (e->phi[kk][jj][ii][eg]*(+e->dphi[k][j][i][0][eg]*n_elem[0][eg]*n_elem[2][eg]
                                                          +e->dphi[k][j][i][1][eg]*n_elem[1][eg]*n_elem[2][eg]
                                                          -e->dphi[k][j][i][2][eg]*(1.-pow(n_elem[2][eg],2)))*pow((dv_elem[c][eg]*u_elem[c][eg]),1))*e->weight[eg]
                    +0.5*e->dphi[k][j][i][2][eg]*pow((dv_elem[c][eg]*u_elem[c][eg]),1)*e->weight[eg];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  
  
  
  
  
  PetscReal velem;
  velem = 0;
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        velem += v_array[ek+k][ej+j][ei+i];
      }
    }
  }
  velem = velem/8.;
  for(l = 0; l < 64; l++){
    if(velem <= 0.2)
      KB_ele[l] += 1e-6;
  }
  
  
  
  ierr = PetscFree6(dv_elem[0],dv_elem[1],dv_elem[2],u_elem[0],u_elem[1],u_elem[2]);CHKERRQ(ierr);
  ierr = PetscFree4(n_elem[0],n_elem[1],n_elem[2],n_mag_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FracFlow_MatD"
extern PetscErrorCode FracFlow_MatD(PetscReal *Kd_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,VFFlowProp flowpropty,PetscReal ****u_array,PetscReal ***v_array)
{
  PetscErrorCode          ierr;
  PetscInt                i,j,k,l,c;
  PetscInt                ii,jj,kk;
  PetscReal               *dv_elem[3],*u_elem[3],*n_elem[3],*n_mag_elem;
  PetscInt                eg;
  PetscReal               mu;
  
  ierr = PetscMalloc6(e->ng,PetscReal,&dv_elem[0],e->ng,PetscReal,&dv_elem[1],e->ng,PetscReal,&dv_elem[2],e->ng,PetscReal,&u_elem[0],e->ng,PetscReal,&u_elem[1],e->ng,PetscReal,&u_elem[2]);CHKERRQ(ierr);
  ierr = PetscMalloc4(e->ng,PetscReal,&n_elem[0],e->ng,PetscReal,&n_elem[1],e->ng,PetscReal,&n_elem[2],e->ng,PetscReal,&n_mag_elem);CHKERRQ(ierr);
  mu     = flowpropty.mu;
  for (eg = 0; eg < e->ng; eg++){
    for(c = 0; c < 3; c++){
      dv_elem[c][eg] = 0;
      u_elem[c][eg] = 0;
      n_elem[c][eg] = 0;
    }
    n_mag_elem[eg] = 0.;
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          for (c = 0; c < 3; c++){
            dv_elem[c][eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
            u_elem[c][eg] += u_array[ek+k][ej+j][ei+i][c] * e->phi[k][j][i][eg];
          }
        }
      }
    }
  }
  for (eg = 0; eg < e->ng; eg++){
    n_mag_elem[eg] = sqrt((pow(dv_elem[0][eg],2))+(pow(dv_elem[1][eg],2))+(pow(dv_elem[2][eg],2)));
    for(c = 0; c < 3; c++){
      n_elem[0][eg] = dv_elem[0][eg]/n_mag_elem[eg];
      n_elem[1][eg] = dv_elem[1][eg]/n_mag_elem[eg];
      n_elem[2][eg] = dv_elem[2][eg]/n_mag_elem[eg];
    }
    if((ierr = PetscIsInfOrNanScalar(n_elem[0][eg])) || (ierr = PetscIsInfOrNanScalar(n_elem[1][eg])) || (ierr = PetscIsInfOrNanScalar(n_elem[2][eg])) )
    {
      n_elem[0][eg] = n_elem[1][eg] = n_elem[2][eg] = n_mag_elem[eg] = 0;
    }
  }
  for (l = 0,k = 0; k < e->nphiz; k++) {
		for (j = 0; j < e->nphiy; j++) {
			for (i = 0; i < e->nphix; i++) {
				for (kk = 0; kk < e->nphiz; kk++) {
					for (jj = 0; jj < e->nphiy; jj++) {
						for (ii = 0; ii < e->nphix; ii++,l++) {
							Kd_ele[l] = 0.;
							for (eg = 0; eg < e->ng; eg++) {
                Kd_ele[l] += ((1.-pow(n_elem[0][eg],2))*e->dphi[kk][jj][ii][0][eg]
                              -e->dphi[kk][jj][ii][1][eg]*n_elem[0][eg]*n_elem[1][eg]
                              -e->dphi[k][j][i][2][eg]*n_elem[0][eg]*n_elem[2][eg])*e->dphi[k][j][i][0][eg]*4*(pow((u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg]),3))*n_mag_elem[eg]/(24.0*mu)*e->weight[eg]
                +(-e->dphi[kk][jj][ii][0][eg]*n_elem[0][eg]*n_elem[1][eg]
                  +(1.-pow(n_elem[1][eg],2))*e->dphi[kk][jj][ii][1][eg]
                  -e->dphi[k][j][i][2][eg]*n_elem[1][eg]*n_elem[2][eg])*e->dphi[k][j][i][1][eg]*4*(pow((u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg]),3))*n_mag_elem[eg]/(24.0*mu)*e->weight[eg]
                +(-e->dphi[kk][jj][ii][0][eg]*n_elem[0][eg]*n_elem[2][eg]
                  -e->dphi[k][j][i][1][eg]*n_elem[1][eg]*n_elem[2][eg]
                  +(1.-pow(n_elem[2][eg],2))*e->dphi[kk][jj][ii][2][eg])*e->dphi[k][j][i][2][eg]*4*(pow((u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg]),3))*n_mag_elem[eg]/(24.0*mu)*e->weight[eg];
              }
						}
					}
				}
			}
		}
	}
  
  
  
  
  PetscReal velem;
  velem = 0;
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        velem += v_array[ek+k][ej+j][ei+i];
      }
    }
  }
  velem = velem/8.;
  for(l = 0; l < 64; l++){
    if(velem <= 0.2)
      Kd_ele[l] += 1e-6;
  }
  
  
  ierr = PetscFree6(dv_elem[0],dv_elem[1],dv_elem[2],u_elem[0],u_elem[1],u_elem[2]);CHKERRQ(ierr);
  ierr = PetscFree4(n_elem[0],n_elem[1],n_elem[2],n_mag_elem);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FracFlow_MatBTranspose"
extern PetscErrorCode FracFlow_MatBTranspose(PetscReal *KB_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscInt dof,VFFlowProp flowpropty,PetscReal ****u_array,PetscReal ***v_array)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,l,c;
  PetscInt       ii,jj,kk;
  PetscReal      *dv_elem[3],*u_elem[3],*n_elem[3],*n_mag_elem;
  PetscInt       eg;
  PetscReal      mu;
  
  PetscFunctionBegin;
  ierr = PetscMalloc6(e->ng,PetscReal,&dv_elem[0],e->ng,PetscReal,&dv_elem[1],e->ng,PetscReal,&dv_elem[2],e->ng,PetscReal,&u_elem[0],e->ng,PetscReal,&u_elem[1],e->ng,PetscReal,&u_elem[2]);CHKERRQ(ierr);
  ierr = PetscMalloc4(e->ng,PetscReal,&n_elem[0],e->ng,PetscReal,&n_elem[1],e->ng,PetscReal,&n_elem[2],e->ng,PetscReal,&n_mag_elem);CHKERRQ(ierr);
  mu     = flowpropty.mu;
  for (eg = 0; eg < e->ng; eg++){
    for(c = 0; c < 3; c++){
      dv_elem[c][eg] = 0;
      u_elem[c][eg] = 0;
      n_elem[c][eg] = 0;
    }
    n_mag_elem[eg] = 0.;
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          for (c = 0; c < 3; c++){
            dv_elem[c][eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
            u_elem[c][eg] += u_array[ek+k][ej+j][ei+i][c] * e->phi[k][j][i][eg];
          }
        }
      }
    }
  }
  for (eg = 0; eg < e->ng; eg++){
    n_mag_elem[eg] = sqrt((pow(dv_elem[0][eg],2))+(pow(dv_elem[1][eg],2))+(pow(dv_elem[2][eg],2)));
    for(c = 0; c < 3; c++){
      n_elem[0][eg] = dv_elem[0][eg]/n_mag_elem[eg];
      n_elem[1][eg] = dv_elem[1][eg]/n_mag_elem[eg];
      n_elem[2][eg] = dv_elem[2][eg]/n_mag_elem[eg];
    }
    if((ierr = PetscIsInfOrNanScalar(n_elem[0][eg])) || (ierr = PetscIsInfOrNanScalar(n_elem[1][eg])) || (ierr = PetscIsInfOrNanScalar(n_elem[2][eg])) )
    {
      n_elem[0][eg] = n_elem[1][eg] = n_elem[2][eg] = n_mag_elem[eg] = 0;
    }
  }
  for (l = 0,k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (kk = 0; kk < e->nphiz; kk++) {
          for (jj = 0; jj < e->nphiy; jj++) {
            for (ii = 0; ii < e->nphix; ii++,l++) {
              KB_ele[l] = 0.;
              if (dof == 0){
                for (eg = 0; eg < e->ng; eg++) {
                  KB_ele[l] += ((1.-pow(n_elem[0][eg],2))*e->dphi[kk][jj][ii][0][eg]
                                -e->dphi[kk][jj][ii][1][eg]*n_elem[0][eg]*n_elem[1][eg]
                                -e->dphi[k][j][i][2][eg]*n_elem[0][eg]*n_elem[2][eg])*e->phi[k][j][i][eg]*4*(pow((u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg]),3))*n_mag_elem[eg]/(24.0*mu)*e->weight[eg];
                }
              }
              if (dof == 1){
                for (eg = 0; eg < e->ng; eg++) {
                  KB_ele[l] += (-e->dphi[kk][jj][ii][0][eg]*n_elem[0][eg]*n_elem[1][eg]
                                +(1.-pow(n_elem[1][eg],2))*e->dphi[kk][jj][ii][1][eg]
                                -e->dphi[k][j][i][2][eg]*n_elem[1][eg]*n_elem[2][eg])*e->phi[k][j][i][eg]*4*(pow((u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg]),3))*n_mag_elem[eg]/(24.0*mu)*e->weight[eg];
                }
              }
              if (dof == 2){
                for (eg = 0; eg < e->ng; eg++) {
                  KB_ele[l] += (-e->dphi[kk][jj][ii][0][eg]*n_elem[0][eg]*n_elem[2][eg]
                                -e->dphi[k][j][i][1][eg]*n_elem[1][eg]*n_elem[2][eg]
                                +(1.-pow(n_elem[2][eg],2))*e->dphi[kk][jj][ii][2][eg])*e->phi[k][j][i][eg]*4*(pow((u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg]),3))*n_mag_elem[eg]/(24.0*mu)*e->weight[eg];
                }
              }
            }
          }
        }
      }
    }
  }
  
  
  PetscReal velem;
  velem = 0;
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        velem += v_array[ek+k][ej+j][ei+i];
      }
    }
  }
  velem = velem/8.;
  for(l = 0; l < 64; l++){
    if(velem <= 0.2)
      KB_ele[l] += 1e-6;
  }
  
  
  
  
  ierr = PetscFree6(dv_elem[0],dv_elem[1],dv_elem[2],u_elem[0],u_elem[1],u_elem[2]);CHKERRQ(ierr);
  ierr = PetscFree4(n_elem[0],n_elem[1],n_elem[2],n_mag_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FracFlow_MatS"
extern PetscErrorCode FracFlow_MatS(PetscReal *S_local,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ****u_array,PetscReal ***v_array)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,l,c;
  PetscInt       ii,jj,kk;
  PetscReal      *dv_elem[3],*u_elem[3];
  PetscInt       eg;
  
  PetscFunctionBegin;
  ierr = PetscMalloc6(e->ng,PetscReal,&dv_elem[0],e->ng,PetscReal,&dv_elem[1],e->ng,PetscReal,&dv_elem[2],e->ng,PetscReal,&u_elem[0],e->ng,PetscReal,&u_elem[1],e->ng,PetscReal,&u_elem[2]);CHKERRQ(ierr);
  for (eg = 0; eg < e->ng; eg++){
    for(c = 0; c < 3; c++){
      dv_elem[c][eg] = 0;
      u_elem[c][eg] = 0;
    }
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          for (c = 0; c < 3; c++){
            dv_elem[c][eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
            u_elem[c][eg] += u_array[ek+k][ej+j][ei+i][c] * e->phi[k][j][i][eg];
          }
        }
      }
    }
  }
  for (l = 0,k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (kk = 0; kk < e->nphiz; kk++) {
          for (jj = 0; jj < e->nphiy; jj++) {
            for (ii = 0; ii < e->nphix; ii++,l++) {
              S_local[l] = 0;
              for (eg = 0; eg < e->ng; eg++) {
                for(c = 0; c < 3; c++){
                  S_local[l] += e->phi[k][j][i][eg]*e->phi[kk][jj][ii][eg]*dv_elem[c][eg]*u_elem[c][eg]*e->weight[eg];
                }
              }
            }
          }
        }
      }
    }
  }
  
  
  PetscReal velem;
  velem = 0;
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        velem += v_array[ek+k][ej+j][ei+i];
      }
    }
  }
  velem = velem/8.;
  for(l = 0; l < 64; l++){
    if(velem <= 0.2)
      S_local[l] += 1e-6;
  }
  
  ierr = PetscFree6(dv_elem[0],dv_elem[1],dv_elem[2],u_elem[0],u_elem[1],u_elem[2]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FracFlow_MatA"
extern PetscErrorCode FracFlow_MatA(PetscReal *A_local,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ****u_array,PetscReal ***v_array)
{
  PetscErrorCode          ierr;
  PetscInt                i,j,k,l,c;
  PetscInt                ii,jj,kk;
  PetscReal               *dv_elem[3],*u_elem[3];
  PetscInt                eg;
  
  PetscFunctionBegin;
  ierr = PetscMalloc6(e->ng,PetscReal,&dv_elem[0],e->ng,PetscReal,&dv_elem[1],e->ng,PetscReal,&dv_elem[2],e->ng,PetscReal,&u_elem[0],e->ng,PetscReal,&u_elem[1],e->ng,PetscReal,&u_elem[2]);CHKERRQ(ierr);
  for (eg = 0; eg < e->ng; eg++){
    for(c = 0; c < 3; c++){
      dv_elem[c][eg] = 0;
      u_elem[c][eg] = 0;
    }
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          for (c = 0; c < 3; c++){
            dv_elem[c][eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
            u_elem[c][eg] += u_array[ek+k][ej+j][ei+i][c] * e->phi[k][j][i][eg];
          }
        }
      }
    }
  }
  for (l = 0,k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (kk = 0; kk < e->nphiz; kk++) {
          for (jj = 0; jj < e->nphiy; jj++) {
            for (ii = 0; ii < e->nphix; ii++,l++) {
              A_local[l] = 0;
              for (eg = 0; eg < e->ng; eg++) {
                for (c = 0; c < 3; c++){
                  A_local[l] += 0.5*e->phi[k][j][i][eg]*e->phi[kk][jj][ii][eg]*dv_elem[c][eg]*u_elem[c][eg]*e->weight[eg];
                }
              }
            }
          }
        }
      }
    }
  }
  
  PetscReal velem;
  velem = 0;
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        velem += v_array[ek+k][ej+j][ei+i];
        }
      }
    }
  velem = velem/8.;
  for(l = 0; l < 64; l++){
    if(velem <= 0.2)
    A_local[l] += 1e-6;
  }
  
  ierr = PetscFree6(dv_elem[0],dv_elem[1],dv_elem[2],u_elem[0],u_elem[1],u_elem[2]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FracFlow_Vecg"
extern PetscErrorCode FracFlow_Vecg(PetscReal *Kg_local,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei, PetscReal ****u_array, PetscReal ***v_array)
{
  PetscErrorCode          ierr;
  PetscInt                i,j,k,l,c;
  PetscReal               *dv_elem[3],*u_elem[3];
  PetscInt                eg;
  
  PetscFunctionBegin;
  ierr = PetscMalloc6(e->ng,PetscReal,&dv_elem[0],e->ng,PetscReal,&dv_elem[1],e->ng,PetscReal,&dv_elem[2],e->ng,PetscReal,&u_elem[0],e->ng,PetscReal,&u_elem[1],e->ng,PetscReal,&u_elem[2]);CHKERRQ(ierr);
  for (eg = 0; eg < e->ng; eg++){
    for(c = 0; c < 3; c++){
      dv_elem[c][eg] = 0;
      u_elem[c][eg] = 0;
    }
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          for (c = 0; c < 3; c++){
            dv_elem[c][eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
            u_elem[c][eg] += u_array[ek+k][ej+j][ei+i][c] * e->phi[k][j][i][eg];
          }
        }
      }
    }
  }
  for (l = 0,k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++,l++) {
        Kg_local[l] = 0.;
        for (eg = 0; eg < e->ng; eg++) {
          for(c = 0; c < 3; c++){
            Kg_local[l] += -1*(dv_elem[c][eg]*u_elem[c][eg]*e->phi[k][j][i][eg]*e->weight[eg]);
          }
        }
      }
    }
  }
  ierr = PetscFree6(dv_elem[0],dv_elem[1],dv_elem[2],u_elem[0],u_elem[1],u_elem[2]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}