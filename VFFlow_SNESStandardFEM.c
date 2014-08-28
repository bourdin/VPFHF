/*
 VFFlow_SNESStandardFEM.c
 
 (c) 2011-2013 C. Chukwudozie, LSU
 */

#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFWell.h"
#include "VFHeat_SNESFEM.h"
#include "VFFlow_KSPMixedFEM.h"
#include "VFFlow_SNESStandardFEM.h"
#include "VFFlow.h"
#include "VFPermfield.h"


/*
 VFFlow_SNESStandardFEM
 */

/*
 ################################################################################################################
 SNES ROUTINE
 ################################################################################################################
 */
#undef __FUNCT__
#define __FUNCT__ "VFFlow_SNESStandardFEMInitialize"
extern PetscErrorCode VFFlow_SNESStandardFEMInitialize(VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode ierr;
  PetscFunctionBegin;
	
	ierr = SNESCreate(PETSC_COMM_WORLD,&ctx->snesP);CHKERRQ(ierr);
  ierr = SNESAppendOptionsPrefix(ctx->snesP,"FlowStSnes_");CHKERRQ(ierr);
  ierr = SNESSetFromOptions(ctx->snesP);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFFlow_SNESStandardFEMFinalize"
extern PetscErrorCode VFFlow_SNESStandardFEMFinalize(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
  ierr = SNESDestroy(&ctx->snesP);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "VF_FormFlowStandardFEMIFunction"
extern PetscErrorCode VF_FormFlowStandardFEMIFunction(SNES snes,Vec Pressure,Vec Func,void *user)
{
	PetscErrorCode ierr;
	VFCtx			 *ctx=(VFCtx*)user;
	Vec         VecRHS;
	PetscReal		theta;
	PetscReal		one_minus_theta;
	
	PetscFunctionBegin;
	ierr = VecSet(Func,0.0);CHKERRQ(ierr);
	theta = ctx->flowprop.theta;
	one_minus_theta = (1.-theta);
	ierr = VecDuplicate(ctx->RHSP,&VecRHS);CHKERRQ(ierr);
	ierr = VecCopy(ctx->RHSP,VecRHS);CHKERRQ(ierr);
	ierr = VecAXPBY(VecRHS,one_minus_theta,theta,ctx->RHSPpre);CHKERRQ(ierr);
	ierr = MatMultAdd(ctx->KPlhs,ctx->pressure_old,VecRHS,VecRHS);CHKERRQ(ierr);
  
  ierr = VecApplyDirichletBC(VecRHS,ctx->PresBCArray,&ctx->bcP[0]);CHKERRQ(ierr);

  
	ierr = MatMult(ctx->KP,Pressure,Func);CHKERRQ(ierr);
  ierr = VecAXPY(Func,-1.0,VecRHS);CHKERRQ(ierr);
	ierr = VecDestroy(&VecRHS);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FractureFlowVelocityCompute_local"
extern PetscErrorCode FractureFlowVelocityCompute_local(PetscReal ****cellvelocityrate_array, PetscReal ***press_array,PetscReal ****u_array, PetscReal ***v_array, VFFlowProp *flowpropty, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e)
{
	PetscErrorCode ierr;
	PetscInt		i, j, k, c;
	PetscInt		eg;
  PetscReal   mu;
  PetscReal		element_vol = 0;
  PetscReal   *gradpress_elem[3];
  PetscReal         *u_elem[3];
  PetscReal         *n_elem[3],*v_mag_elem;

  
	PetscFunctionBegin;
  cellvelocityrate_array[ek][ej][ei][0] = 0;
  cellvelocityrate_array[ek][ej][ei][1] = 0;
  cellvelocityrate_array[ek][ej][ei][2] = 0;
  mu      = flowpropty->mu;
  ierr = PetscMalloc7(e->ng,&n_elem[0],
                      e->ng,&n_elem[1],
                      e->ng,&n_elem[2],
                      e->ng,&u_elem[0],
                      e->ng,&u_elem[1],
                      e->ng,&u_elem[2],
                      e->ng,&v_mag_elem);CHKERRQ(ierr);
  
	ierr = PetscMalloc3(e->ng,&gradpress_elem[0],e->ng,&gradpress_elem[1],e->ng,&gradpress_elem[2]);CHKERRQ(ierr);
  for (eg = 0; eg < e->ng; eg++){
    for (c = 0; c < 3; c++){
      u_elem[c][eg] = 0.;
      n_elem[c][eg] = 0;
      gradpress_elem[c][eg] = 0.;
    }
    v_mag_elem[eg] = 0.;
    element_vol += e->weight[eg];
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          for (c= 0; c < 3; c++){
            n_elem[c][eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
            u_elem[c][eg] += u_array[ek+k][ej+j][ei+i][c]*e->phi[k][j][i][eg];
            gradpress_elem[c][eg] += press_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
          }
        }
      }
    }
  }
  for (eg = 0; eg < e->ng; eg++){
    v_mag_elem[eg] = sqrt((pow(n_elem[0][eg],2))+(pow(n_elem[1][eg],2))+(pow(n_elem[2][eg],2)));
    for(c = 0; c < 3; c++){
      n_elem[c][eg] = n_elem[c][eg]/v_mag_elem[eg];
    }
    if((PetscIsInfOrNanScalar(n_elem[0][eg])) || (PetscIsInfOrNanScalar(n_elem[1][eg])) || (PetscIsInfOrNanScalar(n_elem[2][eg])) )
    {
      n_elem[0][eg] = n_elem[1][eg] = n_elem[2][eg] = v_mag_elem[eg] = 0;
    }
  }
  for (eg = 0; eg < e->ng; eg++) {
    cellvelocityrate_array[ek][ej][ei][0] += -((1.-(pow(n_elem[0][eg],2)))*gradpress_elem[0][eg]
                                              -gradpress_elem[1][eg]*n_elem[0][eg]*n_elem[1][eg]
                                              -gradpress_elem[2][eg]*n_elem[0][eg]*n_elem[2][eg])*4*(pow((u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg]),3))*v_mag_elem[eg]*e->weight[eg]/mu;
    
    cellvelocityrate_array[ek][ej][ei][1] += -(-gradpress_elem[0][eg]*n_elem[0][eg]*n_elem[1][eg]
                                               +(1.-pow(n_elem[1][eg],2))*gradpress_elem[1][eg]
                                               -gradpress_elem[2][eg]*n_elem[1][eg]*n_elem[2][eg])*4*(pow((u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg]),3))*v_mag_elem[eg]*e->weight[eg]/mu;
    
    cellvelocityrate_array[ek][ej][ei][2] += -(-gradpress_elem[0][eg]*n_elem[0][eg]*n_elem[2][eg]
                                               -gradpress_elem[1][eg]*n_elem[1][eg]*n_elem[2][eg]
                                               +(1.-pow(n_elem[2][eg],2))*gradpress_elem[2][eg])*4*(pow((u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg]),3))*v_mag_elem[eg]*e->weight[eg]/mu;
  }
  cellvelocityrate_array[ek][ej][ei][0] = cellvelocityrate_array[ek][ej][ei][0]/element_vol;
  cellvelocityrate_array[ek][ej][ei][1] = cellvelocityrate_array[ek][ej][ei][1]/element_vol;
  cellvelocityrate_array[ek][ej][ei][2] = cellvelocityrate_array[ek][ej][ei][2]/element_vol;
	ierr = PetscFree3(gradpress_elem[0],gradpress_elem[1],gradpress_elem[2]);CHKERRQ(ierr);
  ierr = PetscFree7(n_elem[0],n_elem[1],n_elem[2],u_elem[0],u_elem[1],u_elem[2],v_mag_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FractureFlowVelocityCompute"
extern PetscErrorCode FractureFlowVelocityCompute(VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode  ierr;
	PetscInt		    ek, ej, ei;
	PetscInt		    xs,xm,nx;
	PetscInt		    ys,ym,ny;
	PetscInt		    zs,zm,nz;
	PetscReal		    hx,hy,hz;
	PetscReal       ****coords_array;
  Vec             cellVelocity;
  Vec             cellVelocity_local;
	PetscReal       ****cellVelocity_array;
  PetscReal       ***press_array;
	Vec             press_local;
  PetscReal       ****u_array;
	Vec             u_local;
  PetscReal       ***v_array;
	Vec             v_local;
  
  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daScal,&press_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,fields->pressure,INSERT_VALUES,press_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->pressure,INSERT_VALUES,press_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,press_local,&press_array);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,ctx->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,ctx->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);

  
  
  ierr = DMCreateGlobalVector(ctx->daVectCell,&cellVelocity);CHKERRQ(ierr);
  ierr = VecSet(cellVelocity,0.0);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daVectCell,&cellVelocity_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVectCell,cellVelocity,INSERT_VALUES,cellVelocity_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVectCell,cellVelocity,INSERT_VALUES,cellVelocity_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVectCell,cellVelocity_local,&cellVelocity_array);CHKERRQ(ierr);
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = VFCartFEElement3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
        ierr = FractureFlowVelocityCompute_local(cellVelocity_array, press_array, u_array,v_array, &ctx->flowprop, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);
        
			}
		}
	}
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,press_local,&press_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&press_local);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVectCell,cellVelocity_local,&cellVelocity_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVectCell,&cellVelocity_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daVectCell,cellVelocity_local,ADD_VALUES,cellVelocity);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daVectCell,cellVelocity_local,ADD_VALUES,cellVelocity);CHKERRQ(ierr);
  ierr = VecSet(fields->fracvelocity,0.);CHKERRQ(ierr);
	ierr = CellToNodeInterpolation(fields->fracvelocity,cellVelocity,ctx); CHKERRQ(ierr);
  ierr = VecDestroy(&cellVelocity);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_FlowStandardFEMSNESSolve"
extern PetscErrorCode VF_FlowStandardFEMSNESSolve(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode     ierr;
	SNESConvergedReason reason;
	PetscInt           its;
	PetscReal           Pmin,Pmax;
  
	PetscFunctionBegin;
	ierr = VecCopy(fields->V,ctx->V);CHKERRQ(ierr);
	ierr = VecCopy(fields->U,ctx->U);CHKERRQ(ierr);
  ierr = VF_FormFlowStandardFEMMatricesnVectors(ctx->KP,ctx->KPlhs,ctx->RHSP,fields,ctx);CHKERRQ(ierr);
	ierr = SNESSetFunction(ctx->snesP,ctx->PFunct,VF_FormFlowStandardFEMIFunction,ctx);CHKERRQ(ierr);
  ierr = SNESSetJacobian(ctx->snesP,ctx->JacP,ctx->JacP,VF_FormFlowStandardFEMIJacobian,ctx);CHKERRQ(ierr);
	if (ctx->verbose > 1) {
		ierr = SNESMonitorSet(ctx->snesP,FEMSNESMonitor,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	}
  ierr = SNESSolve(ctx->snesP,PETSC_NULL,fields->pressure);CHKERRQ(ierr);
	ierr = SNESGetConvergedReason(ctx->snesP,&reason);CHKERRQ(ierr);
  ierr = FlowVelocityCompute(ctx,fields);CHKERRQ(ierr);
  if (ctx->FractureFlowCoupling) {
    ierr = FractureFlowVelocityCompute(ctx,fields);CHKERRQ(ierr);
  }
	if (reason < 0) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"[ERROR] snes_StandardFlowSolver diverged with reason %d\n",(int)reason);CHKERRQ(ierr);
	} else {
		ierr = SNESGetIterationNumber(ctx->snesP,&its);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"      snes_StandardFlowSolver converged in %d iterations %d.\n",(int)its,(int)reason);CHKERRQ(ierr);
	}
  ierr = VecMin(fields->pressure,PETSC_NULL,&Pmin);CHKERRQ(ierr);
	ierr = VecMax(fields->pressure,PETSC_NULL,&Pmax);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"      Pressure min / max:     %e %e\n",Pmin,Pmax);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FlowVelocityCompute"
extern PetscErrorCode FlowVelocityCompute(VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode  ierr;
	PetscInt		    ek, ej, ei;
	PetscInt		    xs,xm,nx;
	PetscInt		    ys,ym,ny;
	PetscInt		    zs,zm,nz;
	PetscReal		    hx,hy,hz;
	PetscReal       ****coords_array;
  Vec             cellVelocity;
  Vec             cellVelocity_local;
	PetscReal       ****cellVelocity_array;
  PetscReal       ***press_array;
	Vec             press_local;
  PetscReal       ****perm_array;
	Vec             perm_local;
  PetscReal       ***v_array;
	Vec             v_local;
  
  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ctx->daScal,&press_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,fields->pressure,INSERT_VALUES,press_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->pressure,INSERT_VALUES,press_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,press_local,&press_array);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVFperm,fields->vfperm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVFperm,fields->vfperm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVFperm,perm_local,&perm_array);
  ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx->daVectCell,&cellVelocity);CHKERRQ(ierr);
  ierr = VecSet(cellVelocity,0.0);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daVectCell,&cellVelocity_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVectCell,cellVelocity,INSERT_VALUES,cellVelocity_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVectCell,cellVelocity,INSERT_VALUES,cellVelocity_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVectCell,cellVelocity_local,&cellVelocity_array);CHKERRQ(ierr);
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = VFCartFEElement3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
        ierr = FlowVelocityCompute_local(cellVelocity_array, press_array, perm_array, v_array, &ctx->flowprop, ek, ej, ei, &ctx->e3D);CHKERRQ(ierr);
			}
		}
	}
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,press_local,&press_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&press_local);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVectCell,cellVelocity_local,&cellVelocity_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVectCell,&cellVelocity_local);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daVectCell,cellVelocity_local,ADD_VALUES,cellVelocity);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daVectCell,cellVelocity_local,ADD_VALUES,cellVelocity);CHKERRQ(ierr);
  ierr = VecSet(fields->velocity,0.);CHKERRQ(ierr);
	ierr = CellToNodeInterpolation(fields->velocity,cellVelocity,ctx); CHKERRQ(ierr);
  ierr = VecDestroy(&cellVelocity);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FlowVelocityCompute_local"
extern PetscErrorCode FlowVelocityCompute_local(PetscReal ****cellvelocityrate_array, PetscReal ***press_array,PetscReal ****perm_array, PetscReal ***v_array, VFFlowProp *flowpropty, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e)
{
	PetscErrorCode ierr;
	PetscInt		i, j, k, c;
	PetscInt		eg;
  PetscReal   mu;
  PetscReal		element_vol = 0;
  PetscReal   *gradpress_elem[3];
  PetscReal   *v_elem;
  
	PetscFunctionBegin;
  cellvelocityrate_array[ek][ej][ei][0] = 0;
  cellvelocityrate_array[ek][ej][ei][1] = 0;
  cellvelocityrate_array[ek][ej][ei][2] = 0;
  mu      = flowpropty->mu;
	ierr = PetscMalloc4(e->ng,&gradpress_elem[0],e->ng,&gradpress_elem[1],e->ng,&gradpress_elem[2],e->ng,&v_elem);CHKERRQ(ierr);
	for (eg = 0; eg < e->ng; eg++){
    for(c = 0; c < 3; c++){
      gradpress_elem[c][eg] = 0.;
    }
    v_elem[eg] = 0.;
    element_vol += e->weight[eg];
	}
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          for (c = 0; c < 3; c++){
            gradpress_elem[c][eg] += press_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
          }
          v_elem[eg] += v_array[ek+k][ej+j][ei+i]*e->phi[k][j][i][eg];
        }
      }
    }
  }
  for (eg = 0; eg < e->ng; eg++) {
    cellvelocityrate_array[ek][ej][ei][0] += -pow(v_elem[eg],2)*(perm_array[ek][ej][ei][0]*gradpress_elem[0][eg])*e->weight[eg]/mu;
    cellvelocityrate_array[ek][ej][ei][1] += -pow(v_elem[eg],2)*(perm_array[ek][ej][ei][1]*gradpress_elem[1][eg])*e->weight[eg]/mu;
    cellvelocityrate_array[ek][ej][ei][2] += -pow(v_elem[eg],2)*(perm_array[ek][ej][ei][2]*gradpress_elem[2][eg])*e->weight[eg]/mu;
  }
  cellvelocityrate_array[ek][ej][ei][0] = cellvelocityrate_array[ek][ej][ei][0]/element_vol;
  cellvelocityrate_array[ek][ej][ei][1] = cellvelocityrate_array[ek][ej][ei][1]/element_vol;
  cellvelocityrate_array[ek][ej][ei][2] = cellvelocityrate_array[ek][ej][ei][2]/element_vol;
	ierr = PetscFree4(gradpress_elem[0],gradpress_elem[1],gradpress_elem[2],v_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_FormFlowStandardFEMMatricesnVectors"
extern PetscErrorCode VF_FormFlowStandardFEMMatricesnVectors(Mat K,Mat Krhs,Vec RHS,VFFields * fields,VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       ek,ej,ei;
  PetscInt       i,j,k,l;
  PetscReal      ****perm_array;
  PetscReal      ****coords_array;
  PetscReal      ***RHS_array;
  PetscReal      *RHS_local;
  PetscReal      *RHS1_local;
  Vec            RHS_localVec;
  Vec            perm_local;
  PetscReal      hx,hy,hz;  
  PetscReal      *K1_local,*K2_local,*KS_local,*K3_local,*K4_local,*KD_local,*KDF_local,*KF_local;
  PetscReal      mu;
  PetscReal      theta,timestepsize;
  PetscInt       nrow = ctx->e3D.nphix*ctx->e3D.nphiy*ctx->e3D.nphiz;
  MatStencil     *row;
  PetscReal      ***source_array;
  Vec            source_local;
  PetscReal      ****fluxbc_array;
	Vec            fluxbc_local;
  PetscReal      M_inv;
  FACE           face;
  PetscReal      hwx, hwy, hwz;
  PetscReal      ***v_array;
  Vec            v_local;
  PetscReal      ****u_diff_array;
	Vec            U_diff,u_diff_local;
  PetscReal      ****u_array;
	Vec            u_local;
  PetscReal      ****u_old_array;
	Vec            u_old_local;
  PetscReal      ***v_old_array;
	Vec            v_old_local;
  PetscReal      ***fracflow_array;
  Vec            fracflow_local;
  PetscReal      ***one_array;
  Vec            one_local;
  Vec            Ones;
  PetscReal      alphabiot;
  PetscReal      K_dr;
  PetscReal      ***pressure_diff_array;
	Vec            Pressure_diff,pressure_diff_local;
  PetscReal      ***pmult_array;
  Vec            pmult_local;


  PetscFunctionBegin;  
  alphabiot  = ctx->flowprop.alphabiot;
  K_dr  = ctx->flowprop.K_dr;
  M_inv     = ctx->flowprop.M_inv;
  theta = ctx->flowprop.theta;
  timestepsize = ctx->flowprop.timestepsize;
  mu     = ctx->flowprop.mu;
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = MatZeroEntries(K);CHKERRQ(ierr);
  ierr = MatZeroEntries(Krhs);CHKERRQ(ierr);
  ierr = VecSet(RHS,0.);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daScal,&RHS_localVec);CHKERRQ(ierr);
  ierr = VecSet(RHS_localVec,0.);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,RHS_localVec,&RHS_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ctx->daVect,&fluxbc_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,ctx->VelBCArray,INSERT_VALUES,fluxbc_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,ctx->VelBCArray,INSERT_VALUES,fluxbc_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,fluxbc_local,&fluxbc_array);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ctx->daScal,&source_local);CHKERRQ(ierr);  
  ierr = DMGlobalToLocalBegin(ctx->daScal,ctx->Source,INSERT_VALUES,source_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,ctx->Source,INSERT_VALUES,source_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,source_local,&source_array);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daVFperm,fields->vfperm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daVFperm,fields->vfperm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ctx->daScal,&v_old_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,ctx->V_old,INSERT_VALUES,v_old_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,ctx->V_old,INSERT_VALUES,v_old_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,v_old_local,&v_old_array);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(ctx->daScal,&Pressure_diff);CHKERRQ(ierr);
  ierr = VecSet(Pressure_diff,0.);CHKERRQ(ierr);
  ierr = VecAXPY(Pressure_diff,-1.0,ctx->pressure_old);CHKERRQ(ierr);
  ierr = VecAXPY(Pressure_diff,1.0,fields->pressure);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daScal,&pressure_diff_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,Pressure_diff,INSERT_VALUES,pressure_diff_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,Pressure_diff,INSERT_VALUES,pressure_diff_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,pressure_diff_local,&pressure_diff_array);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(ctx->daVect,&U_diff);CHKERRQ(ierr);
  ierr = VecSet(U_diff,0.);CHKERRQ(ierr);
  ierr = VecAXPY(U_diff,-1.0,ctx->U_old);CHKERRQ(ierr);
  ierr = VecAXPY(U_diff,1.0,ctx->U);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daVect,&u_diff_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,U_diff,INSERT_VALUES,u_diff_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,U_diff,INSERT_VALUES,u_diff_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,u_diff_local,&u_diff_array);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,ctx->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,ctx->U,INSERT_VALUES,u_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ctx->daVect,&u_old_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,ctx->U_old,INSERT_VALUES,u_old_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,ctx->U_old,INSERT_VALUES,u_old_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,u_old_local,&u_old_array);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ctx->daScal,&fracflow_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,ctx->RegFracWellFlowRate,INSERT_VALUES,fracflow_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,ctx->RegFracWellFlowRate,INSERT_VALUES,fracflow_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,fracflow_local,&fracflow_array);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(ctx->daScal,&Ones);CHKERRQ(ierr);
  ierr = VecSet(Ones,1.0);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daScal,&one_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,Ones,INSERT_VALUES,one_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,Ones,INSERT_VALUES,one_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,one_local,&one_array);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ctx->daScalCell,&pmult_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScalCell,fields->pmult,INSERT_VALUES,pmult_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScalCell,fields->pmult,INSERT_VALUES,pmult_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScalCell,pmult_local,&pmult_array);CHKERRQ(ierr);
  ierr = PetscMalloc7(nrow*nrow,&K1_local,
                      nrow*nrow,&K2_local,
                      nrow*nrow,&KS_local,
                      nrow*nrow,&KD_local,
                      nrow,&RHS_local,
                      4,&RHS1_local,
                      nrow,&row);CHKERRQ(ierr);
  ierr = PetscMalloc4(nrow*nrow,&K3_local,nrow*nrow,&K4_local,nrow*nrow,&KDF_local,nrow*nrow,&KF_local);CHKERRQ(ierr);
  
  for (ek = zs; ek < zs+zm; ek++) {
    for (ej = ys; ej < ys+ym; ej++) {
      for (ei = xs; ei < xs+xm; ei++) {
        hx   = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
        hy   = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
        hz   = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
        ierr = VFCartFEElement3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
        
        for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							row[l].i = ei+i;row[l].j = ej+j;row[l].k = ek+k;row[l].c = 0;
						}
					}
				}
        ierr = VF_MatA_local(KS_local,&ctx->e3D,ek,ej,ei,one_array);CHKERRQ(ierr);
        for (l = 0; l < nrow*nrow; l++) {
          KS_local[l] = M_inv*KS_local[l];
        }
        ierr = MatSetValuesStencil(K,nrow,row,nrow,row,KS_local,ADD_VALUES);CHKERRQ(ierr);
        ierr = MatSetValuesStencil(Krhs,nrow,row,nrow,row,KS_local,ADD_VALUES);CHKERRQ(ierr);
        if(ctx->FlowDisplCoupling && ctx->ResFlowMechCoupling == FIXEDSTRESS){
          ierr = VF_MatA_local(KF_local,&ctx->e3D,ek,ej,ei,v_array);CHKERRQ(ierr);
          for (l = 0; l < nrow*nrow; l++) {
            KF_local[l] = pow(alphabiot,2)*KF_local[l]/K_dr;
          }
          ierr = MatSetValuesStencil(K,nrow,row,nrow,row,KF_local,ADD_VALUES);CHKERRQ(ierr);
          ierr = MatSetValuesStencil(Krhs,nrow,row,nrow,row,KF_local,ADD_VALUES);CHKERRQ(ierr);
        }
        ierr = VF_HeatMatK_local(KD_local,&ctx->e3D,ek,ej,ei,perm_array,one_array);CHKERRQ(ierr);
        for (l = 0; l < nrow*nrow; l++) {
					K1_local[l] = theta/mu*timestepsize*KD_local[l];
          K2_local[l] = -1.*(1.-theta)/mu*timestepsize*KD_local[l];
				}
				ierr = MatSetValuesStencil(K,nrow,row,nrow,row,K1_local,ADD_VALUES);CHKERRQ(ierr);
        ierr = MatSetValuesStencil(Krhs,nrow,row,nrow,row,K2_local,ADD_VALUES);CHKERRQ(ierr);
        
        if(ctx->FractureFlowCoupling){
          ierr = VF_MatDFractureFlowCoupling_local(KD_local,&ctx->e3D,ek,ej,ei,u_array,v_array);CHKERRQ(ierr);
          for (l = 0; l < nrow*nrow; l++) {
            K1_local[l] = theta/(12.*mu)*timestepsize*KD_local[l];
            K2_local[l] = -1.*(1.-theta)/(12.*mu)*timestepsize*KD_local[l];
          }
          ierr = MatSetValuesStencil(K,nrow,row,nrow,row,K1_local,ADD_VALUES);CHKERRQ(ierr);
          ierr = MatSetValuesStencil(Krhs,nrow,row,nrow,row,K2_local,ADD_VALUES);CHKERRQ(ierr);
          ierr = VF_RHSFractureFlowCoupling_local(RHS_local,&ctx->e3D,ek,ej,ei,u_array,v_array,u_old_array,v_old_array);
            for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                  RHS_array[ek+k][ej+j][ei+i] += -1*RHS_local[l];
                }
              }
            }
          if(ctx->hasFlowWells){
            ierr = VecApplyFractureWellSource(RHS_local,fracflow_array,&ctx->e3D,ek,ej,ei,ctx,v_array);
//            ierr = VecApplySourceTerms(RHS_local,fracflow_array,&ctx->e3D,ek,ej,ei,ctx,one_array);
            for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                  RHS_array[ek+k][ej+j][ei+i] += -1.0*timestepsize*RHS_local[l];
                }
              }
            }
          }
        }
        ierr = Flow_Vecg(RHS_local,&ctx->e3D,ek,ej,ei,&ctx->flowprop,perm_array,v_array);CHKERRQ(ierr);
        for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
          for (j = 0; j < ctx->e3D.nphiy; j++) {
            for (i = 0; i < ctx->e3D.nphix; i++,l++) {
              RHS_array[ek+k][ej+j][ei+i] += -2.0*timestepsize*RHS_local[l];
            }
          }
        }
        if(ctx->hasFluidSources){
          ierr = VecApplySourceTerms(RHS_local,source_array,&ctx->e3D,ek,ej,ei,ctx,v_array);
          for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
            for (j = 0; j < ctx->e3D.nphiy; j++) {
              for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                RHS_array[ek+k][ej+j][ei+i] += -1.0*timestepsize*RHS_local[l];
              }
            }
          }
        }
        if(ctx->FlowDisplCoupling){
          ierr = VF_RHSFlowMechUCoupling_local(RHS_local,&ctx->e3D,ek,ej,ei,&ctx->flowprop,u_diff_array,v_array);
          for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
            for (j = 0; j < ctx->e3D.nphiy; j++) {
              for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                RHS_array[ek+k][ej+j][ei+i] += -1.0*RHS_local[l];
              }
            }
          }
        }
        if(ctx->FlowDisplCoupling && ctx->ResFlowMechCoupling == FIXEDSTRESS){
          ierr = VF_RHSFlowMechUCouplingFIXSTRESS_local(RHS_local,&ctx->e3D,ek,ej,ei,&ctx->flowprop,ctx->matprop,pressure_diff_array,v_array);
          for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
            for (j = 0; j < ctx->e3D.nphiy; j++) {
              for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                RHS_array[ek+k][ej+j][ei+i] += 1.0*RHS_local[l];
              }
            }
          }
        }
        if (ei == 0) {
/*            					 Face X0  */
					face = X0;
					ierr = VFCartFEElement2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
          if (ctx->bcQ[0].face[face] == FIXED) {
						ierr = VecApplyFluxBC(RHS1_local,fluxbc_array,ek,ej,ei,face,&ctx->e2D,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphix; k++){
							for (j = 0; j < ctx->e2D.nphiy; j++) {
								for (i = 0; i < ctx->e2D.nphiz; i++, l++) {
									RHS_array[ek+k][ej+j][ei+i] +=  -timestepsize*RHS1_local[l];
								}
							}
						}
					}
				}
				if (ei == nx-1) {
/*            					 Face X1  */
					face = X1;
					ierr = VFCartFEElement2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
					if (ctx->bcQ[0].face[face] == FIXED) {
						ierr = VecApplyFluxBC(RHS1_local,fluxbc_array,ek,ej,ei,face,&ctx->e2D,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphix; k++){
							for (j = 0; j < ctx->e2D.nphiy; j++) {
								for (i = 0; i < ctx->e2D.nphiz; i++, l++) {
									RHS_array[ek+k][ej+j][ei+1] +=  -timestepsize*RHS1_local[l];
								}
							}
						}
					}
				}
				if (ej == 0) {
/*            					 Face Y0  */
					face = Y0;
					ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
					if (ctx->bcQ[1].face[face] == FIXED) {
						ierr = VecApplyFluxBC(RHS1_local,fluxbc_array,ek,ej,ei,face,&ctx->e2D,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphiy; k++){
							for (j = 0; j < ctx->e2D.nphiz; j++) {
								for (i = 0; i < ctx->e2D.nphix; i++, l++) {
									RHS_array[ek+k][ej+j][ei+i] +=  -timestepsize*RHS1_local[l];
								}
							}
						}
					}
				}
				if (ej == ny-1) {
/*            					 Face Y1  */
					face = Y1;
					ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
					if (ctx->bcQ[1].face[face] == FIXED) {
						ierr = VecApplyFluxBC(RHS1_local,fluxbc_array,ek,ej,ei,face,&ctx->e2D,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphiy; k++){
							for (j = 0; j < ctx->e2D.nphiz; j++) {
								for (i = 0; i < ctx->e2D.nphix; i++, l++) {
									RHS_array[ek+k][ej+1][ei+i] +=  -timestepsize*RHS1_local[l];
								}
							}
						}
					}
				}
				if (ek == 0) {
/*            					 Face Z0  */
					face = Z0;
					ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
					if (ctx->bcQ[2].face[face] == FIXED) {
						ierr = VecApplyFluxBC(RHS1_local,fluxbc_array,ek,ej,ei,face,&ctx->e2D,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphiz; k++){
							for (j = 0; j < ctx->e2D.nphiy; j++) {
								for (i = 0; i < ctx->e2D.nphix; i++, l++) {
									RHS_array[ek+k][ej+j][ei+i] +=  -timestepsize*RHS1_local[l];
								}
							}
						}
					}
				}
				if (ek == nz-1) {
/*            					 Face Z1  */
					face = Z1;
					ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
					if (ctx->bcQ[2].face[face] == FIXED) {
						ierr = VecApplyFluxBC(RHS1_local,fluxbc_array,ek,ej,ei,face,&ctx->e2D,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphiz; k++){
							for (j = 0; j < ctx->e2D.nphiy; j++) {
								for (i = 0; i < ctx->e2D.nphix; i++, l++) {
									RHS_array[ek+1][ej+j][ei+i] +=  -timestepsize*RHS1_local[l];
								}
							}
						}
					}
				}
      }
    }
  }
  
  PetscInt  w_no = 0;
  PetscInt  w_no1 = 0;
  if(ctx->hasFlowWells){
    while(w_no < ctx->numWells){
      for (ek = zs; ek < zs+zm; ek++) {
        for (ej = ys; ej < ys+ym; ej++) {
          for (ei = xs; ei < xs+xm; ei++) {
            if(
               ((coords_array[ek][ej][ei+1][0] >= ctx->well[w_no].coords[0]) && (coords_array[ek][ej][ei][0] <= ctx->well[w_no].coords[0] ))    &&
               ((coords_array[ek][ej+1][ei][1] >= ctx->well[w_no].coords[1]) && (coords_array[ek][ej][ei][1] <= ctx->well[w_no].coords[1] ))    &&
               ((coords_array[ek+1][ej][ei][2] >= ctx->well[w_no].coords[2]) && (coords_array[ek][ej][ei][2] <= ctx->well[w_no].coords[2] ))
               )
            {
              hx   = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
              hy   = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
              hz   = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
              ierr = VFCartFEElement3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
              hwx = (ctx->well[w_no].coords[0]-coords_array[ek][ej][ei][0])/hx;
              hwy = (ctx->well[w_no].coords[1]-coords_array[ek][ej][ei][1])/hy;
              hwz = (ctx->well[w_no].coords[2]-coords_array[ek][ej][ei][2])/hz;
              if(ctx->well[w_no].condition == RATE){
                ierr = VecApplyWellFlowRate(RHS_local,&ctx->e3D,ctx->well[w_no].Qw,hwx,hwy,hwz,ek,ej,ei,one_array);CHKERRQ(ierr);
                for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
                  for (j = 0; j < ctx->e3D.nphiy; j++) {
                    for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                      if(ctx->well[w_no].type == INJECTOR){
                        RHS_array[ek+k][ej+j][ei+i] += timestepsize*RHS_local[l];
                      }
                      else if(ctx->well[w_no].type == PRODUCER)
                      {
                        RHS_array[ek+k][ej+j][ei+i] += -timestepsize*RHS_local[l];
                      }
                    }
                  }
                }
              }
              else if(ctx->well[w_no].condition == PRESSURE){
              }
              w_no1++;
            }
          }
        }
      }
      ierr = MPI_Allreduce(&w_no1,&w_no,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);

    }
  }

  ierr = MatAssemblyBegin(Krhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Krhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatApplyDirichletBC(K,&ctx->bcP[0]);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Krhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Krhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_old_local,&u_old_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daVect,&u_old_local);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_diff_local,&u_diff_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daVect,&u_diff_local);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScal,pressure_diff_local,&pressure_diff_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&pressure_diff_local);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(ctx->daScal,v_old_local,&v_old_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&v_old_local);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScal,source_local,&source_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&source_local);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,fluxbc_local,&fluxbc_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVect,&fluxbc_local);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScal,RHS_localVec,&RHS_array);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(ctx->daScal,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ctx->daScal,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&RHS_localVec);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(ctx->daScal,fracflow_local,&fracflow_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&fracflow_local);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(ctx->daScal,one_local,&one_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&one_local);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(ctx->daScalCell,pmult_local,&pmult_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScalCell,&pmult_local);CHKERRQ(ierr);
  
  ierr = PetscFree7(K1_local,K2_local,KS_local,KD_local,RHS_local,RHS1_local,row);CHKERRQ(ierr);
  ierr = PetscFree4(K3_local,K4_local,KDF_local,KF_local);CHKERRQ(ierr);
  
  ierr = VecDestroy(&U_diff);CHKERRQ(ierr);
  ierr = VecDestroy(&Pressure_diff);CHKERRQ(ierr);
  ierr = VecDestroy(&Ones);CHKERRQ(ierr);
  
 /*
  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"Matrixf.txt",&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
  ierr = MatView(K,viewer);CHKERRQ(ierr);
  
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"Matrixlhsf.txt",&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
  ierr = MatView(Krhs,viewer);CHKERRQ(ierr);
  
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"RHSf.txt",&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
  ierr = VecView(RHS,viewer);CHKERRQ(ierr);
  Vec rhsc;
  Vec rhsc1;
  ierr = VecDuplicate(RHS,&rhsc);CHKERRQ(ierr);
  ierr = VecDuplicate(RHS,&rhsc1);CHKERRQ(ierr);
  ierr = MatMult(Krhs,ctx->pressure_old,rhsc);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"RHS.txt",&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
  ierr = VecView(rhsc,viewer);CHKERRQ(ierr);
  

  ierr = MatMultAdd(Krhs,ctx->pressure_old ,RHS,rhsc1);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"RHS1.txt",&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
  ierr = VecView(rhsc1,viewer);CHKERRQ(ierr);
 */
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_FormFlowStandardFEMIJacobian"
extern PetscErrorCode VF_FormFlowStandardFEMIJacobian(SNES snes,Vec Pressure,Mat Jac,Mat JacPre,void *user)
{
	PetscErrorCode    ierr;
	VFCtx             *ctx=(VFCtx*)user;
	
	PetscFunctionBegin;
	ierr = MatZeroEntries(JacPre);CHKERRQ(ierr);
	ierr = MatCopy(ctx->KP,JacPre,SAME_NONZERO_PATTERN);
	ierr = MatAssemblyBegin(JacPre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(JacPre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	if (Jac != JacPre) {
		ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}











