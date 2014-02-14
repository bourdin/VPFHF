/*
 VFHeat_SNESFEM.c
 
 (c) 2011-2013 C. Chukwudozie, LSU
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFHeat.h"
#include "VFWell.h"
#include "VFHeat_SNESFEM.h"
#include "VFFlow_KSPMixedFEM.h"
#include "VFFlow_SNESMixedFEM.h"
#include "VFFlow.h"

/*
 VFHeat_SNESFEM
 */

/*
 ################################################################################################################
 SNES ROUTINE
 ################################################################################################################
 */
#undef __FUNCT__
#define __FUNCT__ "VF_FEMSNESHeatSolverInitialize"
extern PetscErrorCode VF_FEMSNESHeatSolverInitialize(VFCtx *ctx, VFFields *fields)
{
	PetscMPIInt    comm_size;
	PetscErrorCode ierr;
  PetscFunctionBegin;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&comm_size);CHKERRQ(ierr);
  ierr = DMCreateMatrix(ctx->daScal,MATAIJ,&ctx->KT);CHKERRQ(ierr);
  ierr = DMCreateMatrix(ctx->daScal,MATAIJ,&ctx->KTlhs);CHKERRQ(ierr);
  ierr = DMCreateMatrix(ctx->daScal,MATAIJ,&ctx->JacT);CHKERRQ(ierr);
	ierr = MatZeroEntries(ctx->KT);CHKERRQ(ierr);
	ierr = MatZeroEntries(ctx->KTlhs);CHKERRQ(ierr);
	ierr = MatZeroEntries(ctx->JacT);CHKERRQ(ierr);
	ierr = MatSetOption(ctx->KT,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
	ierr = MatSetOption(ctx->KTlhs,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
	ierr = MatSetOption(ctx->JacT,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daScal,&ctx->prevT);CHKERRQ(ierr);
	ierr = VecSet(ctx->prevT,0.0);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daScal,&ctx->RHST);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daScal,&ctx->HeatFunct);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)ctx->RHST,"RHS vector of heat equation");CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)ctx->HeatFunct,"RHS of SNES heat solver");CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daScal,&ctx->RHSTpre);CHKERRQ(ierr);
	ierr = VecSet(ctx->RHSTpre,0.);CHKERRQ(ierr);
  
	ierr = SNESCreate(PETSC_COMM_WORLD,&ctx->snesT);CHKERRQ(ierr);
  ierr = SNESAppendOptionsPrefix(ctx->snesT,"HeatSnes_");CHKERRQ(ierr);
  ierr = SNESSetFromOptions(ctx->snesT);CHKERRQ(ierr);
	ierr = BCTInit(&ctx->bcT[0],ctx);
	ierr = BCQTInit(&ctx->bcQT[0],ctx);	
		/*	ierr = GetFlowProp(&ctx->flowprop,ctx->units,ctx->resprop);CHKERRQ(ierr); */
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_FEMSNESHeatSolverFinalize"
extern PetscErrorCode VF_FEMSNESHeatSolverFinalize(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	ierr = MatDestroy(&ctx->KT);CHKERRQ(ierr);
	ierr = MatDestroy(&ctx->KTlhs);CHKERRQ(ierr);
	ierr = MatDestroy(&ctx->JacT);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->prevT);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->RHST);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->HeatFunct);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->RHSTpre);CHKERRQ(ierr);
	ierr = SNESDestroy(&ctx->snesT);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_HeatFEMSNESSolve"
extern PetscErrorCode VF_HeatFEMSNESSolve(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode      ierr;
	SNESConvergedReason reason;	
	PetscInt            its;
	PetscReal           Tmin,Tmax;
	
	PetscFunctionBegin;	
	ierr = FormHeatMatricesnVector(ctx->KT,ctx->KTlhs,ctx->RHST,ctx,fields);CHKERRQ(ierr);
	ierr = SNESSetFunction(ctx->snesT,ctx->HeatFunct,FormSNESHeatIFunction,ctx);CHKERRQ(ierr);
    ierr = SNESSetJacobian(ctx->snesT,ctx->JacT,ctx->JacT,FormSNESHeatIJacobian,ctx);CHKERRQ(ierr);
	if (ctx->verbose > 1) {
		ierr = SNESMonitorSet(ctx->snesT,FEMSNESMonitor,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	}	
    ierr = SNESSolve(ctx->snesT,PETSC_NULL,fields->theta);CHKERRQ(ierr);
	ierr = SNESGetConvergedReason(ctx->snesT,&reason);CHKERRQ(ierr);
	if (reason < 0) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"[ERROR] snesT diverged with reason %d\n",(int)reason);CHKERRQ(ierr);
	} else {
		ierr = SNESGetIterationNumber(ctx->snesT,&its);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"      snesT converged in %d iterations %d.\n",(int)its,(int)reason);CHKERRQ(ierr);
	}
	ierr = VecMin(fields->theta,PETSC_NULL,&Tmin);CHKERRQ(ierr);
	ierr = VecMax(fields->theta,PETSC_NULL,&Tmax);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"      T min / max:     %e %e\n",Tmin,Tmax);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormSNESHeatIJacobian"
extern PetscErrorCode FormSNESHeatIJacobian(SNES snes,Vec T,Mat *Jac,Mat *Jacpre,MatStructure *str,void *user)
{
	PetscErrorCode  ierr;
	VFCtx				   *ctx=(VFCtx*)user;
	
	PetscFunctionBegin;
	*str = DIFFERENT_NONZERO_PATTERN;
	ierr = MatZeroEntries(*Jac);CHKERRQ(ierr);
	ierr = MatCopy(ctx->KT,*Jac,*str);
	ierr = MatAssemblyBegin(*Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	if (*Jac != *Jacpre) {
		ierr = MatAssemblyBegin(*Jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(*Jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormSNESHeatIFunction"
extern PetscErrorCode FormSNESHeatIFunction(SNES snes,Vec T,Vec Func,void *user)
{
	
	PetscErrorCode ierr;
	VFCtx			*ctx=(VFCtx*)user;
	PetscReal		theta,one_minus_theta;
	Vec             VecRHS;
	
	PetscFunctionBegin;
	theta = ctx->flowprop.theta;
	one_minus_theta = (1.-theta);
	
	ierr = VecDuplicate(ctx->RHST,&VecRHS);CHKERRQ(ierr);
	ierr = VecCopy(ctx->RHST,VecRHS);CHKERRQ(ierr);
	ierr = VecSet(Func,0.0);CHKERRQ(ierr);
	ierr = VecAXPBY(VecRHS,one_minus_theta,theta,ctx->RHSTpre);CHKERRQ(ierr);	
	ierr = MatMultAdd(ctx->KTlhs,ctx->prevT,VecRHS,VecRHS);CHKERRQ(ierr);	
	ierr = VecApplyDirichletBC(VecRHS,ctx->TBCArray,&ctx->bcT[0]);CHKERRQ(ierr);
	ierr = MatMult(ctx->KT,T,Func);CHKERRQ(ierr);	
	ierr = VecAXPY(Func,-1.0,VecRHS);CHKERRQ(ierr);
	ierr = VecDestroy(&VecRHS);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormHeatMatricesnVector"
extern PetscErrorCode FormHeatMatricesnVector(Mat K,Mat Klhs,Vec RHS,VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode ierr;
	PetscInt       xs,xm,nx;
	PetscInt       ys,ym,ny;
	PetscInt       zs,zm,nz;
	PetscInt       ek,ej,ei;
	PetscInt       i,j,k,l;
	PetscReal      ****coords_array;
	PetscReal      ***RHS_array;
	PetscReal      *RHS_local;
	PetscReal      *RHS1_local;
	Vec            RHS_localVec;
	PetscReal      hx,hy,hz;
	PetscReal      theta,timestepsize;
	PetscInt       nrow = ctx->e3D.nphix*ctx->e3D.nphiy*ctx->e3D.nphiz;
	MatStencil     *row;
	FACE           face;
	MatStructure   flg;
	PetscReal      *KM_local,*KC_local,*KD_local,*K1_local,*K2_local,*KN_local,*KS_local;
	Vec            diffsvty_local;
	PetscReal      ****diffsvty_array;
	PetscReal      ***heatsource_array;
	Vec            heatsource_local;
	Vec            vel_local;
	PetscReal      ****vel_array;
	PetscReal      ****fluxbc_array;
	Vec            fluxbc_local;
	/*	Vec				rhoCp_eff_local,rho_liq_local,Cp_liq_local,Cp_sol_local;	*/
	PetscReal      rhoCp_eff_array;	
	PetscReal      rho_liq_array;
	PetscReal      Cp_liq_array;
	PetscReal      rho_sol_array;
	PetscReal      Cp_sol_array;
  PetscReal      ***v_array;
	Vec            v_local;
	
	PetscFunctionBegin;
	flg = SAME_NONZERO_PATTERN;
	timestepsize = ctx->flowprop.timestepsize;	
	theta = ctx->flowprop.theta;
	rho_liq_array = ctx->flowprop.rho;
	Cp_liq_array = ctx->flowprop.Cp;
	rho_sol_array = ctx->matprop[0].rho;
	Cp_sol_array = ctx->matprop[0].Cp;
	ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = MatZeroEntries(K);CHKERRQ(ierr);
	ierr = MatZeroEntries(Klhs);CHKERRQ(ierr);
	ierr = VecSet(RHS,0.);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScal,&RHS_localVec);CHKERRQ(ierr);
	ierr = VecSet(RHS_localVec,0.);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,RHS_localVec,&RHS_array);CHKERRQ(ierr);	
	
	/*
	 ierr = DMGetLocalVector(ctx->daScal,&rhoCp_eff_local);CHKERRQ(ierr);
	 ierr = DMGlobalToLocalBegin(ctx->daScal,,INSERT_VALUES,rhoCp_eff_local);CHKERRQ(ierr);
	 ierr = DMGlobalToLocalEnd(ctx->daScal,,INSERT_VALUES,rhoCp_eff_local);CHKERRQ(ierr);
	 ierr = DMDAVecGetArray(ctx->daScal,rhoCp_eff_local,&rhoCp_eff_array);CHKERRQ(ierr);
	 
	 ierr = DMGetLocalVector(ctx->daScal,&rho_liq_local);CHKERRQ(ierr);
	 ierr = DMGlobalToLocalBegin(ctx->daScal,,INSERT_VALUES,rho_liq_local);CHKERRQ(ierr);
	 ierr = DMGlobalToLocalEnd(ctx->daScal,,INSERT_VALUES,rho_liq_local);CHKERRQ(ierr);
	 ierr = DMDAVecGetArray(ctx->daScal,rho_liq_local,&rho_liq_array);CHKERRQ(ierr);
	 
	 ierr = DMGetLocalVector(ctx->daScal,&Cp_liq_local);CHKERRQ(ierr);
	 ierr = DMGlobalToLocalBegin(ctx->daScal,,INSERT_VALUES,Cp_liq_local);CHKERRQ(ierr);
	 ierr = DMGlobalToLocalEnd(ctx->daScal,,INSERT_VALUES,Cp_liq_local);CHKERRQ(ierr);
	 ierr = DMDAVecGetArray(ctx->daScal,Cp_liq_local,&Cp_liq_array);CHKERRQ(ierr);
	 
	 ierr = DMGetLocalVector(ctx->daScal,&rho_sol_local);CHKERRQ(ierr);
	 ierr = DMGlobalToLocalBegin(ctx->daScal,,INSERT_VALUES,rho_sol_local);CHKERRQ(
	 ierr = DMGlobalToLocalEnd(ctx->daScal,,INSERT_VALUES,rho_sol_local);CHKERRQ(ierr);
	 ierr = DMDAVecGetArray(ctx->daScal,rho_sol_local,&rho_sol_array);CHKERRQ(ierr);
	 
	 ierr = DMGetLocalVector(ctx->daScal,&Cp_sol_local);CHKERRQ(ierr);
	 ierr = DMGlobalToLocalBegin(ctx->daScal,,INSERT_VALUES,Cp_sol_local);CHKERRQ(ierr);
	 ierr = DMGlobalToLocalEnd(ctx->daScal,,INSERT_VALUES,Cp_sol_local);CHKERRQ(ierr);
	 ierr = DMDAVecGetArray(ctx->daScal,Cp_sol_local,&Cp_sol_array);CHKERRQ(ierr);
	 */
	ierr = DMGetLocalVector(ctx->daVect,&fluxbc_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,ctx->HeatFluxBCArray,INSERT_VALUES,fluxbc_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,ctx->HeatFluxBCArray,INSERT_VALUES,fluxbc_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,fluxbc_local,&fluxbc_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScal,&heatsource_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,ctx->HeatSource,INSERT_VALUES,heatsource_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,ctx->HeatSource,INSERT_VALUES,heatsource_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,heatsource_local,&heatsource_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daVFperm,&diffsvty_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVFperm,ctx->Cond,INSERT_VALUES,diffsvty_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVFperm,ctx->Cond,INSERT_VALUES,diffsvty_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVFperm,diffsvty_local,&diffsvty_array);CHKERRQ(ierr);	
	ierr = DMGetLocalVector(ctx->daVect,&vel_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,fields->velocity,INSERT_VALUES,vel_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,fields->velocity,INSERT_VALUES,vel_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,vel_local,&vel_array);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = PetscMalloc5(nrow*nrow,PetscReal,&KM_local,
						nrow*nrow,PetscReal,&KC_local,
						nrow*nrow,PetscReal,&KD_local,
						nrow*nrow,PetscReal,&K1_local,
						nrow*nrow,PetscReal,&K2_local);CHKERRQ(ierr);
	ierr = PetscMalloc5(nrow*nrow,PetscReal,&KN_local,
						nrow*nrow,PetscReal,&KS_local,
						nrow,PetscReal,&RHS_local,
						nrow,PetscReal,&RHS1_local,
						nrow,MatStencil,&row);CHKERRQ(ierr);	
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx   = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy   = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz   = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				/*	Assembling the sub-Matrices	*/
				rhoCp_eff_array=rho_liq_array*Cp_liq_array+rho_sol_array*Cp_sol_array;
				ierr = VF_MatA_local(KM_local,&ctx->e3D,ek,ej,ei,v_array);CHKERRQ(ierr);
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							row[l].i = ei+i;row[l].j = ej+j;row[l].k = ek+k;row[l].c = 0;
						}
					}
				}
				ierr = MatSetValuesStencil(K,nrow,row,nrow,row,KM_local,ADD_VALUES);CHKERRQ(ierr);
				ierr = MatSetValuesStencil(Klhs,nrow,row,nrow,row,KM_local,ADD_VALUES);CHKERRQ(ierr);
				ierr = VF_HeatMatK_local(KD_local,&ctx->e3D,ek,ej,ei,diffsvty_array,v_array);CHKERRQ(ierr);

				for (l = 0; l < nrow*nrow; l++) {
					K1_local[l] = theta*timestepsize*(1./rhoCp_eff_array)*KD_local[l];
					K2_local[l] = -1.*(1.-theta)*timestepsize*(1./rhoCp_eff_array)*KD_local[l];
				}
				ierr = MatSetValuesStencil(K,nrow,row,nrow,row,K1_local,ADD_VALUES);CHKERRQ(ierr);
				ierr = MatSetValuesStencil(Klhs,nrow,row,nrow,row,K2_local,ADD_VALUES);CHKERRQ(ierr);
				ierr = VF_HeatMatC_local(KC_local,&ctx->e3D,ek,ej,ei,vel_array,v_array);CHKERRQ(ierr);

				for (l = 0; l < nrow*nrow; l++) {
					K1_local[l] = theta*timestepsize*(rho_liq_array*Cp_liq_array/rhoCp_eff_array)*KC_local[l];
/*					K1_local[l] = timestepsize*KC_local[l];   */
					K2_local[l] = -1.0*(1.-theta)*timestepsize*(rho_liq_array*Cp_liq_array/rhoCp_eff_array)*KC_local[l];
				}
				ierr = MatSetValuesStencil(K,nrow,row,nrow,row,K1_local,ADD_VALUES);CHKERRQ(ierr);
				ierr = MatSetValuesStencil(Klhs,nrow,row,nrow,row,K2_local,ADD_VALUES);CHKERRQ(ierr);
				ierr = VF_HeatMatN_local(KN_local,&ctx->e3D,ek,ej,ei,vel_array,v_array);CHKERRQ(ierr);

				for (l = 0; l < nrow*nrow; l++) {		
					K1_local[l] = 1./2.*theta*timestepsize*timestepsize*(rho_liq_array*Cp_liq_array/rhoCp_eff_array)*(rho_liq_array*Cp_liq_array/rhoCp_eff_array)*KN_local[l];
					K2_local[l] = -1./2.*(1-theta)*timestepsize*timestepsize*(rho_liq_array*Cp_liq_array/rhoCp_eff_array)*(rho_liq_array*Cp_liq_array/rhoCp_eff_array)*KN_local[l];
          
        }

/*        ierr = MatSetValuesStencil(K,nrow,row,nrow,row,K1_local,ADD_VALUES);CHKERRQ(ierr);
				ierr = MatSetValuesStencil(Klhs,nrow,row,nrow,row,K2_local,ADD_VALUES);CHKERRQ(ierr);   */
        /*	Assembling the righthand side vector	*/
				if(ctx->hasHeatSources){
					ierr = VecApplyHeatSourceTerms(RHS_local,RHS1_local,heatsource_array,&ctx->e3D,ek,ej,ei,ctx,vel_array,v_array);
					for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
						for (j = 0; j < ctx->e3D.nphiy; j++) {
							for (i = 0; i < ctx->e3D.nphix; i++,l++) {
								RHS_array[ek+k][ej+j][ei+i] += timestepsize*(1./rhoCp_eff_array)*RHS_local[l];
							}
						}
					}
				}
				/*	Assembling contribution from flux boundary terms	*/
				if (ei == 0) {
					/*					 Face X0			*/
					face = X0;	
					ierr = CartFE_Element2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
					if (ctx->bcQT[0].face[face] == FIXED) {
						ierr = VecApplyFluxBC(RHS_local,fluxbc_array,ek,ej,ei,face,&ctx->e2D,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphix; k++){
							for (j = 0; j < ctx->e2D.nphiy; j++) {
								for (i = 0; i < ctx->e2D.nphiz; i++, l++) {
									RHS_array[ek+k][ej+j][ei+i] +=  -timestepsize*(1./rhoCp_eff_array)*RHS_local[l];
								}
							}
						}
					}
				}
				if (ei == nx-1) {
					/*					 Face X1		*/
					face = X1;
					ierr = CartFE_Element2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
					if (ctx->bcQT[0].face[face] == FIXED) {
						ierr = VecApplyFluxBC(RHS_local,fluxbc_array,ek,ej,ei,face,&ctx->e2D,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphix; k++){
							for (j = 0; j < ctx->e2D.nphiy; j++) {
								for (i = 0; i < ctx->e2D.nphiz; i++, l++) {
									RHS_array[ek+k][ej+j][ei+1] +=  -timestepsize*(1./rhoCp_eff_array)*RHS_local[l];
								}
							}
						}
					}
				}				
				if (ej == 0) {
					/*					 Face Y0		*/
					face = Y0;
					ierr = CartFE_Element2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
					if (ctx->bcQT[1].face[face] == FIXED) {
						ierr = VecApplyFluxBC(RHS_local,fluxbc_array,ek,ej,ei,face,&ctx->e2D,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphiy; k++){
							for (j = 0; j < ctx->e2D.nphiz; j++) {
								for (i = 0; i < ctx->e2D.nphix; i++, l++) {
									RHS_array[ek+k][ej+j][ei+i] +=  -timestepsize*(1./rhoCp_eff_array)*RHS_local[l];
								}
							}
						}
					}
				}
				if (ej == ny-1) {
					/*					 Face Y1		*/
					face = Y1;
					ierr = CartFE_Element2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
					if (ctx->bcQT[1].face[face] == FIXED) {
						ierr = VecApplyFluxBC(RHS_local,fluxbc_array,ek,ej,ei,face,&ctx->e2D,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphiy; k++){
							for (j = 0; j < ctx->e2D.nphiz; j++) {
								for (i = 0; i < ctx->e2D.nphix; i++, l++) {
									RHS_array[ek+k][ej+1][ei+i] +=  -timestepsize*(1./rhoCp_eff_array)*RHS_local[l];
								}
							}
						}
					}
				}
				if (ek == 0) {
					/*					 Face Z0		*/
					face = Z0;
					ierr = CartFE_Element2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
					if (ctx->bcQT[2].face[face] == FIXED) {
						ierr = VecApplyFluxBC(RHS_local,fluxbc_array,ek,ej,ei,face,&ctx->e2D,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphiz; k++){
							for (j = 0; j < ctx->e2D.nphiy; j++) {
								for (i = 0; i < ctx->e2D.nphix; i++, l++) {
									RHS_array[ek+k][ej+j][ei+i] +=  -timestepsize*(1./rhoCp_eff_array)*RHS_local[l];
								}
							}
						}
					}
				}
				if (ek == nz-1) {
					/*					 Face Z1		*/
					face = Z1;
					ierr = CartFE_Element2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
					if (ctx->bcQT[2].face[face] == FIXED) {
						ierr = VecApplyFluxBC(RHS_local,fluxbc_array,ek,ej,ei,face,&ctx->e2D,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphiz; k++){
							for (j = 0; j < ctx->e2D.nphiy; j++) {
								for (i = 0; i < ctx->e2D.nphix; i++, l++) {
									RHS_array[ek+1][ej+j][ei+i] +=  -timestepsize*(1./rhoCp_eff_array)*RHS_local[l];
								}
							}
						}
					}
				}
				
			}
		}
	}
	ierr = MatAssemblyBegin(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);	
	ierr = MatApplyDirichletBC(K,ctx->daScal,&ctx->bcT[0]);CHKERRQ(ierr);	
	ierr = MatAssemblyBegin(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,heatsource_local,&heatsource_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&heatsource_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,diffsvty_local,&diffsvty_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVFperm,&diffsvty_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,vel_local,&vel_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVect,&vel_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,RHS_localVec,&RHS_array);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScal,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScal,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&RHS_localVec);CHKERRQ(ierr);	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	/*		
	 ierr = DMDAVecRestoreArray(ctx->daScal,rhoCp_eff_local,&rhoCp_eff_array);CHKERRQ(ierr);
	 ierr = DMRestoreLocalVector(ctx->daScal,&rhoCp_eff_local);CHKERRQ(ierr);
	 ierr = DMDAVecRestoreArray(ctx->daScal,rho_liq_local,&rho_liq_array);CHKERRQ(ierr);
	 ierr = DMRestoreLocalVector(ctx->daScal,&rho_liq_local);CHKERRQ(ierr);
	 ierr = DMDAVecRestoreArray(ctx->daScal,Cp_liq_local,&Cp_liq_array);CHKERRQ(ierr);
	 ierr = DMRestoreLocalVector(ctx->daScal,&Cp_liq_local);CHKERRQ(ierr);
	 ierr = DMDAVecRestoreArray(ctx->daScal,rho_sol_local,&rho_sol_array);CHKERRQ(ierr);
	 ierr = DMRestoreLocalVector(ctx->daScal,&rho_sol_local);CHKERRQ(ierr);
	 ierr = DMDAVecRestoreArray(ctx->daScal,Cp_sol_local,&Cp_sol_array);CHKERRQ(ierr);
	 ierr = DMRestoreLocalVector(ctx->daScal,&Cp_sol_local);CHKERRQ(ierr);
	 */
  

   PetscViewer viewer;
   ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"Matrix.txt",&viewer);CHKERRQ(ierr);
   ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
   ierr = MatView(K,viewer);CHKERRQ(ierr);
   
   ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"Matrixlhs.txt",&viewer);CHKERRQ(ierr);
   ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
   ierr = MatView(Klhs,viewer);CHKERRQ(ierr);
   
   ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"RHS.txt",&viewer);CHKERRQ(ierr);
   ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
   ierr = VecView(RHS,viewer);CHKERRQ(ierr);
 
	ierr = DMDAVecRestoreArray(ctx->daVect,fluxbc_local,&fluxbc_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVect,&fluxbc_local);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = PetscFree5(KM_local,KC_local,KD_local,K1_local,K2_local);CHKERRQ(ierr);
	ierr = PetscFree5(KN_local,KS_local,RHS_local,RHS1_local,row);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "VecApplyFluxBC"
extern PetscErrorCode VecApplyFluxBC(PetscReal *RHS_local,PetscReal ****flux_array,PetscInt ek,PetscInt ej,PetscInt ei,FACE face,CartFE_Element2D *e,PetscReal ***v_array)
{
	PetscErrorCode ierr;
	PetscInt       i,j,k,l,g;
	PetscReal      *flux_elem,*v_elem;
	
	PetscFunctionBegin;
  ierr = PetscMalloc2(e->ng,PetscReal,&flux_elem,e->ng,PetscReal,&v_elem);CHKERRQ(ierr);
	for (g = 0; g < e->ng; g++) {
    v_elem[g] = 0;
    flux_elem[g] = 0;
	}
	switch (face) {
		case X0:
			for (k = 0; k < e->nphix; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphiz; i++) {
						for (g = 0; g < e->ng; g++) {
							flux_elem[g] += e->phi[i][j][k][g]*flux_array[ek+k][ej+j][ei+i][0];
							v_elem[g] += e->phi[i][j][k][g]*v_array[ek+k][ej+j][ei+i];
						}
					}
				}
			}
			/*			Accumulate		*/
			for (l=0,k = 0; k < e->nphix; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphiz; i++, l++) {
						RHS_local[l] = 0.;
						for (g = 0; g < e->ng; g++) {
							RHS_local[l] += -1*e->weight[g]*e->phi[i][j][k][g]*(pow(v_elem[g],2))*flux_elem[g];
						}
					}
				}
			}
			break;
		case X1:
			for (k = 0; k < e->nphix; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphiz; i++) {
						for (g = 0; g < e->ng; g++) {
							flux_elem[g] += e->phi[i][j][k][g]*flux_array[ek+k][ej+j][ei+1][0];
							v_elem[g] += e->phi[i][j][k][g]*v_array[ek+k][ej+j][ei+1];
						}
					}
				}
			}
			/*			Accumulate		*/
			for (l=0,k = 0; k < e->nphix; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphiz; i++, l++) {
						RHS_local[l] = 0.;
						for (g = 0; g < e->ng; g++) {
							RHS_local[l] += e->weight[g]*e->phi[i][j][k][g]*(pow(v_elem[g],2))*flux_elem[g];
						}
					}
				}
			}
			break;
		case Y0:
			for (k = 0; k < e->nphiy; k++) {
				for (j = 0; j < e->nphiz; j++) {
					for (i = 0; i < e->nphix; i++) {
						for (g = 0; g < e->ng; g++) {
							flux_elem[g] += e->phi[j][k][i][g]*flux_array[ek+k][ej+j][ei+i][1];
							v_elem[g] += e->phi[j][k][i][g]*v_array[ek+k][ej+j][ei+i];
						}
					}
				}
			}
			/*			Accumulate		*/
			for (l=0,k = 0; k < e->nphiy; k++) {
				for (j = 0; j < e->nphiz; j++) {
					for (i = 0; i < e->nphix; i++, l++) {
						RHS_local[l] = 0.;
						for (g = 0; g < e->ng; g++) {
							RHS_local[l] += -1*e->weight[g]*e->phi[j][k][i][g]*(pow(v_elem[g],2))*flux_elem[g];
						}
					}
				}
			}
			break;
		case Y1:
			for (k = 0; k < e->nphiy; k++) {
				for (j = 0; j < e->nphiz; j++) {
					for (i = 0; i < e->nphix; i++) {
						for (g = 0; g < e->ng; g++) {
							flux_elem[g] += e->phi[j][k][i][g]*flux_array[ek+k][ej+1][ei+i][1];
							v_elem[g] += e->phi[j][k][i][g]*v_array[ek+k][ej+1][ei+i];
						}
					}
				}
			}
			/*			Accumulate		*/
			for (l=0,k = 0; k < e->nphiy; k++) {
				for (j = 0; j < e->nphiz; j++) {
					for (i = 0; i < e->nphix; i++, l++) {
						RHS_local[l] = 0.;
						for (g = 0; g < e->ng; g++) {
							RHS_local[l] += e->weight[g]*e->phi[j][k][i][g]*(pow(v_elem[g],2))*flux_elem[g];
						}
					}
				}
			}
			break;
		case Z0:
			for (k = 0; k < e->nphiz; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphix; i++) {
						for (g = 0; g < e->ng; g++) {
							flux_elem[g] += e->phi[k][j][i][g]*flux_array[ek][ej+j][ei+i][2];
							v_elem[g] += e->phi[k][j][i][g]*v_array[ek][ej+j][ei+i];
						}
					}
				}
			}
			/*			Accumulate		*/
			for (l=0,k = 0; k < e->nphiz; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphix; i++, l++) {
						RHS_local[l] = 0.;
						for (g = 0; g < e->ng; g++) {
							RHS_local[l] += -1*e->weight[g]*e->phi[k][j][i][g]*(pow(v_elem[g],2))*flux_elem[g];
						}
					}
				}
			}
			break;
		case Z1:
			for (k = 0; k < e->nphiz; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphix; i++) {
						for (g = 0; g < e->ng; g++) {
							flux_elem[g] += e->phi[k][j][i][g]*flux_array[ek+1][ej+j][ei+i][2];
							v_elem[g] += e->phi[k][j][i][g]*v_array[ek+1][ej+j][ei+i];
						}
					}
				}
			}
			/*			Accumulate		*/
			for (l=0,k = 0; k < e->nphiz; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphix; i++, l++) {
						RHS_local[l] = 0.;
						for (g = 0; g < e->ng; g++) {
							RHS_local[l] += e->weight[g]*e->phi[k][j][i][g]*(pow(v_elem[g],2))*flux_elem[g];
						}
					}
				}
			}
			break;
	}  
  ierr = PetscFree2(flux_elem,v_elem);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_HeatMatK_local"
extern PetscErrorCode VF_HeatMatK_local(PetscReal *Kd_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ****diffsvty_array,PetscReal ***v_array)
{
  PetscErrorCode  ierr;
	PetscInt        i,j,k,l;
	PetscInt        ii,jj,kk;
	PetscInt        eg;
	PetscReal       kx,ky,kz,kxy,kxz,kyz;
  PetscReal		   *v_elem;
  
	PetscFunctionBegin;
  ierr = PetscMalloc(e->ng*sizeof(PetscReal),&v_elem);CHKERRQ(ierr);
	kx     = diffsvty_array[ek][ej][ei][0];
	ky     = diffsvty_array[ek][ej][ei][1];
	kz     = diffsvty_array[ek][ej][ei][2];
	kxy    = diffsvty_array[ek][ej][ei][3];
	kxz    = diffsvty_array[ek][ej][ei][4];
	kyz    = diffsvty_array[ek][ej][ei][5];
  for (eg = 0; eg < e->ng; eg++){
		v_elem[eg] = 0.;
	}
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          v_elem[eg] += v_array[ek+k][ej+j][ei+i] * e->phi[k][j][i][eg];
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
							Kd_ele[l] = 0.;
							for (eg = 0; eg < e->ng; eg++) {
								Kd_ele[l] += ((kx*e->dphi[kk][jj][ii][0][eg]+kxy*e->dphi[kk][jj][ii][1][eg]+kxz*e->dphi[kk][jj][ii][2][eg])*e->dphi[k][j][i][0][eg]
                              +(kxy*e->dphi[kk][jj][ii][0][eg]+ky*e->dphi[kk][jj][ii][1][eg]+kyz*e->dphi[kk][jj][ii][2][eg])*e->dphi[k][j][i][1][eg]
                              +(kxz*e->dphi[kk][jj][ii][0][eg]+kyz*e->dphi[kk][jj][ii][1][eg]+kz*e->dphi[kk][jj][ii][2][eg])*e->dphi[k][j][i][2][eg])*(pow(v_elem[eg],2))*e->weight[eg];
     
/*                The line below highlights the fact that v^2 was removed from above
                +(kxz*e->dphi[kk][jj][ii][0][eg]+kyz*e->dphi[kk][jj][ii][1][eg]+kz*e->dphi[kk][jj][ii][2][eg])*e->dphi[k][j][i][2][eg])*(pow(v_elem[eg],2))*e->weight[eg];  */

							}
						}
					}
				}
			}
		}
	}
  ierr = PetscFree(v_elem);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VF_HeatMatN_local"
extern PetscErrorCode VF_HeatMatN_local(PetscReal *KN_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei, PetscReal ****vel_array,PetscReal ***v_array)
{
	PetscErrorCode ierr;
	PetscInt i,j,k,l;
	PetscInt ii,jj,kk;
	PetscInt eg;
	PetscReal		*velx_loc;
	PetscReal		*vely_loc;
	PetscReal		*velz_loc;  
	PetscReal		*dvelx_dx_loc;  
	PetscReal		*dvely_dy_loc;  
	PetscReal		*dvelz_dz_loc;
	PetscReal		*v_elem;

	PetscFunctionBegin;
	ierr = PetscMalloc3(e->ng,PetscReal,&velx_loc,e->ng,PetscReal,&vely_loc,e->ng,PetscReal,&velz_loc);CHKERRQ(ierr);
	ierr = PetscMalloc4(e->ng,PetscReal,&dvelx_dx_loc,e->ng,PetscReal,&dvely_dy_loc,e->ng,PetscReal,&dvelz_dz_loc,e->ng,PetscReal,&v_elem);CHKERRQ(ierr);
	for (eg = 0; eg < e->ng; eg++){
		velx_loc[eg] = 0.;
		vely_loc[eg] = 0.;
		velz_loc[eg] = 0.;
		
		dvelx_dx_loc[eg] = 0.;
		dvely_dy_loc[eg] = 0.;
		dvelz_dz_loc[eg] = 0.;
		v_elem[eg] = 0.;
		
	}
	for (eg = 0; eg < e->ng; eg++) {
		for (k = 0; k < e->nphiz; k++) {
			for (j = 0; j < e->nphiy; j++) {
				for (i = 0; i < e->nphix; i++) {
					velx_loc[eg] += vel_array[ek+k][ej+j][ei+i][0] * e->phi[k][j][i][eg];
					vely_loc[eg] += vel_array[ek+k][ej+j][ei+i][1] * e->phi[k][j][i][eg];
					velz_loc[eg] += vel_array[ek+k][ej+j][ei+i][2] * e->phi[k][j][i][eg];
					dvelx_dx_loc[eg] += vel_array[ek+k][ej+j][ei+i][0] * e->dphi[k][j][i][0][eg];
					dvely_dy_loc[eg] += vel_array[ek+k][ej+j][ei+i][1] * e->dphi[k][j][i][1][eg];
					dvelz_dz_loc[eg] += vel_array[ek+k][ej+j][ei+i][2] * e->dphi[k][j][i][2][eg];
					
          v_elem[eg] += v_array[ek+k][ej+j][ei+i] * e->phi[k][j][i][eg];

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
							KN_ele[l] = 0.;
							for (eg = 0; eg < e->ng; eg++) {
								KN_ele[l] += ( e->phi[k][j][i][eg]*(dvelx_dx_loc[eg]+dvely_dy_loc[eg]+dvelz_dz_loc[eg]) + (velx_loc[eg]*e->dphi[k][j][i][0][eg]
																														   +vely_loc[eg]*e->dphi[k][j][i][1][eg]+velz_loc[eg]*e->dphi[k][j][i][2][eg]) )
								*(velx_loc[eg]*e->dphi[kk][jj][ii][0][eg]
								  +vely_loc[eg]*e->dphi[kk][jj][ii][1][eg]+velz_loc[eg]*e->dphi[kk][jj][ii][2][eg])*(pow(v_elem[eg],2))*e->weight[eg];
							}
						}
					}
				}
			}
		}
	}
	ierr = PetscFree3(dvelx_dx_loc,dvely_dy_loc,dvelz_dz_loc);CHKERRQ(ierr);
	ierr = PetscFree4(velx_loc,vely_loc,velz_loc,v_elem);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_HeatMatC_local"
extern PetscErrorCode VF_HeatMatC_local(PetscReal *KC_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei, PetscReal ****vel_array,PetscReal ***v_array)
{
  
	PetscErrorCode ierr;
	PetscInt i,j,k,l;
	PetscInt ii,jj,kk;
	PetscInt eg;
	PetscReal		*velx_loc;
	PetscReal		*vely_loc;
	PetscReal		*velz_loc;
	PetscReal		*v_elem;
  
	PetscFunctionBegin;
	ierr = PetscMalloc4(e->ng,PetscReal,&velx_loc,e->ng,PetscReal,&vely_loc,e->ng,PetscReal,&velz_loc,e->ng,PetscReal,&v_elem);CHKERRQ(ierr);
	for (eg = 0; eg < e->ng; eg++){
		velx_loc[eg] = 0.;
		vely_loc[eg] = 0.;
		velz_loc[eg] = 0.;
		v_elem[eg] = 0.;
	}
	for (eg = 0; eg < e->ng; eg++) {
		for (k = 0; k < e->nphiz; k++) {
			for (j = 0; j < e->nphiy; j++) {
				for (i = 0; i < e->nphix; i++) {
					velx_loc[eg] += vel_array[ek+k][ej+j][ei+i][0] * e->phi[k][j][i][eg];
					vely_loc[eg] += vel_array[ek+k][ej+j][ei+i][1] * e->phi[k][j][i][eg];
					velz_loc[eg] += vel_array[ek+k][ej+j][ei+i][2] * e->phi[k][j][i][eg];
					v_elem[eg] += v_array[ek+k][ej+j][ei+i] * e->phi[k][j][i][eg];
				}
			}
		}
/*    ierr = PetscPrintf(PETSC_COMM_WORLD,"eg = %d \t velx = %g \t vely = %g velz = %g \n",eg, velx_loc[eg],vely_loc[eg],velz_loc[eg]);CHKERRQ(ierr);   */
	}
	for (l = 0,k = 0; k < e->nphiz; k++) {
		for (j = 0; j < e->nphiy; j++) {
			for (i = 0; i < e->nphix; i++) {
				for (kk = 0; kk < e->nphiz; kk++) {
					for (jj = 0; jj < e->nphiy; jj++) {
						for (ii = 0; ii < e->nphix; ii++,l++) {
							KC_ele[l] = 0.;
							for (eg = 0; eg < e->ng; eg++) {
								KC_ele[l] += e->phi[k][j][i][eg]*(velx_loc[eg]*e->dphi[kk][jj][ii][0][eg]
																  +vely_loc[eg]*e->dphi[kk][jj][ii][1][eg]+velz_loc[eg]*e->dphi[kk][jj][ii][2][eg])*(pow(v_elem[eg],2))*e->weight[eg];
                
                
/*                KC_ele[l] += e->phi[kk][jj][ii][eg]*(velx_loc[eg]*e->dphi[k][j][i][0][eg])*(pow(v_elem[eg],2))*e->weight[eg]; */
                
							}
						}
					}
				}
			}
		}
	}
	ierr = PetscFree4(velx_loc,vely_loc,velz_loc,v_elem);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecApplyHeatSourceTerms"
extern PetscErrorCode VecApplyHeatSourceTerms(PetscReal *K1_local,PetscReal *K2_local,PetscReal ***source_array,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,VFCtx *ctx, PetscReal ****vel_array,PetscReal ***v_array)
{
	PetscErrorCode ierr;
	PetscInt       i,j,k,l;
	PetscInt       eg;
	PetscReal      *loc_source;
	PetscReal		*dsource_dx_loc;
	PetscReal		*dsource_dy_loc;
	PetscReal		*dsource_dz_loc;
	PetscReal		*velx_loc;
	PetscReal		*vely_loc;
	PetscReal		*velz_loc;
	PetscReal		*v_elem;
	
	PetscFunctionBegin;
	ierr = PetscMalloc4(e->ng,PetscReal,&velx_loc,e->ng,PetscReal,&vely_loc,e->ng,PetscReal,&velz_loc,e->ng,PetscReal,&v_elem);CHKERRQ(ierr);
	ierr = PetscMalloc4(e->ng,PetscReal,&loc_source,e->ng,PetscReal,&dsource_dx_loc,e->ng,PetscReal,&dsource_dy_loc,e->ng,PetscReal,&dsource_dz_loc);CHKERRQ(ierr);
	for (eg = 0; eg < e->ng; eg++) {
		loc_source[eg] = 0.;
		dsource_dx_loc[eg] = 0.;
		dsource_dy_loc[eg] = 0.;
		dsource_dz_loc[eg] = 0.;
		velx_loc[eg] = 0.;
		vely_loc[eg] = 0.;
		velz_loc[eg] = 0.;
		v_elem[eg] = 0.;
	}
	for (k = 0; k < e->nphiz; k++) {
		for (j = 0; j < e->nphiy; j++) {
			for (i = 0; i < e->nphix; i++) {
				for (eg = 0; eg < e->ng; eg++) {
					loc_source[eg] += source_array[ek+k][ej+j][ei+i]*e->phi[k][j][i][eg];
					dsource_dx_loc[eg] += source_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][0][eg];
					dsource_dx_loc[eg] += source_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][1][eg];
					dsource_dx_loc[eg] += source_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][2][eg];
					velx_loc[eg] += vel_array[ek+k][ej+j][ei+i][0] * e->phi[k][j][i][eg];
					vely_loc[eg] += vel_array[ek+k][ej+j][ei+i][1] * e->phi[k][j][i][eg];
					velz_loc[eg] += vel_array[ek+k][ej+j][ei+i][2] * e->phi[k][j][i][eg];
          v_elem[eg] += v_array[ek+k][ej+j][ei+i] * e->phi[k][j][i][eg];
				}
			}
		}
	}
	for (l = 0,k = 0; k < e->nphiz; k++) {
		for (j = 0; j < e->nphiy; j++) {
			for (i = 0; i < e->nphix; i++,l++) {
				K1_local[l] = 0.;
				K2_local[l] = 0.;
				for (eg = 0; eg < e->ng; eg++) {
					K1_local[l] += loc_source[eg]*e->phi[k][j][i][eg]*(pow(v_elem[eg],2))*e->weight[eg];
					K2_local[l] += (velx_loc[eg]*dsource_dx_loc[eg]+vely_loc[eg]*dsource_dy_loc[eg]+velz_loc[eg]*dsource_dz_loc[eg])*e->phi[k][j][i][eg]*(pow(v_elem[eg],2))*e->weight[eg];
				}
			}
		}
	}
	ierr = PetscFree4(velx_loc,vely_loc,velz_loc,v_elem);CHKERRQ(ierr);
	ierr = PetscFree4(loc_source,dsource_dx_loc,dsource_dy_loc,dsource_dz_loc);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


















