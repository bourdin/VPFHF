/*
 VFFlow_MixedFEM.c
 A mixed finite elements Darcy solver based on the method in
 Masud, A. and Hughes, T. J. (2002). A stabilized mixed finite element method for
 Darcy flow. Computer Methods in Applied Mechanics and Engineering, 191(3940):43414370.
 
 (c) 2011-2012 C. Chukwudozie, LSU
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFFlow.h"
/* #include "PetscFixes.h" */
#include "VFFlow_KSPMixedFEM.h"

/*
 MixedFlowFEMKSPSolve
 */


#undef __FUNCT__
#define __FUNCT__ "MixedFEMFlowSolverInitialize"
extern PetscErrorCode MixedFEMFlowSolverInitialize(VFCtx *ctx, VFFields *fields)
{
  PetscErrorCode      ierr;
  
  PetscFunctionBegin;
  
  ierr = KSPCreate(PETSC_COMM_WORLD,&ctx->kspVelP);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ctx->kspVelP,1.e-6,1.e-6,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = KSPSetOperators(ctx->kspVelP,ctx->KVelP,ctx->KVelP,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ctx->kspVelP,PETSC_TRUE);CHKERRQ(ierr);
  ierr = KSPAppendOptionsPrefix(ctx->kspVelP,"Flowksp_");CHKERRQ(ierr);
  ierr = KSPSetType(ctx->kspVelP,KSPBCGSL);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ctx->kspVelP);CHKERRQ(ierr);
  ierr = KSPGetPC(ctx->kspVelP,&ctx->pcVelP);CHKERRQ(ierr);
  ierr = PCSetType(ctx->pcVelP,PCJACOBI);CHKERRQ(ierr);
  ierr = PCSetFromOptions(ctx->pcVelP);CHKERRQ(ierr);
  
  /*
   ierr = KSPCreate(PETSC_COMM_WORLD,&ctx->kspVelP);CHKERRQ(ierr);
   ierr = KSPSetTolerances(ctx->kspVelP,1.e-6,1.e-6,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
   ierr = KSPSetOperators(ctx->kspVelP,ctx->KVelP,ctx->KVelP,flg);CHKERRQ(ierr);
   ierr = KSPSetInitialGuessNonzero(ctx->kspVelP,PETSC_TRUE);CHKERRQ(ierr);
   ierr = KSPAppendOptionsPrefix(ctx->kspVelP,"VelP_");CHKERRQ(ierr);
   ierr = KSPSetType(ctx->kspVelP,KSPFGMRES);CHKERRQ(ierr);
   ierr = KSPGetPC(ctx->kspVelP,&ctx->pcVelP);CHKERRQ(ierr);
   const PetscInt ufields[] = {0,1,2},pfields[] = {3};
   ierr = PCSetType(ctx->pcVelP,PCFIELDSPLIT);CHKERRQ(ierr);
   ierr = PCFieldSplitSetBlockSize(ctx->pcVelP,4);
   ierr = PCFieldSplitSetFields(ctx->pcVelP,"u",3,ufields,ufields);CHKERRQ(ierr);
   ierr = PCFieldSplitSetFields(ctx->pcVelP,"p",1,pfields,pfields);CHKERRQ(ierr);
   ierr = PetscOptionsSetValue("-VelP_fieldsplit_u_pc_type","mg");CHKERRQ(ierr); //mg is a type of SOR (see -VelP_ksp_view)
   ierr = PetscOptionsSetValue("-VelP_fieldsplit_p_pc_type","none");CHKERRQ(ierr);
   ierr = PetscOptionsSetValue("-VelP_fieldsplit_u_ksp_type","bcgsl");CHKERRQ(ierr);
   ierr = PetscOptionsSetValue("-VelP_fieldsplit_p_ksp_type","fgmres");CHKERRQ(ierr);
   ierr = PetscOptionsSetValue("-VelP_pc_fieldsplit_type","schur");CHKERRQ(ierr);
   ierr = PetscOptionsSetValue("-VelP_ksp_atol","5e-9");CHKERRQ(ierr);
   ierr = PetscOptionsSetValue("-VelP_ksp_rtol","5e-9");CHKERRQ(ierr);
   //     ierr = PetscOptionsSetValue("-VelP_pc_ﬁeldsplit_schur_precondition","self");CHKERRQ(ierr);
   ierr = PetscOptionsSetValue("-VelP_pc_fieldsplit_detect_saddle_point","true");CHKERRQ(ierr);
   ierr = PCSetFromOptions(ctx->pcVelP);CHKERRQ(ierr);
   ierr = KSPSetFromOptions(ctx->kspVelP);CHKERRQ(ierr);
   */
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MixedFEMFlowSolverFinalize"
extern PetscErrorCode MixedFEMFlowSolverFinalize(VFCtx *ctx,VFFields *fields)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = KSPDestroy(&ctx->kspVelP);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MixedFlowFEMKSPSolve"
extern PetscErrorCode MixedFlowFEMKSPSolve(VFCtx *ctx,VFFields *fields)
{
  PetscErrorCode     ierr;
  PetscInt           xs,xm,ys,ym,zs,zm;
  PetscInt           i,j,k,c,veldof = 3;
  PetscInt           its;
  KSPConvergedReason reason;
  PetscReal          ****VelnPress_array;
  PetscReal          ***Press_array;
  Vec                VecRHS;
  PetscReal          ****vel_array;
  PetscReal          theta,one_minus_theta;
  Vec                vec;
	PetscReal           Velmin,Velmax;
	PetscReal           Pmin,Pmax;
  
  PetscFunctionBegin;
  theta = ctx->flowprop.theta;
  one_minus_theta = (1.-theta);
  ierr = VecDuplicate(ctx->RHSVelP,&VecRHS);CHKERRQ(ierr);
  ierr = VecDuplicate(ctx->RHSVelP,&vec);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daFlow,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = VecCopy(fields->V,ctx->V);CHKERRQ(ierr);
	ierr = VecCopy(fields->U,ctx->U);CHKERRQ(ierr);
  ierr = FlowMatnVecAssemble(ctx->KVelP,ctx->KVelPlhs,ctx->RHSVelP,fields,ctx);CHKERRQ(ierr);
  ierr = VecCopy(ctx->RHSVelP,VecRHS);CHKERRQ(ierr);
  ierr = VecAXPBY(VecRHS,one_minus_theta,theta,ctx->RHSVelPpre);CHKERRQ(ierr);
  ierr = MatMultAdd(ctx->KVelPlhs,ctx->PreFlowFields,VecRHS,VecRHS);CHKERRQ(ierr);
  ierr = VecApplyVelocityBC(VecRHS,ctx->VelBCArray,&ctx->bcQ[0],ctx);CHKERRQ(ierr);
  if (ctx->verbose > 1) {
    ierr = KSPMonitorSet(ctx->kspVelP,MixedFEMKSPMonitor,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  }
  ierr = KSPSolve(ctx->kspVelP,VecRHS,fields->VelnPress);CHKERRQ(ierr);
  
/*
  PetscViewer viewer;
   ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"RHSvec.txt",&viewer);CHKERRQ(ierr);
   ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
   ierr = VecView(vec,viewer);CHKERRQ(ierr);
   
   
   ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"Matrix.txt",&viewer);CHKERRQ(ierr);
   ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
   ierr = MatView(ctx->KVelP,viewer);CHKERRQ(ierr);
   
   ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"Matrixlhs.txt",&viewer);CHKERRQ(ierr);
   ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
   ierr = MatView(ctx->KVelPlhs,viewer);CHKERRQ(ierr);
   
   ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"RHS.txt",&viewer);CHKERRQ(ierr);
   ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
   ierr = VecView(VecRHS,viewer);CHKERRQ(ierr);
   
   ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"RHSor.txt",&viewer);CHKERRQ(ierr);
   ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
   ierr = VecView(ctx->RHSVelP,viewer);CHKERRQ(ierr);
   
   
   ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"Solution.txt",&viewer);CHKERRQ(ierr);
   ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
   ierr = VecView(fields->VelnPress,viewer);CHKERRQ(ierr);
*/
  ierr = KSPGetConvergedReason(ctx->kspVelP,&reason);CHKERRQ(ierr);
  if (reason < 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"[ERROR] ksp_MixedFlowSolver diverged with reason %d\n",(int)reason);CHKERRQ(ierr);
  } else {
    ierr = KSPGetIterationNumber(ctx->kspVelP,&its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"      ksp_MixedFlowSolver converged in %d iterations %d.\n",(int)its,(int)reason);CHKERRQ(ierr);
  }
  ierr = DMDAVecGetArray(ctx->daScal,fields->pressure,&Press_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daVect,fields->velocity,&vel_array);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daFlow,fields->VelnPress,&VelnPress_array);CHKERRQ(ierr);
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
  ierr = DMDAVecRestoreArray(ctx->daScal,fields->pressure,&Press_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,fields->velocity,&vel_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daFlow,fields->VelnPress,&VelnPress_array);CHKERRQ(ierr);
  ierr = VecDestroy(&VecRHS);CHKERRQ(ierr);
  ierr = VecDestroy(&vec);CHKERRQ(ierr);
  ierr = VecMin(fields->velocity,PETSC_NULL,&Velmin);CHKERRQ(ierr);
	ierr = VecMax(fields->velocity,PETSC_NULL,&Velmax);CHKERRQ(ierr);
  ierr = VecMin(fields->pressure,PETSC_NULL,&Pmin);CHKERRQ(ierr);
	ierr = VecMax(fields->pressure,PETSC_NULL,&Pmax);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"      Velocity min / max:     %e %e\n",Velmin,Velmax);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"      Pressure min / max:     %e %e\n",Pmin,Pmax);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MixedFEMKSPMonitor"
extern PetscErrorCode MixedFEMKSPMonitor(KSP ksp,PetscInt its,PetscReal fnorm,void* ptr)
{
  PetscErrorCode ierr;
  PetscReal      norm,vmax,vmin;
  MPI_Comm       comm;
  Vec                             solution;
  
  PetscFunctionBegin;
  ierr = KSPBuildSolution(ksp,PETSC_NULL,&solution);CHKERRQ(ierr);
  ierr = VecNorm(solution,NORM_1,&norm);CHKERRQ(ierr);
  ierr = VecMax(solution,PETSC_NULL,&vmax);CHKERRQ(ierr);
  ierr = VecMin(solution,PETSC_NULL,&vmin);CHKERRQ(ierr);
  ierr = PetscObjectGetComm((PetscObject)ksp,&comm);CHKERRQ(ierr);
  ierr = PetscPrintf(comm,"ksp_iter_step %D :solution norm = %G, max sol. value  = %G, min sol. value = %G\n",its,norm,vmax,vmin);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecApplyFractureWellSource"
extern PetscErrorCode VecApplyFractureWellSource(PetscReal *Ks_local,PetscReal ***source_array,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,VFCtx *ctx,PetscReal ***v_array)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,l,c;
  PetscReal      alpha_c;
  PetscReal      mu;
  PetscReal      *loc_source,*dv_elem[3],*v_elem,*v_mag_elem;
  PetscInt       eg;
  
  PetscFunctionBegin;
  alpha_c = ctx->flowprop.alpha;
  mu     = ctx->flowprop.mu;
  ierr = PetscMalloc6(e->ng,PetscReal,&loc_source,e->ng,PetscReal,&v_elem,e->ng,PetscReal,&dv_elem[0],e->ng,PetscReal,&dv_elem[1],e->ng,PetscReal,&dv_elem[2],e->ng,PetscReal,&v_mag_elem);CHKERRQ(ierr);
  for (eg = 0; eg < e->ng; eg++){
    loc_source[eg] = 0.;
    v_elem[eg] = 0.;
    v_mag_elem[eg] = 0.;
    for(c = 0; c < 3; c++){
      dv_elem[c][eg] = 0.;
    }
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          loc_source[eg] += source_array[ek+k][ej+j][ei+i]*e->phi[k][j][i][eg];
          v_elem[eg] += v_array[ek+k][ej+j][ei+i] * e->phi[k][j][i][eg];
          for(c = 0; c < 3; c++){
            dv_elem[c][eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
          }
        }
      }
    }
  }
  for (eg = 0; eg < e->ng; eg++){
    v_mag_elem[eg] = sqrt((pow(dv_elem[0][eg],2))+(pow(dv_elem[1][eg],2))+(pow(dv_elem[2][eg],2)));
  }
  for (l = 0,k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++,l++) {
        Ks_local[l] = 0.;
        for (eg = 0; eg < e->ng; eg++) {
          Ks_local[l] += -loc_source[eg]*e->phi[k][j][i][eg]*v_mag_elem[eg]*e->weight[eg]/alpha_c;
        }
      }
    }
  }
  /*
  for (l = 0,k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++,l++) {
        if(v_array[ek+k][ej+j][ei+i] > 0.05){
          Ks_local[l] = 0.;
        }
      }
    }
  }
  */
  ierr = PetscFree6(loc_source,v_elem,dv_elem[0],dv_elem[1],dv_elem[2],v_mag_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecApplySourceTerms"
extern PetscErrorCode VecApplySourceTerms(PetscReal *Ks_local,PetscReal ***source_array,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,VFCtx *ctx,PetscReal ***v_array)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,l;
  PetscReal      alpha_c;
  PetscReal      mu;
  PetscReal      *loc_source, *v_elem;
  PetscInt       eg;
  
  PetscFunctionBegin;
  alpha_c = ctx->flowprop.alpha;
  mu     = ctx->flowprop.mu;
  ierr = PetscMalloc2(e->ng,PetscReal,&loc_source,e->ng,PetscReal,&v_elem);CHKERRQ(ierr);
  for (eg = 0; eg < e->ng; eg++){
    loc_source[eg] = 0.;
    v_elem[eg] = 0.;
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          loc_source[eg] += source_array[ek+k][ej+j][ei+i]*e->phi[k][j][i][eg];
          v_elem[eg] += v_array[ek+k][ej+j][ei+i] * e->phi[k][j][i][eg];
        }
      }
    }
  }
  for (l = 0,k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++,l++) {
        Ks_local[l] = 0.;
        for (eg = 0; eg < e->ng; eg++) {
          Ks_local[l] += -loc_source[eg]*e->phi[k][j][i][eg]*(pow(v_elem[eg],2))*e->weight[eg]/alpha_c;
        }
      }
    }
  }
  /*
  for (l = 0,k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++,l++) {
        if(v_array[ek+k][ej+j][ei+i] > 0.05 ){
          Ks_local[l] = 0.;
        }
      }
    }
  }
  */
  ierr = PetscFree2(loc_source,v_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecApplyVelocityBC"
extern PetscErrorCode VecApplyVelocityBC(Vec RHS,Vec BCV, BC *bcQ,VFCtx *ctx)
{
	PetscErrorCode ierr;
	PetscInt       xs,xm,nx;
	PetscInt       ys,ym,ny;
	PetscInt       zs,zm,nz;
	PetscInt       dim,dof=3;
	PetscInt       i,j,k,c;
	PetscReal		****velbc_array;
	PetscReal		****RHS_array;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(ctx->daFlow,&dim,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daFlow,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,BCV,&velbc_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daFlow,RHS,&RHS_array);CHKERRQ(ierr);
	for (c = 0; c < dof; c++) {
		if (xs == 0) {
			i = 0;
			if (bcQ[c].face[X0] == FIXED) {
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
					}
				}
			}
		}
		if (xs+xm == nx) {
			i = nx-1;
			if (bcQ[c].face[X1] == FIXED) {
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
					}
				}
			}
		}
		if (ys == 0) {
			j = 0;
			if (bcQ[c].face[Y0] == FIXED) {
				for (k = zs; k < zs+zm; k++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
					}
				}
			}
		}
		if (ys+ym == ny) {
			j = ny-1;
			if (bcQ[c].face[Y1] == FIXED) {
				for (k = zs; k < zs+zm; k++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
						
					}
				}
			}
		}
		if (zs == 0) {
			k = 0;
			if (bcQ[c].face[Z0] == FIXED) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
					}
				}
			}
		}
		if (zs+zm == nz) {
			k = nz-1;
			if (bcQ[c].face[Z1] == FIXED) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
					}
				}
			}
		}
		if (xs == 0 && zs == 0) {
			k = 0;i = 0;
			if (bcQ[c].edge[X0Z0] == FIXED) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && zs == 0) {
			k = 0;i = nx-1;
			if (bcQ[c].edge[X1Z0] == FIXED) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (ys == 0 && zs == 0) {
			k = 0;j = 0;
			if (bcQ[c].edge[Y0Z0] == FIXED) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (ys+ym == ny && zs == 0) {
			k = 0;j = 0;
			if (bcQ[c].edge[Y1Z0] == FIXED) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && zs+zm == nz) {
			k = nz-1;i = 0;
			if (bcQ[c].edge[X0Z1] == FIXED) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && zs+zm == nz) {
			k = nz-1;i = nx-1;
			if (bcQ[c].edge[X1Z1] == FIXED) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (ys == 0 && zs+zm == nz) {
			k = nz-1;j = 0;
			if (bcQ[c].edge[Y0Z1] == FIXED) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (ys+ym == ny && zs+zm == nz) {
			k = nz-1;j = ny-1;
			if (bcQ[c].edge[Y1Z1] == FIXED) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && ys == 0) {
			j = 0;i = 0;
			if (bcQ[c].edge[X0Y0] == FIXED) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && ys+ym == ny) {
			j = ny-1;i = 0;
			if (bcQ[c].edge[X0Y1] == FIXED) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && ys == 0) {
			j = 0;i = nx-1;
			if (bcQ[c].edge[X1Y0] == FIXED) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && ys+ym == ny) {
			j = ny-1;i = nx-1;
			if (bcQ[c].edge[X1Y1] == FIXED) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && ys == 0 && zs == 0) {
			k = 0;j = 0;i = 0;
			if (bcQ[c].vertex[X0Y0Z0] == FIXED) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys == 0 && zs == 0) {
			k = 0;j = 0;i = nx-1;
			if (bcQ[c].vertex[X1Y0Z0] == FIXED) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs == 0 && ys+ym == ny && zs == 0) {
			k = 0;j = ny-1;i = 0;
			if (bcQ[c].vertex[X0Y1Z0] == FIXED) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys+ym == ny && zs == 0) {
			k = 0;j = ny-1;i = nx-1;
			if (bcQ[c].vertex[X1Y1Z0] == FIXED) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs == 0 && ys == 0 && zs+zm == nz) {
			k = nz-1;j = 0;i = 0;
			if (bcQ[c].vertex[X0Y0Z1] == FIXED) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys == 0 && zs+zm == nz) {
			k = nz-1;j = 0;i = nx-1;
			if (bcQ[c].vertex[X1Y0Z1] == FIXED) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs == 0 && ys+ym == ny && zs+zm == nz) {
			k = nz-1;j = ny-1;i = 0;
			
			if (bcQ[c].vertex[X0Y1Z1] == FIXED) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys+ym == ny && zs+zm == nz) {
			k = nz-1;j = ny-1;i = nx-1;
			if (bcQ[c].vertex[X1Y1Z1] == FIXED) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
	}
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,BCV,&velbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daFlow,RHS,&RHS_array);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FlowMatnVecAssemble"
extern PetscErrorCode FlowMatnVecAssemble(Mat K,Mat Krhs,Vec RHS,VFFields *fields,VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       ek,ej,ei;
  PetscInt       i,j,k,l;
  PetscInt       veldof = 3;
  PetscInt       c;
  PetscReal      ****perm_array;
  PetscReal      ****coords_array;
  PetscReal      ****RHS_array;
  PetscReal      *RHS_local;
  Vec            RHS_localVec;
  Vec            perm_local;
  PetscReal      hx,hy,hz;
  PetscReal      *KA_local,*KB_local,*KD_local,*KBTrans_local,*KS_local;
  PetscReal      *KArhs_local,*KBrhs_local,*KDrhs_local,*KBTransrhs_local;
  PetscReal      *KL_local,*KLrhs_local;
  PetscReal      *KAf_local,*KBf_local,*KDf_local,*KBfTrans_local;
  PetscReal      *KAfrhs_local,*KBfrhs_local,*KDfrhs_local,*KBfTransrhs_local;
  PetscReal      beta_c,alpha_c,mu,gx,gy,gz;
  PetscReal      theta,timestepsize;
  PetscInt       nrow = ctx->e3D.nphix*ctx->e3D.nphiy*ctx->e3D.nphiz;
  MatStencil     *row,*row1;
  PetscReal      ***source_array;
  Vec            source_local;
  PetscReal      M_inv;
  PetscReal      one_minus_theta;
  FACE           face;
  PetscReal      ***prebc_array;
  Vec            prebc_local;
  MatStructure   flg;
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
  PetscReal      ***one_array;
  Vec            one_local;
  Vec            Ones;
  PetscReal      ***fracflow_array;
  Vec            fracflow_local;
  PetscReal      alphabiot;
  PetscReal      K_dr;
  PetscReal      ***pressure_diff_array;
	Vec            Pressure_diff,pressure_diff_local;
  
  
  PetscFunctionBegin;
  flg = SAME_NONZERO_PATTERN;
  M_inv     = ctx->flowprop.M_inv;
  beta_c = ctx->flowprop.beta;
  alpha_c = ctx->flowprop.alpha;
  theta = ctx->flowprop.theta;
  timestepsize = ctx->flowprop.timestepsize;
  one_minus_theta = -1.*(1.-theta);
  alphabiot  = ctx->flowprop.alphabiot;
  K_dr  = ctx->flowprop.K_dr;
  mu     = ctx->flowprop.mu;
  gx     = ctx->flowprop.g[0];
  gy     = ctx->flowprop.g[1];
  gz     = ctx->flowprop.g[2];
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = MatZeroEntries(K);CHKERRQ(ierr);
  ierr = MatZeroEntries(Krhs);CHKERRQ(ierr);
  ierr = VecSet(RHS,0.);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ctx->daFlow,&RHS_localVec);CHKERRQ(ierr);
  ierr = VecSet(RHS_localVec,0.);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daFlow,RHS_localVec,&RHS_array);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ctx->daScal,&source_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,ctx->Source,INSERT_VALUES,source_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,ctx->Source,INSERT_VALUES,source_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,source_local,&source_array);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daVFperm,fields->vfperm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daVFperm,fields->vfperm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daScal,&prebc_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,ctx->PresBCArray,INSERT_VALUES,prebc_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,ctx->PresBCArray,INSERT_VALUES,prebc_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,prebc_local,&prebc_array);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
  
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

  ierr = DMGetLocalVector(ctx->daScal,&v_old_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,ctx->V_old,INSERT_VALUES,v_old_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,ctx->V_old,INSERT_VALUES,v_old_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,v_old_local,&v_old_array);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(ctx->daScal,&Ones);CHKERRQ(ierr);
  ierr = VecSet(Ones,1.0);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daScal,&one_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,Ones,INSERT_VALUES,one_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,Ones,INSERT_VALUES,one_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,one_local,&one_array);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ctx->daScal,&fracflow_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,ctx->RegFracWellFlowRate,INSERT_VALUES,fracflow_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,ctx->RegFracWellFlowRate,INSERT_VALUES,fracflow_local);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,fracflow_local,&fracflow_array);CHKERRQ(ierr);
  
  ierr = PetscMalloc5(nrow*nrow,PetscReal,&KA_local,
                      nrow*nrow,PetscReal,&KB_local,
                      nrow*nrow,PetscReal,&KD_local,
                      nrow*nrow,PetscReal,&KBTrans_local,
                      nrow*nrow,PetscReal,&KS_local);CHKERRQ(ierr);
  ierr = PetscMalloc4(nrow*nrow,PetscReal,&KArhs_local,
                      nrow*nrow,PetscReal,&KBrhs_local,
                      nrow*nrow,PetscReal,&KDrhs_local,
                      nrow*nrow,PetscReal,&KBTransrhs_local);CHKERRQ(ierr);
  
  ierr = PetscMalloc5(nrow*nrow,PetscReal,&KAf_local,
                      nrow*nrow,PetscReal,&KBf_local,
                      nrow*nrow,PetscReal,&KDf_local,
                      nrow*nrow,PetscReal,&KBfTrans_local,
                      nrow*nrow,PetscReal,&KL_local);CHKERRQ(ierr);
  
  ierr = PetscMalloc5(nrow*nrow,PetscReal,&KAfrhs_local,
                      nrow*nrow,PetscReal,&KBfrhs_local,
                      nrow*nrow,PetscReal,&KDfrhs_local,
                      nrow*nrow,PetscReal,&KBfTransrhs_local,
                      nrow*nrow,PetscReal,&KLrhs_local);CHKERRQ(ierr);
  
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
        /*This computes the local contribution of the global A matrix*/
        ierr = VF_MatA_local(KS_local,&ctx->e3D,ek,ej,ei,v_array);CHKERRQ(ierr);
        for (l = 0; l < nrow*nrow; l++) {
          if(ctx->FlowDisplCoupling && ctx->ResFlowMechCoupling == FIXEDSTRESS){
            KS_local[l] = -1.*(M_inv+alphabiot*alphabiot/K_dr)*KS_local[l]/alpha_c;
          }
          else{
            KS_local[l] = -1.*M_inv*KS_local[l]/alpha_c;
          }
        }
        for (c = 0; c < veldof; c++) {
          ierr = Flow_MatA(KA_local,&ctx->e3D,ek,ej,ei,c,ctx->flowprop,perm_array,one_array);CHKERRQ(ierr);
//          ierr = Flow_MatB(KB_local,&ctx->e3D,ek,ej,ei,c,v_array);CHKERRQ(ierr);
//          ierr = Flow_MatBTranspose(KBTrans_local,&ctx->e3D,ek,ej,ei,c,v_array);CHKERRQ(ierr);
          ierr = Flow_MatB(KB_local,&ctx->e3D,ek,ej,ei,c,one_array);CHKERRQ(ierr);
          ierr = Flow_MatBTranspose(KBTrans_local,&ctx->e3D,ek,ej,ei,c,one_array);CHKERRQ(ierr);
          for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
            for (j = 0; j < ctx->e3D.nphiy; j++) {
              for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                row[l].i  = ei+i;row[l].j = ej+j;row[l].k = ek+k;row[l].c = c;
                row1[l].i = ei+i;row1[l].j = ej+j;row1[l].k = ek+k;row1[l].c = 3;
              }
            }
          }
          for (l = 0; l < nrow*nrow; l++) {
            KArhs_local[l] = one_minus_theta*KA_local[l];
            KA_local[l] = theta*KA_local[l];
            KBrhs_local[l] = one_minus_theta*KB_local[l];
            KB_local[l] = theta*KB_local[l];
            KBTransrhs_local[l] = timestepsize*one_minus_theta*KBTrans_local[l];
            KBTrans_local[l] = timestepsize*theta*KBTrans_local[l];
          }
          ierr = MatSetValuesStencil(K,nrow,row,nrow,row,KA_local,ADD_VALUES);CHKERRQ(ierr);
          ierr = MatSetValuesStencil(Krhs,nrow,row,nrow,row,KArhs_local,ADD_VALUES);CHKERRQ(ierr);
          
          ierr = MatSetValuesStencil(K,nrow,row,nrow,row1,KB_local,ADD_VALUES);CHKERRQ(ierr);
          ierr = MatSetValuesStencil(Krhs,nrow,row,nrow,row1,KBrhs_local,ADD_VALUES);CHKERRQ(ierr);
          
          ierr = MatSetValuesStencil(K,nrow,row1,nrow,row,KBTrans_local,ADD_VALUES);CHKERRQ(ierr);
          ierr = MatSetValuesStencil(Krhs,nrow,row1,nrow,row,KBTransrhs_local,ADD_VALUES);CHKERRQ(ierr);
        }
        ierr = Flow_MatD(KD_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,perm_array,one_array);CHKERRQ(ierr);
        for (l = 0; l < nrow*nrow; l++) {
          KDrhs_local[l] = timestepsize*one_minus_theta*KD_local[l];
          KD_local[l] = timestepsize*theta*KD_local[l];
        }
        
        ierr = MatSetValuesStencil(K,nrow,row1,nrow,row1,KD_local,ADD_VALUES);CHKERRQ(ierr);
        ierr = MatSetValuesStencil(Krhs,nrow,row1,nrow,row1,KDrhs_local,ADD_VALUES);CHKERRQ(ierr);
        
        ierr = MatSetValuesStencil(K,nrow,row1,nrow,row1,KS_local,ADD_VALUES);CHKERRQ(ierr);
        ierr = MatSetValuesStencil(Krhs,nrow,row1,nrow,row1,KS_local,ADD_VALUES);CHKERRQ(ierr);
        
        
        
        
        
        
        
   
        if(ctx->FractureFlowCoupling){
          for (c = 0; c < veldof; c++) {
            ierr = VF_MatAFractureFlowCoupling_local(KAf_local,&ctx->e3D,ek,ej,ei,u_array,v_array);CHKERRQ(ierr);
            ierr = VF_MatBTFractureFlowCoupling_local(KBfTrans_local,&ctx->e3D,ek,ej,ei,c,u_array,v_array);CHKERRQ(ierr);
            ierr = VF_MatBFractureFlowCoupling_local(KBf_local,&ctx->e3D,ek,ej,ei,c,u_array,v_array);CHKERRQ(ierr);
            ierr = VF_MatLeakOff_local(KL_local,&ctx->e3D,ek,ej,ei,c,v_array);CHKERRQ(ierr);
            for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                  row[l].i  = ei+i;row[l].j = ej+j;row[l].k = ek+k;row[l].c = c;
                  row1[l].i = ei+i;row1[l].j = ej+j;row1[l].k = ek+k;row1[l].c = 3;
                }
              }
            }
            for (l = 0; l < nrow*nrow; l++) {
              KAfrhs_local[l] = mu*one_minus_theta*KAf_local[l];
              KAf_local[l] = mu*theta*KAf_local[l];
              
              KBfrhs_local[l] = one_minus_theta*KBf_local[l]/12.0;
              KBf_local[l] = theta*KBf_local[l]/12.0;

              KBfTransrhs_local[l] = timestepsize*one_minus_theta*KBfTrans_local[l];
              KBfTrans_local[l] = timestepsize*theta*KBfTrans_local[l];

              KLrhs_local[l] = -1.0*timestepsize*one_minus_theta*KL_local[l];
              KL_local[l] = -1.0*timestepsize*theta*KL_local[l];
            }
            ierr = MatSetValuesStencil(K,nrow,row,nrow,row,KAf_local,ADD_VALUES);CHKERRQ(ierr);
            ierr = MatSetValuesStencil(Krhs,nrow,row,nrow,row,KAfrhs_local,ADD_VALUES);CHKERRQ(ierr);
            
  
            ierr = MatSetValuesStencil(K,nrow,row,nrow,row1,KBf_local,ADD_VALUES);CHKERRQ(ierr);
            ierr = MatSetValuesStencil(Krhs,nrow,row,nrow,row1,KBfrhs_local,ADD_VALUES);CHKERRQ(ierr);
 
            ierr = MatSetValuesStencil(K,nrow,row1,nrow,row,KBfTrans_local,ADD_VALUES);CHKERRQ(ierr);
            ierr = MatSetValuesStencil(Krhs,nrow,row1,nrow,row,KBfTransrhs_local,ADD_VALUES);CHKERRQ(ierr);

            ierr = MatSetValuesStencil(K,nrow,row1,nrow,row,KL_local,ADD_VALUES);CHKERRQ(ierr);
            ierr = MatSetValuesStencil(Krhs,nrow,row1,nrow,row,KLrhs_local,ADD_VALUES);CHKERRQ(ierr);
   
          }
          
          ierr = VF_MatDFractureFlowCoupling_local(KDf_local,&ctx->e3D,ek,ej,ei,u_array,v_array);CHKERRQ(ierr);
          for (l = 0; l < nrow*nrow; l++) {
            KDfrhs_local[l] = -timestepsize*one_minus_theta/(24.*mu)*KDf_local[l];
            KDf_local[l] = -timestepsize*theta/(24.*mu)*KDf_local[l];
            
          }
          ierr = MatSetValuesStencil(K,nrow,row1,nrow,row1,KDf_local,ADD_VALUES);CHKERRQ(ierr);
          ierr = MatSetValuesStencil(Krhs,nrow,row1,nrow,row1,KDfrhs_local,ADD_VALUES);CHKERRQ(ierr);
  
          if(ctx->hasFlowWells){
            ierr = VecApplyFractureWellSource(RHS_local,fracflow_array,&ctx->e3D,ek,ej,ei,ctx,v_array);
            for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                  RHS_array[ek+k][ej+j][ei+i][3] += timestepsize*RHS_local[l];
                }
              }
            }
          }
          
          ierr = VF_RHSFractureFlowCoupling_local(RHS_local,&ctx->e3D,ek,ej,ei,u_array,v_array,u_old_array,v_old_array);
          for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
            for (j = 0; j < ctx->e3D.nphiy; j++) {
              for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                RHS_array[ek+k][ej+j][ei+i][3] += RHS_local[l];
              }
            }
          }
        }
        
        
        
        
        
        
        
        
        /*Assembling the righthand side vector f*/
        for (c = 0; c < veldof; c++) {
          ierr = Flow_Vecf(RHS_local,&ctx->e3D,ek,ej,ei,c,ctx->flowprop,v_array);CHKERRQ(ierr);
          for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
            for (j = 0; j < ctx->e3D.nphiy; j++) {
              for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                RHS_array[ek+k][ej+j][ei+i][c] += RHS_local[l];
              }
            }
          }
        }
        /*Assembling the righthand side vector g*/
        ierr = Flow_Vecg(RHS_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,perm_array,one_array);CHKERRQ(ierr);
        for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
          for (j = 0; j < ctx->e3D.nphiy; j++) {
            for (i = 0; i < ctx->e3D.nphix; i++,l++) {
              RHS_array[ek+k][ej+j][ei+i][3] += timestepsize*RHS_local[l];
            }
          }
        }
        if(ctx->hasFluidSources){
          ierr = VecApplySourceTerms(RHS_local,source_array,&ctx->e3D,ek,ej,ei,ctx,v_array);
          for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
            for (j = 0; j < ctx->e3D.nphiy; j++) {
              for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                RHS_array[ek+k][ej+j][ei+i][3] += timestepsize*RHS_local[l];
              }
            }
          }
        }
        if(ctx->FlowDisplCoupling){
          ierr = VF_RHSFlowMechUCoupling_local(RHS_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,u_diff_array,v_array);
          for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
            for (j = 0; j < ctx->e3D.nphiy; j++) {
              for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                RHS_array[ek+k][ej+j][ei+i][3] += RHS_local[l];
              }
            }
          }
        }
        if(ctx->FlowDisplCoupling && ctx->ResFlowMechCoupling == FIXEDSTRESS){
          ierr = VF_RHSFlowMechUCouplingFIXSTRESS_local(RHS_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,ctx->matprop,pressure_diff_array,v_array);
          for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
            for (j = 0; j < ctx->e3D.nphiy; j++) {
              for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                RHS_array[ek+k][ej+j][ei+i][3] += -1*RHS_local[l];
              }
            }
          }
        }
        if (ei == 0) {
          /*                                       Face X0                        */
          face = X0;
          ierr = CartFE_Element2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
          if (ctx->bcP[0].face[face] == FIXED) {
            ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array,v_array);CHKERRQ(ierr);
            for (l=0,k = 0; k < ctx->e2D.nphix; k++){
              for (j = 0; j < ctx->e2D.nphiy; j++) {
                for (i = 0; i < ctx->e2D.nphiz; i++, l++) {
                  RHS_array[ek+k][ej+j][ei+i][0] -= RHS_local[l];
                }
              }
            }
          }
        }
        if (ei == nx-1) {
          /*                                       Face X1                */
          face = X1;
          ierr = CartFE_Element2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
          if (ctx->bcP[0].face[face] == FIXED) {
            ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array,v_array);CHKERRQ(ierr);
            for (l=0,k = 0; k < ctx->e2D.nphix; k++){
              for (j = 0; j < ctx->e2D.nphiy; j++) {
                for (i = 0; i < ctx->e2D.nphiz; i++, l++) {
                  RHS_array[ek+k][ej+j][ei+1][0] += RHS_local[l];
                }
              }
            }
          }
        }
        if (ej == 0) {
          /*                                       Face Y0                */
          face = Y0;
          ierr = CartFE_Element2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
          if (ctx->bcP[0].face[face] == FIXED) {
            ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array,v_array);CHKERRQ(ierr);
            for (l=0,k = 0; k < ctx->e2D.nphiy; k++){
              for (j = 0; j < ctx->e2D.nphiz; j++) {
                for (i = 0; i < ctx->e2D.nphix; i++, l++) {
                  RHS_array[ek+k][ej+j][ei+i][1] -= RHS_local[l];
                }
              }
            }
          }
        }
        if (ej == ny-1) {
          /*                                       Face Y1                */
          face = Y1;
          ierr = CartFE_Element2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
          if (ctx->bcP[0].face[face] == FIXED) {
            ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array,v_array);CHKERRQ(ierr);
            for (l=0,k = 0; k < ctx->e2D.nphiy; k++){
              for (j = 0; j < ctx->e2D.nphiz; j++) {
                for (i = 0; i < ctx->e2D.nphix; i++, l++) {
                  RHS_array[ek+k][ej+1][ei+i][1] += RHS_local[l];
                }
              }
            }
          }
        }
        if (ek == 0) {
          /*                                       Face Z0                */
          face = Z0;
          ierr = CartFE_Element2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
          if (ctx->bcP[0].face[face] == FIXED) {
            ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array,v_array);CHKERRQ(ierr);
            for (l=0,k = 0; k < ctx->e2D.nphiz; k++){
              for (j = 0; j < ctx->e2D.nphiy; j++) {
                for (i = 0; i < ctx->e2D.nphix; i++, l++) {
                  RHS_array[ek+k][ej+j][ei+i][2] -= RHS_local[l];
                }
              }
            }
          }
        }
        if (ek == nz-1) {
          /*                                       Face Z1                */
          face = Z1;
          ierr = CartFE_Element2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
          if (ctx->bcP[0].face[face] == FIXED) {
            ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array,v_array);CHKERRQ(ierr);
            for (l=0,k = 0; k < ctx->e2D.nphiz; k++){
              for (j = 0; j < ctx->e2D.nphiy; j++) {
                for (i = 0; i < ctx->e2D.nphix; i++, l++) {
                  RHS_array[ek+1][ej+j][ei+i][2] += RHS_local[l];
                }
              }
            }
          }
        }
      }
    }
  }
  PetscInt  w_no = 0;
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
              ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
              hwx = (ctx->well[w_no].coords[0]-coords_array[ek][ej][ei][0])/hx;
              hwy = (ctx->well[w_no].coords[1]-coords_array[ek][ej][ei][1])/hy;
              hwz = (ctx->well[w_no].coords[2]-coords_array[ek][ej][ei][2])/hz;
              if(ctx->well[w_no].condition == RATE){
                ierr = VecApplyWellFlowRate(RHS_local,&ctx->e3D,ctx->well[w_no].Qw,hwx,hwy,hwz,ek,ej,ei,v_array);CHKERRQ(ierr);
                for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
                  for (j = 0; j < ctx->e3D.nphiy; j++) {
                    for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                      if(ctx->well[w_no].type == INJECTOR){
                        RHS_array[ek+k][ej+j][ei+i][3] -= timestepsize*RHS_local[l];
                      }
                      else if(ctx->well[w_no].type == PRODUCER)
                      {
                        RHS_array[ek+k][ej+j][ei+i][3] -= -timestepsize*RHS_local[l];
                      }
                    }
                  }
                }
              }
              else if(ctx->well[w_no].condition == PRESSURE){
              }
              w_no++;
            }
          }
        }
      }
    }
    
 /*
    w_no = 0;
    while(w_no < ctx->numfracWells){
      for (ek = zs; ek < zs+zm; ek++) {
        for (ej = ys; ej < ys+ym; ej++) {
          for (ei = xs; ei < xs+xm; ei++) {
            if(
               ((coords_array[ek][ej][ei+1][0] >= ctx->fracwell[w_no].coords[0]) && (coords_array[ek][ej][ei][0] <= ctx->fracwell[w_no].coords[0] ))    &&
               ((coords_array[ek][ej+1][ei][1] >= ctx->fracwell[w_no].coords[1]) && (coords_array[ek][ej][ei][1] <= ctx->fracwell[w_no].coords[1] ))    &&
               ((coords_array[ek+1][ej][ei][2] >= ctx->fracwell[w_no].coords[2]) && (coords_array[ek][ej][ei][2] <= ctx->fracwell[w_no].coords[2] ))
               )
            {
              hx   = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
              hy   = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
              hz   = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
              ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
              hwx = (ctx->fracwell[w_no].coords[0]-coords_array[ek][ej][ei][0])/hx;
              hwy = (ctx->fracwell[w_no].coords[1]-coords_array[ek][ej][ei][1])/hy;
              hwz = (ctx->fracwell[w_no].coords[2]-coords_array[ek][ej][ei][2])/hz;
              if(ctx->fracwell[w_no].condition == RATE){
                ierr = VecApplyFractureWellFlowRate(RHS_local,&ctx->e3D,ctx->fracwell[w_no].Qw,hwx,hwy,hwz,ek,ej,ei,v_array);CHKERRQ(ierr);
                for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
                  for (j = 0; j < ctx->e3D.nphiy; j++) {
                    for (i = 0; i < ctx->e3D.nphix; i++,l++) {
                      if(ctx->fracwell[w_no].type == INJECTOR){
                        RHS_array[ek+k][ej+j][ei+i][3] -= timestepsize*RHS_local[l];
                      }
                      else if(ctx->well[w_no].type == PRODUCER)
                      {
                        RHS_array[ek+k][ej+j][ei+i][3] -= -timestepsize*RHS_local[l];
                      }
                    }
                  }
                }
              }
              else if(ctx->fracwell[w_no].condition == PRESSURE){
              }
            }
          }
        }
      }
      w_no++;
    }
*/
    
    
    
  }
  ierr = MatAssemblyBegin(Krhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Krhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatApplyKSPVelocityBC(K,Krhs,&ctx->bcQ[0]);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Krhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Krhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,source_local,&source_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&source_local);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daFlow,RHS_localVec,&RHS_array);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(ctx->daFlow,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ctx->daFlow,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daFlow,&RHS_localVec);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,prebc_local,&prebc_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&prebc_local);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScal,pressure_diff_local,&pressure_diff_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&pressure_diff_local);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_diff_local,&u_diff_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daVect,&u_diff_local);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_old_local,&u_old_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daVect,&u_old_local);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_local,&u_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daVect,&u_local);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(ctx->daScal,v_old_local,&v_old_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&v_old_local);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScal,one_local,&one_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&one_local);CHKERRQ(ierr);

  ierr = PetscFree5(KAf_local,KBf_local,KDf_local,KBfTrans_local,KL_local);CHKERRQ(ierr);
  ierr = PetscFree5(KAfrhs_local,KBfrhs_local,KDfrhs_local,KBfTransrhs_local,KLrhs_local);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(ctx->daScal,fracflow_local,&fracflow_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&fracflow_local);CHKERRQ(ierr);
  
  ierr = PetscFree5(KA_local,KB_local,KD_local,KBTrans_local,KS_local);CHKERRQ(ierr);
  ierr = PetscFree4(KArhs_local,KBrhs_local,KDrhs_local,KBTransrhs_local);CHKERRQ(ierr);
  ierr = PetscFree3(RHS_local,row,row1);CHKERRQ(ierr);
  
  ierr = VecDestroy(&U_diff);CHKERRQ(ierr);
  ierr = VecDestroy(&Pressure_diff);CHKERRQ(ierr);
  ierr = VecDestroy(&Ones);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecApplyWellFlowRate"
extern PetscErrorCode VecApplyWellFlowRate(PetscReal *RHS_local,CartFE_Element3D *e,PetscReal Q,PetscReal hwx,PetscReal hwy,PetscReal hwz,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ***v_array)
{
  PetscReal       phi[2][2][2];
  PetscInt        i,j,k,l;
  PetscReal             v_elem_pt = 0.;
  
  PetscFunctionBegin;
  phi[0][0][0] = (1.-hwz)*(1.-hwy)*(1.-hwx);
  phi[0][0][1] = (1.-hwz)*(1.-hwy)*hwx;
  phi[0][1][0] = (1.-hwz)*hwy*(1.-hwx);
  phi[0][1][1] = (1.-hwz)*hwy*hwx;
  phi[1][0][0] = hwz*(1.-hwy)*(1.-hwx);
  phi[1][0][1] = hwz*(1.-hwy)*hwx;
  phi[1][1][0] = hwz*hwy*(1.-hwx);
  phi[1][1][1] = hwz*hwy*hwx;
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        v_elem_pt += v_array[ek+k][ej+j][ei+i] * phi[k][j][i];
      }
    }
  }
  for(l = 0, k = 0; k < 2; k++){
    for(j = 0; j < 2; j++){
      for(i = 0; i < 2; i++, l++){
        RHS_local[l] = 0;
        RHS_local[l] = Q*phi[k][j][i]*(pow(v_elem_pt,2));
      }
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecApplyFractureWellFlowRate"
extern PetscErrorCode VecApplyFractureWellFlowRate(PetscReal *RHS_local,CartFE_Element3D *e,PetscReal Q,PetscReal hwx,PetscReal hwy,PetscReal hwz,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ***v_array)
{
  PetscReal       phi[2][2][2];
  PetscInt        i,j,k,l;
  PetscReal             v_elem_pt = 0.;
  
  PetscFunctionBegin;
  phi[0][0][0] = (1.-hwz)*(1.-hwy)*(1.-hwx);
  phi[0][0][1] = (1.-hwz)*(1.-hwy)*hwx;
  phi[0][1][0] = (1.-hwz)*hwy*(1.-hwx);
  phi[0][1][1] = (1.-hwz)*hwy*hwx;
  phi[1][0][0] = hwz*(1.-hwy)*(1.-hwx);
  phi[1][0][1] = hwz*(1.-hwy)*hwx;
  phi[1][1][0] = hwz*hwy*(1.-hwx);
  phi[1][1][1] = hwz*hwy*hwx;
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        v_elem_pt += v_array[ek+k][ej+j][ei+i] * phi[k][j][i];
      }
    }
  }
  for(l = 0, k = 0; k < 2; k++){
    for(j = 0; j < 2; j++){
      for(i = 0; i < 2; i++, l++){
        RHS_local[l] = 0;
        RHS_local[l] = Q*phi[k][j][i]*(1.0-(pow(v_elem_pt,2)));
      }
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecApplyRegularizedFractureWellFlowRate"
extern PetscErrorCode VecApplyRegularizedFractureWellFlowRate(PetscReal *RHS_local,CartFE_Element3D *e,PetscReal Q,PetscReal hwx,PetscReal hwy,PetscReal hwz,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ***v_array)
{
  PetscReal       phi[2][2][2];
  PetscInt        i,j,k,l;
  PetscReal             v_elem_pt = 0.;
  
  PetscFunctionBegin;
  phi[0][0][0] = (1.-hwz)*(1.-hwy)*(1.-hwx);
  phi[0][0][1] = (1.-hwz)*(1.-hwy)*hwx;
  phi[0][1][0] = (1.-hwz)*hwy*(1.-hwx);
  phi[0][1][1] = (1.-hwz)*hwy*hwx;
  phi[1][0][0] = hwz*(1.-hwy)*(1.-hwx);
  phi[1][0][1] = hwz*(1.-hwy)*hwx;
  phi[1][1][0] = hwz*hwy*(1.-hwx);
  phi[1][1][1] = hwz*hwy*hwx;
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        v_elem_pt += v_array[ek+k][ej+j][ei+i] * phi[k][j][i];
      }
    }
  }
  for(l = 0, k = 0; k < 2; k++){
    for(j = 0; j < 2; j++){
      for(i = 0; i < 2; i++, l++){
        RHS_local[l] = 0;
        RHS_local[l] = Q*phi[k][j][i]*(1.0-(pow(v_elem_pt,2)));
      }
    }
  }
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "VecApplyPressureBC"
extern PetscErrorCode VecApplyPressureBC(PetscReal *RHS_local,PetscReal ***pre_array,PetscInt ek,PetscInt ej,PetscInt ei,FACE face,CartFE_Element2D *e,VFFlowProp flowpropty,PetscReal ****perm_array,PetscReal ***v_array)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,l,g;
  PetscReal      *pre_elem,*v_elem;
  PetscReal      beta_c,mu;
  PetscReal      kx,ky,kz;
  
  PetscFunctionBegin;
  ierr = PetscMalloc2(e->ng,PetscReal,&pre_elem,e->ng,PetscReal,&v_elem);CHKERRQ(ierr);
  beta_c  = flowpropty.beta;
  mu      = flowpropty.mu;
  kx = perm_array[ek][ej][ei][0];
  ky = perm_array[ek][ej][ei][1];
  kz = perm_array[ek][ej][ei][2];
    // Initialize pre_Elem
  for (g = 0; g < e->ng; g++) {
    pre_elem[g] = 0;
    v_elem[g] = 0.;
  }
  switch (face) {
    case X0:
      for (k = 0; k < e->nphix; k++) {
        for (j = 0; j < e->nphiy; j++) {
          for (i = 0; i < e->nphiz; i++) {
            for (g = 0; g < e->ng; g++) {
              pre_elem[g] += e->phi[i][j][k][g]*pre_array[ek+k][ej+j][ei+i];
              v_elem[g] += e->phi[i][j][k][g]*v_array[ek+k][ej+j][ei+i];
            }
          }
        }
      }
      /*                      Accumulate              */
      for (l=0,k = 0; k < e->nphix; k++) {
        for (j = 0; j < e->nphiy; j++) {
          for (i = 0; i < e->nphiz; i++, l++) {
            RHS_local[l] = 0.;
            for (g = 0; g < e->ng; g++) {
              RHS_local[l] -= e->weight[g]*e->phi[i][j][k][g]*(pow(v_elem[g],2))*pre_elem[g];
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
              pre_elem[g] += e->phi[i][j][k][g]*pre_array[ek+k][ej+j][ei+1];
              v_elem[g] += e->phi[i][j][k][g]*v_array[ek+k][ej+j][ei+1];
            }
          }
        }
      }
      /*                      Accumulate              */
      for (l=0,k = 0; k < e->nphix; k++) {
        for (j = 0; j < e->nphiy; j++) {
          for (i = 0; i < e->nphiz; i++, l++) {
            RHS_local[l] = 0.;
            for (g = 0; g < e->ng; g++) {
              RHS_local[l] -= e->weight[g]*e->phi[i][j][k][g]*(pow(v_elem[g],2))*pre_elem[g];
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
              pre_elem[g] += e->phi[j][k][i][g]*pre_array[ek+k][ej+j][ei+i];
              v_elem[g] += e->phi[j][k][i][g]*v_array[ek+k][ej+j][ei+i];
            }
          }
        }
      }
      /*                      Accumulate              */
      for (l=0,k = 0; k < e->nphiy; k++) {
        for (j = 0; j < e->nphiz; j++) {
          for (i = 0; i < e->nphix; i++, l++) {
            RHS_local[l] = 0.;
            for (g = 0; g < e->ng; g++) {
              RHS_local[l] -= e->weight[g]*e->phi[j][k][i][g]*(pow(v_elem[g],2))*pre_elem[g];
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
              pre_elem[g] += e->phi[j][k][i][g]*pre_array[ek+k][ej+1][ei+i];
              v_elem[g] += e->phi[j][k][i][g]*v_array[ek+k][ej+1][ei+i];
            }
          }
        }
      }
      /*                      Accumulate              */
      for (l=0,k = 0; k < e->nphiy; k++) {
        for (j = 0; j < e->nphiz; j++) {
          for (i = 0; i < e->nphix; i++, l++) {
            RHS_local[l] = 0.;
            for (g = 0; g < e->ng; g++) {
              RHS_local[l] -= e->weight[g]*e->phi[j][k][i][g]*(pow(v_elem[g],2))*pre_elem[g];
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
              pre_elem[g] += e->phi[k][j][i][g]*pre_array[ek][ej+j][ei+i];
              v_elem[g] += e->phi[k][j][i][g]*v_array[ek][ej+j][ei+i];
            }
          }
        }
      }
      /*                      Accumulate              */
      for (l=0,k = 0; k < e->nphiz; k++) {
        for (j = 0; j < e->nphiy; j++) {
          for (i = 0; i < e->nphix; i++, l++) {
            RHS_local[l] = 0.;
            for (g = 0; g < e->ng; g++) {
              RHS_local[l] -= e->weight[g]*e->phi[k][j][i][g]*(pow(v_elem[g],2))*pre_elem[g];
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
              pre_elem[g] += e->phi[k][j][i][g]*pre_array[ek+1][ej+j][ei+i];
              v_elem[g] += e->phi[k][j][i][g]*v_array[ek+1][ej+j][ei+i];
            }
          }
        }
      }
      /*                      Accumulate              */
      for (l=0,k = 0; k < e->nphiz; k++) {
        for (j = 0; j < e->nphiy; j++) {
          for (i = 0; i < e->nphix; i++, l++) {
            RHS_local[l] = 0.;
            for (g = 0; g < e->ng; g++) {
              RHS_local[l] -= e->weight[g]*e->phi[k][j][i][g]*(pow(v_elem[g],2))*pre_elem[g];
            }
          }
        }
      }
      break;
  }
  ierr = PetscFree2(pre_elem,v_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatApplyKSPVelocityBC"
extern PetscErrorCode MatApplyKSPVelocityBC(Mat K,Mat Klhs,BC *bcQ)
{
  PetscErrorCode ierr;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       i,j,k,c;
  MatStencil    *row;
  PetscReal      one=1.;
  PetscInt       numBC=0,l=0;
  PetscInt       dim,dof=3;
  DM                              da;
  
  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject) K,"DM",(PetscObject *) &da); CHKERRQ(ierr);
  if (!da) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG," Matrix not generated from a DA");
  
  ierr = DMDAGetInfo(da,&dim,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  /*
   Compute the number of boundary nodes on each processor.
   Edges and corners are counted multiple times (2 and 3 resp)
   */
  for (c = 0; c < dof; c++){
    if (xs == 0       && bcQ[c].face[X0] == FIXED)             numBC += ym * zm;
    if (xs + xm == nx && bcQ[c].face[X1] == FIXED)             numBC += ym * zm;
    if (ys == 0       && bcQ[c].face[Y0] == FIXED)             numBC += xm * zm;
    if (ys + ym == ny && bcQ[c].face[Y1] == FIXED)             numBC += xm * zm;
    if (zs == 0       && bcQ[c].face[Z0] == FIXED && dim == 3) numBC += xm * ym;
    if (zs + zm == nz && bcQ[c].face[Z1] == FIXED && dim == 3) numBC += xm * ym;
    if (xs == 0       && ys == 0       && zs == 0       && bcQ[c].vertex[X0Y0Z0] == FIXED) numBC++;
    if (xs == 0       && ys + ym == ny && zs == 0       && bcQ[c].vertex[X0Y1Z0] == FIXED) numBC++;
    if (xs + xm == nx && ys == 0       && zs == 0       && bcQ[c].vertex[X1Y0Z0] == FIXED) numBC++;
    if (xs + xm == nx && ys + ym == ny && zs == 0       && bcQ[c].vertex[X1Y1Z0] == FIXED) numBC++;
    if (xs == 0       && ys == 0       && zs + zm == nz && bcQ[c].vertex[X0Y0Z1] == FIXED && dim == 3) numBC++;
    if (xs == 0       && ys + ym == ny && zs + zm == nz && bcQ[c].vertex[X0Y1Z1] == FIXED && dim == 3) numBC++;
    if (xs + xm == nx && ys == 0       && zs + zm == nz && bcQ[c].vertex[X1Y0Z1] == FIXED && dim == 3) numBC++;
    if (xs + xm == nx && ys + ym == ny && zs + zm == nz && bcQ[c].vertex[X1Y1Z1] == FIXED && dim == 3) numBC++;
  }
  ierr = PetscMalloc(numBC * sizeof(MatStencil),&row);CHKERRQ(ierr);
  /*
   Create an array of rows to be zeroed out
   */
  /*
   i == 0
   */
  for (c = 0; c < dof; c++) {
    if (xs == 0 && bcQ[c].face[X0] == FIXED) {
      for (k = zs; k < zs + zm; k++) {
        for (j = ys; j < ys + ym; j++) {
          row[l].i = 0; row[l].j = j; row[l].k = k; row[l].c = c;
          l++;
        }
      }
    }
    /*
     i == nx-1
     */
    if (xs + xm == nx && bcQ[c].face[X1] == FIXED) {
      for (k = zs; k < zs + zm; k++) {
        for (j = ys; j < ys + ym; j++) {
          row[l].i = nx-1; row[l].j = j; row[l].k = k; row[l].c = c;
          l++;
        }
      }
    }
    /*
     y == 0
     */
    if (ys == 0 && bcQ[c].face[Y0] == FIXED) {
      for (k = zs; k < zs + zm; k++) {
        for (i = xs; i < xs + xm; i++) {
          row[l].i = i; row[l].j = 0; row[l].k = k; row[l].c = c;
          l++;
        }
      }
    }
    /*
     y == ny-1
     */
    if (ys + ym == ny && bcQ[c].face[Y1] == FIXED) {
      for (k = zs; k < zs + zm; k++) {
        for (i = xs; i < xs + xm; i++) {
          row[l].i = i; row[l].j = ny-1; row[l].k = k; row[l].c = c;
          l++;
        }
      }
    }
    if (dim==3){
      /*
       z == 0
       */
      if (zs == 0 && bcQ[c].face[Z0] == FIXED) {
        for (j = ys; j < ys + ym; j++) {
          for (i = xs; i < xs + xm; i++) {
            row[l].i = i; row[l].j = j; row[l].k = 0; row[l].c = c;
            l++;
          }
        }
      }
      /*
       z == nz-1
       */
      if (zs + zm == nz && bcQ[c].face[Z1] == FIXED) {
        for (j = ys; j < ys + ym; j++) {
          for (i = xs; i < xs + xm; i++) {
            row[l].i = i; row[l].j = j; row[l].k = nz-1; row[l].c = c;
            l++;
          }
        }
      }
    }
    if (xs == 0       && ys == 0       && zs == 0       && bcQ[c].vertex[X0Y0Z0] == FIXED) {
      row[l].i = 0; row[l].j = 0; row[l].k = 0; row[l].c = c;
      l++;
    }
    if (xs == 0       && ys == 0       && zs + zm == nz && bcQ[c].vertex[X0Y0Z1] == FIXED && dim ==3) {
      row[l].i = 0; row[l].j = 0; row[l].k = nz-1; row[l].c = c;
      l++;
    }
    if (xs == 0       && ys + ym == ny && zs == 0       && bcQ[c].vertex[X0Y1Z0] == FIXED) {
      row[l].i = 0; row[l].j = ny-1; row[l].k = 0; row[l].c = c;
      l++;
    }
    if (xs == 0       && ys + ym == ny && zs + zm == nz && bcQ[c].vertex[X0Y1Z1] == FIXED && dim ==3) {
      row[l].i = 0; row[l].j = ny-1; row[l].k = nz-1; row[l].c = c;
      l++;
    }
    if (xs + xm == nx && ys == 0       && zs == 0       && bcQ[c].vertex[X1Y0Z0] == FIXED) {
      row[l].i = nx-1; row[l].j = 0; row[l].k = 0; row[l].c = c;
      l++;
    }
    if (xs + xm == nx && ys == 0       && zs + zm == nz && bcQ[c].vertex[X1Y0Z1] == FIXED && dim ==3) {
      row[l].i = nx-1; row[l].j = 0; row[l].k = nz-1; row[l].c = c;
      l++;
    }
    if (xs + xm == nx && ys + ym == ny && zs == 0       && bcQ[c].vertex[X1Y1Z0] == FIXED) {
      row[l].i = nx-1; row[l].j = ny-1; row[l].k = 0; row[l].c = c;
      l++;
    }
    if (xs + xm == nx && ys + ym == ny && zs + zm == nz && bcQ[c].vertex[X1Y1Z1] == FIXED && dim ==3) {
      row[l].i = nx=1; row[l].j = ny-1; row[l].k = nz-1; row[l].c = c;
      l++;
    }
    
  }
  ierr = MatZeroRowsStencil(K,numBC,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
    //      ierr = MatZeroRowsStencil(Klhs,numBC,row,zero,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscFree(row);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Flow_Vecg"
extern PetscErrorCode Flow_Vecg(PetscReal *Kg_local,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,VFFlowProp flowpropty,PetscReal ****perm_array,PetscReal ***v_array)
{
  PetscErrorCode    ierr;
  PetscInt          i,j,k,l;
  PetscInt          eg;
  PetscReal         beta_c,mu,rho,gamma_c,gx,gy,gz;
  PetscReal         kx,ky,kz,kxy,kxz,kyz;
  PetscReal         *v_elem;
  
  PetscFunctionBegin;
  ierr = PetscMalloc(e->ng*sizeof(PetscReal),&v_elem);CHKERRQ(ierr);
  beta_c  = flowpropty.beta;
  gamma_c = flowpropty.gamma;
  rho     = flowpropty.rho;
  mu      = flowpropty.mu;
  gx      = flowpropty.g[0];
  gy      = flowpropty.g[1];
  gz      = flowpropty.g[2];
  kx  = perm_array[ek][ej][ei][0];
  ky  = perm_array[ek][ej][ei][1];
  kz  = perm_array[ek][ej][ei][2];
  kxy = perm_array[ek][ej][ei][3];
  kxz = perm_array[ek][ej][ei][4];
  kyz = perm_array[ek][ej][ei][5];
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
      for (i = 0; i < e->nphix; i++,l++) {
        Kg_local[l] = 0.;
        for (eg = 0; eg < e->ng; eg++) {
          Kg_local[l] += -0.5*rho*gamma_c*beta_c/mu*((kx*gx+kxy*gy+kxz*gz)*e->dphi[k][j][i][0][eg]
                                                     +(kxy*gx+ky*gy+kyz*gz)*e->dphi[k][j][i][1][eg]
                                                     +(kxz*gx+kyz*gy+kz*gz)*e->dphi[k][j][i][2][eg])*(pow(v_elem[eg],2))*e->weight[eg];
        }
      }
    }
  }
  ierr = PetscFree(v_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Flow_Vecf"
extern PetscErrorCode Flow_Vecf(PetscReal *Kf_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscInt c,VFFlowProp flowpropty,PetscReal ***v_array)
{
  PetscErrorCode    ierr;
  PetscInt          i,j,k,l;
  PetscInt          eg;
  PetscReal         rho,gamma_c;
  PetscReal         grav;
  PetscReal         *v_elem;
  
  PetscFunctionBegin;
  ierr = PetscMalloc(e->ng*sizeof(PetscReal),&v_elem);CHKERRQ(ierr);
  gamma_c = flowpropty.gamma;
  rho     = flowpropty.rho;
  grav      = flowpropty.g[c];
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
      for (i = 0; i < e->nphix; i++,l++) {
        Kf_ele[l] = 0.;
        for (eg = 0; eg < e->ng; eg++) {
          Kf_ele[l] += 0.5*grav*gamma_c*rho*e->phi[k][j][i][eg]*(pow(v_elem[eg],2))*e->weight[eg];
        }
      }
    }
  }
  ierr = PetscFree(v_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "Flow_MatD"
extern PetscErrorCode Flow_MatD(PetscReal *Kd_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,VFFlowProp flowpropty,PetscReal ****perm_array,PetscReal ***v_array)
{
  PetscErrorCode    ierr;
  PetscInt          i,j,k,l;
  PetscInt          ii,jj,kk;
  PetscInt          eg;
  PetscReal         beta_c,mu;
  PetscReal         kx,ky,kz,kxy,kxz,kyz;
  PetscReal         *v_elem;
  
  PetscFunctionBegin;
  ierr = PetscMalloc(e->ng*sizeof(PetscReal),&v_elem);CHKERRQ(ierr);
  beta_c = flowpropty.beta;
  mu     = flowpropty.mu;
  kx     = perm_array[ek][ej][ei][0];
  ky     = perm_array[ek][ej][ei][1];
  kz     = perm_array[ek][ej][ei][2];
  kxy    = perm_array[ek][ej][ei][3];
  kxz    = perm_array[ek][ej][ei][4];
  kyz    = perm_array[ek][ej][ei][5];
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
                Kd_ele[l] += -0.5*beta_c/mu*
                ((kx*e->dphi[k][j][i][0][eg]+kxy*e->dphi[k][j][i][1][eg]+kxz*e->dphi[k][j][i][2][eg])*e->dphi[kk][jj][ii][0][eg]
                 +(kxy*e->dphi[k][j][i][0][eg]+ky*e->dphi[k][j][i][1][eg]+kyz*e->dphi[k][j][i][2][eg])*e->dphi[kk][jj][ii][1][eg]
                 +(kxz*e->dphi[k][j][i][0][eg]+kyz*e->dphi[k][j][i][1][eg]+kz*e->dphi[k][j][i][2][eg])*e->dphi[kk][jj][ii][2][eg])*(pow(v_elem[eg],2))*e->weight[eg];
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
#define __FUNCT__ "Flow_MatBTranspose"
extern PetscErrorCode Flow_MatBTranspose(PetscReal *KB_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscInt c,PetscReal ***v_array)
{
  PetscErrorCode    ierr;
  PetscInt          i,j,k,l;
  PetscInt          ii,jj,kk;
  PetscInt          eg;
  PetscReal         *v_elem;
  
  PetscFunctionBegin;
  ierr = PetscMalloc(e->ng*sizeof(PetscReal),&v_elem);CHKERRQ(ierr);
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
              KB_ele[l] = 0.;
              for (eg = 0; eg < e->ng; eg++) {
                KB_ele[l] += -(e->phi[k][j][i][eg]*e->dphi[kk][jj][ii][c][eg]
                               +0.5*e->dphi[k][j][i][c][eg]*e->phi[kk][jj][ii][eg])*(pow(v_elem[eg],2))*e->weight[eg];
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
#define __FUNCT__ "Flow_MatB"
extern PetscErrorCode Flow_MatB(PetscReal *KB_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscInt c,PetscReal ***v_array)
{
  PetscErrorCode    ierr;
  PetscInt          i,j,k,l;
  PetscInt          ii,jj,kk;
  PetscInt          eg;
  PetscReal         *v_elem;
  
  PetscFunctionBegin;
  ierr = PetscMalloc(e->ng*sizeof(PetscReal),&v_elem);CHKERRQ(ierr);
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
              KB_ele[l] = 0.;
              for (eg = 0; eg < e->ng; eg++) {
                KB_ele[l] += -(e->phi[kk][jj][ii][eg]*e->dphi[k][j][i][c][eg]
                               +0.5*e->dphi[kk][jj][ii][c][eg]*e->phi[k][j][i][eg])*(pow(v_elem[eg],2))*e->weight[eg];
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
#define __FUNCT__ "Flow_MatA"
extern PetscErrorCode Flow_MatA(PetscReal *A_local,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscInt c,VFFlowProp flowpropty,PetscReal ****perm_array,PetscReal ***v_array)
{
  PetscErrorCode    ierr;
  PetscInt          i,j,k,l;
  PetscInt          ii,jj,kk;
  PetscInt          eg;
  PetscReal         *v_elem;
  PetscReal         perm_inv = 0;
  PetscReal         beta_c,mu;
  
  PetscFunctionBegin;
  ierr = PetscMalloc(e->ng*sizeof(PetscReal),&v_elem);CHKERRQ(ierr);
  beta_c = flowpropty.beta;
  mu     = flowpropty.mu;
  perm_inv     = 1./perm_array[ek][ej][ei][c];
  if((ierr = PetscIsInfOrNanScalar(perm_inv)))
  {
    perm_inv = 0;
  }
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
              A_local[l] = 0;
              for (eg = 0; eg < e->ng; eg++) {
                A_local[l] += 0.5*mu*perm_inv/(beta_c)*e->phi[k][j][i][eg]*e->phi[kk][jj][ii][eg]*(pow(v_elem[eg],2))*e->weight[eg];
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
#define __FUNCT__ "VF_MatA_local"
extern PetscErrorCode VF_MatA_local(PetscReal *A_local,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ***v_array)
{
  PetscErrorCode    ierr;
  PetscInt          i,j,k,l;
  PetscInt          ii,jj,kk;
  PetscInt          eg;
  PetscReal         *v_elem;
  
  PetscFunctionBegin;
  ierr = PetscMalloc(e->ng*sizeof(PetscReal),&v_elem);CHKERRQ(ierr);
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
              A_local[l] = 0;
              for (eg = 0; eg < e->ng; eg++) {
                A_local[l] += e->phi[k][j][i][eg]*e->phi[kk][jj][ii][eg]*(pow(v_elem[eg],2))*e->weight[eg];
//                A_local[l] += e->phi[k][j][i][eg]*e->phi[kk][jj][ii][eg]*e->weight[eg];
              }
            }
          }
        }
      }
    }
  }
  /*
  for (l = 0,k = 0; k < e->nphiz; k++) {
		for (j = 0; j < e->nphiy; j++) {
			for (i = 0; i < e->nphix; i++) {
				for (kk = 0; kk < e->nphiz; kk++) {
					for (jj = 0; jj < e->nphiy; jj++) {
						for (ii = 0; ii < e->nphix; ii++,l++) {
              if(v_array[ek+k][ej+j][ei+i] <= 0.05 || v_array[ek+kk][ej+jj][ei+ii] <= 0.05 ){
                A_local[l] = 0.;
              }
            }
          }
        }
      }
    }
  }
   */
  ierr = PetscFree(v_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_RHSFlowMechUCoupling_local"
extern PetscErrorCode VF_RHSFlowMechUCoupling_local(PetscReal *K_local,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,VFFlowProp flowpropty,PetscReal ****u_diff_array,PetscReal ***v_array)
{
  PetscErrorCode    ierr;
  PetscInt          i,j,k,l,c;
  PetscInt          eg;
  PetscReal         *v_elem;
  PetscReal         *du_elem[3],alphabiot;
  
  PetscFunctionBegin;
  alphabiot  = flowpropty.alphabiot;
  ierr = PetscMalloc4(e->ng,PetscReal,&v_elem,
                      e->ng,PetscReal,&du_elem[0],
                      e->ng,PetscReal,&du_elem[1],
                      e->ng,PetscReal,&du_elem[2]);CHKERRQ(ierr);
  for (eg = 0; eg < e->ng; eg++){
    v_elem[eg] = 0.;
    for (c = 0; c < 3; c++){
      du_elem[c][eg] = 0.;
    }
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          v_elem[eg] += v_array[ek+k][ej+j][ei+i] * e->phi[k][j][i][eg];
          for (c=0; c<3; c++){
            du_elem[c][eg] += u_diff_array[ek+k][ej+j][ei+i][c]*e->dphi[k][j][i][c][eg];
          }
        }
      }
    }
  }
  for (l = 0,k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++,l++) {
        K_local[l] = 0.;
        for (eg = 0; eg < e->ng; eg++) {
          for(c = 0; c < 3; c++){
            K_local[l] += alphabiot*du_elem[c][eg]*e->phi[k][j][i][eg]*(pow(v_elem[eg],2))*e->weight[eg];
          }
        }
      }
    }
  }
  ierr = PetscFree4(v_elem,du_elem[0],du_elem[1],du_elem[2]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "VF_RHSFlowMechUCouplingFIXSTRESS_local"
extern PetscErrorCode VF_RHSFlowMechUCouplingFIXSTRESS_local(PetscReal *K_local,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,VFFlowProp flowpropty,VFMatProp *matprop,PetscReal ***pressure_diff_array,PetscReal ***v_array)
{
  PetscErrorCode    ierr;
  PetscInt          i,j,k,l,c;
  PetscInt          eg;
  PetscReal         *v_elem;
  PetscReal         *p_diff_elem,alphabiot,K_dr,nu,E;
  
  PetscFunctionBegin;
  alphabiot  = flowpropty.alphabiot;
  K_dr  = flowpropty.K_dr;
  nu       = matprop->nu;
  E       = matprop->E;
  ierr = PetscMalloc2(e->ng,PetscReal,&v_elem,e->ng,PetscReal,&p_diff_elem);CHKERRQ(ierr);
  for (eg = 0; eg < e->ng; eg++){
    v_elem[eg] = 0.;
    p_diff_elem[eg] = 0.;
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          v_elem[eg] += v_array[ek+k][ej+j][ei+i] * e->phi[k][j][i][eg];
          p_diff_elem[eg] += pressure_diff_array[ek+k][ej+j][ei+i] * e->phi[k][j][i][eg];
        }
      }
    }
  }
  for (l = 0,k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++,l++) {
        K_local[l] = 0.;
        for (eg = 0; eg < e->ng; eg++) {
          K_local[l] += alphabiot*alphabiot/K_dr*p_diff_elem[eg]*e->phi[k][j][i][eg]*(pow(v_elem[eg],2))*e->weight[eg];
        }
      }
    }
  }
  
  ierr = PetscFree2(v_elem,p_diff_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}













#undef __FUNCT__
#define __FUNCT__ "VF_RHSFractureFlowCoupling_local"
extern PetscErrorCode VF_RHSFractureFlowCoupling_local(PetscReal *K_local,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ****u_array,PetscReal ***v_array,PetscReal ****u_old_array,PetscReal ***v_old_array)
{
  PetscErrorCode    ierr;
  PetscInt          i,j,k,l,c;
  PetscInt          eg;
  PetscReal         *u_elem[3],*dv_elem[3];
  PetscReal         *u_old_elem[3],*dv_old_elem[3];
  PetscReal         *n_elem[3],*v_mag_elem;
  PetscReal         *n_old_elem[3],*v_old_mag_elem;
  
  PetscFunctionBegin;
  ierr = PetscMalloc6(e->ng,PetscReal,&dv_elem[0],
                      e->ng,PetscReal,&dv_elem[1],
                      e->ng,PetscReal,&dv_elem[2],
                      e->ng,PetscReal,&u_elem[0],
                      e->ng,PetscReal,&u_elem[1],
                      e->ng,PetscReal,&u_elem[2]);CHKERRQ(ierr);

  ierr = PetscMalloc6(e->ng,PetscReal,&dv_old_elem[0],
                      e->ng,PetscReal,&dv_old_elem[1],
                      e->ng,PetscReal,&dv_old_elem[2],
                      e->ng,PetscReal,&u_old_elem[0],
                      e->ng,PetscReal,&u_old_elem[1],
                      e->ng,PetscReal,&u_old_elem[2]);CHKERRQ(ierr);
  ierr = PetscMalloc4(e->ng,PetscReal,&n_elem[0],e->ng,PetscReal,&n_elem[1],e->ng,PetscReal,&n_elem[2],e->ng,PetscReal,&v_mag_elem);CHKERRQ(ierr);
  ierr = PetscMalloc4(e->ng,PetscReal,&n_old_elem[0],e->ng,PetscReal,&n_old_elem[1],e->ng,PetscReal,&n_old_elem[2],e->ng,PetscReal,&v_old_mag_elem);CHKERRQ(ierr);

  for (eg = 0; eg < e->ng; eg++){
    for (c = 0; c < 3; c++){
      dv_elem[c][eg] = 0.;
      dv_old_elem[c][eg] = 0.;
      u_elem[c][eg] = 0.;
      u_old_elem[c][eg] = 0.;
      n_elem[c][eg] = 0;
      n_old_elem[c][eg] = 0;
    }
    v_mag_elem[eg] = 0.;
    v_old_mag_elem[eg] = 0.;
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          for (c= 0; c < 3; c++){
            dv_elem[c][eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
            dv_old_elem[c][eg] += v_old_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
            u_elem[c][eg] += u_array[ek+k][ej+j][ei+i][c]*e->phi[k][j][i][eg];
            u_old_elem[c][eg] += u_old_array[ek+k][ej+j][ei+i][c]*e->phi[k][j][i][eg];
          }
        }
      }
    }
  }
  for (eg = 0; eg < e->ng; eg++){
    v_mag_elem[eg] = sqrt((pow(dv_elem[0][eg],2))+(pow(dv_elem[1][eg],2))+(pow(dv_elem[2][eg],2)));
    v_old_mag_elem[eg] = sqrt((pow(dv_old_elem[0][eg],2))+(pow(dv_old_elem[1][eg],2))+(pow(dv_old_elem[2][eg],2)));
    for(c = 0; c < 3; c++){
      n_elem[c][eg] = dv_elem[c][eg]/v_mag_elem[eg];
      n_old_elem[c][eg] = dv_old_elem[c][eg]/v_old_mag_elem[eg];
    }
    if((ierr = PetscIsInfOrNanScalar(n_elem[0][eg])) || (ierr = PetscIsInfOrNanScalar(n_elem[1][eg])) || (ierr = PetscIsInfOrNanScalar(n_elem[2][eg])) )
    {
      n_elem[0][eg] = n_elem[1][eg] = n_elem[2][eg] = v_mag_elem[eg] = 0;
    }
    if((ierr = PetscIsInfOrNanScalar(n_old_elem[0][eg])) || (ierr = PetscIsInfOrNanScalar(n_old_elem[1][eg])) || (ierr = PetscIsInfOrNanScalar(n_old_elem[2][eg])) )
    {
      n_old_elem[0][eg] = n_old_elem[1][eg] = n_old_elem[2][eg] = v_old_mag_elem[eg] = 0;
    }
  }
  for (l = 0,k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++,l++) {
        K_local[l] = 0.;
        for (eg = 0; eg < e->ng; eg++) {
          for(c = 0; c < 3; c++){
            K_local[l] += (dv_elem[c][eg]*u_elem[c][eg]-dv_old_elem[c][eg]*u_old_elem[c][eg])*e->phi[k][j][i][eg]*e->weight[eg];
          }
        }
      }
    }
  }
  
  
  for (l = 0,k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++,l++) {
        if(v_array[ek+k][ej+j][ei+i] > 0.05 ){
          K_local[l] = 0.;
        }
      }
    }
  }
  
  
  ierr = PetscFree6(dv_elem[0],dv_elem[1],dv_elem[2],u_elem[0],u_elem[1],u_elem[2]);CHKERRQ(ierr);
  ierr = PetscFree6(dv_old_elem[0],dv_old_elem[1],dv_old_elem[2],u_old_elem[0],u_old_elem[1],u_old_elem[2]);CHKERRQ(ierr);
  ierr = PetscFree4(n_elem[0],n_elem[1],n_elem[2],v_mag_elem);CHKERRQ(ierr);
  ierr = PetscFree4(n_old_elem[0],n_old_elem[1],n_old_elem[2],v_old_mag_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}















#undef __FUNCT__
#define __FUNCT__ "VF_MatAFractureFlowCoupling_local"
extern PetscErrorCode VF_MatAFractureFlowCoupling_local(PetscReal *Kd_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ****u_array,PetscReal ***v_array)
{
  PetscErrorCode          ierr;
  PetscInt                i,j,k,l,c;
  PetscInt                ii,jj,kk;
  PetscReal               *dv_elem[3],*u_elem[3];
  PetscInt                eg;
  
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
							Kd_ele[l] = 0.;
							for (eg = 0; eg < e->ng; eg++) {
                for(c = 0; c < 3; c++){
                  Kd_ele[l] += e->phi[k][j][i][eg]*e->phi[kk][jj][ii][eg]*pow((dv_elem[c][eg]*u_elem[c][eg]),1)*e->weight[eg];
                }
                
              }
						}
					}
				}
			}
		}
	}
  ierr = PetscFree6(dv_elem[0],dv_elem[1],dv_elem[2],u_elem[0],u_elem[1],u_elem[2]);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}





#undef __FUNCT__
#define __FUNCT__ "VF_MatLeakOff_local"
extern PetscErrorCode VF_MatLeakOff_local(PetscReal *Kd_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscInt dof,PetscReal ***v_array)
{
  PetscErrorCode          ierr;
  PetscInt                i,j,k,l,c;
  PetscInt                ii,jj,kk;
  PetscReal               *dv_elem;
  PetscInt                eg;
  
  ierr = PetscMalloc(e->ng*sizeof(PetscReal),&dv_elem);CHKERRQ(ierr);
  for (eg = 0; eg < e->ng; eg++){
    dv_elem[eg] = 0;
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          dv_elem[eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][dof][eg];
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
                  Kd_ele[l] += -e->phi[kk][jj][ii][eg]*dv_elem[eg]*e->phi[k][j][i][eg]*e->weight[eg];
                }
						}
					}
				}
			}
		}
	}
  ierr = PetscFree(dv_elem);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VF_MatBTFractureFlowCoupling_local"
extern PetscErrorCode VF_MatBTFractureFlowCoupling_local(PetscReal *Kd_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscInt dof,PetscReal ****u_array,PetscReal ***v_array)
{
  PetscErrorCode          ierr;
  PetscInt                i,j,k,l,c;
  PetscInt                ii,jj,kk;
  PetscReal               *dv_elem[3],*u_elem[3],*n_elem[3],*v_mag_elem;
  PetscInt                eg;
  
  ierr = PetscMalloc6(e->ng,PetscReal,&dv_elem[0],e->ng,PetscReal,&dv_elem[1],e->ng,PetscReal,&dv_elem[2],e->ng,PetscReal,&u_elem[0],e->ng,PetscReal,&u_elem[1],e->ng,PetscReal,&u_elem[2]);CHKERRQ(ierr);
  ierr = PetscMalloc4(e->ng,PetscReal,&n_elem[0],e->ng,PetscReal,&n_elem[1],e->ng,PetscReal,&n_elem[2],e->ng,PetscReal,&v_mag_elem);CHKERRQ(ierr);
  for (eg = 0; eg < e->ng; eg++){
    for(c = 0; c < 3; c++){
      dv_elem[c][eg] = 0;
      n_elem[c][eg] = 0;
      u_elem[c][eg] = 0;
    }
    v_mag_elem[eg] = 0.;
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
    v_mag_elem[eg] = sqrt((pow(dv_elem[0][eg],2))+(pow(dv_elem[1][eg],2))+(pow(dv_elem[2][eg],2)));
    for(c = 0; c < 3; c++){
      n_elem[c][eg] = dv_elem[c][eg]/v_mag_elem[eg];
    }
    if((ierr = PetscIsInfOrNanScalar(n_elem[0][eg])) || (ierr = PetscIsInfOrNanScalar(n_elem[1][eg])) || (ierr = PetscIsInfOrNanScalar(n_elem[2][eg])) )
    {
      n_elem[0][eg] = n_elem[1][eg] = n_elem[2][eg] = v_mag_elem[eg] = 0;
    }
  }
  for (l = 0,k = 0; k < e->nphiz; k++) {
		for (j = 0; j < e->nphiy; j++) {
			for (i = 0; i < e->nphix; i++) {
				for (kk = 0; kk < e->nphiz; kk++) {
					for (jj = 0; jj < e->nphiy; jj++) {
						for (ii = 0; ii < e->nphix; ii++,l++) {
							Kd_ele[l] = 0.;
              
              if (dof == 0){
                for (eg = 0; eg < e->ng; eg++) {
                    Kd_ele[l] += ((1.-(pow(n_elem[0][eg],2)))*e->dphi[k][j][i][0][eg]
                                  -e->dphi[k][j][i][1][eg]*n_elem[0][eg]*n_elem[1][eg]
                                  -e->dphi[k][j][i][2][eg]*n_elem[0][eg]*n_elem[2][eg])*e->phi[kk][jj][ii][eg]*(u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg])*v_mag_elem[eg]*e->weight[eg];
                }
              }
              if (dof == 1){
                for (eg = 0; eg < e->ng; eg++) { 
                    Kd_ele[l] += (-e->dphi[k][j][i][0][eg]*n_elem[0][eg]*n_elem[1][eg]
                                  +(1.-pow(n_elem[1][eg],2))*e->dphi[k][j][i][1][eg]
                                  -e->dphi[k][j][i][2][eg]*n_elem[1][eg]*n_elem[2][eg])*e->phi[kk][jj][ii][eg]*(u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg])*v_mag_elem[eg]*e->weight[eg];
                }
              }
              if (dof == 2){
                for (eg = 0; eg < e->ng; eg++) {
                    Kd_ele[l] += (-e->dphi[k][j][i][0][eg]*n_elem[0][eg]*n_elem[2][eg]
                                -e->dphi[k][j][i][1][eg]*n_elem[1][eg]*n_elem[2][eg]
                                +(1.-pow(n_elem[2][eg],2))*e->dphi[k][j][i][2][eg])*e->phi[kk][jj][ii][eg]*(u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg])*v_mag_elem[eg]*e->weight[eg];
                }
              }
						}
					}
				}
			}
		}
	}
  ierr = PetscFree6(dv_elem[0],dv_elem[1],dv_elem[2],u_elem[0],u_elem[1],u_elem[2]);CHKERRQ(ierr);
  ierr = PetscFree4(n_elem[0],n_elem[1],n_elem[2],v_mag_elem);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}






#undef __FUNCT__
#define __FUNCT__ "VF_MatBFractureFlowCoupling_local"
extern PetscErrorCode VF_MatBFractureFlowCoupling_local(PetscReal *Kd_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscInt dof,PetscReal ****u_array,PetscReal ***v_array)
{
  PetscErrorCode          ierr;
  PetscInt                i,j,k,l,c;
  PetscInt                ii,jj,kk;
  PetscReal               *dv_elem[3],*u_elem[3],*n_elem[3],*v_mag_elem;
  PetscInt                eg;
  
  ierr = PetscMalloc6(e->ng,PetscReal,&dv_elem[0],e->ng,PetscReal,&dv_elem[1],e->ng,PetscReal,&dv_elem[2],e->ng,PetscReal,&u_elem[0],e->ng,PetscReal,&u_elem[1],e->ng,PetscReal,&u_elem[2]);CHKERRQ(ierr);
  ierr = PetscMalloc4(e->ng,PetscReal,&n_elem[0],e->ng,PetscReal,&n_elem[1],e->ng,PetscReal,&n_elem[2],e->ng,PetscReal,&v_mag_elem);CHKERRQ(ierr);
  for (eg = 0; eg < e->ng; eg++){
    for(c = 0; c < 3; c++){
      dv_elem[c][eg] = 0;
      n_elem[c][eg] = 0;
      u_elem[c][eg] = 0;
    }
    v_mag_elem[eg] = 0.;
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
    v_mag_elem[eg] = sqrt((pow(dv_elem[0][eg],2))+(pow(dv_elem[1][eg],2))+(pow(dv_elem[2][eg],2)));
    for(c = 0; c < 3; c++){
      n_elem[c][eg] = dv_elem[c][eg]/v_mag_elem[eg];
    }
    if((PetscIsInfOrNanScalar(n_elem[0][eg])) || (PetscIsInfOrNanScalar(n_elem[1][eg])) || (PetscIsInfOrNanScalar(n_elem[2][eg])) )
    {
      n_elem[0][eg] = n_elem[1][eg] = n_elem[2][eg] = v_mag_elem[eg] = 0;
    }
  }
  for (l = 0,k = 0; k < e->nphiz; k++) {
		for (j = 0; j < e->nphiy; j++) {
			for (i = 0; i < e->nphix; i++) {
				for (kk = 0; kk < e->nphiz; kk++) {
					for (jj = 0; jj < e->nphiy; jj++) {
						for (ii = 0; ii < e->nphix; ii++,l++) {
							Kd_ele[l] = 0.;
              if (dof == 0){
                for (eg = 0; eg < e->ng; eg++) {
                    Kd_ele[l] += ((1.-(pow(n_elem[0][eg],2)))*e->dphi[kk][jj][ii][0][eg]
                                  -e->dphi[kk][jj][ii][1][eg]*n_elem[0][eg]*n_elem[1][eg]
                                  -e->dphi[kk][jj][ii][2][eg]*n_elem[0][eg]*n_elem[2][eg])*e->phi[k][j][i][eg]*4*(pow((u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg]),3))*v_mag_elem[eg]*e->weight[eg];
                  }
                }
              if (dof == 1){
                for (eg = 0; eg < e->ng; eg++) {
                    Kd_ele[l] += (-e->dphi[kk][jj][ii][0][eg]*n_elem[0][eg]*n_elem[1][eg]
                              +(1.-pow(n_elem[1][eg],2))*e->dphi[kk][jj][ii][1][eg]
                              -e->dphi[kk][jj][ii][2][eg]*n_elem[1][eg]*n_elem[2][eg])*e->phi[k][j][i][eg]*4*(pow((u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg]),3))*v_mag_elem[eg]*e->weight[eg];
                  }
                }
        
              if (dof == 2){
                for (eg = 0; eg < e->ng; eg++) {
                    Kd_ele[l] += (-e->dphi[kk][jj][ii][0][eg]*n_elem[0][eg]*n_elem[2][eg]
                                  -e->dphi[kk][jj][ii][1][eg]*n_elem[1][eg]*n_elem[2][eg]
                                  +(1.-pow(n_elem[2][eg],2))*e->dphi[kk][jj][ii][2][eg])*e->phi[k][j][i][eg]*4*(pow((u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg]),3))*v_mag_elem[eg]*e->weight[eg];
                }
              }
						}
					}
				}
			}
    }
  }
  ierr = PetscFree6(dv_elem[0],dv_elem[1],dv_elem[2],u_elem[0],u_elem[1],u_elem[2]);CHKERRQ(ierr);
  ierr = PetscFree4(n_elem[0],n_elem[1],n_elem[2],v_mag_elem);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VF_MatDFractureFlowCoupling_local"
extern PetscErrorCode VF_MatDFractureFlowCoupling_local(PetscReal *Kd_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ****u_array,PetscReal ***v_array)
{
  PetscErrorCode          ierr;
  PetscInt                i,j,k,l,c;
  PetscInt                ii,jj,kk;
  PetscReal               *dv_elem[3],*u_elem[3],*n_elem[3],*v_mag_elem,*v_elem;
  PetscInt                eg;
  
  ierr = PetscMalloc6(e->ng,PetscReal,&dv_elem[0],e->ng,PetscReal,&dv_elem[1],e->ng,PetscReal,&dv_elem[2],e->ng,PetscReal,&u_elem[0],e->ng,PetscReal,&u_elem[1],e->ng,PetscReal,&u_elem[2]);CHKERRQ(ierr);
  ierr = PetscMalloc5(e->ng,PetscReal,&n_elem[0],e->ng,PetscReal,&n_elem[1],e->ng,PetscReal,&n_elem[2],e->ng,PetscReal,&v_mag_elem,e->ng,PetscReal,&v_elem);CHKERRQ(ierr);
  for (eg = 0; eg < e->ng; eg++){
    for(c = 0; c < 3; c++){
      dv_elem[c][eg] = 0;
      u_elem[c][eg] = 0;
      n_elem[c][eg] = 0;
    }
    v_mag_elem[eg] = 0.;
    v_elem[eg] = 0.;
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (eg = 0; eg < e->ng; eg++) {
          v_elem[eg] += v_array[ek+k][ej+j][ei+i] * e->phi[k][j][i][eg];
          for (c = 0; c < 3; c++){
            dv_elem[c][eg] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][c][eg];
            u_elem[c][eg] += u_array[ek+k][ej+j][ei+i][c] * e->phi[k][j][i][eg];
          }
        }
      }
    }
  }
  for (eg = 0; eg < e->ng; eg++){
    v_mag_elem[eg] = sqrt((pow(dv_elem[0][eg],2))+(pow(dv_elem[1][eg],2))+(pow(dv_elem[2][eg],2)));
    for(c = 0; c < 3; c++){
      n_elem[c][eg] = dv_elem[c][eg]/v_mag_elem[eg];
    }
    if((ierr = PetscIsInfOrNanScalar(n_elem[0][eg])) || (ierr = PetscIsInfOrNanScalar(n_elem[1][eg])) || (ierr = PetscIsInfOrNanScalar(n_elem[2][eg])) )
    {
      n_elem[0][eg] = n_elem[1][eg] = n_elem[2][eg] = v_mag_elem[eg] = 0;
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
                 
                Kd_ele[l] += ((1.-(pow(n_elem[0][eg],2)))*e->dphi[k][j][i][0][eg]
                              -e->dphi[k][j][i][1][eg]*n_elem[0][eg]*n_elem[1][eg]
                              -e->dphi[k][j][i][2][eg]*n_elem[0][eg]*n_elem[2][eg])*
                ((1.-(pow(n_elem[0][eg],2)))*e->dphi[kk][jj][ii][0][eg]
                 -e->dphi[kk][jj][ii][1][eg]*n_elem[0][eg]*n_elem[1][eg]
                 -e->dphi[kk][jj][ii][2][eg]*n_elem[0][eg]*n_elem[2][eg])
                *4*(pow((u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg]),3))*(pow(v_mag_elem[eg],1))*e->weight[eg];
                
                Kd_ele[l] += (-e->dphi[k][j][i][0][eg]*n_elem[0][eg]*n_elem[1][eg]
                              +(1.-pow(n_elem[1][eg],2))*e->dphi[k][j][i][1][eg]
                              -e->dphi[k][j][i][2][eg]*n_elem[1][eg]*n_elem[2][eg])
                *(-e->dphi[kk][jj][ii][0][eg]*n_elem[0][eg]*n_elem[1][eg]
                  +(1.-pow(n_elem[1][eg],2))*e->dphi[kk][jj][ii][1][eg]
                  -e->dphi[kk][jj][ii][2][eg]*n_elem[1][eg]*n_elem[2][eg])
                *4*(pow((u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg]),3))*(pow(v_mag_elem[eg],1))*e->weight[eg];
                
                Kd_ele[l] += (-e->dphi[k][j][i][0][eg]*n_elem[0][eg]*n_elem[2][eg]
                              -e->dphi[k][j][i][1][eg]*n_elem[1][eg]*n_elem[2][eg]
                              +(1.-pow(n_elem[2][eg],2))*e->dphi[k][j][i][2][eg])
                *(-e->dphi[kk][jj][ii][0][eg]*n_elem[0][eg]*n_elem[2][eg]
                  -e->dphi[kk][jj][ii][1][eg]*n_elem[1][eg]*n_elem[2][eg]
                  +(1.-pow(n_elem[2][eg],2))*e->dphi[kk][jj][ii][2][eg])
                *4*(pow((u_elem[0][eg]*n_elem[0][eg] + u_elem[1][eg]*n_elem[1][eg] + u_elem[2][eg]*n_elem[2][eg]),3))*(pow(v_mag_elem[eg],1))*e->weight[eg];
              }
						}
					}
				}
			}
		}
	}
  ierr = PetscFree6(dv_elem[0],dv_elem[1],dv_elem[2],u_elem[0],u_elem[1],u_elem[2]);CHKERRQ(ierr);
  ierr = PetscFree5(n_elem[0],n_elem[1],n_elem[2],v_mag_elem,v_elem);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}





#undef __FUNCT__
#define __FUNCT__ "VF_MatFluidCompreStiffMatrix_local"
extern PetscErrorCode VF_MatFluidCompreStiffMatrix_local(PetscReal *Kd_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ****u_array,PetscReal ***v_array)
{
  PetscErrorCode          ierr;
  PetscInt                i,j,k,l,c;
  PetscInt                ii,jj,kk;
  PetscReal               *dv_elem[3],*u_elem[3];
  PetscInt                eg;
  
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
							Kd_ele[l] = 0.;
							for (eg = 0; eg < e->ng; eg++) {
                for(c = 0; c < 3; c++){
                  Kd_ele[l] += (u_elem[c][eg]*dv_elem[c][eg])*e->phi[k][j][i][eg]*e->weight[eg];
                }
              }
						}
					}
				}
			}
		}
	}
  ierr = PetscFree6(dv_elem[0],dv_elem[1],dv_elem[2],u_elem[0],u_elem[1],u_elem[2]);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}









#undef __FUNCT__
#define __FUNCT__ "VF_FlowStiffnessMatrixEps_local"
extern PetscErrorCode VF_FlowStiffnessMatrixEps_local(PetscReal *Kd_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ***v_array)
{
  PetscErrorCode  ierr;
	PetscInt        i,j,k,l;
	PetscInt        ii,jj,kk;
	PetscInt        eg;
	PetscReal       kx_ep,ky_ep,kz_ep;
  PetscReal		   *v_elem;
  PetscReal		   val = 1e-15;
  
	PetscFunctionBegin;
  kx_ep = val;
  ky_ep = val;
  kz_ep = val;
  ierr = PetscMalloc(e->ng*sizeof(PetscReal),&v_elem);CHKERRQ(ierr);
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
                Kd_ele[l] += ((kx_ep*e->dphi[kk][jj][ii][0][eg]*e->dphi[k][j][i][0][eg])
                              +(ky_ep*e->dphi[kk][jj][ii][1][eg]*e->dphi[k][j][i][1][eg])
                              +(kz_ep*e->dphi[kk][jj][ii][2][eg]*e->dphi[k][j][i][2][eg]))*(1-v_elem[eg])*(1-v_elem[eg])*e->weight[eg];
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