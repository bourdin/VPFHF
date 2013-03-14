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
#define __FUNCT__ "MixedFlowFEMKSPSolve"
extern PetscErrorCode MixedFlowFEMKSPSolve(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode     ierr;
	PetscViewer        viewer;
	PetscInt           xs,xm,ys,ym,zs,zm;
	PetscInt           i,j,k,c,veldof = 3;
	PetscInt           its;
	KSPConvergedReason reason;
	PetscReal          ****VelnPress_array;
	PetscReal          ***Press_array;
	Vec                VecRHS;
	PetscReal          ****vel_array;
	PetscReal		   ****velbc_array;
	Vec					velbc_local;
	PetscReal			theta,one_minus_theta,timestepsize;
	PetscReal			dt_dot_theta,dt_dot_one_minus_theta;
	Vec                vec;
	PetscReal           VelPmin,VelPmax;

	PetscFunctionBegin;
	timestepsize = ctx->flowprop.timestepsize;
	theta = ctx->flowprop.theta;
	one_minus_theta = (1.-theta);
	dt_dot_theta = timestepsize*theta;
	dt_dot_one_minus_theta = timestepsize*(1.-theta);
	ierr = VecDuplicate(ctx->RHSVelP,&VecRHS);CHKERRQ(ierr);
	ierr = VecDuplicate(ctx->RHSVelP,&vec);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daFlow,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = FlowMatnVecAssemble(ctx->KVelP,ctx->KVelPlhs,ctx->RHSVelP,fields,ctx);CHKERRQ(ierr);
	ierr = VecCopy(ctx->RHSVelP,VecRHS);CHKERRQ(ierr);
	ierr = VecAXPBY(VecRHS,dt_dot_one_minus_theta,dt_dot_theta,ctx->RHSVelPpre);CHKERRQ(ierr);
	ierr = VecCopy(ctx->RHSVelP,ctx->RHSVelPpre);CHKERRQ(ierr);
	ierr = MatMultAdd(ctx->KVelPlhs,fields->VelnPress,VecRHS,VecRHS);CHKERRQ(ierr);
	
	ierr = DMGetLocalVector(ctx->daVect,&velbc_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVect,ctx->VelBCArray,INSERT_VALUES,velbc_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVect,ctx->VelBCArray,INSERT_VALUES,velbc_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,velbc_local,&velbc_array);CHKERRQ(ierr); 
	ierr = VecApplyVelocityBC(VecRHS,&ctx->bcQ[0],ctx,velbc_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,velbc_local,&velbc_array);CHKERRQ(ierr); 
	ierr = DMRestoreLocalVector(ctx->daVect,&velbc_local);CHKERRQ(ierr);
	ierr = KSPMonitorSet(ctx->kspVelP,MixedFEMKSPMonitor,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = KSPSolve(ctx->kspVelP,VecRHS,fields->VelnPress);CHKERRQ(ierr);
	
	/*
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
		ierr = PetscPrintf(PETSC_COMM_WORLD,"[ERROR] kspVelP diverged with reason %d\n",(int)reason);CHKERRQ(ierr);
	} else {
		ierr = KSPGetIterationNumber(ctx->kspVelP,&its);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"      kspVelP converged in %d iterations %d.\n",(int)its,(int)reason);CHKERRQ(ierr);
	}
	/*The next few lines equate the values of pressure calculated from the flow solver, to the pressure defined in da=daScal*/
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
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MixedFEMKSPMonitor"
extern PetscErrorCode MixedFEMKSPMonitor(KSP ksp,PetscInt its,PetscReal fnorm,void* ptr)
{
	PetscErrorCode ierr;
	PetscReal      norm,vmax,vmin;
	MPI_Comm       comm;
	Vec				solution;
	
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
#define __FUNCT__ "ReSETFlowBC"
extern PetscErrorCode ReSETFlowBC(BC *bcP,BC *bcQ, FlowCases flowcase)
{
	PetscInt i,c;
	
	PetscFunctionBegin;
	for (i = 0; i < 6; i++) {
		bcP[0].face[i] = NONE;
		for (c = 0; c < 3; c++) {
			bcQ[c].face[i] = NONE;
		}
	}
	for (i = 0; i < 12; i++) {
		bcP[0].face[i] = NONE;
		for (c = 0; c < 3; c++) {
			bcQ[c].edge[i] = NONE;
		}
	}
	for (i = 0; i < 8; i++) {
		bcP[0].face[i] = NONE;
		for (c = 0; c < 3; c++) {
			bcQ[c].vertex[i] = NONE;
		}
	}
	switch (flowcase) {
		case ALLPRESSUREBC:
			for (i = 0; i < 6; i++) {
				bcP[0].face[i] = VALUE;
			}
			break;
		case ALLNORMALFLOWBC:			
			bcQ[0].face[X0] = VALUE;
			bcQ[0].face[X1] = VALUE;
			bcQ[1].face[Y0] = VALUE;
			bcQ[1].face[Y1] = VALUE;
			bcQ[2].face[Z0] = VALUE;
			bcQ[2].face[Z1] = VALUE;
			break;
		default:
			SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: [%s] unknown FLOWCASE %i .\n",__FUNCT__,flowcase);
			break;
	}
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "VecApplySourceTerms"
extern PetscErrorCode VecApplySourceTerms(PetscReal *Ks_local,PetscReal ***source_array,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,VFCtx *ctx)
{
	PetscErrorCode ierr;
	PetscInt       i,j,k,l;
	PetscReal      alpha_c;
	PetscReal      mu;
	PetscReal      *loc_source;
	PetscInt       eg;
	
	PetscFunctionBegin;
	alpha_c = ctx->flowprop.alpha;
	mu     = ctx->flowprop.mu;
	ierr   = PetscMalloc(e->ng*sizeof(PetscReal),&loc_source);CHKERRQ(ierr);
	
	for (eg = 0; eg < e->ng; eg++) loc_source[eg] = 0.;
	for (k = 0; k < e->nphiz; k++) {
		for (j = 0; j < e->nphiy; j++) {
			for (i = 0; i < e->nphix; i++) {
				for (eg = 0; eg < e->ng; eg++) {
					loc_source[eg] += source_array[ek+k][ej+j][ei+i]*e->phi[k][j][i][eg];
				}
			}
		}
	}
	for (eg = 0; eg < e->ng; eg++)
		for (l = 0,k = 0; k < e->nphiz; k++) {
			for (j = 0; j < e->nphiy; j++) {
				for (i = 0; i < e->nphix; i++,l++) {
					Ks_local[l] = 0.;
					for (eg = 0; eg < e->ng; eg++) {
						Ks_local[l] += -loc_source[eg]*e->phi[k][j][i][eg]*e->weight[eg]/alpha_c;
					}
				}
			}
		}
	PetscFree(loc_source);
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "VecApplyVelocityBC"
extern PetscErrorCode VecApplyVelocityBC(Vec RHS,BC *bcQ,VFCtx *ctx, PetscReal ****velbc_array)
{
	PetscErrorCode ierr;
	PetscInt       xs,xm,nx;
	PetscInt       ys,ym,ny;
	PetscInt       zs,zm,nz;
	PetscInt       dim,dof=3;
	PetscInt       i,j,k,c;
	DM             da;
	PetscReal      ****RHS_array;
	PetscReal      hx,hy,hz;
	
	PetscFunctionBegin;
	ierr = PetscObjectQuery((PetscObject)RHS,"DM",(PetscObject*)&da);CHKERRQ(ierr);
	if (!da) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Vector not generated from a DMDA");
	
	ierr = DMDAGetInfo(da,&dim,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(da,RHS,&RHS_array);CHKERRQ(ierr);
	hx   = 1./(nx-1);hy = 1./(ny-1);hz = 1./(nz-1);
	for (c = 0; c < dof; c++) {
		if (xs == 0) {
			i = 0;
			if (bcQ[c].face[X0] == VALUE) {
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
					}
				}
			}
		}
		if (xs+xm == nx) {
			i = nx-1;
			if (bcQ[c].face[X1] == VALUE) {
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
					}
				}
			}
		}
		if (ys == 0) {
			j = 0;
			if (bcQ[c].face[Y0] == VALUE) {
				for (k = zs; k < zs+zm; k++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
					}
				}
			}
		}
		if (ys+ym == ny) {
			j = ny-1;
			if (bcQ[c].face[Y1] == VALUE) {
				for (k = zs; k < zs+zm; k++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
						
					}
				}
			}
		}
		if (zs == 0) {
			k = 0;
			if (bcQ[c].face[Z0] == VALUE) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
					}
				}
			}
		}
		if (zs+zm == nz) {
			k = nz-1;
			if (bcQ[c].face[Z1] == VALUE) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
					}
				}
			}
		}
		if (xs == 0 && zs == 0) {
			k = 0;i = 0;
			if (bcQ[c].edge[X0Z0] == VALUE) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && zs == 0) {
			k = 0;i = nx-1;
			if (bcQ[c].edge[X1Z0] == VALUE) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (ys == 0 && zs == 0) {
			k = 0;j = 0;
			if (bcQ[c].edge[Y0Z0] == VALUE) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (ys+ym == ny && zs == 0) {
			k = 0;j = 0;
			if (bcQ[c].edge[Y1Z0] == VALUE) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && zs+zm == nz) {
			k = nz-1;i = 0;
			if (bcQ[c].edge[X0Z1] == VALUE) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && zs+zm == nz) {
			k = nz-1;i = nx-1;
			if (bcQ[c].edge[X1Z1] == VALUE) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (ys == 0 && zs+zm == nz) {
			k = nz-1;j = 0;
			if (bcQ[c].edge[Y0Z1] == VALUE) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (ys+ym == ny && zs+zm == nz) {
			k = nz-1;j = ny-1;
			if (bcQ[c].edge[Y1Z1] == VALUE) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && ys == 0) {
			j = 0;i = 0;
			if (bcQ[c].edge[X0Y0] == VALUE) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && ys+ym == ny) {
			j = ny-1;i = 0;
			if (bcQ[c].edge[X0Y1] == VALUE) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && ys == 0) {
			j = 0;i = nx-1;
			if (bcQ[c].edge[X1Y0] == VALUE) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && ys+ym == ny) {
			j = ny-1;i = nx-1;
			if (bcQ[c].edge[X1Y1] == VALUE) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && ys == 0 && zs == 0) {
			k = 0;j = 0;i = 0;
			if (bcQ[c].vertex[X0Y0Z0] == VALUE) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys == 0 && zs == 0) {
			k = 0;j = 0;i = nx-1;
			if (bcQ[c].vertex[X1Y0Z0] == VALUE) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs == 0 && ys+ym == ny && zs == 0) {
			k = 0;j = ny-1;i = 0;
			if (bcQ[c].vertex[X0Y1Z0] == VALUE) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys+ym == ny && zs == 0) {
			k = 0;j = ny-1;i = nx-1;
			if (bcQ[c].vertex[X1Y1Z0] == VALUE) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs == 0 && ys == 0 && zs+zm == nz) {
			k = nz-1;j = 0;i = 0;
			if (bcQ[c].vertex[X0Y0Z1] == VALUE) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys == 0 && zs+zm == nz) {
			k = nz-1;j = 0;i = nx-1;
			if (bcQ[c].vertex[X1Y0Z1] == VALUE) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs == 0 && ys+ym == ny && zs+zm == nz) {
			k = nz-1;j = ny-1;i = 0;
			if (bcQ[c].vertex[X0Y1Z1] == VALUE) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys+ym == ny && zs+zm == nz) {
			k = nz-1;j = ny-1;i = nx-1;
			if (bcQ[c].vertex[X1Y1Z1] == VALUE) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
	}
	ierr = DMDAVecRestoreArrayDOF(da,RHS,&RHS_array);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ReSETSourceTerms"
extern PetscErrorCode ReSETSourceTerms(Vec Src,FlowProp flowpropty)
{
	PetscErrorCode ierr;
	PetscReal      pi;
	PetscReal      ***source_array;
	PetscInt       xs,xm,nx;
	PetscInt       ys,ym,ny;
	PetscInt       zs,zm,nz;
	PetscInt       dim,dof;
	PetscInt       ei,ej,ek,c;
	DM             da;
	PetscReal      mu,beta_c;
	PetscReal      hx,hy,hz;
	
	PetscFunctionBegin;
	beta_c = flowpropty.beta;
	mu     = flowpropty.mu;
	ierr   = PetscObjectQuery((PetscObject)Src,"DM",(PetscObject*)&da);CHKERRQ(ierr);
	if (!da) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Vector not generated from a DA");
	
	ierr = DMDAGetInfo(da,&dim,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   &dof,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da,Src,&source_array);CHKERRQ(ierr);
	pi   = 6.*asin(0.5);
	hx   = 1./(nx-1);
	hy   = 1./(nx-1);
	hz   = 1./(nz-1);
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				source_array[ek][ej][ei] = 0.; 
			}
		}
	}
	ierr = DMDAVecRestoreArray(da,Src,&source_array);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FlowMatnVecAssemble"
extern PetscErrorCode FlowMatnVecAssemble(Mat K,Mat Krhs,Vec RHS,VFFields * fields,VFCtx *ctx)
{
	PetscErrorCode ierr;
	PetscInt       xs,xm,nx;
	PetscInt       ys,ym,ny;
	PetscInt       zs,zm,nz;
	PetscInt       ek,ej,ei;
	PetscInt       i,j,k,l,ii;
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
	PetscReal      beta_c,alpha_c,mu,gx,gy,gz;
	PetscReal	   theta,timestepsize;
	PetscInt       nrow = ctx->e3D.nphix*ctx->e3D.nphiy*ctx->e3D.nphiz;
	MatStencil     *row,*row1;
	PetscReal      ***source_array;
	Vec            source_local;
	PetscReal      M_inv;
	Mat			   K1,K2;
	PetscReal		time_theta,time_one_minus_theta;
	FACE			face;
	PetscReal		***prebc_array;
	Vec				prebc_local;
	MatStructure	flg;
	PetscReal		hwx, hwy, hwz;	
	
	PetscFunctionBegin;
	flg = SAME_NONZERO_PATTERN;
	M_inv     = ctx->flowprop.M_inv;
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
	ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	/* This line ensures that the number of cells is one less than the number of nodes. Force processing of cells to stop once the second to the last node is processed */
	ierr = MatZeroEntries(K);CHKERRQ(ierr);
	ierr = MatZeroEntries(Krhs);CHKERRQ(ierr);
	ierr = MatDuplicate(K,MAT_COPY_VALUES,&K1);CHKERRQ(ierr);
	ierr = MatDuplicate(K,MAT_COPY_VALUES,&K2);CHKERRQ(ierr);
	ierr = VecSet(RHS,0.);CHKERRQ(ierr);
	/* Get coordinates from daVect since ctx->coordinates was created as an object in daVect */
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
	ierr = DMGlobalToLocalEnd(ctx->daScal,fields->PresBCArray,INSERT_VALUES,prebc_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,prebc_local,&prebc_array);CHKERRQ(ierr); 
	ierr = PetscMalloc5(nrow*nrow,PetscReal,&KA_local,
						nrow*nrow,PetscReal,&KB_local,
						nrow*nrow,PetscReal,&KD_local,
						nrow*nrow,PetscReal,&KBTrans_local,
						nrow*nrow,PetscReal,&KS_local);CHKERRQ(ierr);	
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
				ierr = FLow_MatA(KA_local,&ctx->e3D,ek,ej,ei);CHKERRQ(ierr);
				for (l = 0; l < nrow*nrow; l++) {
					KS_local[l] = -2.*M_inv*KA_local[l]/alpha_c;
				}
				for (c = 0; c < veldof; c++) {
					ierr = FLow_MatB(KB_local,&ctx->e3D,ek,ej,ei,c);CHKERRQ(ierr);
					ierr = FLow_MatBTranspose(KBTrans_local,&ctx->e3D,ek,ej,ei,c,ctx->flowprop,perm_array);CHKERRQ(ierr);
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
				ierr = FLow_MatD(KD_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,perm_array);CHKERRQ(ierr);
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							row[l].i = ei+i;row[l].j = ej+j;row[l].k = ek+k;row[l].c = 3;
						}
					}
				}
				ierr = MatSetValuesStencil(K1,nrow,row,nrow,row,KD_local,ADD_VALUES);CHKERRQ(ierr);
				ierr = MatSetValuesStencil(K2,nrow,row,nrow,row,KS_local,ADD_VALUES);CHKERRQ(ierr);
				/*Assembling the righthand side vector f*/
				for (c = 0; c < veldof; c++) {
					ierr = FLow_Vecf(RHS_local,&ctx->e3D,ek,ej,ei,c,ctx->flowprop,perm_array);CHKERRQ(ierr);
					for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
						for (j = 0; j < ctx->e3D.nphiy; j++) {
							for (i = 0; i < ctx->e3D.nphix; i++,l++) {
								RHS_array[ek+k][ej+j][ei+i][c] += RHS_local[l];
							}
						}
					}
				}
				/*Assembling the righthand side vector g*/
				ierr = FLow_Vecg(RHS_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,perm_array);CHKERRQ(ierr);
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							RHS_array[ek+k][ej+j][ei+i][3] += RHS_local[l];
						}
					}
				}
				if(ctx->hasFluidSources){
					ierr = VecApplySourceTerms(RHS_local,source_array,&ctx->e3D,ek,ej,ei,ctx);
					for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
						for (j = 0; j < ctx->e3D.nphiy; j++) {
							for (i = 0; i < ctx->e3D.nphix; i++,l++) {
								RHS_array[ek+k][ej+j][ei+i][3] += RHS_local[l];
							}
						}
					}
				}
				if (ei == 0) {
					/*					 Face X0			*/
					face = X0;	
					ierr = CartFE_Element2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == VALUE) {
						ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array);CHKERRQ(ierr);
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
					/*					 Face X1		*/
					face = X1;
					ierr = CartFE_Element2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == VALUE) {
						ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array);CHKERRQ(ierr);
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
					/*					 Face Y0		*/
					face = Y0;
					ierr = CartFE_Element2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == VALUE) {
						ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array);CHKERRQ(ierr);
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
					/*					 Face Y1		*/
					face = Y1;
					ierr = CartFE_Element2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == VALUE) {
						ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array);CHKERRQ(ierr);
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
					/*					 Face Z0		*/
					face = Z0;
					ierr = CartFE_Element2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == VALUE) {
						ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array);CHKERRQ(ierr);
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
					/*					 Face Z1		*/
					face = Z1;
					ierr = CartFE_Element2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == VALUE) {
						ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array);CHKERRQ(ierr);
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
	if(ctx->hasFlowWells){
		for(ii = 0; ii < ctx->numWells; ii++){
			for (ek = zs; ek < zs+zm; ek++) {
				for (ej = ys; ej < ys+ym; ej++) {
					for (ei = xs; ei < xs+xm; ei++) {
						if(  
						   ((coords_array[ek][ej][ei+1][0] >= ctx->well[ii].coords[0]) && (coords_array[ek][ej][ei][0] <= ctx->well[ii].coords[0] ))	&&
						   ((coords_array[ek][ej+1][ei][1] >= ctx->well[ii].coords[1]) && (coords_array[ek][ej][ei][1] <= ctx->well[ii].coords[1] ))	&&
						   ((coords_array[ek+1][ej][ei][2] >= ctx->well[ii].coords[2]) && (coords_array[ek][ej][ei][2] <= ctx->well[ii].coords[2] ))
						   )
						{
							hwx = (ctx->well[ii].coords[0]-coords_array[ek][ej][ei][0])/(coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0]); 
							hwy = (ctx->well[ii].coords[1]-coords_array[ek][ej][ei][1])/(coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1]); 
							hwz = (ctx->well[ii].coords[2]-coords_array[ek][ej][ei][2])/(coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2]); 
							if(ctx->well[ii].condition == RATE){
								ierr = VecApplyWEllFlowRate(RHS_local,ctx->well[ii].Qw,hwx,hwy,hwz);CHKERRQ(ierr);
								for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
									for (j = 0; j < ctx->e3D.nphiy; j++) {
										for (i = 0; i < ctx->e3D.nphix; i++,l++) {
											if(ctx->well[ii].type == INJECTOR){
												RHS_array[ek+k][ej+j][ei+i][3] -= RHS_local[l];
											}
											else if(ctx->well[ii].type == PRODUCER)
											{
												RHS_array[ek+k][ej+j][ei+i][3] -= -1*RHS_local[l];
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
	ierr = MatAXPY(Krhs,time_one_minus_theta,K1,flg);CHKERRQ(ierr);
	ierr = MatAXPY(Krhs,1.0,K2,flg);CHKERRQ(ierr);
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
	ierr = PetscFree5(KA_local,KB_local,KD_local,KBTrans_local,KS_local);CHKERRQ(ierr);
	ierr = PetscFree3(RHS_local,row,row1);CHKERRQ(ierr);	
	ierr = MatDestroy(&K1);CHKERRQ(ierr);
	ierr = MatDestroy(&K2);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecApplyWEllFlowRate"
extern PetscErrorCode VecApplyWEllFlowRate(PetscReal *RHS_local,PetscReal Q,PetscReal hwx,PetscReal hwy,PetscReal hwz)
{
	PetscReal	phi[2][2][2];
	PetscInt	i,j,k,l;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	phi[0][0][0] = (1.-hwz)*(1.-hwy)*(1.-hwx);
	phi[0][0][1] = (1.-hwz)*(1.-hwy)*hwx;
	phi[0][1][0] = (1.-hwz)*hwy*(1.-hwx);
	phi[0][1][1] = (1.-hwz)*hwy*hwx;
	phi[1][0][0] = hwz*(1.-hwy)*(1.-hwx);
	phi[1][0][1] = hwz*(1.-hwy)*hwx;
	phi[1][1][0] = hwz*hwy*(1.-hwx);
	phi[1][1][1] = hwz*hwy*hwx;
	
	for(l = 0, k = 0; k < 2; k++){
		for(j = 0; j < 2; j++){
			for(i = 0; i < 2; i++, l++){
				RHS_local[l] = 0;
				RHS_local[l] = Q*phi[k][j][i];
			}
		}
	}
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecApplyPressureBC"
extern PetscErrorCode VecApplyPressureBC(PetscReal *RHS_local,PetscReal ***pre_array,PetscInt ek,PetscInt ej,PetscInt ei,FACE face,CartFE_Element2D *e,FlowProp flowpropty,PetscReal ****perm_array)
{
	PetscErrorCode ierr;
	PetscInt       i,j,k,l,g;
	PetscReal      *pre_elem;
	PetscReal		beta_c,mu;
	PetscReal		kx,ky,kz;
	
	PetscFunctionBegin;
	beta_c  = flowpropty.beta;
	mu      = flowpropty.mu;
	ierr = PetscMalloc(e->ng*sizeof(PetscReal),&pre_elem);CHKERRQ(ierr);
		// Initialize pre_Elem	 
	for (g = 0; g < e->ng; g++) {
		pre_elem[g] = 0;
	}
	switch (face) {
		case X0:
			kx = perm_array[ek][ej][ei][0];
			for (k = 0; k < e->nphix; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphiz; i++) {
						for (g = 0; g < e->ng; g++) {
							pre_elem[g] += e->phi[i][j][k][g]*pre_array[ek+k][ej+j][ei+i];
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
							RHS_local[l] -= kx*beta_c/mu * e->weight[g]*e->phi[i][j][k][g]*pre_elem[g];
						}
					}
				}
			}
			break;
		case X1:
			kx = perm_array[ek][ej][ei][0];
			for (k = 0; k < e->nphix; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphiz; i++) {
						for (g = 0; g < e->ng; g++) {
							pre_elem[g] += e->phi[i][j][k][g]*pre_array[ek+k][ej+j][ei+1];
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
							RHS_local[l] -= kx*beta_c/mu * e->weight[g]*e->phi[i][j][k][g]*pre_elem[g];
						}
					}
				}
			}
			break;
		case Y0:
			ky = perm_array[ek][ej][ei][1];
			for (k = 0; k < e->nphiy; k++) {
				for (j = 0; j < e->nphiz; j++) {
					for (i = 0; i < e->nphix; i++) {
						for (g = 0; g < e->ng; g++) {
							pre_elem[g] += e->phi[j][k][i][g]*pre_array[ek+k][ej+j][ei+i];
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
							RHS_local[l] -= ky*beta_c/mu * e->weight[g]*e->phi[j][k][i][g]*pre_elem[g];
						}
					}
				}
			}
			break;
		case Y1:
			ky = perm_array[ek][ej][ei][1];
			for (k = 0; k < e->nphiy; k++) {
				for (j = 0; j < e->nphiz; j++) {
					for (i = 0; i < e->nphix; i++) {
						for (g = 0; g < e->ng; g++) {
							pre_elem[g] += e->phi[j][k][i][g]*pre_array[ek+k][ej+1][ei+i];
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
							RHS_local[l] -= ky*beta_c/mu * e->weight[g]*e->phi[j][k][i][g]*pre_elem[g];
						}
					}
				}
			}
			break;
		case Z0:
			kz = perm_array[ek][ej][ei][2];
			for (k = 0; k < e->nphiz; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphix; i++) {
						for (g = 0; g < e->ng; g++) {
							pre_elem[g] += e->phi[k][j][i][g]*pre_array[ek][ej+j][ei+i];
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
							RHS_local[l] -= kz*beta_c/mu * e->weight[g]*e->phi[k][j][i][g]*pre_elem[g];
						}
					}
				}
			}
			break;
		case Z1:
			kz = perm_array[ek][ej][ei][2];
			for (k = 0; k < e->nphiz; k++) {
				for (j = 0; j < e->nphiy; j++) {
					for (i = 0; i < e->nphix; i++) {
						for (g = 0; g < e->ng; g++) {
							pre_elem[g] += e->phi[k][j][i][g]*pre_array[ek+1][ej+j][ei+i];
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
							RHS_local[l] -= kz*beta_c/mu * e->weight[g]*e->phi[k][j][i][g]*pre_elem[g];
						}
					}
				}
			}
			break;
	}
	ierr = PetscFree(pre_elem);CHKERRQ(ierr);
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
	DM				da;
	PetscReal      zero=0.0;
	
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
		if (xs == 0       && bcQ[c].face[X0] == VALUE)             numBC += ym * zm;
		if (xs + xm == nx && bcQ[c].face[X1] == VALUE)             numBC += ym * zm;
		if (ys == 0       && bcQ[c].face[Y0] == VALUE)             numBC += xm * zm;
		if (ys + ym == ny && bcQ[c].face[Y1] == VALUE)             numBC += xm * zm;
		if (zs == 0       && bcQ[c].face[Z0] == VALUE && dim == 3) numBC += xm * ym;
		if (zs + zm == nz && bcQ[c].face[Z1] == VALUE && dim == 3) numBC += xm * ym;
		if (xs == 0       && ys == 0       && zs == 0       && bcQ[c].vertex[X0Y0Z0] == VALUE) numBC++;
		if (xs == 0       && ys + ym == ny && zs == 0       && bcQ[c].vertex[X0Y1Z0] == VALUE) numBC++;
		if (xs + xm == nx && ys == 0       && zs == 0       && bcQ[c].vertex[X1Y0Z0] == VALUE) numBC++;
		if (xs + xm == nx && ys + ym == ny && zs == 0       && bcQ[c].vertex[X1Y1Z0] == VALUE) numBC++;
		if (xs == 0       && ys == 0       && zs + zm == nz && bcQ[c].vertex[X0Y0Z1] == VALUE && dim == 3) numBC++;
		if (xs == 0       && ys + ym == ny && zs + zm == nz && bcQ[c].vertex[X0Y1Z1] == VALUE && dim == 3) numBC++;
		if (xs + xm == nx && ys == 0       && zs + zm == nz && bcQ[c].vertex[X1Y0Z1] == VALUE && dim == 3) numBC++;
		if (xs + xm == nx && ys + ym == ny && zs + zm == nz && bcQ[c].vertex[X1Y1Z1] == VALUE && dim == 3) numBC++;
	}
	ierr = PetscMalloc(numBC * sizeof(MatStencil),&row);CHKERRQ(ierr);
	/*
	 Create an array of rows to be zeroed out
	 */
	/*
	 i == 0
	 */
	for (c = 0; c < dof; c++) {
		if (xs == 0 && bcQ[c].face[X0] == VALUE) {
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
		if (xs + xm == nx && bcQ[c].face[X1] == VALUE) {
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
		if (ys == 0 && bcQ[c].face[Y0] == VALUE) {
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
		if (ys + ym == ny && bcQ[c].face[Y1] == VALUE) {
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
			if (zs == 0 && bcQ[c].face[Z0] == VALUE) {
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
			if (zs + zm == nz && bcQ[c].face[Z1] == VALUE) {
				for (j = ys; j < ys + ym; j++) {
					for (i = xs; i < xs + xm; i++) {
						row[l].i = i; row[l].j = j; row[l].k = nz-1; row[l].c = c; 
						l++;
					}
				}
			}
		}
		if (xs == 0       && ys == 0       && zs == 0       && bcQ[c].vertex[X0Y0Z0] == VALUE) { 
			row[l].i = 0; row[l].j = 0; row[l].k = 0; row[l].c = c; 
			l++;
		}
		if (xs == 0       && ys == 0       && zs + zm == nz && bcQ[c].vertex[X0Y0Z1] == VALUE && dim ==3) { 
			row[l].i = 0; row[l].j = 0; row[l].k = nz-1; row[l].c = c; 
			l++;
		}
		if (xs == 0       && ys + ym == ny && zs == 0       && bcQ[c].vertex[X0Y1Z0] == VALUE) { 
			row[l].i = 0; row[l].j = ny-1; row[l].k = 0; row[l].c = c; 
			l++;
		}
		if (xs == 0       && ys + ym == ny && zs + zm == nz && bcQ[c].vertex[X0Y1Z1] == VALUE && dim ==3) { 
			row[l].i = 0; row[l].j = ny-1; row[l].k = nz-1; row[l].c = c; 
			l++;
		}
		if (xs + xm == nx && ys == 0       && zs == 0       && bcQ[c].vertex[X1Y0Z0] == VALUE) { 
			row[l].i = nx-1; row[l].j = 0; row[l].k = 0; row[l].c = c; 
			l++;
		}
		if (xs + xm == nx && ys == 0       && zs + zm == nz && bcQ[c].vertex[X1Y0Z1] == VALUE && dim ==3) { 
			row[l].i = nx-1; row[l].j = 0; row[l].k = nz-1; row[l].c = c; 
			l++;
		}
		if (xs + xm == nx && ys + ym == ny && zs == 0       && bcQ[c].vertex[X1Y1Z0] == VALUE) { 
			row[l].i = nx-1; row[l].j = ny-1; row[l].k = 0; row[l].c = c; 
			l++;
		}
		if (xs + xm == nx && ys + ym == ny && zs + zm == nz && bcQ[c].vertex[X1Y1Z1] == VALUE && dim ==3) { 
			row[l].i = nx=1; row[l].j = ny-1; row[l].k = nz-1; row[l].c = c; 
			l++;
		}
		
	}
	ierr = MatZeroRowsStencil(K,numBC,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
		//	ierr = MatZeroRowsStencil(Klhs,numBC,row,zero,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscFree(row);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FLow_Vecg"
extern PetscErrorCode FLow_Vecg(PetscReal *Kg_local,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,FlowProp flowpropty,PetscReal ****perm_array)
{
	PetscInt  i,j,k,l;
	PetscInt  eg;
	PetscReal beta_c,mu,rho,gamma_c,gx,gy,gz;
	PetscReal kx,ky,kz,kxy,kxz,kyz;
	
	PetscFunctionBegin;
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
	
	for (l = 0,k = 0; k < e->nphiz; k++) {
		for (j = 0; j < e->nphiy; j++) {
			for (i = 0; i < e->nphix; i++,l++) {
				Kg_local[l] = 0.;
				for (eg = 0; eg < e->ng; eg++) {
					/* Need to multiply by the permability when it is available*/
					Kg_local[l] += -0.5*rho*gamma_c*beta_c/mu*((kx*gx+kxy*gy+kxz*gz)*e->dphi[k][j][i][0][eg]
															   +(kxy*gx+ky*gy+kyz*gz)*e->dphi[k][j][i][1][eg]
															   +(kxz*gx+kyz*gy+kz*gz)*e->dphi[k][j][i][2][eg])*e->weight[eg];
				}
			}
		}
	}
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "FLow_Vecf"
extern PetscErrorCode FLow_Vecf(PetscReal *Kf_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscInt c,FlowProp flowpropty,PetscReal ****perm_array)
{
	PetscInt  i,j,k,l;
	PetscInt  eg;
	PetscReal beta_c,mu,rho,gamma_c;
	PetscReal k1,k2,k3;
	PetscReal gx,gy,gz;
	
	PetscFunctionBegin;
	beta_c  = flowpropty.beta;
	gamma_c = flowpropty.gamma;
	rho     = flowpropty.rho;
	mu      = flowpropty.mu;
	gx      = flowpropty.g[0];
	gy      = flowpropty.g[1];
	gz      = flowpropty.g[2];
	if (c == 0) {
		k1 = perm_array[ek][ej][ei][0];
		k2 = perm_array[ek][ej][ei][3];
		k3 = perm_array[ek][ej][ei][4];
	}
	if (c == 1) {
		k1 = perm_array[ek][ej][ei][3];
		k2 = perm_array[ek][ej][ei][1];
		k3 = perm_array[ek][ej][ei][5];
	}
	if (c == 2) {
		k1 = perm_array[ek][ej][ei][4];
		k2 = perm_array[ek][ej][ei][5];
		k3 = perm_array[ek][ej][ei][2];
	}
	for (l = 0,k = 0; k < e->nphiz; k++) {
		for (j = 0; j < e->nphiy; j++) {
			for (i = 0; i < e->nphix; i++,l++) {
				Kf_ele[l] = 0.;
				for (eg = 0; eg < e->ng; eg++) {
					/* Need to multiply by the permability when it is available*/
					Kf_ele[l] += 0.5*(k1*gx+k2*gy+k3*gz)*gamma_c*beta_c*rho/mu*e->phi[k][j][i][eg]*e->weight[eg];
				}
			}
		}
	}
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "FLow_MatD"
extern PetscErrorCode FLow_MatD(PetscReal *Kd_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,FlowProp flowpropty,PetscReal ****perm_array)
{
	PetscInt  i,j,k,l;
	PetscInt  ii,jj,kk;
	PetscInt  eg;
	PetscReal beta_c,mu;
	PetscReal kx,ky,kz,kxy,kxz,kyz;
	
	PetscFunctionBegin;
	beta_c = flowpropty.beta;
	mu     = flowpropty.mu;
	kx     = perm_array[ek][ej][ei][0];
	ky     = perm_array[ek][ej][ei][1];
	kz     = perm_array[ek][ej][ei][2];
	kxy    = perm_array[ek][ej][ei][3];
	kxz    = perm_array[ek][ej][ei][4];
	kyz    = perm_array[ek][ej][ei][5];
	
	for (l = 0,k = 0; k < e->nphiz; k++) {
		for (j = 0; j < e->nphiy; j++) {
			for (i = 0; i < e->nphix; i++) {
				for (kk = 0; kk < e->nphiz; kk++) {
					for (jj = 0; jj < e->nphiy; jj++) {
						for (ii = 0; ii < e->nphix; ii++,l++) {
							Kd_ele[l] = 0.;
							for (eg = 0; eg < e->ng; eg++) {
								/* Need to multiply by the permability when it is available*/
								Kd_ele[l] += -0.5*beta_c/mu*
								((kx*e->dphi[k][j][i][0][eg]+kxy*e->dphi[k][j][i][1][eg]+kxz*e->dphi[k][j][i][2][eg])*e->dphi[kk][jj][ii][0][eg]
								 +(kxy*e->dphi[k][j][i][0][eg]+ky*e->dphi[k][j][i][1][eg]+kyz*e->dphi[k][j][i][2][eg])*e->dphi[kk][jj][ii][1][eg]
								 +(kxz*e->dphi[k][j][i][0][eg]+kyz*e->dphi[k][j][i][1][eg]+kz*e->dphi[k][j][i][2][eg])*e->dphi[kk][jj][ii][2][eg])*e->weight[eg];
							}
						}
					}
				}
			}
		}
	}
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "FLow_MatB"
extern PetscErrorCode FLow_MatB(PetscReal *KB_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscInt c)
{
	PetscInt i,j,k,l;
	PetscInt ii,jj,kk;
	PetscInt eg;
	PetscFunctionBegin;
	for (l = 0,k = 0; k < e->nphiz; k++) {
		for (j = 0; j < e->nphiy; j++) {
			for (i = 0; i < e->nphix; i++) {
				for (kk = 0; kk < e->nphiz; kk++) {
					for (jj = 0; jj < e->nphiy; jj++) {
						for (ii = 0; ii < e->nphix; ii++,l++) {
							KB_ele[l] = 0.;
							for (eg = 0; eg < e->ng; eg++) {
								KB_ele[l] += -(e->phi[k][j][i][eg]*e->dphi[kk][jj][ii][c][eg]
											   +0.5*e->dphi[k][j][i][c][eg]*e->phi[kk][jj][ii][eg])*e->weight[eg];
							}
						}
					}
				}
			}
		}
	}
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "FLow_MatBTranspose"
extern PetscErrorCode FLow_MatBTranspose(PetscReal *KB_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscInt c,FlowProp flowpropty,PetscReal ****perm_array)
{
	PetscInt  i,j,k,l;
	PetscInt  ii,jj,kk;
	PetscInt  eg;
	PetscReal beta_c,mu;
	PetscReal k1,k2,k3;
	
	PetscFunctionBegin;
	beta_c = flowpropty.beta;
	mu     = flowpropty.mu;	
	if (c == 0) {
		k1 = perm_array[ek][ej][ei][0];
		k2 = perm_array[ek][ej][ei][3];
		k3 = perm_array[ek][ej][ei][4];
	}
	if (c == 1) {
		k1 = perm_array[ek][ej][ei][3];
		k2 = perm_array[ek][ej][ei][1];
		k3 = perm_array[ek][ej][ei][5];
	}
	if (c == 2) {
		k1 = perm_array[ek][ej][ei][4];
		k2 = perm_array[ek][ej][ei][5];
		k3 = perm_array[ek][ej][ei][2];
	}
	for (l = 0,k = 0; k < e->nphiz; k++) {
		for (j = 0; j < e->nphiy; j++) {
			for (i = 0; i < e->nphix; i++) {
				for (kk = 0; kk < e->nphiz; kk++) {
					for (jj = 0; jj < e->nphiy; jj++) {
						for (ii = 0; ii < e->nphix; ii++,l++) {
							KB_ele[l] = 0.;
							for (eg = 0; eg < e->ng; eg++) {
								KB_ele[l] += -beta_c/mu*(e->phi[kk][jj][ii][eg]*(k1*e->dphi[k][j][i][0][eg]+k2*e->dphi[k][j][i][1][eg]+k3*e->dphi[k][j][i][2][eg])
														 +0.5*(k1*e->dphi[kk][jj][ii][0][eg]+k2*e->dphi[kk][jj][ii][1][eg]+k3*e->dphi[kk][jj][ii][2][eg])*e->phi[k][j][i][eg])*e->weight[eg];
							}
						}
					}
				}
			}
		}
	}
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "FLow_MatA"
extern PetscErrorCode FLow_MatA(PetscReal *A_local,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei)
{
	PetscInt i,j,k,l;
	PetscInt ii,jj,kk;
	PetscInt eg;
	
	PetscFunctionBegin;
	for (l = 0,k = 0; k < e->nphiz; k++) {
		for (j = 0; j < e->nphiy; j++) {
			for (i = 0; i < e->nphix; i++) {
				for (kk = 0; kk < e->nphiz; kk++) {
					for (jj = 0; jj < e->nphiy; jj++) {
						for (ii = 0; ii < e->nphix; ii++,l++) {
							A_local[l] = 0;
							for (eg = 0; eg < e->ng; eg++) {
								A_local[l] += 0.5*e->phi[k][j][i][eg]*e->phi[kk][jj][ii][eg]*e->weight[eg];
							}
						}
					}
				}
			}
		}
	}
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "ReSETBoundaryTerms"
extern PetscErrorCode ReSETBoundaryTerms(VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode ierr;
	PetscReal		****perm_array;
	PetscReal		****vel_array;
	PetscReal		***press_array;
	PetscInt		xs,xm,nx;
	PetscInt		ys,ym,ny;
	PetscInt		zs,zm,nz;
	PetscInt		dim, dof;
	PetscInt		ei,ej,ek,c;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVFperm,fields->vfperm,&perm_array);CHKERRQ(ierr); 
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->VelBCArray,&vel_array);CHKERRQ(ierr); 
	ierr = DMDAVecGetArray(ctx->daScal,ctx->PresBCArray,&press_array);CHKERRQ(ierr); 
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				for(c = 3; c < 6; c++){
					perm_array[ek][ej][ei][c] = 0.;
				}
			}
		}
	}
	ierr = DMDAGetInfo(ctx->daFlow,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daFlow,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	for (ek = zs; ek < zs + zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				vel_array[ek][ej][ei][0] = 0.;
				vel_array[ek][ej][ei][1] = 0.;
				vel_array[ek][ej][ei][2] = 0.;
				press_array[ek][ej][ei] = 0.;
			}
		}
	}
	ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,fields->vfperm,&perm_array);CHKERRQ(ierr);	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->VelBCArray,&vel_array);CHKERRQ(ierr);	
	ierr = DMDAVecRestoreArray(ctx->daScal,ctx->PresBCArray,&press_array);CHKERRQ(ierr);	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MixedFEMFlowSolverInitialize"
extern PetscErrorCode MixedFEMFlowSolverInitialize(VFCtx *ctx, VFFields *fields)
{
	PetscMPIInt    comm_size;
	PetscErrorCode ierr;
	MatStructure flg;
	
	PetscFunctionBegin;
	ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"","");CHKERRQ(ierr);
	{
		ctx->units    = UnitaryUnits;
		ierr          = PetscOptionsEnum("-flowunits","\n\tFlow solver","",FlowUnitName,(PetscEnum)ctx->units,(PetscEnum*)&ctx->units,PETSC_NULL);CHKERRQ(ierr);
		/*	ctx->flowcase = ALLPRESSUREBC; */
			//		ctx->flowcase = ALLNORMALFLOWBC;
		ierr          = PetscOptionsEnum("-flow boundary conditions","\n\tFlow solver","",FlowBC_Case,(PetscEnum)ctx->flowcase,(PetscEnum*)&ctx->flowcase,PETSC_NULL);CHKERRQ(ierr);
	}
	ierr = PetscOptionsEnd();CHKERRQ(ierr);	
	/* MixedFEM ksp flow solver initialize*/
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&comm_size);CHKERRQ(ierr);
	if (comm_size == 1) {
		ierr = DMCreateMatrix(ctx->daFlow,MATSEQAIJ,&ctx->KVelP);CHKERRQ(ierr);
		ierr = DMCreateMatrix(ctx->daFlow,MATSEQAIJ,&ctx->KVelPlhs);CHKERRQ(ierr);
	} else {
		ierr = DMCreateMatrix(ctx->daFlow,MATMPIAIJ,&ctx->KVelP);CHKERRQ(ierr);
		ierr = DMCreateMatrix(ctx->daFlow,MATMPIAIJ,&ctx->KVelPlhs);CHKERRQ(ierr);
	}
	ierr = MatSetOption(ctx->KVelP,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
	ierr = MatSetOption(ctx->KVelPlhs,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daFlow,&ctx->RHSVelP);CHKERRQ(ierr);
	ierr = VecSet(ctx->RHSVelP,0.);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)ctx->RHSVelP,"RHS of flow solver");CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daFlow,&ctx->RHSVelPpre);CHKERRQ(ierr);
	ierr = VecSet(ctx->RHSVelPpre,0.);CHKERRQ(ierr);
	ierr = KSPCreate(PETSC_COMM_WORLD,&ctx->kspVelP);CHKERRQ(ierr);	
	ierr = KSPSetTolerances(ctx->kspVelP,1.e-6,1.e-6,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSetOperators(ctx->kspVelP,ctx->KVelP,ctx->KVelP,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = KSPSetInitialGuessNonzero(ctx->kspVelP,PETSC_TRUE);CHKERRQ(ierr);
	ierr = KSPAppendOptionsPrefix(ctx->kspVelP,"VelP_");CHKERRQ(ierr);
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
	 //	ierr = PetscOptionsSetValue("-VelP_pc_ﬁeldsplit_schur_precondition","self");CHKERRQ(ierr);
	 ierr = PetscOptionsSetValue("-VelP_pc_fieldsplit_detect_saddle_point","true");CHKERRQ(ierr);
	 ierr = PCSetFromOptions(ctx->pcVelP);CHKERRQ(ierr);
	 ierr = KSPSetFromOptions(ctx->kspVelP);CHKERRQ(ierr);
	 */
	ierr = BCPInit(&ctx->bcP[0],ctx);
	ierr = BCQInit(&ctx->bcQ[0],ctx);
	ierr = GetFlowProp(&ctx->flowprop,ctx->units,ctx->resprop);CHKERRQ(ierr);
		//	ierr = ReSETFlowBC(&ctx->bcP[0],&ctx->bcQ[0],ctx->flowcase);CHKERRQ(ierr);	
	ierr = ReSETSourceTerms(ctx->Source,ctx->flowprop);		
	ierr = ReSETBoundaryTerms(ctx,fields);CHKERRQ(ierr);
		//	ierr = FormInitialSolution(fields->VelnPress,fields->FlowBCArray,&ctx->bcFlow[0],ctx);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MixedFEMFlowSolverFinalize"
extern PetscErrorCode MixedFEMFlowSolverFinalize(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	ierr = KSPDestroy(&ctx->kspVelP);CHKERRQ(ierr);
	ierr = MatDestroy(&ctx->KVelP);CHKERRQ(ierr);
	ierr = MatDestroy(&ctx->KVelPlhs);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->RHSVelP);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->RHSVelPpre);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "GetFlowProp"
extern PetscErrorCode GetFlowProp(FlowProp *flowprop,FlowUnit flowunit,ResProp resprop)
{
	PetscFunctionBegin;
	flowprop->theta    = 1.;						  /*Time paramter*/
	flowprop->timestepsize = 1;				    	/*Time step size	*/
	flowprop->M_inv = 0.;
	switch (flowunit) {
		case UnitaryUnits:
			flowprop->mu    = 1.;                     /*viscosity in cp*/
			flowprop->rho   = 1.;                     /*density in lb/ft^3*/
			flowprop->cf    = 1.;                     /*compressibility in psi^{-1}*/
			flowprop->beta  = 1;                      /*flow rate conversion constant*/
			flowprop->gamma = 1;                      /*pressure conversion constant*/
			flowprop->alpha = 1;                      /*volume conversion constatnt*/
			flowprop->g[0]  = 0.;                     /*x-component of gravity. unit is ft/s^2*/
			flowprop->g[1]  = 0.;                     /*y-component of gravity. unit is ft/s^2*/
			flowprop->g[2]  = 0.;                     /* 32.17;									/ *z-component of gravity. unit is ft/s^2* / */
			flowprop->Cp = 1.;                      /*Liquid specific heat capacity*/
			break;
		case FieldUnits: 
			flowprop->mu = resprop.visc;				  /* viscosity in cp */ 
			flowprop->rho = 62.428*resprop.fdens;	/* density in lb/ft^3 */ 
			flowprop->cf = resprop.rock_comp;						/* compressibility in psi^{-1} */ 
			flowprop->beta = 1.127;								/* flow rate conversion constant */ 
			flowprop->gamma = 2.158e-4;						/* pressue conversion constant */ 
			flowprop->alpha = 5.615;							/* volume conversion constatnt*/ 
			flowprop->g[0] = 0.;									/* x-componenet of gravity. unit is ft/s^2 */ 
			flowprop->g[1] = 0.;									/* y-component of gravity. unit is ft/s^2 */ 
			flowprop->g[2] = 0.;//32.17;					/* z-component of gravity. unit is ft/s^2 */ 
			flowprop->Cp = 1.;                      /*Liquid specific heat capacity*/
			break; 
		case MetricUnits:
			flowprop->mu    = 0.001*resprop.visc;     /*viscosity in Pa.s*/
			flowprop->rho   = 1000*resprop.fdens;     /*density in kg/m^3*/
			flowprop->cf    = 1.450e-4*resprop.wat_comp;    /*compressibility in Pa^{-1}*/
			flowprop->beta  = 86.4e-6;                /*flow rate conversion constant*/
			flowprop->gamma = 1e-3;                   /*pressue conversion constant*/
			flowprop->alpha = 1;                      /*volume conversion constatnt*/
			flowprop->g[0]  = 0.;                     /*x-component of gravity. unit is m/s^2*/
			flowprop->g[1]  = 0.;                     /*y-component of gravity. unit is m/s^2*/
			flowprop->g[2]  = 9.81;                   /*z-component of gravity. unit is m/s^2*/
			flowprop->Cp = 1.;                      /*Liquid specific heat capacity*/
			break;
		default:
			SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: [%s] unknown FLOWCASE %i .\n",__FUNCT__,flowunit);
			break;
	}
	PetscFunctionReturn(0);
}
