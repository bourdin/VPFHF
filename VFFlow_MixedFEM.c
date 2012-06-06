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
/* #include "PetscFixes.h" */
#include "VFFlow_MixedFEM.h"

/*
 VFFlow_DarcyMixedFEMSteadyState
*/
#undef __FUNCT__
#define __FUNCT__ "VFFlow_DarcyMixedFEMSteadyState"
extern PetscErrorCode VFFlow_DarcyMixedFEMSteadyState(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode     ierr;
	PetscViewer        viewer;
	PetscInt           xs,xm,ys,ym,zs,zm;
	PetscInt           i,j,k,c,veldof = 3;
	PetscInt           its;
	KSPConvergedReason reason;
	PetscReal          ****VelnPress_array;
	PetscReal          ***Press_array;
	Vec                vec;
	PetscReal          ****vel_array;
	
	PetscFunctionBegin;
	ierr = VecDuplicate(ctx->RHSVelP,&vec);CHKERRQ(ierr);
	
	ierr = DMDAGetCorners(ctx->daFlow,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = FlowMatnVecAssemble(ctx->KVelP,ctx->RHSVelP,fields,ctx);CHKERRQ(ierr);
	ierr = KSPSolve(ctx->kspVelP,ctx->RHSVelP,fields->VelnPress);CHKERRQ(ierr);
	
	/*
	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"Matrix.txt",&viewer);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
	ierr = MatView(ctx->KVelP,viewer);CHKERRQ(ierr);
	
	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"RHS.txt",&viewer);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
	ierr = VecView(ctx->RHSVelP,viewer);CHKERRQ(ierr);
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
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "SETFlowBC"
extern PetscErrorCode SETFlowBC(FLOWBC *BC,FlowCases flowcase)
{
	PetscInt i,c;
	
	PetscFunctionBegin;
	for (i = 0; i < 6; i++) {
		for (c = 0; c < 4; c++) {
			BC[c].face[i] = NOBC;
		}
	}
	for (i = 0; i < 12; i++) {
		for (c = 0; c < 4; c++) {
			BC[c].edge[i] = NOBC;
		}
	}
	for (i = 0; i < 8; i++) {
		for (c = 0; c < 4; c++) {
			BC[c].vertex[i] = NOBC;
		}
	}
	switch (flowcase) {
		case ALLPRESSUREBC:
			for (i = 0; i < 6; i++) {
				BC[3].face[i] = PRESSURE;
			}
			break;
		case ALLNORMALFLOWBC:			
			BC[0].face[X0] = VELOCITY;
			BC[0].face[X1] = VELOCITY;
			BC[1].face[Y0] = VELOCITY;
			BC[1].face[Y1] = VELOCITY;
			BC[2].face[Z0] = VELOCITY;
			BC[2].face[Z1] = VELOCITY;
			break;
		default:
			SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: [%s] unknown FLOWCASE %i .\n",__FUNCT__,flowcase);
			break;
	}
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "VecApplyWellFlowBC"
extern PetscErrorCode VecApplyWellFlowBC(PetscReal *Ks_local,PetscReal ***source_array,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,VFCtx *ctx)
{
	PetscErrorCode ierr;
	PetscInt       i,j,k,l;
	PetscReal      beta_c;
	PetscReal      mu;
	PetscReal      *loc_source;
	PetscInt       eg;
	
	PetscFunctionBegin;
	beta_c = ctx->flowprop.beta;
	mu     = ctx->flowprop.mu;
	ierr   = PetscMalloc(e->ng*sizeof(PetscReal),&loc_source);CHKERRQ(ierr);
	
	for (eg = 0; eg < e->ng; eg++) loc_source[eg] = 0.;
	for (k = 0; k < e->nphiz; k++) {
		for (j = 0; j < e->nphiy; j++) {
			for (i = 0; i < e->nphix; i++) {
				for (eg = 0; eg < e->ng; eg++) {
					loc_source[eg] += source_array[ek+k][ej+j][ei+i]*e->phi[k][j][i][eg];
				}
				;
			}
		}
	}
	for (eg = 0; eg < e->ng; eg++)
		for (l = 0,k = 0; k < e->nphiz; k++) {
			for (j = 0; j < e->nphiy; j++) {
				for (i = 0; i < e->nphix; i++,l++) {
					Ks_local[l] = 0.;
					for (eg = 0; eg < e->ng; eg++) {
						Ks_local[l] += -loc_source[eg]*e->phi[k][j][i][eg]*e->weight[eg];
					}
					;
				}
			}
		}
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "VecApplyFlowBC"
extern PetscErrorCode VecApplyFlowBC(Vec RHS,FLOWBC *BC,VFCtx *ctx, PetscReal ****UnPre_array)
{
	PetscErrorCode ierr;
	PetscInt       xs,xm,nx;
	PetscInt       ys,ym,ny;
	PetscInt       zs,zm,nz;
	PetscInt       dim,dof;
	PetscInt       i,j,k,c;
	DM             da;
	PetscReal      ****RHS_array;
	PetscReal      press;
	PetscReal      vel;
	PetscReal      hx,hy,hz;
	
	PetscFunctionBegin;
	ierr = PetscObjectQuery((PetscObject)RHS,"DM",(PetscObject*)&da);CHKERRQ(ierr);
	if (!da) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Vector not generated from a DMDA");
	
	ierr = DMDAGetInfo(da,&dim,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,&dof,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(da,RHS,&RHS_array);CHKERRQ(ierr);
	hx   = 1./(nx-1);hy = 1./(ny-1);hz = 1./(nz-1);
	for (c = 0; c < dof; c++) {
		if (xs == 0) {
			i = 0;
			if (BC[c].face[X0] == PRESSURE) {
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
			if (BC[c].face[X0] == VELOCITY) {
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
		}
		if (xs+xm == nx) {
			i = nx-1;
			if (BC[c].face[X1] == PRESSURE) {
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
			if (BC[c].face[X1] == VELOCITY) {
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
		}
		if (ys == 0) {
			j = 0;
			if (BC[c].face[Y0] == PRESSURE) {
				for (k = zs; k < zs+zm; k++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
			if (BC[c].face[Y0] == VELOCITY) {
				for (k = zs; k < zs+zm; k++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
		}
		if (ys+ym == ny) {
			j = ny-1;
			if (BC[c].face[Y1] == PRESSURE) {
				for (k = zs; k < zs+zm; k++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
			if (BC[c].face[Y1] == VELOCITY) {
				for (k = zs; k < zs+zm; k++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];

					}
				}
			}
		}
		if (zs == 0) {
			k = 0;
			if (BC[c].face[Z0] == PRESSURE) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
			if (BC[c].face[Z0] == VELOCITY) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
		}
		if (zs+zm == nz) {
			k = nz-1;
			if (BC[c].face[Z1] == PRESSURE) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
			if (BC[c].face[Z1] == VELOCITY) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
		}
		if (xs == 0 && zs == 0) {
			k = 0;i = 0;
			if (BC[c].edge[X0Z0] == PRESSURE) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
			if (BC[c].edge[X0Z0] == VELOCITY) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && zs == 0) {
			k = 0;i = nx-1;
			if (BC[c].edge[X1Z0] == PRESSURE) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
			if (BC[c].edge[X1Z0] == VELOCITY) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
		if (ys == 0 && zs == 0) {
			k = 0;j = 0;
			if (BC[c].edge[Y0Z0] == PRESSURE) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
			if (BC[c].edge[Y0Z0] == VELOCITY) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
		if (ys+ym == ny && zs == 0) {
			k = 0;j = 0;
			if (BC[c].edge[Y1Z0] == PRESSURE) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
			if (BC[c].edge[Y1Z0] == VELOCITY) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && zs+zm == nz) {
			k = nz-1;i = 0;
			if (BC[c].edge[X0Z1] == PRESSURE) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
			if (BC[c].edge[X0Z1] == VELOCITY) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && zs+zm == nz) {
			k = nz-1;i = nx-1;
			if (BC[c].edge[X1Z1] == PRESSURE) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
			if (BC[c].edge[X1Z1] == VELOCITY) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
		if (ys == 0 && zs+zm == nz) {
			k = nz-1;j = 0;
			if (BC[c].edge[Y0Z1] == PRESSURE) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
			if (BC[c].edge[Y0Z1] == VELOCITY) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
		if (ys+ym == ny && zs+zm == nz) {
			k = nz-1;j = ny-1;
			if (BC[c].edge[Y1Z1] == PRESSURE) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
			if (BC[c].edge[Y1Z1] == VELOCITY) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && ys == 0) {
			j = 0;i = 0;
			if (BC[c].edge[X0Y0] == PRESSURE) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
			if (BC[c].edge[X0Y0] == VELOCITY) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && ys+ym == ny) {
			j = ny-1;i = 0;
			if (BC[c].edge[X0Y1] == PRESSURE) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
			if (BC[c].edge[X0Y1] == VELOCITY) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && ys == 0) {
			j = 0;i = nx-1;
			if (BC[c].edge[X1Y0] == PRESSURE) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
			if (BC[c].edge[X1Y0] == VELOCITY) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && ys+ym == ny) {
			j = ny-1;i = nx-1;
			if (BC[c].edge[X1Y1] == PRESSURE) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
			if (BC[c].edge[X1Y1] == VELOCITY) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && ys == 0 && zs == 0) {
			k = 0;j = 0;i = 0;
			if (BC[c].vertex[X0Y0Z0] == PRESSURE) {
				RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
			}
			if (BC[c].vertex[X0Y0Z0] == VELOCITY) {
				RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys == 0 && zs == 0) {
			k = 0;j = 0;i = nx-1;
			if (BC[c].vertex[X1Y0Z0] == PRESSURE) {
				RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
			}
			if (BC[c].vertex[X1Y0Z0] == VELOCITY) {
				RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
			}
		}
		if (xs == 0 && ys+ym == ny && zs == 0) {
			k = 0;j = ny-1;i = 0;
			if (BC[c].vertex[X0Y1Z0] == PRESSURE) {
				RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
			}
			if (BC[c].vertex[X0Y1Z0] == VELOCITY) {
				RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys+ym == ny && zs == 0) {
			k = 0;j = ny-1;i = nx-1;
			if (BC[c].vertex[X1Y1Z0] == PRESSURE) {
				RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
			}
			if (BC[c].vertex[X1Y1Z0] == VELOCITY) {
				RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
			}
		}
		if (xs == 0 && ys == 0 && zs+zm == nz) {
			k = nz-1;j = 0;i = 0;
			if (BC[c].vertex[X0Y0Z1] == PRESSURE) {
				RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
			}
			if (BC[c].vertex[X0Y0Z1] == VELOCITY) {
				RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys == 0 && zs+zm == nz) {
			k = nz-1;j = 0;i = nx-1;
			if (BC[c].vertex[X1Y0Z1] == PRESSURE) {
				RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
			}
			if (BC[c].vertex[X1Y0Z1] == VELOCITY) {
				RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
			}
		}
		if (xs == 0 && ys+ym == ny && zs+zm == nz) {
			k = nz-1;j = ny-1;i = 0;
			if (BC[c].vertex[X0Y1Z1] == PRESSURE) {
				RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
			}
			if (BC[c].vertex[X0Y1Z1] == VELOCITY) {
				RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys+ym == ny && zs+zm == nz) {
			k = nz-1;j = ny-1;i = nx-1;
			if (BC[c].vertex[X1Y1Z1] == PRESSURE) {
				RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
			}
			if (BC[c].vertex[X1Y1Z1] == VELOCITY) {
				RHS_array[k][j][i][c] = UnPre_array[k][j][i][c];
			}
		}
	}
	ierr = DMDAVecRestoreArrayDOF(da,RHS,&RHS_array);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatApplyFlowBC"
extern PetscErrorCode MatApplyFlowBC(Mat K,DM da,FLOWBC *BC)
{
	PetscErrorCode ierr;
	PetscInt       xs,xm,nx;
	PetscInt       ys,ym,ny;
	PetscInt       zs,zm,nz;
	PetscInt       dim,dof;
	PetscReal      one = 1.;
	PetscInt       i,j,k,c;
	MatStencil     *row;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetInfo(da,&dim,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,&dof,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = PetscMalloc(1*sizeof(MatStencil),&row);CHKERRQ(ierr);
	
	for (c = 0; c < dof; c++) {
		/* 1. */
		if (xs == 0 && BC[c].face[X0] != NOBC) {
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					row[0].i = 0;row[0].j = j;row[0].k = k;row[0].c = c;
					ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
				}
			}
		}
		/* 2. */
		if (xs+xm == nx && BC[c].face[X1] != NOBC) {
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					row[0].i = nx-1;row[0].j = j;row[0].k = k;row[0].c = c;
					ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
				}
			}
		}
		/* 3 */
		if (ys == 0 && BC[c].face[Y0] != NOBC) {
			for (k = zs; k < zs+zm; k++) {
				for (i = xs; i < xs+xm; i++) {
					row[0].i = i;row[0].j = 0;row[0].k = k;row[0].c = c;
					ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
				}
			}
		}
		/* 4 */
		if (ys+ym == ny && BC[c].face[Y1] != NOBC) {
			for (k = zs; k < zs+zm; k++) {
				for (i = xs; i < xs+xm; i++) {
					row[0].i = i;row[0].j = ny-1;row[0].k = k;row[0].c = c;
					ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
				}
			}
		}
		/* 5 */
		if (zs == 0 && BC[c].face[Z0] != NOBC) {
			for (j = ys; j < ys+ym; j++) {
				for (i = xs; i < xs+xm; i++) {
					row[0].i = i;row[0].j = j;row[0].k = 0;row[0].c = c;
					ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
				}
			}
		}
		/* 6 */
		if (zs+zm == nz && BC[c].face[Z1] != NOBC) {
			for (j = ys; j < ys+ym; j++) {
				for (i = xs; i < xs+xm; i++) {
					row[0].i = i;row[0].j = j;row[0].k = nz-1;row[0].c = c;
					ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
				}
			}
		}
		/* 1. */
		if (xs == 0 && zs == 0 && BC[c].edge[X0Z0] != NOBC) {
			for (j = ys; j < ys+ym; j++) {
				row[0].i = 0;row[0].j = j;row[0].k = 0;row[0].c = c;
				ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
			}
		}
		/* 2. */
		if (xs+xm == nx && zs == 0 && BC[c].edge[X1Z0] != NOBC) {
			for (j = ys; j < ys+ym; j++) {
				row[0].i = nx-1;row[0].j = j;row[0].k = 0;row[0].c = c;
				ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
			}
		}
		/* 3.		Realignment of partitioning to avoid corner */
		if (ys == 0 && zs == 0 && BC[c].edge[Y0Z0] != NOBC) {
			for (i = xs; i < xs+xm; i++) {
				row[0].i = i;row[0].j = 0;row[0].k = 0;row[0].c = c;
				ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
			}
		}
		/* 4.		Realignment of partitioning to avoid corner */
		if (ys+ym == ny && zs == 0 && BC[c].edge[Y1Z0] != NOBC) {
			for (i = xs; i < xs+xm; i++) {
				row[0].i = i;row[0].j = ny-1;row[0].k = 0;row[0].c = c;
				ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
			}
		}
		/* 5. */
		if (xs == 0 && zs+zm == nz && BC[c].edge[X0Z1] != NOBC) {
			for (j = ys; j < ys+ym; j++) {
				row[0].i = 0;row[0].j = j;row[0].k = nz-1;row[0].c = c;
				ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
			}
		}
		/* 6. */
		if (xs+xm == nx && zs+zm == nz && BC[c].edge[X1Z1] != NOBC) {
			for (j = ys; j < ys+ym; j++) {
				row[0].i = nx-1;row[0].j = j;row[0].k = nz-1;row[0].c = c;
				ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
			}
		}
		/* 7. */
		if (ys == 0 && zs+zm == nz && BC[c].edge[Y0Z1] != NOBC) {
			for (i = xs; i < xs+xm; i++) {
				row[0].i = i;row[0].j = 0;row[0].k = nz-1;row[0].c = c;
				ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
			}
		}
		/* 8. */
		if (ys+ym == ny && zs+zm == nz && BC[c].edge[Y1Z1] != NOBC) {
			for (i = xs; i < xs+xm; i++) {
				row[0].i = i;row[0].j = ny-1;row[0].k = nz-1;row[0].c = c;
				ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
			}
		}
		/* 9. */
		if (xs == 0 && ys == 0 && BC[c].edge[X0Y0] != NOBC) {
			for (k = zs; k < zs+zm; k++) {
				row[0].i = 0;row[0].j = 0;row[0].k = k;row[0].c = c;
				ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
			}
		}
		/* 10. */
		if (xs == 0 && ys+ym == ny && BC[c].edge[X0Y1] != NOBC) {
			for (k = zs; k < zs+zm; k++) {
				row[0].i = 0;row[0].j = ny-1;row[0].k = k;row[0].c = c;
				ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
			}
		}
		/* 11. */
		if (xs+xm == nx && ys == 0 && BC[c].edge[X1Y0] != NOBC) {
			for (k = zs; k < zs+zm; k++) {
				row[0].i = nx-1;row[0].j = 0;row[0].k = k;row[0].c = c;
				ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
			}
		}
		/* 12. */
		if (xs+xm == nx && ys+ym == ny && BC[c].edge[X1Y1] != NOBC) {
			for (k = zs; k < zs+zm; k++) {
				row[0].i = nx-1;row[0].j = ny-1;row[0].k = k;row[0].c = c;
				ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
			}
		}
		if (xs == 0 && ys == 0 && zs == 0 && BC[c].vertex[X0Y0Z0] != NOBC) {
			row[0].i = 0;row[0].j = 0;row[0].k = 0;row[0].c = c;
			ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
		}
		if (xs == 0 && ys == 0 && zs+zm == nz && BC[c].vertex[X0Y0Z1] != NOBC) {
			row[0].i = 0;row[0].j = 0;row[0].k = nz-1;row[0].c = c;
			ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
		}
		if (xs == 0 && ys+ym == ny && zs == 0 && BC[c].vertex[X0Y1Z0] != NOBC) {
			row[0].i = 0;row[0].j = ny-1;row[0].k = 0;row[0].c = c;
			ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
		}
		if (xs == 0 && ys+ym == ny && zs+zm == nz && BC[c].vertex[X0Y1Z1] != NOBC) {
			row[0].i = 0;row[0].j = ny-1;row[0].k = nz-1;row[0].c = c;
			ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
		}
		if (xs+xm == nx && ys == 0 && zs == 0 && BC[c].vertex[X1Y0Z0] != NOBC) {
			row[0].i = nx-1;row[0].j = 0;row[0].k = 0;row[0].c = c;
			ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
		}
		if (xs+xm == nx && ys == 0 && zs+zm == nz && BC[c].vertex[X1Y0Z1] != NOBC) {
			row[0].i = nx-1;row[0].j = 0;row[0].k = nz-1;row[0].c = c;
			ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
		}
		if (xs+xm == nx && ys+ym == ny && zs == 0 && BC[c].vertex[X1Y1Z0] != NOBC) {
			row[0].i = nx-1;row[0].j = ny-1;row[0].k = 0;row[0].c = c;
			ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
		}
		if (xs+xm == nx && ys+ym == ny && zs+zm == nz && BC[c].vertex[X1Y1Z1] != NOBC) {
			row[0].i = nx-1;row[0].j = ny-1;row[0].k = nz-1;row[0].c = c;
			ierr     = MatZeroRowsStencil(K,1,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
		}
	}
	
	ierr = PetscFree(row);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "SETSourceTerms"
extern PetscErrorCode SETSourceTerms(Vec Src,FlowProp flowpropty)
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
				/*		source_array[ek][ej][ei] = beta_c/mu*cos(pi*ek*hz)*cos(pi*ej*hy)*cos(pi*ei*hx); */
				source_array[ek][ej][ei] = 4.*pi*pi*3.*beta_c/mu*sin(2.*pi*ek*hz)*sin(2.*pi*ej*hy)*sin(2.*pi*ei*hx);
			}
		}
	}
	ierr = DMDAVecRestoreArray(da,Src,&source_array);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "FlowMatnVecAssemble"
extern PetscErrorCode FlowMatnVecAssemble(Mat K,Vec RHS,VFFields * fields,VFCtx *ctx)
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
	PetscReal      *KA_local,*KB_local,*KD_local,*KBTrans_local;
	PetscReal      beta_c,mu,gx,gy,gz;
	PetscInt       nrow = ctx->e3D.nphix*ctx->e3D.nphiy*ctx->e3D.nphiz;
	MatStencil     *row,*row1;
	PetscReal      ***source_array;
	Vec            source_local;
	PetscReal      pi;
	PetscReal		 ****velnpre_array;
	Vec            velnpre_local;
	
	PetscFunctionBegin;
	beta_c = ctx->flowprop.beta;
	mu     = ctx->flowprop.mu;
	gx     = ctx->flowprop.g[0];
	gy     = ctx->flowprop.g[1];
	gz     = ctx->flowprop.g[2];
	
	ierr = DMDAGetInfo(ctx->daFlow,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daFlow,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	/* This line ensures that the number of cells is one less than the number of nodes. Force processing of cells to stop once the second to the last node is processed */
	ierr = MatZeroEntries(K);CHKERRQ(ierr);
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
	
	ierr = DMGetLocalVector(ctx->daFlow,&velnpre_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daFlow,fields->FlowBCArray,INSERT_VALUES,velnpre_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daFlow,fields->FlowBCArray,INSERT_VALUES,velnpre_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daFlow,velnpre_local,&velnpre_array);CHKERRQ(ierr); 
	
	ierr = PetscMalloc7(nrow*nrow,PetscReal,&KA_local,
						nrow*nrow,PetscReal,&KB_local,
						nrow*nrow,PetscReal,&KD_local,
						nrow*nrow,PetscReal,&KBTrans_local,
						nrow,PetscReal,&RHS_local,
						nrow,MatStencil,&row,
						nrow,MatStencil,&row1);CHKERRQ(ierr);
	
	if (xs+xm == nx) xm--;
	if (ys+ym == ny) ym--;
	if (zs+zm == nz) zm--;
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx   = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy   = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz   = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				/*This computes the local contribution of the global A matrix*/
				ierr = FLow_MatA(KA_local,&ctx->e3D,ek,ej,ei);CHKERRQ(ierr);
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
					ierr = MatSetValuesStencil(K,nrow,row,nrow,row,KA_local,ADD_VALUES);CHKERRQ(ierr);
					ierr = MatSetValuesStencil(K,nrow,row1,nrow,row,KB_local,ADD_VALUES);CHKERRQ(ierr);
					ierr = MatSetValuesStencil(K,nrow,row,nrow,row1,KBTrans_local,ADD_VALUES);CHKERRQ(ierr);
				}
				ierr = FLow_MatD(KD_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,perm_array);CHKERRQ(ierr);
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							row[l].i = ei+i;row[l].j = ej+j;row[l].k = ek+k;row[l].c = 3;
						}
					}
				}
				ierr = MatSetValuesStencil(K,nrow,row,nrow,row,KD_local,ADD_VALUES);CHKERRQ(ierr);
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
				ierr = VecApplyWellFlowBC(RHS_local,source_array,&ctx->e3D,ek,ej,ei,ctx);
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							RHS_array[ek+k][ej+j][ei+i][3] += RHS_local[l];
						}
					}
				}
			}
		}
	}
	if (xs+xm == nx) xm++;
	if (ys+ym == ny) ym++;
	if (zs+zm == nz) zm++;
	
	ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatApplyFlowBC(K,ctx->daFlow,&ctx->bcFlow[0]);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(ctx->daScal,source_local,&source_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScal,&source_local);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArrayDOF(ctx->daFlow,RHS_localVec,&RHS_array);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daFlow,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daFlow,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daFlow,&RHS_localVec);CHKERRQ(ierr);
	/*
	PetscViewer				viewer;
	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"2RHS.txt",&viewer);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
	ierr = VecView(RHS,viewer);CHKERRQ(ierr);
	*/
	ierr = VecApplyFlowBC(RHS,&ctx->bcFlow[0],ctx,velnpre_array);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArrayDOF(ctx->daFlow,velnpre_local,&velnpre_array);CHKERRQ(ierr); 
	ierr = DMRestoreLocalVector(ctx->daFlow,&velnpre_local);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = PetscFree7(KA_local,KB_local,KD_local,KBTrans_local,RHS_local,row,row1);CHKERRQ(ierr);
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
#define __FUNCT__ "SETBoundaryTerms"
extern PetscErrorCode SETBoundaryTerms(VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode ierr;
	PetscReal		****perm_array;
	PetscReal		****velnpress_array;
	PetscInt		xs,xm,nx;
	PetscInt		ys,ym,ny;
	PetscInt		zs,zm,nz;
	PetscInt		dim, dof;
	PetscInt		ei,ej,ek,c;
	PetscReal		beta_c,gamma_c,rho,mu,pi,gx,gy,gz;
	PetscReal		hi,hj,hk;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(ctx->daFlow,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daFlow,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVFperm,fields->vfperm,&perm_array);CHKERRQ(ierr); 
	ierr = DMDAVecGetArrayDOF(ctx->daFlow,fields->FlowBCArray,&velnpress_array);CHKERRQ(ierr); 
	pi      = 6.*asin(0.5);
	beta_c = ctx->flowprop.beta;
	gamma_c = ctx->flowprop.gamma;
	rho     = ctx->flowprop.rho;
	mu     = ctx->flowprop.mu;
	gx = ctx->flowprop.g[0];
	gy = ctx->flowprop.g[1];
	gz = ctx->flowprop.g[2];
	hi = 1./(nx-1.);
	hj = 1./(ny-1.);
	hk = 1./(nz-1.);
	for (ek = zs; ek < zs + zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				for(c = 3; c < 6; c++){
					perm_array[ek][ej][ei][c] = 0.;
				}
			}
		}
	}
	for (ek = zs; ek < zs + zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				velnpress_array[ek][ej][ei][0] = -beta_c/mu*(2.*pi*cos(2.*pi*ei*hi)*sin(2.*pi*ej*hj)*sin(2.*pi*ek*hk)-gamma_c*rho*gx);
				velnpress_array[ek][ej][ei][1] = -beta_c/mu*(2.*pi*sin(2.*pi*ei*hi)*cos(2.*pi*ej*hj)*sin(2.*pi*ek*hk)-gamma_c*rho*gy);
				velnpress_array[ek][ej][ei][2] = -beta_c/mu*(2.*pi*sin(2.*pi*ei*hi)*sin(2.*pi*ej*hj)*cos(2.*pi*ek*hk)-gamma_c*rho*gz);
				velnpress_array[ek][ej][ei][3] = sin(2.*pi*ek*hk)*sin(2.*pi*ej*hj)*sin(2.*pi*ei*hi);
					//	printf("\nbc[2] = %f\n", velnpress_array[ek][ej][ei][2]);
			}
		}
	}
	ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,fields->vfperm,&perm_array);CHKERRQ(ierr);	
	ierr = DMDAVecRestoreArrayDOF(ctx->daFlow,fields->FlowBCArray,&velnpress_array);CHKERRQ(ierr);	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MixedFEMFlowSolverInitialize"
extern PetscErrorCode MixedFEMFlowSolverInitialize(VFCtx *ctx, VFFields *fields)
{
	PetscMPIInt    comm_size;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"","");CHKERRQ(ierr);
	{
		ctx->units    = UnitaryUnits;
		ierr          = PetscOptionsEnum("-flowunits","\n\tFlow solver","",FlowUnitName,(PetscEnum)ctx->units,(PetscEnum*)&ctx->units,PETSC_NULL);CHKERRQ(ierr);
		/*	ctx->flowcase = ALLPRESSUREBC; */
		ctx->flowcase = ALLNORMALFLOWBC;
		ierr          = PetscOptionsEnum("-flow boundary conditions","\n\tFlow solver","",FlowBC_Case,(PetscEnum)ctx->flowcase,(PetscEnum*)&ctx->flowcase,PETSC_NULL);CHKERRQ(ierr);
	}
	ierr = PetscOptionsEnd();CHKERRQ(ierr);
	
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&comm_size);CHKERRQ(ierr);
	if (comm_size == 1) {
		ierr = DMCreateMatrix(ctx->daFlow,MATSEQAIJ,&ctx->KVelP);CHKERRQ(ierr);
	} else {
		ierr = DMCreateMatrix(ctx->daFlow,MATMPIAIJ,&ctx->KVelP);CHKERRQ(ierr);
	}
	ierr = MatSetOption(ctx->KVelP,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daFlow,&ctx->RHSVelP);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)ctx->RHSVelP,"RHS of flow solver");CHKERRQ(ierr);
	
	ierr = KSPCreate(PETSC_COMM_WORLD,&ctx->kspVelP);CHKERRQ(ierr);
	
	ierr = KSPSetTolerances(ctx->kspVelP,1.e-8,1.e-8,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSetOperators(ctx->kspVelP,ctx->KVelP,ctx->KVelP,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = KSPSetInitialGuessNonzero(ctx->kspVelP,PETSC_TRUE);CHKERRQ(ierr);
	ierr = KSPAppendOptionsPrefix(ctx->kspVelP,"VelP_");CHKERRQ(ierr);
	ierr = KSPSetType(ctx->kspVelP,KSPBICG);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ctx->kspVelP);CHKERRQ(ierr);
	ierr = KSPGetPC(ctx->kspVelP,&ctx->pcVelP);CHKERRQ(ierr);
	ierr = PCSetType(ctx->pcVelP,PCJACOBI);CHKERRQ(ierr);
	ierr = PCSetFromOptions(ctx->pcVelP);CHKERRQ(ierr);
	
	ierr = GetFlowProp(&ctx->flowprop,ctx->units,ctx->resprop);CHKERRQ(ierr);
	ierr = SETFlowBC(&ctx->bcFlow[0],ctx->flowcase);CHKERRQ(ierr);
	ierr = SETSourceTerms(ctx->Source,ctx->flowprop);
	
	ierr = SETBoundaryTerms(ctx,fields);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MixedFEMFlowSolverFinalize"
extern PetscErrorCode MixedFEMFlowSolverFinalize(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = DMDestroy(&ctx->daFlow);CHKERRQ(ierr);
	ierr = VecDestroy(&fields->VelnPress);CHKERRQ(ierr);
	ierr = KSPDestroy(&ctx->kspVelP);CHKERRQ(ierr);
	ierr = MatDestroy(&ctx->KVelP);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->RHSVelP);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->Source);CHKERRQ(ierr);
	ierr = VecDestroy(&fields->FlowBCArray);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "GetFlowProp"
extern PetscErrorCode GetFlowProp(FlowProp *flowprop,FlowUnit flowunit,ResProp resprop)
{
	PetscFunctionBegin;
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
			break;
			
		default:
			SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_USER,"ERROR: [%s] unknown FLOWCASE %i .\n",__FUNCT__,flowunit);
			break;
	}
	PetscFunctionReturn(0);
}
