/*
 VFFlow_TSFEM.c
 Finite element implementation of  Darcy flow 
 (c) 2010-2013 K. Yoshioka
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFFlow.h"
#include "VFFlow_TSFEM.h"
#include "VFFlow_KSPMixedFEM.h"

#undef __FUNCT__
#define __FUNCT__ "FEMTSFlowSolverInitialize"
extern PetscErrorCode FEMTSFlowSolverInitialize(VFCtx *ctx, VFFields *fields)
{
	PetscMPIInt    comm_size;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"","");CHKERRQ(ierr);
	{
		ctx->units    = UnitaryUnits;
		ierr          = PetscOptionsEnum("-flowunits","\n\tFlow solver","",FlowUnitName,(PetscEnum)ctx->units,(PetscEnum*)&ctx->units,PETSC_NULL);CHKERRQ(ierr);
	}
	ierr = PetscOptionsEnd();CHKERRQ(ierr);	
	
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&comm_size);CHKERRQ(ierr);
	if (comm_size == 1) {
		ierr = DMCreateMatrix(ctx->daScal,MATSEQAIJ,&ctx->KP);CHKERRQ(ierr);
		ierr = DMCreateMatrix(ctx->daScal,MATSEQAIJ,&ctx->KPlhs);CHKERRQ(ierr);		
		ierr = DMCreateMatrix(ctx->daScal,MATSEQAIJ,&ctx->JacP);CHKERRQ(ierr);
	} else {
		ierr = DMCreateMatrix(ctx->daScal,MATMPIAIJ,&ctx->KP);CHKERRQ(ierr);
		ierr = DMCreateMatrix(ctx->daScal,MATMPIAIJ,&ctx->KPlhs);CHKERRQ(ierr);
		ierr = DMCreateMatrix(ctx->daScal,MATMPIAIJ,&ctx->JacP);CHKERRQ(ierr);
	}
    ierr = MatZeroEntries(ctx->JacP);CHKERRQ(ierr);
    ierr = MatSetOption(ctx->KP,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
    ierr = MatSetOption(ctx->KPlhs,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
	ierr = MatSetOption(ctx->JacP,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(ctx->daScal,&ctx->RHSP);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) ctx->RHSP,"RHS of FEM Pressure");CHKERRQ(ierr);	
	ierr = DMCreateGlobalVector(ctx->daScal,&ctx->PFunct);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)ctx->PFunct,"RHS of FEM TS flow solver");CHKERRQ(ierr);

	ierr = TSCreate(PETSC_COMM_WORLD,&ctx->tsP);CHKERRQ(ierr);
	ierr = TSSetDM(ctx->tsP,ctx->daScal);CHKERRQ(ierr);
	ierr = TSSetProblemType(ctx->tsP,TS_LINEAR);CHKERRQ(ierr);
	ierr = TSSetType(ctx->tsP,TSBEULER);CHKERRQ(ierr);
		
	ierr = BCPInit(&ctx->bcP[0],ctx);

	ierr = GetFlowProp(&ctx->flowprop,ctx->units,ctx->resprop);CHKERRQ(ierr);
//	ierr = SETFlowBC(&ctx->bcP[0],&ctx->bcQ[0],ctx->flowcase);CHKERRQ(ierr);	// Currently BCpres is a PetscOption received from command line. Also done in test35
	ierr = SETSourceTerms(ctx->Source,ctx->flowprop);
	ierr = SETBoundaryTerms_P(ctx,fields);CHKERRQ(ierr); // Set fields.pressure according to BC & IC
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FEMTSFlowSolverFinalize"
extern PetscErrorCode FEMTSFlowSolverFinalize(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = MatDestroy(&ctx->KP);CHKERRQ(ierr);
	ierr = MatDestroy(&ctx->JacP);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->RHSP);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->PFunct);CHKERRQ(ierr);
	ierr = TSDestroy(&ctx->tsP);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FlowFEMTSSolve"
extern PetscErrorCode FlowFEMTSSolve(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode     ierr;
	KSPConvergedReason reason;	
	PetscReal          ***Press_array;
	PetscInt           i,j,k,c;
	PetscInt           xs,xm,ys,ym,zs,zm;
	PetscInt           its;

	PetscFunctionBegin;	
	ierr = DMCreateGlobalVector(ctx->daVFperm,&ctx->Perm);CHKERRQ(ierr);
	ierr = VecSet(ctx->Perm,0.0);CHKERRQ(ierr);
	ierr = VecCopy(fields->vfperm,ctx->Perm);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(ctx->daScal,&ctx->PresBC);CHKERRQ(ierr);
	ierr = VecSet(ctx->PresBC,0.0);CHKERRQ(ierr);
	ierr = VecCopy(fields->PresBCArray,ctx->PresBC);CHKERRQ(ierr);
	
	ierr = TSSetIFunction(ctx->tsP,PETSC_NULL,FormIFunction_P,ctx);CHKERRQ(ierr);
    ierr = TSSetIJacobian(ctx->tsP,ctx->JacP,ctx->JacP,FormIJacobian_P,ctx);CHKERRQ(ierr);
	ierr = TSSetRHSFunction(ctx->tsP,PETSC_NULL,FormFunction_P,ctx);CHKERRQ(ierr);
  
	ierr = TSSetSolution(ctx->tsP,fields->pressure);CHKERRQ(ierr);
	ierr = TSSetInitialTimeStep(ctx->tsP,0.0,ctx->timevalue);CHKERRQ(ierr);
    ierr = TSSetDuration(ctx->tsP,ctx->maxtimestep,ctx->maxtimevalue);CHKERRQ(ierr);
   
	ierr = TSMonitorSet(ctx->tsP,FEMTSMonitor,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = TSSetFromOptions(ctx->tsP);CHKERRQ(ierr);	

    ierr = TSSolve(ctx->tsP,fields->pressure,PETSC_NULL);CHKERRQ(ierr);
    ierr = TSGetTimeStepNumber(ctx->tsP,&ctx->timestep);CHKERRQ(ierr);
	
	ierr = VecDestroy(&ctx->Perm);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->PresBC);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "FEMTSMonitor"
extern PetscErrorCode FEMTSMonitor(TS ts,PetscInt timestep,PetscReal timevalue,Vec pressure,void* ptr)
{
	PetscErrorCode ierr;
	PetscReal      norm,vmax,vmin;
	MPI_Comm       comm;
	PetscViewer    viewer;
	PetscBool      drawcontours;


	PetscFunctionBegin;
/*	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"Solution.txt",&viewer);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
	ierr = VecView(pressure,viewer);CHKERRQ(ierr); */
	ierr = VecNorm(pressure,NORM_1,&norm);CHKERRQ(ierr);
	ierr = VecMax(pressure,PETSC_NULL,&vmax);CHKERRQ(ierr);
	ierr = VecMin(pressure,PETSC_NULL,&vmin);CHKERRQ(ierr);
	ierr = PetscObjectGetComm((PetscObject)ts,&comm);CHKERRQ(ierr);
	ierr = PetscPrintf(comm,"timestep %D: time %G, solution norm %G, max %G, min %G\n",timestep,timevalue,norm,vmax,vmin);CHKERRQ(ierr);
		
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormIFunction_P"
extern PetscErrorCode FormIFunction_P(TS ts,PetscReal t,Vec pressure,Vec pressuredot,Vec Func,void *user)
{
	PetscErrorCode ierr;
	VFCtx			*ctx=(VFCtx*)user;
	PetscViewer     viewer;
	
	PetscFunctionBegin;
	ierr = FormTSMatricesnVector_P(ctx->KP,ctx->KPlhs,ctx->RHSP,ctx);CHKERRQ(ierr);
/*	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"RHSVector.txt",&viewer);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
	ierr = VecView(ctx->RHSP,viewer);CHKERRQ(ierr);	*/
	ierr = VecSet(Func,0.0);CHKERRQ(ierr);
	ierr = MatMult(ctx->KP,pressure,Func);CHKERRQ(ierr);	
	ierr = MatMultAdd(ctx->KPlhs,pressuredot,Func,Func);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormFunction_P"
extern PetscErrorCode FormFunction_P(TS ts,PetscReal t,Vec vec1,Vec Func,void *user)
{
	PetscErrorCode ierr;
	VFCtx			*ctx=(VFCtx*)user;
	
	PetscFunctionBegin;
	ierr = VecSet(Func,0.0);CHKERRQ(ierr);
	ierr = VecCopy(ctx->RHSP,Func);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormIJacobian_P"
extern PetscErrorCode FormIJacobian_P(TS ts,PetscReal t,Vec pressure,Vec pressuredot,PetscReal shift,Mat *Jac,Mat *Jacpre,MatStructure *str,void *user)
{
	PetscErrorCode ierr;
	VFCtx          *ctx=(VFCtx*)user;
	PetscViewer    viewer;

	PetscFunctionBegin;
	*str = DIFFERENT_NONZERO_PATTERN;
	ierr = MatZeroEntries(*Jac);CHKERRQ(ierr);
	ierr = MatCopy(ctx->KP,*Jac,*str);
	ierr = MatAXPY(*Jac,shift,ctx->KPlhs,*str);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(*Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	if (*Jac != *Jacpre) {
		ierr = MatAssemblyBegin(*Jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(*Jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}	
	
/*	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"Jacobian.txt",&viewer);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
	ierr = MatView(*Jac,viewer);CHKERRQ(ierr);
*/
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormTSMatricesnVector_P"
extern PetscErrorCode FormTSMatricesnVector_P(Mat K,Mat Klhs,Vec RHS,VFCtx *ctx)
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
	Vec            RHS_localVec;	
	PetscReal      *RHS_local;
	Vec            perm_local;
	PetscReal      hx,hy,hz;
	PetscReal      *KA_local,*KD_local,*Klhs_local;
	PetscReal      beta_c,alpha_c,mu,gx,gy,gz;
	PetscInt       nrow = ctx->e3D.nphix*ctx->e3D.nphiy*ctx->e3D.nphiz;
	MatStencil     *row;
	PetscReal      ***source_array;
	Vec            source_local;
	PetscReal      M_inv;
	FACE           face;
	PetscReal      ***presbc_array;
	Vec            presbc_local;

	PetscFunctionBegin;
	M_inv   = ctx->flowprop.M_inv;
	beta_c  = ctx->flowprop.beta;
	alpha_c = ctx->flowprop.alpha;
	mu      = ctx->flowprop.mu;
	gx      = ctx->flowprop.g[0];
	gy      = ctx->flowprop.g[1];
	gz      = ctx->flowprop.g[2];	
	ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = MatZeroEntries(K);CHKERRQ(ierr);
	ierr = MatZeroEntries(Klhs);CHKERRQ(ierr);
	ierr = VecSet(RHS,0.);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScal,&RHS_localVec);CHKERRQ(ierr);
	ierr = VecSet(RHS_localVec,0.);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,RHS_localVec,&RHS_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScal,&source_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,ctx->Source,INSERT_VALUES,source_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,ctx->Source,INSERT_VALUES,source_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,source_local,&source_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVFperm,ctx->Perm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVFperm,ctx->Perm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);	

	ierr = DMGetLocalVector(ctx->daScal,&presbc_local);CHKERRQ(ierr);

	ierr = DMGlobalToLocalBegin(ctx->daScal,ctx->PresBC,INSERT_VALUES,presbc_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,ctx->PresBC,INSERT_VALUES,presbc_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daScal,presbc_local,&presbc_array);CHKERRQ(ierr); 
	
	ierr = PetscMalloc3(nrow*nrow,PetscReal,&KA_local,nrow*nrow,PetscReal,&KD_local,nrow*nrow,PetscReal,&Klhs_local);CHKERRQ(ierr);
	ierr = PetscMalloc2(nrow,PetscReal,&RHS_local,nrow,MatStencil,&row);CHKERRQ(ierr);
	
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx   = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy   = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz   = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				//This computes the local contribution of the global A matrix
				ierr = FLow_MatA(KA_local,&ctx->e3D,ek,ej,ei);CHKERRQ(ierr);
				for (l = 0; l < nrow*nrow; l++) {
					Klhs_local[l] = -2*M_inv*KA_local[l]/alpha_c;
				}
				ierr = FLow_MatD(KD_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,perm_array);CHKERRQ(ierr);
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							row[l].i = ei+i;row[l].j = ej+j;row[l].k = ek+k;row[l].c = 0;
						}
					}
				}
				ierr = MatSetValuesStencil(K,nrow,row,nrow,row,KD_local,ADD_VALUES);CHKERRQ(ierr);
				ierr = MatSetValuesStencil(Klhs,nrow,row,nrow,row,Klhs_local,ADD_VALUES);CHKERRQ(ierr);
				/*Assembling the right hand side vector g*/
				ierr = FLow_Vecg(RHS_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,perm_array);CHKERRQ(ierr);
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							RHS_array[ek+k][ej+j][ei+i] += RHS_local[l];
						}
					}
				}
				/* Adding source term */
				ierr = VecApplySourceTerms(RHS_local,source_array,&ctx->e3D,ek,ej,ei,ctx);
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							RHS_array[ek+k][ej+j][ei+i] += RHS_local[l];
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
    ierr = MatApplyTSPressureBC(K,Klhs,&ctx->bcP[0]);CHKERRQ(ierr);

    ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		
	ierr = DMDAVecRestoreArray(ctx->daScal,source_local,&source_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&source_local);CHKERRQ(ierr);	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,RHS_localVec,&RHS_array);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScal,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScal,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&RHS_localVec);CHKERRQ(ierr);	
	
    ierr = VecApplyTSPressureBC(RHS,ctx->PresBC,&ctx->bcP[0]);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daScal,presbc_local,&presbc_array);CHKERRQ(ierr); 
	ierr = DMRestoreLocalVector(ctx->daScal,&presbc_local);CHKERRQ(ierr);

	ierr = PetscFree3(KA_local,KD_local,Klhs_local);CHKERRQ(ierr);
	ierr = PetscFree2(RHS_local,row);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatApplyTSPressureBC"
extern PetscErrorCode MatApplyTSPressureBC(Mat K,Mat Klhs,BC *bcP)
{
	PetscErrorCode ierr;
	PetscInt       xs,xm,nx;
	PetscInt       ys,ym,ny;
	PetscInt       zs,zm,nz;
	PetscInt       i,j,k;
	MatStencil    *row;
	PetscReal      zero=0.0;
	PetscReal      one=1.;
	PetscInt       numBC=0,l=0,dim;
	DM             da;
	
	PetscFunctionBegin;
	ierr = PetscObjectQuery((PetscObject) K,"DM",(PetscObject *) &da); CHKERRQ(ierr);
	if (!da) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG," Matrix not generated from a DA");
	ierr = DMDAGetInfo(da,&dim,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

    if (xs == 0       && bcP[0].face[X0] == VALUE)             numBC += ym * zm;
    if (xs + xm == nx && bcP[0].face[X1] == VALUE)             numBC += ym * zm;
    if (ys == 0       && bcP[0].face[Y0] == VALUE)             numBC += xm * zm;
    if (ys + ym == ny && bcP[0].face[Y1] == VALUE)             numBC += xm * zm;
    if (zs == 0       && bcP[0].face[Z0] == VALUE && dim == 3) numBC += xm * ym;
    if (zs + zm == nz && bcP[0].face[Z1] == VALUE && dim == 3) numBC += xm * ym;

	ierr = PetscMalloc(numBC * sizeof(MatStencil),&row);CHKERRQ(ierr);
	/*
	 Create an array of rows to be zeroed out
	 */
	/*
	 i == 0
	 */
    if (xs == 0 && bcP[0].face[X0] == VALUE) {
        for (k = zs; k < zs + zm; k++) {
            for (j = ys; j < ys + ym; j++) {
                row[l].i = 0; row[l].j = j; row[l].k = k; row[l].c = 0; 
                l++;
            }
        }
    }
    /* 
     i == nx-1
    */
    if (xs + xm == nx && bcP[0].face[X1] == VALUE) {
        for (k = zs; k < zs + zm; k++) {
            for (j = ys; j < ys + ym; j++) {
                row[l].i = nx-1; row[l].j = j; row[l].k = k; row[l].c = 0; 
                l++;
            }
        }
    }
    /*
     y == 0
    */
    if (ys == 0 && bcP[0].face[Y0] == VALUE) {
        for (k = zs; k < zs + zm; k++) {
            for (i = xs; i < xs + xm; i++) {
                row[l].i = i; row[l].j = 0; row[l].k = k; row[l].c = 0; 
                l++;
            }
        }
    }
    /*
     y == ny-1
    */
    if (ys + ym == ny && bcP[0].face[Y1] == VALUE) {
        for (k = zs; k < zs + zm; k++) {
            for (i = xs; i < xs + xm; i++) {
                row[l].i = i; row[l].j = ny-1; row[l].k = k; row[l].c = 0; 
                l++;
            }
        }
    }
    /*
     z == 0
    */
    if (zs == 0 && bcP[0].face[Z0] == VALUE) {
        for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
                row[l].i = i; row[l].j = j; row[l].k = 0; row[l].c = 0; 
                l++;
            }
        }
    }
    /*
     z == nz-1
    */
    if (zs + zm == nz && bcP[0].face[Z1] == VALUE) {
        for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
                row[l].i = i; row[l].j = j; row[l].k = nz-1; row[l].c = 0; 
                l++;
            }
        }
    }
	
	ierr = MatZeroRowsStencil(K,numBC,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = MatZeroRowsStencil(Klhs,numBC,row,zero,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscFree(row);CHKERRQ(ierr);	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecApplyTSPressureBC"
extern PetscErrorCode VecApplyTSPressureBC(Vec RHS,Vec BCF,BC *BC)
{
  PetscErrorCode ierr;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       i,j,k,c;
  DM             da;
  PetscReal  ****RHS_array;
  PetscReal  ****BCF_array;
  PetscInt       dim,dof;
  
  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject) RHS,"DM",(PetscObject *) &da); CHKERRQ(ierr);
    if (!da) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Vector not generated from a DMDA");
  
  ierr = DMDAGetInfo(da,&dim,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                    &dof,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  
  if (dim == 2) {
    ierr = PetscMalloc(sizeof(PetscReal ***),&RHS_array);CHKERRQ(ierr);
    ierr = PetscMalloc(sizeof(PetscReal ***),&BCF_array);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(da,RHS,&RHS_array[0]);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(da,BCF,&BCF_array[0]);CHKERRQ(ierr);
  } else {
    ierr = DMDAVecGetArrayDOF(da,RHS,&RHS_array);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(da,BCF,&BCF_array);CHKERRQ(ierr);
  }
  
  for (c = 0;c < dof;c++) {
    /* 
      Faces
    */
    if (xs == 0) {
      /*
        x == 0
      */
      i = 0;
      if (BC[c].face[X0] == VALUE) {
        for (k = zs; k < zs + zm; k++) {
          for (j = ys; j < ys + ym; j++) {
            RHS_array[k][j][i][c] = BCF_array[k][j][i][c];
          }
        }
      }  
    }
    if (xs + xm == nx) {
      /*
        x == nx-1
      */
      i = nx-1;
      if (BC[c].face[X1] == VALUE) {
        for (k = zs; k < zs + zm; k++) {
          for (j = ys; j < ys + ym; j++) {
            RHS_array[k][j][i][c] = BCF_array[k][j][i][c];
          }
        }
      }  
    }
    if (ys == 0) {
      /*
        y == 0
      */
      j = 0;
      if (BC[c].face[Y0] == VALUE) {
        for (k = zs; k < zs + zm; k++) {
          for (i = xs; i < xs + xm; i++) {
            RHS_array[k][j][i][c] = BCF_array[k][j][i][c];
          }
        }
      }
    }
    if (ys + ym == ny) {
      /*
        y == ny-1
      */
      j = ny-1;
      if (BC[c].face[Y1] == VALUE) {
        for (k = zs; k < zs + zm; k++) {
          for (i = xs; i < xs + xm; i++) {
            RHS_array[k][j][i][c] = BCF_array[k][j][i][c];
          }
        }
      }
    }
    if (dim == 3) {
      if (zs == 0) {
        /*
          z == 0
        */
        k = 0;
        if (BC[c].face[Z0] == VALUE) {
          for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
              RHS_array[k][j][i][c] = BCF_array[k][j][i][c];
            }
          }
        }
      }
      if (zs + zm == nz) {
        /*
          z == nz-1
        */
        k = nz-1;
        if (BC[c].face[Z1] == VALUE) {
          for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
              RHS_array[k][j][i][c] = BCF_array[k][j][i][c];
            }
          }
        }
      }
    }
  }
  
  if (dim == 2) {
    ierr = DMDAVecRestoreArrayDOF(da,RHS,&RHS_array[0]);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(da,BCF,&BCF_array[0]);CHKERRQ(ierr);
    ierr = PetscFree(RHS_array);CHKERRQ(ierr);
    ierr = PetscFree(BCF_array);CHKERRQ(ierr);
  } else {
    ierr = DMDAVecRestoreArrayDOF(da,RHS,&RHS_array);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(da,BCF,&BCF_array);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
