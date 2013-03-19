/*
   VFFlow.c
   Generic interface to flow solvers

     (c) 2010-2012 Blaise Bourdin, LSU. bourdin@lsu.edu
                    Keita Yoshioka, Chevron ETC
*/
#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFFlow.h"
#include "VFFlow_FEM.h"
#include "VFFlow_Fake.h"
#include "VFFlow_KSPMixedFEM.h"
#include "VFFlow_SNESMixedFEM.h"
#include "VFFlow_TSMixedFEM.h"
#include "VFHeat_SNESFEM.h"



#undef __FUNCT__
#define __FUNCT__ "FlowSolverFinalize"
/* 
   FlowSolverFinalize
   Keita Yoshioka yoshk@chevron.com
*/
extern PetscErrorCode FlowSolverFinalize(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	switch (ctx->flowsolver) {
		case FLOWSOLVER_KSPMIXEDFEM:       
			ierr = MixedFEMFlowSolverFinalize(ctx,fields);CHKERRQ(ierr);
			break;
		case FLOWSOLVER_TSMIXEDFEM:       
			ierr = MixedFEMTSFlowSolverFinalize(ctx,fields);CHKERRQ(ierr);
			break;
		case FLOWSOLVER_SNESMIXEDFEM:       
			ierr = MixedFEMSNESFlowSolverFinalize(ctx,fields);CHKERRQ(ierr);
			break;
		case FLOWSOLVER_FEM:
			break; 
		case FLOWSOLVER_TS:
		    ierr = FEMTSFlowSolverFinalize(ctx,fields);CHKERRQ(ierr);
		    break;
		case FLOWSOLVER_SNES:
		    ierr = FEMSNESFlowSolverFinalize(ctx,fields);CHKERRQ(ierr);
		    break;
		case FLOWSOLVER_FAKE:
			break; 
		case FLOWSOLVER_READFROMFILES:
			break;
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FlowSolverInitialize"
extern PetscErrorCode FlowSolverInitialize(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	switch (ctx->flowsolver) {
		case FLOWSOLVER_KSPMIXEDFEM:       
			ierr = MixedFEMFlowSolverInitialize(ctx,fields);CHKERRQ(ierr);
			break;
		case FLOWSOLVER_TSMIXEDFEM:       
			ierr = MixedFEMTSFlowSolverInitialize(ctx,fields);CHKERRQ(ierr);
			break;
		case FLOWSOLVER_SNESMIXEDFEM:       
			ierr = MixedFEMSNESFlowSolverInitialize(ctx,fields);CHKERRQ(ierr);
			break;
		case FLOWSOLVER_FEM:
			break;
		case FLOWSOLVER_TS:
		    ierr = FEMTSFlowSolverInitialize(ctx,fields);CHKERRQ(ierr);
		    break; 
		case FLOWSOLVER_SNES:
		    ierr = FEMSNESFlowSolverInitialize(ctx,fields);CHKERRQ(ierr);
		    break;
		case FLOWSOLVER_FAKE:
			break; 
		case FLOWSOLVER_READFROMFILES:
			break;
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BCQInit"
extern PetscErrorCode BCQInit(BC *BCQ,VFCtx *ctx)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = BCInit(BCQ,3);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BCPInit"
extern PetscErrorCode BCPInit(BC *BCP,VFCtx *ctx)
{
  PetscErrorCode ierr;
 
  PetscFunctionBegin;
  ierr = BCInit(BCP,1);CHKERRQ(ierr);


//  When positive value is given, boundary condiction has a numerical value
/*
  if(ctx->BCpres[0] > -1.e-8) BCP[0].face[X0] = VALUE;
  if(ctx->BCpres[1] > -1.e-8) BCP[0].face[X1] = VALUE;
  if(ctx->BCpres[2] > -1.e-8) BCP[0].face[Y0] = VALUE;
  if(ctx->BCpres[3] > -1.e-8) BCP[0].face[Y1] = VALUE;
  if(ctx->BCpres[4] > -1.e-8) BCP[0].face[Z0] = VALUE;
  if(ctx->BCpres[5] > -1.e-8) BCP[0].face[Z1] = VALUE;
 */
  PetscFunctionReturn(0);
}  
  
#undef __FUNCT__
#define __FUNCT__ "SETBoundaryTerms_P"
extern PetscErrorCode SETBoundaryTerms_P(VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode ierr;
	PetscReal		***pressure_array;
	PetscInt		xs,xm,nx;
	PetscInt		ys,ym,ny;
	PetscInt		zs,zm,nz;
	PetscInt		dim, dof;
	PetscInt		i,j,k;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(ctx->daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,ctx->PresBCArray,&pressure_array);CHKERRQ(ierr); 

    for (k = zs; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
            for (i = xs; i < xs+xm; i++) {
		        // x == 0
                if (i == 0 && ctx->bcP[0].face[X0] == VALUE) {
                    pressure_array[k][j][i] = ctx->BCpres[X0];			
                }
		        // x == nx-1
                else if (i == nx-1 && ctx->bcP[0].face[X1] == VALUE) {
                    pressure_array[k][j][i] = ctx->BCpres[X1];
                }
		        // y == 0
                else if (j == 0 && ctx->bcP[0].face[Y0] == VALUE) {
                    pressure_array[k][j][i] = ctx->BCpres[Y0];
                }
                //  y == ny-1
                else if (j == ny-1 && ctx->bcP[0].face[Y1] == VALUE) {
                    pressure_array[k][j][i] = ctx->BCpres[Y1];
                }
                //  z == 0
                else if (k == 0 && ctx->bcP[0].face[Z0] == VALUE) {
                    pressure_array[k][j][i] = ctx->BCpres[Z0];
                }
                //  z == nz-1
                else if (k == nz-1 && ctx->bcP[0].face[Z1] == VALUE) {
                    pressure_array[k][j][i] = ctx->BCpres[Z1];
                }else  {
                    pressure_array[k][j][i] = ctx->resprop.Pinit;
                }
            }
        }
    }
	ierr = DMDAVecRestoreArray(ctx->daScal,ctx->PresBCArray,&pressure_array);CHKERRQ(ierr);	

	PetscFunctionReturn(0);
}

/*
  VFFlowTimeStep: Does one time step of the flow solver selected in ctx.flowsolver
*/
#undef __FUNCT__
#define __FUNCT__ "VFFlowTimeStep"
extern PetscErrorCode VFFlowTimeStep(VFCtx *ctx,VFFields *fields)
{
  char           filename[FILENAME_MAX];
  PetscViewer    viewer;
  PetscErrorCode     ierr;

  PetscFunctionBegin;
  switch (ctx->flowsolver) {
    case FLOWSOLVER_TS:
	  ierr = FlowFEMTSSolve(ctx,fields);CHKERRQ(ierr);
	  break;
    case FLOWSOLVER_FEM:       
//      ierr = VFFlow_FEM(ctx,fields);CHKERRQ(ierr);
      break;
	case FLOWSOLVER_SNES:
	  ierr = FlowFEMSNESSolve(ctx,fields);CHKERRQ(ierr);
      break;
    case FLOWSOLVER_KSPMIXEDFEM:
      ierr = MixedFlowFEMKSPSolve(ctx,fields);CHKERRQ(ierr);
      break;
	case FLOWSOLVER_TSMIXEDFEM:
		ierr = MixedFlowFEMTSSolve(ctx,fields);CHKERRQ(ierr);
		break;
	case FLOWSOLVER_SNESMIXEDFEM:
		ierr = MixedFlowFEMSNESSolve(ctx,fields);CHKERRQ(ierr);
		break;
	case FLOWSOLVER_FAKE:
		  ierr = VFFlow_Fake(ctx,fields);CHKERRQ(ierr);
		  break; 
    case FLOWSOLVER_READFROMFILES:
      //ierr = PetscLogStagePush(ctx->vflog.VF_IOStage);CHKERRQ(ierr);
      switch (ctx->fileformat) {
        case FILEFORMAT_HDF5:
#ifdef PETSC_HAVE_HDF5
          ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
          /*
            ierr = VecLoad(viewer,fields->theta);CHKERRQ(ierr);
          */
          ierr = VecLoad(fields->pressure,viewer);CHKERRQ(ierr);
          ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);    
#endif
          break;
        case FILEFORMAT_BIN:
          SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Reading from binary files not implemented yet");
          break;
      }
      //ierr = PetscLogStagePop();CHKERRQ(ierr);
	// eventually replace FLOWSOLVER_FEM with this after confirming it provides the same results

  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecApplyPressureBC_FEM"
extern PetscErrorCode VecApplyPressureBC_FEM(Vec RHS,Vec BCF,BC *BC)
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

#undef __FUNCT__
#define __FUNCT__ "MatApplyPressureBC_FEM"
extern PetscErrorCode MatApplyPressureBC_FEM(Mat K,Mat M,BC *bcP)
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
	ierr = MatZeroRowsStencil(M,numBC,row,zero,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscFree(row);CHKERRQ(ierr);	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFFlow_FEM_MatKPAssembly3D_local"
/*
   VFFlow_FEM_MatPAssembly3D_local
*/

extern PetscErrorCode VFFlow_FEM_MatKPAssembly3D_local(PetscReal *Mat_local,FlowProp *flowprop,PetscReal ****perm_array, PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e)
{
  PetscInt       g,i1,i2,j1,j2,k1,k2,l;
  PetscReal      rho,relk,mu;
  PetscReal      kxx,kyy,kzz,kxy,kxz,kyz;
  PetscReal      DCoef_P;
  PetscErrorCode ierr;

  PetscFunctionBegin;
/*
  The following properties should be changed to a function of pressure and temperature (and saturation for multi-phase)
*/
  rho     = flowprop->rho;
  mu      = flowprop->mu;
  DCoef_P = -1*rho/mu;
/*
  Permeability should be associated with each element later
*/
  kxx = perm_array[ek][ej][ei][0];
  kyy = perm_array[ek][ej][ei][1];
  kzz = perm_array[ek][ej][ei][2];
  kxy = perm_array[ek][ej][ei][3];
  kxz = perm_array[ek][ej][ei][4];
  kyz = perm_array[ek][ej][ei][5];

  for (l = 0,k1 = 0; k1 < e->nphiz; k1++) {
    for (j1 = 0; j1 < e->nphiy; j1++) {
      for (i1 = 0; i1 < e->nphix; i1++) {
        for (k2 = 0; k2 < e->nphiz; k2++) {
          for (j2 = 0; j2 < e->nphiy; j2++) {
            for (i2 = 0; i2 < e->nphix; i2++,l++) {
              Mat_local[l] = 0.;
              for (g = 0; g < e->ng; g++) {
                Mat_local[l] += e->weight[g]*DCoef_P*(kxx*e->dphi[k1][j1][i1][0][g]*e->dphi[k2][j2][i2][0][g]
                                                      +kyy*e->dphi[k1][j1][i1][1][g]*e->dphi[k2][j2][i2][1][g]
                                                      +kzz*e->dphi[k1][j1][i1][2][g]*e->dphi[k2][j2][i2][2][g]
                                                      +kxy*(e->dphi[k1][j1][i1][1][g]*e->dphi[k2][j2][i2][0][g]+e->dphi[k1][j1][i1][0][g]*e->dphi[k2][j2][i2][1][g])
                                                      +kxz*(e->dphi[k1][j1][i1][2][g]*e->dphi[k2][j2][i2][0][g]+e->dphi[k1][j1][i1][0][g]*e->dphi[k2][j2][i2][2][g])
                                                      +kyz*(e->dphi[k1][j1][i1][2][g]*e->dphi[k2][j2][i2][1][g]+e->dphi[k1][j1][i1][1][g]*e->dphi[k2][j2][i2][2][g]));
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
#define __FUNCT__ "VFFlow_FEM_MatMPAssembly3D_local"
/*
   VFFlow_FEM_MassMatPAssembly3D_local
*/

extern PetscErrorCode VFFlow_FEM_MatMPAssembly3D_local(PetscReal *Mat_local,FlowProp *flowprop,PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e)
{
  PetscInt       g,i1,i2,j1,j2,k1,k2,l;
  PetscReal      ACoef_P;
  PetscErrorCode ierr;

  PetscFunctionBegin;
/*
  The following properties should be changed to a function of pressure and temperature (and saturation for multi-phase)
*/
  ACoef_P   = flowprop->M_inv;

  for (l = 0,k1 = 0; k1 < e->nphiz; k1++) {
    for (j1 = 0; j1 < e->nphiy; j1++) {
      for (i1 = 0; i1 < e->nphix; i1++) {
        for (k2 = 0; k2 < e->nphiz; k2++) {
          for (j2 = 0; j2 < e->nphiy; j2++) {
            for (i2 = 0; i2 < e->nphix; i2++,l++) {
              Mat_local[l] = 0.;
              for (g = 0; g < e->ng; g++) {
                Mat_local[l] += e->weight[g]*ACoef_P*e->phi[k1][j1][i1][g]*e->phi[k2][j2][i2][g];
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
#define __FUNCT__ "FEMSNESMonitor"
extern PetscErrorCode FEMSNESMonitor(SNES snes,PetscInt its,PetscReal fnorm,void* ptr)
{
	PetscErrorCode ierr;
	PetscReal      norm,vmax,vmin;
	MPI_Comm       comm;
	Vec            solution;


	PetscFunctionBegin;
    ierr = SNESGetSolution(snes, &solution);CHKERRQ(ierr);
	ierr = VecNorm(solution,NORM_1,&norm);CHKERRQ(ierr);
	ierr = VecMax(solution,PETSC_NULL,&vmax);CHKERRQ(ierr);
	ierr = VecMin(solution,PETSC_NULL,&vmin);CHKERRQ(ierr);
	ierr = PetscObjectGetComm((PetscObject)snes,&comm);CHKERRQ(ierr);
	ierr = PetscPrintf(comm,"snes_iter_step %D :solution norm = %G, max sol. value  = %G, min sol. value = %G\n",its,norm,vmax,vmin);CHKERRQ(ierr);
		
	PetscFunctionReturn(0);
}