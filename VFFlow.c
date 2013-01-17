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
		case FLOWSOLVER_DARCYMIXEDFEMSTEADYSTATE:       
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
		case FLOWSOLVER_DARCYMIXEDFEMSTEADYSTATE:       
			ierr = MixedFEMFlowSolverInitialize(ctx,fields);CHKERRQ(ierr);
			break;
		case FLOWSOLVER_TSMIXEDFEM:       
			ierr = MixedFEMTSFlowSolverInitialize(ctx,fields);CHKERRQ(ierr);
			break;
		case FLOWSOLVER_SNESMIXEDFEM:       
			ierr = MixedFEMSNESFlowSolverInitialize(ctx,fields);CHKERRQ(ierr);
			break;
		case FLOWSOLVER_TS:
		    ierr = FEMTSFlowSolverInitialize(ctx,fields);CHKERRQ(ierr);
		    break;
		case FLOWSOLVER_FEM:
			break; 
		case FLOWSOLVER_SNES:
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

/*
  When positive value is given, boundary condiction has a numerical value
*/
  if(ctx->BCpres[0] > -1.e-8) BCP[0].face[X0] = VALUE;
  if(ctx->BCpres[1] > -1.e-8) BCP[0].face[X1] = VALUE;
  if(ctx->BCpres[2] > -1.e-8) BCP[0].face[Y0] = VALUE;
  if(ctx->BCpres[3] > -1.e-8) BCP[0].face[Y1] = VALUE;
  if(ctx->BCpres[4] > -1.e-8) BCP[0].face[Z0] = VALUE;
  if(ctx->BCpres[5] > -1.e-8) BCP[0].face[Z1] = VALUE;

  BCP[0].face[X0] = VALUE;
  BCP[0].face[X1] = VALUE;
  BCP[0].face[Y0] = VALUE;
  BCP[0].face[Y1] = VALUE;
  BCP[0].face[Z0] = VALUE;
  BCP[0].face[Z1] = VALUE;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BCTInit"
/*
  BCTInit
  Keita Yoshioka yoshk@chevron.com
*/
extern PetscErrorCode BCTInit(BC *BCT,VFCtx *ctx)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = BCInit(BCT,1);CHKERRQ(ierr);
  
/*
  When positive value is given, boundary condiction has a numecial value
*/
  if(ctx->BCtheta[0] > -1.e-8) BCT[0].face[X0] = VALUE;
  if(ctx->BCtheta[1] > -1.e-8) BCT[0].face[X1] = VALUE;
  if(ctx->BCtheta[2] > -1.e-8) BCT[0].face[Y0] = VALUE;
  if(ctx->BCtheta[3] > -1.e-8) BCT[0].face[Y1] = VALUE;
  if(ctx->BCtheta[4] > -1.e-8) BCT[0].face[Z0] = VALUE;
  if(ctx->BCtheta[5] > -1.e-8) BCT[0].face[Z1] = VALUE;

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
	ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,fields->PresBCArray,&pressure_array);CHKERRQ(ierr); 

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
	ierr = DMDAVecRestoreArray(ctx->daScal,fields->PresBCArray,&pressure_array);CHKERRQ(ierr);	
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
      ierr = VFFlow_FEM(ctx,fields);CHKERRQ(ierr);
      break;
	case FLOWSOLVER_SNES:
//	  ierr = VFFlow_FEM(ctx,fields);CHKERRQ(ierr);
	  ierr = VFFlow_SNES_FEM(ctx,fields);CHKERRQ(ierr);
      break;
    case FLOWSOLVER_DARCYMIXEDFEMSTEADYSTATE:
      ierr = VFFlow_DarcyMixedFEMSteadyState(ctx,fields);CHKERRQ(ierr);
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



  
