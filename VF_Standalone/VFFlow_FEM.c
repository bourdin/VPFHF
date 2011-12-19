/*
   VFFlow_FEM.c
   A direct solver for the Stokes equation 
    
   (c) 2011 K. Yoshioka, CHEVRON ETC
*/

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"


#undef __FUNCT__
#define __FUNCT__ "VFFlow_FEM_MatPAssembly3D_local"
/* 
   VFFlow_FEM_MatPAssembly3D_local
*/

extern PetscErrorCode VFFlow_FEM_MatPAssembly3D_local(PetscReal *Mat_local,ResProp *resprop,PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e)
{
  PetscInt       g,i1,i2,j1,j2,k1,k2,l;
  PetscReal      fdens,relk,visc;
  PetscReal      kxx,kyy,kzz,kxy,kxz,kyz;
  PetscReal      DCoef_P;
  PetscErrorCode ierr;

  PetscFunctionBegin;
/*
  The following properties should be changed to a function of pressure and temperature (and saturation for multi-phase) 
*/
  fdens = resprop->fdens;
  relk = resprop->relk;
  visc = resprop->visc;
  DCoef_P = fdens*relk/visc;
/*
  Permeability should be associated with each element later
*/
  kxx = resprop->perm;
  kyy = resprop->perm;
  kzz = resprop->perm;
  kxy = 0.*resprop->perm;
  kxz = 0.*resprop->perm;
  kyz = 0.*resprop->perm;
  
  for (l = 0,k1 = 0; k1 < e->nphiz; k1++) {
    for (j1 = 0; j1 < e->nphiy; j1++) {
      for (i1 = 0; i1 < e->nphix; i1++) {
        for (k2 = 0; k2 < e->nphiz; k2++) {
          for (j2 = 0; j2 < e->nphiy; j2++) {
            for (i2 = 0; i2 < e->nphix; i2++,l++) {
              Mat_local[l] = 0.;
              for (g = 0; g < e->ng; g++) {
                Mat_local[l] += e->weight[g] * DCoef_P * ( kxx * e->dphi[k1][j1][i1][0][g] * e->dphi[k2][j2][i2][0][g]
                                                         + kyy * e->dphi[k1][j1][i1][1][g] * e->dphi[k2][j2][i2][1][g]
                                                         + kzz * e->dphi[k1][j1][i1][2][g] * e->dphi[k2][j2][i2][2][g]
                                                         + kxy * (e->dphi[k1][j1][i1][1][g] * e->dphi[k2][j2][i2][0][g] + e->dphi[k1][j1][i1][0][g]*e->dphi[k2][j2][i2][1][g])
                                                         + kxz * (e->dphi[k1][j1][i1][2][g] * e->dphi[k2][j2][i2][0][g] + e->dphi[k1][j1][i1][0][g]*e->dphi[k2][j2][i2][2][g])
                                                         + kyz * (e->dphi[k1][j1][i1][2][g] * e->dphi[k2][j2][i2][1][g] + e->dphi[k1][j1][i1][1][g]*e->dphi[k2][j2][i2][2][g]) );
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
#define __FUNCT__ "VFFlow_FEM_MatPAssembly3D"
/* 
   VFFlow_FEM_MatPAssembly3D
*/
extern PetscErrorCode VFFlow_FEM_MatPAssembly3D(Mat K,Vec RHS,VFFields *fields,VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       ei,ej,ek,i,j,k,l;
  PetscInt       nrow = ctx->e3D.nphix * ctx->e3D.nphiy * ctx->e3D.nphiz;
  Vec            RHS_localVec;
  PetscReal      ***RHS_array;
  PetscReal      *RHS_local;
  PetscReal      *K_local;
  MatStencil     *row;
  PetscReal      hx,hy,hz;
  PetscReal      ****coords_array;
 
  PetscFunctionBegin;
  ierr = PetscLogStagePush(ctx->vflog.VF_PAssemblyStage);CHKERRQ(ierr);
  ierr = DAGetInfo(ctx->daVect,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DAGetCorners(ctx->daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  if (xs+xm == nx) xm--;
  if (ys+ym == ny) ym--;
  if (zs+zm == nz) zm--;

  ierr = MatZeroEntries(K);CHKERRQ(ierr);
  ierr = VecSet(RHS,0.);CHKERRQ(ierr);

  /* 
    Get coordinates
  */
  ierr = DAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  /*
   Get local mat and RHS
  */
  ierr = PetscMalloc(nrow * nrow *sizeof(PetscReal), &K_local);CHKERRQ(ierr);
  ierr = PetscMalloc(nrow * sizeof(MatStencil),&row);CHKERRQ(ierr);
  ierr = DAGetLocalVector(ctx->daScal,&RHS_localVec);CHKERRQ(ierr);
  ierr = VecSet(RHS_localVec,0.);CHKERRQ(ierr);
  ierr = DAVecGetArray(ctx->daScal,RHS_localVec,&RHS_array);CHKERRQ(ierr);

  ierr = PetscMalloc(nrow * sizeof(PetscReal), &RHS_local);CHKERRQ(ierr);

  /*
    loop through all elements
  */
  for (ek = zs; ek < zs+zm; ek++) {
    for (ej = ys; ej < ys+ym; ej++) {
      for (ei = xs; ei < xs+xm; ei++) {
        hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
        hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
        hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
        ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
        /*
          Accumulate stiffness matrix
        */
        for (l = 0; l < nrow * nrow; l++) K_local[l] = 0.;
        ierr = PetscLogEventBegin(ctx->vflog.VF_MatPLocalEvent,0,0,0,0);CHKERRQ(ierr);
        ierr = VFFlow_FEM_MatPAssembly3D_local(K_local,&ctx->resprop,ek,ej,ei,&ctx->e3D);
        ierr = PetscLogEventEnd(ctx->vflog.VF_MatPLocalEvent,0,0,0,0);CHKERRQ(ierr);
        
        for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
          for (j = 0; j < ctx->e3D.nphiy; j++) {
            for (i = 0; i < ctx->e3D.nphix; i++,l++) {
              row[l].i = ei + i; row[l].j = ej + j; row[l].k = ek + k; row[l].c = 0;
            }  
          }
        } 
        
 
        ierr = MatSetValuesStencil(K,nrow,row,nrow,row,K_local,ADD_VALUES);CHKERRQ(ierr);

        /*
         Jump to next element
        */
      }  
    }   
  }  

  /*
    Source term defined at a node
  */

  if(xs+xm == nx-1) xm++;
  if(ys+ym == ny-1) ym++;
  if(zs+zm == nz-1) zm++;

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if(i == ctx->SrcLoc[0] && j == ctx->SrcLoc[1] && k == ctx->SrcLoc[2]) {
          RHS_array[k][j][i] = ctx->SrcRate;
        }
      }
    }
  }


  /*
    Global Assembly
  */
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatApplyDirichletBC(K,ctx->daScal,&ctx->bcP[0]);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
 
  

  ierr = DAVecRestoreArrayDOF(ctx->daScal,RHS_localVec,&RHS_array);CHKERRQ(ierr); 
  ierr = DALocalToGlobalBegin(ctx->daScal,RHS_localVec,RHS);CHKERRQ(ierr);
  ierr = DALocalToGlobalEnd(ctx->daScal,RHS_localVec,RHS);CHKERRQ(ierr);
  ierr = VecApplyDirichletFlowBC(RHS,fields->pressure,&ctx->bcP[0],ctx->BCpres);CHKERRQ(ierr);
  /*
   Cleanup
  */
  ierr = DAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = PetscFree3(RHS_local,K_local,row);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFFlow_FEM_MatTAssembly3D_local"
/*
  VFFlow_FEM_MatTAssembly3D_local
  
  Keita Yoshioak yoshk@chevron.com
*/
extern PetscErrorCode VFFlow_FEM_MatTAssembly3D_local(PetscReal *Mat_local,ResProp *resprop,PetscInt ek, PetscInt ej,PetscInt ei,CartFE_Element3D *e)
{
  PetscInt        g,i1,i2,j1,j2,k1,k2,l;
  PetscReal       TCond_X,TCond_Y,TCond_Z;
  PetscErrorCode  ierr;
  
  PetscFunctionBegin;
  
/*
  Thermal properties should be an element property 
*/
  TCond_X = resprop->TCond_X;
  TCond_Y = resprop->TCond_Y;
  TCond_Z = resprop->TCond_Z;
  
  
  for (l = 0,k1 = 0; k1 < e->nphiz; k1++) {
    for (j1 = 0; j1 < e->nphiy; j1++) {
	  for (i1 = 0; i1 < e->nphix; i1++) {
	    for (k2 = 0; k2 < e->nphiz; k2++) {
		  for (j2 = 0; j2 < e->nphiy; j2++) {
		    for (i2 = 0; i2 < e->nphix; i2++,l++) {
			  Mat_local[l] = 0.;
			    for (g = 0; g < e->ng; g++) {
				  Mat_local[l] += e->weight[g] * (TCond_X * e->dphi[k1][j1][i1][0][g] * e->dphi[k2][j2][i2][0][g]
			                                      + TCond_Y * e->dphi[k1][j1][i1][1][g] * e->dphi[k2][j2][i2][1][g]
												  + TCond_Z * e->dphi[k1][j1][i1][2][g] * e->dphi[k2][j2][i2][2][g]); 
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
#define __FUNCT__ "VFFlow_FEM_MatTAssembly3D"
/*
  VFFlow_FEM_MatTAssembly3D

  Keita Yoshioka yoshk@chevron.com
*/
extern PetscErrorCode VFFlow_FEM_MatTAssembly3D(Mat K,Vec RHS,VFFields *fields,VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       ei,ej,ek,i,j,k,l;
  PetscInt       nrow = ctx->e3D.nphix * ctx->e3D.nphiy * ctx->e3D.nphiz;
  Vec            RHS_localVec;
  PetscReal      ***RHS_array;
  PetscReal      *RHS_local;
  PetscReal      *K_local;
  MatStencil     *row;
  PetscReal      hx,hy,hz;
  PetscReal      ****coords_array;
  
  PetscFunctionBegin;
  ierr = PetscLogStagePush(ctx->vflog.VF_TAssemblyStage);CHKERRQ(ierr);
  ierr = DAGetInfo(ctx->daVect,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DAGetCorners(ctx->daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  if (xs+xm == nx) xm--;
  if (ys+ym == ny) ym--;
  if (zs+zm == nz) zm--;
  
  ierr = MatZeroEntries(K);CHKERRQ(ierr);
  ierr = VecSet(RHS,0.);CHKERRQ(ierr);
  
  /*
    Get coordinates
  */
  ierr = DAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  /*
    Get local mat and RHS
  */
  ierr = PetscMalloc(nrow * nrow * sizeof(PetscReal), &K_local);CHKERRQ(ierr);
  ierr = PetscMalloc(nrow * sizeof(MatStencil),&row);CHKERRQ(ierr);
  ierr = DAGetLocalVector(ctx->daScal,&RHS_localVec);CHKERRQ(ierr);
  ierr = VecSet(RHS_localVec,0.);CHKERRQ(ierr);
  ierr = DAVecGetArray(ctx->daScal,RHS_localVec,&RHS_array);CHKERRQ(ierr);
  
  ierr = PetscMalloc(nrow * sizeof(PetscReal), &RHS_local);CHKERRQ(ierr);

  /*
    loop through all elements
  */
  for (ek = zs; ek < zs+zm; ek++) {
    for (ej = ys; ej < ys+ym; ej++) {
      for (ei = xs; ei < xs+xm; ei++) {
	    hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
		hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
		hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
		ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
		
		for (l = 0; l < nrow * nrow; l++) K_local[l] = 0.;
		ierr = PetscLogEventBegin(ctx->vflog.VF_MatTLocalEvent,0,0,0,0);CHKERRQ(ierr);
		ierr = VFFlow_FEM_MatTAssembly3D_local(K_local,&ctx->resprop,ek,ej,ei,&ctx->e3D);
		ierr = PetscLogEventEnd(ctx->vflog.VF_MatTLocalEvent,0,0,0,0);CHKERRQ(ierr);
		
		// *row: MatStencil
		for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
		  for (j = 0; j < ctx->e3D.nphiy; j++) {
            for (i = 0; i < ctx->e3D.nphix; i++,l++) {
              row[l].i = ei + i; row[l].j = ej + j; row[l].k = ek + k; row[l].c = 0;
			}  
          }
		}
		
		ierr = MatSetValuesStencil(K,nrow,row,nrow,row,K_local,ADD_VALUES);CHKERRQ(ierr);
		
		// Jump to next element
	  }
    }
  }  


  /*
    Global Assembly
  */
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatApplyDirichletBC(K,ctx->daScal,&ctx->bcT[0]);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  ierr = DAVecRestoreArrayDOF(ctx->daScal,RHS_localVec,&RHS_array);CHKERRQ(ierr);
  ierr = DALocalToGlobalBegin(ctx->daScal,RHS_localVec,RHS);CHKERRQ(ierr);
  ierr = DALocalToGlobalEnd(ctx->daScal,RHS_localVec,RHS);CHKERRQ(ierr);
  ierr = VecApplyDirichletFlowBC(RHS,fields->theta,&ctx->bcT[0],ctx->BCtheta);CHKERRQ(ierr);
  
  /*
    Clean up
  */
  ierr = DAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = PetscFree3(RHS_local,K_local,row);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFFlow_FEM"
/* 
   Linear flow for testing purpose
*/

extern PetscErrorCode VFFlow_FEM(VFCtx *ctx, VFFields *fields)
{
  PetscErrorCode     ierr;
  KSPConvergedReason reasonP,reasonT;
  PetscInt           itsP,itsT;
  PetscReal          Pmin,Pmax;
  PetscReal          Tmin,Tmax;

  PetscFunctionBegin;

  if (ctx->verbose >0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Entering flow solver %s implemented in %s\n",__FUNCT__,__FILE__);CHKERRQ(ierr);
  }

  ierr = VFFlow_FEM_MatPAssembly3D(ctx->KP,ctx->RHSP,fields,ctx);CHKERRQ(ierr);
  ierr = VFFlow_FEM_MatTAssembly3D(ctx->KT,ctx->RHST,fields,ctx);CHKERRQ(ierr);
  if (ctx->verbose > 1) {
    ierr = MatView(ctx->KP,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	ierr = MatView(ctx->KT,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = VecView(ctx->RHSP,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	ierr = VecView(ctx->RHST,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  } 

  ierr = PetscLogStagePush(ctx->vflog.VF_PSolverStage);CHKERRQ(ierr);
  ierr = KSPSolve(ctx->kspP,ctx->RHSP,fields->pressure);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);
  
  ierr = PetscLogStagePush(ctx->vflog.VF_TSolverStage);CHKERRQ(ierr);
  ierr = KSPSolve(ctx->kspT,ctx->RHST,fields->theta);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  if (ctx->verbose >1) {
    ierr = VecView(fields->pressure,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	ierr = VecView(fields->theta,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  if (ctx->verbose > 0) {
    ierr = VecMin(fields->pressure,PETSC_NULL,&Pmin);CHKERRQ(ierr);
    ierr = VecMax(fields->pressure,PETSC_NULL,&Pmax);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Pmin = %g, Pmax = %g\n",Pmin,Pmax);CHKERRQ(ierr);
	ierr = VecMin(fields->theta,PETSC_NULL,&Tmin);CHKERRQ(ierr);
	ierr = VecMax(fields->theta,PETSC_NULL,&Tmax);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Tmin = %g, Tmax = %g\n",Tmin,Tmax);CHKERRQ(ierr);
  }

  ierr = KSPGetConvergedReason(ctx->kspP,&reasonP);CHKERRQ(ierr);
  ierr = KSPGetConvergedReason(ctx->kspT,&reasonT);CHKERRQ(ierr);
  if (reasonP < 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"[ERROR] kspP diverged with reason %d\n", (int)reasonP);CHKERRQ(ierr);
  } else {
    ierr = KSPGetIterationNumber(ctx->kspP,&itsP);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"      kspP converged in %d iterations %d. \n",(int)itsP,(int)reasonP);CHKERRQ(ierr);
  } 
  if(reasonT < 0){
	ierr = PetscPrintf(PETSC_COMM_WORLD,"[ERROR] kspT diverged with reason %d\n", (int)reasonT);CHKERRQ(ierr);
  } else {
	ierr = KSPGetIterationNumber(ctx->kspT,&itsT);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"      kspT converged in %d iterations %d. \n",(int)itsT,(int)reasonT);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
/*
#undef __FUNCT__
#define __FUNCT__ "VFFormFunction_Flow_local"

extern PetscErrorCode VFFormFunction_Flow_local(PetscReal *func_local,ResProp *resprop, PetscInt ek, PetscInt ej, PetscInt ei,CartFE_Element3D *e)
{
  PetscErrorCode ierr;
  PetscInt       g,i1,i2,j1,j2,k1,k2,i,j,k,l;
  
  PetscFunctionBegin;
  
  for(g=0; g < e->ng; g+=) {
    for(l=0,k=0; k < e->nphiz; k++) {
	  for(j=0; j < e->nphiy; j++) {
	    for(i=0; i < e->nphix; i++) {
		  pressur[
		  
  
  PetscFunctionReturn(0);
}
  
  */

#undef __FUNCT__
#define __FUNCT__ "VFFormIBCondition_Flow"
/*
  Set initial and Dirichlet boundary condition
*/
extern PetscErrorCode VFFormIBCondition_Flow(VFCtx *ctx, VFFields *fields)
{
  PetscInt        xs,xm,nx;
  PetscInt        ys,ym,ny;
  PetscInt        zs,zm,nz;
  PetscInt        i,j,k;
  PetscErrorCode  ierr;
  PetscReal       ***pressure_array;
  
  PetscFunctionBegin;
  
  ierr = DAGetInfo(ctx->daVect,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DAGetCorners(ctx->daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);  
  ierr = DAVecGetArray(ctx->daScal,fields->pressure,&pressure_array);
  
  for (k = zs; k < zs+zm; k++) {
    for (j = ys; j < ys+ym; j++) {
      for (i = xs; i < xs+xm; i++) {
	    /*
		  x == 0
		*/
	    if (i == 0 && ctx->bcP[0].face[X0] == VALUE) {
		  pressure_array[k][j][i] = ctx->BCpres[X0];
		}
		/*
		  x == nx-1
		*/
		else if (i == nx-1 && ctx->bcP[0].face[X1] == VALUE) {
		  pressure_array[k][j][i] = ctx->BCpres[X1];
		}
		/*
		  y == 0
		*/
		else if (j == 0 && ctx->bcP[0].face[Y0] == VALUE) {
		  pressure_array[k][j][i] = ctx->BCpres[Y0];
		}
		/*
		  y == ny-1
		*/
		else if (j == ny-1 && ctx->bcP[0].face[Y1] == VALUE) {
		  pressure_array[k][j][i] = ctx->BCpres[Y1];
		}
		/*
		  z == 0
		*/
		else if (k == 0 && ctx->bcP[0].face[Z0] == VALUE) {
		  pressure_array[k][j][i] = ctx->BCpres[Z0];
		}
		/*
		  z == nz-1
		*/
		else if (k == nz-1 && ctx->bcP[0].face[Z1] == VALUE) {
		  pressure_array[k][j][i] = ctx->BCpres[Z1];
		}
		else {
		  pressure_array[k][j][i] = ctx->resprop.Pinit;
		}
	  }
	}
  }  
  PetscFunctionReturn(0);
}
  

#undef __FUNCT__
#define __FUNCT__ "VFFormFunction_Flow"

extern PetscErrorCode VFFormFunction_Flow(SNES snes, Vec pressure_Vec, Vec F, void *voidctx)
{
  VFCtx           *ctx = (VFCtx*)voidctx;
  PetscErrorCode  ierr;
  PetscInt        xs,xm,nx;
  PetscInt        ys,ym,ny;
  PetscInt        zs,zm,nz;
  PetscInt        ei,ej,ek;
  PetscInt        i1,j1,k1,i2,j2,k2,g;
  PetscInt        nrow = ctx->e3D.nphix * ctx->e3D.nphiy * ctx->e3D.nphiz;
  Vec             pressure_localVec, func_localVec;
  PetscReal       ***pressure_array,***func_array;
  PetscReal       ****coords_array;
  PetscReal       hx,hy,hz;
  PetscReal       fdens,relk,visc;
  PetscReal       kxx,kyy,kzz,kxy,kxz,kyz;
  PetscReal       DCoef_P;
  
  PetscFunctionBegin;
  
  fdens = ctx->resprop.fdens;
  relk = ctx->resprop.relk;
  visc = ctx->resprop.visc;
  kxx = ctx->resprop.perm;
  kyy = ctx->resprop.perm;
  kzz = ctx->resprop.perm;
  DCoef_P = fdens*relk/visc;


  
  ierr = DAGetInfo(ctx->daVect,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DAGetCorners(ctx->daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  if (xs+xm == nx) xm--;
  if (ys+ym == ny) ym--;
  if (zs+zm == nz) zm--;
  
  /*
    get coordinates
  */
  ierr = DAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  
  /*
    get pressure and function array
  */  
  ierr = DAGetLocalVector(ctx->daScal,&pressure_localVec);CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(ctx->daScal,pressure_Vec,INSERT_VALUES,pressure_localVec);CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(ctx->daScal,pressure_Vec,INSERT_VALUES,pressure_localVec);CHKERRQ(ierr);
  ierr = DAVecGetArray(ctx->daScal,pressure_localVec,&pressure_array);CHKERRQ(ierr);
//  ierr = DAVecGetArray(ctx->daScal,pressure_Vec,&pressure_array);CHKERRQ(ierr);
  
  ierr = DAGetLocalVector(ctx->daScal,&func_localVec);CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(ctx->daScal,F,INSERT_VALUES,func_localVec);CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(ctx->daScal,F,INSERT_VALUES,func_localVec);CHKERRQ(ierr);
  ierr = DAVecGetArray(ctx->daScal,func_localVec,&func_array);CHKERRQ(ierr);
//  ierr = DAVecGetArray(ctx->daScal,F,&func_array);CHKERRQ(ierr);
 
  /*
    loop through all the elements
  */
  for (ek = zs; ek < zs+zm; ek++) {
    for (ej = ys; ej < ys+ym; ej++) {
	  for (ei = xs; ei < xs+xm; ei++) {
        hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
        hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
        hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
        ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
//		func_array[ek][ej][ei] = pressure_array[ek][ej][ei]*pressure_array[ek][ej][ei];
//		func_array[ek][ej][ei] = 1.2345;
//		pressure_array[ek][ej][ei] = 0.99;
		for(g = 0; g < ctx->e3D.ng; g++) {
		  for (k1 = 0; k1 < ctx->e3D.nphiz; k1++) {
		    for (j1 = 0; j1 < ctx->e3D.nphiy; j1++) {
			  for (i1 = 0; i1 < ctx->e3D.nphix; i1++) {
			    for(k2 = 0; k2 < ctx->e3D.nphiz; k2++) {
				  for(j2 = 0; j2 < ctx->e3D.nphiy; j2++) {
				    for(i2 = 0; i2 < ctx->e3D.nphix; i2++) {					 
					  /*
					    Dirichlet Boundary Condition
					  */
					  if ((ei+i1 == 0 && ctx->bcP[0].face[X0] == VALUE) || (ei+i1 == nx-1 && ctx->bcP[0].face[X1] == VALUE) \
					  ||  (ej+j1 == 0 && ctx->bcP[0].face[Y0] == VALUE) || (ej+j1 == ny-1 && ctx->bcP[0].face[Y1] == VALUE) \
					  ||  (ek+k1 == 0 && ctx->bcP[0].face[Z0] == VALUE) || (ek+k1 == nz-1 && ctx->bcP[0].face[Z1] == VALUE)) {
					    func_array[ek+k1][ej+j1][ei+i1] = pressure_array[ek+k1][ej+j1][ei+i1];
					  }
                      else {
					    func_array[ek+k1][ej+j1][ei+i1] += ctx->e3D.weight[g] * pressure_array[ek+k2][ej+j2][ei+i2] *DCoef_P* pressure_array[ek+k2][ej+j2][ei+i2] *
					    ( kxx * ctx->e3D.dphi[k1][j1][i1][0][g] * ctx->e3D.dphi[k2][j2][i2][0][g]
					    + kyy * ctx->e3D.dphi[k1][j1][i1][1][g] * ctx->e3D.dphi[k2][j2][i2][1][g]
					    + kzz * ctx->e3D.dphi[k1][j1][i1][2][g] * ctx->e3D.dphi[k2][j2][i2][2][g] );
					  } 
					}
				  }
				}
			  }
			}
		  }
		}
		PetscPrintf(PETSC_COMM_WORLD,"func_array[%d][%d][%d] = %e \n",ek,ej,ei,func_array[ek][ej][ei]);
	  }
	}
  }
  
  /*
    (Dirichlet) Boundary Condition
  */
	
  /*	
    Clean up
  */
  ierr = DAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);

  ierr = DAVecRestoreArray(ctx->daScal,pressure_localVec,&pressure_array);CHKERRQ(ierr);
  ierr = DARestoreLocalVector(ctx->daScal,&pressure_localVec);CHKERRQ(ierr);
 
  ierr = DAVecRestoreArray(ctx->daScal,func_localVec,&func_array);CHKERRQ(ierr);
//  ierr = DALocalToGlobal(ctx->daScal,func_localVec,F);CHKERRQ(ierr);
  ierr = DALocalToGlobalBegin(ctx->daScal,func_localVec,F);CHKERRQ(ierr);
  ierr = DALocalToGlobalEnd(ctx->daScal,func_localVec,F);CHKERRQ(ierr);
//  ierr = DARestoreLocalVector(ctx->daScal,&func_localVec);CHKERRQ(ierr);
	
//	ierr = DAVecRestoreArray(ctx->daScal,pressure_Vec,&pressure_array);CHKERRQ(ierr);
	
	
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VFFormJacobian_Flow"

extern PetscErrorCode VFFormJacobian_Flow(SNES snes, Vec pressure_Vec, Mat *J, Mat *B, MatStructure *flag, void *voidctx)
{
  VFCtx           *ctx = (VFCtx*)voidctx;
  PetscErrorCode ierr;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       ei,ej,ek,i,j,k,l;
  PetscInt       nrow = ctx->e3D.nphix * ctx->e3D.nphiy * ctx->e3D.nphiz;
  PetscReal      *J_local;
  MatStencil     *row;
  PetscReal      hx,hy,hz;
  PetscReal      ****coords_array;
 
  PetscFunctionBegin;
  ierr = PetscLogStagePush(ctx->vflog.VF_PAssemblyStage);CHKERRQ(ierr);
  ierr = DAGetInfo(ctx->daVect,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DAGetCorners(ctx->daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  if (xs+xm == nx) xm--;
  if (ys+ym == ny) ym--;
  if (zs+zm == nz) zm--;

  ierr = MatZeroEntries(J);CHKERRQ(ierr);

  /* 
    Get coordinates
  */
  ierr = DAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  /*
   Get local Jacobian
  */
  ierr = PetscMalloc(nrow * nrow *sizeof(PetscReal), &J_local);CHKERRQ(ierr);
  ierr = PetscMalloc(nrow * sizeof(MatStencil),&row);CHKERRQ(ierr);

  /*
    loop through all elements
  */
  for (ek = zs; ek < zs+zm; ek++) {
    for (ej = ys; ej < ys+ym; ej++) {
      for (ei = xs; ei < xs+xm; ei++) {
        hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
        hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
        hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
        ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
        /*
          Accumulate Jacobian matrix
        */
        for (l = 0; l < nrow * nrow; l++) J_local[l] = 0.;
        ierr = PetscLogEventBegin(ctx->vflog.VF_MatPLocalEvent,0,0,0,0);CHKERRQ(ierr);
        ierr = VFFlow_FEM_MatPAssembly3D_local(J_local,&ctx->resprop,ek,ej,ei,&ctx->e3D);
        ierr = PetscLogEventEnd(ctx->vflog.VF_MatPLocalEvent,0,0,0,0);CHKERRQ(ierr);
        
        for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
          for (j = 0; j < ctx->e3D.nphiy; j++) {
            for (i = 0; i < ctx->e3D.nphix; i++,l++) {
              row[l].i = ei + i; row[l].j = ej + j; row[l].k = ek + k; row[l].c = 0;
            }  
          }
        } 
        
 
        ierr = MatSetValuesStencil(J,nrow,row,nrow,row,J_local,ADD_VALUES);CHKERRQ(ierr);

        /*
         Jump to next element
        */
      }  
    }   
  }  

  /*
    Global Assembly
  */
  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatApplyDirichletBC(J,ctx->daScal,&ctx->bcP[0]);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /*
   Cleanup
  */
  ierr = DAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = PetscFree2(J_local,row);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);
  
  *flag = SAME_NONZERO_PATTERN;
  
  PetscFunctionReturn(0);
}
  

#undef __FUNCT__
#define __FUNCT__ "VFFlow_SNES_FEM"

extern PetscErrorCode VFFlow_SNES_FEM(VFCtx *ctx, VFFields *fields)
{
  PetscErrorCode       ierr;
  ISColoring           iscoloring;
  PetscInt             itsP;
  SNESConvergedReason  reasonP;
  PetscReal            Pmin,Pmax;
  Vec                  r;
  Mat                  J;
  MatFDColoring        matfdcoloring;
  SNES                 snes;
  
  PetscFunctionBegin;
  
  //ierr = PetscPrintf(PETSC_COMM_WORLD,"No function is defined yet \n");CHKERRQ(ierr);
  
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = VecDuplicate(fields->pressure,&r);CHKERRQ(ierr);
  
  ierr = SNESSetFunction(snes,r,VFFormFunction_Flow,ctx);CHKERRQ(ierr);
  
  // create if-then based on Jacobian computation choise (numerical or analytical) later
/*  ierr = DAGetColoring(ctx->daScal,IS_COLORING_GLOBAL,MATAIJ,&iscoloring);CHKERRQ(ierr);
  ierr = DAGetMatrix(ctx->daScal,MATAIJ,&J);CHKERRQ(ierr);
  ierr = MatFDColoringCreate(J,iscoloring,&matfdcoloring);CHKERRQ(ierr);
  ierr = ISColoringDestroy(iscoloring);CHKERRQ(ierr);
  ierr = MatFDColoringSetFunction(matfdcoloring,(PetscErrorCode (*)(void))VFFormFunction_Flow,&ctx);CHKERRQ(ierr);
  ierr = MatFDColoringSetFromOptions(matfdcoloring);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,J,J,SNESDefaultComputeJacobianColor,matfdcoloring);CHKERRQ(ierr);
*/	  
  ierr = SNESSetJacobian(snes,J,J,VFFormJacobian_Flow,&ctx);CHKERRQ(ierr);  

	  // move to VFCommon later?
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  
  ierr = VFFormIBCondition_Flow(ctx,fields); CHKERRQ(ierr);
  
  ierr = PetscLogStagePush(ctx->vflog.VF_PSolverStage);CHKERRQ(ierr);
  ierr = SNESSolve(snes,PETSC_NULL,fields->pressure);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);
  
  // TEST - explicitly calcualte the residual
/*  PetscReal fnorm;
  ierr = VFFormFunction_Flow(snes,fields->pressure,r,ctx);CHKERRQ(ierr);
  ierr = VecNorm(r,NORM_2,&fnorm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"fnorm %g \n",fnorm);CHKERRQ(ierr);
  ierr = VecView(r,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = VecView(fields->pressure,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
*/  
  if (ctx->verbose >1) {
    ierr = VecView(fields->pressure,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//	ierr = VecView(fields->theta,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  if (ctx->verbose > 0) {
    ierr = VecMin(fields->pressure,PETSC_NULL,&Pmin);CHKERRQ(ierr);
    ierr = VecMax(fields->pressure,PETSC_NULL,&Pmax);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Pmin = %g, Pmax = %g\n",Pmin,Pmax);CHKERRQ(ierr);
/*	ierr = VecMin(fields->theta,PETSC_NULL,&Tmin);CHKERRQ(ierr);
	ierr = VecMax(fields->theta,PETSC_NULL,&Tmax);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Tmin = %g, Tmax = %g\n",Tmin,Tmax);CHKERRQ(ierr); */
  } 
  ierr = SNESGetConvergedReason(snes, &reasonP);CHKERRQ(ierr);
  
  if (reasonP < 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"[ERROR] snesF diverged with reason %d\n", (int)reasonP);CHKERRQ(ierr);
  } else {
    ierr = SNESGetIterationNumber(snes,&itsP);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"      snesF converged in %d iterations with reason %d. \n",(int)itsP,(int)reasonP);CHKERRQ(ierr);
  } 
  
  /*
    clean up
  */
  ierr = MatDestroy(J);CHKERRQ(ierr);
//  ierr = MatFDColoringDestroy(matfdcoloring);CHKERRQ(ierr);
  ierr = VecDestroy(r);CHKERRQ(ierr);
  ierr = SNESDestroy(snes);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}











































