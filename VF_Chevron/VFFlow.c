#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"


#undef __FUNCT__
#define __FUNCT__ "BCPInit"
/* 
   BCPInit
*/
extern PetscErrorCode BCPInit(BC *BC,VFCtx *ctx)
{
  PetscErrorCode ierr;
 
  PetscFunctionBegin;
  ierr = BCInit(BC,1);CHKERRQ(ierr);

/*
  if given value is positive, this applies 
*/
  if(ctx->BCpres[0] > -1.e-8) BC[0].face[X0] = VALUE;
  if(ctx->BCpres[1] > -1.e-8) BC[0].face[X1] = VALUE;
  if(ctx->BCpres[2] > -1.e-8) BC[0].face[Y0] = VALUE;
  if(ctx->BCpres[3] > -1.e-8) BC[0].face[Y1] = VALUE;
  if(ctx->BCpres[4] > -1.e-8) BC[0].face[Z0] = VALUE;
  if(ctx->BCpres[5] > -1.e-8) BC[0].face[Z1] = VALUE;

  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "VF_MatP3D_local"
/* 
   VF_MatP3D_local
*/

extern PetscErrorCode VF_MatP3D_local(PetscReal *Mat_local,ResProp *resprop,PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e)
{
  PetscInt       g,i1,i2,j1,j2,k1,k2,l;
  PetscReal      fdens,relk,perm,visc;
  PetscReal      DCoef_P;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  fdens = resprop->fdens;
  relk = resprop->relk;
  perm = resprop->perm;
  visc = resprop->visc;
  DCoef_P = fdens*relk*perm/visc;
  
  for (l = 0,k1 = 0; k1 < e->nphiz; k1++) {
    for (j1 = 0; j1 < e->nphiy; j1++) {
      for (i1 = 0; i1 < e->nphix; i1++) {
        for (k2 = 0; k2 < e->nphiz; k2++) {
          for (j2 = 0; j2 < e->nphiy; j2++) {
            for (i2 = 0; i2 < e->nphix; i2++,l++) {
              Mat_local[l] = 0.;
              for (g = 0; g < e->ng; g++) {
                Mat_local[l] += e->weight[g] * DCoef_P * ( e->dphi[k1][j1][i1][0][g] * e->dphi[k2][j2][i2][0][g]
                                                         + e->dphi[k1][j1][i1][1][g] * e->dphi[k2][j2][i2][1][g]
                                                         + e->dphi[k1][j1][i1][2][g] * e->dphi[k2][j2][i2][2][g]);
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
#define __FUNCT__ "VF_PAssembly3D"
/* 
   VF_PAssembly3D
*/
extern PetscErrorCode VF_PAssembly3D(Mat K,Vec RHS,VFFields *fields,VFCtx *ctx)
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
  ierr = DAGetLocalVector(ctx->daScal,&RHS_localVec);CHKERRQ(ierr)
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
        ierr = VF_MatP3D_local(K_local,&ctx->resprop,ek,ej,ei,&ctx->e3D);
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
          RHS - add source terms
        */

        /*
         Jump to next element
        */
      }  
    }   
  }  

  /*
    Global Assembly and Source/Sink terms (Boundary Conditions ?)
  */
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatApplyDirichletBC(K,ctx->daScal,&ctx->bcP[0]);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = DAVecRestoreArrayDOF(ctx->daScal,RHS_localVec,&RHS_array);CHKERRQ(ierr); 
  ierr = DALocalToGlobalBegin(ctx->daScal,RHS_localVec,RHS);CHKERRQ(ierr);
  ierr = DALocalToGlobalEnd(ctx->daScal,RHS_localVec,RHS);CHKERRQ(ierr);
  ierr = VecApplyDirichletFlowBC(RHS,fields->pressure,&ctx->bcP[0],&ctx->BCpres);CHKERRQ(ierr);
  /*
   Cleanup
  */
  ierr = DAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = PetscFree3(RHS_local,K_local,row);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VF_FakeFlow"
/* 
   Fake flow solver for VF_Chevron.c test
*/

extern PetscErrorCode VF_FakeFlow(VFCtx *ctx, VFFields *fields)
{
  PetscErrorCode ierr;
  PetscReal      time;
  PetscReal      pres,temp;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       i,j,k;
  Vec            theta_localVec;
  Vec            pressure_localVec;
  PetscReal      ***theta_array;
  PetscReal      ***pressure_array;
  PetscReal      ****coords_array;
  PetscReal      x,y,z,r;
  PetscReal      Tinit,Pinit,P1,P2,T1,T2;

  PetscFunctionBegin;

  ierr = DAGetInfo(ctx->daVect,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DAGetCorners(ctx->daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  /*
     Get coordinate, temperature array, and pressure array
  */
  ierr = DAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DAVecGetArray(ctx->daScal,fields->theta,&theta_array);CHKERRQ(ierr);
  ierr = DAVecGetArray(ctx->daScal,fields->pressure,&pressure_array);CHKERRQ(ierr);

  time = ctx->timevalue;
  Tinit = ctx->resprop.Tinit;
  Pinit = ctx->resprop.Pinit;
  P1 = 1.;
  P2 = 1.e6;
  T1 = 2.;
  T2 = 1.e6;

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        x = coords_array[k][j][i][0]+0.5;
        y = coords_array[k][j][i][1]+0.5;
        z = coords_array[k][j][i][2]+0.5;
        r = sqrt(x*x+y*y);
        theta_array[k][j][i] = Tinit - T1*log(T2*time/(r*r));
        pressure_array[k][j][i] = Pinit + P1*log(P2*time/(r*r));
      }
    }
  } 
  
  /*
     Cleanup
  */
  ierr = DAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(ctx->daScal,fields->theta,&theta_array);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(ctx->daScal,fields->pressure,&pressure_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_LinFlow"
/* 
   Linear flow for testing purpose
*/

extern PetscErrorCode VF_LinFlow(VFCtx *ctx, VFFields *fields)
{
  PetscErrorCode     ierr;
  KSPConvergedReason reason;
  PetscInt           its;
  PetscReal          Pmin,Pmax;

  PetscFunctionBegin;
  ierr = VF_PAssembly3D(ctx->KP,ctx->RHSP,fields,ctx);CHKERRQ(ierr);
  if (ctx->verbose > 1) {
    ierr = MatView(ctx->KP,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = VecView(ctx->RHSP,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  } 

  ierr = PetscLogStagePush(ctx->vflog.VF_PSolverStage);CHKERRQ(ierr);
  ierr = KSPSolve(ctx->kspP,ctx->RHSP,fields->pressure);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  if (ctx->verbose >1) {
    ierr = VecView(fields->pressure,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  if (ctx->verbose > 0) {
    ierr = VecMin(fields->pressure,PETSC_NULL,&Pmin);CHKERRQ(ierr);
    ierr = VecMax(fields->pressure,PETSC_NULL,&Pmax);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Pmin = %g, Pmax = %g\n",Pmin,Pmax);CHKERRQ(ierr);
  }

  ierr = KSPGetConvergedReason(ctx->kspP,&reason);CHKERRQ(ierr);
  if (reason < 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"[ERROR] kspP diverged with reason %d\n", (int)reason);CHKERRQ(ierr);
  } else {
    ierr = KSPGetIterationNumber(ctx->kspP,&its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"      kspP converged in %d iterations %d. \n",(int)its,(int)reason);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VFFlowTimeStep"
extern PetscErrorCode VFFlowTimeStep(VFCtx *ctx,VFFields *fields)
{
  PetscErrorCode     ierr;

  PetscFunctionBegin;
  
//  ierr = VF_FakeFlow(ctx,fields); 
  ierr = VF_LinFlow(ctx,fields);
  PetscFunctionReturn(0);
}



  
