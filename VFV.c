#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFU.h"
#include "VFV.h"


#undef __FUNCT__
#define __FUNCT__ "BCVInit"
/*
  BCVInit

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode BCVInit(BC *BC,VFPreset preset)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = BCInit(BC,1);CHKERRQ(ierr);
  switch (preset) {
      /*
        Preventing fracture through the lower, upper, and non-symmetry planes, when BC are imposed
      */
    case SYMXY:
    case SYMX:
    case SYMY:
    case NOSYM:
       BC[0].face[Z1] = ONE; 
      break;
    case TEST_MANUAL:
      ierr = BCGet(BC,"V",1);
      break;
    default:
      SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_USER,"ERROR: [%s] unknown preset %i.\n",__FUNCT__,preset);
      break;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BCVUpdate"
/*
  BCVUpdate

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode BCVUpdate(BC *BC,VFPreset preset)
{
  PetscFunctionBegin;
  switch (preset) {
      /*
        Preventing fracture through the lower, upper, and non-symmetry planes, when BC are imposed
      */
    case SYMXY:
       BC[0].face[X1] = ONE; 
       BC[0].face[Y1] = ONE; 
       BC[0].face[Z0] = ONE; 
       BC[0].face[Z1] = ONE; 
       break;
     case SYMX:
       BC[0].face[X1] = ONE; 
       BC[0].face[Y0] = ONE; 
       BC[0].face[Y1] = ONE; 
       BC[0].face[Z0] = ONE; 
       BC[0].face[Z1] = ONE; 
       break;
     case SYMY:
       BC[0].face[X0] = ONE; 
       BC[0].face[Y0] = ONE; 
       BC[0].face[Z0] = ONE; 
       BC[0].face[Z1] = ONE; 
       break;
     case NOSYM:
       BC[0].face[X0] = ONE; 
       BC[0].face[X1] = ONE; 
       BC[0].face[Y0] = ONE; 
       BC[0].face[Y1] = ONE; 
       BC[0].face[Z0] = ONE; 
       BC[0].face[Z1] = ONE; 
       break;
    case TEST_MANUAL:
      break;
    default:
      SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_USER,"ERROR: [%s] unknown preset %i.\n",__FUNCT__,preset);
      break;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_MatVAT2Surface3D_local"
/*
  VF_MatVAT2Surface3D_local

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_MatVAT2Surface3D_local(PetscReal *Mat_local,MatProp *matprop,VFProp *vfprop,CartFE_Element3D *e,PetscReal Gc)
{
  PetscInt       g,i1,i2,j1,j2,k1,k2,l;
  PetscReal      coef = Gc / vfprop->atCv * .5;

  PetscFunctionBegin;
  for (l = 0,k1 = 0; k1 < e->nphiz; k1++) {
    for (j1 = 0; j1 < e->nphiy; j1++) {
      for (i1 = 0; i1 < e->nphix; i1++) {
        for (k2 = 0; k2 < e->nphiz; k2++) {
          for (j2 = 0; j2 < e->nphiy; j2++) {
            for (i2 = 0; i2 < e->nphix; i2++,l++) {
              for (g = 0; g < e->ng; g++) {
                Mat_local[l] += coef * e->weight[g] * (e->phi[k1][j1][i1][g] *     e->phi[k2][j2][i2][g] / vfprop->epsilon +
                                                      (e->dphi[k1][j1][i1][0][g] * e->dphi[k2][j2][i2][0][g] 
                                                     + e->dphi[k1][j1][i1][1][g] * e->dphi[k2][j2][i2][1][g]
                                                     + e->dphi[k1][j1][i1][2][g] * e->dphi[k2][j2][i2][2][g]) * vfprop->epsilon);
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
#define __FUNCT__ "VF_MatVCoupling3D_local"
/*
  VF_MatVCoupling3D_local

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_MatVCoupling3D_local(PetscReal *Mat_local,PetscReal ****U_array,
                                              PetscReal ***theta_array,PetscReal ***thetaRef_array,
                                              PetscReal ***pressure_array,PetscReal ***pressureRef_array,
                                              MatProp *matprop,VFProp *vfprop,
                                              PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e)
{
  PetscErrorCode ierr;
  PetscInt       g,i1,i2,j1,j2,k1,k2,l;
  PetscReal      *ElasticEnergyDensity_local;
  
  PetscFunctionBegin;
  ierr = PetscMalloc(e->ng * sizeof(PetscReal),&ElasticEnergyDensity_local);CHKERRQ(ierr);
  for (g = 0; g < e->ng; g++) ElasticEnergyDensity_local[g] = 0;
  ierr = ElasticEnergyDensity3D_local(ElasticEnergyDensity_local,U_array,
                                      theta_array,thetaRef_array,
                                      pressure_array,pressureRef_array,
                                      matprop,ek,ej,ei,e);CHKERRQ(ierr);

  for (l = 0,k1 = 0; k1 < e->nphiz; k1++) {
    for (j1 = 0; j1 < e->nphiy; j1++) {
      for (i1 = 0; i1 < e->nphix; i1++) {
        for (k2 = 0; k2 < e->nphiz; k2++) {
          for (j2 = 0; j2 < e->nphiy; j2++) {
            for (i2 = 0; i2 < e->nphix; i2++,l++) {
              for (g = 0; g < e->ng; g++) {
                Mat_local[l] += e->weight[g] * e->phi[k1][j1][i1][g] * e->phi[k2][j2][i2][g] * ElasticEnergyDensity_local[g] * 2.;
              }
            }
          }
        }
      }
    }
  }
  ierr = PetscFree(ElasticEnergyDensity_local);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VF_MatVCouplingShearOnly3D_local"
/*
  VF_MatVCouplingShearOnly3D_local

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_MatVCouplingShearOnly3D_local(PetscReal *Mat_local,PetscReal ****U_array,
                                                       PetscReal ***theta_array,PetscReal ***thetaRef_array,
                                                       PetscReal ***pressure_array,PetscReal ***pressureRef_array,
                                                       MatProp *matprop,VFProp *vfprop,
                                                       PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e)
{
  PetscErrorCode ierr;
  PetscInt       g,i1,i2,j1,j2,k1,k2,l;
  PetscReal      *ElasticEnergyDensityS_local,*ElasticEnergyDensityD_local;
  
  PetscFunctionBegin;
  ierr = PetscMalloc2(e->ng,PetscReal,&ElasticEnergyDensityS_local,e->ng,PetscReal,&ElasticEnergyDensityD_local);CHKERRQ(ierr);
  ierr = ElasticEnergyDensitySphericalDeviatoric3D_local(ElasticEnergyDensityS_local,ElasticEnergyDensityD_local,
                                                         U_array,theta_array,thetaRef_array,
                                                         pressure_array,pressureRef_array,
                                                         matprop,ek,ej,ei,e);CHKERRQ(ierr);
  for (l = 0,k1 = 0; k1 < e->nphiz; k1++) {
    for (j1 = 0; j1 < e->nphiy; j1++) {
      for (i1 = 0; i1 < e->nphix; i1++) {
        for (k2 = 0; k2 < e->nphiz; k2++) {
          for (j2 = 0; j2 < e->nphiy; j2++) {
            for (i2 = 0; i2 < e->nphix; i2++,l++) {
              for (g = 0; g < e->ng; g++) {
                Mat_local[l] += e->weight[g] * e->phi[k1][j1][i1][g] * e->phi[k2][j2][i2][g] * ElasticEnergyDensityD_local[g] * 2.;
              }
            }
          }
        }
      }
    }
  }
  ierr = PetscFree2(ElasticEnergyDensityS_local,ElasticEnergyDensityD_local);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VF_RHSVPressure3D_local"
/*
  VF_RHSVPressure3D_local

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_RHSVPressure3D_local(PetscReal *RHS_local,PetscReal ****u_array,
                                              PetscReal ***pressure_array,PetscReal ***pressureRef_array,
                                              MatProp *matprop,VFProp *vfprop,
                                              PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e)
{
  PetscErrorCode ierr;
  PetscInt       l,i,j,k,g,c;
  PetscReal      *pressure_elem,*u_elem[3];

  PetscFunctionBegin;
  ierr = PetscMalloc4(e->ng,PetscReal,&pressure_elem,
                      e->ng,PetscReal,&u_elem[0],
                      e->ng,PetscReal,&u_elem[1],
                      e->ng,PetscReal,&u_elem[2]);CHKERRQ(ierr);

  /*
    Compute the projection of the fields in the local base functions basis
  */
  for (g = 0; g < e->ng; g++) {
    pressure_elem[g] = 0;
    for (c=0; c<3; c++) u_elem[c][g] = 0.;
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (g = 0; g < e->ng; g++) {
          pressure_elem[g] += e->phi[k][j][i][g] 
                              * (pressure_array[ek+k][ej+j][ei+i] - pressureRef_array[ek+k][ej+j][ei+i]);
          for (c=0; c<3; c++) u_elem[c][g] += e->phi[k][j][i][g] * u_array[ek+k][ej+j][ei+i][c];
        }
      }
    }
  }  
  ierr = PetscLogFlops(9 * e->ng * e->nphix * e->nphiy * e->nphiz);CHKERRQ(ierr);

  /*
    Accumulate the contribution of the current element to the local
    version of the RHS
  */
  for (l=0,k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++,l++) {
        for (c = 0; c < 3; c++) {
          for (g = 0; g < e->ng; g++) {
            RHS_local[l] += e->weight[g] * pressure_elem[g] 
                          * u_elem[c][g]            
                          * e->dphi[k][j][i][c][g];  
          }
        }
      }
    }
  }
  ierr = PetscLogFlops(12 * e->ng * e->nphix * e->nphiy * e->nphiz);CHKERRQ(ierr);
  /*
    Clean up
  */
  ierr = PetscFree4(pressure_elem,u_elem[0],u_elem[1],u_elem[2]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_RHSVAT2Surface3D_local"
/*
  VF_RHSVAT2Surface3D_local

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_RHSVAT2Surface3D_local(PetscReal *RHS_local,MatProp *matprop,VFProp *vfprop,CartFE_Element3D *e, PetscReal Gc)
{
  PetscInt       g,i,j,k,l;
  PetscReal      coef = Gc / vfprop->atCv / vfprop->epsilon *.5;

  PetscFunctionBegin;
  for (l=0,k=0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++,l++) {
        for (g = 0; g < e->ng; g++) {
          RHS_local[l] += e->weight[g] * e->phi[k][j][i][g] * coef;
        }
      }
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_VAssembly3D"
/*
  VF_VAssembly3D

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_VAssembly3D(Mat K,Vec RHS,VFFields *fields,VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       ei,ej,ek,i,j,k,l;
  PetscInt       nrow = ctx->e3D.nphix * ctx->e3D.nphiy * ctx->e3D.nphiz;
  Vec            RHS_localVec,U_localVec;
  Vec            theta_localVec,thetaRef_localVec;
  Vec            pressure_localVec,pressureRef_localVec;
  Vec            Gc_localVec;
  PetscReal      ***RHS_array,****U_array;
  PetscReal      ***theta_array,***thetaRef_array;
  PetscReal      ***pressure_array,***pressureRef_array;
  PetscReal      ***Gc_array;
  PetscReal      *RHS_local;
  PetscReal      *K_local;
  MatStencil     *row;
  PetscReal      hx,hy,hz;
  PetscReal      ****coords_array;
  
  PetscFunctionBegin;
  /* 
    Get global number of vertices along each coordinate axis on the ENTIRE mesh
  */
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  /*
    Get informations on LOCAL slice (i.e. subdomain)
    xs, ys, ym = 1st index in x,y,z direction
    xm, ym, zm = number of vertices in the x, y, z directions
  */
  ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  ierr = MatZeroEntries(K);CHKERRQ(ierr);
  ierr = VecSet(RHS,0.);CHKERRQ(ierr);
  /*
    Get coordinates
  */
  ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  
  /*
    get U_array
  */
  ierr = DMGetLocalVector(ctx->daVect,&U_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daVect,fields->U,INSERT_VALUES,U_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daVect,fields->U,INSERT_VALUES,U_localVec);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daVect,U_localVec,&U_array);CHKERRQ(ierr);    
  /*
    get theta_array, thetaRef_array
  */
  ierr = DMGetLocalVector(ctx->daScal,&theta_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,fields->theta,INSERT_VALUES,theta_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,fields->theta,INSERT_VALUES,theta_localVec);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,theta_localVec,&theta_array);CHKERRQ(ierr);    

  ierr = DMGetLocalVector(ctx->daScal,&thetaRef_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,fields->thetaRef,INSERT_VALUES,thetaRef_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,fields->thetaRef,INSERT_VALUES,thetaRef_localVec);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,thetaRef_localVec,&thetaRef_array);CHKERRQ(ierr);    
  /*
    get pressure_array, pressureRef_array
  */
  ierr = DMGetLocalVector(ctx->daScal,&pressure_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,fields->pressure,INSERT_VALUES,pressure_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,fields->pressure,INSERT_VALUES,pressure_localVec);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,pressure_localVec,&pressure_array);CHKERRQ(ierr);    

  ierr = DMGetLocalVector(ctx->daScal,&pressureRef_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,fields->pressureRef,INSERT_VALUES,pressureRef_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,fields->pressureRef,INSERT_VALUES,pressureRef_localVec);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,pressureRef_localVec,&pressureRef_array);CHKERRQ(ierr);  
  /*
    get Gc_array
  */
  ierr = DMGetLocalVector(ctx->daScalCell,&Gc_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScalCell,ctx->matprop->VecGc,INSERT_VALUES,Gc_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScalCell,ctx->matprop->VecGc,INSERT_VALUES,Gc_localVec);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScalCell,Gc_localVec,&Gc_array);CHKERRQ(ierr);    
  /*
    get local mat and RHS
  */
  ierr = PetscMalloc3(nrow,PetscReal,&RHS_local,
                      nrow * nrow,PetscReal,&K_local,
                      nrow,MatStencil,&row);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ctx->daScal,&RHS_localVec);CHKERRQ(ierr);
  ierr = VecSet(RHS_localVec,0.);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,RHS_localVec,&RHS_array);CHKERRQ(ierr);    

  /*
    loop through all elements (ei,ej)
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
        for (l = 0; l < nrow * nrow; l++) {
          K_local[l] = 0.;
        }
        ierr = VF_MatVAT2Surface3D_local(K_local,&ctx->matprop[ctx->layer[ek]],&ctx->vfprop,&ctx->e3D,Gc_array[ek][ej][ei]);CHKERRQ(ierr);
        switch (ctx->unilateral) {
          case UNILATERAL_NONE:
            ierr = VF_MatVCoupling3D_local(K_local,U_array,theta_array,thetaRef_array,
                                           pressure_array,pressureRef_array,
                                           &ctx->matprop[ctx->layer[ek]],&ctx->vfprop,ek,ej,ei,
                                           &ctx->e3D);CHKERRQ(ierr);
            break;
          case UNILATERAL_SHEARONLY:
            ierr = VF_MatVCouplingShearOnly3D_local(K_local,U_array,theta_array,thetaRef_array,
                                                    pressure_array,pressureRef_array,
                                                    &ctx->matprop[ctx->layer[ek]],&ctx->vfprop,ek,ej,ei,
                                                    &ctx->e3D);CHKERRQ(ierr);
            break;
        }
        /*
         Generate array of grid indices in the linear system's ordering.
         i.e. tells MatSetValuesStencil where to store values from  K_local
         if the element is indexed by (ek, ej, ei), the associated degrees of freedom
         have indices (ek ... ek + nphiz), (ej .. ej + nphiy), (ei ... ei + nphix)
        */
        
        for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
          for (j = 0; j < ctx->e3D.nphiy; j++) {
            for (i = 0; i < ctx->e3D.nphix; i++,l++) {
              row[l].i = ei + i; row[l].j = ej + j; row[l].k = ek + k; row[l].c = 0;
            }
          }
        }
        
        /*
          Add local stiffness matrix to global stiffness natrix
        */
        ierr = MatSetValuesStencil(K,nrow,row,nrow,row,K_local,ADD_VALUES);CHKERRQ(ierr);
        /*
          Accumulate Surface energy contribution to RHS
        */
        for (l = 0; l < nrow; l++) RHS_local[l] = 0.;
        ierr = VF_RHSVAT2Surface3D_local(RHS_local,&ctx->matprop[ctx->layer[ek]],&ctx->vfprop,&ctx->e3D,Gc_array[ek][ej][ei]);CHKERRQ(ierr);
        if (ctx->hasCrackPressure) {
          ierr = VF_RHSVPressure3D_local(RHS_local,U_array,pressure_array,pressureRef_array,
                                        &ctx->matprop[ctx->layer[ek]],&ctx->vfprop,ek,ej,ei,
                                        &ctx->e3D);CHKERRQ(ierr);
        }
        
        for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
          for (j = 0; j < ctx->e3D.nphiy; j++) {
            for (i = 0; i < ctx->e3D.nphix; i++,l++) {
              RHS_array[ek+k][ej+j][ei+i] += RHS_local[l];
            }
          }
        }
        /*
          Jump to next element
        */
      }
    }
  }
  /*
    Global Assembly and Boundary Conditions
  */
  /*
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatApplyDirichletBC(K,ctx->daScal,&ctx->bcV[0]);CHKERRQ(ierr);
  */
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(ctx->daScal,RHS_localVec,&RHS_array);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(ctx->daScal,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ctx->daScal,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&RHS_localVec);CHKERRQ(ierr);
  /*
  ierr = VecApplyDirichletBC(RHS,fields->V,&ctx->bcV[0]);CHKERRQ(ierr);
  */
  /*
    Cleanup
  */
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,U_localVec,&U_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daVect,&U_localVec);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,theta_localVec,&theta_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&theta_localVec);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,thetaRef_localVec,&thetaRef_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&thetaRef_localVec);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,pressure_localVec,&pressure_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&pressure_localVec);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,pressureRef_localVec,&pressureRef_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&pressureRef_localVec);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(ctx->daScalCell,Gc_localVec,&Gc_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScalCell,&Gc_localVec);CHKERRQ(ierr);
  
  ierr = PetscFree3(RHS_local,K_local,row);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_SurfaceEnergy3D_local"
/*
  VF_SurfaceEnergy3D_local: 

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_SurfaceEnergy3D_local(PetscReal *SurfaceEnergy_local,PetscReal ***v_array,MatProp *matprop,VFProp *vfprop,PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e,PetscReal Gc)
{
  PetscInt       g,i,j,k;
  PetscReal      *v_elem,*gradv_elem[3];
  PetscErrorCode ierr;
  PetscReal      coef = Gc / vfprop->atCv * .25;

  PetscFunctionBegin;
  ierr = PetscMalloc4(e->ng,PetscReal,&v_elem,e->ng,PetscReal,&gradv_elem[0],e->ng,PetscReal,&gradv_elem[1],e->ng,PetscReal,&gradv_elem[2]);CHKERRQ(ierr);
  
  for (g = 0; g < e->ng; g++) {
    v_elem[g] = 0.;
    gradv_elem[0][g] = 0.;
    gradv_elem[1][g] = 0.;
    gradv_elem[2][g] = 0.;
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (g = 0; g < e->ng; g++) {
          v_elem[g] += v_array[ek+k][ej+j][ei+i] * e->phi[k][j][i][g];
          gradv_elem[0][g] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][0][g];
          gradv_elem[1][g] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][1][g];
          gradv_elem[2][g] += v_array[ek+k][ej+j][ei+i] * e->dphi[k][j][i][2][g];
        }
      }
    }
  }
  for (g = 0; g < e->ng; g++) {
    *SurfaceEnergy_local += e->weight[g] * ((1. - v_elem[g]) * (1. - v_elem[g]) / vfprop->epsilon 
      + (gradv_elem[0][g] * gradv_elem[0][g] + gradv_elem[1][g] * gradv_elem[1][g] + gradv_elem[2][g] * gradv_elem[2][g]) * vfprop->epsilon) * coef;
  }
  ierr = PetscFree4(v_elem,gradv_elem[0],gradv_elem[1],gradv_elem[2]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_VEnergy3D"
/*
  VF_VEnergy3D

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_VEnergy3D(PetscReal *SurfaceEnergy,VFFields *fields,VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       ei,ej,ek;
  Vec            v_localVec;
  PetscReal      ***v_array;
  PetscReal      mySurfaceEnergy=0.;
  PetscReal      hx,hy,hz;
  PetscReal      ****coords_array;
  Vec            Gc_localVec;
  PetscReal      ***Gc_array;

  PetscFunctionBegin;
  
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  
  /*
    Get coordinates, if necessary
  */
  ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ctx->daScal,&v_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_localVec);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,v_localVec,&v_array);CHKERRQ(ierr);    
  ierr = DMGetLocalVector(ctx->daScalCell,&Gc_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScalCell,ctx->matprop->VecGc,INSERT_VALUES,Gc_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScalCell,ctx->matprop->VecGc,INSERT_VALUES,Gc_localVec);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScalCell,Gc_localVec,&Gc_array);CHKERRQ(ierr);        
  
  for (ek = zs; ek < zs + zm; ek++) {
    for (ej = ys; ej < ys + ym; ej++) {
      for (ei = xs; ei < xs + xm; ei++) {
        hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
        hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
        hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
        ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
        ierr = VF_SurfaceEnergy3D_local(&mySurfaceEnergy,v_array,&ctx->matprop[ctx->layer[ek]],&ctx->vfprop,ek,ej,ei,&ctx->e3D,Gc_array[ek][ej][ei]);CHKERRQ(ierr);
      }
    }
  }
  *SurfaceEnergy = 0.;
  ierr = MPI_Reduce(&mySurfaceEnergy,SurfaceEnergy,1,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(ctx->daScal,v_localVec,&v_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&v_localVec);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScalCell,Gc_localVec,&Gc_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScalCell,&Gc_localVec);CHKERRQ(ierr);  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_IrrevApplyEQ"
/*
  VF_IrrevApplyEQ: Apply irreversibility conditions using truncation and equality constraints

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_IrrevApplyEQ(Mat K,Vec RHS,Vec V,Vec VIrrev,VFProp *vfprop,VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscInt       myirrevnum = 0;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       i,j,k,l = 0;
  PetscReal      ***V_array,***VIrrev_array,***RHS_array;
  PetscReal      one = 1.;
  MatStencil     *row;
  PetscFunctionBegin;
  
  ierr = DMDAGetInfo(ctx->daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(ctx->daScal,V,&V_array);CHKERRQ(ierr);        
  ierr = DMDAVecGetArray(ctx->daScal,VIrrev,&VIrrev_array);CHKERRQ(ierr);        
  ierr = DMDAVecGetArray(ctx->daScal,RHS,&RHS_array);CHKERRQ(ierr);        
  
  /* 
    first pass: RHS, V, and count
  */
  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {    
        if (VIrrev_array[k][j][i] <= vfprop->irrevtol) {
          myirrevnum++;
          RHS_array[k][j][i] = 0.;
          V_array[k][j][i] = 0.;
        }
      }
    }
  }
  ierr = PetscMalloc(myirrevnum * sizeof(MatStencil),&row);CHKERRQ(ierr);
  /* 
    second pass: Matrix
  */
  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {    
        if (VIrrev_array[k][j][i] <= vfprop->irrevtol) {
          row[l].i = i; row[l].j = j; row[l].k = k; row[l].c = 0; 
          l++; 
        }
      }
    }
  }
  ierr = MatZeroRowsStencil(K,myirrevnum,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscFree(row);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(ctx->daScal,V,&V_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,VIrrev,&VIrrev_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,RHS,&RHS_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_StepV"
/*
  VF_StepV

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_StepV(VFFields *fields,VFCtx *ctx)
{
  PetscErrorCode      ierr;
  SNESConvergedReason reason;
  PetscInt            its;
  PetscReal           Vmin,Vmax;
  
  PetscFunctionBegin;
  ierr = VF_VAssembly3D(ctx->KV,ctx->RHSV,fields,ctx);CHKERRQ(ierr);
  /*
    Take care of irreversibility
  */
  ierr = VF_IrrevApplyEQ(ctx->KV,ctx->RHSV,fields->V,fields->VIrrev,&ctx->vfprop,ctx);CHKERRQ(ierr);
  
  if (ctx->verbose > 1) {
    ierr = VecView(ctx->RHSV,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = VecView(fields->U,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = MatView(ctx->KV,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  ierr = SNESSetFunction(ctx->snesV,ctx->VResidual,VF_VResidual,ctx);CHKERRQ(ierr);
  ierr = SNESSetJacobian(ctx->snesV,ctx->JacV,ctx->JacV,VF_VIJacobian,ctx);CHKERRQ(ierr);
  if (ctx->verbose > 1) {
    ierr = SNESMonitorSet(ctx->snesV,VF_VSNESMonitor,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  }
  ierr = SNESSolve(ctx->snesV,PETSC_NULL,fields->V);CHKERRQ(ierr);
  if (ctx->verbose > 1) {
    ierr = VecView(fields->V,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  
  ierr = SNESGetConvergedReason(ctx->snesV,&reason);CHKERRQ(ierr);
  if (reason < 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"[ERROR] snesV diverged with reason %d\n",(int)reason);CHKERRQ(ierr);
  } else {
    ierr = SNESGetIterationNumber(ctx->snesV,&its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"      snesV converged in %d iterations %d.\n",(int)its,(int)reason);CHKERRQ(ierr);
  }
  /* Get Min / Max of V */
  ierr = VecMin(fields->V,PETSC_NULL,&Vmin);CHKERRQ(ierr);
  ierr = VecMax(fields->V,PETSC_NULL,&Vmax);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"      V min / max:     %e %e\n",Vmin,Vmax);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_VResidual"
extern PetscErrorCode VF_VResidual(SNES snes,Vec V,Vec residual,void *user)
{
  PetscErrorCode ierr;
  VFCtx     *ctx=(VFCtx*)user;
  
  PetscFunctionBegin;
  ierr = VecSet(residual,0.0);CHKERRQ(ierr);
  ierr = MatMult(ctx->KV,V,residual);CHKERRQ(ierr); 
  ierr = VecAXPY(residual,-1.0,ctx->RHSV);CHKERRQ(ierr);
  ierr = ResidualApplyDirichletBC(residual,V,V,ctx->bcV);CHKERRQ(ierr);
  /*
    UGLY HACK:
    We can pass V for the BC, because we know that the only BC used in the V problem are 0 and 1, so 
    the 3rd argument is never used
  */
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_VIJacobian"
extern PetscErrorCode VF_VIJacobian(SNES snes,Vec V,Mat *Jac,Mat *Jac1,MatStructure *str,void *user)
{
  PetscErrorCode ierr;
  VFCtx       *ctx=(VFCtx*)user;
  
  PetscFunctionBegin;
  *str = SAME_NONZERO_PATTERN;
  ierr = MatCopy(ctx->KV,*Jac,*str);
  ierr = MatAssemblyBegin(*Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (*Jac != *Jac1) {
    ierr = MatAssemblyBegin(*Jac1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*Jac1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
  ierr = MatApplyDirichletBCRowCol(*Jac1,ctx->daScal,ctx->bcV);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "VF_VSNESMonitor"
extern PetscErrorCode VF_VSNESMonitor(SNES snes,PetscInt its,PetscReal fnorm,void* ptr)
{
  PetscErrorCode ierr;
  PetscReal      norm,vmax,vmin;
  MPI_Comm       comm;
  Vec       V;
  
  PetscFunctionBegin;
  ierr = SNESGetSolution(snes,&V);CHKERRQ(ierr);  
  ierr = VecNorm(V,NORM_1,&norm);CHKERRQ(ierr);
  ierr = VecMax(V,PETSC_NULL,&vmax);CHKERRQ(ierr);
  ierr = VecMin(V,PETSC_NULL,&vmin);CHKERRQ(ierr);
  ierr = PetscObjectGetComm((PetscObject)snes,&comm);CHKERRQ(ierr);
  ierr = PetscPrintf(comm,"V_snes_iter_step %D :solution norm = %G, max sol. value  = %G, min sol. value = %G\n",its,norm,vmax,vmin);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}





