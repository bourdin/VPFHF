#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFV.h"
#include "VFU.h"


#undef __FUNCT__
#define __FUNCT__ "BCUInit"
/*
  BCUInit

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/

extern PetscErrorCode BCUInit(BC *BC,VFPreset preset)
{
  PetscInt       c,loc;

  PetscFunctionBegin;
  for (loc = 0; loc < 6; loc++) {
    for (c=0; c < 3; c++) {
      BC[c].face[loc] = NONE;
    }
  }
  for (loc = 0; loc < 12; loc++) {
    for (c=0; c < 3; c++) {
      BC[c].edge[loc] = NONE;
    }
  }
  for (loc = 0; loc < 8; loc++) {
    for (c=0; c < 3; c++) {
      BC[c].vertex[loc] = NONE;
    }
  }
  switch (preset) {
    case SYMXY:
      /*
        blocking vertical displacement on the plane z=z_min, 
        + symmetry with respect to the x=0 and y=0 planes
      */
      BC[0].face[X0] = ZERO;
      BC[1].face[Y0] = ZERO;
      BC[2].face[Z1] = ZERO;
      break;
    case SYMX:
      /*
        blocking vertical displacement on the plane z=z_max, 
        + symmetry with respect to the x=0 plane 
      */
      BC[0].face[X0]       = ZERO;
      BC[1].vertex[X0Y0Z1] = ZERO;
      BC[2].face[Z1]       = ZERO;
    case SYMY:
      /*
        blocking vertical displacement on the plane z=z_max, 
        + symmetry with respect to the y=0 plane 
      */
      BC[0].vertex[X0Y0Z1] = ZERO;
      BC[1].face[Y0]       = ZERO;
      BC[2].face[Z1]       = ZERO;
    case NOSYM: 
      /*
        blocking vertical displacement on the plane z=z_max, 
        + rigid motions
      */
      BC[0].vertex[X0Y0Z1] = ZERO;
      BC[1].vertex[X0Y0Z1] = ZERO;
      BC[1].vertex[X1Y0Z1] = ZERO;
      BC[2].face[Z1]       = ZERO;
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
#define __FUNCT__ "BCUUpdate"
/*
  BCUUpdate:  Update boundary condition flag for U when prescribed displacements are applied
              Typically, this is used to account for insitu stresses by first computing the 
              boundary displacement using the reference pressure and porosity.

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode BCUUpdate(BC *BC,VFPreset preset)
{
  PetscInt       c;

  PetscFunctionBegin;
  switch (preset) {
     case SYMXY:
       /*
         symmetry with respect to the x=0 and y=0 planes, all other faces from data dile
       */
       BC[0].face[X0] = ZERO;
       BC[0].face[X1] = FIXED;   BC[1].face[X1] = FIXED;   BC[2].face[X1] = FIXED; 
       BC[1].face[Y0] = ZERO;
       BC[0].face[Y1] = FIXED;   BC[1].face[Y1] = FIXED;   BC[2].face[Y1] = FIXED; 
       BC[0].face[Z0] = FIXED;   BC[1].face[Z0] = FIXED;   BC[2].face[Z0] = FIXED; 
       BC[0].face[Z1] = FIXED;   BC[1].face[Z1] = FIXED;   BC[2].face[Z1] = FIXED; 
       break;
    case SYMX:
       /*
         symmetry with respect to the x=0 plane, all other faces from data dile
       */
       BC[0].face[X0] = ZERO;
       BC[0].face[X1] = FIXED;   BC[1].face[X1] = FIXED;   BC[2].face[X1] = FIXED; 
       BC[0].face[Y0] = FIXED;   BC[1].face[Y0] = FIXED;   BC[2].face[Y0] = FIXED; 
       BC[0].face[Y1] = FIXED;   BC[1].face[Y1] = FIXED;   BC[2].face[Y1] = FIXED; 
       BC[0].face[Z0] = FIXED;   BC[1].face[Z0] = FIXED;   BC[2].face[Z0] = FIXED; 
       BC[0].face[Z1] = FIXED;   BC[1].face[Z1] = FIXED;   BC[2].face[Z1] = FIXED; 
       break;
    case SYMY:
       /*
         symmetry with respect to the x=0 plane, all other faces from data dile
       */
       BC[0].face[X0] = FIXED;   BC[1].face[X0] = FIXED;   BC[2].face[X0] = FIXED; 
       BC[0].face[X1] = FIXED;   BC[1].face[X1] = FIXED;   BC[2].face[X1] = FIXED; 
       BC[0].face[Y0] = ZERO;
       BC[0].face[Y1] = FIXED;   BC[1].face[Y1] = FIXED;   BC[2].face[Y1] = FIXED; 
       BC[0].face[Z0] = FIXED;   BC[1].face[Z0] = FIXED;   BC[2].face[Z0] = FIXED; 
       BC[0].face[Z1] = FIXED;   BC[1].face[Z1] = FIXED;   BC[2].face[Z1] = FIXED; 
       break;
    case NOSYM: 
       /* 
          Reading all boundary displacements from file
       */
       for (c = 0; c < 3; c++) {
         BC[c].face[X0] = FIXED;
         BC[c].face[X1] = FIXED;
         BC[c].face[Y0] = FIXED;
         BC[c].face[Y1] = FIXED;
         BC[c].face[Z0] = FIXED;
         BC[c].face[Z1] = FIXED;
       }
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
#define __FUNCT__ "ElasticEnergyDensity3D_local"
/*
  ElasticEnergyDensity3D_local: computes the local contribution of the elastic energy desity
                                at each Gauss point.
  
  The elastic energy density is defiend as 
      1/2 A \epsilon : \epsilon 
  where epsilon is the inelastic strain:
       \epsilon = e(u) - \alpha (\theta-\theta_0) Id - \beta (p-p_0) A^{-1}Id
                = e(u) - \alpha (\theta-\theta_0) Id - \beta (p-p_0) / (3\lambda + 2\u) Id
  and 
     * ":" denotes tensor contraction, i.e. M:N = tr M^tN
     * A is the Hooke's law, \lambda and \mu are the lam\'e coefficients
     * \alpha is the thermal expansion coefficient
     * \beta is the Biot coefficient
     * the index "0" denotes a reference field (temperature or pressure)

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode ElasticEnergyDensity3D_local(PetscReal *ElasticEnergyDensity_local,
                                                   PetscReal ****u_array,
                                                   PetscReal ***theta_array,PetscReal ***thetaRef_array,
                                                   PetscReal ***pressure_array,PetscReal ***pressureRef_array,
                                                   MatProp *matprop,PetscInt ek,PetscInt ej,PetscInt ei,
                                                   CartFE_Element3D *e)
{
  PetscErrorCode ierr;
  PetscReal      *epsilon11_elem,*epsilon22_elem,*epsilon33_elem,*epsilon12_elem,*epsilon23_elem,*epsilon13_elem;
  PetscReal      *sigma11_elem,*sigma22_elem,*sigma33_elem,*sigma12_elem,*sigma23_elem,*sigma13_elem;
  PetscInt       i,j,k,g;
  PetscReal      lambda,mu,alpha,coefbeta;
  /*
  PetscReal      myElasticEnergyDensity = 0;
  */

  PetscFunctionBegin;
  lambda   = matprop->lambda;
  mu       = matprop->mu;
  alpha    = matprop->alpha;
  coefbeta = matprop->beta / (3.*lambda + 2.*mu);
    
  ierr = PetscMalloc3(e->ng,PetscReal,&sigma11_elem,e->ng,PetscReal,&sigma22_elem,e->ng,PetscReal,&sigma33_elem);CHKERRQ(ierr);
  ierr = PetscMalloc3(e->ng,PetscReal,&sigma12_elem,e->ng,PetscReal,&sigma23_elem,e->ng,PetscReal,&sigma13_elem);CHKERRQ(ierr);
  ierr = PetscMalloc3(e->ng,PetscReal,&epsilon11_elem,e->ng,PetscReal,&epsilon22_elem,e->ng,PetscReal,&epsilon33_elem);CHKERRQ(ierr);
  ierr = PetscMalloc3(e->ng,PetscReal,&epsilon12_elem,e->ng,PetscReal,&epsilon23_elem,e->ng,PetscReal,&epsilon13_elem);CHKERRQ(ierr);

  for (g = 0; g < e->ng; g++) {
    epsilon11_elem[g] = 0;
    epsilon22_elem[g] = 0;
    epsilon33_elem[g] = 0;
    epsilon12_elem[g] = 0;
    epsilon23_elem[g] = 0;
    epsilon13_elem[g] = 0;
    ElasticEnergyDensity_local[g] = 0.;
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (g = 0; g < e->ng; g++) {
          epsilon11_elem[g] +=  e->dphi[k][j][i][0][g] * u_array[ek+k][ej+j][ei+i][0] 
                              - alpha * e->phi[k][j][i][g] 
                                * (theta_array[ek+k][ej+j][ei+i]-thetaRef_array[ek+k][ej+j][ei+i])
                              - coefbeta * e->phi[k][j][i][g] 
                                * (pressure_array[ek+k][ej+j][ei+i]-pressureRef_array[ek+k][ej+j][ei+i]);
          epsilon22_elem[g] +=  e->dphi[k][j][i][1][g] * u_array[ek+k][ej+j][ei+i][1] 
                              - alpha * e->phi[k][j][i][g] 
                                * (theta_array[ek+k][ej+j][ei+i]-thetaRef_array[ek+k][ej+j][ei+i])
                              - coefbeta * e->phi[k][j][i][g] 
                                * (pressure_array[ek+k][ej+j][ei+i]-pressureRef_array[ek+k][ej+j][ei+i]);
          epsilon33_elem[g] +=  e->dphi[k][j][i][2][g] * u_array[ek+k][ej+j][ei+i][2] 
                              - alpha * e->phi[k][j][i][g] 
                                * (theta_array[ek+k][ej+j][ei+i]-thetaRef_array[ek+k][ej+j][ei+i])
                              - coefbeta * e->phi[k][j][i][g] 
                                * (pressure_array[ek+k][ej+j][ei+i]-pressureRef_array[ek+k][ej+j][ei+i]);

          epsilon12_elem[g] += (e->dphi[k][j][i][0][g] * u_array[ek+k][ej+j][ei+i][1] 
                              + e->dphi[k][j][i][1][g] * u_array[ek+k][ej+j][ei+i][0]) * .5;
          epsilon23_elem[g] += (e->dphi[k][j][i][1][g] * u_array[ek+k][ej+j][ei+i][2] 
                              + e->dphi[k][j][i][2][g] * u_array[ek+k][ej+j][ei+i][1]) * .5;
          epsilon13_elem[g] += (e->dphi[k][j][i][0][g] * u_array[ek+k][ej+j][ei+i][2] 
                              + e->dphi[k][j][i][2][g] * u_array[ek+k][ej+j][ei+i][0]) * .5;
        }
      }
    }  
  }
  ierr = PetscLogFlops(33 * e->ng * e->nphix * e->nphiy * e->nphiz);CHKERRQ(ierr);
  
  for (g = 0; g < e->ng; g++) {
    sigma11_elem[g] = (lambda + 2.*mu) * epsilon11_elem[g] + lambda * epsilon22_elem[g] + lambda * epsilon33_elem[g];
    sigma22_elem[g] = lambda * epsilon11_elem[g] + (lambda + 2.*mu) * epsilon22_elem[g] + lambda * epsilon33_elem[g];
    sigma33_elem[g] = lambda * epsilon11_elem[g] + lambda * epsilon22_elem[g] + (lambda + 2.*mu) * epsilon33_elem[g];
    sigma12_elem[g] = 2. * mu * epsilon12_elem[g];
    sigma23_elem[g] = 2. * mu * epsilon23_elem[g];
    sigma13_elem[g] = 2. * mu * epsilon13_elem[g];
    ElasticEnergyDensity_local[g] = (sigma11_elem[g] * epsilon11_elem[g] 
                                   + sigma22_elem[g] * epsilon22_elem[g] 
                                   + sigma33_elem[g] * epsilon33_elem[g]) * .5 
                                   + sigma12_elem[g] * epsilon12_elem[g]
                                   + sigma23_elem[g] * epsilon23_elem[g]
                                   + sigma13_elem[g] * epsilon13_elem[g];
  }
  ierr = PetscLogFlops(46 * e->ng);CHKERRQ(ierr);

  ierr = PetscFree3(sigma11_elem,sigma22_elem,sigma33_elem);CHKERRQ(ierr);
  ierr = PetscFree3(sigma12_elem,sigma23_elem,sigma13_elem);CHKERRQ(ierr);
  ierr = PetscFree3(epsilon11_elem,epsilon22_elem,epsilon33_elem);CHKERRQ(ierr);
  ierr = PetscFree3(epsilon12_elem,epsilon23_elem,epsilon13_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ElasticEnergyDensitySphericalDeviatoric3D_local"
/*
  ElasticEnergyDensitySphericalDeviatoric3D_local
  Compute the spherical and deviatoric parts of the elastic energy density in an element
  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode ElasticEnergyDensitySphericalDeviatoric3D_local(PetscReal *ElasticEnergyDensityS_local,PetscReal *ElasticEnergyDensityD_local,
                                                                      PetscReal ****u_array,
                                                                      PetscReal ***theta_array,PetscReal ***thetaRef_array,
                                                                      PetscReal ***pressure_array,PetscReal ***pressureRef_array,
                                                                      MatProp *matprop,PetscInt ek,PetscInt ej,PetscInt ei,
                                                                      CartFE_Element3D *e)
{
  PetscErrorCode ierr;
  PetscReal      *epsilon11_elem,*epsilon22_elem,*epsilon33_elem,*epsilon12_elem,*epsilon23_elem,*epsilon13_elem;
  PetscReal      *sigma11_elem,*sigma22_elem,*sigma33_elem,*sigma12_elem,*sigma23_elem,*sigma13_elem;
  PetscInt       i,j,k,g;
  PetscReal      lambda,mu,alpha,kappa,coefbeta;
  PetscReal      *ElasticEnergyDensity_local;
  
  PetscFunctionBegin;
  lambda = matprop->lambda;
  mu     = matprop->mu;
  alpha  = matprop->alpha;
  kappa  = lambda + 2.* mu / 3.;
  coefbeta = matprop->beta / (3.*lambda + 2.*mu);
  
  ierr = PetscMalloc3(e->ng,PetscReal,&sigma11_elem,e->ng,PetscReal,&sigma22_elem,e->ng,PetscReal,&sigma33_elem);CHKERRQ(ierr);
  ierr = PetscMalloc3(e->ng,PetscReal,&sigma12_elem,e->ng,PetscReal,&sigma23_elem,e->ng,PetscReal,&sigma13_elem);CHKERRQ(ierr);
  ierr = PetscMalloc3(e->ng,PetscReal,&epsilon11_elem,e->ng,PetscReal,&epsilon22_elem,e->ng,PetscReal,&epsilon33_elem);CHKERRQ(ierr);
  ierr = PetscMalloc3(e->ng,PetscReal,&epsilon12_elem,e->ng,PetscReal,&epsilon23_elem,e->ng,PetscReal,&epsilon13_elem);CHKERRQ(ierr);
   
  ierr = PetscMalloc(e->ng * sizeof(PetscReal),&ElasticEnergyDensity_local);CHKERRQ(ierr);
  for (g = 0; g < e->ng; g++) {
    epsilon11_elem[g] = 0;
    epsilon22_elem[g] = 0;
    epsilon33_elem[g] = 0;
    epsilon12_elem[g] = 0;
    epsilon23_elem[g] = 0;
    epsilon13_elem[g] = 0;
    ElasticEnergyDensity_local[g]  = 0.;
    ElasticEnergyDensityS_local[g] = 0.;
    ElasticEnergyDensityD_local[g] = 0.;
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (g = 0; g < e->ng; g++) {
          epsilon11_elem[g] +=  e->dphi[k][j][i][0][g] * u_array[ek+k][ej+j][ei+i][0] 
                              - alpha * e->phi[k][j][i][g] 
                                * (theta_array[ek+k][ej+j][ei+i]-thetaRef_array[ek+k][ej+j][ei+i])
                              - coefbeta * e->phi[k][j][i][g] 
                                * (pressure_array[ek+k][ej+j][ei+i]-pressureRef_array[ek+k][ej+j][ei+i]);
          epsilon22_elem[g] +=  e->dphi[k][j][i][1][g] * u_array[ek+k][ej+j][ei+i][1] 
                              - alpha * e->phi[k][j][i][g] 
                                * (theta_array[ek+k][ej+j][ei+i]-thetaRef_array[ek+k][ej+j][ei+i])
                              - coefbeta * e->phi[k][j][i][g] 
                                * (pressure_array[ek+k][ej+j][ei+i]-pressureRef_array[ek+k][ej+j][ei+i]);
          epsilon33_elem[g] +=  e->dphi[k][j][i][2][g] * u_array[ek+k][ej+j][ei+i][2] 
                              - alpha * e->phi[k][j][i][g] 
                                * (theta_array[ek+k][ej+j][ei+i]-thetaRef_array[ek+k][ej+j][ei+i])
                              - coefbeta * e->phi[k][j][i][g] 
                                * (pressure_array[ek+k][ej+j][ei+i]-pressureRef_array[ek+k][ej+j][ei+i]);

          epsilon12_elem[g] += (e->dphi[k][j][i][0][g] * u_array[ek+k][ej+j][ei+i][1] 
                              + e->dphi[k][j][i][1][g] * u_array[ek+k][ej+j][ei+i][0]) * .5;
          epsilon23_elem[g] += (e->dphi[k][j][i][1][g] * u_array[ek+k][ej+j][ei+i][2] 
                              + e->dphi[k][j][i][2][g] * u_array[ek+k][ej+j][ei+i][1]) * .5;
          epsilon13_elem[g] += (e->dphi[k][j][i][0][g] * u_array[ek+k][ej+j][ei+i][2] 
                              + e->dphi[k][j][i][2][g] * u_array[ek+k][ej+j][ei+i][0]) * .5;

        }
      }
    }  
  }
  ierr = PetscLogFlops(30 * e->ng * e->nphix * e->nphiy * e->nphiz);CHKERRQ(ierr);
  
  for (g = 0; g < e->ng; g++) {
    sigma11_elem[g] = (lambda + 2.*mu) * epsilon11_elem[g] + lambda * epsilon22_elem[g] + lambda * epsilon33_elem[g];
    sigma22_elem[g] = lambda * epsilon11_elem[g] + (lambda + 2.*mu) * epsilon22_elem[g] + lambda * epsilon33_elem[g];
    sigma33_elem[g] = lambda * epsilon11_elem[g] + lambda * epsilon22_elem[g] + (lambda + 2.*mu) * epsilon33_elem[g];
    sigma12_elem[g] = 2.*mu * epsilon12_elem[g];
    sigma23_elem[g] = 2.*mu * epsilon23_elem[g];
    sigma13_elem[g] = 2.*mu * epsilon13_elem[g];
    ElasticEnergyDensity_local[g]  = (sigma11_elem[g] * epsilon11_elem[g] 
                                    + sigma22_elem[g] * epsilon22_elem[g] 
                                    + sigma33_elem[g] * epsilon33_elem[g]) * .5 
                                    + sigma12_elem[g] * epsilon12_elem[g]
                                    + sigma23_elem[g] * epsilon23_elem[g]
                                    + sigma13_elem[g] * epsilon13_elem[g];
    ElasticEnergyDensityS_local[g] = kappa * (epsilon11_elem[g] + epsilon22_elem[g] + epsilon33_elem[g]) 
                                           * (epsilon11_elem[g] + epsilon22_elem[g] + epsilon33_elem[g]) * .5;
    ElasticEnergyDensityD_local[g] = ElasticEnergyDensity_local[g] - ElasticEnergyDensityS_local[g];
  }
  ierr = PetscLogFlops(47 * e->ng);CHKERRQ(ierr);
  
  ierr = PetscFree3(sigma11_elem,sigma22_elem,sigma33_elem);CHKERRQ(ierr);
  ierr = PetscFree3(sigma12_elem,sigma23_elem,sigma13_elem);CHKERRQ(ierr);
  ierr = PetscFree3(epsilon11_elem,epsilon22_elem,epsilon33_elem);CHKERRQ(ierr);
  ierr = PetscFree3(epsilon12_elem,epsilon23_elem,epsilon13_elem);CHKERRQ(ierr);
  ierr = PetscFree(ElasticEnergyDensity_local);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_MatU3D_local"
/*
  VF_MatU3D_local

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_MatU3D_local(PetscReal *Mat_local,PetscReal ***v_array,MatProp *matprop,VFProp *vfprop,PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e)
{
  PetscInt       g,i1,i2,j1,j2,k1,k2,c1,c2,l;
  PetscReal      lambda,mu;
  PetscReal      *v_elem;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscMalloc(e->ng * sizeof(PetscReal),&v_elem);CHKERRQ(ierr);
  lambda = matprop->lambda;
  mu     = matprop->mu;
  
  for (g = 0; g < e->ng; g++) v_elem[g] = 0.;
  for (k1 = 0; k1 < e->nphiz; k1++) {
    for (j1 = 0; j1 < e->nphiy; j1++) {
      for (i1 = 0; i1 < e->nphix; i1++) {
        for (g = 0; g < e->ng; g++) {
          v_elem[g] += v_array[ek+k1][ej+j1][ei+i1] * e->phi[k1][j1][i1][g];
        }
      }
    }
  }
  ierr = PetscLogFlops(2 * e->ng * e->nphix * e->nphiy * e->nphiz);CHKERRQ(ierr);

  for (l = 0,k1 = 0; k1 < e->nphiz; k1++) {
    for (j1 = 0; j1 < e->nphiy; j1++) {
      for (i1 = 0; i1 < e->nphix; i1++) {
        for (c1 = 0; c1 < e->dim; c1++) {
          for (k2 = 0; k2 < e->nphiz; k2++) {
            for (j2 = 0; j2 < e->nphiy; j2++) {
              for (i2 = 0; i2 < e->nphix; i2++) {
                for (c2 = 0; c2 < e->dim; c2++,l++) {
                  Mat_local[l] = 0; 
                  for (g = 0; g < e->ng; g++) {
                    if (c1 == c2) {
                      Mat_local[l] += e->weight[g] * ( mu * ( e->dphi[k1][j1][i1][0][g] * e->dphi[k2][j2][i2][0][g] 
                                                            + e->dphi[k1][j1][i1][1][g] * e->dphi[k2][j2][i2][1][g]
                                                            + e->dphi[k1][j1][i1][2][g] * e->dphi[k2][j2][i2][2][g])
                                                     + (lambda + mu) * e->dphi[k1][j1][i1][c1][g] * e->dphi[k2][j2][i2][c2][g])
                                                   * (v_elem[g] * v_elem[g] + vfprop->eta);
                      
                      ierr = PetscLogFlops(15);CHKERRQ(ierr);
                    } else {
                      Mat_local[l] += e->weight[g] * (lambda * e->dphi[k1][j1][i1][c1][g] * e->dphi[k2][j2][i2][c2][g] 
                                                        + mu * e->dphi[k1][j1][i1][c2][g] * e->dphi[k2][j2][i2][c1][g])
                                                   * (v_elem[g] * v_elem[g] + vfprop->eta);
                      ierr = PetscLogFlops(10);CHKERRQ(ierr);
                      
                    }
                  }
                }                                                            
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
#define __FUNCT__ "VF_MatUShearOnly3D_local"
/*
  VF_MatUShearOnly3D_local

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_MatUShearOnly3D_local(PetscReal *Mat_local,PetscReal ***v_array,MatProp *matprop,VFProp *vfprop,PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e)
{
  PetscInt       g,i1,i2,j1,j2,k1,k2,c1,c2,l;
  PetscReal      lambda,mu,kappa;
  PetscReal      *v_elem;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscMalloc(e->ng * sizeof(PetscReal),&v_elem);CHKERRQ(ierr);
  lambda = matprop->lambda;
  mu     = matprop->mu;
  kappa  = lambda + 2. * mu / 3.;
  
  for (g = 0; g < e->ng; g++) v_elem[g] = 0.;
  for (k1 = 0; k1 < e->nphiz; k1++) {
    for (j1 = 0; j1 < e->nphiy; j1++) {
      for (i1 = 0; i1 < e->nphix; i1++) {
        for (g = 0; g < e->ng; g++) {
          v_elem[g] += v_array[ek+k1][ej+j1][ei+i1] * e->phi[k1][j1][i1][g];
        }
      }
    }
  }
  ierr = PetscLogFlops(2 * e->ng * e->nphix * e->nphiy * e->nphiz);CHKERRQ(ierr);

  for (l = 0,k1 = 0; k1 < e->nphiz; k1++) {
    for (j1 = 0; j1 < e->nphiy; j1++) {
      for (i1 = 0; i1 < e->nphix; i1++) {
        for (c1 = 0; c1 < e->dim; c1++) {
          for (k2 = 0; k2 < e->nphiz; k2++) {
            for (j2 = 0; j2 < e->nphiy; j2++) {
              for (i2 = 0; i2 < e->nphix; i2++) {
                for (c2 = 0; c2 < e->dim; c2++,l++) {
                  Mat_local[l] = 0; 
                  for (g = 0; g < e->ng; g++) {
                    if (c1 == c2) {
                      Mat_local[l] += e->weight[g] * ( mu * ( e->dphi[k1][j1][i1][0][g] * e->dphi[k2][j2][i2][0][g] 
                                                            + e->dphi[k1][j1][i1][1][g] * e->dphi[k2][j2][i2][1][g]
                                                            + e->dphi[k1][j1][i1][2][g] * e->dphi[k2][j2][i2][2][g]
                                                            + e->dphi[k1][j1][i1][c1][g] * e->dphi[k2][j2][i2][c2][g] / 3.)
                                                          * (v_elem[g] * v_elem[g] + vfprop->eta)
                                                      + kappa * e->dphi[k1][j1][i1][c1][g] * e->dphi[k2][j2][i2][c2][g]);
                      
                      ierr = PetscLogFlops(17);CHKERRQ(ierr);
                    } else {
                      Mat_local[l] += e->weight[g] * ( mu * ( e->dphi[k1][j1][i1][c2][g] * e->dphi[k2][j2][i2][c1][g] 
                                                            - 2. * e->dphi[k1][j1][i1][c1][g] * e->dphi[k2][j2][i2][c2][g] / 3.)
                                                          * (v_elem[g] * v_elem[g] + vfprop->eta)
                                                      + kappa * e->dphi[k1][j1][i1][c1][g] * e->dphi[k2][j2][i2][c2][g]);
                      ierr = PetscLogFlops(14);CHKERRQ(ierr);
                      
                    }
                  }
                }                                                            
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
#define __FUNCT__ "VF_RHSUCoupling3D_local"
/*
  VF_RHSUCoupling3D_local

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_RHSUCoupling3D_local(PetscReal *RHS_local,PetscReal ***v_array,
                                              PetscReal ***theta_array,PetscReal ***thetaRef_array,
                                              PetscReal ***pressure_array,PetscReal ***pressureRef_array,
                                              MatProp *matprop,VFProp *vfprop,
                                              PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e)
{
  PetscErrorCode ierr;
  PetscInt       l,i,j,k,g,c;
  PetscInt       dim=3;
  PetscReal      *theta_elem,*pressure_elem,*v_elem;
  PetscReal      coefalpha,beta;
  
  PetscFunctionBegin;
  coefalpha = (3.* matprop->lambda + 2. * matprop->mu) * matprop->alpha;
  beta      = matprop->beta;
  
  ierr = PetscMalloc3(e->ng,PetscReal,&theta_elem,e->ng,PetscReal,&pressure_elem,e->ng,PetscReal,&v_elem);CHKERRQ(ierr);
  /*
    Initialize pressure_Elem, theta_Elem and v_elem
  */
  for (g = 0; g < e->ng; g++) {
    theta_elem[g] = 0;
    pressure_elem[g] = 0;
    v_elem[g] = 0.;
  }
  /*
    Compute theta_elem
  */
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (g = 0; g < e->ng; g++) {
          theta_elem[g]    += e->phi[k][j][i][g] 
                              * (theta_array[ek+k][ej+j][ei+i] - thetaRef_array[ek+k][ej+j][ei+i]);
          pressure_elem[g] += e->phi[k][j][i][g] 
                              * (pressure_array[ek+k][ej+j][ei+i] - pressureRef_array[ek+k][ej+j][ei+i]);
          v_elem[g]        += e->phi[k][j][i][g] * v_array[ek+k][ej+j][ei+i];
        }
      }
    }
  }  
  ierr = PetscLogFlops(8 * e->ng * e->nphix * e->nphiy * e->nphiz);CHKERRQ(ierr);
  /*
    Accumulate the contribution of the current element to the local
    version of the RHS
  */
  for (l=0,k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (c = 0; c < dim; c++,l++) {
          for (g = 0; g < e->ng; g++) {
            RHS_local[l] += e->weight[g] * e->dphi[k][j][i][c][g] 
                            * (coefalpha * theta_elem[g] + beta * pressure_elem[g])
                            * (v_elem[g] * v_elem[g] + vfprop->eta);  
          }
        }
      }
    }
  }
  ierr = PetscLogFlops(9 * dim * e->ng * e->nphix * e->nphiy * e->nphiz);CHKERRQ(ierr);
  /*
    Clean up
  */
  ierr = PetscFree3(theta_elem,pressure_elem,v_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VF_RHSUCouplingShearOnly3D_local"
/*
  VF_RHSUCouplingShearOnly3D_local

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_RHSUCouplingShearOnly3D_local(PetscReal *RHS_local,PetscReal ***v_array,
                                                       PetscReal ***theta_array,PetscReal ***thetaRef_array,
                                                       PetscReal ***pressure_array,PetscReal ***pressureRef_array,
                                                       MatProp *matprop,VFProp *vfprop,
                                                       PetscInt ek,PetscInt ej,PetscInt ei,
                                                       CartFE_Element3D *e)
{
  PetscErrorCode ierr;
  PetscInt       l,i,j,k,g,c;
  PetscInt       dim=3;
  PetscReal      *theta_elem,*pressure_elem;
  PetscReal      coefalpha,beta;
  
  PetscFunctionBegin;
  coefalpha = (3.* matprop->lambda + 2. * matprop->mu) * matprop->alpha;
  beta      = matprop->beta;
  ierr = PetscMalloc2(e->ng,PetscReal,&theta_elem,e->ng,PetscReal,&pressure_elem);CHKERRQ(ierr);
  /*
    Initialize theta_Elem and v_elem
  */
  for (g = 0; g < e->ng; g++) {
    theta_elem[g] = 0;
    pressure_elem[g] = 0;
  }
  /*
    Compute theta_elem
  */
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (g = 0; g < e->ng; g++) {
          theta_elem[g]    += e->phi[k][j][i][g] 
                              * (theta_array[ek+k][ej+j][ei+i] - thetaRef_array[ek+k][ej+j][ei+i]);
          pressure_elem[g] += e->phi[k][j][i][g] 
                              * (pressure_array[ek+k][ej+j][ei+i] - pressureRef_array[ek+k][ej+j][ei+i]);
        }
      }
    }
  }  
  ierr = PetscLogFlops(6 * e->ng * e->nphix * e->nphiy * e->nphiz);CHKERRQ(ierr);
  /*
    Accumulate the contribution of the current element to the local
    version of the RHS.
    Note that only sperical terms appear in the RHS so there is no V
  */
  for (l=0,k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (c = 0; c < dim; c++,l++) {
          for (g = 0; g < e->ng; g++) {
            RHS_local[l] += e->weight[g] * e->dphi[k][j][i][c][g] 
                            * (theta_elem[g] * coefalpha + pressure_elem[g] * beta);  
          }
        }
      }
    }
  }
  ierr = PetscLogFlops(6 * dim * e->ng * e->nphix * e->nphiy * e->nphiz);CHKERRQ(ierr);
  /*
    Clean up
  */
  ierr = PetscFree2(theta_elem,pressure_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_RHSUPressure3D_local"
/*
  VF_RHSUPressure3D_local

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_RHSUPressure3D_local(PetscReal *RHS_local,PetscReal ***v_array,
                                              PetscReal ***pressure_array,PetscReal ***pressureRef_array,
                                              MatProp *matprop,VFProp *vfprop,
                                              PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e)
{
  PetscErrorCode ierr;
  PetscInt       l,i,j,k,g,c;
  PetscReal      *pressure_elem,*gradv_elem[3];

  PetscFunctionBegin;
  ierr = PetscMalloc4(e->ng,PetscReal,&pressure_elem,
                      e->ng,PetscReal,&gradv_elem[0],
                      e->ng,PetscReal,&gradv_elem[1],
                      e->ng,PetscReal,&gradv_elem[2]);CHKERRQ(ierr);

  /*
    Compute the projection of the fields in the local base functions basis
  */
  for (g = 0; g < e->ng; g++) {
    pressure_elem[g] = 0;
    for (c=0; c<3; c++) gradv_elem[c][g] = 0.;
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (g = 0; g < e->ng; g++) {
          pressure_elem[g] += e->phi[k][j][i][g] 
                              * (pressure_array[ek+k][ej+j][ei+i] - pressureRef_array[ek+k][ej+j][ei+i]);
          for (c=0; c<3; c++) gradv_elem[c][g] += e->dphi[k][j][i][c][g] * v_array[ek+k][ej+j][ei+i];
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
      for (i = 0; i < e->nphix; i++) {
        for (c = 0; c < 3; c++,l++) {
          for (g = 0; g < e->ng; g++) {
            RHS_local[l] += e->weight[g] * pressure_elem[g] 
                          * gradv_elem[c][g]            
                          * e->phi[k][j][i][g];  
          }
        }
      }
    }
  }
  ierr = PetscLogFlops(12 * e->ng * e->nphix * e->nphiy * e->nphiz);CHKERRQ(ierr);
  /*
    Clean up
  */
  ierr = PetscFree4(pressure_elem,gradv_elem[0],gradv_elem[1],gradv_elem[2]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_RHSUInSituStresses3D_local"
/*
  VF_RHSUInSituStresses3D_local: Accumulates the contribution of surface forces along the face of an element

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_RHSUInSituStresses3D_local(PetscReal *RHS_local,PetscReal ****f_array,PetscInt ek,PetscInt ej,PetscInt ei,FACE face,CartFE_Element3D *e)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,l,c,g;
  PetscReal      *f_elem[3];

  PetscFunctionBegin;
  ierr = PetscMalloc3(e->ng,PetscReal,&f_elem[0],e->ng,PetscReal,&f_elem[1],e->ng,PetscReal,&f_elem[2]);CHKERRQ(ierr);
  /*
    Initialize f_Elem
  */
  for (c = 0; c < e->dim; c++) {
    for (g = 0; g < e->ng; g++) {
      f_elem[c][g] = 0;
    }
  }
  /*
    Compute f_elem
  */
  switch (face) {
    case X0:
      for (k = 0; k < e->nphiz; k++) {
        for (j = 0; j < e->nphiy; j++) {
          for (i = 0; i < e->nphix; i++) {
            for (c = 0; c < e->dim; c++) {
              for (g = 0; g < e->ng; g++) {
                f_elem[c][g] += e->phi[k][j][i][g] * f_array[ek+k][ej+j][ei][c] / e->lx;
              }
            }
          }
        }
      }
      break;
    case X1:
      for (k = 0; k < e->nphiz; k++) {
        for (j = 0; j < e->nphiy; j++) {
          for (i = 0; i < e->nphix; i++) {
            for (c = 0; c < e->dim; c++) {
              for (g = 0; g < e->ng; g++) {
                f_elem[c][g] += e->phi[k][j][i][g] * f_array[ek+k][ej+j][ei+e->nphix-1][c] / e->lx;
              }
            }
          }
        }
      }
      break;
    case Y0:
      for (k = 0; k < e->nphiz; k++) {
        for (j = 0; j < e->nphiy; j++) {
          for (i = 0; i < e->nphix; i++) {
            for (c = 0; c < e->dim; c++) {
              for (g = 0; g < e->ng; g++) {
                f_elem[c][g] += e->phi[k][j][i][g] * f_array[ek+k][ej][ei+i][c] / e->ly;
              }
            }
          }
        }
      }
      break;
    case Y1:
      for (k = 0; k < e->nphiz; k++) {
        for (j = 0; j < e->nphiy; j++) {
          for (i = 0; i < e->nphix; i++) {
            for (c = 0; c < e->dim; c++) {
              for (g = 0; g < e->ng; g++) {
                f_elem[c][g] += e->phi[k][j][i][g] * f_array[ek+k][ej+e->nphiy-1][ei+i][c] / e->ly;
              }
            }
          }
        }
      }
      break;
    case Z0:
      for (k = 0; k < e->nphiz; k++) {
        for (j = 0; j < e->nphiy; j++) {
          for (i = 0; i < e->nphix; i++) {
            for (c = 0; c < e->dim; c++) {
              for (g = 0; g < e->ng; g++) {
                f_elem[c][g] += e->phi[k][j][i][g] * f_array[ek][ej+j][ei+i][c] / e->lz;
              }
            }
          }
        }
      }
      break;
    case Z1:
      for (k = 0; k < e->nphiz; k++) {
        for (j = 0; j < e->nphiy; j++) {
          for (i = 0; i < e->nphix; i++) {
            for (c = 0; c < e->dim; c++) {
              for (g = 0; g < e->ng; g++) {
                f_elem[c][g] += e->phi[k][j][i][g] * f_array[ek+e->nphiz-1][ej+j][ei+i][c] / e->lz;
              }
            }
          }
        }
      }
      break;
  }
  /* 
    Accumulate
  */
  for (l=0,k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (c = 0; c < e->dim; c++,l++) {
          for (g = 0; g < e->ng; g++) {
            RHS_local[l] += e->weight[g] * e->phi[k][j][i][g] * f_elem[c][g];
          }
        }
      }
    }
  }
  /*
    Clean up
  */
  ierr = PetscFree3(f_elem[0],f_elem[1],f_elem[2]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VF_UAssembly3D"
/*
  VF_UAssembly3D

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_UAssembly3D(Mat K,Vec RHS,VFFields *fields,VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       ei,ej,ek,i,j,k,c,l;
  PetscInt       dim=3;
  PetscInt       nrow = dim * ctx->e3D.nphix * ctx->e3D.nphiy * ctx->e3D.nphiz;
  Vec            RHS_localVec,V_localVec;
  Vec            theta_localVec,thetaRef_localVec;
  Vec            pressure_localVec,pressureRef_localVec;
  PetscReal      ****RHS_array,***v_array;
  PetscReal      ***theta_array,***thetaRef_array;
  PetscReal      ***pressure_array,***pressureRef_array;
  PetscReal      *RHS_local;
  PetscReal      *K_local;
  MatStencil     *row;
  PetscReal      hx,hy,hz;
  PetscReal      ****coords_array;
  PetscReal      ****f_array;
  Vec            f_localVec;
  FACE           face;
  PetscReal      z;
  int            stresscomp[3];
  PetscReal      stressdir[3];
  PetscReal      stressmag;
  PetscReal      BBmin[3],BBmax[3];
  
  PetscFunctionBegin;
  //ierr = PetscLogStagePush(ctx->vflog.VF_UAssemblyStage);CHKERRQ(ierr);
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  ierr = MatZeroEntries(K);CHKERRQ(ierr);
  ierr = VecSet(RHS,0.);CHKERRQ(ierr);
  /*
    Get coordinates
  */
  ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  /*
    Get bounding box from petsc DA
  */
  ierr = DMDAGetBoundingBox(ctx->daVect,BBmin,BBmax);CHKERRQ(ierr);
  /*
    get v_array
  */
  ierr = DMGetLocalVector(ctx->daScal,&V_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,V_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,V_localVec);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,V_localVec,&v_array);CHKERRQ(ierr);    
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
    Allocating multi-dimensional vectors in C is a pain, so for the in-situ stresses / surface stresses, I get a 
    3 component Vec for the external forces, then get a full array. This is a bit wastefull since we never use the 
    internal values...
    Note that we cannot fully initialize it since we need the normal direction to each face.
  */
  if (ctx->hasInsitu) {
    ierr = DMGetLocalVector(ctx->daVect,&f_localVec);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(ctx->daVect,f_localVec,&f_array);CHKERRQ(ierr); 
  }
  /*
    get local mat and RHS
  */
  ierr = PetscMalloc3(nrow,PetscReal,&RHS_local,
                      nrow * nrow,PetscReal,&K_local,
                      nrow,MatStencil,&row);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ctx->daVect,&RHS_localVec);CHKERRQ(ierr);
  ierr = VecSet(RHS_localVec,0.);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daVect,RHS_localVec,&RHS_array);CHKERRQ(ierr);    

  /*
    loop through all elements (ei,ej)
  */
  for (ek = zs; ek < zs + zm; ek++) {
    for (ej = ys; ej < ys+ym; ej++) {
      for (ei = xs; ei < xs+xm; ei++) {
        hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
        hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
        hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
        ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
        /* 
          Compute and accumulate the contribution of the local stiffness matrix to the global stiffness matrix
        */
        for (l = 0; l < nrow * nrow; l++) K_local[l] = 0.;
        //ierr = PetscLogEventBegin(ctx->vflog.VF_MatULocalEvent,0,0,0,0);CHKERRQ(ierr);
        switch (ctx->unilateral) {
          case UNILATERAL_NONE:
            ierr = VF_MatU3D_local(K_local,v_array,&ctx->matprop[ctx->layer[ek]],&ctx->vfprop,
                                   ek,ej,ei,&ctx->e3D);
            break;
          case UNILATERAL_SHEARONLY:
            ierr = VF_MatUShearOnly3D_local(K_local,v_array,&ctx->matprop[ctx->layer[ek]],&ctx->vfprop,
                                            ek,ej,ei,&ctx->e3D);
            break;
        }
        //ierr = PetscLogEventEnd(ctx->vflog.VF_MatULocalEvent,0,0,0,0);CHKERRQ(ierr);
      
        for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
          for (j = 0; j < ctx->e3D.nphiy; j++) {
            for (i = 0; i < ctx->e3D.nphix; i++) {
              for (c = 0; c < dim; c++,l++) {
                row[l].i = ei + i; row[l].j = ej + j; row[l].k = ek + k; row[l].c = c;
              }
            }
          }
        }
        ierr = MatSetValuesStencil(K,nrow,row,nrow,row,K_local,ADD_VALUES);CHKERRQ(ierr);
        /*
          Compute and accumulate the local contribution of the effective strain contribution to the global RHS
        */
        for (l = 0; l < nrow; l++) RHS_local[l] = 0.;
        //ierr = PetscLogEventBegin(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
        switch (ctx->unilateral) {
          case UNILATERAL_NONE:
            ierr = VF_RHSUCoupling3D_local(RHS_local,v_array,theta_array,thetaRef_array,
                                           pressure_array,pressureRef_array,
                                           &ctx->matprop[ctx->layer[ek]],&ctx->vfprop,
                                           ek,ej,ei,&ctx->e3D);CHKERRQ(ierr);
            break;
          case UNILATERAL_SHEARONLY:
            ierr = VF_RHSUCouplingShearOnly3D_local(RHS_local,v_array,theta_array,thetaRef_array,
                                                    pressure_array,pressureRef_array,
                                                    &ctx->matprop[ctx->layer[ek]],&ctx->vfprop,
                                                    ek,ej,ei,&ctx->e3D);CHKERRQ(ierr);
            break;
        }
        if (ctx->hasCrackPressure) {
          ierr = VF_RHSUPressure3D_local(RHS_local,v_array,pressure_array,pressureRef_array,
                                           &ctx->matprop[ctx->layer[ek]],&ctx->vfprop,
                                           ek,ej,ei,&ctx->e3D);CHKERRQ(ierr);
        }
        //ierr = PetscLogEventEnd(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
        for (l=0,k = 0; k < ctx->e3D.nphiz; k++) {
          for (j = 0; j < ctx->e3D.nphiy; j++) {
            for (i = 0; i < ctx->e3D.nphix; i++){
              for (c = 0; c < dim; c++,l++) {
                RHS_array[ek+k][ej+j][ei+i][c] += RHS_local[l];
                RHS_local[l]=0;
              }
            }
          }
        }
        /*
          Compute and accumulate the local contribution of the insitu stresses to the global RHS
        */
        if (ctx->hasInsitu) {
          if (ek == 0 ) {
            /* 
              Face Z0
              sigma.(0,0,-1) = (s_13,s_23,-s_33) = (S4,S3,-S2)
            */
            face = Z0;
            stresscomp[0] = 4; stressdir[0] = 1.;
            stresscomp[1] = 3; stressdir[1] = 1.;
            stresscomp[2] = 2; stressdir[2] = -1.;
            for (c = 0; c < 3; c++) {
              if (ctx->bcU[c].face[face] == NONE) {
                for (k = 0; k < ctx->e3D.nphiz; k++){
                  z = coords_array[ek+k][ej][ei][2];
                  stressmag = stressdir[c] * 
                           (ctx->insitumin[stresscomp[c]] + (z - BBmin[2]) / (BBmax[2] - BBmin[2])
                            * (ctx->insitumax[stresscomp[c]] - ctx->insitumin[stresscomp[c]]));
                  for (j = 0; j < ctx->e3D.nphiy; j++) {
                    for (i = 0; i < ctx->e3D.nphix; i++) {
                      f_array[ek+k][ej+j][ei+i][c] = stressmag;
                    }
                  }
                }
              }
            }
            //ierr = PetscLogEventBegin(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            ierr = VF_RHSUInSituStresses3D_local(RHS_local,f_array,ek,ej,ei,face,&ctx->e3D);CHKERRQ(ierr);
            //ierr = PetscLogEventEnd(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            for (l=0,k = 0; k < ctx->e3D.nphiz; k++){
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++) {
                  for (c = 0; c < dim; c++,l++) {
                    RHS_array[0][ej+j][ei+i][c] += RHS_local[l];
                    RHS_local[l]=0;
                  }
                }
              }
            }
          }
          if (ek == nz-2) {
            /* 
              Face Z1
              sigma.(0,0,1) = (s_13,s_23,s_33) = (S4,S3,S2)
            */
            face = Z1;
            stresscomp[0] = 4; stressdir[0] = 1.;
            stresscomp[1] = 3; stressdir[1] = 1.;
            stresscomp[2] = 2; stressdir[2] = 1.;
            for (k = 0; k < ctx->e3D.nphiz; k++){
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++) {
                    z = coords_array[ek+k][ej+j][ei+i][2]; 
                    for (c = 0; c < 3; c++) {
                       if (ctx->bcU[c].face[face] == NONE) {
                          f_array[ek+k][ej+j][ei+i][c] = stressdir[c] * 
                                                         (ctx->insitumin[stresscomp[c]] + 
                                                          (z - BBmin[2]) / (BBmax[2] - BBmin[2]) * 
                                                          (ctx->insitumax[stresscomp[c]] - ctx->insitumin[stresscomp[c]]));
                       }
                    }
                }
              }
            }
            //ierr = PetscLogEventBegin(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            ierr = VF_RHSUInSituStresses3D_local(RHS_local,f_array,ek,ej,ei,face,&ctx->e3D);CHKERRQ(ierr);
            //ierr = PetscLogEventEnd(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            for (l=0,k = 0; k < ctx->e3D.nphiz; k++){
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++) {
                  for (c = 0; c < dim; c++,l++) {
                    RHS_array[nz-1][ej+j][ei+i][c] += RHS_local[l];
                    RHS_local[l]=0;
                  }
                }
              }
            }
          }
        
          if (ej == 0) {
            /* 
              Face Y0
              sigma.(0,-1,0) = (s_12,-s_22,-s_23) = (S5,-S1,-S3)
              the negative sign in the 3rd component comes from thaty the z-axis is pointing down
            */
            face = Y0;
            stresscomp[0] = 5; stressdir[0] = 1.;
            stresscomp[1] = 1; stressdir[1] = -1.;
            stresscomp[2] = 3; stressdir[2] = -1.;
            for (k = 0; k < ctx->e3D.nphiz; k++){
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++) {
                    z = coords_array[ek+k][ej+j][ei+i][2];
                    for (c = 0; c < 3; c++) {
                       if (ctx->bcU[c].face[face] == NONE) {
                          f_array[ek+k][ej+j][ei+i][c] = stressdir[c] * 
                                                         (ctx->insitumin[stresscomp[c]] + 
                                                          (z - BBmin[2]) / (BBmax[2] - BBmin[2]) * 
                                                          (ctx->insitumax[stresscomp[c]] - ctx->insitumin[stresscomp[c]]));
                       }
                    }
                }
              }
            }
            //ierr = PetscLogEventBegin(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            ierr = VF_RHSUInSituStresses3D_local(RHS_local,f_array,ek,ej,ei,face,&ctx->e3D);CHKERRQ(ierr);
            //ierr = PetscLogEventEnd(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            for (l=0,k = 0; k < ctx->e3D.nphiz; k++){
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++) {
                  for (c = 0; c < dim; c++,l++) {
                    RHS_array[ek+k][0][ei+i][c] += RHS_local[l];
                    RHS_local[l]=0;
                  }
                }
              }
            }
          }
          if (ej == ny-2) {
            /* 
              Face Y1
              sigma.(0,1,0) = (s_12,s_22,-s_23) = (S5,S1,-S3)
              the negative sign in the 3rd component comes from thaty the z-axis is pointing down
            */
            face = Y1;
            stresscomp[0] = 5; stressdir[0] = 1.;
            stresscomp[1] = 1; stressdir[1] = 1.;
            stresscomp[2] = 3; stressdir[2] = -1.;
            for (k = 0; k < ctx->e3D.nphiz; k++){
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++) {
                    z = coords_array[ek+k][ej+j][ei+i][2];
                    for (c = 0; c < 3; c++) {
                       if (ctx->bcU[c].face[face] == NONE) {
                          f_array[ek+k][ej+j][ei+i][c] = stressdir[c] * 
                                                         (ctx->insitumin[stresscomp[c]] + 
                                                          (z - BBmin[2]) / (BBmax[2] - BBmin[2]) * 
                                                          (ctx->insitumax[stresscomp[c]] - ctx->insitumin[stresscomp[c]]));
                       }
                    }
                }
              }
            }
            //ierr = PetscLogEventBegin(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            ierr = VF_RHSUInSituStresses3D_local(RHS_local,f_array,ek,ej,ei,face,&ctx->e3D);CHKERRQ(ierr);
            //ierr = PetscLogEventEnd(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            for (l=0,k = 0; k < ctx->e3D.nphiz; k++){
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++) {
                  for (c = 0; c < dim; c++,l++) {
                    RHS_array[ek+k][ny-1][ei+i][c] += RHS_local[l];
                    RHS_local[l]=0;
                  }
                }
              }
            }
          }
        
          if (ei == 0) {
            /* 
              Face X0
              sigma.(-1,0,0) = (-s_11,s_12,-s_13) = (-S0,S5,-S4)
              the negative sign in the 3rd component comes from thaty the z-axis is pointing down
            */
            face = X0;
            stresscomp[0] = 0; stressdir[0] = -1.;
            stresscomp[1] = 5; stressdir[1] = 1.;
            stresscomp[2] = 4; stressdir[2] = -1.;
            for (k = 0; k < ctx->e3D.nphiz; k++){
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++) {
                    z = coords_array[ek+k][ej+j][ei+i][2];
                    for (c = 0; c < 3; c++) {
                       if (ctx->bcU[c].face[face] == NONE) {
                          f_array[ek+k][ej+j][ei+i][c] = stressdir[c] * 
                                                         (ctx->insitumin[stresscomp[c]] + 
                                                          (z - BBmin[2]) / (BBmax[2] - BBmin[2]) * 
                                                          (ctx->insitumax[stresscomp[c]] - ctx->insitumin[stresscomp[c]]));
                       }
                    }
                }
              }
            }
            //ierr = PetscLogEventBegin(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            ierr = VF_RHSUInSituStresses3D_local(RHS_local,f_array,ek,ej,ei,face,&ctx->e3D);CHKERRQ(ierr);
            //ierr = PetscLogEventEnd(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            for (l=0,k = 0; k < ctx->e3D.nphiz; k++){
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++) {
                  for (c = 0; c < dim; c++,l++) {
                    RHS_array[ek+k][ej+j][0][c] += RHS_local[l];
                    RHS_local[l]=0;
                  }
                }
              }
            }
          }
          if (ei == nx-2) {
            /* 
              Face X1
              sigma.(1,0,0) = (s_11,s_12,-s_13) = (S0,S5,-S4)
              the negative sign in the 3rd component comes from thaty the z-axis is pointing down
            */
            face = X1;
            stresscomp[0] = 0; stressdir[0] = 1.;
            stresscomp[1] = 5; stressdir[1] = 1.;
            stresscomp[2] = 4; stressdir[2] = -1.;
            for (k = 0; k < ctx->e3D.nphiz; k++){
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++) {
                  z = coords_array[ek+k][ej+j][ei+i][2];
                   for (c = 0; c < 3; c++) {
                     if (ctx->bcU[c].face[face] == NONE) {
                       f_array[ek+k][ej+j][ei+i][c] = stressdir[c] * 
                                                      (ctx->insitumin[stresscomp[c]] + 
                                                        (z - BBmin[2]) / (BBmax[2] - BBmin[2]) * 
                                                        (ctx->insitumax[stresscomp[c]] - ctx->insitumin[stresscomp[c]]));
                     }
                   }
                }
              }
            }
            //ierr = PetscLogEventBegin(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            ierr = VF_RHSUInSituStresses3D_local(RHS_local,f_array,ek,ej,ei,face,&ctx->e3D);CHKERRQ(ierr);
            //ierr = PetscLogEventEnd(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            for (l=0,k = 0; k < ctx->e3D.nphiz; k++){
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++) {
                  for (c = 0; c < dim; c++,l++) {
                    RHS_array[ek+k][ej+j][nx-1][c] += RHS_local[l];
                    RHS_local[l]=0;
                  }
                }
              }
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
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  ierr = MatApplyDirichletBC(K,ctx->daVect,&ctx->bcU[0]);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
 
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,RHS_localVec,&RHS_array);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(ctx->daVect,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ctx->daVect,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daVect,&RHS_localVec);CHKERRQ(ierr);
  ierr = VecApplyDirichletBC(RHS,fields->BCU,&ctx->bcU[0]);CHKERRQ(ierr);
  
/*
    Cleanup
  */
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,V_localVec,&v_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&V_localVec);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,theta_localVec,&theta_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&theta_localVec);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,thetaRef_localVec,&thetaRef_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&thetaRef_localVec);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,pressure_localVec,&pressure_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&pressure_localVec);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,pressureRef_localVec,&pressureRef_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&pressureRef_localVec);CHKERRQ(ierr);
  if (ctx->hasInsitu) {
    ierr = DMDAVecRestoreArrayDOF(ctx->daVect,f_localVec,&f_array);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(ctx->daVect,&f_localVec);CHKERRQ(ierr);
  }
  ierr = PetscFree3(RHS_local,K_local,row);CHKERRQ(ierr);
  //ierr = PetscLogStagePop();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_ElasticEnergy3D_local"
/*
  VF_ElasticEnergy3D_local: 

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_ElasticEnergy3D_local(PetscReal *ElasticEnergy_local,
                                               PetscReal ****u_array,PetscReal ***v_array,
                                               PetscReal ***theta_array,PetscReal ***thetaRef_array,
                                               PetscReal ***pressure_array,PetscReal ***pressureRef_array,
                                               MatProp *matprop,VFProp *vfprop,
                                               PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e)
{
  PetscInt       g,i,j,k;
  PetscReal      *epsilon11_elem,*epsilon22_elem,*epsilon33_elem,*epsilon12_elem,*epsilon23_elem,*epsilon13_elem;
  PetscReal      *sigma11_elem,*sigma22_elem,*sigma33_elem,*sigma12_elem,*sigma23_elem,*sigma13_elem;
  PetscReal      lambda,mu,alpha,coefbeta;
  PetscReal      *v_elem;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscMalloc3(e->ng,PetscReal,&sigma11_elem,e->ng,PetscReal,&sigma22_elem,e->ng,PetscReal,&sigma33_elem);CHKERRQ(ierr);
  ierr = PetscMalloc3(e->ng,PetscReal,&sigma12_elem,e->ng,PetscReal,&sigma23_elem,e->ng,PetscReal,&sigma13_elem);CHKERRQ(ierr);
  ierr = PetscMalloc3(e->ng,PetscReal,&epsilon11_elem,e->ng,PetscReal,&epsilon22_elem,e->ng,PetscReal,&epsilon33_elem);CHKERRQ(ierr);
  ierr = PetscMalloc3(e->ng,PetscReal,&epsilon12_elem,e->ng,PetscReal,&epsilon23_elem,e->ng,PetscReal,&epsilon13_elem);CHKERRQ(ierr);
  ierr = PetscMalloc(e->ng * sizeof(PetscReal),&v_elem);CHKERRQ(ierr);
  lambda   = matprop->lambda;
  mu       = matprop->mu;
  alpha    = matprop->alpha;
  coefbeta = matprop->beta / (3.*lambda + 2. *mu);
  
  for (g = 0; g < e->ng; g++) v_elem[g] = 0.;
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (g = 0; g < e->ng; g++) {
          v_elem[g] += v_array[ek+k][ej+j][ei+i] * e->phi[k][j][i][g];
        }
      }
    }
  }
  /*
    epsilon is the inelastic strain
  */
  for (g = 0; g < e->ng; g++) {
    epsilon11_elem[g] = 0;
    epsilon22_elem[g] = 0;
    epsilon33_elem[g] = 0;
    epsilon12_elem[g] = 0;
    epsilon23_elem[g] = 0;
    epsilon13_elem[g] = 0;
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (g = 0; g < e->ng; g++) {
          epsilon11_elem[g] +=  e->dphi[k][j][i][0][g] * u_array[ek+k][ej+j][ei+i][0] 
                              - e->phi[k][j][i][g] * ( 
                                  coefbeta * (pressure_array[ek+k][ej+j][ei+i] - pressureRef_array[ek+k][ej+j][ei+i])
                                + alpha * (theta_array[ek+k][ej+j][ei+i] - thetaRef_array[ek+k][ej+j][ei+i]));
          epsilon22_elem[g] +=  e->dphi[k][j][i][1][g] * u_array[ek+k][ej+j][ei+i][1] 
                              - e->phi[k][j][i][g] * ( 
                                  coefbeta * (pressure_array[ek+k][ej+j][ei+i] - pressureRef_array[ek+k][ej+j][ei+i])
                                + alpha * (theta_array[ek+k][ej+j][ei+i] - thetaRef_array[ek+k][ej+j][ei+i]));
          epsilon33_elem[g] +=  e->dphi[k][j][i][2][g] * u_array[ek+k][ej+j][ei+i][2] 
                              - e->phi[k][j][i][g] * ( 
                                  coefbeta * (pressure_array[ek+k][ej+j][ei+i] - pressureRef_array[ek+k][ej+j][ei+i])
                                + alpha * (theta_array[ek+k][ej+j][ei+i] - thetaRef_array[ek+k][ej+j][ei+i]));

          epsilon12_elem[g] += (e->dphi[k][j][i][0][g] * u_array[ek+k][ej+j][ei+i][1] 
                              + e->dphi[k][j][i][1][g] * u_array[ek+k][ej+j][ei+i][0]) * .5;
          epsilon23_elem[g] += (e->dphi[k][j][i][1][g] * u_array[ek+k][ej+j][ei+i][2] 
                              + e->dphi[k][j][i][2][g] * u_array[ek+k][ej+j][ei+i][1]) * .5;
          epsilon13_elem[g] += (e->dphi[k][j][i][0][g] * u_array[ek+k][ej+j][ei+i][2] 
                              + e->dphi[k][j][i][2][g] * u_array[ek+k][ej+j][ei+i][0]) * .5;
        }
      }
    }  
  }

  *ElasticEnergy_local = 0.;
  for (g = 0; g < e->ng; g++) {
    sigma11_elem[g] = (lambda + 2.*mu) * epsilon11_elem[g] + lambda * epsilon22_elem[g] + lambda * epsilon33_elem[g];
    sigma22_elem[g] = lambda * epsilon11_elem[g] + (lambda + 2.*mu) * epsilon22_elem[g] + lambda * epsilon33_elem[g];
    sigma33_elem[g] = lambda * epsilon11_elem[g] + lambda * epsilon22_elem[g] + (lambda + 2.*mu) * epsilon33_elem[g];
    sigma12_elem[g] = 2. * mu * epsilon12_elem[g];
    sigma23_elem[g] = 2. * mu * epsilon23_elem[g];
    sigma13_elem[g] = 2. * mu * epsilon13_elem[g];
    *ElasticEnergy_local += ((    sigma11_elem[g] * epsilon11_elem[g] 
                                + sigma22_elem[g] * epsilon22_elem[g] 
                                + sigma33_elem[g] * epsilon33_elem[g]) * .5
                              + sigma12_elem[g] * epsilon12_elem[g]
                              + sigma23_elem[g] * epsilon23_elem[g]
                              + sigma13_elem[g] * epsilon13_elem[g] ) * (v_elem[g] * v_elem[g] + vfprop->eta) * e->weight[g];
  }
  
  /*
  ierr = PetscPrintf(PETSC_COMM_SELF,"%s Elast E[%i,%i,%i]: %e\n",__FUNCT__,ei,ej,ek,*ElasticEnergy_local);
  */
  ierr = PetscFree3(sigma11_elem,sigma22_elem,sigma33_elem);CHKERRQ(ierr);
  ierr = PetscFree3(sigma12_elem,sigma23_elem,sigma13_elem);CHKERRQ(ierr);
  ierr = PetscFree3(epsilon11_elem,epsilon22_elem,epsilon33_elem);CHKERRQ(ierr);
  ierr = PetscFree3(epsilon12_elem,epsilon23_elem,epsilon13_elem);CHKERRQ(ierr);
  ierr = PetscFree(v_elem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_PressureWork3D_local"
/*
  VF_PressureWork3D_local: Compute the contribution of an element to the work of the pressure forces along the crack walls,
  given by
    \int_e p(x)\nabla v(x) \cdot u(x) \, dx

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_PressureWork3D_local(PetscReal *PressureWork_local,PetscReal ****u_array,PetscReal ***v_array,
                                              PetscReal ***pressure_array,PetscReal ***pressureRef_array,
                                              MatProp *matprop,VFProp *vfprop,
                                              PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,g,c;
  PetscReal      *pressure_elem,*gradv_elem[3],*u_elem[3];

  PetscFunctionBegin;
  ierr = PetscMalloc4(e->ng,PetscReal,&pressure_elem,
                      e->ng,PetscReal,&gradv_elem[0],
                      e->ng,PetscReal,&gradv_elem[1],
                      e->ng,PetscReal,&gradv_elem[2]);CHKERRQ(ierr);
  ierr = PetscMalloc3(e->ng,PetscReal,&u_elem[0],
                      e->ng,PetscReal,&u_elem[1],
                      e->ng,PetscReal,&u_elem[2]);CHKERRQ(ierr);

  /*
    Compute the projection of the fields in the local base functions basis
  */
  for (g = 0; g < e->ng; g++) {
    pressure_elem[g] = 0;
    for (c=0; c<3; c++) {
      gradv_elem[c][g] = 0.;
      u_elem[c][g] = 0.;
    }
  }
  for (k = 0; k < e->nphiz; k++) {
    for (j = 0; j < e->nphiy; j++) {
      for (i = 0; i < e->nphix; i++) {
        for (g = 0; g < e->ng; g++) {
          pressure_elem[g] += e->phi[k][j][i][g] 
                              * (pressure_array[ek+k][ej+j][ei+i] - pressureRef_array[ek+k][ej+j][ei+i]);
          for (c = 0; c < 3; c++) {
            gradv_elem[c][g] += e->dphi[k][j][i][c][g] * v_array[ek+k][ej+j][ei+i];
            u_elem[c][g]     += e->phi[k][j][i][g]     * u_array[ek+k][ej+j][ei+i][c];
          }
        }
      }
    }
  }  
  ierr = PetscLogFlops(15 * e->ng * e->nphix * e->nphiy * e->nphiz);CHKERRQ(ierr);

  /*
    Accumulate the contribution of the current element
  */
  for (g = 0; g < e->ng; g++) {
    *PressureWork_local += e->weight[g] * pressure_elem[g] 
                         * (u_elem[0][g] * gradv_elem[0][g] + u_elem[1][g] * gradv_elem[1][g] + u_elem[2][g] * gradv_elem[2][g]);  
  }
  ierr = PetscLogFlops(e->ng * e->nphix * e->nphiy * e->nphiz);CHKERRQ(ierr);
  /*
    Clean up
  */
  ierr = PetscFree4(pressure_elem,gradv_elem[0],gradv_elem[1],gradv_elem[2]);CHKERRQ(ierr);
  ierr = PetscFree3(u_elem[0],u_elem[1],u_elem[2]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_InSituStressWork3D_local"
/*
  VF_InSituStressWork3D_local

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_InSituStressWork3D_local(PetscReal *Work_local,PetscReal ****u_array,PetscReal ****f_array,PetscInt ek,PetscInt ej,PetscInt ei,FACE face,CartFE_Element3D *e)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,c,g;
  PetscInt       dim=3;
  PetscReal      *u_elem[3],*f_elem[3];

  PetscFunctionBegin;
  /*
    Initialize
  */
  ierr = PetscMalloc3(e->ng,PetscReal,&u_elem[0],e->ng,PetscReal,&u_elem[1],e->ng,PetscReal,&u_elem[2]);CHKERRQ(ierr);
  ierr = PetscMalloc3(e->ng,PetscReal,&f_elem[0],e->ng,PetscReal,&f_elem[1],e->ng,PetscReal,&f_elem[2]);CHKERRQ(ierr);
  for (c = 0; c < dim; c++) {
    for (g = 0; g < e->ng; g++) {
      u_elem[c][g] = 0;
      f_elem[c][g] = 0;
    }
  }
  /*
    Compute
  */
  switch(face) {
    case X0:
      for (k = 0; k < e->nphiz; k++) {
        for (j = 0; j < e->nphiy; j++) {
          for (i = 0; i < e->nphix; i++){
            for (c = 0; c < dim; c++) {
              for (g = 0; g < e->ng; g++) {
                u_elem[c][g] += e->phi[k][j][i][g] * u_array[ek+k][ej+j][ei][c] / e->lx;
                f_elem[c][g] += e->phi[k][j][i][g] * f_array[ek+k][ej+j][ei][c];
              }
            }
          }
        }
      }
      break;
    case X1:
      for (k = 0; k < e->nphiz; k++) {
        for (j = 0; j < e->nphiy; j++) {
          for (i = 0; i < e->nphix; i++){
            for (c = 0; c < dim; c++) {
              for (g = 0; g < e->ng; g++) {
                u_elem[c][g] += e->phi[k][j][i][g] * u_array[ek+k][ej+j][ei+e->nphix-1][c] / e->lx;
                f_elem[c][g] += e->phi[k][j][i][g] * f_array[ek+k][ej+j][ei+e->nphix-1][c];
              }
            }
          }
        }
      }
      break;
    case Y0:
      for (k = 0; k < e->nphiz; k++) {
        for (j = 0; j < e->nphiy; j++) {
          for (i = 0; i < e->nphix; i++){
            for (c = 0; c < dim; c++) {
              for (g = 0; g < e->ng; g++) {
                u_elem[c][g] += e->phi[k][j][i][g] * u_array[ek+k][ej][ei+i][c] / e->ly;
                f_elem[c][g] += e->phi[k][j][i][g] * f_array[ek+k][ej][ei+i][c];
              }
            }
          }
        }
      }
      break;
    case Y1:
      for (k = 0; k < e->nphiz; k++) {
        for (j = 0; j < e->nphiy; j++) {
          for (i = 0; i < e->nphix; i++){
            for (c = 0; c < dim; c++) {
              for (g = 0; g < e->ng; g++) {
                u_elem[c][g] += e->phi[k][j][i][g] * u_array[ek+k][ej+e->nphiy-1][ei+i][c] / e->ly;
                f_elem[c][g] += e->phi[k][j][i][g] * f_array[ek+k][ej+e->nphiy-1][ei+i][c];
              }
            }
          }
        }
      }
      break;
    case Z0:
      for (k = 0; k < e->nphiz; k++) {
        for (j = 0; j < e->nphiy; j++) {
          for (i = 0; i < e->nphix; i++){
            for (c = 0; c < dim; c++) {
              for (g = 0; g < e->ng; g++) {
                u_elem[c][g] += e->phi[k][j][i][g] * u_array[ek][ej+j][ei+i][c] / e->lz;
                f_elem[c][g] += e->phi[k][j][i][g] * f_array[ek][ej+j][ei+i][c];
              }
            }
          }
        }
      }
      break;
    case Z1:
      for (k = 0; k < e->nphiz; k++) {
        for (j = 0; j < e->nphiy; j++) {
          for (i = 0; i < e->nphix; i++){
            for (c = 0; c < dim; c++) {
              for (g = 0; g < e->ng; g++) {
                u_elem[c][g] += e->phi[k][j][i][g] * u_array[ek+e->nphiz-1][ej+j][ei+i][c] / e->lz;
                f_elem[c][g] += e->phi[k][j][i][g] * f_array[ek+e->nphiz-1][ej+j][ei+i][c];
              }
            }
          }
        }
      }
      break;
  }
  /*
    Accumulate
  */
  for (c = 0; c < dim; c++) {
    for (g = 0; g < e->ng; g++) {
      *Work_local += e->weight[g] * u_elem[c][g] * f_elem[c][g];
    }
  }
  ierr = PetscFree3(u_elem[0],u_elem[1],u_elem[2]);CHKERRQ(ierr);
  ierr = PetscFree3(f_elem[0],f_elem[1],f_elem[2]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_UEnergy3D"
/*
  VF_UEnergy3D

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_UEnergy3D(PetscReal *ElasticEnergy,PetscReal *InsituWork,PetscReal *PressureWork,VFFields *fields,VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscInt        xs,xm,nx;
  PetscInt        ys,ym,ny;
  PetscInt        zs,zm,nz;
  PetscInt        ei,ej,ek;
  PetscInt        i,j,k,c;
  Vec             u_localVec,v_localVec;
  Vec             theta_localVec,thetaRef_localVec;
  Vec             pressure_localVec,pressureRef_localVec;
  PetscReal   ****u_array;
  PetscReal    ***v_array;
  PetscReal    ***theta_array;
  PetscReal    ***thetaRef_array;
  PetscReal    ***pressure_array;
  PetscReal    ***pressureRef_array;
  PetscReal       myInsituWork=0.,myPressureWork=0.;
  PetscReal       myElasticEnergy=0.,myElasticEnergyLocal=0.;
  PetscReal       hx,hy,hz;
  PetscReal   ****coords_array;
  PetscReal   ****f_array;
  Vec             f_localVec;
  FACE            face;
  PetscReal       z;
  int             stresscomp[3];
  PetscReal       stressdir[3];
  PetscReal       stressmag;
  PetscReal       BBmin[3],BBmax[3];

  PetscFunctionBegin;
  //ierr = PetscLogStagePush(ctx->vflog.VF_EnergyStage);CHKERRQ(ierr);
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  /*
    Get coordinates
  */
  ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  /*
    Get bounding box from petsc DA
  */
  ierr = DMDAGetBoundingBox(ctx->daVect,BBmin,BBmax);CHKERRQ(ierr);
  /*
    get U_array
  */
  ierr = DMGetLocalVector(ctx->daVect,&u_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daVect,fields->U,INSERT_VALUES,u_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daVect,fields->U,INSERT_VALUES,u_localVec);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(ctx->daVect,u_localVec,&u_array);CHKERRQ(ierr);
  /*
    get v_array
  */
  ierr = DMGetLocalVector(ctx->daScal,&v_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->daScal,fields->V,INSERT_VALUES,v_localVec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->daScal,fields->V,INSERT_VALUES,v_localVec);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,v_localVec,&v_array);CHKERRQ(ierr);    
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
    Allocating multi-dimensional vectors in C is a pain, so for the in-situ stresses / surface stresses, 
    I get a 3 component Vec for the external forces, 
    then get a full array. This is a bit wastefull since we never use the internal values...
    Note that we cannot fully initialize it since we need the normal direction to each face.
  */
  if (ctx->hasInsitu) {
    ierr = DMGetLocalVector(ctx->daVect,&f_localVec);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(ctx->daVect,f_localVec,&f_array);CHKERRQ(ierr); 
  }
  
  *ElasticEnergy = 0;
  *InsituWork = 0;
  *PressureWork = 0.;
  for (ek = zs; ek < zs + zm; ek++) {
    for (ej = ys; ej < ys+ym; ej++) {
      for (ei = xs; ei < xs+xm; ei++) {
        hx = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
        hy = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
        hz = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
        ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
        /*
          Elastic Energy is trivial
        */
        //ierr = PetscLogEventBegin(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
        ierr = VF_ElasticEnergy3D_local(&myElasticEnergyLocal,u_array,v_array,
                                        theta_array,thetaRef_array,
                                        pressure_array,pressureRef_array,
                                        &ctx->matprop[ctx->layer[ek]],&ctx->vfprop,
                                        ek,ej,ei,&ctx->e3D);CHKERRQ(ierr); 
        myElasticEnergy += myElasticEnergyLocal;
        if (ctx->hasCrackPressure) {
          ierr = VF_PressureWork3D_local(&myPressureWork,u_array,v_array,pressure_array,pressureRef_array,
                                           &ctx->matprop[ctx->layer[ek]],&ctx->vfprop,
                                           ek,ej,ei,&ctx->e3D);CHKERRQ(ierr);
        }
        //ierr = PetscLogEventEnd(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
        
        if (ctx->hasInsitu) {
          /*
            We need to reconstruct the external forces before computing their work. 
            This take a bit more effort
          */
          if (ek == 0) {
            /* 
              Face Z0
              sigma.(0,0,1) = (s_13,s_23,s_33) = (S4,S3,S2)
            */
            face = Z0;
            stresscomp[0] = 4; stressdir[0] = 1.;
            stresscomp[1] = 3; stressdir[1] = 1.;
            stresscomp[2] = 2; stressdir[2] = 1.;
            for (c = 0; c < 3; c++) {
              if (ctx->bcU[c].face[face] == NONE) {
                for (k = 0; k < ctx->e3D.nphiz; k++){
                  z = coords_array[ek+k][ej][ei][2];
                  stressmag = stressdir[c] * 
                           (ctx->insitumin[stresscomp[c]] + (z - BBmin[2]) / (BBmax[2] - BBmin[2])
                            * (ctx->insitumax[stresscomp[c]] - ctx->insitumin[stresscomp[c]]));
                  for (j = 0; j < ctx->e3D.nphiy; j++) {
                    for (i = 0; i < ctx->e3D.nphix; i++) {
                      f_array[ek+k][ej+j][ei+i][c] = stressmag;
                    }
                  }
                }
              }
            }
            //ierr = PetscLogEventBegin(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            ierr = VF_InSituStressWork3D_local(&myInsituWork,u_array,f_array,ek,ej,ei,face,&ctx->e3D);CHKERRQ(ierr); 
            //ierr = PetscLogEventEnd(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
          }
        
          if (ek == nz-2) {
            /* 
              Face Z1
              sigma.(0,0,-1) = (s_13,s_23,-s_33) = (S4,S3,-S2)
            */
            face = Z1;
            stresscomp[0] = 4; stressdir[0] = 1.;
            stresscomp[1] = 3; stressdir[1] = 1.;
            stresscomp[2] = 2; stressdir[2] = -1.;
            for (c = 0; c < 3; c++) {
              if (ctx->bcU[c].face[face] == NONE) {
                for (k = 0; k < ctx->e3D.nphiz; k++){
                  z = coords_array[ek+k][ej][ei][2];
                  stressmag = stressdir[c] * 
                             (ctx->insitumin[stresscomp[c]] + (z - BBmin[2]) / (BBmax[2] - BBmin[2])
                              * (ctx->insitumax[stresscomp[c]] - ctx->insitumin[stresscomp[c]]));
                  for (j = 0; j < ctx->e3D.nphiy; j++) {
                    for (i = 0; i < ctx->e3D.nphix; i++) {
                       f_array[ek+k][ej+j][ei+i][c] = stressmag;
                     }
                   }
                }
              }
            }
            //ierr = PetscLogEventBegin(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            ierr = VF_InSituStressWork3D_local(&myInsituWork,u_array,f_array,ek,ej,ei,face,&ctx->e3D);CHKERRQ(ierr); 
            //ierr = PetscLogEventEnd(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
          }
        
          if (ej == 0) {
            /* 
              Face Y0
              sigma.(0,1,0) = (s_12,s_22,s_23) = (S5,S1,S3)
            */
            face = Y0;
            stresscomp[0] = 5; stressdir[0] = 1.;
            stresscomp[1] = 1; stressdir[1] = 1.;
            stresscomp[2] = 3; stressdir[2] = 1.;
            for (k = 0; k < ctx->e3D.nphiz; k++){
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++) {
                   z = coords_array[ek+k][ej+j][ei+i][2];
                   for (c = 0; c < 3; c++) {
                     if (ctx->bcU[0].face[face] == NONE) {
                       f_array[ek+k][ej+j][ei+i][c] = stressdir[c] * 
                                                      (ctx->insitumin[stresscomp[c]] + 
                                                        (z - BBmin[2]) / (BBmax[2] - BBmin[2]) 
                                                        * (ctx->insitumax[stresscomp[c]] - ctx->insitumin[stresscomp[c]]));
                     }
                   }
                }
              }
            }
            //ierr = PetscLogEventBegin(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            ierr = VF_InSituStressWork3D_local(&myInsituWork,u_array,f_array,ek,ej,ei,face,&ctx->e3D);CHKERRQ(ierr); 
            //ierr = PetscLogEventEnd(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
          }
          
          if (ej == ny-2) {
            /* 
              Face Y1
              sigma.(0,-1,0) = (s_12,-s_22,s_23) = (S5,-S1,S3)
            */
            face = Y1;
            stresscomp[0] = 5; stressdir[0] =  1.;
            stresscomp[1] = 1; stressdir[1] = -1.;
            stresscomp[2] = 3; stressdir[2] =  1.;
            for (k = 0; k < ctx->e3D.nphiz; k++){
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++) {
                   z = coords_array[ek+k][ej+j][ei+i][2];
                   for (c = 0; c < 3; c++) {
                     if (ctx->bcU[0].face[face] == NONE) {
                       f_array[ek+k][ej+j][ei+i][c] = stressdir[c] * 
                                                      (ctx->insitumin[stresscomp[c]] + 
                                                        (z - BBmin[2]) / (BBmax[2] - BBmin[2]) 
                                                        * (ctx->insitumax[stresscomp[c]] - ctx->insitumin[stresscomp[c]]));
                     }
                   }
                }
              }
            }
            //ierr = PetscLogEventBegin(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            ierr = VF_InSituStressWork3D_local(&myInsituWork,u_array,f_array,ek,ej,ei,face,&ctx->e3D);CHKERRQ(ierr); 
            //ierr = PetscLogEventEnd(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
          }
  
          if (ei == 0) {
            /* 
              Face X0
              sigma.(1,0,0) = (s_11,s_12,s_13) = (S0,S5,S4)
            */
            face = X0;
            stresscomp[0] = 0; stressdir[0] = 1.;
            stresscomp[1] = 5; stressdir[1] = 1.;
            stresscomp[2] = 4; stressdir[2] = 1.;
            for (k = 0; k < ctx->e3D.nphiz; k++){
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++) {
                   z = coords_array[ek+k][ej+j][ei+i][2];
                   for (c = 0; c < 3; c++) {
                     if (ctx->bcU[0].face[face] == NONE) {
                       f_array[ek+k][ej+j][ei+i][c] = stressdir[c] * 
                                                      (ctx->insitumin[stresscomp[c]] + 
                                                        (z - BBmin[2]) / (BBmax[2] - BBmin[2]) 
                                                        * (ctx->insitumax[stresscomp[c]] - ctx->insitumin[stresscomp[c]]));
                     }
                   }
                }
              }
            }
            //ierr = PetscLogEventBegin(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            ierr = VF_InSituStressWork3D_local(&myInsituWork,u_array,f_array,ek,ej,ei,face,&ctx->e3D);CHKERRQ(ierr); 
            //ierr = PetscLogEventEnd(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
          }
  
          if (ei == nx-2) {
            /* 
              Face X1
              sigma.(-1,0,0) = (-s_11,s_12,s_13) = (-S0,S5,S4)
            */
            face = X1;
            stresscomp[0] = 0; stressdir[0] = -1.;
            stresscomp[1] = 5; stressdir[1] =  1.;
            stresscomp[2] = 4; stressdir[2] =  1.;
            for (k = 0; k < ctx->e3D.nphiz; k++){
              for (j = 0; j < ctx->e3D.nphiy; j++) {
                for (i = 0; i < ctx->e3D.nphix; i++) {
                   z = coords_array[ek+k][ej+j][ei+i][2];
                   for (c = 0; c < 3; c++) {
                     if (ctx->bcU[0].face[face] != FIXED) {
                       f_array[ek+k][ej+j][ei+i][c] = stressdir[c] * 
                                                      (ctx->insitumin[stresscomp[c]] + 
                                                        (z - BBmin[2]) / (BBmax[2] - BBmin[2]) 
                                                        * (ctx->insitumax[stresscomp[c]] - ctx->insitumin[stresscomp[c]]));
                     }
                   }
                }
              }
            }
            //ierr = PetscLogEventBegin(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
            ierr = VF_InSituStressWork3D_local(&myInsituWork,u_array,f_array,ek,ej,ei,face,&ctx->e3D);CHKERRQ(ierr); 
            //ierr = PetscLogEventEnd(ctx->vflog.VF_VecULocalEvent,0,0,0,0);CHKERRQ(ierr);
          }
        }
      }
    }
  }
  ierr = MPI_Reduce(&myElasticEnergy,ElasticEnergy,1,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Reduce(&myPressureWork,PressureWork,1,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Reduce(&myInsituWork,InsituWork,1,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD);CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,u_localVec,&u_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daVect,&u_localVec);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,v_localVec,&v_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&v_localVec);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(ctx->daScal,theta_localVec,&theta_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&theta_localVec);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,thetaRef_localVec,&thetaRef_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&thetaRef_localVec);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,pressure_localVec,&pressure_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&pressure_localVec);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(ctx->daScal,pressureRef_localVec,&pressureRef_array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->daScal,&pressureRef_localVec);CHKERRQ(ierr);
  if (ctx->hasInsitu) {
    ierr = DMDAVecRestoreArrayDOF(ctx->daVect,f_localVec,&f_array);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(ctx->daVect,&f_localVec);CHKERRQ(ierr);
  }

  //ierr = PetscLogStagePop();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_StepU"
/*
  VF_StepU

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_StepU(VFFields *fields,VFCtx *ctx)
{
  PetscErrorCode     ierr;
//  KSPConvergedReason reason;
	SNESConvergedReason reason;
  PetscInt           its;
  PetscReal          Umin,Umax;
  
  PetscFunctionBegin;
  ierr = VF_UAssembly3D(ctx->KU,ctx->RHSU,fields,ctx);CHKERRQ(ierr);
  if (ctx->hasInsitu) {
    ierr = VecApplyDirichletBC(fields->U,fields->BCU,&ctx->bcU[0]);CHKERRQ(ierr);
  }

  if (ctx->verbose > 1) {
    ierr = MatView(ctx->KU,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);
    ierr = VecView(ctx->RHSU,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  
  //ierr = PetscLogStagePush(ctx->vflog.VF_USolverStage);CHKERRQ(ierr);
  /* 
    REMOVE THIS
  */
  //PetscBool flg = PETSC_FALSE;
  //ierr = MatNullSpaceTest(ctx->nullspaceU,ctx->KU,&flg);CHKERRQ(ierr);
  //ierr = MatNullSpaceRemove(ctx->nullspaceU,ctx->RHSU,PETSC_NULL);CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n\n Null space test returned %i %i\n\n\n",flg,PETSC_TRUE);CHKERRQ(ierr);
	ierr = SNESSetFunction(ctx->snesU,ctx->UFunct,VF_UFunction,ctx);CHKERRQ(ierr);
    ierr = SNESSetJacobian(ctx->snesU,ctx->JacU,ctx->JacU,VF_UIJacobian,ctx);CHKERRQ(ierr);
	if (ctx->verbose > 1) {
		ierr = SNESMonitorSet(ctx->snesVelP,VF_USNESMonitor,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	}
    ierr = SNESSolve(ctx->snesU,PETSC_NULL,fields->U);CHKERRQ(ierr);
//ierr = PetscLogStagePop();CHKERRQ(ierr);
  if (ctx->verbose > 1) {
    ierr = VecView(fields->U,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  if (ctx->verbose > 0) {
    ierr = VecMin(fields->U,PETSC_NULL,&Umin);CHKERRQ(ierr);
    ierr = VecMax(fields->U,PETSC_NULL,&Umax);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Umin = %g, Umax = %g\n",Umin,Umax);CHKERRQ(ierr);
  }
  
  ierr = SNESGetConvergedReason(ctx->snesU,&reason);CHKERRQ(ierr);
  if (reason < 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"[ERROR] snesU diverged with reason %d\n",(int)reason);CHKERRQ(ierr);
  } else {
    ierr = SNESGetIterationNumber(ctx->snesU,&its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"      snesU converged in %d iterations %d.\n",(int)its,(int)reason);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_ComputeBCU"
/*
  VF_ComputeBCU: Compute the boundary displacement by performing an elastic solve with insitu stresses
                 The sets ctx->hasInsitu to PETSC_FALSE and update the BC flags

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VF_ComputeBCU(VFFields *fields,VFCtx *ctx)
{
  PetscErrorCode     ierr;
  SNESConvergedReason reason;
  PetscInt           its;
  PetscReal          Umin,Umax;
  
  PetscFunctionBegin;
  ierr = VF_UAssembly3D(ctx->KU,ctx->RHSU,fields,ctx);CHKERRQ(ierr);

  if (ctx->verbose > 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"%s: computing boundary displacements\n",__FUNCT__);CHKERRQ(ierr);
  }
  if (ctx->verbose > 1) {
    ierr = MatView(ctx->KU,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);
    ierr = VecView(ctx->RHSU,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  
  //ierr = PetscLogStagePush(ctx->vflog.VF_USolverStage);CHKERRQ(ierr);
	ierr = SNESSetFunction(ctx->snesU,ctx->UFunct,VF_UFunction,ctx);CHKERRQ(ierr);
    ierr = SNESSetJacobian(ctx->snesU,ctx->JacU,ctx->JacU,VF_UIJacobian,ctx);CHKERRQ(ierr);
	if (ctx->verbose > 1) {
		ierr = SNESMonitorSet(ctx->snesVelP,VF_USNESMonitor,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	}
    ierr = SNESSolve(ctx->snesU,PETSC_NULL,fields->BCU);CHKERRQ(ierr);
  //ierr = PetscLogStagePop();CHKERRQ(ierr);

  if (ctx->verbose > 1) {
    ierr = VecView(fields->BCU,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  if (ctx->verbose > 0) {
    ierr = VecMin(fields->BCU,PETSC_NULL,&Umin);CHKERRQ(ierr);
    ierr = VecMax(fields->BCU,PETSC_NULL,&Umax);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Umin = %g, Umax = %g\n",Umin,Umax);CHKERRQ(ierr);
  }
  
  ierr = SNESGetConvergedReason(ctx->snesU,&reason);CHKERRQ(ierr);
  if (reason < 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"[ERROR] snesU diverged with reason %d\n",(int)reason);CHKERRQ(ierr);
  } else {
    ierr = SNESGetIterationNumber(ctx->snesU,&its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"      snesU converged in %d iterations %d.\n",(int)its,(int)reason);CHKERRQ(ierr);
  }

  /*
    Copy elastic solution to U, since it is going to be the solution of the first time step
  */
  ierr = VecCopy(fields->BCU,fields->U);CHKERRQ(ierr);

  /* 
    Update boundary condition flags and others
  */
  ctx->hasInsitu = PETSC_FALSE;
  ierr = BCUUpdate(&ctx->bcU[0],ctx->preset);CHKERRQ(ierr);
  if (ctx->verbose > 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"BCU: after update\n");CHKERRQ(ierr);
    ierr = BCView(&ctx->bcU[0],PETSC_VIEWER_STDOUT_WORLD,3);CHKERRQ(ierr);
  }


  ierr = BCVUpdate(&ctx->bcV[0],ctx->preset);CHKERRQ(ierr);
  if (ctx->verbose > 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"BCV: after update\n");CHKERRQ(ierr);
    ierr = BCView(&ctx->bcV[0],PETSC_VIEWER_STDOUT_WORLD,1);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_UFunction"
extern PetscErrorCode VF_UFunction(SNES snes,Vec U,Vec Func,void *user)
{
	PetscErrorCode ierr;
	VFCtx			*ctx=(VFCtx*)user;
	
	PetscFunctionBegin;
	ierr = VecSet(Func,0.0);CHKERRQ(ierr);
	ierr = MatMult(ctx->KU,U,Func);CHKERRQ(ierr);	
	ierr = VecAXPY(Func,-1.0,ctx->RHSU);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_UIJacobian"
extern PetscErrorCode VF_UIJacobian(SNES snes,Vec U,Mat *Jac,Mat *Jac1,MatStructure *str,void *user)
{
	PetscErrorCode ierr;
	VFCtx				*ctx=(VFCtx*)user;
	
	PetscFunctionBegin;
	*str = DIFFERENT_NONZERO_PATTERN;
	ierr = MatZeroEntries(*Jac);CHKERRQ(ierr);
	ierr = MatCopy(ctx->KU,*Jac,*str);
	ierr = MatAssemblyBegin(*Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	if (*Jac != *Jac1) {
		ierr = MatAssemblyBegin(*Jac1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(*Jac1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "VF_USNESMonitor"
extern PetscErrorCode VF_USNESMonitor(SNES snes,PetscInt U_its,PetscReal fnorm,void* ptr)
{
	PetscErrorCode ierr;
	PetscReal      norm,vmax,vmin;
	MPI_Comm       comm;
	Vec				U;
	
	PetscFunctionBegin;
	ierr = SNESGetSolution(snes,&U);CHKERRQ(ierr);	
	ierr = VecNorm(U,NORM_1,&norm);CHKERRQ(ierr);
	ierr = VecMax(U,PETSC_NULL,&vmax);CHKERRQ(ierr);
	ierr = VecMin(U,PETSC_NULL,&vmin);CHKERRQ(ierr);
	ierr = PetscObjectGetComm((PetscObject)snes,&comm);CHKERRQ(ierr);
	ierr = PetscPrintf(comm,"U_snes_iter_step %D :solution norm = %G, max sol. value  = %G, min sol. value = %G\n",U_its,norm,vmax,vmin);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}







 