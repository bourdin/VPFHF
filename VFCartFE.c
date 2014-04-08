#define VFCartFEC
#include "petsc.h"
#include "VFCartFE.h"

static const char *BCTYPE_NAME[] = {
  "NONE",
  "ZERO",
  "ONE",
  "FIXED",
  "BCTYPE_NAME","",0
};
/*
 Type of boundary condition:
 NONE  = Natural boundary condition
 ZERO  = homogeneous Dirichlet
 ONE   = inhomogenous Dirichlet fixed to one
 FIXED = insert from external field (inhomogeneous Dirichlet BC)
*/

static const char *FACE_NAME[] = {
	"X0","X1",
	"Y0","Y1",
	"Z0","Z1",
	"FACE_NAME","",0
};

/*
static const char *EDGE_NAME[] = {
	"X0Z0","X1Z0","Y0Z0","Y1Z0",
	"X0Z1","X1Z1","Y0Z1","Y1Z1",
	"X0Y0","X0Y1","X1Y0","X1Y1",
	"EDGE_NAME","",0
};
*/
static const char *VERTEX_NAME[] = {
	"X0Y0Z0","X1Y0Z0","X0Y1Z0","X1Y1Z0",
	"X0Y0Z1","X1Y0Z1","X0Y1Z1","X1Y1Z1",
	"VERTEX_NAME","",0};


#undef __FUNCT__
#define __FUNCT__ "VFCartFEInit"
/*
  VFCartFEInit

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFCartFEInit()
{
  /*
  PetscErrorCode ierr;
  */
  
  PetscFunctionBegin;
  /*
  ierr = PetscClassIdRegister("VFCartFEElement",&VFCartFEElement);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("VFCartFEElem Init",VFCartFEElement,&VFCartFEElementInitEvent);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("VFCartFEElem Init",&VFCartFEElementInitStage);CHKERRQ(ierr);
  */
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VFCartFEElement1DCreate"
/*
  VFCartFEElement1DCreate: Current element structures are static, so it does not really allocates, but instead 
  initializes the sizes

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFCartFEElement1DCreate(VFCartFEElement1D *e)
{
  PetscFunctionBegin;
  e->dim       = 1;
  e->ng        = 3;
  e->nphix     = 2;
  e->nphiy     = 1;
  e->nphiz     = 1;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFCartFEElement1DInit"
/*
  VFCartFEElement1DInit

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFCartFEElement1DInit(VFCartFEElement1D *e,PetscReal lx)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  e->lx        = lx;
  e->ly        = 0.;
  e->lz        = 0.;
  
  e->weight[0] = 5.*e->lx / 18.;
  e->weight[1] = 8.*e->lx / 18.;
  e->weight[2] = 5.*e->lx / 18.;
  /* 
    Value of the basis function at the integration points 
  */
  e->phi[0][0][0][0] = (1. + sqrt(3./5.)) * .5; e->phi[0][0][0][1] = .5; e->phi[0][0][0][2] = (1. - sqrt(3./5.)) * .5; 
  e->phi[0][0][1][0] = (1. - sqrt(3./5.)) * .5; e->phi[0][0][1][1] = .5; e->phi[0][0][1][2] = (1. + sqrt(3./5.)) * .5; 
  
  /* 
    Value of the derivative of the basis functions at the integration points
  */
  e->dphi[0][0][0][0] = -1. / e->lx; e->dphi[0][0][0][1] = -1. / e->lx; e->dphi[0][0][0][2] = -1. / e->lx; 
  e->dphi[0][0][1][0] =  1. / e->lx; e->dphi[0][0][1][1] =  1. / e->lx; e->dphi[0][0][1][2] =  1. / e->lx;
  ierr = PetscLogFlops(24);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFCartFEElement2DCreate"
/*
  VFCartFEElement2DCreate: Current element structures are static, so it does not really allocates, but instead 
  initializes the sizes

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFCartFEElement2DCreate(VFCartFEElement2D *e)
{
  PetscFunctionBegin;
  e->dim       = 2;
  e->ng        = 9;
  e->nphix     = 2;
  e->nphiy     = 2;
  e->nphiz     = 1;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFCartFEElement2DInit"
/*
  VFCartFEElement2DInit

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFCartFEElement2DInit(VFCartFEElement2D *e,PetscReal lx,PetscReal ly)
{
  PetscErrorCode     ierr;
  VFCartFEElement1D   ex,ey;
  PetscInt           gi,gj,g;   /* local and lexicographic indices of the integration points */
  PetscInt           i,j;
  
  PetscFunctionBegin;
  ierr = VFCartFEElement1DCreate(&ex);CHKERRQ(ierr);
  ierr = VFCartFEElement1DInit(&ex,lx);CHKERRQ(ierr);
  ierr = VFCartFEElement1DCreate(&ey);CHKERRQ(ierr);
  ierr = VFCartFEElement1DInit(&ey,ly);CHKERRQ(ierr);
  
  ierr = VFCartFEElement2DCreate(e);CHKERRQ(ierr);
  e->lx    = lx; 
  e->ly    = ly;
  e->lz    = 0;
  
  for (g = 0,gj = 0; gj < ey.ng; gj++) {
    for (gi = 0; gi < ex.ng; gi++,g++) {
      e->weight[g] = ey.weight[gj] * ex.weight[gi];
    }
  }
  for (j = 0; j < e->nphiy; j++) {
    for (i = 0; i < e->nphix; i++) {
      for (g = 0,gj = 0; gj < ey.ng; gj++) {
        for (gi = 0; gi < ex.ng; gi++,g++) {
          /* 
            Value of the basis functions and their derivatives at the integration points 
            Using the conventions
            phi[j][i](x,y) = phi[j](x) . phi[i](y)
            \partial phi[j][i] / \partial x = phi[j]'(x) . phi[i](y)
            \partial phi[j][i] / \partial y = phi[j](x) . phi[i]'(y)
          */
          
          e->phi[0][j][i][g]     = ey.phi[0][0][j][gj] * ex.phi[0][0][i][gi];
          e->dphi[0][j][i][0][g] = ey.phi[0][0][j][gj] * ex.dphi[0][0][i][gi];          
          e->dphi[0][j][i][1][g] = ey.dphi[0][0][j][gj] * ex.phi[0][0][i][gi];          
        }
      }
    }
  }
  ierr = PetscLogFlops(1 + e->ng + 3 * e->ng * e->nphiy * e->nphix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFCartFEElement3DCreate"
/*
  VFCartFEElement3DCreate: Current element structures are static, so it does not really allocates, but instead 
  initializes the sizes

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFCartFEElement3DCreate(VFCartFEElement3D *e)
{
  PetscFunctionBegin;
  e->dim       = 3;
  e->ng        = 27;
  e->nphix     = 2;
  e->nphiy     = 2;
  e->nphiz     = 2;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFCartFEElement3DInit"
/*
  VFCartFEElement3DInit

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFCartFEElement3DInit(VFCartFEElement3D *e,PetscReal lx,PetscReal ly,PetscReal lz)
{
  PetscErrorCode     ierr;
  VFCartFEElement1D   ex,ey,ez;
  PetscInt           gi,gj,gk,g;   /* local and lexicographic indices of the integration points */
  PetscInt           i,j,k;
  
  PetscFunctionBegin;
  ierr = VFCartFEElement1DCreate(&ex);CHKERRQ(ierr);
  ierr = VFCartFEElement1DInit(&ex,lx);CHKERRQ(ierr);
  ierr = VFCartFEElement1DCreate(&ey);CHKERRQ(ierr);
  ierr = VFCartFEElement1DInit(&ey,ly);CHKERRQ(ierr);
  ierr = VFCartFEElement1DCreate(&ez);CHKERRQ(ierr);
  ierr = VFCartFEElement1DInit(&ez,lz);CHKERRQ(ierr);
  
  ierr = VFCartFEElement3DCreate(e);CHKERRQ(ierr);
  e->lx    = lx; 
  e->ly    = ly;
  e->lz    = lz;
  
  for (g = 0,gk = 0; gk < ez.ng; gk++) {
    for (gj = 0; gj < ey.ng; gj++) {
      for (gi = 0; gi < ex.ng; gi++,g++) {
        e->weight[g] = ez.weight[gk] * ey.weight[gj] * ex.weight[gi];
      }
    }
  }
  for (g = 0,gk = 0; gk < ez.ng; gk++) {
    for (gj = 0; gj < ey.ng; gj++) {
      for (gi = 0; gi < ex.ng; gi++,g++) {
        for (k = 0; k < e->nphiz; k++) {
          for (j = 0; j < e->nphiy; j++) {
            for (i = 0; i < e->nphix; i++) { 
              /* 
                Value of the basis functions and their derivatives at the integration points 
                Using the conventions 
                phi[i][j][k](x,y,z) = phi[i](x) . phi[j](y) . phi[k](z)
                \partial phi[k][j][i] / \partial x = phi[i]'(x) . phi[j](y) . phi[k](z) 
                \partial phi[k][j][i] / \partial y = phi[i](x) . phi[j]'(y) . phi[k](z) 
                \partial phi[k][j][i] / \partial z = phi[i](x) . phi[j](y) . phi[k]'(z) 
              */
              e->phi[k][j][i][g]     = ez.phi[0][0][k][gk] * ey.phi[0][0][j][gj] * ex.phi[0][0][i][gi];
              e->dphi[k][j][i][0][g] = ez.phi[0][0][k][gk] * ey.phi[0][0][j][gj] * ex.dphi[0][0][i][gi];          
              e->dphi[k][j][i][1][g] = ez.phi[0][0][k][gk] * ey.dphi[0][0][j][gj] * ex.phi[0][0][i][gi];          
              e->dphi[k][j][i][2][g] = ez.dphi[0][0][k][gk] * ey.phi[0][0][j][gj] * ex.phi[0][0][i][gi];          
            }
          }
        }
      }
    }
  }
  ierr = PetscLogFlops(2 + 2 * e->ng + 8 * e->ng * e->nphiz * e->nphiy * e->nphix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecSetFromBC"
/*
  
  VecSetFromBC:
  
  (c) 2014 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VecSetFromBC(Vec BCVec,VFBC *BC)
{
  PetscErrorCode ierr;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       i,j,k,c;
  DM             da;
  PetscReal  ****BC_array,BCValue;
  PetscInt       dof;
  
  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject) BCVec,"DM",(PetscObject *) &da); CHKERRQ(ierr);
  if (!da) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Vector not generated from a DMDA");
  
  ierr = DMDAGetInfo(da,NULL,&nx,&ny,&nz,NULL,NULL,NULL,
                    &dof,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  
  ierr = DMDAVecGetArrayDOF(da,BCVec,&BC_array);CHKERRQ(ierr);
  
  for (c = 0;c < dof;c++) {
    /* 
      Faces
    */
    if (xs == 0 && BC[c].face[X0] != NONE) {
      /*
        x == 0
      */
      switch(BC[c].face[X0]) {
        case ZERO:
          BCValue = 0.;
          break;
        case ONE:
          BCValue = 1.;
        case FIXED:
          BCValue = BC[c].faceValue[X0];
        case NONE:
          break;
      }      
      for (k = zs; k < zs + zm; k++) {
        for (j = ys; j < ys + ym; j++) {
          BC_array[k][j][0][c] = BCValue;
        }
      }
    }

    if (xs + xm == nx && BC[c].face[X1] != NONE) {
      /*
        x == nx-1
      */
      switch(BC[c].face[X1]) {
        case ZERO:
          BCValue = 0.;
          break;
        case ONE:
          BCValue = 1.;
          break;
        case FIXED:
          BCValue = BC[c].faceValue[X1];
          break;
        case NONE:
          break;
      }      
      for (k = zs; k < zs + zm; k++) {
        for (j = ys; j < ys + ym; j++) {
          BC_array[k][j][nx-1][c] = BCValue;
        }
      }
    }

    if (ys == 0 && BC[c].face[Y0] != NONE) {
      /*
        y == 0
      */
      switch(BC[c].face[Y0]) {
        case ZERO:
          BCValue = 0.;
          break;
        case ONE:
          BCValue = 1.;
          break;
        case FIXED:
          BCValue = BC[c].faceValue[Y0];
          break;
        case NONE:
          break;
      }      
      for (k = zs; k < zs + zm; k++) {
        for (i = xs; i < xs + xm; i++) {
          BC_array[k][0][i][c] = BCValue;
        }
      }
    }
    if (ys + ym == ny && BC[c].face[Y1] != NONE) {
      /*
        y == ny-1
      */
      switch(BC[c].face[Y1]) {
        case ZERO:
          BCValue = 0.;
          break;
        case ONE:
          BCValue = 1.;
          break;
        case FIXED:
          BCValue = BC[c].faceValue[Y0];
          break;
        case NONE:
          break;
      }      
      for (k = zs; k < zs + zm; k++) {
        for (i = xs; i < xs + xm; i++) {
          BC_array[k][ny-1][i][c] = BCValue;
        }
      }
    }
    if (zs == 0 && BC[c].face[Z0] != NONE) {
        /*
          z == 0
        */
      switch(BC[c].face[Z0]) {
        case ZERO:
          BCValue = 0.;
          break;
        case ONE:
          BCValue = 1.;
          break;
        case FIXED:
          BCValue = BC[c].faceValue[Y0];
          break;
        case NONE:
          break;
      }      
      for (j = ys; j < ys + ym; j++) {
        for (i = xs; i < xs + xm; i++) {
          BC_array[0][j][i][c] = BCValue;
        }
      }
    }
    if (zs + zm == nz && BC[c].face[Z1] != NONE) {
      /*
        z == nz-1
      */
      switch(BC[c].face[Z0]) {
        case ZERO:
          BCValue = 0.;
          break;
        case ONE:
          BCValue = 1.;
          break;
        case FIXED:
          BCValue = BC[c].faceValue[Y0];
          break;
        case NONE:
          break;
      }      
      for (j = ys; j < ys + ym; j++) {
        for (i = xs; i < xs + xm; i++) {
          BC_array[nz-1][j][i][c] = BCValue;
        }
      }
    }
    /*
      edges
    */
    /* 
      vertices
    */
    if (xs == 0       && ys == 0       && zs == 0       && BC[c].vertex[X0Y0Z0] == FIXED) BC_array[0][0][0][c]          = BC[c].vertexValue[X0Y0Z0];
    if (xs == 0       && ys == 0       && zs == 0       && BC[c].vertex[X0Y0Z0] == ONE)   BC_array[0][0][0][c]          = 1.;
    if (xs == 0       && ys == 0       && zs == 0       && BC[c].vertex[X0Y0Z0] == ZERO)  BC_array[0][0][0][c]          = 0.;

    if (xs == 0       && ys == 0       && zs + zm == nz && BC[c].vertex[X0Y0Z1] == FIXED) BC_array[nz-1][0][0][c]       = BC[c].vertexValue[X0Y0Z1];
    if (xs == 0       && ys == 0       && zs + zm == nz && BC[c].vertex[X0Y0Z1] == ONE)   BC_array[nz-1][0][0][c]       = 1.;
    if (xs == 0       && ys == 0       && zs + zm == nz && BC[c].vertex[X0Y0Z1] == ZERO)  BC_array[nz-1][0][0][c]       = 0.;

    if (xs == 0       && ys + ym == ny && zs == 0       && BC[c].vertex[X0Y1Z0] == FIXED) BC_array[0][ny-1][0][c]       = BC[c].vertexValue[X0Y1Z0];
    if (xs == 0       && ys + ym == ny && zs == 0       && BC[c].vertex[X0Y1Z0] == ONE)   BC_array[0][ny-1][0][c]       = 1.;
    if (xs == 0       && ys + ym == ny && zs == 0       && BC[c].vertex[X0Y1Z0] == ZERO)  BC_array[0][ny-1][0][c]       = 0.;

    if (xs == 0       && ys + ym == ny && zs + zm == nz && BC[c].vertex[X0Y1Z1] == FIXED) BC_array[nz-1][ny-1][0][c]    = BC[c].vertexValue[X0Y1Z1];
    if (xs == 0       && ys + ym == ny && zs + zm == nz && BC[c].vertex[X0Y1Z1] == ONE)   BC_array[nz-1][ny-1][0][c]    = 1.;
    if (xs == 0       && ys + ym == ny && zs + zm == nz && BC[c].vertex[X0Y1Z1] == ZERO)  BC_array[nz-1][ny-1][0][c]    = 0.;

    if (xs + xm == nx && ys == 0       && zs == 0       && BC[c].vertex[X1Y0Z0] == FIXED) BC_array[0][0][nx-1][c]       = BC[c].vertexValue[X1Y0Z0];
    if (xs + xm == nx && ys == 0       && zs == 0       && BC[c].vertex[X1Y0Z0] == ONE)   BC_array[0][0][nx-1][c]       = 1.;
    if (xs + xm == nx && ys == 0       && zs == 0       && BC[c].vertex[X1Y0Z0] == ZERO)  BC_array[0][0][nx-1][c]       = 0.;

    if (xs + xm == nx && ys == 0       && zs + zm == nz && BC[c].vertex[X1Y0Z1] == FIXED) BC_array[nz-1][0][nx-1][c]    = BC[c].vertexValue[X1Y0Z1];
    if (xs + xm == nx && ys == 0       && zs + zm == nz && BC[c].vertex[X1Y0Z1] == ONE)   BC_array[nz-1][0][nx-1][c]    = 1.;
    if (xs + xm == nx && ys == 0       && zs + zm == nz && BC[c].vertex[X1Y0Z1] == ZERO)  BC_array[nz-1][0][nx-1][c]    = 0.;

    if (xs + xm == nx && ys + ym == ny && zs == 0       && BC[c].vertex[X1Y1Z0] == FIXED) BC_array[0][ny-1][nx-1][c]    = BC[c].vertexValue[X1Y1Z0];
    if (xs + xm == nx && ys + ym == ny && zs == 0       && BC[c].vertex[X1Y1Z0] == ONE)   BC_array[0][ny-1][nx-1][c]    = 1.;
    if (xs + xm == nx && ys + ym == ny && zs == 0       && BC[c].vertex[X1Y1Z0] == ZERO)  BC_array[0][ny-1][nx-1][c]    = 0.;

    if (xs + xm == nx && ys + ym == ny && zs + zm == nz && BC[c].vertex[X1Y1Z1] == FIXED) BC_array[nz-1][ny-1][nx-1][c] = BC[c].vertexValue[X1Y1Z1];
    if (xs + xm == nx && ys + ym == ny && zs + zm == nz && BC[c].vertex[X1Y1Z1] == ONE)   BC_array[nz-1][ny-1][nx-1][c] = 1.;
    if (xs + xm == nx && ys + ym == ny && zs + zm == nz && BC[c].vertex[X1Y1Z1] == ZERO)  BC_array[nz-1][ny-1][nx-1][c] = 0.;
  }
  
  ierr = DMDAVecRestoreArrayDOF(da,BCVec,&BC_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
  

#undef __FUNCT__
#define __FUNCT__ "VecApplyDirichletBC"
/*
  VecApplyDirichletBC

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VecApplyDirichletBC(Vec RHS,Vec BCU,VFBC *BC)
{
  PetscErrorCode ierr;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       i,j,k,c;
  DM             da;
  PetscReal  ****RHS_array;
  PetscReal  ****BCU_array;
  PetscInt       dim,dof;
  
  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject) RHS,"DM",(PetscObject *) &da); CHKERRQ(ierr);
    if (!da) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Vector not generated from a DMDA");
  
  ierr = DMDAGetInfo(da,&dim,&nx,&ny,&nz,NULL,NULL,NULL,
                    &dof,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  
  if (dim == 2) {
    ierr = PetscMalloc(sizeof(PetscReal ***),&RHS_array);CHKERRQ(ierr);
    ierr = PetscMalloc(sizeof(PetscReal ***),&BCU_array);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(da,RHS,&RHS_array[0]);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(da,BCU,&BCU_array[0]);CHKERRQ(ierr);
  } else {
    ierr = DMDAVecGetArrayDOF(da,RHS,&RHS_array);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(da,BCU,&BCU_array);CHKERRQ(ierr);
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
      if (BC[c].face[X0] == FIXED) {
        for (k = zs; k < zs + zm; k++) {
          for (j = ys; j < ys + ym; j++) {
            RHS_array[k][j][i][c] = BCU_array[k][j][i][c];
          }
        }
      }
      if (BC[c].face[X0] == ONE) {
        for (k = zs; k < zs + zm; k++) {
          for (j = ys; j < ys + ym; j++) {
            RHS_array[k][j][i][c] = 1.;
          }
        }
      }
      if (BC[c].face[X0] == ZERO) {
        for (k = zs; k < zs + zm; k++) {
          for (j = ys; j < ys + ym; j++) {
            RHS_array[k][j][i][c] = 0.;
          }
        }
      }
    }
    if (xs + xm == nx) {
      /*
        x == nx-1
      */
      i = nx-1;
      if (BC[c].face[X1] == FIXED) {
        for (k = zs; k < zs + zm; k++) {
          for (j = ys; j < ys + ym; j++) {
            RHS_array[k][j][i][c] = BCU_array[k][j][i][c];
          }
        }
      }
      if (BC[c].face[X1] == ONE) {
        for (k = zs; k < zs + zm; k++) {
          for (j = ys; j < ys + ym; j++) {
            RHS_array[k][j][i][c] = 1.;
          }
        }
      }
      if (BC[c].face[X1] == ZERO) {
        for (k = zs; k < zs + zm; k++) {
          for (j = ys; j < ys + ym; j++) {
            RHS_array[k][j][i][c] = 0.;
          }
        }
      }
    }
    if (ys == 0) {
      /*
        y == 0
      */
      j = 0;
      if (BC[c].face[Y0] == FIXED) {
        for (k = zs; k < zs + zm; k++) {
          for (i = xs; i < xs + xm; i++) {
            RHS_array[k][j][i][c] = BCU_array[k][j][i][c];
          }
        }
      }
      if (BC[c].face[Y0] == ONE) {
        for (k = zs; k < zs + zm; k++) {
          for (i = xs; i < xs + xm; i++) {
            RHS_array[k][j][i][c] = 1.;
          }
        }
      }
      if (BC[c].face[Y0] == ZERO) {
        for (k = zs; k < zs + zm; k++) {
          for (i = xs; i < xs + xm; i++) {
            RHS_array[k][j][i][c] = 0.;
          }
        }
      }
    }
    if (ys + ym == ny) {
      /*
        y == ny-1
      */
      j = ny-1;
      if (BC[c].face[Y1] == FIXED) {
        for (k = zs; k < zs + zm; k++) {
          for (i = xs; i < xs + xm; i++) {
            RHS_array[k][j][i][c] = BCU_array[k][j][i][c];
          }
        }
      }
      if (BC[c].face[Y1] == ONE) {
        for (k = zs; k < zs + zm; k++) {
          for (i = xs; i < xs + xm; i++) {
            RHS_array[k][j][i][c] = 1.;
          }
        }
      }
      if (BC[c].face[Y1] == ZERO) {
        for (k = zs; k < zs + zm; k++) {
          for (i = xs; i < xs + xm; i++) {
            RHS_array[k][j][i][c] = 0.;
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
        if (BC[c].face[Z0] == FIXED) {
          for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
              RHS_array[k][j][i][c] = BCU_array[k][j][i][c];
            }
          }
        }
        if (BC[c].face[Z0] == ONE) {
          for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
              RHS_array[k][j][i][c] = 1.;
            }
          }
        }
        if (BC[c].face[Z0] == ZERO) {
          for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
              RHS_array[k][j][i][c] = 0.;
            }
          }
        }
      }
      if (zs + zm == nz) {
        /*
          z == nz-1
        */
        k = nz-1;
        if (BC[c].face[Z1] == FIXED) {
          for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
              RHS_array[k][j][i][c] = BCU_array[k][j][i][c];
            }
          }
        }
        if (BC[c].face[Z1] == ONE) {
          for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
              RHS_array[k][j][i][c] = 1.;
            }
          }
        }
        if (BC[c].face[Z1] == ZERO) {
          for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
              RHS_array[k][j][i][c] = 0.;
            }
          }
        }
      }
    }
    /*
      edges
    */
    /* 
      vertices
    */
    if (xs == 0       && ys == 0       && zs == 0       && BC[c].vertex[X0Y0Z0] == FIXED) RHS_array[0][0][0][c]          = BCU_array[0][0][0][c];
    if (xs == 0       && ys == 0       && zs == 0       && BC[c].vertex[X0Y0Z0] == ONE)   RHS_array[0][0][0][c]          = 1.;
    if (xs == 0       && ys == 0       && zs == 0       && BC[c].vertex[X0Y0Z0] == ZERO)  RHS_array[0][0][0][c]          = 0.;

    if (xs == 0       && ys == 0       && zs + zm == nz && BC[c].vertex[X0Y0Z1] == FIXED) RHS_array[nz-1][0][0][c]       = BCU_array[nz-1][0][0][c];
    if (xs == 0       && ys == 0       && zs + zm == nz && BC[c].vertex[X0Y0Z1] == ONE)   RHS_array[nz-1][0][0][c]       = 1.;
    if (xs == 0       && ys == 0       && zs + zm == nz && BC[c].vertex[X0Y0Z1] == ZERO)  RHS_array[nz-1][0][0][c]       = 0.;

    if (xs == 0       && ys + ym == ny && zs == 0       && BC[c].vertex[X0Y1Z0] == FIXED) RHS_array[0][ny-1][0][c]       = BCU_array[0][ny-1][0][c];
    if (xs == 0       && ys + ym == ny && zs == 0       && BC[c].vertex[X0Y1Z0] == ONE)   RHS_array[0][ny-1][0][c]       = 1.;
    if (xs == 0       && ys + ym == ny && zs == 0       && BC[c].vertex[X0Y1Z0] == ZERO)  RHS_array[0][ny-1][0][c]       = 0.;

    if (xs == 0       && ys + ym == ny && zs + zm == nz && BC[c].vertex[X0Y1Z1] == FIXED) RHS_array[nz-1][ny-1][0][c]    = BCU_array[nz-1][ny-1][0][c];
    if (xs == 0       && ys + ym == ny && zs + zm == nz && BC[c].vertex[X0Y1Z1] == ONE)   RHS_array[nz-1][ny-1][0][c]    = 1.;
    if (xs == 0       && ys + ym == ny && zs + zm == nz && BC[c].vertex[X0Y1Z1] == ZERO)  RHS_array[nz-1][ny-1][0][c]    = 0.;

    if (xs + xm == nx && ys == 0       && zs == 0       && BC[c].vertex[X1Y0Z0] == FIXED) RHS_array[0][0][nx-1][c]       = BCU_array[0][0][nx-1][c];
    if (xs + xm == nx && ys == 0       && zs == 0       && BC[c].vertex[X1Y0Z0] == ONE)   RHS_array[0][0][nx-1][c]       = 1.;
    if (xs + xm == nx && ys == 0       && zs == 0       && BC[c].vertex[X1Y0Z0] == ZERO)  RHS_array[0][0][nx-1][c]       = 0.;

    if (xs + xm == nx && ys == 0       && zs + zm == nz && BC[c].vertex[X1Y0Z1] == FIXED) RHS_array[nz-1][0][nx-1][c]    = BCU_array[nz-1][0][nx-1][c];
    if (xs + xm == nx && ys == 0       && zs + zm == nz && BC[c].vertex[X1Y0Z1] == ONE)   RHS_array[nz-1][0][nx-1][c]    = 1.;
    if (xs + xm == nx && ys == 0       && zs + zm == nz && BC[c].vertex[X1Y0Z1] == ZERO)  RHS_array[nz-1][0][nx-1][c]    = 0.;

    if (xs + xm == nx && ys + ym == ny && zs == 0       && BC[c].vertex[X1Y1Z0] == FIXED) RHS_array[0][ny-1][nx-1][c]    = BCU_array[0][ny-1][nx-1][c];
    if (xs + xm == nx && ys + ym == ny && zs == 0       && BC[c].vertex[X1Y1Z0] == ONE)   RHS_array[0][ny-1][nx-1][c]    = 1.;
    if (xs + xm == nx && ys + ym == ny && zs == 0       && BC[c].vertex[X1Y1Z0] == ZERO)  RHS_array[0][ny-1][nx-1][c]    = 0.;

    if (xs + xm == nx && ys + ym == ny && zs + zm == nz && BC[c].vertex[X1Y1Z1] == FIXED) RHS_array[nz-1][ny-1][nx-1][c] = BCU_array[nz-1][ny-1][nx-1][c];
    if (xs + xm == nx && ys + ym == ny && zs + zm == nz && BC[c].vertex[X1Y1Z1] == ONE)   RHS_array[nz-1][ny-1][nx-1][c] = 1.;
    if (xs + xm == nx && ys + ym == ny && zs + zm == nz && BC[c].vertex[X1Y1Z1] == ZERO)  RHS_array[nz-1][ny-1][nx-1][c] = 0.;
  }
  
  if (dim == 2) {
    ierr = DMDAVecRestoreArrayDOF(da,RHS,&RHS_array[0]);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(da,BCU,&BCU_array[0]);CHKERRQ(ierr);
    ierr = PetscFree(RHS_array);CHKERRQ(ierr);
    ierr = PetscFree(BCU_array);CHKERRQ(ierr);
  } else {
    ierr = DMDAVecRestoreArrayDOF(da,RHS,&RHS_array);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(da,BCU,&BCU_array);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ResidualApplyDirichletBC"
/*
  ResidualApplyDirichletBC

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode ResidualApplyDirichletBC(Vec residual,Vec U,Vec BCU,VFBC *BC)
{
  PetscErrorCode ierr;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       i,j,k,c;
  DM             da;
  PetscReal  ****residual_array;
  PetscReal  ****BCU_array,****U_array;
  PetscInt       dim,dof;
  
  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject) residual,"DM",(PetscObject *) &da); CHKERRQ(ierr);
    if (!da) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Vector not generated from a DMDA");
  
  ierr = DMDAGetInfo(da,&dim,&nx,&ny,&nz,NULL,NULL,NULL,
                    &dof,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  
  if (dim == 2) {
    ierr = PetscMalloc(sizeof(PetscReal ***),&residual_array);CHKERRQ(ierr);
    ierr = PetscMalloc(sizeof(PetscReal ***),&BCU_array);CHKERRQ(ierr);
    ierr = PetscMalloc(sizeof(PetscReal ***),&U_array);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(da,residual,&residual_array[0]);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(da,BCU,&BCU_array[0]);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(da,U,&U_array[0]);CHKERRQ(ierr);
  } else {
    ierr = DMDAVecGetArrayDOF(da,residual,&residual_array);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(da,BCU,&BCU_array);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(da,U,&U_array);CHKERRQ(ierr);
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
      if (BC[c].face[X0] == FIXED) {
        for (k = zs; k < zs + zm; k++) {
          for (j = ys; j < ys + ym; j++) {
            residual_array[k][j][i][c] = U_array[k][j][i][c]-BCU_array[k][j][i][c];
          }
        }
      }
      if (BC[c].face[X0] == ONE) {
        for (k = zs; k < zs + zm; k++) {
          for (j = ys; j < ys + ym; j++) {
            residual_array[k][j][i][c] = U_array[k][j][i][c]-1.;
          }
        }
      }
      if (BC[c].face[X0] == ZERO) {
        for (k = zs; k < zs + zm; k++) {
          for (j = ys; j < ys + ym; j++) {
            residual_array[k][j][i][c] = U_array[k][j][i][c];
          }
        }
      }
    }
    if (xs + xm == nx) {
      /*
        x == nx-1
      */
      i = nx-1;
      if (BC[c].face[X1] == FIXED) {
        for (k = zs; k < zs + zm; k++) {
          for (j = ys; j < ys + ym; j++) {
            residual_array[k][j][i][c] = U_array[k][j][i][c]-BCU_array[k][j][i][c];
          }
        }
      }
      if (BC[c].face[X1] == ONE) {
        for (k = zs; k < zs + zm; k++) {
          for (j = ys; j < ys + ym; j++) {
            residual_array[k][j][i][c] = U_array[k][j][i][c]-1.;
          }
        }
      }
      if (BC[c].face[X1] == ZERO) {
        for (k = zs; k < zs + zm; k++) {
          for (j = ys; j < ys + ym; j++) {
            residual_array[k][j][i][c] = U_array[k][j][i][c];
          }
        }
      }
    }
    if (ys == 0) {
      /*
        y == 0
      */
      j = 0;
      if (BC[c].face[Y0] == FIXED) {
        for (k = zs; k < zs + zm; k++) {
          for (i = xs; i < xs + xm; i++) {
            residual_array[k][j][i][c] = U_array[k][j][i][c]-BCU_array[k][j][i][c];
          }
        }
      }
      if (BC[c].face[Y0] == ONE) {
        for (k = zs; k < zs + zm; k++) {
          for (i = xs; i < xs + xm; i++) {
            residual_array[k][j][i][c] = U_array[k][j][i][c]-1;
          }
        }
      }
      if (BC[c].face[Y0] == ZERO) {
        for (k = zs; k < zs + zm; k++) {
          for (i = xs; i < xs + xm; i++) {
            residual_array[k][j][i][c] = U_array[k][j][i][c];
          }
        }
      }
    }
    if (ys + ym == ny) {
      /*
        y == ny-1
      */
      j = ny-1;
      if (BC[c].face[Y1] == FIXED) {
        for (k = zs; k < zs + zm; k++) {
          for (i = xs; i < xs + xm; i++) {
            residual_array[k][j][i][c] = U_array[k][j][i][c]-BCU_array[k][j][i][c];
          }
        }
      }
      if (BC[c].face[Y1] == ONE) {
        for (k = zs; k < zs + zm; k++) {
          for (i = xs; i < xs + xm; i++) {
            residual_array[k][j][i][c] = U_array[k][j][i][c]-1.;
          }
        }
      }
      if (BC[c].face[Y1] == ZERO) {
        for (k = zs; k < zs + zm; k++) {
          for (i = xs; i < xs + xm; i++) {
            residual_array[k][j][i][c] = U_array[k][j][i][c];
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
        if (BC[c].face[Z0] == FIXED) {
          for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
              residual_array[k][j][i][c] = U_array[k][j][i][c]-BCU_array[k][j][i][c];
            }
          }
        }
        if (BC[c].face[Z0] == ONE) {
          for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
              residual_array[k][j][i][c] = U_array[k][j][i][c]-1.;
            }
          }
        }
        if (BC[c].face[Z0] == ZERO) {
          for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
              residual_array[k][j][i][c] = U_array[k][j][i][c];
            }
          }
        }
      }
      if (zs + zm == nz) {
        /*
          z == nz-1
        */
        k = nz-1;
        if (BC[c].face[Z1] == FIXED) {
          for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
              residual_array[k][j][i][c] = U_array[k][j][i][c]-BCU_array[k][j][i][c];
            }
          }
        }
        if (BC[c].face[Z1] == ONE) {
          for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
              residual_array[k][j][i][c] = U_array[k][j][i][c]-1.;
            }
          }
        }
        if (BC[c].face[Z1] == ZERO) {
          for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
              residual_array[k][j][i][c] = U_array[k][j][i][c];
            }
          }
        }
      }
    }
    /*
      edges
    */
    /* 
      vertices
    */
    if (xs == 0       && ys == 0       && zs == 0       && BC[c].vertex[X0Y0Z0] == FIXED) residual_array[0][0][0][c]          = U_array[0][0][0][c]-BCU_array[0][0][0][c];
    if (xs == 0       && ys == 0       && zs == 0       && BC[c].vertex[X0Y0Z0] == ONE)   residual_array[0][0][0][c]          = U_array[0][0][0][c]-1.;
    if (xs == 0       && ys == 0       && zs == 0       && BC[c].vertex[X0Y0Z0] == ZERO)  residual_array[0][0][0][c]          = U_array[0][0][0][c];

    if (xs == 0       && ys == 0       && zs + zm == nz && BC[c].vertex[X0Y0Z1] == FIXED) residual_array[nz-1][0][0][c]       = U_array[nz-1][0][0][c]-BCU_array[nz-1][0][0][c];
    if (xs == 0       && ys == 0       && zs + zm == nz && BC[c].vertex[X0Y0Z1] == ONE)   residual_array[nz-1][0][0][c]       = U_array[nz-1][0][0][c]-1.;
    if (xs == 0       && ys == 0       && zs + zm == nz && BC[c].vertex[X0Y0Z1] == ZERO)  residual_array[nz-1][0][0][c]       = U_array[nz-1][0][0][c];

    if (xs == 0       && ys + ym == ny && zs == 0       && BC[c].vertex[X0Y1Z0] == FIXED) residual_array[0][ny-1][0][c]       = U_array[0][ny-1][0][c]-BCU_array[0][ny-1][0][c];
    if (xs == 0       && ys + ym == ny && zs == 0       && BC[c].vertex[X0Y1Z0] == ONE)   residual_array[0][ny-1][0][c]       = U_array[0][ny-1][0][c]-1.;
    if (xs == 0       && ys + ym == ny && zs == 0       && BC[c].vertex[X0Y1Z0] == ZERO)  residual_array[0][ny-1][0][c]       = U_array[0][ny-1][0][c];

    if (xs == 0       && ys + ym == ny && zs + zm == nz && BC[c].vertex[X0Y1Z1] == FIXED) residual_array[nz-1][ny-1][0][c]    = U_array[nz-1][ny-1][0][c]-BCU_array[nz-1][ny-1][0][c];
    if (xs == 0       && ys + ym == ny && zs + zm == nz && BC[c].vertex[X0Y1Z1] == ONE)   residual_array[nz-1][ny-1][0][c]    = U_array[nz-1][ny-1][0][c]-1.;
    if (xs == 0       && ys + ym == ny && zs + zm == nz && BC[c].vertex[X0Y1Z1] == ZERO)  residual_array[nz-1][ny-1][0][c]    = U_array[nz-1][ny-1][0][c];

    if (xs + xm == nx && ys == 0       && zs == 0       && BC[c].vertex[X1Y0Z0] == FIXED) residual_array[0][0][nx-1][c]       = U_array[0][0][nx-1][c]-BCU_array[0][0][nx-1][c];
    if (xs + xm == nx && ys == 0       && zs == 0       && BC[c].vertex[X1Y0Z0] == ONE)   residual_array[0][0][nx-1][c]       = U_array[0][0][nx-1][c]-1.;
    if (xs + xm == nx && ys == 0       && zs == 0       && BC[c].vertex[X1Y0Z0] == ZERO)  residual_array[0][0][nx-1][c]       = U_array[0][0][nx-1][c];

    if (xs + xm == nx && ys == 0       && zs + zm == nz && BC[c].vertex[X1Y0Z1] == FIXED) residual_array[nz-1][0][nx-1][c]    = U_array[nz-1][0][nx-1][c]-BCU_array[nz-1][0][nx-1][c];
    if (xs + xm == nx && ys == 0       && zs + zm == nz && BC[c].vertex[X1Y0Z1] == ONE)   residual_array[nz-1][0][nx-1][c]    = U_array[nz-1][0][nx-1][c]-1.;
    if (xs + xm == nx && ys == 0       && zs + zm == nz && BC[c].vertex[X1Y0Z1] == ZERO)  residual_array[nz-1][0][nx-1][c]    = U_array[nz-1][0][nx-1][c];

    if (xs + xm == nx && ys + ym == ny && zs == 0       && BC[c].vertex[X1Y1Z0] == FIXED) residual_array[0][ny-1][nx-1][c]    = U_array[0][ny-1][nx-1][c]-BCU_array[0][ny-1][nx-1][c];
    if (xs + xm == nx && ys + ym == ny && zs == 0       && BC[c].vertex[X1Y1Z0] == ONE)   residual_array[0][ny-1][nx-1][c]    = U_array[0][ny-1][nx-1][c]-1.;
    if (xs + xm == nx && ys + ym == ny && zs == 0       && BC[c].vertex[X1Y1Z0] == ZERO)  residual_array[0][ny-1][nx-1][c]    = U_array[0][ny-1][nx-1][c];

    if (xs + xm == nx && ys + ym == ny && zs + zm == nz && BC[c].vertex[X1Y1Z1] == FIXED) residual_array[nz-1][ny-1][nx-1][c] = U_array[nz-1][ny-1][nx-1][c]-BCU_array[nz-1][ny-1][nx-1][c];
    if (xs + xm == nx && ys + ym == ny && zs + zm == nz && BC[c].vertex[X1Y1Z1] == ONE)   residual_array[nz-1][ny-1][nx-1][c] = U_array[nz-1][ny-1][nx-1][c]-1.;
    if (xs + xm == nx && ys + ym == ny && zs + zm == nz && BC[c].vertex[X1Y1Z1] == ZERO)  residual_array[nz-1][ny-1][nx-1][c] = U_array[nz-1][ny-1][nx-1][c];
  }
  
  if (dim == 2) {
    ierr = DMDAVecRestoreArrayDOF(da,residual,&residual_array[0]);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(da,BCU,&BCU_array[0]);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(da,U,&U_array[0]);CHKERRQ(ierr);
    ierr = PetscFree(residual_array);CHKERRQ(ierr);
    ierr = PetscFree(BCU_array);CHKERRQ(ierr);
    ierr = PetscFree(U_array);CHKERRQ(ierr);
  } else {
    ierr = DMDAVecRestoreArrayDOF(da,residual,&residual_array);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(da,BCU,&BCU_array);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(da,U,&U_array);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatApplyDirichletBC"
/*
  MatApplyDirichletBC

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode MatApplyDirichletBC(Mat K,DM da,VFBC *BC)
{
  PetscErrorCode ierr;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       i,j,k,c;
  MatStencil    *row;
  PetscReal      one=1.;
  PetscInt       numBC=0,l=0;
  PetscInt       dim,dof;

  PetscFunctionBegin;
  
  /*
    This is only implemented in petsc-dev (as of petsc-3.1 days)
  */
  /*
    ierr = PetscObjectQuery((PetscObject) K,"DM",(PetscObject *) &da); CHKERRQ(ierr);
    if (!da) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG," Matrix not generated from a DA");
  */
  ierr = DMDAGetInfo(da,&dim,&nx,&ny,&nz,NULL,NULL,NULL,
                    &dof,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  /*
    Compute the number of boundary nodes on each processor. 
    Edges and corners are counted multiple times (2 and 3 resp)
  */
  for (c = 0; c < dof; c++){
    if (xs == 0       && BC[c].face[X0] != NONE)             numBC += ym * zm;
    if (xs + xm == nx && BC[c].face[X1] != NONE)             numBC += ym * zm;
    if (ys == 0       && BC[c].face[Y0] != NONE)             numBC += xm * zm;
    if (ys + ym == ny && BC[c].face[Y1] != NONE)             numBC += xm * zm;
    if (zs == 0       && BC[c].face[Z0] != NONE && dim == 3) numBC += xm * ym;
    if (zs + zm == nz && BC[c].face[Z1] != NONE && dim == 3) numBC += xm * ym;
    if (xs == 0       && ys == 0       && zs == 0       && BC[c].vertex[X0Y0Z0] != NONE) numBC++;
    if (xs == 0       && ys + ym == ny && zs == 0       && BC[c].vertex[X0Y1Z0] != NONE) numBC++;
    if (xs + xm == nx && ys == 0       && zs == 0       && BC[c].vertex[X1Y0Z0] != NONE) numBC++;
    if (xs + xm == nx && ys + ym == ny && zs == 0       && BC[c].vertex[X1Y1Z0] != NONE) numBC++;
    if (xs == 0       && ys == 0       && zs + zm == nz && BC[c].vertex[X0Y0Z1] != NONE && dim == 3) numBC++;
    if (xs == 0       && ys + ym == ny && zs + zm == nz && BC[c].vertex[X0Y1Z1] != NONE && dim == 3) numBC++;
    if (xs + xm == nx && ys == 0       && zs + zm == nz && BC[c].vertex[X1Y0Z1] != NONE && dim == 3) numBC++;
    if (xs + xm == nx && ys + ym == ny && zs + zm == nz && BC[c].vertex[X1Y1Z1] != NONE && dim == 3) numBC++;
  }
  ierr = PetscMalloc(numBC * sizeof(MatStencil),&row);CHKERRQ(ierr);
  /*
    Create an array of rows to be zeroed out
  */
  /*
    i == 0
  */
  for (c = 0; c < dof; c++) {
    if (xs == 0 && BC[c].face[X0] != NONE) {
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
    if (xs + xm == nx && BC[c].face[X1] != NONE) {
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
    if (ys == 0 && BC[c].face[Y0] != NONE) {
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
    if (ys + ym == ny && BC[c].face[Y1] != NONE) {
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
      if (zs == 0 && BC[c].face[Z0] != NONE) {
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
      if (zs + zm == nz && BC[c].face[Z1] != NONE) {
        for (j = ys; j < ys + ym; j++) {
          for (i = xs; i < xs + xm; i++) {
            row[l].i = i; row[l].j = j; row[l].k = nz-1; row[l].c = c; 
            l++;
          }
        }
      }
    }
    if (xs == 0       && ys == 0       && zs == 0       && BC[c].vertex[X0Y0Z0] != NONE) { 
      row[l].i = 0; row[l].j = 0; row[l].k = 0; row[l].c = c; 
      l++;
    }
    if (xs == 0       && ys == 0       && zs + zm == nz && BC[c].vertex[X0Y0Z1] != NONE && dim ==3) { 
      row[l].i = 0; row[l].j = 0; row[l].k = nz-1; row[l].c = c; 
      l++;
    }
    if (xs == 0       && ys + ym == ny && zs == 0       && BC[c].vertex[X0Y1Z0] != NONE) { 
      row[l].i = 0; row[l].j = ny-1; row[l].k = 0; row[l].c = c; 
      l++;
    }
    if (xs == 0       && ys + ym == ny && zs + zm == nz && BC[c].vertex[X0Y1Z1] != NONE && dim ==3) { 
      row[l].i = 0; row[l].j = ny-1; row[l].k = nz-1; row[l].c = c; 
      l++;
    }
    if (xs + xm == nx && ys == 0       && zs == 0       && BC[c].vertex[X1Y0Z0] != NONE) { 
      row[l].i = nx-1; row[l].j = 0; row[l].k = 0; row[l].c = c; 
      l++;
    }
    if (xs + xm == nx && ys == 0       && zs + zm == nz && BC[c].vertex[X1Y0Z1] != NONE && dim ==3) { 
      row[l].i = nx-1; row[l].j = 0; row[l].k = nz-1; row[l].c = c; 
      l++;
    }
    if (xs + xm == nx && ys + ym == ny && zs == 0       && BC[c].vertex[X1Y1Z0] != NONE) { 
      row[l].i = nx-1; row[l].j = ny-1; row[l].k = 0; row[l].c = c; 
      l++;
    }
    if (xs + xm == nx && ys + ym == ny && zs + zm == nz && BC[c].vertex[X1Y1Z1] != NONE && dim ==3) { 
      row[l].i = nx-1; row[l].j = ny-1; row[l].k = nz-1; row[l].c = c; 
      l++;
    }
    
  }
  ierr = MatZeroRowsStencil(K,numBC,row,one,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscFree(row);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatApplyDirichletBCRowCol"
/*
  MatApplyDirichletBCRowCol

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode MatApplyDirichletBCRowCol(Mat K,DM da,VFBC *BC)
{
  PetscErrorCode ierr;
  PetscInt       xs,xm,nx;
  PetscInt       ys,ym,ny;
  PetscInt       zs,zm,nz;
  PetscInt       i,j,k,c;
  MatStencil    *row;
  PetscReal      one=1.;
  PetscInt       numBC=0,l=0;
  PetscInt       dim,dof;

  PetscFunctionBegin;
  
  /*
    This is only implemented in petsc-dev (as of petsc-3.1 days)
  */
  /*
    ierr = PetscObjectQuery((PetscObject) K,"DM",(PetscObject *) &da); CHKERRQ(ierr);
    if (!da) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG," Matrix not generated from a DA");
  */
  ierr = DMDAGetInfo(da,&dim,&nx,&ny,&nz,NULL,NULL,NULL,
                    &dof,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  /*
    Compute the number of boundary nodes on each processor. 
    Edges and corners are counted multiple times (2 and 3 resp)
  */
  for (c = 0; c < dof; c++){
    if (xs == 0       && BC[c].face[X0] != NONE)             numBC += ym * zm;
    if (xs + xm == nx && BC[c].face[X1] != NONE)             numBC += ym * zm;
    if (ys == 0       && BC[c].face[Y0] != NONE)             numBC += xm * zm;
    if (ys + ym == ny && BC[c].face[Y1] != NONE)             numBC += xm * zm;
    if (zs == 0       && BC[c].face[Z0] != NONE && dim == 3) numBC += xm * ym;
    if (zs + zm == nz && BC[c].face[Z1] != NONE && dim == 3) numBC += xm * ym;
    if (xs == 0       && ys == 0       && zs == 0       && BC[c].vertex[X0Y0Z0] != NONE) numBC++;
    if (xs == 0       && ys + ym == ny && zs == 0       && BC[c].vertex[X0Y1Z0] != NONE) numBC++;
    if (xs + xm == nx && ys == 0       && zs == 0       && BC[c].vertex[X1Y0Z0] != NONE) numBC++;
    if (xs + xm == nx && ys + ym == ny && zs == 0       && BC[c].vertex[X1Y1Z0] != NONE) numBC++;
    if (xs == 0       && ys == 0       && zs + zm == nz && BC[c].vertex[X0Y0Z1] != NONE && dim == 3) numBC++;
    if (xs == 0       && ys + ym == ny && zs + zm == nz && BC[c].vertex[X0Y1Z1] != NONE && dim == 3) numBC++;
    if (xs + xm == nx && ys == 0       && zs + zm == nz && BC[c].vertex[X1Y0Z1] != NONE && dim == 3) numBC++;
    if (xs + xm == nx && ys + ym == ny && zs + zm == nz && BC[c].vertex[X1Y1Z1] != NONE && dim == 3) numBC++;
  }
  ierr = PetscMalloc(numBC * sizeof(MatStencil),&row);CHKERRQ(ierr);
  /*
    Create an array of rows to be zeroed out
  */
  /*
    i == 0
  */
  for (c = 0; c < dof; c++) {
    if (xs == 0 && BC[c].face[X0] != NONE) {
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
    if (xs + xm == nx && BC[c].face[X1] != NONE) {
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
    if (ys == 0 && BC[c].face[Y0] != NONE) {
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
    if (ys + ym == ny && BC[c].face[Y1] != NONE) {
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
      if (zs == 0 && BC[c].face[Z0] != NONE) {
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
      if (zs + zm == nz && BC[c].face[Z1] != NONE) {
        for (j = ys; j < ys + ym; j++) {
          for (i = xs; i < xs + xm; i++) {
            row[l].i = i; row[l].j = j; row[l].k = nz-1; row[l].c = c; 
            l++;
          }
        }
      }
    }
    if (xs == 0       && ys == 0       && zs == 0       && BC[c].vertex[X0Y0Z0] != NONE) { 
      row[l].i = 0; row[l].j = 0; row[l].k = 0; row[l].c = c; 
      l++;
    }
    if (xs == 0       && ys == 0       && zs + zm == nz && BC[c].vertex[X0Y0Z1] != NONE && dim ==3) { 
      row[l].i = 0; row[l].j = 0; row[l].k = nz-1; row[l].c = c; 
      l++;
    }
    if (xs == 0       && ys + ym == ny && zs == 0       && BC[c].vertex[X0Y1Z0] != NONE) { 
      row[l].i = 0; row[l].j = ny-1; row[l].k = 0; row[l].c = c; 
      l++;
    }
    if (xs == 0       && ys + ym == ny && zs + zm == nz && BC[c].vertex[X0Y1Z1] != NONE && dim ==3) { 
      row[l].i = 0; row[l].j = ny-1; row[l].k = nz-1; row[l].c = c; 
      l++;
    }
    if (xs + xm == nx && ys == 0       && zs == 0       && BC[c].vertex[X1Y0Z0] != NONE) { 
      row[l].i = nx-1; row[l].j = 0; row[l].k = 0; row[l].c = c; 
      l++;
    }
    if (xs + xm == nx && ys == 0       && zs + zm == nz && BC[c].vertex[X1Y0Z1] != NONE && dim ==3) { 
      row[l].i = nx-1; row[l].j = 0; row[l].k = nz-1; row[l].c = c; 
      l++;
    }
    if (xs + xm == nx && ys + ym == ny && zs == 0       && BC[c].vertex[X1Y1Z0] != NONE) { 
      row[l].i = nx-1; row[l].j = ny-1; row[l].k = 0; row[l].c = c; 
      l++;
    }
    if (xs + xm == nx && ys + ym == ny && zs + zm == nz && BC[c].vertex[X1Y1Z1] != NONE && dim ==3) { 
      row[l].i = nx-1; row[l].j = ny-1; row[l].k = nz-1; row[l].c = c; 
      l++;
    }
    
  }
  ierr = MatZeroRowsColumnsStencil(K,numBC,row,one,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscFree(row);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAReadCoordinatesHDF5"
/*
  DAReadCoordinatesHDF5

  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode DAReadCoordinatesHDF5(DM da,const char filename[])
{
#ifdef PETSC_HAVE_HDF5
  PetscErrorCode ierr;
  Vec            Coords;
  PetscViewer    HDF5Viewer;
  
  PetscFunctionBegin;
  ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&HDF5Viewer);
  ierr = DMGetGlobalVector(da,&Coords);CHKERRQ(ierr);
  ierr = VecLoad(Coords,HDF5Viewer);CHKERRQ(ierr);
  ierr = DMDASetCoordinates(da,Coords);CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(da,&Coords);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&HDF5Viewer);CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFBCCreate"
/*
  VFBCCreate: Create an array of BC and set all boundary conditions to FREE

  (c) 2010-2014 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFBCCreate(VFBC *bc,PetscInt dof)
{
  PetscInt       i,c;

  PetscFunctionBegin;
  /*
    faces 
  */
  for (i = 0; i < 6; i++) {
    for (c = 0; c < dof; c++) {
      bc[c].face[i] = NONE;
      bc[c].faceValue[i] = 0.;
    }
  }
  /*
    edges 
  */
  for (i = 0; i < 12; i++) {
    for (c = 0; c < dof; c++) {
      bc[c].edge[i] = NONE;
      bc[c].edgeValue[i] = 0.;
    }
  }
  /*
    vertices 
  */
  for (i = 0; i < 8; i++) {
    for (c = 0; c < dof; c++) {
      bc[c].vertex[i] = NONE;
      bc[c].vertexValue[i] = 0.;
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFBCSetFromOptions"
/*
  VFBCSetFromOptions: Get boundary condition flag for each face, edge, vertex of the domain.
  The option names are -<prefix>_<name> where name is the name of a geometric entity

  (c) 2010-2014 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFBCSetFromOptions(VFBC *bc,const char prefix[],PetscInt dof)
{
  PetscErrorCode ierr;
  PetscInt       i,c;
  char           optstr[256],optstr2[256],helpstr[256],helpstr2[256];

  PetscFunctionBegin;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"\n\nVFCartFE: Boundary Conditions:",__FILE__);CHKERRQ(ierr);
  /*
    faces 
  */
  for (i = 0; i < 6; i++) {
    for (c = 0; c < dof; c++) {
      if (dof == 1) {
        ierr = PetscSNPrintf(optstr,256,"-%s_%s_BC",prefix,FACE_NAME[i]);
        ierr = PetscSNPrintf(helpstr,256,"\t%s face %s boundary condition type",prefix,FACE_NAME[i]);
        ierr = PetscSNPrintf(optstr2,256,"-%s_%s",prefix,FACE_NAME[i]);
        ierr = PetscSNPrintf(helpstr2,256,"\t%s face %s boundary value",prefix,FACE_NAME[i]);
      } else {
        ierr = PetscSNPrintf(optstr,256,"-%s_%s_BC_%i",prefix,FACE_NAME[i],c);
        ierr = PetscSNPrintf(helpstr,256,"\t%s face %s component %i boundary condition type ",prefix,FACE_NAME[i],c);
        ierr = PetscSNPrintf(optstr2,256,"-%s_%s_%i",prefix,FACE_NAME[i],c);
        ierr = PetscSNPrintf(helpstr2,256,"\t%s face %s component %i boundary value ",prefix,FACE_NAME[i],c);
      }
      ierr = PetscOptionsEnum(optstr,helpstr,"",BCTYPE_NAME,(PetscEnum)bc[c].face[i],(PetscEnum*)&bc[c].face[i],NULL);CHKERRQ(ierr);
      ierr = PetscOptionsReal(optstr2,helpstr2,"",bc[c].faceValue[i],&bc[c].faceValue[i],NULL);CHKERRQ(ierr);
    }
  }
  /*
    edges
  */
  /* 
    Edge BC handling is incomplete, so I am removing the parsing for now
  
  for (i = 0; i < 12; i++) {
    for (c = 0; c < dof; c++) {
      if (dof == 1) {
        ierr = PetscSNPrintf(optstr,256,"-%s_%s_BC",prefix,EDGE_NAME[i]);
        ierr = PetscSNPrintf(helpstr,256,"\tedge %s boundary condition type",EDGE_NAME[i]);
      } else {
        ierr = PetscSNPrintf(optstr,256,"-%s_%s_BC_%i",prefix,EDGE_NAME[i],c);
        ierr = PetscSNPrintf(helpstr,256,"\tedge %s component %i boundary condition type ",EDGE_NAME[i],c);
      }
      ierr = PetscOptionsEnum(optstr,helpstr,"",BCTYPE_NAME,(PetscEnum)bc[c].edge[i],(PetscEnum*)&bc[c].edge[i],NULL);CHKERRQ(ierr);
    }
  }
  */
  /*
    vertex
  */
  for (i = 0; i < 8; i++) {
    for (c = 0; c < dof; c++) {
      if (dof == 1) {
        ierr = PetscSNPrintf(optstr,256,"-%s_%s_BC",prefix,VERTEX_NAME[i]);
        ierr = PetscSNPrintf(helpstr,256,"\t%s vertex %s boundary condition type",prefix,VERTEX_NAME[i]);
        ierr = PetscSNPrintf(optstr2,256,"-%s_%s",prefix,VERTEX_NAME[i]);
        ierr = PetscSNPrintf(helpstr2,256,"\t%s vertex %s boundary value",prefix,VERTEX_NAME[i]);
      } else {
        ierr = PetscSNPrintf(optstr,256,"-%s_%s_BC_%i",prefix,VERTEX_NAME[i],c);
        ierr = PetscSNPrintf(helpstr,256,"\t%s vertex %s component %i boundary condition type ",prefix,VERTEX_NAME[i],c);
        ierr = PetscSNPrintf(optstr2,256,"-%s_%s_%i",prefix,VERTEX_NAME[i],c);
        ierr = PetscSNPrintf(helpstr2,256,"\t%s vertex %s component %i boundary value",prefix,VERTEX_NAME[i],c);
      }
      ierr = PetscOptionsEnum(optstr,helpstr,"",BCTYPE_NAME,(PetscEnum)bc[c].vertex[i],(PetscEnum*)&bc[c].vertex[i],NULL);CHKERRQ(ierr);
    }
  }
  ierr = PetscOptionsEnd();
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFBCView"
/*
  VFBCView

  (c) 2010-2014 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFBCView(VFBC *bc,PetscViewer viewer,PetscInt dof)
{
  PetscErrorCode ierr;
  PetscInt       i,c;
  
  PetscFunctionBegin;
  /*
    faces 
  */
  ierr = PetscViewerASCIIPrintf(viewer,"faces:\n");CHKERRQ(ierr);
  for (i = 0; i < 6; i++) {
    for (c = 0; c < dof; c++) {
      ierr = PetscViewerASCIIPrintf(viewer,"    %s[%i]=%s:%f ",FACE_NAME[i],c,BCTYPE_NAME[bc[c].face[i]],bc[c].faceValue[i]);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
  }
  /*
    edges 
  */
  /*
  ierr = PetscViewerASCIIPrintf(viewer,"edges:\n");CHKERRQ(ierr);
  for (i = 0; i < 12; i++) {
    for (c = 0; c < dof; c++) {
      ierr = PetscViewerASCIIPrintf(viewer,"  %s[%i]=%s ",EDGE_NAME[i],c,BCTYPE_NAME[ bc[c].edge[i] ]);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
  }
  */
  /*
    vertices 
  */
  ierr = PetscViewerASCIIPrintf(viewer,"vertices:\n");CHKERRQ(ierr);
  for (i = 0; i < 8; i++) {
    for (c = 0; c < dof; c++) {
      ierr = PetscViewerASCIIPrintf(viewer,"%s[%i]=%s:%f ",VERTEX_NAME[i],c,BCTYPE_NAME[bc[c].vertex[i]],bc[c].vertexValue[i]);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecApplyDirichletFlowBC"
/*
  VecApplyDirichletBC
  (c) 2011 K. Yoshioka, CHEVRON ETC
*/
extern PetscErrorCode VecApplyDirichletFlowBC(Vec RHS,Vec BCF,VFBC *BC,PetscReal *BCflow)
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
  
  ierr = DMDAGetInfo(da,&dim,&nx,&ny,&nz,NULL,NULL,NULL,
                    &dof,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
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
      if (BC[c].face[X0] == FIXED) {
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
      if (BC[c].face[X1] == FIXED) {
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
      if (BC[c].face[Y0] == FIXED) {
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
      if (BC[c].face[Y1] == FIXED) {
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
        if (BC[c].face[Z0] == FIXED) {
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
        if (BC[c].face[Z1] == FIXED) {
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

