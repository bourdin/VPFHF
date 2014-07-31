#ifndef VFVFCartFEH
#define VFVFCartFEH
/*
  Implement a simple P1 Lagrange element on a rectangle aligned with axis directions.
  (c) 2010-2012 B. Bourdin, bourdin@lsu.edu
  
  Uses 3 points / 5th order quadrature from 
    A. Ern  and J.-L. Guermond "Theory and Practice of Finite Elements"
    Table 8.1 p. 359
*/

typedef struct {
  PetscInt     dim;                  /* dimension of the space */
  PetscInt     ng;                   /* number of integration points */
  PetscInt     nphix;                /* number of basis functions along the x axis */
  PetscInt     nphiy;                /* number of basis functions along the y axis */
  PetscInt     nphiz;                /* number of basis functions along the z axis */
  PetscReal    lx;                   /* length of the element */     
  PetscReal    ly;                   /* length of the element */     
  PetscReal    lz;                   /* length of the element */     
  /*PetscReal    *weight;*/              /* integration weight */
  PetscReal    weight[3];            /* integration weight */
  PetscReal    phi[1][1][2][3];      /* phi[i][g] = value of i^th basis function at g^th integration point */
  PetscReal    dphi[1][1][2][3];     /* dphi[i][g] = value of the derivative of the i^th basis function at g^th integration point */
} VFCartFEElement1D;

typedef struct {
  PetscInt     dim;                  /* dimension of the space */
  PetscInt     ng;                   /* number of integration points */
  PetscInt     nphix;                /* number of basis functions along the x axis */
  PetscInt     nphiy;                /* number of basis functions along the y axis */
  PetscInt     nphiz;                /* number of basis functions along the z axis */
  PetscReal    lx;                   /* length of the element */     
  PetscReal    ly;                   /* length of the element */     
  PetscReal    lz;                   /* length of the element */     
  /*PetscReal    *weight;*/              /* integration weight */
  PetscReal    weight[9];            /* integration weight */
  PetscReal    phi[1][2][2][9];      /* phi[j][i][g] = value of (i,j)^th basis function at g^th integration point */
  PetscReal    dphi[1][2][2][2][9];  /* phi[j][i][l][g] = value of the derivative w.r.t. x_l of the (i,j)^th basis function at g^th integration point */
} VFCartFEElement2D;

typedef struct {
  PetscInt     dim;                  /* dimension of the space */
  PetscInt     ng;                   /* number of integration points */
  PetscInt     nphix;                /* number of basis functions along the x axis */
  PetscInt     nphiy;                /* number of basis functions along the y axis */
  PetscInt     nphiz;                /* number of basis functions along the z axis */
  PetscReal    lx;                   /* length of the element */     
  PetscReal    ly;                   /* length of the element */     
  PetscReal    lz;                   /* length of the element */     
  /*PetscReal    *weight;*/              /* integration weight */
  PetscReal    weight[27];           /* integration weight */
  PetscReal    phi[2][2][2][27];     /* phi[k][j][i][g] = value of (i,j,k)^th basis function at g^th integration point */
  PetscReal    dphi[2][2][2][3][27]; /* phi[k][j][i][l][g] = value of the derivative w.r.t. x_l of the (i,j,k)^th basis function at g^th integration point */
} VFCartFEElement3D;

typedef struct {
  PetscInt     dim;                  /* dimension of the space */
  PetscInt     ng;                   /* number of integration points */
  PetscInt     nphix;                /* number of basis functions along the x axis */
  PetscInt     nphiy;                /* number of basis functions along the y axis */
  PetscInt     nphiz;                /* number of basis functions along the z axis */
  PetscReal    lx;                   /* length of the element */     
  PetscReal    ly;                   /* length of the element */     
  PetscReal    lz;                   /* length of the element */     
  PetscReal    *weight;              /* integration weight */
  PetscReal    ****phi;              /* phi[k][j][i][g] = value of (i,j,k)^th basis function at g^th integration point */
  PetscReal    *****dphi;            /* phi[k][j][i][l][g] = value of the derivative w.r.t. x_l of the (i,j,k)^th basis function at g^th integration point */
} CartesianElement;

typedef struct {
  PetscReal    x;
  PetscReal    y;
  PetscReal    z;
} coord3d;

typedef enum {
  NONE,
  ZERO,
  ONE,
  FIXED,
} BCTYPE;

typedef enum {
  X0,
  X1,
  Y0,
  Y1,
  Z0,
  Z1
} FACE;

typedef enum {
  X0Z0,
  X1Z0,
  Y0Z0,
  Y1Z0,
  X0Z1,
  X1Z1,
  Y0Z1,
  Y1Z1,
  X0Y0,
  X0Y1,
  X1Y0,
  X1Y1
} EDGE;

typedef enum {
  X0Y0Z0,
  X1Y0Z0,
  X0Y1Z0,
  X1Y1Z0,
  X0Y0Z1,
  X1Y0Z1,
  X0Y1Z1,
  X1Y1Z1
} VERTEX;

typedef struct {
  BCTYPE        face[6];
  PetscReal     faceValue[6];
  BCTYPE        edge[12];
  PetscReal     edgeValue[12];
  BCTYPE        vertex[8];   
  PetscReal     vertexValue[8];   
} VFBC;

#ifndef VFCartFEC
extern PetscErrorCode VFCartFEInit();
extern PetscErrorCode VFCartFEElement1DCreate1D(VFCartFEElement1D *e);
extern PetscErrorCode VFCartFEElement1DInit(VFCartFEElement1D *e,PetscReal lx);
extern PetscErrorCode VFCartFEElement2DCreate(VFCartFEElement2D *e);
extern PetscErrorCode VFCartFEElement2DInit(VFCartFEElement2D *e,PetscReal lx,PetscReal ly);
extern PetscErrorCode VFCartFEElement3DCreate(VFCartFEElement3D *e);
extern PetscErrorCode VFCartFEElement3DInit(VFCartFEElement3D *e,PetscReal lx,PetscReal ly,PetscReal lz);

extern PetscErrorCode DAReadCoordinatesHDF5(DM da,const char filename[]);

extern PetscErrorCode VFBCCreate(VFBC *bc,PetscInt dof);
extern PetscErrorCode VFBCSetFromOptions(VFBC *bc,const char prefix[],PetscInt dof);
extern PetscErrorCode VFBCView(VFBC *bc,PetscViewer viewer,PetscInt dof);
extern PetscErrorCode VecSetFromBC(Vec BCVec,VFBC *BC);

extern PetscErrorCode VecApplyDirichletBC(Vec RHS,Vec BCU,VFBC *BC);
extern PetscErrorCode ResidualApplyDirichletBC(Vec Residual,Vec U,Vec BCU,VFBC *BC);
extern PetscErrorCode GradientApplyDirichletBC(Vec gradient,VFBC *BC);
extern PetscErrorCode MatApplyDirichletBC(Mat K,VFBC *BC);
extern PetscErrorCode MatApplyDirichletBCRowCol(Mat K,VFBC *BC);
extern PetscErrorCode VecApplyDirichletFlowBC(Vec RHS,Vec BCU,VFBC *BC,PetscReal *BCpres);

#endif
#endif /* VFCartFEH */

