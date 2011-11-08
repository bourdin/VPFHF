#ifndef CARTFE_H
#define CARTFE_H
/* 
  Implement a simple P1 Lagrange element on a rectangle aligned with axis directions.
  (c) 2010-2011 B. Bourdin, bourdin@lsu.edu
  
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
} CartFE_Element1D;

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
} CartFE_Element2D;

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
} CartFE_Element3D;

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
  VALUE
} BCTYPE;

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

typedef enum {
  X0,
  X1,
  Y0,
  Y1,
  Z0,
  Z1
} FACE;

static const char *FACE_NAME[] = {
	"X0","X1",
	"Y0","Y1",
	"Z0","Z1",
	"FACE_NAME","",0
};

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

static const char *EDGE_NAME[] = {
	"X0Z0","X1Z0","Y0Z0","Y1Z0",
	"X0Z1","X1Z1","Y0Z1","Y1Z1",
	"X0Y0","X0Y1","X1Y0","X1Y1",
	"EDGE_NAME","",0
};

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

static const char *VERTEX_NAME[] = {
	"X0Y0Z0","X1Y0Z0","X0Y1Z0","X1Y1Z0",
	"X0Y0Z1","X1Y0Z1","X0Y1Z1","X1Y1Z1",
	"VERTEX_NAME","",0};

typedef struct {
  BCTYPE     face[6];
  BCTYPE     edge[12];
  BCTYPE     vertex[8];   
} BC;

#ifndef CARTFE_C
extern PetscErrorCode CartFE_Init();
extern PetscErrorCode CartFE_Element1DCreate1D(CartFE_Element1D *e);
extern PetscErrorCode CartFE_Element1DInit(CartFE_Element1D *e,PetscReal lx);
extern PetscErrorCode CartFE_Element2DCreate(CartFE_Element2D *e);
extern PetscErrorCode CartFE_Element2DInit(CartFE_Element2D *e,PetscReal lx,PetscReal ly);
extern PetscErrorCode CartFE_Element3DCreate(CartFE_Element3D *e);
extern PetscErrorCode CartFE_Element3DInit(CartFE_Element3D *e,PetscReal lx,PetscReal ly,PetscReal lz);

extern PetscErrorCode DAReadCoordinatesHDF5(DA da,const char filename[]);
extern PetscErrorCode BCInit(BC *bc,PetscInt dof);
extern PetscErrorCode BCGet(BC *bc,const char prefix[],PetscInt dof);
extern PetscErrorCode BCView(BC *bc,PetscViewer viewer,PetscInt dof);

extern PetscErrorCode VecApplyDirichletBC(Vec RHS,Vec BCU,BC *BC);
extern PetscErrorCode MatApplyDirichletBC(Mat K,DA da,BC *BC);
extern PetscErrorCode VecApplyDirichletFlowBC(Vec RHS,Vec BCU,BC *BC,PetscReal *BCpres);
#endif
#endif /* CARTFE_H */

