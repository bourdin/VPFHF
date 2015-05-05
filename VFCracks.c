#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFCracks.h"


#undef __FUNCT__
#define __FUNCT__ "VFPennyCrackGet"
/*


 VFPennyCrackGet (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFPennyCrackGet(const char prefix[],VFPennyCrack *PennyCrack)
{
  PetscErrorCode ierr;
  PetscInt       nval=3;
  PetscFunctionBegin;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,prefix,"\n\nVF: Penny-shaped crack description:","");CHKERRQ(ierr);
  {
    ierr = PetscOptionsString("-name","\n\tPenny-shaped crack name","",PennyCrack->name,PennyCrack->name,sizeof(PennyCrack->name),PETSC_NULL);CHKERRQ(ierr);
    nval = 3;

    ierr = PetscOptionsRealArray("-center","\n\tPenny-shaped crack center of center  (comma separated).","",PennyCrack->center,&nval,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-r","\n\t Penny-shaped crack radius","",PennyCrack->r,&PennyCrack->r,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-theta","\n\t Penny-shaped crack polar angle (in degrees)","",PennyCrack->theta,&PennyCrack->theta,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-phi","\n\t Penny-shaped crack co-latitude (in degrees)","",PennyCrack->phi,&PennyCrack->phi,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-thickness","\n\t Penny-shaped crack thickness","",PennyCrack->thickness,&PennyCrack->thickness,PETSC_NULL);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFPennyCrackCreate"
/*
 VFPennyCrackCreate: Allocates a PennyCrack data structure

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFPennyCrackCreate(VFPennyCrack *PennyCrack)
{
  PetscErrorCode ierr;
  int            i;

  PetscFunctionBegin;
  ierr = PetscStrcpy(PennyCrack->name,"PennyCrack");CHKERRQ(ierr);
  for (i=0; i<3; i++) PennyCrack->center[i] = 0.;
  PennyCrack->r         = 0.;
  PennyCrack->theta     = 0.;
  PennyCrack->phi       = 0.;
  PennyCrack->thickness = 0.;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFPennyCrackView"
/*
 VFPennyCrackView

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFPennyCrackView(VFPennyCrack *PennyCrack,PetscViewer viewer)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(viewer,"PennyCrack object \"%s\":\n",PennyCrack->name);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"center:     \t%e \t%e \t%e\n",PennyCrack->center[0],PennyCrack->center[1],PennyCrack->center[2]);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"radius:     \t%e\n",PennyCrack->r);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"polar angle:\t%e\n",PennyCrack->theta);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"co-latitude:\t%e\n",PennyCrack->phi);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"thickness:\t%e\n",PennyCrack->thickness);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFPennyCrackSetName"
/*
 VFPennyCrackSetName

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFPennyCrackSetName(VFPennyCrack *PennyCrack,const char name[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscStrcpy(PennyCrack->name,name);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFDistanceToPennyCrack"
/*
 VFDistanceToPennyCrack: Computes the distance between a point with coordinates x and a PennyCrack described by PennyCrack


 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFDistanceToPennyCrack(PetscReal *d,PetscReal *x,VFPennyCrack *PennyCrack)
{
  PetscReal n[3],tau[3],l,xdotn;

  PetscFunctionBegin;
  /*
   n: normal vector to the disk
   */
  n[0] = cos(PennyCrack->theta * PETSC_PI/180.)*sin(PennyCrack->phi * PETSC_PI/180.);
  n[1] = sin(PennyCrack->theta * PETSC_PI/180.)*sin(PennyCrack->phi * PETSC_PI/180.);
  n[2] = cos(PennyCrack->phi * PETSC_PI/180.);
  /*
   tau: projection onto the disk plane
   */
  xdotn = (x[0]-PennyCrack->center[0])*n[0] +
          (x[1]-PennyCrack->center[1])*n[1] +
          (x[2]-PennyCrack->center[2])*n[2];
  tau[0] = (x[0]-PennyCrack->center[0]) - xdotn*n[0];
  tau[1] = (x[1]-PennyCrack->center[1]) - xdotn*n[1];
  tau[2] = (x[2]-PennyCrack->center[2]) - xdotn*n[2];
  l      = PetscSqrtReal(tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2]);

  if (l <= PennyCrack->r) *d = PetscSqrtReal(xdotn*xdotn);
  else *d = PetscSqrtReal(pow(xdotn,2)+pow((l-PennyCrack->r),2));
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFPennyCrackBuildVAT2"
/*
 VFPennyCrackBuildVAT2:  Build the V-field associated with the array of penny-shaped cracks
 following the construction in Bourdin-Francfort-Marigo '08.

 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFPennyCrackBuildVAT2(Vec V,VFPennyCrack *crack,VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal      ****coords_array;
  PetscReal      ***v_array;
  PetscReal      x[3];
  PetscReal      dist;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);

  ierr = VecSet(V,1.0);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  for (k = zs; k < zs+zm; k++) {
    for (j = ys; j < ys+ym; j++)
      for (i = xs; i < xs+xm; i++) {
        x[2] = coords_array[k][j][i][2];
        x[1] = coords_array[k][j][i][1];
        x[0] = coords_array[k][j][i][0];
        if (crack->r == 0) v_array[k][j][i] = 1.;
        else {
          ierr = VFDistanceToPennyCrack(&dist,x,crack);CHKERRQ(ierr);
          dist -= crack->thickness/2.;
          if (dist <= 0) {
            v_array[k][j][i] = 0.;
          } else {
            v_array[k][j][i] = 1.-exp(-dist/2/ctx->vfprop.epsilon);
          }
        }
      }
  }
  ierr = DMDAVecRestoreArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFPennyCrackBuildVAT1"
/*
 VFPennyCrackBuildVAT1:  Build the V-field associated with the array of penny-shaped cracks
 following the construction in Bourdin-Francfort-Marigo '08.

 (c) 2013 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFPennyCrackBuildVAT1(Vec V,VFPennyCrack *crack,VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal      ****coords_array;
  PetscReal      ***v_array;
  PetscReal      x[3];
  PetscReal      dist;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);

  ierr = VecSet(V,1.0);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  for (k = zs; k < zs+zm; k++) {
    for (j = ys; j < ys+ym; j++)
      for (i = xs; i < xs+xm; i++) {
        x[2] = coords_array[k][j][i][2];
        x[1] = coords_array[k][j][i][1];
        x[0] = coords_array[k][j][i][0];
        if (crack->r == 0) {
          v_array[k][j][i] = 1.;
        } else {
          ierr = VFDistanceToPennyCrack(&dist,x,crack);CHKERRQ(ierr);
          dist -= crack->thickness/2.;
          if (dist <= 0) {
            v_array[k][j][i] = 0.;
          } else if (dist < 2. * ctx->vfprop.epsilon) {
            v_array[k][j][i] = dist/ctx->vfprop.epsilon * (1.- .25*dist/ctx->vfprop.epsilon);
          } else {
            v_array[k][j][i] = 1.;
          }
        }
      }
  }
  ierr = DMDAVecRestoreArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFRectangularCrackGet"
/*
 VFRectangularCrackGet (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFRectangularCrackGet(const char prefix[],VFRectangularCrack *RectangularCrack)
{
  PetscErrorCode ierr;
  PetscInt       nval=3;
  PetscFunctionBegin;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,prefix,"\n\nVF: Rectangular-shaped crack description:","");CHKERRQ(ierr);
  {
    ierr = PetscOptionsString("-name","\n\tRectangular-shaped crack name","",RectangularCrack->name,RectangularCrack->name,sizeof(RectangularCrack->name),PETSC_NULL);CHKERRQ(ierr);
    nval = 9;
    ierr = PetscOptionsRealArray("-corners","\n\tRectangular-shaped crack corners coordinates (x0,y0,z0, x1,y1,z1, x2,y2,z2)  (comma separated).","",RectangularCrack->corners,&nval,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-thickness","\n\tRectangular-shaped crack thickness","",RectangularCrack->thickness,&RectangularCrack->thickness,PETSC_NULL);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "VFRectangularCrackCreate"
/*
 VFRectangularCrackCreate: Allocates a RectangularCrack data structure
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFRectangularCrackCreate(VFRectangularCrack *RectangularCrack)
{
  PetscErrorCode ierr;
  int            i;

  PetscFunctionBegin;
  ierr = PetscStrcpy(RectangularCrack->name,"RectangularCrack");CHKERRQ(ierr);
  for (i=0; i<9; i++) RectangularCrack->corners[i] = 0.;
  RectangularCrack->thickness = 0.;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFRectangularCrackView"
/*
 VFRectangularCrackView
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFRectangularCrackView(VFRectangularCrack *RectangularCrack,PetscViewer viewer)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(viewer,"RectangularCrack object \"%s\":\n",RectangularCrack->name);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"corners:     \t%e \t%e \t%e\n\t\t%e \t%e \t%e\n\t\t%e \t%e \t%e\n",
                                RectangularCrack->corners[0],RectangularCrack->corners[1],RectangularCrack->corners[2],
                                RectangularCrack->corners[3],RectangularCrack->corners[4],RectangularCrack->corners[5],
                                RectangularCrack->corners[6],RectangularCrack->corners[7],RectangularCrack->corners[8]);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"thickness:\t%e\n",RectangularCrack->thickness);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFRectangularCrackSetName"
/*
 VFRectangularCrackSetName
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFRectangularCrackSetName(VFRectangularCrack *RectangularCrack,const char name[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscStrcpy(RectangularCrack->name,name);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFDistanceToRectangularCrack"
/*
 VFDistanceToRectangularCrack: Computes the distance between a point with coordinates x and a RectangularCrack described by RectangularCrack
 (c) 2010-2012 Chukwudi Chukwudozie cchukw1@tigers.lsu.edu
 */
extern PetscErrorCode VFDistanceToRectangularCrack(PetscReal *d,PetscReal *x,VFRectangularCrack *Crack)
{
  PetscReal n_1[3],n_2[3],n[3],n_q[3],l,xdotn;
  PetscReal a[3],b[3],center[3],corner4[3],n_norm;
  PetscReal length_to_side;
  PetscReal p[3],q[3],tau[3];
  PetscReal dist_a,dist_b,theta_b, theta, theta_side, theta_top;
  PetscFunctionBegin;
  /*
   Corner allocation:

   pt3;[6,7,8]           pt4;corner4
          +-----------------------+
          |         Side_D        |
          |                       |
          |                       |
          |         center        |
   Side_A |           o           |    Side_B
          |                       |
          |                       |
          |                       |
          |         Side_C        |
          +-----------------------+
   pt1;[0,1,2]           pt2;[3,4,5]
   */
  center[0] = (Crack->corners[6]+Crack->corners[3])/2.;
  center[1] = (Crack->corners[7]+Crack->corners[4])/2.;
  center[2] = (Crack->corners[8]+Crack->corners[5])/2.;

  corner4[0] = 2.*center[0]-Crack->corners[0];
  corner4[1] = 2.*center[1]-Crack->corners[1];
  corner4[2] = 2.*center[2]-Crack->corners[2];
  /*
   n: normal vector to the plane
   */
  a[0] = Crack->corners[3]-Crack->corners[0];
  a[1] = Crack->corners[4]-Crack->corners[1];
  a[2] = Crack->corners[5]-Crack->corners[2];

  b[0] = Crack->corners[6]-Crack->corners[0];
  b[1] = Crack->corners[7]-Crack->corners[1];
  b[2] = Crack->corners[8]-Crack->corners[2];

  n[0]   = a[1]*b[2] - a[2]*b[1];
  n[1]   = a[2]*b[0] - a[0]*b[2];
  n[2]   = a[0]*b[1] - a[1]*b[0];
  n_norm = PetscSqrtReal(pow(n[0],2)+pow(n[1],2)+pow(n[2],2));
  n[0]   = n[0]/n_norm;n[1] = n[1]/n_norm;n[2] = n[2]/n_norm;
  /*
   tau: projection onto the rectangular plane
   */
  xdotn = (x[0]-center[0])*n[0] +
          (x[1]-center[1])*n[1] +
          (x[2]-center[2])*n[2];

  tau[0] = (x[0]-center[0]) - xdotn*n[0];
  tau[1] = (x[1]-center[1]) - xdotn*n[1];
  tau[2] = (x[2]-center[2]) - xdotn*n[2];

  p[0] = center[0]+xdotn*n[0];p[1] = center[1]+xdotn*n[1];p[2] = center[2]+xdotn*n[2];
  q[0] = x[0]+center[0]-p[0];q[1] = x[1]+center[1]-p[1];q[2] = x[2]+center[2]-p[2];
  /**Check which side of polygon point is*/
  /*Vector in the direction of center and corner3*/
  n_1[0] = Crack->corners[6]-center[0];
  n_1[1] = Crack->corners[7]-center[1];
  n_1[2] = Crack->corners[8]-center[2];
  n_norm = PetscSqrtReal(pow(n_1[0],2)+pow(n_1[1],2)+pow(n_1[2],2));
  n_1[0] = n_1[0]/n_norm;n_1[1] = n_1[1]/n_norm;n_1[2] = n_1[2]/n_norm;
  /*Vector in the direction of center and corner4*/
  n_2[0] = corner4[0]-center[0];
  n_2[1] = corner4[1]-center[1];
  n_2[2] = corner4[2]-center[2];
  n_norm = PetscSqrtReal(pow(n_2[0],2)+pow(n_2[1],2)+pow(n_2[2],2));
  n_2[0] = n_2[0]/n_norm;n_2[1] = n_2[1]/n_norm;n_2[2] = n_2[2]/n_norm;

  theta_top  = acos(n_2[0]*n_1[0]+n_2[1]*n_1[1]+n_2[2]*n_1[2]);
  theta_side = PETSC_PI-theta_top;
  /*Vector in the direction of q*/
  n_q[0] = q[0]-center[0];
  n_q[1] = q[1]-center[1];
  n_q[2] = q[2]-center[2];
  n_norm = PetscSqrtReal(pow(n_q[0],2)+pow(n_q[1],2)+pow(n_q[2],2));
  n_q[0] = n_q[0]/n_norm;n_q[1] = n_q[1]/n_norm;n_q[2] = n_q[2]/n_norm;
  theta  = acos(n_2[0]*n_q[0]+n_2[1]*n_q[1]+n_2[2]*n_q[2]) + acos(-n_1[0]*n_q[0]-n_1[1]*n_q[1]-n_1[2]*n_q[2]);

  if ((theta-theta_side) < 1e-6 || (2.*PETSC_PI-theta_side-theta) < 1e-6) {
    n_1[0]         = corner4[0]-Crack->corners[3];
    n_1[1]         = corner4[1]-Crack->corners[4];
    n_1[2]         = corner4[2]-Crack->corners[5];
    n_norm         = PetscSqrtReal(pow(n_1[0],2)+pow(n_1[1],2)+pow(n_1[2],2));
    n_1[0]         = n_1[0]/n_norm;n_1[1] = n_1[1]/n_norm;n_1[2] = n_1[2]/n_norm;
    dist_a         = PetscSqrtReal(pow((corner4[0]-center[0]),2)+pow((corner4[1]-center[1]),2)+pow((corner4[2]-center[2]),2));
    dist_b         = (corner4[0]-center[0])*n_1[0]+(corner4[1]-center[1])*n_1[1]+(corner4[2]-center[2])*n_1[2];
    dist_b         = PetscAbs(dist_b);
    dist_b         = PetscSqrtReal(pow(dist_a,2)-pow(dist_b,2));
    theta_b        = acos(n_1[0]*n_q[0]+n_1[1]*n_q[1]+n_1[2]*n_q[2]);
    length_to_side = dist_b/sin(theta_b);
  }else {
    n_1[0]         = corner4[0]-Crack->corners[6];
    n_1[1]         = corner4[1]-Crack->corners[7];
    n_1[2]         = corner4[2]-Crack->corners[8];
    n_norm         = PetscSqrtReal(pow(n_1[0],2)+pow(n_1[1],2)+pow(n_1[2],2));
    n_1[0]         = n_1[0]/n_norm;n_1[1] = n_1[1]/n_norm;n_1[2] = n_1[2]/n_norm;
    dist_a         = PetscSqrtReal(pow((corner4[0]-center[0]),2)+pow((corner4[1]-center[1]),2)+pow((corner4[2]-center[2]),2));
    dist_b         = (corner4[0]-center[0])*n_1[0]+(corner4[1]-center[1])*n_1[1]+(corner4[2]-center[2])*n_1[2];
    dist_b         = PetscAbs(dist_b);
    dist_b         = PetscSqrtReal(pow(dist_a,2)-pow(dist_b,2));
    theta_b        = acos(n_1[0]*n_q[0]+n_1[1]*n_q[1]+n_1[2]*n_q[2]);
    length_to_side = dist_b/sin(theta_b);
  }
  l = PetscSqrtReal(tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2]);
  if (l > length_to_side) *d = PetscSqrtReal(pow((l-length_to_side),2) + pow(xdotn,2));
  else *d = PetscSqrtReal(xdotn*xdotn);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFRectangularCrackBuildVAT2"
/*
 VFRectangularCrackBuildVAT2:  Build the V-field associated with the array of rectangular cracks
 following the construction in Bourdin-Francfort-Marigo '08.

 (c) 2010-2012 Chukwudi Chukwuzodie cchukw1@tigers.lsu.edu
 */
extern PetscErrorCode VFRectangularCrackBuildVAT2(Vec V,VFRectangularCrack *crack,VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal      ****coords_array;
  PetscReal      ***v_array;
  PetscReal      x[3];
  PetscReal      dist;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);

  ierr = VecSet(V,1.0);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  for (k = zs; k < zs+zm; k++) {
    for (j = ys; j < ys+ym; j++)
      for (i = xs; i < xs+xm; i++) {
        x[2] = coords_array[k][j][i][2];
        x[1] = coords_array[k][j][i][1];
        x[0] = coords_array[k][j][i][0];
        ierr = VFDistanceToRectangularCrack(&dist,x,crack);CHKERRQ(ierr);
        if (dist <= crack->thickness) v_array[k][j][i] = 0.;
        else v_array[k][j][i] = 1.-exp(-dist/2/ctx->vfprop.epsilon);
      }
  }
  ierr = DMDAVecRestoreArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFRectangularCrackBuildVAT1"
/*
 VFRectangularCrackBuildVAT1:  Build the V-field associated with the array of rectangular cracks
 following the construction in Bourdin-Francfort-Marigo '08.

 (c) 2010-2013 Chukwudi Chukwuzodie cchukw1@tigers.lsu.edu
               Blaise Bourdin       bourdin@lsu.edu
 */
extern PetscErrorCode VFRectangularCrackBuildVAT1(Vec V,VFRectangularCrack *crack,VFCtx *ctx)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal      ****coords_array;
  PetscReal      ***v_array;
  PetscReal      x[3];
  PetscReal      dist;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);

  ierr = VecSet(V,1.0);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  for (k = zs; k < zs+zm; k++) {
    for (j = ys; j < ys+ym; j++)
      for (i = xs; i < xs+xm; i++) {
        x[2] = coords_array[k][j][i][2];
        x[1] = coords_array[k][j][i][1];
        x[0] = coords_array[k][j][i][0];
        ierr = VFDistanceToRectangularCrack(&dist,x,crack);CHKERRQ(ierr);
        dist -= crack->thickness/2.;
        if (dist <= 0) {
          v_array[k][j][i] = 0.;
        } else if (dist < 2. * ctx->vfprop.epsilon) {
          v_array[k][j][i] = dist/ctx->vfprop.epsilon * (1.- .25*dist/ctx->vfprop.epsilon);
        } else {
          v_array[k][j][i] = 1.;
        }
      }
  }
  ierr = DMDAVecRestoreArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

