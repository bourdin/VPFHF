#include "petsc.h"
#include "CartFE.h"
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
  PennyCrack->r        = 0.;
  PennyCrack->theta    = 0.;
  PennyCrack->phi      = 0.;
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
  ierr = PetscViewerASCIIPrintf(viewer,"center:   \t%e \t%e \t%e\n",PennyCrack->center[0],PennyCrack->center[1],PennyCrack->center[2]);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"radius:     \t%e\n",PennyCrack->r);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"polar angle:\t%e\n",PennyCrack->theta);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"co-latitude:\t%e\n",PennyCrack->phi);CHKERRQ(ierr);
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
  PetscErrorCode      ierr;
  PetscReal           n[3],tau[3],l,xdotn, taumag;

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
  xdotn  = (x[0]-PennyCrack->center[0])*n[0] + 
           (x[1]-PennyCrack->center[1])*n[1] + 
           (x[2]-PennyCrack->center[2])*n[2];
  tau[0] = (x[0]-PennyCrack->center[0]) - xdotn*n[0];
  tau[1] = (x[1]-PennyCrack->center[1]) - xdotn*n[1];
  tau[2] = (x[2]-PennyCrack->center[2]) - xdotn*n[2];
  l = sqrt(tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2]);
  if (PennyCrack->r == 0.) {
	  *d = 0.;
  } else {
    if (l <= PennyCrack->r) {
      *d = sqrt(xdotn*xdotn);
    } else {		
	  *d = sqrt(pow(xdotn,2)+pow((l-PennyCrack->r),2));
    } 
  }
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
  PetscErrorCode      ierr;
  PetscInt            i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal       ****coords_array;
  PetscReal        ***v_array;
  PetscReal           x[3];
  PetscReal           dist;  

  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);

	ierr = VecSet(V,1.0);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  for (k = zs; k < zs+zm; k++) {
    for (j = ys; j < ys+ym; j++) {
      for (i = xs; i < xs+xm; i++) { 
        x[2] = coords_array[k][j][i][2];
        x[1] = coords_array[k][j][i][1];
        x[0] = coords_array[k][j][i][0];
        ierr = VFDistanceToPennyCrack(&dist,x,crack);CHKERRQ(ierr);
		  if(crack->r == 0){
			  v_array[k][j][i] = 1.;
		  }
		  if(dist == 0){
			  v_array[k][j][i] = 0.;
		  }
		  else{
			  v_array[k][j][i] = 1.-exp(-dist/2/ctx->vfprop.epsilon);
		  }
      }
    }
  }      
	ierr = DMDAVecRestoreArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
//
// #undef __FUNCT__
// #define __FUNCT__ "VFRectangularCrackGet"
// /*
//   
// 
//   VFRectangularCrackGet (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
// */
// extern PetscErrorCode VFRectangularCrackGet(const char prefix[],VFRectangularCrack *RectangularCrack)
// {
//   PetscErrorCode ierr;
//   PetscInt       nval=3;
//   PetscFunctionBegin;
//   ierr = PetscOptionsBegin(PETSC_COMM_WORLD,prefix,"\n\nVF: Rectangular-shaped crack description:","");CHKERRQ(ierr);
//   {
//     ierr = PetscOptionsString("-name","\n\tRectangular-shaped crack name","",RectangularCrack->name,RectangularCrack->name,sizeof(RectangularCrack->name),PETSC_NULL);CHKERRQ(ierr);
//     nval = 9;    
//     ierr = PetscOptionsRealArray("-corners","\n\tRectangular-shaped crack corners coordinates (x0,y0,z0, x1,y1,z1, x2,y2,z2)  (comma separated).","",RectangularCrack->corner,&nval,PETSC_NULL);CHKERRQ(ierr);
//   }
//   ierr = PetscOptionsEnd();CHKERRQ(ierr);
//   PetscFunctionReturn(0);
// }
// 
// #undef __FUNCT__
// #define __FUNCT__ "VFRectangularCrackCreate"
// /*
//   VFRectangularCrackCreate: Allocates a RectangularCrack data structure
//   
//   (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
// */
// extern PetscErrorCode VFRectangularCrackCreate(VFRectangularCrack *RectangularCrack)
// {
//   PetscErrorCode ierr;
//   int            i;
//   
//   PetscFunctionBegin;
//   ierr = PetscStrcpy(RectangularCrack->name,"RectangularCrack");CHKERRQ(ierr);
//   for (i=0; i<9; i++) RectangularCrack->corners[i] = 0.;
//   PetscFunctionReturn(0);
// }
// 
// #undef __FUNCT__
// #define __FUNCT__ "VFRectangularCrackView"
// /*
//   VFRectangularCrackView 
//   
//   (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
// */
// extern PetscErrorCode VFRectangularCrackView(VFRectangularCrack *RectangularCrack,PetscViewer viewer)
// {
//   PetscErrorCode ierr;
//   
//   PetscFunctionBegin;
//   ierr = PetscViewerASCIIPrintf(viewer,"RectangularCrack object \"%s\":\n",RectangularCrack->name);CHKERRQ(ierr);
//   ierr = PetscViewerASCIIPrintf(viewer,"corners:     \t%e \t%e \t%e \t%e\n\t\t%e \t%e \t%e \t%e\n\t\t%e \t%e \t%e \t%e\n",
//                                RectangularCrack->corners[0],RectangularCrack->corners[1],RectangularCrack->corners[2],
//                                RectangularCrack->corners[3],RectangularCrack->corners[4],RectangularCrack->corners[5],
//                                RectangularCrack->corners[6],RectangularCrack->corners[7],RectangularCrack->corners[8]);CHKERRQ(ierr);
//   PetscFunctionReturn(0);
// }
// 
// #undef __FUNCT__
// #define __FUNCT__ "VFRectangularCrackSetName"
// /*
//   VFRectangularCrackSetName 
//   
//   (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
// */
// extern PetscErrorCode VFRectangularCrackSetName(VFRectangularCrack *RectangularCrack,const char name[])
// {
//   PetscErrorCode ierr;
//   
//   PetscFunctionBegin;
//   ierr = PetscStrcpy(RectangularCrack->name,name);CHKERRQ(ierr);
//   PetscFunctionReturn(0);
// }
// 
// #undef __FUNCT__
// #define __FUNCT__ "VFDistanceToRectangularCrack"
// /*
//   VFDistanceToRectangularCrack: Computes the distance between a point with coordinates x and a RectangularCrack described by RectangularCrack
//   
//   
//   (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
// */
// extern PetscErrorCode VFDistanceToRectangularCrack(PetscReal *d,PetscReal *x,VFRectangularCrack *Crack)
// {
//   PetscErrorCode      ierr;
//   PetscReal           n[3],tau[3],l,xdotn;
// 
//   PetscFunctionBegin;
//   /*
//     n: normal vector to the disk
//   */
//   n[0] = cos(Crack->theta * PETSC_PI/180.)*sin(Crack->phi * PETSC_PI/180.);
//   n[1] = sin(Crack->theta * PETSC_PI/180.)*sin(Crack->phi * PETSC_PI/180.);
//   n[2] = cos(Crack->phi * PETSC_PI/180.);
//   
//   /*
//     tau: projection onto the disk plane
//   */
//   xdotn  = (x[0]-PennyCrack->corner[0])*n[0] + 
//            (x[1]-PennyCrack->corner[1])*n[1] + 
//            (x[2]-PennyCrack->corner[2])*n[2];
//   tau[0] = (x[0]-PennyCrack->corner[0]) - xdotn*n[0];
//   tau[1] = (x[1]-PennyCrack->corner[1]) - xdotn*n[1];
//   tau[2] = (x[2]-PennyCrack->corner[2]) - xdotn*n[2];
//   l = sqrt(tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2]);
//   if (PennyCrack->r == 0.) {
//     *d = sqrt(pow(x[0] - PennyCrack->center[0],2) +  
//               pow(x[1] - PennyCrack->center[1],2) + 
//               pow(x[2] - PennyCrack->center[2],2));
//   } else {
//     if (l < PennyCrack->r) {
//       *d = sqrt(xdotn*xdotn);
//     } else {
//       *d = sqrt(pow(x[0] - PennyCrack->center[0] - tau[0] / l * PennyCrack->r,2) +  
//                 pow(x[1] - PennyCrack->center[1] - tau[1] / l * PennyCrack->r,2) + 
//                 pow(x[2] - PennyCrack->center[2] - tau[2] / l * PennyCrack->r,2));
//     } 
//   }
//   PetscFunctionReturn(0);
// }

