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
  //ierr = PetscMalloc(3*sizeof(PetscReal),&(PennyCrack->center));CHKERRQ(ierr);
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
  PetscReal           n[3],tau[3],l,xdotn;

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
    *d = sqrt(pow(x[0] - PennyCrack->center[0],2) +  
              pow(x[1] - PennyCrack->center[1],2) + 
              pow(x[2] - PennyCrack->center[2],2));
  } else {
    if (l < PennyCrack->r) {
      *d = sqrt(xdotn*xdotn);
    } else {
      *d = sqrt(pow(x[0] - PennyCrack->center[0] - tau[0] / l * PennyCrack->r,2) +  
                pow(x[1] - PennyCrack->center[1] - tau[1] / l * PennyCrack->r,2) + 
                pow(x[2] - PennyCrack->center[2] - tau[2] / l * PennyCrack->r,2));
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
  ierr = DAGetInfo(ctx->daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DAGetCorners(ctx->daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

	ierr = DAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);

	ierr = VecSet(V,1.0);CHKERRQ(ierr);
	ierr = DAVecGetArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  for (k = zs; k < zs+zm; k++) {
    for (j = ys; j < ys+ym; j++) {
      for (i = xs; i < xs+xm; i++) { 
        x[2] = coords_array[k][j][i][2];
        x[1] = coords_array[k][j][i][1];
        x[0] = coords_array[k][j][i][0];
        ierr = VFDistanceToPennyCrack(&dist,x,crack);CHKERRQ(ierr);
        v_array[k][j][i] = 1.-exp(-dist/2/ctx->vfprop.epsilon);
      }
    }
  }      
	ierr = DAVecRestoreArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  ierr = DAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}