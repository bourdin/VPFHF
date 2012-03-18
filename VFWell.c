#include "petsc.h"
#include "VFWell.h"

#undef __FUNCT__
#define __FUNCT__ "VFWellGet"
/*
  

  VFWellGet (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFWellGet(const char prefix[],VFWell *well)
{
  PetscErrorCode ierr;
  PetscInt       nval=3;
  PetscFunctionBegin;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,prefix,"\n\nVF: well description:","");CHKERRQ(ierr);
  {
    ierr = PetscOptionsString("-name","\n\twell name","",well->name,well->name,sizeof(well->name),PETSC_NULL);CHKERRQ(ierr);
    nval = 3;
    ierr = PetscOptionsRealArray("-top","\n\twell top coordinates (comma separated).","",well->top,&nval,PETSC_NULL);CHKERRQ(ierr);
    nval = 3;
    ierr = PetscOptionsRealArray("-bottom","\n\t well bottom coordinates  (comma separated).","",well->bottom,&nval,PETSC_NULL);CHKERRQ(ierr);
    /*
    ierr = PetscOptionsEnum("-bcV","\n\t boundary condition on V field","",BCTYPE_NAME,(PetscEnum)well->BCV,(PetscEnum*)&(well->BCV),PETSC_NULL);CHKERRQ(ierr);
    */
}
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFWellCreate"
/*
  VFWellCreate: Allocates a well data structure
  
  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFWellCreate(VFWell *well)
{
  PetscErrorCode ierr;
  int            i;
  
  PetscFunctionBegin;
  ierr = PetscStrcpy(well->name,"well");CHKERRQ(ierr);
  for (i=0; i<3; i++) well->top[i] = 0.;
  for (i=0; i<3; i++) well->bottom[i] = 0.;
  /*
  well->BCV = NONE;
  */
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFWellView"
/*
  VFWellView 
  
  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFWellView(VFWell *well,PetscViewer viewer)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(viewer,"Well object \"%s\":\n",well->name);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"top:    \t%e \t%e \t%e\n",well->top[0],well->top[1],well->top[2]);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"bottom: \t%e \t%e \t%e\n",well->bottom[0],well->bottom[1],well->bottom[2]);CHKERRQ(ierr);
  /*
  ierr = PetscViewerASCIIPrintf(viewer,"BCV:    \t%s\n",BCTYPE_NAME[well->BCV]);CHKERRQ(ierr);
  */
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFWellSetName"
/*
  VFWellSetName 
  
  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFWellSetName(VFWell *well,const char name[])
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscStrcpy(well->name,name);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFDistanceToWell"
/*
  VFDistanceToWell: Computes the distance between a point with coordinates x and a well described by well
  
  
  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFDistanceToWell(PetscReal *d,PetscReal *x,VFWell *well)
{
  PetscErrorCode ierr;
  PetscReal      xx0,x1x0,xx0x1x0,l;
  PetscReal      x0[3],x1[3];
  PetscInt       i;

  PetscFunctionBegin;
  for (i = 0; i < 3; i++) x0[i] = well->top[i];
  for (i = 0; i < 3; i++) x1[i] = well->bottom[i];

  /*
    l is the projection of x along (x0,x1)
    l = (x-x_0)\cdot(x_1-x_0) / \|x1-x_0\|^2
    x is between the planes normal to (x0,x1) if and only if 0 \le l \le 1
  */
  x1x0 =    (x1[0] - x0[0]) * (x1[0] - x0[0]) +
            (x1[1] - x0[1]) * (x1[1] - x0[1]) +
            (x1[2] - x0[2]) * (x1[2] - x0[2]);
  l = ( (x[0] - x0[0]) * (x1[0] - x0[0]) +
        (x[1] - x0[1]) * (x1[1] - x0[1]) +
        (x[2] - x0[2]) * (x1[2] - x0[2]) ) / x1x0;
  if (l < 0.) {
    /*
      d = \|x-x_0\|
    */
    *d = sqrt( (x[0] - x0[0]) * (x[0] - x0[0]) +
               (x[1] - x0[1]) * (x[1] - x0[1]) +
               (x[2] - x0[2]) * (x[2] - x0[2]) );
    ierr = PetscLogFlops(24);CHKERRQ(ierr);
  } else if ( l < 1.) {
    /*
      d^2 = \|x-x_0\|^2 - [(x-x_0).(x_1-x_0)] / \|x_1-x_0\|^2
    */
    xx0 =     (x[0] - x0[0]) * (x[0] - x0[0]) +
              (x[1] - x0[1]) * (x[1] - x0[1]) +
              (x[2] - x0[2]) * (x[2] - x0[2]);
    xx0x1x0 = (x[0] - x0[0]) * (x1[0] - x0[0]) +
              (x[1] - x0[1]) * (x1[1] - x0[1]) +
              (x[2] - x0[2]) * (x1[2] - x0[2]);
    *d = sqrt( xx0 - xx0x1x0 * xx0x1x0 / x1x0 );
    ierr = PetscLogFlops(47);CHKERRQ(ierr);
  } else {
    /*
      d = \|x-x_1\|
    */
    *d = sqrt( (x[0] - x1[0]) * (x[0] - x1[0]) +
               (x[1] - x1[1]) * (x[1] - x1[1]) +
               (x[2] - x1[2]) * (x[2] - x1[2]) );
    ierr = PetscLogFlops(24);CHKERRQ(ierr);
  }
  /*
    ierr = PetscPrintf(PETSC_COMM_SELF,"d(%f,%f,%f) = %e\n",x[0],x[1],x[2],*d);
  
  */
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFWellBuildVAT2"
/*
  VFPennyCrackBuildVAT2:  Build the V-field associated with an array of wells 
                          following the construction in Bourdin-Francfort-Marigo '08.
  
  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFWellBuildVAT2(Vec V,VFWell *well,VFCtx *ctx)
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
        ierr = VFDistanceToWell(&dist,x,well);CHKERRQ(ierr);
        v_array[k][j][i] = 1.-exp(-dist/2/ctx->vfprop.epsilon);
      }
    }
  }      
	ierr = DMDAVecRestoreArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}