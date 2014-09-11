#include "petsc.h"
#include "VFWell.h"

static const char *WellConstraint_Name[] = {
	"PRESSURE",
	"RATE",
	"WellConstraint_Name",
	"",
	0
};

static const char *WellType_Name[] = {
	"INJECTOR",
	"PRODUCER",
	"WellType_Name",
	"",
	0
};

#undef __FUNCT__
#define __FUNCT__ "VFWellGet"
/*
 
 
 VFWellGet (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
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
    nval = 3;
    ierr = PetscOptionsRealArray("-coords","\n\t well x, y & z coordinates  (comma separated).","",well->coords,&nval,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-Qw","\n\t well flow rate","",well->Qw,&well->Qw,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-Pw","\n\t well bottomhole pressure","",well->Pw,&well->Pw,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-rw","\n\t well radius","",well->rw,&well->rw,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsEnum("-constraint","\n\t\n\t well constraint type","",WellConstraint_Name,(PetscEnum)well->condition,(PetscEnum*)&well->condition,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsEnum("-type","\n\t\n\t well type","",WellType_Name,(PetscEnum)well->type,(PetscEnum*)&well->type,PETSC_NULL);CHKERRQ(ierr);
    
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
 
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFWellCreate(VFWell *well)
{
  PetscErrorCode ierr;
  int            i;
  
  PetscFunctionBegin;
  ierr = PetscStrcpy(well->name,"well");CHKERRQ(ierr);
  for (i=0; i<3; i++) well->top[i] = 0.;
  for (i=0; i<3; i++) well->bottom[i] = 0.;
  for (i=0; i<3; i++) well->coords[i] = 0.;
  well->condition = RATE;
  well->type = PRODUCER;
  well->Qw = 0.;
  well->rw = 0.;
  well->Pw = 0.;
  /*
   well->BCV = NONE;
   */
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFWellView"
/*
 VFWellView
 
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFWellView(VFWell *well,PetscViewer viewer)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(viewer,"Well object \"%s\":\n",well->name);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"top:    \t%e \t%e \t%e\n",well->top[0],well->top[1],well->top[2]);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"bottom: \t%e \t%e \t%e\n",well->bottom[0],well->bottom[1],well->bottom[2]);CHKERRQ(ierr);
	
	ierr = PetscViewerASCIIPrintf(viewer,"wellnodes: \t%e \t%e \t%e\n",well->coords[0],well->coords[1],well->coords[2]);CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(viewer,"well rate:     \t%e\n",well->Qw);CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(viewer,"well pressure:     \t%e\n",well->Pw);CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(viewer,"well radius:     \t%e\n",well->rw);CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(viewer,"Well Type:	\"%s\" well under \"%s\" condition\n",WellType_Name[well->type],WellConstraint_Name[well->condition]);CHKERRQ(ierr);
	
	/*
   ierr = PetscViewerASCIIPrintf(viewer,"BCV:    \t%s\n",BCTYPE_NAME[well->BCV]);CHKERRQ(ierr);
   */
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFWellSetName"
/*
 VFWellSetName
 
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
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
 
 
 (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
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
#define __FUNCT__ "VFWellBuildVAT1"
/*
 VFWellBuildVAT1:  Build the V-field associated with an array of wells
 
 (c) 2010-2014 Blaise Bourdin bourdin@lsu.edu
 */
extern PetscErrorCode VFWellBuildVAT1(Vec V,VFWell *well,VFCtx *ctx)
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
        if (well->rw <= 0.) {
          v_array[k][j][i] = 1.;
        } else {
          ierr = VFDistanceToWell(&dist,x,well);CHKERRQ(ierr);
          dist -= well->rw;
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
  }
	ierr = DMDAVecRestoreArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFWellBuildVAT2"
/*
 VFWellBuildVAT2:  Build the V-field associated with an array of wells
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
        if(dist <= ctx->well->rw*2.0){
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
#define __FUNCT__ "VFFakeWellBuildVAT2"
/*
 (no fracture toughness zone)
 */
extern PetscErrorCode VFFakeWellBuildVAT2(VFWell *well,VFCtx *ctx)
{
  PetscErrorCode      ierr;
  PetscInt            ei,ej,ek,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal       ****coords_array;
  PetscReal        ***Gc_array;
  PetscReal           x[3];
  PetscReal           dist;
  
  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  
  ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  
  ierr = DMDAVecGetArray(ctx->daScalCell,ctx->matprop->VecGc,&Gc_array);CHKERRQ(ierr);
  for (ek = zs; ek < zs+zm; ek++) {
    for (ej = ys; ej < ys+ym; ej++) {
      for (ei = xs; ei < xs+xm; ei++) {
        x[2] = coords_array[ek][ej][ei][2];
        x[1] = coords_array[ek][ej][ei][1];
        x[0] = coords_array[ek][ej][ei][0];
        ierr = VFDistanceToWell(&dist,x,well);CHKERRQ(ierr);
        if(dist <= ctx->well->rw*2.0){
          Gc_array[ek][ej][ei] = 1.e-9;
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(ctx->daScalCell,ctx->matprop->VecGc,&Gc_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFRegDiracDeltaFunction2"
extern PetscErrorCode VFRegDiracDeltaFunction2(Vec RegV,VFWell *well,VFPennyCrack *crack,VFCtx *ctx, Vec V)
{
  PetscErrorCode      ierr;
  PetscInt            i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal           ****coords_array;
  PetscReal           ***v_array;
  PetscReal           ***regv_array;
  PetscReal           x[3],x0[3];
  PetscReal           dist;
  PetscReal           epsilon;
  
  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  epsilon = ctx->vfprop.epsilon;
  for (i = 0; i < 3; i++) x0[i] = well->coords[i];
  
	ierr = DMDAVecGetArray(ctx->daScal,RegV,&regv_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  for (k = zs; k < zs+zm; k++) {
    for (j = ys; j < ys+ym; j++) {
      for (i = xs; i < xs+xm; i++) {
        x[2] = coords_array[k][j][i][2];
        x[1] = coords_array[k][j][i][1];
        x[0] = coords_array[k][j][i][0];
        dist =    sqrt((x[0] - x0[0]) * (x[0] - x0[0]) +
                       (x[1] - x0[1]) * (x[1] - x0[1]) +
                       (x[2] - x0[2]) * (x[2] - x0[2]));
        if(dist <= crack->thickness/2+epsilon/3){
          regv_array[k][j][i] += well->Qw;
        }
        else {
          regv_array[k][j][i] += 0;
        }
      }
    }
  }
	ierr = DMDAVecRestoreArray(ctx->daScal,RegV,&regv_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VFRegDiracDeltaFunction1"
extern PetscErrorCode VFRegDiracDeltaFunction1(Vec V,VFWell *well,VFPennyCrack *crack,VFCtx *ctx)
{
  PetscErrorCode      ierr;
  PetscInt            i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal           ****coords_array;
  PetscReal           ***v_array;
  PetscReal           x[3],x0[3];
  PetscReal           dist;
  
  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  
  for (i = 0; i < 3; i++) x0[i] = well->coords[i];
  
	ierr = DMDAVecGetArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  for (k = zs; k < zs+zm; k++) {
    for (j = ys; j < ys+ym; j++) {
      for (i = xs; i < xs+xm; i++) {
        x[2] = coords_array[k][j][i][2];
        x[1] = coords_array[k][j][i][1];
        x[0] = coords_array[k][j][i][0];
        dist =    sqrt((x[0] - x0[0]) * (x[0] - x0[0]) +
                       (x[1] - x0[1]) * (x[1] - x0[1]) +
                       (x[2] - x0[2]) * (x[2] - x0[2]));
        if(dist <= crack->thickness/2.){
          v_array[k][j][i] = well->Qw;
        }
        else if (dist < 2. * ctx->vfprop.epsilon+crack->thickness/2.) {
          v_array[k][j][i] = well->Qw*(1.-dist/ctx->vfprop.epsilon * (1.- .25*dist/ctx->vfprop.epsilon));
        }
        else {
          v_array[k][j][i] = 0.;
        }
      }
    }
  }
	ierr = DMDAVecRestoreArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "VFRegDiracDeltaFunction"
extern PetscErrorCode VFRegDiracDeltaFunction(Vec RegV,VFWell *well,VFPennyCrack *crack,VFCtx *ctx, Vec V)
{
  PetscErrorCode      ierr;
  PetscInt            i,j,k,nx,ny,nz,xs,xm,ys,ym,zs,zm;
  PetscReal           ****coords_array;
  PetscReal           ***v_array,***regv_array;
  PetscReal           x[3],x0[3];
  PetscReal           dist;
  PetscFunctionBegin;
  ierr = DMDAGetInfo(ctx->daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  
  ierr = DMDAVecGetArray(ctx->daScal,RegV,&regv_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  for (i = 0; i < 3; i++) x0[i] = well->coords[i];
	ierr = DMDAVecGetArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  for (k = zs; k < zs+zm; k++) {
    for (j = ys; j < ys+ym; j++) {
      for (i = xs; i < xs+xm; i++) {
        x[2] = coords_array[k][j][i][2];
        x[1] = coords_array[k][j][i][1];
        x[0] = coords_array[k][j][i][0];
        dist =    sqrt((x[0] - x0[0]) * (x[0] - x0[0]) +
                       (x[1] - x0[1]) * (x[1] - x0[1]) +
                       (x[2] - x0[2]) * (x[2] - x0[2]));
        if(dist <= crack->thickness/2.){
          regv_array[k][j][i] += well->Qw*PETSC_PI/(4.*pow(ctx->vfprop.epsilon,1));
        }
        else {
          regv_array[k][j][i] += well->Qw*exp(-dist/ctx->vfprop.epsilon)*PETSC_PI/(4.*pow(ctx->vfprop.epsilon,1));
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(ctx->daScal,RegV,&regv_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,V,&v_array);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
