static char help[] = "Demonstrates using the VFWell Object\n\n";
#include "petsc.h"
#include "../VF/VFWell.h"
#include "../Utils/xdmf.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode      ierr;
  VFWell              *well;
  PetscInt            numwell=0;
  char                prefix[256]="";
  PetscInt            *n,nval=3;
  PetscReal           *l;
  DA                  da,da3;
  Vec                 coord,d,v;
  PetscReal           ****coord_array,***d_array,***v_array;
  PetscInt            xs,xm,ys,ym,zs,zm,nx,ny,nz,k,j,i,w,c;
  PetscReal           d0,d1,*x;
  PetscViewer         h5viewer,XDMFviewer,binviewer;
  PetscReal           epsilon;
  
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);CHKERRQ(ierr);
  
  /* 
    Get geometry informations
  */
  ierr = PetscMalloc(3 * sizeof(PetscInt),&n);CHKERRQ(ierr);
  ierr = PetscMalloc(3 * sizeof(PetscReal),&l);CHKERRQ(ierr);
  for (i = 0; i < 3; i++) n[i] = 10;
  for (i = 0; i < 3; i++) l[i] = 1.;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"\n\ntestwell: test program for well functions","");CHKERRQ(ierr);
  {
    nval = 3;
    ierr = PetscOptionsIntArray("-n","\n\tnumber of cells (default 10), comma separated","",n,&nval,PETSC_NULL);CHKERRQ(ierr);
    nval = 3;
    ierr = PetscOptionsRealArray("-l","\n\tDomain dimensions (default 1.), comma separated","",l,&nval,PETSC_NULL);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  
  /*
    Initializes geometry
  */
  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,n[0],n[1],n[2],
                    PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,
                    PETSC_NULL,PETSC_NULL,PETSC_NULL,&da);CHKERRQ(ierr);
  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,n[0],n[1],n[2],
                    PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,3,1,
                    PETSC_NULL,PETSC_NULL,PETSC_NULL,&da3);CHKERRQ(ierr);
  ierr = DASetUniformCoordinates(da,0.,l[0],0.,l[1],0.,l[2]);CHKERRQ(ierr);
  ierr = DASetUniformCoordinates(da3,0.,l[0],0.,l[1],0.,l[2]);CHKERRQ(ierr);

  ierr = DACreateGlobalVector(da,&d);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) d,"Distance to wells");CHKERRQ(ierr);
  ierr = VecSet(d,0.);CHKERRQ(ierr);
  ierr = DACreateGlobalVector(da,&v);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) v,"V field");CHKERRQ(ierr);
  ierr = VecSet(v,0.);CHKERRQ(ierr);
  ierr = DAGetCoordinates(da3,&coord);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) coord,"Coordinates");CHKERRQ(ierr);

  epsilon = (l[0]/n[0] + l[1]/n[1] + l[2]/n[2])/3.;
  ierr = PetscOptionsGetReal("","-epsilon",&epsilon,PETSC_NULL);CHKERRQ(ierr);
  
  /* 
    Get wells informations
  */
  ierr = PetscOptionsGetInt("","-numwell",&numwell,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscMalloc(numwell*sizeof(VFWell),&well);CHKERRQ(ierr)
  for (i=0; i < numwell; i++) {
    ierr = PetscSNPrintf(prefix,255,"well%i_",i);CHKERRQ(ierr);
    ierr = VFWellCreate(&well[i]);CHKERRQ(ierr)
    ierr = VFWellGet(prefix,&well[i]);CHKERRQ(ierr);
    ierr = VFWellView(&well[i],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  
  /*
    Compute distance to wells and fake v field
  */
  ierr = PetscMalloc(3 * sizeof(PetscReal),&x);CHKERRQ(ierr);
  ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DAVecGetArrayDOF(da3,coord,&coord_array);CHKERRQ(ierr);
  ierr = DAVecGetArray(da,d,&d_array);CHKERRQ(ierr);
  ierr = DAVecGetArray(da,v,&v_array);CHKERRQ(ierr);

  for (k = zs; k < zs+zm; k++) {
    for (j = ys; j < ys+ym; j++) {
      for (i = xs; i < xs+xm; i++) {
        d0 = 1e+30;
        for (c = 0; c < 3; c++) x[c] = coord_array[k][j][i][c];
        d_array[k][j][i] = 1e+30;
        for (w = 0; w < numwell; w++) {
          ierr = VFDistanceToWell(&d1,x,&well[w]);CHKERRQ(ierr);
          if (d1 < d_array[k][j][i]) d_array[k][j][i] = d1;
        }
        v_array[k][j][i] = PetscExpScalar(-d_array[k][j][i]/epsilon);
      }
    }
  }

  ierr = DAVecRestoreArrayDOF(da3,coord,&coord_array);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(da,d,&d_array);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(da,v,&v_array);CHKERRQ(ierr);


  /* 
    save fields in an hdf5 file
  */
  ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,"TEST.h5",FILE_MODE_WRITE,&h5viewer);
  ierr = VecView(coord,h5viewer);CHKERRQ(ierr);
  ierr = VecView(d,h5viewer);CHKERRQ(ierr);
  ierr = VecView(v,h5viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(h5viewer);CHKERRQ(ierr);
  /*
    Write the XDMF Description
  */
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"TEST.xmf",&XDMFviewer);CHKERRQ(ierr);
  ierr = XDMFuniformgridInitialize(XDMFviewer,(PetscReal) 1.0,"TEST.h5");CHKERRQ(ierr); 
  ierr = XDMFtopologyAdd (XDMFviewer,n[0],n[1],n[2],"TEST.h5","Coordinates");CHKERRQ(ierr);
  ierr = XDMFattributeAdd(XDMFviewer,n[0],n[1],n[2],1,"Scalar","Node","TEST.h5","Distance to wells");CHKERRQ(ierr);
  ierr = XDMFattributeAdd(XDMFviewer,n[0],n[1],n[2],1,"Scalar","Node","TEST.h5","V field");CHKERRQ(ierr);
  ierr = XDMFuniformgridFinalize(XDMFviewer);CHKERRQ(ierr); 
  ierr = PetscViewerDestroy(XDMFviewer);CHKERRQ(ierr);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"TEST.bin",FILE_MODE_WRITE,&binviewer);CHKERRQ(ierr);
	ierr = DAView(da,binviewer);CHKERRQ(ierr);
	ierr = VecView(coord,binviewer);CHKERRQ(ierr);
ierr = PetscViewerDestroy(binviewer);CHKERRQ(ierr);

  /*
    Clean up
  */
  ierr = PetscFree(n);CHKERRQ(ierr);
  ierr = PetscFree(l);CHKERRQ(ierr);
  ierr = VecDestroy(d);CHKERRQ(ierr);
  ierr = VecDestroy(v);CHKERRQ(ierr);
  ierr = VecDestroy(coord);CHKERRQ(ierr);
  ierr = DADestroy(da);CHKERRQ(ierr);
  ierr = DADestroy(da3);CHKERRQ(ierr);
  for (i=0; i < numwell; i++) {
    ierr = VFWellDestroy(&well[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(x);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return(0);
}
