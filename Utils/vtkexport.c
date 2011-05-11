/*
  vtkexport.c
  (c) 2010 Blaise Bourdin bourdin@lsu.edu
*/

static const char banner[] = "vtkexport:\nconvert petsc binary files into vtk files\n(c) 2010 Blaise Bourdin Louisiana State University bourdin@lsu.edu\n\n";

#include "petsc.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  char                vtkfilename[FILENAME_MAX],petscfilename[FILENAME_MAX];
  char                prefix[FILENAME_MAX];
  PetscViewer         viewer,vtkviewer;
  DA                  daVect,daScal;
  Vec                 coordinates,U,V,temp,pres,pmult;
  PetscInt            nx,ny,nz,step,maxstep;
  PetscErrorCode      ierr;
  FILE                *file;
  
  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"\nh5xfer: ","");CHKERRQ(ierr);
  {
    maxstep = 1000;
    ierr = PetscOptionsInt("-n","maximum number of time steps to convert\t","",maxstep,&maxstep,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscSNPrintf(prefix,FILENAME_MAX,"TEST");CHKERRQ(ierr);
    ierr = PetscOptionsString("-p","file prefix","",prefix,prefix,PETSC_MAX_PATH_LEN-1,PETSC_NULL);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  /*
    Load DA from files
    As of version 3.1, there is a bug in petsc preventing to save 2 DA with different number of degrees of freedoms
    per node in a single file, so we read the Scalar DA, then recreate the Scal DA from the saved informations.
    Also, for some reason coordinates vectors obtained using DAGetCoordinates don't inherit from the proper layout 
    when saved in an hdf5 file. So instead of simply using the coordinate vector obtained from DAScal, we will 
    read it from the file a bit later
  */
  ierr = PetscSNPrintf(petscfilename,FILENAME_MAX,"%s.bin",prefix);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,petscfilename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = DALoad(viewer,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,&daScal);CHKERRQ(ierr);
  
  ierr = DAGetInfo(daScal,0,&nx,&ny,&nz,0,0,0, 0,0,0,0);CHKERRQ(ierr);
  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,nx,ny,nz,
                    PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,3,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                    &(daVect));CHKERRQ(ierr);
  /* 
    Create Vecs from DA
  */
  ierr = DACreateGlobalVector(daScal,&pres);CHKERRQ(ierr);
  ierr = DACreateGlobalVector(daScal,&pmult);CHKERRQ(ierr);
  ierr = DACreateGlobalVector(daScal,&temp);CHKERRQ(ierr);
  ierr = DACreateGlobalVector(daVect,&U);CHKERRQ(ierr);
  ierr = DACreateGlobalVector(daScal,&V);CHKERRQ(ierr);
  ierr = DACreateGlobalVector(daVect,&coordinates);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) pres,"Pressure");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) pmult,"Permeability");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) temp,"Temperature");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) U,"Displacement");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) V,"Fracture");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) coordinates,"Coordinates");CHKERRQ(ierr);
  
  /*
    Get coordinates from petsc binary files, set DA coordinates, and save coordinates in h5 file
  */
  ierr = VecLoadIntoVector(viewer,coordinates);CHKERRQ(ierr);
  ierr = DASetCoordinates(daVect,coordinates);CHKERRQ(ierr);
  ierr = DASetCoordinates(daScal,coordinates);CHKERRQ(ierr);
  
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  /* 
    Open multistep XDMF file
  */  
  for (step = 1; step < maxstep+1; step++) {
    ierr = PetscSNPrintf(petscfilename,FILENAME_MAX,"%s.%.5i.bin",prefix,step);CHKERRQ(ierr);
    if ( (file = fopen(petscfilename,"r")) ) {
      fclose(file);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Processing scalar fields for time step %i\n",step);CHKERRQ(ierr);

      /*
        Open vtk viewer
      */
      ierr = PetscSNPrintf(vtkfilename,FILENAME_MAX,"%s.%.5i.vtk",prefix,step);CHKERRQ(ierr);
      ierr = PetscViewerCreate(PETSC_COMM_WORLD,&vtkviewer);CHKERRQ(ierr);
      ierr = PetscViewerSetType(vtkviewer,PETSC_VIEWER_ASCII);CHKERRQ(ierr);
      ierr = PetscViewerFileSetName(vtkviewer,vtkfilename);CHKERRQ(ierr);
      ierr = PetscViewerSetFormat(vtkviewer,PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
      ierr = DAView(daScal,vtkviewer);CHKERRQ(ierr);
      /* 
        Read and save fields.
        Fields NEED to be read in th eorder they were saved in.
      */
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,petscfilename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  
      ierr = VecLoadIntoVector(viewer,U);CHKERRQ(ierr);
      ierr = VecView(U,vtkviewer);CHKERRQ(ierr);
  
      ierr = VecLoadIntoVector(viewer,V);CHKERRQ(ierr);
      ierr = VecView(V,vtkviewer);CHKERRQ(ierr);
  
      ierr = VecLoadIntoVector(viewer,pmult);CHKERRQ(ierr);
      ierr = VecView(pmult,vtkviewer);CHKERRQ(ierr);
  
      ierr = VecLoadIntoVector(viewer,temp);CHKERRQ(ierr);
      ierr = VecView(temp,vtkviewer);CHKERRQ(ierr);
  
      ierr = VecLoadIntoVector(viewer,pres);CHKERRQ(ierr);
      ierr = VecView(pres,vtkviewer);CHKERRQ(ierr);
  
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(vtkviewer);CHKERRQ(ierr);
    }
  }
  
  ierr = VecDestroy(pres);CHKERRQ(ierr);
  ierr = VecDestroy(pmult);CHKERRQ(ierr);
  ierr = VecDestroy(temp);CHKERRQ(ierr);
  ierr = VecDestroy(U);CHKERRQ(ierr);
  ierr = VecDestroy(V);CHKERRQ(ierr);
  ierr = VecDestroy(coordinates);CHKERRQ(ierr);
  ierr = PetscFinalize();
}
