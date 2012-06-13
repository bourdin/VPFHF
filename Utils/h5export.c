/*
  h5export.c
  (c) 2010-2012 Blaise Bourdin bourdin@lsu.edu
*/

static const char banner[] = "h5export:\nconvert petsc binary files into hdf5 / xmf files\n(c) 2010-2012 Blaise Bourdin Louisiana State University bourdin@lsu.edu\n\n";

#include "petsc.h"
#include "../xdmf.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  char                h5filename[FILENAME_MAX],petscfilename[FILENAME_MAX],XDMFfilename[FILENAME_MAX];
  char                h5coordfilename[FILENAME_MAX],prefix[FILENAME_MAX];
  PetscViewer         viewer,h5viewer,XDMFviewer,XDMFviewer2;
  DM                  daVect,daScal;
  Vec                 coordinates,U,V,temp,pres,pmult,FVel,FVCellndof;
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

  ierr = PetscSNPrintf(petscfilename,FILENAME_MAX,"%s.bin",prefix);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,petscfilename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = DMCreate(PETSC_COMM_WORLD,&daScal);CHKERRQ(ierr);
  ierr = DMLoad(daScal,viewer);CHKERRQ(ierr);
  ierr = DMView(daScal,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  
  //ierr = DMDAGetReducedDA(daScal,1,&daVect);CHKERRQ(ierr);
  ierr = DMDAGetInfo(daScal,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,
						DMDA_STENCIL_BOX,nx,ny,nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,3,1,
						PETSC_NULL,PETSC_NULL,PETSC_NULL,&daVect);CHKERRQ(ierr);
  /* 
    Create Vecs from DA
  */
  ierr = VecCreate(PETSC_COMM_WORLD,&pres);CHKERRQ(ierr);
  ierr = VecSetBlockSize(pres,1);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&pmult);CHKERRQ(ierr);
  ierr = VecSetBlockSize(pmult,1);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&temp);CHKERRQ(ierr);
  ierr = VecSetBlockSize(temp,1);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&U);CHKERRQ(ierr);
  ierr = VecSetBlockSize(U,3);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&V);CHKERRQ(ierr);
  ierr = VecSetBlockSize(V,1);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&FVel);CHKERRQ(ierr);
  ierr = VecSetBlockSize(FVel,3);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&FVCellndof);CHKERRQ(ierr);
  ierr = VecSetBlockSize(FVCellndof,1);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&coordinates);CHKERRQ(ierr);
  ierr = VecSetBlockSize(coordinates,3);CHKERRQ(ierr);

  ierr = PetscObjectSetName((PetscObject) pres,"Pressure");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) pmult,"Permeability");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) temp,"Temperature");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) U,"Displacement");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) V,"Fracture");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) FVel,"Fluid_Velocity");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) FVCellndof,"FVCellndof");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) coordinates,"Coordinates");CHKERRQ(ierr);
  
  /*
    Get coordinates from petsc binary files, set DA coordinates, and save coordinates in h5 file
  */
  //ierr = PetscSNPrintf(h5coordfilename,FILENAME_MAX,"%s.bin",prefix);CHKERRQ(ierr);
  //ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,h5coordfilename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(coordinates,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = DMDASetCoordinates(daVect,coordinates);CHKERRQ(ierr);
  ierr = DMDASetCoordinates(daScal,coordinates);CHKERRQ(ierr);
  
  ierr = PetscSNPrintf(h5filename,FILENAME_MAX,"%s.h5",prefix);CHKERRQ(ierr);
  ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,h5filename,FILE_MODE_WRITE,&h5viewer);
  ierr = VecView(coordinates,h5viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&h5viewer);CHKERRQ(ierr);
  
  /* 
    Open multistep XDMF file
  */
  ierr = PetscSNPrintf(XDMFfilename,FILENAME_MAX,"%s.xmf",prefix);CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,XDMFfilename,&XDMFviewer2);CHKERRQ(ierr);
  ierr = XDMFmultistepInitialize(XDMFviewer2);CHKERRQ(ierr);
  
  for (step = 1; step < maxstep+1; step++) {
    /*
      Tests if file exists
    */
    ierr = PetscSNPrintf(petscfilename,FILENAME_MAX,"%s.%.5i.bin",prefix,step);CHKERRQ(ierr);
    if ( (file = fopen(petscfilename,"r")) ) {
      fclose(file);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Processing fields for time step %i\n",step);CHKERRQ(ierr);
      /* 
        Read and save fields.
        Fields NEED to be read in the order they were saved in.
      */
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,petscfilename,FILE_MODE_READ,&viewer);
  
      /*
        Open h5 viewer
      */
      ierr = PetscSNPrintf(h5filename,FILENAME_MAX,"%s.%.5i.h5",prefix,step);CHKERRQ(ierr);
      ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,h5filename,FILE_MODE_WRITE,&h5viewer);
  
      ierr = VecLoad(U,viewer);CHKERRQ(ierr);
      ierr = VecView(U,h5viewer);CHKERRQ(ierr);
  
      ierr = VecLoad(FVel,viewer);CHKERRQ(ierr);
      ierr = VecView(FVel,h5viewer);CHKERRQ(ierr);
  
      ierr = VecLoad(V,viewer);CHKERRQ(ierr);
      ierr = VecView(V,h5viewer);CHKERRQ(ierr);
  
      ierr = VecLoad(pmult,viewer);CHKERRQ(ierr);
      ierr = VecView(pmult,h5viewer);CHKERRQ(ierr);
  
      ierr = VecLoad(temp,viewer);CHKERRQ(ierr);
      ierr = VecView(temp,h5viewer);CHKERRQ(ierr);
  
      ierr = VecLoad(pres,viewer);CHKERRQ(ierr);
      ierr = VecView(pres,h5viewer);CHKERRQ(ierr);

      ierr = VecLoad(FVCellndof,viewer);CHKERRQ(ierr);
      ierr = VecView(FVCellndof,h5viewer);CHKERRQ(ierr);
	  
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&h5viewer);CHKERRQ(ierr);
      /*
        Write the XDMF Description
      */
      ierr = PetscSNPrintf(XDMFfilename,FILENAME_MAX,"%s.%.5i.xmf",prefix,step);CHKERRQ(ierr);
      ierr = PetscSNPrintf(h5coordfilename,FILENAME_MAX,"%s.h5",prefix);CHKERRQ(ierr);
  
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,XDMFfilename,&XDMFviewer);CHKERRQ(ierr);
      ierr = XDMFuniformgridInitialize(XDMFviewer,(PetscReal) step,h5filename);CHKERRQ(ierr); 
      ierr = XDMFtopologyAdd (XDMFviewer,nx,ny,nz,h5coordfilename,"Coordinates");CHKERRQ(ierr);
      ierr = XDMFattributeAdd(XDMFviewer,nx,ny,nz,3,"Vector","Node",h5filename,"Displacement");CHKERRQ(ierr);
      ierr = XDMFattributeAdd(XDMFviewer,nx,ny,nz,3,"Vector","Node",h5filename,"Fluid_Velocity");CHKERRQ(ierr);
      ierr = XDMFattributeAdd(XDMFviewer,nx,ny,nz,1,"Scalar","Node",h5filename,"Fracture");CHKERRQ(ierr);
      ierr = XDMFattributeAdd(XDMFviewer,nx,ny,nz,1,"Scalar","Node",h5filename,"Permeability");CHKERRQ(ierr);
      ierr = XDMFattributeAdd(XDMFviewer,nx,ny,nz,1,"Scalar","Node",h5filename,"Temperature");CHKERRQ(ierr);
      ierr = XDMFattributeAdd(XDMFviewer,nx,ny,nz,1,"Scalar","Node",h5filename,"Pressure");CHKERRQ(ierr);
      ierr = XDMFattributeAdd(XDMFviewer,nx,ny,nz,1,"Scalar","Node",h5filename,"FVCellndof");CHKERRQ(ierr);	  
      ierr = XDMFuniformgridFinalize(XDMFviewer);CHKERRQ(ierr); 
      ierr = PetscViewerDestroy(&XDMFviewer);CHKERRQ(ierr);
      
      ierr = XDMFmultistepAddstep(XDMFviewer2,XDMFfilename);CHKERRQ(ierr);
    }
  }  

  ierr = XDMFmultistepFinalize(XDMFviewer2);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&XDMFviewer2);CHKERRQ(ierr);

  ierr = VecDestroy(&pres);CHKERRQ(ierr);
  ierr = VecDestroy(&pmult);CHKERRQ(ierr);
  ierr = VecDestroy(&temp);CHKERRQ(ierr);
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&V);CHKERRQ(ierr);
  ierr = VecDestroy(&FVel);CHKERRQ(ierr);
  ierr = VecDestroy(&FVCellndof);CHKERRQ(ierr);  
  ierr = VecDestroy(&coordinates);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}
