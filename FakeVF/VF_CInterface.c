#include "petsc.h"
#include "VF_CInterface.h"
#include "../Utils/xdmf.h"
#include "../Utils/VFGMRS.h"

VFCtx user;

#undef __FUNCT__
#define __FUNCT__ "VFInitialize"
/*
  VFInitialize: Initialize the VF code. Called by the fortran implementation of VIADAT

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VFInitialize(PetscInt nx,PetscInt ny,PetscInt nz,PetscReal *dx,PetscReal *dy,PetscReal *dz)
{
  PetscErrorCode ierr;
  PetscReal      ****coords_array;
  PetscInt       xs,xm,ys,ym,zs,zm;
  PetscInt       i,j,k;
  PetscReal      *X,*Y,*Z;
  char           filename[FILENAME_MAX];
  PetscViewer    viewer,h5viewer;  
  PetscTruth     flg;

  PetscFunctionBegin;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"\nFakeVF: ","");CHKERRQ(ierr);
  {
    ierr = PetscSNPrintf(user.prefix,FILENAME_MAX,"TEST");CHKERRQ(ierr);
    user.fileformat = FILEFORMAT_BIN;
    ierr = PetscOptionsEnum("-format","FileFormat","",VFFileFormatName,(PetscEnum)user.fileformat,(PetscEnum*)&user.fileformat,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsString("-p","file prefix","",user.prefix,user.prefix,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  
  ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.log",user.prefix);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n log file is %s\n\n",filename);CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&user.logviewer);CHKERRQ(ierr);  
  ierr = PetscViewerASCIIPrintf(user.logviewer,"This is %s\n",__FUNCT__);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"This is the PETSc / C implementation of %s\n",__FUNCT__);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"nx      = %i\n",nx);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"ny      = %i\n",ny);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"nz      = %i\n",nz);CHKERRQ(ierr);

  user.nx = nx-1;
  user.ny = ny-1;
  user.nz = nz-1;

  /*
    Create a DA
  */
  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,user.nx+1,user.ny+1,user.nz+1,
                    PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,3,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                    &(user.daVect));CHKERRQ(ierr);
  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,user.nx+1,user.ny+1,user.nz+1,
                    PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                    &(user.daScal));CHKERRQ(ierr);

  /*
    Set its coordinates from dx, dy, dz
  */
  ierr = DACreateGlobalVector(user.daVect,&(user.coordinates));CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) user.coordinates,"Coordinates");CHKERRQ(ierr);
  ierr = DAVecGetArrayDOF(user.daVect,user.coordinates,&coords_array);CHKERRQ(ierr);

  ierr = DAGetCorners(user.daVect,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  ierr = PetscMalloc3(nx+1,PetscReal,&X,ny+1,PetscReal,&Y,nz+1,PetscReal,&Z);CHKERRQ(ierr);
  Z[0] = 0.;
  Y[0] = 0.;
  X[0] = 0.;
  for (k = 1; k < nz+1; k++) Z[k] = Z[k-1] + dz[k-1];
  for (j = 1; j < ny+1; j++) Y[j] = Y[j-1] + dy[j-1];
  for (i = 1; i < nx+1; i++) X[i] = X[i-1] + dx[i-1];
  
  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        coords_array[k][j][i][2] = Z[k];
        coords_array[k][j][i][1] = Y[j];
        coords_array[k][j][i][0] = X[i]; 
      }
    }
  }
  ierr = PetscFree3(X,Y,Z);CHKERRQ(ierr);

  ierr = DAVecRestoreArrayDOF(user.daVect,user.coordinates,&coords_array);CHKERRQ(ierr);

  ierr = DASetCoordinates(user.daVect,user.coordinates);CHKERRQ(ierr);
  ierr = DASetCoordinates(user.daScal,user.coordinates);CHKERRQ(ierr);

  ierr = DACreateGlobalVector(user.daScal,&(user.pmult));CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) user.pmult,"Permeability");CHKERRQ(ierr);
  ierr = VecSet(user.pmult,0.);CHKERRQ(ierr);
  ierr = DACreateGlobalVector(user.daScal,&(user.tempr));CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) user.tempr,"Temperature");CHKERRQ(ierr);
  ierr = DACreateGlobalVector(user.daScal,&(user.pres));CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) user.pres,"Pressure");CHKERRQ(ierr);
  ierr = DACreateGlobalVector(user.daVect,&(user.U));CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) user.U,"Displacement");CHKERRQ(ierr);
  ierr = VecSet(user.U,0.);CHKERRQ(ierr);
  ierr = DACreateGlobalVector(user.daScal,&(user.V));CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) user.V,"Fracture");CHKERRQ(ierr);
  ierr = VecSet(user.V,1.);CHKERRQ(ierr);
  
  /* 
    Saves DA informations
  */
	ierr = DAView(user.daScal,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  switch (user.fileformat) {
    case FILEFORMAT_BIN:
      ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.bin",user.prefix);CHKERRQ(ierr);
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = DAView(user.daScal,viewer);CHKERRQ(ierr);
      ierr = VecView(user.coordinates,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
      break;
    case FILEFORMAT_HDF5:
#ifdef PETSC_HAVE_HDF5
      /*
        Write headers in multistep XDMF file
      */
      ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.xmf",user.prefix);CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&user.XDMFviewer);
      ierr = XDMFmultistepInitialize(user.XDMFviewer);CHKERRQ(ierr);

      ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.h5",user.prefix);CHKERRQ(ierr);
      ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&h5viewer);CHKERRQ(ierr);
      ierr = VecView(user.coordinates,h5viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(h5viewer);CHKERRQ(ierr);
#endif
    break;
  }
      
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "vendit"
/*
  vendit:

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode vendit(PetscInt rank,PetscReal tim,PetscInt nstep,PetscInt nfout,PetscInt nfbug)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(user.logviewer,"This is %s\n",__FUNCT__);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(user.logviewer,"\tTime=%g\tnstep=%i\n",tim,nstep);CHKERRQ(ierr);

  /*
    Closes XDMF files
  */
#ifdef PETSC_HAVE_HDF5
  if (user.fileformat == FILEFORMAT_HDF5) {
    ierr = XDMFuniformgridFinalize(user.XDMFviewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(user.XDMFviewer);CHKERRQ(ierr);
  }
#endif
  ierr = PetscViewerDestroy(user.logviewer);CHKERRQ(ierr); 
  ierr = PetscFinalize();
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "vtdata"
/*
  vtdata: 

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode vtdata(PetscInt rank,PetscReal tim,PetscInt nstep,PetscInt nfout,PetscInt nfbug)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(user.logviewer,"\nThis is %s\n",__FUNCT__);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(user.logviewer,"\tTime=\t%g\n\tnstep=\t%i\n",tim,nstep);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nThis is the PETSc / C implementation of %s\n",__FUNCT__);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"   TIME = %15.5E\n",tim);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"   NSTEP = %10i\n",nstep);CHKERRQ(ierr);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"   user.nx = %i\n",user.nx);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"   user.ny = %i\n",user.ny);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"   user.nz = %i\n",user.nz);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "vperm"
/*
  vperm: This is where the permeability is updated from the pressure + temperature, i.e. the main hook for the fracture computation

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode vperm(PetscInt rank,PetscReal *pres,PetscReal *tempr,PetscReal *pmult,PetscReal tim,PetscInt nstep,PetscInt newt,PetscInt nfout,PetscInt nfbug,PetscInt n)
{
  PetscReal      fieldmin,fieldmax;   
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(user.logviewer,"\nThis is %s\n",__FUNCT__);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(user.logviewer,"\tTime=%g\tnstep=%i\tnewt=%i\n",tim,nstep,newt);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nThis is the PETSc / C implementation of %s\n",__FUNCT__);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"   TIME  = %15.5E\n",tim);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"   NSTEP = %10i\n",nstep);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"   NEWT  = %10i\n",newt);CHKERRQ(ierr);

  /*
    Get Pressure and temperature from GMRS
  */
  ierr = VecGMRSToPetsc(pres,user.pres);CHKERRQ(ierr);  
  ierr = VecGMRSToPetsc(tempr,user.tempr);CHKERRQ(ierr);  

  ierr = VecMin(user.tempr,PETSC_NULL,&fieldmin);CHKERRQ(ierr);
  ierr = VecMax(user.tempr,PETSC_NULL,&fieldmax);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Temperature range: %g - %g\n",fieldmin,fieldmax);CHKERRQ(ierr);
  ierr = VecMin(user.pres,PETSC_NULL,&fieldmin);CHKERRQ(ierr);
  ierr = VecMax(user.pres,PETSC_NULL,&fieldmax);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Pressure range:    %g - %g\n",fieldmin,fieldmax);CHKERRQ(ierr);
  /* 
    Compute Displacement and Fracture, then derive permeability multipliers
  */
  /*
    Send permeability multiplier back to GMRS
  */
  /*
  p1 = .5;
  p2 = .75;
  ierr = VecCopy(user.tempr,user.pmult);CHKERRQ(ierr);
  ierr = VecScale(user.pmult,-1.0*(p2-p1)/(tempmin-tempmax));CHKERRQ(ierr);
  ierr = VecShift(user.pmult,(p2*tempmin-p1*tempmax)/(tempmin-tempmax));CHKERRQ(ierr);
  ierr = VecPetscToGMRS(user.pmult,pmult);CHKERRQ(ierr); 
  */
  /*
    Calling my vstdout since this is where I write fields
  */
  /*
  ierr = vstdout(rank,tim,nstep,nfout,nfbug);CHKERRQ(ierr);
  */
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "vstdout"
/*
  vstdout: 

  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode vstdout(PetscInt rank,PetscReal tim,PetscInt nstep,PetscInt nfout,PetscInt nfbug)
{
  char           filename[FILENAME_MAX],h5filename[FILENAME_MAX],h5coordfilename[FILENAME_MAX],XDMFfilename[FILENAME_MAX];  
  PetscViewer    viewer,XDMFviewer,h5viewer;
  PetscErrorCode ierr;
   
  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(user.logviewer,"\nThis is %s\n",__FUNCT__);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(user.logviewer,"\tTime=%g\tnstep=%i\n",tim,nstep);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nThis is the PETSc / C implementation of %s\n",__FUNCT__);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"   TIME  = %15.5g\n",tim);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"   NSTEP = %10i\n",nstep);CHKERRQ(ierr);

  switch (user.fileformat) {
    case FILEFORMAT_BIN:
      /*
        Save all fields in the same petsc binary file
      */
      ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.%.5i.bin",user.prefix,nstep);CHKERRQ(ierr);
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = VecView(user.U,viewer);CHKERRQ(ierr);
      ierr = VecView(user.V,viewer);CHKERRQ(ierr);
      ierr = VecView(user.pmult,viewer);CHKERRQ(ierr);
      ierr = VecView(user.tempr,viewer);CHKERRQ(ierr);
      ierr = VecView(user.pres,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
      break;
    case FILEFORMAT_HDF5:
      /*
        Save the actual data
      */
#ifdef PETSC_HAVE_HDF5
      ierr = PetscSNPrintf(h5filename,FILENAME_MAX,"%s.%.5i.h5",user.prefix,nstep);CHKERRQ(ierr);
      ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,h5filename,FILE_MODE_WRITE,&h5viewer);
      ierr = VecView(user.U,h5viewer);CHKERRQ(ierr);
      ierr = VecView(user.V,h5viewer);CHKERRQ(ierr);
      ierr = VecView(user.pmult,h5viewer);CHKERRQ(ierr);
      ierr = VecView(user.tempr,h5viewer);CHKERRQ(ierr);
      ierr = VecView(user.pres,h5viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(h5viewer);CHKERRQ(ierr);
      /*
        Write the XDMF Description
      */
      ierr = PetscSNPrintf(XDMFfilename,FILENAME_MAX,"%s.%.5i.xmf",user.prefix,nstep);CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,XDMFfilename,&XDMFviewer);CHKERRQ(ierr);
      ierr = XDMFuniformgridInitialize(XDMFviewer,tim,h5filename);CHKERRQ(ierr); 
      ierr = PetscSNPrintf(h5coordfilename,FILENAME_MAX,"%s.h5",user.prefix);CHKERRQ(ierr);
      ierr = XDMFtopologyAdd(XDMFviewer,user.nx+1,user.ny+1,user.nz+1,h5coordfilename,"Coordinates");CHKERRQ(ierr);
      ierr = XDMFattributeAdd(XDMFviewer,user.nx+1,user.ny+1,user.nz+1,3,"Vector","Node",h5filename,"Displacement");CHKERRQ(ierr);
      ierr = XDMFattributeAdd(XDMFviewer,user.nx+1,user.ny+1,user.nz+1,1,"Scalar","Node",h5filename,"Fracture");CHKERRQ(ierr);
      ierr = XDMFattributeAdd(XDMFviewer,user.nx+1,user.ny+1,user.nz+1,1,"Scalar","Node",h5filename,"Permeability");CHKERRQ(ierr);
      ierr = XDMFattributeAdd(XDMFviewer,user.nx+1,user.ny+1,user.nz+1,1,"Scalar","Node",h5filename,"Temperature");CHKERRQ(ierr);
      ierr = XDMFattributeAdd(XDMFviewer,user.nx+1,user.ny+1,user.nz+1,1,"Scalar","Node",h5filename,"Pressure");CHKERRQ(ierr);
      ierr = XDMFuniformgridFinalize(XDMFviewer);CHKERRQ(ierr); 
      ierr = PetscViewerDestroy(XDMFviewer);CHKERRQ(ierr);

      /*
        Add an entry in the multistep XDMF file
      */
      ierr = XDMFmultistepAddstep(user.XDMFviewer,XDMFfilename);CHKERRQ(ierr);
#endif
    break;
  }
 PetscFunctionReturn(0);
}
