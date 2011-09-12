/*
  VFracture.c
  (c) 2010-2011 Blaise Bourdin bourdin@lsu.edu
*/

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFV.h"
#include "VFU.h"
#include "VFFlow.h"

VFCtx               ctx;
VFFields            fields;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode      ierr;
  
  /*
    Iteration
  */
  PetscReal           ElasticEnergy = 0;
  PetscReal           InsituWork = 0;
  PetscReal           SurfaceEnergy = 0;
  PetscInt            DebugMode = 0;
  FILE                *file;
  char                H5filename[FILENAME_MAX],filename[FILENAME_MAX];
  PetscInt            *n,i,nval=3;
  PetscReal           *l,*dx,*dy,*dz;
  PetscTruth          flg;
  /*
    This is essentially VIADAT
  */
  ierr = PetscInitialize(&argc,&argv,(char*)0,banner);CHKERRQ(ierr);
  
  ierr = PetscOptionsGetInt("","-debug",&DebugMode,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"\n\nVF-Chevron: uncoupled variational fracture model specific options:","");CHKERRQ(ierr);
  {
    ierr = PetscMalloc(3 * sizeof(PetscInt),&n);CHKERRQ(ierr);
    for (i = 0; i < 3; i++) n[i] = 10;
    ierr = PetscOptionsIntArray("-n","\n\tnumber of cells (default 10), comma separated","",n,&nval,PETSC_NULL);CHKERRQ(ierr);
    if (nval != 3 && nval != 0) {
      SETERRQ4(PETSC_ERR_USER,"ERROR: Expecting 3 values for option %s, got only %i in %s\n",n,"-n",nval,__FUNCT__);
    }

    ierr = PetscMalloc(3 * sizeof(PetscReal),&l);CHKERRQ(ierr);
    for (i = 0; i < 3; i++) l[i] = 1.;
    ierr = PetscOptionsRealArray("-l","\n\tDomain dimensions (default 1.), comma separated","",l,&nval,PETSC_NULL);CHKERRQ(ierr);
    if (nval != 3 && nval != 0) {
      SETERRQ4(PETSC_ERR_USER,"ERROR: Expecting 3 values for option %s, got only %i in %s\n",n,"-l",nval,__FUNCT__);
    }
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  ierr = PetscMalloc3(n[0],PetscReal,&dx,n[1],PetscReal,&dy,n[2],PetscReal,&dz);CHKERRQ(ierr);


  for (i = 0; i < n[0]; i++) dx[i] = l[0]/n[0];
  for (i = 0; i < n[1]; i++) dy[i] = l[1]/n[1];
  for (i = 0; i < n[2]; i++) dz[i] = l[2]/n[2];

  ierr = VFInitialize(n[0],n[1],n[2],dx,dy,dz);CHKERRQ(ierr);

 /* end VIADAT */

  /* start of time step */
  ctx.timestep = 1 ;

  switch (ctx.fileformat) {
   case FILEFORMAT_HDF5:       
     ierr = FieldsH5Write(&ctx,&fields);
   break;
   case FILEFORMAT_BIN:
     ierr = FieldsBinaryWrite(&ctx,&fields);
   break; 
  } 
  ctx.timestep++; 
  ierr = PetscSNPrintf(H5filename,FILENAME_MAX,"%s.%.5i.h5",ctx.prefix,ctx.timestep);CHKERRQ(ierr);
  while ( file = fopen(H5filename,"w") ) {
    fclose(file);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nProcessing step %i (%s).\n",ctx.timestep,H5filename);CHKERRQ(ierr);
    ctx.timevalue = (PetscReal) ctx.timestep;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\ntime value %f \n",ctx.timevalue);CHKERRQ(ierr);
    /*
      This is called in VTDATA
    */
    ierr = VFTimeStepPrepare(&ctx,&fields);CHKERRQ(ierr);
    /*
      End VTDATA
    */
   
    /*
      Update pressure and theta here
    */

    ierr = VFFlowTimeStep(&ctx,&fields);CHKERRQ(ierr);
    

    /*
      This is essentially VPERM (minus pmult update, and transfer from GMRS)
    */
    switch (ctx.mode) {
      case ELASTICITY:
        ierr = VFElasticityTimeStep(&ctx,&fields);CHKERRQ(ierr);
      break;
      case FRACTURE:
        ierr = VFFractureTimeStep(&ctx,&fields);CHKERRQ(ierr);
      break;
    }

    /*
      This is essentially VSTOUT
    */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic Energy:            %e\n",ctx.ElasticEnergy);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Work of surface forces:    %e\n",ctx.InsituWork);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Surface energy:            %e\n",ctx.SurfaceEnergy);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Total energy:              %e\n",ctx.ElasticEnergy+SurfaceEnergy-InsituWork);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(ctx.energyviewer,"%i   \t%e   \t%e   \t%e   \t%e\n",ctx.timestep,ElasticEnergy,InsituWork,SurfaceEnergy,
                                  ElasticEnergy - InsituWork + SurfaceEnergy);CHKERRQ(ierr); 

     switch (ctx.fileformat) {
      case FILEFORMAT_HDF5:       
        ierr = FieldsH5Write(&ctx,&fields);
      break;
      case FILEFORMAT_BIN:
        ierr = FieldsBinaryWrite(&ctx,&fields);
      break; 
    } 
    ierr = PetscSNPrintf(filename,FILENAME_MAX,"%s.log",ctx.prefix);CHKERRQ(ierr);
    ierr = PetscLogPrintSummary(PETSC_COMM_WORLD,filename);CHKERRQ(ierr);

    /* 
      flushes out statistics about the current run
    */
    
    /*
      build filename of next step and loops
    */
    ctx.timestep++;
    if((ctx.timestep > ctx.maxtimestep) || (ctx.timevalue > ctx.maxtimevalue)) break;
    ierr = PetscSNPrintf(H5filename,FILENAME_MAX,"%s.%.5i.h5",ctx.prefix,ctx.timestep);CHKERRQ(ierr);
  }
  /* end of time step */

  ierr = VFFinalize(&ctx,&fields);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}
  
