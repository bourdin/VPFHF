#include "petsc.h"

#undef __FUNCT__
#define __FUNCT__ "VecGMRSToPetsc"
/*
  VecGMRSToPetsc: convert a GMRS 1d array sent from CPU 0 into a distributed Petsc Vec

  (c) 2010 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VecGMRSToPetsc(PetscReal *GMRS_array,Vec x)
{
  Vec				     natural,zero_Vec;
  PetscReal      *zero_array;
  MPI_Comm       comm;
  VecScatter		 tozero;
  PetscMPIInt	   rank;
  DA				     da;
  PetscInt	     i,j,k,l,mx,my,mz,n,N;
  PetscErrorCode ierr;

  
  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject) x,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  
  ierr = PetscObjectQuery((PetscObject) x,"DA",(PetscObject *) &da);CHKERRQ(ierr);
  if (!da) SETERRQ(PETSC_ERR_ARG_WRONG,"Vector not generated from a DA");
  ierr = DAGetInfo(da,0,&mx,&my,&mz,0,0,0, 0,0,0,0);CHKERRQ(ierr);
  N = mx * my * mz;
  if (!rank) {
    n = N;
  } else {
    n = 0;
  }
  
  ierr = DACreateNaturalVector(da,&natural);CHKERRQ(ierr);
  ierr = VecScatterCreateToZero(natural,&tozero,&zero_Vec);CHKERRQ(ierr);
  
  /*
    This would work if the dof ordering scheme for petsc and GMRs would coincide. 
    ierr = VecCreateMPIWithArray(comm,n,N,zero_array,&GMRS_Vec);CHKERRQ(ierr);
  */
  ierr = VecGetArray(zero_Vec,&zero_array);CHKERRQ(ierr);
  if (!rank) {
    for (l=0,i = 0; i < mx; i++) {
      for (j = 0; j < my; j++) {
        for (k = 0; k < mz; k++,l++) {
          zero_array[i + j*mx + k*mx*my] = GMRS_array[l];
        }
      }
    }
  }
  ierr = VecRestoreArray(zero_Vec,&zero_array);CHKERRQ(ierr);
  
  ierr = VecScatterBegin(tozero,zero_Vec,natural,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
  ierr = VecScatterEnd  (tozero,zero_Vec,natural,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
  ierr = VecScatterDestroy(tozero); CHKERRQ(ierr);
  ierr = VecDestroy(zero_Vec);CHKERRQ(ierr);
  
  ierr = DANaturalToGlobalBegin(da, natural, INSERT_VALUES, x);CHKERRQ(ierr);
  ierr = DANaturalToGlobalEnd  (da, natural, INSERT_VALUES, x);CHKERRQ(ierr);
  ierr = VecDestroy(natural);CHKERRQ(ierr);  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecPetscToGMRS"
/*
  VecPetscToGMRS: convert a distributed Petsc Vec into GMRS 1d array with all values on CPU 0

  (c) 2010 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode VecPetscToGMRS(Vec x,PetscReal *GMRS_array)
{
  Vec				     natural,zero_Vec;
  PetscReal      *zero_array;
  MPI_Comm       comm;
  VecScatter		 tozero;
  PetscMPIInt	   rank;
  DA				     da;
  PetscInt	     i,j,k,l,mx,my,mz,n,N;
  PetscErrorCode ierr;

  
  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject) x,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  
  ierr = PetscObjectQuery((PetscObject) x,"DA",(PetscObject *) &da);CHKERRQ(ierr);
  if (!da) SETERRQ(PETSC_ERR_ARG_WRONG,"Vector not generated from a DA");
  ierr = DAGetInfo(da,0,&mx,&my,&mz,0,0,0, 0,0,0,0);CHKERRQ(ierr);
  N = mx * my * mz;
  if (!rank) {
    n = N;
  } else {
    n = 0;
  }
  
  ierr = DACreateNaturalVector(da,&natural);CHKERRQ(ierr);
  ierr = DAGlobalToNaturalBegin(da,x,INSERT_VALUES,natural);CHKERRQ(ierr);
  ierr = DAGlobalToNaturalEnd  (da,x,INSERT_VALUES,natural);CHKERRQ(ierr);

  ierr = VecScatterCreateToZero(natural,&tozero,&zero_Vec);CHKERRQ(ierr);
  ierr = VecScatterBegin(tozero,natural,zero_Vec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd  (tozero,natural,zero_Vec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterDestroy(tozero); CHKERRQ(ierr);
  /*
    This would work if the dof ordering scheme for petsc and GMRs would coincide. 
    ierr = VecCreateMPIWithArray(comm,n,N,zero_array,&GMRS_Vec);CHKERRQ(ierr);
  */

  ierr = VecGetArray(zero_Vec,&zero_array);CHKERRQ(ierr);
  if (!rank) {
    for (l=0,i = 0; i < mx; i++) {
      for (j = 0; j < my; j++) {
        for (k = 0; k < mz; k++,l++) {
          GMRS_array[l] = zero_array[i + j*mx + k*mx*my];
        }
      }
    }
  }
  ierr = VecRestoreArray(zero_Vec,&zero_array);CHKERRQ(ierr);  
  ierr = VecDestroy(natural);CHKERRQ(ierr);  
  ierr = VecDestroy(zero_Vec);CHKERRQ(ierr);  
  PetscFunctionReturn(0);
}
