#include "petscsys.h"
#include "petscfix.h"
#include "petscvec.h"
#include "petscda.h"


#ifdef PETSC_USE_POINTER_CONVERSION
#if defined(__cplusplus)
extern "C" { 
#endif 
extern void *PetscToPointer(void*);
extern int PetscFromPointer(void *);
extern void PetscRmPointer(void*);
#if defined(__cplusplus)
} 
#endif 

#else

#define PetscToPointer(a) (*(long *)(a))
#define PetscFromPointer(a) (long)(a)
#define PetscRmPointer(a)
#endif

#include "VF_CInterface.h"

#ifdef PETSC_HAVE_FORTRAN_CAPS
#define vfinitialize_ VFINITIALIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define vfinitialize_ vfinitialize
#endif

#ifdef PETSC_HAVE_FORTRAN_CAPS
#define vendit_ VENDIT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define vendit_ vendit
#endif

#ifdef PETSC_HAVE_FORTRAN_CAPS
#define vtdata_ VTDATA
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define vtdata_ vtdata
#endif

#ifdef PETSC_HAVE_FORTRAN_CAPS
#define vperm_ VPERM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define vperm_ vperm
#endif

#ifdef PETSC_HAVE_FORTRAN_CAPS
#define vstdout_ VSTDOUT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define vstdout_ vstdout
#endif

#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL vfinitialize_(PetscInt *nx,PetscInt *ny,PetscInt *nz,PetscReal *dx,PetscReal *dy,PetscReal *dz,PetscErrorCode *ierr){
  *ierr = VFInitialize(*nx,*ny,*nz,dx,dy,dz);
}
#if defined(__cplusplus)
}
#endif

#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL vendit_(PetscInt *rank,PetscReal *tim,PetscInt *nstep,PetscInt *nfout,PetscInt *nfbug){
  PetscErrorCode ierr;
  ierr = vendit(*rank,*tim,*nstep,*nfout,*nfbug);
}
#if defined(__cplusplus)
}
#endif

#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL vtdata_(PetscInt *rank,PetscReal *tim,PetscInt *nstep,int *nfout,int *nfbug){
  PetscErrorCode ierr;
  ierr = vtdata(*rank,*tim,*nstep,*nfout,*nfbug);
}
#if defined(__cplusplus)
}
#endif

#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL vperm_(PetscInt* rank,PetscReal *pres,PetscReal *tempr,PetscReal *pmult,PetscReal *tim,PetscInt *nstep,PetscInt *newt,PetscInt *nfout,PetscInt *nfbug,PetscInt *n){
  PetscErrorCode ierr;
  ierr = vperm(*rank,pres,tempr,pmult,*tim,*nstep,*newt,*nfout,*nfbug,*n);
}
#if defined(__cplusplus)
}
#endif

#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL vstdout_(PetscInt *rank,PetscReal *tim,PetscInt *nstep,int *nfout,int *nfbug){
  PetscErrorCode ierr;
  ierr = vstdout(*rank,*tim,*nstep,*nfout,*nfbug);
}
#if defined(__cplusplus)
}
#endif
