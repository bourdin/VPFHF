Subroutine VIADAT(myprc,numprc,nx,ny,nz,dx,dy,dz,nfout,nfbug)
!*********************************************************************
!  THIS ROUTINE IS CALLED AT INITIALIZATION BY ALL PROCESSORS
!  MYPRC = PROCESSOR NUMBER (MASTER PROCESSOR IS 0)
!  NUMPRC = NUMBER OF PROCESSORS FOR PARALLEL RUN
!  NX = NUMBER OF CELLS IN X-DIRECTION
!  NY = NUMBER OF CELLS IN Y-DIRECTION
!  NZ = NUMBER OF CELLS IN Z-DIRECTION
!  DX = GRID SPACING IN X-DIRECTION (FT)
!  DY = GRID SPACING IN Y-DIRECTION (FT)
!  DZ = GRID SPACING IN Z-DIRECTION (FT)  
!  NFOUT = UNIT NUMBER FOR PRINTING ON MYPRC = 0
!  NFBUG = UNIT NUMBER FOR PRINTING FOR PROCESSOR MYPRC /= 0
!          REQUIRES BUGKEY(I)=TRUE FOR SOME I /= 6 OR 7
!          OUTPUT FILE WILL BE debug##.JOBID.out WHERE ## IS PROCESSOR NUMBER
!          AND JOBID IS DETERMINED BY THE SYSTEM AT RUNTIME
!*********************************************************************
#include "finclude/petscdef.h"
   Use petsc
   Implicit NONE
   Integer, Intent(IN) :: myprc,numprc,nx,ny,nz,nfout,nfbug
   Real*8, Intent(IN)  :: dx(nx),dy(ny),dz(nz)
   PetscInt            :: ierr
   integer             :: junk

   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)
   
   Call VFInitialize(nx,ny,nz,dx,dy,dz,ierr);
   If (ierr /= 0) Then
      Call MPI_Finalize(ierr)
      STOP
   End If
   ! Getting rid of an annoying warning.
   junk=myprc
   junk=numprc
   junk=nfout
   junk=nfbug
End Subroutine VIADAT

