      MODULE TRANSIENTS
      CONTAINS
         DOUBLE PRECISION FUNCTION TEMP(THETAF,THETAR,R,L,X)
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(IN)  :: THETAF,THETAR,R,L,X
          TEMP=MIN(THETAR,MAX(THETAF,THETAF+(THETAR-THETAF)*(X-R+L)/L))
C          TEMP = THETAF
         END FUNCTION TEMP 
         
         DOUBLE PRECISION FUNCTION PRESSURE(PF,PR,R,L,X)
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(IN)  :: PF,PR,R,L,X

          DOUBLE PRECISION              :: PMIN, PMAX
          PMAX = MAX(PR,PF)
          PMIN = MIN(PR,PF)
          PRESSURE=MIN(PMAX,MAX(PMIN,PF+(PR-PF)/L*(X-R+L)))
C          PRESSURE = PF
         END FUNCTION PRESSURE 
      END MODULE TRANSIENTS
      
      PROGRAM FAKEGMRS
      USE TRANSIENTS
      INCLUDE 'mpif.h'
      
C     SAMPLE PARAMETERS SIMILAR TO THE BOUND EXAMPLE KEITA ONCE GAVE ME      
      INTEGER, PARAMETER          :: NDAYS = 88
      DOUBLE PRECISION, PARAMETER :: LX = 1000
      DOUBLE PRECISION, PARAMETER :: LY = 1000
      DOUBLE PRECISION, PARAMETER :: LZ = 200
      DOUBLE PRECISION, PARAMETER :: FLOWRATE = 20000000.0
      DOUBLE PRECISION, PARAMETER :: TEMP_W = 80.0
      DOUBLE PRECISION, PARAMETER :: TEMP_R = 180.0
      DOUBLE PRECISION, PARAMETER :: PRES_W = 3200.0
      DOUBLE PRECISION, PARAMETER :: PRES_R = 2400.0
      DOUBLE PRECISION, PARAMETER :: LDIFF = 300

      INTEGER, PARAMETER          :: NX = 51!51!75
      INTEGER, PARAMETER          :: NY = 51!51!40
      INTEGER, PARAMETER          :: NZ = 13!51!72
      INTEGER, PARAMETER          :: NSTEP = 11

      INTEGER                     :: IERR, RANK, NUMPROC
      DOUBLE PRECISION            :: DX(NX-1)
      DOUBLE PRECISION            :: DY(NY-1)
      DOUBLE PRECISION            :: DZ(NZ-1)
      INTEGER                     :: NFOUT,NFBUG
      REAL*8                      :: TIM
      INTEGER                     :: NEWT
      INTEGER                     :: N
      
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PRES,TEMPR,PMULT
      DOUBLE PRECISION                            :: TMPVAL
            
      INTEGER                     :: I,J,K,STEP
      DOUBLE PRECISION            :: FLUIDVOL, FLUIDRAD, R
      DOUBLE PRECISION            :: X(NX)
      DOUBLE PRECISION            :: Y(NY)
      DOUBLE PRECISION            :: Z(NZ)
      
      DOUBLE PRECISION            :: FMIN,FMAX
      
      CALL MPI_INIT(IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUMPROC, IERR)
      
C     DO SOME BOGUS INIT:
C      DX(1:NX/2-1) = 2 * LX / NX / 3.
C      DX(NX/2:NX)  = 4 * LX / NX / 3.
C
C      DY(1:NY/2-1)  = 2 * LY / NY / 3. 
C      DY(NY/2:NY) = 4 * LY / NY / 3.
C
C      DZ(1:NZ/4) = 4 * LZ / NZ / 3.
C      DZ(NZ/4+1:3*NZ/4-1) = 2 * LZ / NZ  / 3.
C      DZ(3*NZ/4:NZ) = 4 * LZ / NZ / 3.

      DX = LX / (NX-1)
      DY = LY / (NY-1)
      DZ = LZ / (NZ-1)
      
      X(1) = DX(1)/2.
      DO I = 1, NX-1
         X(I+1) = X(I) + DX(I+1)
      END DO

      Y(1) = DY(1)/2.
      DO J = 1, NY-1
         Y(J+1) = Y(J) + DY(J+1)
      END DO
      
      Z(1) = DZ(1)/2.
      DO K = 1, NZ-1
         Z(K+1) = Z(K) + DZ(K+1)
      END DO
        IF (RANK .EQ. 0) THEN
            WRITE(*,*) 'Reservoir extends: ', SUM(DX), SUM(DY), SUM(DZ)
            WRITE(*,*) '                X: ', X(1), X(NX)
            WRITE(*,*) '                Y: ', Y(1), Y(NY)
            WRITE(*,*) '                Z: ', Z(1), Z(NZ)
            WRITE(*,*) 'Reservoir volume: ', SUM(DX) * SUM(DY) * SUM(DZ)
      END IF
      
      NFOUT = 6
      NFBUG = 7

      CALL VIADAT(RANK,NUMPROC,NX,NY,NZ,DX,DY,DZ,NFOUT,NFBUG)
            
      IF (RANK .EQ. 0) THEN
         N = NX * NY * NZ
      ELSE
         N = 0
      END IF
      
      ALLOCATE(PRES(N))
      ALLOCATE(TEMPR(N))
      ALLOCATE(PMULT(N))

      TEMPR=0.
      PRES=0.
      DO STEP = 1, NSTEP
         TIM   = NDAYS * (STEP - 1.) / (NSTEP - 1.)
         FLUIDVOL = TIM * FLOWRATE
         FLUIDRAD = (FLUIDVOL * 3. / 4. / 3.14159) ** (1./3.)
         IF (RANK .EQ. 0) THEN
            WRITE(*,*) 'FLUIDRAD = ', FLUIDRAD
            DO I = 1, NX
               DO J = 1, NY
                  DO K = 1, NZ
                     R = SQRT(X(I)**2 + Y(J)**2 + 16.0*(Z(K)-LZ/2)**2)
                     TMPVAL = TEMP(TEMP_W,TEMP_R,FLUIDRAD,LDIFF,R)
                     TEMPR(K+(J-1)*NZ+(I-1)*NY*NZ)=TMPVAL
                     TMPVAL=PRESSURE(PRES_W,PRES_R,FLUIDRAD,LDIFF,R)
                     PRES(K+(J-1)*NZ+(I-1)*NY*NZ)=TMPVAL
                  END DO
               END DO
            END DO
         END IF
         NEWT = 1
         CALL VTDATA(RANK,TIM,STEP,NFOUT,NFBUG)
         CALL VPERM(RANK,PRES,TEMPR,PMULT,TIM,STEP,NEWT,NFOUT,NFBUG,N)
         CALL VSTDOUT(RANK,TIM,STEP,NFOUT,NFBUG)
      END DO
      
      DEALLOCATE(PRES)
      DEALLOCATE(TEMPR)
      DEALLOCATE(PMULT)
      
      CALL VENDIT(RANK,TIM,STEP,NFOUT,NFBUG)
      CALL MPI_FINALIZE(IERR)
      END PROGRAM FAKEGMRS
