SUBROUTINE CART_GRIDHF
! PERFORM CARTESIAN GRID-BASED HF CALCULATION

   USE CONSTANTS
   USE CONTROL
   USE STRUCTURE
   USE INTEGRAL
   USE BASISSET
   USE DFT
   USE CARTESIANGRID

   IMPLICIT NONE
   INTEGER,PARAMETER :: MAXTRL = 100
   DOUBLE PRECISION,ALLOCATABLE :: TRL1(:,:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: PRD1(:,:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: TMP(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: H(:,:),VL(:,:),VR(:,:),ER(:),EI(:),WK(:)
   DOUBLE PRECISION,ALLOCATABLE :: R1(:)
   LOGICAL,ALLOCATABLE :: DONE(:)
   DOUBLE PRECISION :: RESIDUAL
   INTEGER :: NTRL
   INTEGER :: INFO
   CHARACTER(LEN=20) :: FILENAME
   INTEGER :: I,J
   DOUBLE PRECISION :: A
integer ix,iy

   MAXN=IOPTN(100)
   MAXN3=(2*MAXN+1)**3
   DELTAH=DOPTN(101)

   WRITE(6,'(A)') "CARTESIAN-GRID-BASED HF CALCULATION"
   WRITE(6,'(A,I3)') "NUMBER OF GRID POINTS = ",2*MAXN+1
   WRITE(6,'(A,F20.15,A)') "GRID SPACING = ",DELTAH, " BOHR"

   ALLOCATE(TRL1(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,MAXTRL))
   ALLOCATE(PRD1(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,MAXTRL))
   ALLOCATE(TMP(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN))
   ALLOCATE(OCC(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,IOCC))
   ALLOCATE(ORB(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,1:IALL(0)-IVIRTCORE))
   ALLOCATE(DEL2(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN))
   ALLOCATE(IPIV(MAXN3))
   ALLOCATE(H(MAXTRL,MAXTRL),VL(MAXTRL,MAXTRL),VR(MAXTRL,MAXTRL),ER(MAXTRL),EI(MAXTRL),WK(4*MAXTRL),R1(MAXTRL),DONE(MAXTRL))

write(*,*) "construct orbitals"
   CALL CART_ORBITALS
write(*,*) "construct J Del squared"
!  CALL CART_DEL2
write(*,*) "LU decomp J Del squared"
!  CALL DGETRF(MAXN3,MAXN3,DEL2,MAXN3,IPIV,INFO)
write(*,*) "INFO = ",INFO

   NTRL=IOCC
   DO I=1,IOCC
    TRL1(:,:,:,I)=ORB(:,:,:,I)
do ix=-maxn,maxn
 do iy=-maxn,maxn
  write(40+i,*) ix,iy,orb(ix,iy,0,i)
 enddo
enddo
    DONE(I)=.FALSE.
    CALL CART_NORMALIZE(TRL1(:,:,:,I))
    IF (I > 1) THEN
     DO J=1,I-1
      CALL CART_SCHMIDT(TRL1(:,:,:,I),TRL1(:,:,:,J))
     ENDDO
     CALL CART_NORMALIZE(TRL1(:,:,:,I))
    ENDIF
    OCC(:,:,:,I)=TRL1(:,:,:,I)
do ix=-maxn,maxn
 do iy=-maxn,maxn
  write(50+i,*) ix,iy,occ(ix,iy,0,i)
 enddo
enddo
   ENDDO
stop

   RESIDUAL=1.0D10
   DO WHILE(RESIDUAL > DOPTN(67))
    DO I=1,NTRL
     IF (.NOT.DONE(I)) THEN
      CALL CART_ZERO(PRD1(:,:,:,I))
      CALL CART_T(TRL1(:,:,:,I),PRD1(:,:,:,I))
      CALL CART_N(TRL1(:,:,:,I),PRD1(:,:,:,I))
      CALL CART_J(TRL1(:,:,:,I),PRD1(:,:,:,I))
      CALL CART_K(TRL1(:,:,:,I),PRD1(:,:,:,I))
      DONE(I)=.TRUE.
     ENDIF
    ENDDO
do ix=-maxn,maxn
 do iy=-maxn,maxn
  write(60+i,*) ix,iy,prd1(ix,iy,0,i)
 enddo
enddo
    DO I=1,NTRL
     DO J=1,NTRL
      CALL CART_PRD(TRL1(:,:,:,I),PRD1(:,:,:,J),H(I,J))
     ENDDO
    ENDDO
call dump10(h,ntrl,maxtrl)
    CALL DGEEV('N','V',NTRL,H,MAXTRL,ER,EI,VL,MAXTRL,VR,MAXTRL,WK,4*MAXTRL,INFO)
write(*,*) "info=",info
write(*,*) (er(i),i=1,ntrl)
    DO I=1,IOCC
     CALL CART_LC(TRL1,OCC(:,:,:,I),VR(:,I),NTRL)
     CALL CART_LC(PRD1,TMP,VR(:,I),NTRL)
     TMP=TMP-OCC(:,:,:,I)*ER(I)
do ix=-maxn,maxn
 do iy=-maxn,maxn
  write(70+i,*) ix,iy,occ(ix,iy,0,i),tmp(ix,iy,0)
 enddo
enddo
do ix=-maxn,maxn
 do iy=-maxn,maxn
  write(80+i,*) ix,iy,occ(ix,iy,0,i),tmp(ix,iy,0)
 enddo
enddo
     CALL CART_PRD(TMP,TMP,A)
write(*,*) "residual =",a
    ENDDO
stop
   ENDDO

   DEALLOCATE(DONE,R1,WK,EI,ER,VR,VL,H,IPIV,DEL2,ORB,OCC,TMP,PRD1,TRL1)

   RETURN
END SUBROUTINE



SUBROUTINE CART_SCHMIDT(X,Y)
! SCHMIDT ORTHOGONALIZE X AGAINST Y

   USE CARTESIANGRID

   IMPLICIT NONE
   DOUBLE PRECISION :: X(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   DOUBLE PRECISION :: Y(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   INTEGER :: IX,IY,IZ
   DOUBLE PRECISION :: A

   A=0.0D0
   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      A=A+Y(IX,IY,IZ)*Y(IX,IY,IZ)*DELTAH**3
     ENDDO
    ENDDO
   ENDDO
   A=0.0D0
   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      A=A+X(IX,IY,IZ)*Y(IX,IY,IZ)*DELTAH**3
     ENDDO
    ENDDO
   ENDDO
   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      X(IX,IY,IZ)=X(IX,IY,IZ)-A*Y(IX,IY,IZ)
     ENDDO
    ENDDO
   ENDDO

   RETURN
END SUBROUTINE



SUBROUTINE CART_NORMALIZE(X)
! NORMALIZE X

   USE CARTESIANGRID

   IMPLICIT NONE
   DOUBLE PRECISION :: X(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   INTEGER :: IX,IY,IZ
   DOUBLE PRECISION :: A

   A=0.0D0
   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      A=A+X(IX,IY,IZ)*X(IX,IY,IZ)*DELTAH**3
     ENDDO
    ENDDO
   ENDDO
   IF (A==0.0D0) CALL PABORT('NULL TRIAL VECTOR')
   A=1.0D0/DSQRT(A)
   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      X(IX,IY,IZ)=X(IX,IY,IZ)*A
     ENDDO
    ENDDO
   ENDDO

   RETURN
END SUBROUTINE



SUBROUTINE CART_ZERO(X)
! ZERO SCRATCH X

   USE CARTESIANGRID

   IMPLICIT NONE
   DOUBLE PRECISION :: X(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   INTEGER :: IX,IY,IZ

   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      X(IX,IY,IZ)=0.0D0
     ENDDO
    ENDDO
   ENDDO

   RETURN
END SUBROUTINE



SUBROUTINE CART_T(X,Y)
! ACT THE KINETIC-ENERGY OPERATOR ON X AND ADD THE RESULT TO Y

   USE CARTESIANGRID

   IMPLICIT NONE
   DOUBLE PRECISION :: X(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   DOUBLE PRECISION :: Y(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   INTEGER :: IX,IY,IZ

   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      IF ((IX >= -MAXN+2).AND.(IX <= MAXN-2)) THEN
       Y(IX,IY,IZ)=Y(IX,IY,IZ)-(0.5D0/12.0D0/DELTAH**2)*( &
                            -1.0D0*X(IX-2,IY,IZ) &
                           +16.0D0*X(IX-1,IY,IZ) &
                           -30.0D0*X(IX,  IY,IZ) &
                           +16.0D0*X(IX+1,IY,IZ) &
                            -1.0D0*X(IX+2,IY,IZ))
      ELSE IF (IX == -MAXN+1) THEN
       Y(IX,IY,IZ)=Y(IX,IY,IZ)-(0.5D0/12.0D0/DELTAH**2)*( &
                           +11.0D0*X(IX-1,IY,IZ) &
                           -20.0D0*X(IX,  IY,IZ) &
                            +6.0D0*X(IX+1,IY,IZ) &
                            +4.0D0*X(IX+2,IY,IZ) &
                            -1.0D0*X(IX+3,IY,IZ))
      ELSE IF (IX == MAXN-1) THEN
       Y(IX,IY,IZ)=Y(IX,IY,IZ)-(0.5D0/12.0D0/DELTAH**2)*( &
                            -1.0D0*X(IX-3,IY,IZ) &
                            +4.0D0*X(IX-2,IY,IZ) &
                            +6.0D0*X(IX-1,IY,IZ) &
                           -20.0D0*X(IX,  IY,IZ) &
                           +11.0D0*X(IX+1,IY,IZ))
      ELSE IF (IX == -MAXN) THEN
       Y(IX,IY,IZ)=Y(IX,IY,IZ)-(0.5D0/12.0D0/DELTAH**2)*( &
                           +35.0D0*X(IX,  IY,IZ) &
                          -104.0D0*X(IX+1,IY,IZ) &
                          +114.0D0*X(IX+2,IY,IZ) &
                           -56.0D0*X(IX+3,IY,IZ) &
                           +11.0D0*X(IX+4,IY,IZ))
      ELSE IF (IX == MAXN) THEN
       Y(IX,IY,IZ)=Y(IX,IY,IZ)-(0.5D0/12.0D0/DELTAH**2)*( &
                           +11.0D0*X(IX-4,IY,IZ) &
                           -56.0D0*X(IX-3,IY,IZ) &
                          +114.0D0*X(IX-2,IY,IZ) &
                          -104.0D0*X(IX-1,IY,IZ) &
                           +35.0D0*X(IX,  IY,IZ))
      ENDIF
      IF ((IY >= -MAXN+2).AND.(IY <= MAXN-2)) THEN
       Y(IX,IY,IZ)=Y(IX,IY,IZ)-(0.5D0/12.0D0/DELTAH**2)*( &
                            -1.0D0*X(IX,IY-2,IZ) &
                           +16.0D0*X(IX,IY-1,IZ) &
                           -30.0D0*X(IX,IY,  IZ) &
                           +16.0D0*X(IX,IY+1,IZ) &
                            -1.0D0*X(IX,IY+2,IZ))
      ELSE IF (IY == -MAXN+1) THEN
       Y(IX,IY,IZ)=Y(IX,IY,IZ)-(0.5D0/12.0D0/DELTAH**2)*( &
                           +11.0D0*X(IX,IY-1,IZ) &
                           -20.0D0*X(IX,IY,  IZ) &
                            +6.0D0*X(IX,IY+1,IZ) &
                            +4.0D0*X(IX,IY+2,IZ) &
                            -1.0D0*X(IX,IY+3,IZ))
      ELSE IF (IY == MAXN-1) THEN
       Y(IX,IY,IZ)=Y(IX,IY,IZ)-(0.5D0/12.0D0/DELTAH**2)*( &
                            -1.0D0*X(IX,IY-3,IZ) &
                            +4.0D0*X(IX,IY-2,IZ) &
                            +6.0D0*X(IX,IY-1,IZ) &
                           -20.0D0*X(IX,IY,  IZ) &
                           +11.0D0*X(IX,IY+1,IZ))
      ELSE IF (IY == -MAXN) THEN
       Y(IX,IY,IZ)=Y(IX,IY,IZ)-(0.5D0/12.0D0/DELTAH**2)*( &
                           +35.0D0*X(IX,IY,  IZ) &
                          -104.0D0*X(IX,IY+1,IZ) &
                          +114.0D0*X(IX,IY+2,IZ) &
                           -56.0D0*X(IX,IY+3,IZ) &
                           +11.0D0*X(IX,IY+4,IZ))
      ELSE IF (IY == MAXN) THEN
       Y(IX,IY,IZ)=Y(IX,IY,IZ)-(0.5D0/12.0D0/DELTAH**2)*( &
                           +11.0D0*X(IX,IY-4,IZ) &
                           -56.0D0*X(IX,IY-3,IZ) &
                          +114.0D0*X(IX,IY-2,IZ) &
                          -104.0D0*X(IX,IY-1,IZ) &
                           +35.0D0*X(IX,IY,  IZ))
      ENDIF
      IF ((IZ >= -MAXN+2).AND.(IZ <= MAXN-2)) THEN
       Y(IX,IY,IZ)=Y(IX,IY,IZ)-(0.5D0/12.0D0/DELTAH**2)*( &
                            -1.0D0*X(IX,IY,IZ-2) &
                           +16.0D0*X(IX,IY,IZ-1) &
                           -30.0D0*X(IX,IY,IZ  ) &
                           +16.0D0*X(IX,IY,IZ+1) &
                            -1.0D0*X(IX,IY,IZ+2))
      ELSE IF (IZ == -MAXN+1) THEN
       Y(IX,IY,IZ)=Y(IX,IY,IZ)-(0.5D0/12.0D0/DELTAH**2)*( &
                           +11.0D0*X(IX,IY,IZ-1) &
                           -20.0D0*X(IX,IY,IZ  ) &
                            +6.0D0*X(IX,IY,IZ+1) &
                            +4.0D0*X(IX,IY,IZ+2) &
                            -1.0D0*X(IX,IY,IZ+3))
      ELSE IF (IZ == MAXN-1) THEN
       Y(IX,IY,IZ)=Y(IX,IY,IZ)-(0.5D0/12.0D0/DELTAH**2)*( &
                            -1.0D0*X(IX,IY,IZ-3) &
                            +4.0D0*X(IX,IY,IZ-2) &
                            +6.0D0*X(IX,IY,IZ-1) &
                           -20.0D0*X(IX,IY,IZ  ) &
                           +11.0D0*X(IX,IY,IZ+1))
      ELSE IF (IZ == -MAXN) THEN
       Y(IX,IY,IZ)=Y(IX,IY,IZ)-(0.5D0/12.0D0/DELTAH**2)*( &
                           +35.0D0*X(IX,IY,IZ  ) &
                          -104.0D0*X(IX,IY,IZ+1) &
                          +114.0D0*X(IX,IY,IZ+2) &
                           -56.0D0*X(IX,IY,IZ+3) &
                           +11.0D0*X(IX,IY,IZ+4))
      ELSE IF (IZ == MAXN) THEN
       Y(IX,IY,IZ)=Y(IX,IY,IZ)-(0.5D0/12.0D0/DELTAH**2)*( &
                           +11.0D0*X(IX,IY,IZ-4) &
                           -56.0D0*X(IX,IY,IZ-3) &
                          +114.0D0*X(IX,IY,IZ-2) &
                          -104.0D0*X(IX,IY,IZ-1) &
                           +35.0D0*X(IX,IY,IZ  ))
      ENDIF
     ENDDO
    ENDDO
   ENDDO

   RETURN
END SUBROUTINE



SUBROUTINE CART_N(X,Y)
! ACT THE NUCLEAR-ATTRACTION POTENTIAL ON X AND ADD THE RESULT TO Y

   USE STRUCTURE
   USE INTEGRAL
   USE CARTESIANGRID

   IMPLICIT NONE
   DOUBLE PRECISION :: X(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   DOUBLE PRECISION :: Y(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   INTEGER :: IX,IY,IZ
   DOUBLE PRECISION :: X1,Y1,Z1,X2,Y2,Z2,R
   INTEGER :: IA

   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      X1=DFLOAT(IX)*DELTAH
      Y1=DFLOAT(IY)*DELTAH
      Z1=DFLOAT(IZ)*DELTAH
      DO IA=1,NATOM
       X2=ATOMX(IA)
       Y2=ATOMY(IA)
       Z2=ATOMZ(IA)
       R=(X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2
       IF (R > 0.0D0) R=1.0D0/SQRT(R)
       Y(IX,IY,IZ)=Y(IX,IY,IZ)-DFLOAT(IATOM(IA))*R*X(IX,IY,IZ)
      ENDDO
     ENDDO
    ENDDO
   ENDDO

   RETURN
END SUBROUTINE



SUBROUTINE CART_J(X,Y)
! ACT THE COULOMB POTENTIAL ON X AND ADD THE RESULT TO Y

   USE CONSTANTS
   USE STRUCTURE
   USE INTEGRAL
   USE CARTESIANGRID

   IMPLICIT NONE
   DOUBLE PRECISION :: X(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   DOUBLE PRECISION :: Y(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   DOUBLE PRECISION :: Z(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   INTEGER :: IX,IY,IZ
   DOUBLE PRECISION :: X1,Y1,Z1,R,C
   INTEGER :: IA
   INTEGER :: INFO

   C=0.0D0
   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      DO IA=1,IOCC
       C=C+2.0D0*OCC(IX,IY,IZ,IA)**2*DELTAH**3
      ENDDO
     ENDDO
    ENDDO
   ENDDO
!  IF (C /= DFLOAT(2*IOCC)) CALL WARNING("ERROR IN COULOMB CHARGE")

   Z=0.0D0
   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      DO IA=1,IOCC
       Z(IX,IY,IZ)=Z(IX,IY,IZ)-8.0D0*PI*OCC(IX,IY,IZ,IA)**2
      ENDDO
      X1=DFLOAT(IX)*DELTAH
      Y1=DFLOAT(IY)*DELTAH
      Z1=DFLOAT(IZ)*DELTAH
      R=X1**2+Y1**2+Z1**2
      IF (R > 0.0D0) R=C/DSQRT(R)
      IF (IX == -MAXN+1) Z(IX,IY,IZ)=Z(IX,IY,IZ) +1.0D0/12.0D0/DELTAH**2*R
      IF (IX ==  MAXN-1) Z(IX,IY,IZ)=Z(IX,IY,IZ) +1.0D0/12.0D0/DELTAH**2*R
      IF (IX == -MAXN  ) Z(IX,IY,IZ)=Z(IX,IY,IZ)-11.0D0/12.0D0/DELTAH**2*R
      IF (IX ==  MAXN  ) Z(IX,IY,IZ)=Z(IX,IY,IZ)-11.0D0/12.0D0/DELTAH**2*R
      IF (IY == -MAXN+1) Z(IX,IY,IZ)=Z(IX,IY,IZ) +1.0D0/12.0D0/DELTAH**2*R
      IF (IY ==  MAXN-1) Z(IX,IY,IZ)=Z(IX,IY,IZ) +1.0D0/12.0D0/DELTAH**2*R
      IF (IY == -MAXN  ) Z(IX,IY,IZ)=Z(IX,IY,IZ)-11.0D0/12.0D0/DELTAH**2*R
      IF (IY ==  MAXN  ) Z(IX,IY,IZ)=Z(IX,IY,IZ)-11.0D0/12.0D0/DELTAH**2*R
      IF (IZ == -MAXN+1) Z(IX,IY,IZ)=Z(IX,IY,IZ) +1.0D0/12.0D0/DELTAH**2*R
      IF (IZ ==  MAXN-1) Z(IX,IY,IZ)=Z(IX,IY,IZ) +1.0D0/12.0D0/DELTAH**2*R
      IF (IZ == -MAXN  ) Z(IX,IY,IZ)=Z(IX,IY,IZ)-11.0D0/12.0D0/DELTAH**2*R
      IF (IZ ==  MAXN  ) Z(IX,IY,IZ)=Z(IX,IY,IZ)-11.0D0/12.0D0/DELTAH**2*R
     ENDDO
    ENDDO
   ENDDO

   CALL DGETRS('N',MAXN3,1,DEL2,MAXN3,IPIV,Z,MAXN3,INFO)

   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      Y(IX,IY,IZ)=Y(IX,IY,IZ)+Z(IX,IY,IZ)*X(IX,IY,IZ)
     ENDDO
    ENDDO
   ENDDO

   RETURN
END SUBROUTINE




SUBROUTINE CART_K(X,Y)
! ACT THE EXCHANGE OPERATOR ON X AND ADD THE RESULT TO Y

   USE CONSTANTS
   USE STRUCTURE
   USE INTEGRAL
   USE CARTESIANGRID

   IMPLICIT NONE
   DOUBLE PRECISION :: X(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   DOUBLE PRECISION :: Y(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   DOUBLE PRECISION :: Z(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   INTEGER :: IX,IY,IZ
   DOUBLE PRECISION :: X1,Y1,Z1,R,C
   INTEGER :: IA
   INTEGER :: INFO

   DO IA=1,IOCC

    C=0.0D0
    DO IX=-MAXN,MAXN
     DO IY=-MAXN,MAXN
      DO IZ=-MAXN,MAXN
       C=C+OCC(IX,IY,IZ,IA)*X(IX,IY,IZ)*DELTAH**3
      ENDDO
     ENDDO
    ENDDO

    Z=0.0D0
    DO IX=-MAXN,MAXN
     DO IY=-MAXN,MAXN
      DO IZ=-MAXN,MAXN
       Z(IX,IY,IZ)=Z(IX,IY,IZ)-4.0D0*PI*OCC(IX,IY,IZ,IA)*X(IX,IY,IZ)
       X1=DFLOAT(IX)*DELTAH
       Y1=DFLOAT(IY)*DELTAH
       Z1=DFLOAT(IZ)*DELTAH
       R=X1**2+Y1**2+Z1**2
       IF (R > 0.0D0) R=C/DSQRT(R)
       IF (IX == -MAXN+1) Z(IX,IY,IZ)=Z(IX,IY,IZ) +1.0D0/12.0D0/DELTAH**2*R
       IF (IX ==  MAXN-1) Z(IX,IY,IZ)=Z(IX,IY,IZ) +1.0D0/12.0D0/DELTAH**2*R
       IF (IX == -MAXN  ) Z(IX,IY,IZ)=Z(IX,IY,IZ)-11.0D0/12.0D0/DELTAH**2*R
       IF (IX ==  MAXN  ) Z(IX,IY,IZ)=Z(IX,IY,IZ)-11.0D0/12.0D0/DELTAH**2*R
       IF (IY == -MAXN+1) Z(IX,IY,IZ)=Z(IX,IY,IZ) +1.0D0/12.0D0/DELTAH**2*R
       IF (IY ==  MAXN-1) Z(IX,IY,IZ)=Z(IX,IY,IZ) +1.0D0/12.0D0/DELTAH**2*R
       IF (IY == -MAXN  ) Z(IX,IY,IZ)=Z(IX,IY,IZ)-11.0D0/12.0D0/DELTAH**2*R
       IF (IY ==  MAXN  ) Z(IX,IY,IZ)=Z(IX,IY,IZ)-11.0D0/12.0D0/DELTAH**2*R
       IF (IZ == -MAXN+1) Z(IX,IY,IZ)=Z(IX,IY,IZ) +1.0D0/12.0D0/DELTAH**2*R
       IF (IZ ==  MAXN-1) Z(IX,IY,IZ)=Z(IX,IY,IZ) +1.0D0/12.0D0/DELTAH**2*R
       IF (IZ == -MAXN  ) Z(IX,IY,IZ)=Z(IX,IY,IZ)-11.0D0/12.0D0/DELTAH**2*R
       IF (IZ ==  MAXN  ) Z(IX,IY,IZ)=Z(IX,IY,IZ)-11.0D0/12.0D0/DELTAH**2*R
      ENDDO
     ENDDO
    ENDDO
    CALL DGETRS('N',MAXN3,1,DEL2,MAXN3,IPIV,Z,MAXN3,INFO)
    DO IX=-MAXN,MAXN
     DO IY=-MAXN,MAXN
      DO IZ=-MAXN,MAXN
       Y(IX,IY,IZ)=Y(IX,IY,IZ)-Z(IX,IY,IZ)*X(IX,IY,IZ)
      ENDDO
     ENDDO
    ENDDO

   ENDDO

   RETURN
END SUBROUTINE




SUBROUTINE CART_PRD(X,Y,A)
! RETURN AS A THE VECTOR INNER PRODUCT OF X AND Y

   USE CARTESIANGRID

   IMPLICIT NONE
   DOUBLE PRECISION :: X(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   DOUBLE PRECISION :: Y(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   INTEGER :: IX,IY,IZ
   DOUBLE PRECISION :: A

   A=0.0D0
   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      A=A+X(IX,IY,IZ)*Y(IX,IY,IZ)*DELTAH**3
     ENDDO
    ENDDO
   ENDDO

   RETURN
END SUBROUTINE




SUBROUTINE CART_LC(X,Y,C,N)
! RETURN A LINEAR COMBINATION, Y = C X

   USE CARTESIANGRID

   IMPLICIT NONE
   INTEGER :: N
   DOUBLE PRECISION :: X(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,N)
   DOUBLE PRECISION :: Y(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   DOUBLE PRECISION :: C(N)
   INTEGER :: IX,IY,IZ,I

   Y=0.0D0
   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      DO I=1,N
       Y(IX,IY,IZ)=Y(IX,IY,IZ)+X(IX,IY,IZ,N)*C(N)
      ENDDO
     ENDDO
    ENDDO
   ENDDO

   RETURN
END SUBROUTINE
