SUBROUTINE CART_GRIDMP2
! PERFORM CARTESIAN GRID-BASED MP2-R12 CALCULATION BY SOLVING SINANOGLU EQUATION.

   USE CONSTANTS
   USE CONTROL
   USE STRUCTURE
   USE INTEGRAL
   USE BASISSET
   USE DFT
   USE CARTESIANGRID

   IMPLICIT NONE
   INTEGER :: I,J
   DOUBLE PRECISION,ALLOCATABLE :: RHS(:,:,:,:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: TRL1(:,:,:,:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: TRL2(:,:,:,:,:,:)
   INTEGER :: NTRL
   INTEGER :: INFO
   CHARACTER(LEN=20) :: FILENAME
   DOUBLE PRECISION :: RESIDUAL  ! RESIDUAL; 'CICONV' DETERMINES THE CONVERGENCE CRITERION
   DOUBLE PRECISION :: E_SINGLET ! SINGLET PAIR ENERGY
integer :: jx,jy

   MAXN=IOPTN(100)
   MAXN3=(2*MAXN+1)**3
   DELTAH=DOPTN(101)

   WRITE(6,'(A)') "CARTESIAN-GRID-BASED MP2 CALCULATION"
   WRITE(6,'(A,I3)') "NUMBER OF GRID POINTS = ",2*MAXN+1
   WRITE(6,'(A,F20.15,A)') "GRID SPACING = ",DELTAH, " BOHR"

   ALLOCATE(RHS(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN))
   ALLOCATE(TRL1(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN))
   ALLOCATE(TRL2(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN))
   ALLOCATE(ORB(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,1:IALL(0)-IVIRTCORE))
   ALLOCATE(DEL2(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN))
   ALLOCATE(IPIV(MAXN3))

write(*,*) "construct orbitals"
   CALL CART_ORBITALS
write(*,*) "construct J Del squared"
   CALL CART_DEL2
write(*,*) "LU decomp J Del squared"
   CALL DGETRF(MAXN3,MAXN3,DEL2,MAXN3,IPIV,INFO)
write(*,*) "INFO = ",INFO

write(*,*) "singlet pairs"
   DO I=1,IOCC
    DO J=I,IOCC
write(*,*) "occ-occ pair",i,j
write(*,*) "construct rhs"
     NTRL=0
     RESIDUAL=1.0D10
     CALL CART_SINGLET(RHS,I,J,I,J,.TRUE.)
     CALL PROJECTOR(RHS,.TRUE.)
write(*,*) "an initial guess"
     CALL CART_SINGLET(TRL1,I,J,IOCC+1,IOCC+1,.FALSE.)

     DO WHILE(RESIDUAL > DOPTN(67))
      NTRL=NTRL+1
      IF (NTRL > MAXNTRL) CALL PABORT('MAXNTRL REACHED')
      WRITE(FILENAME,'(A,".",I2.2)') TRIM(COPTN(1)),NTRL
write(*,*) 100+NTRL,filename
      OPEN(100+NTRL,FILE=FILENAME,FORM='UNFORMATTED')
      REWIND(100+NTRL)
      WRITE(100+NTRL) TRL1
do jx=-maxn,maxn
 do jy=-maxn,maxn
  write(60,*) jx,jy,trl1(jx,jy,0,0,0,0)
 enddo
enddo
      TRL2=0.0D0
      CALL CART_FOCK(TRL1,TRL2,I,J,.TRUE.)
do jx=-maxn,maxn
 do jy=-maxn,maxn
  write(61,*) jx,jy,trl2(jx,jy,0,0,0,0)
 enddo
enddo
      CALL CART_KINETIC(TRL1,TRL2,I,J,.TRUE.)
do jx=-maxn,maxn
 do jy=-maxn,maxn
  write(62,*) jx,jy,trl2(jx,jy,0,0,0,0)
 enddo
enddo
stop
      CALL PROJECTOR(TRL2,.TRUE.)
      REWIND(200+NTRL)
      WRITE(200+NTRL) TRL2
      CALL CART_RESIDUAL(I,J,RHS,TRL2,RESIDUAL)
write(*,*) "residual = ",residual
     ENDDO

    ENDDO
   ENDDO

   DEALLOCATE(IPIV,DEL2,ORB,TRL2,TRL1,RHS)

   RETURN
END SUBROUTINE



SUBROUTINE CART_SINGLET(X,I,J,A,B,TOGGLE)
! CONSTRUCT A SINGLET PAIR FUNCTION ON THE CARTESIAN GRID.

   USE STRUCTURE
   USE BASISSET
   USE INTEGRAL
   USE CARTESIANGRID

   IMPLICIT NONE
   DOUBLE PRECISION :: X(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   INTEGER :: I,J,A,B
   LOGICAL :: TOGGLE ! TRUE FOR 1/R12; FALSE FOR 1/(EI+EJ-EA-EB)
   INTEGER :: IX,IY,IZ,JX,JY,JZ
   DOUBLE PRECISION :: X1,Y1,Z1,X2,Y2,Z2,R

   X=0.0D0
   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      X1=DFLOAT(IX)*DELTAH
      Y1=DFLOAT(IY)*DELTAH
      Z1=DFLOAT(IZ)*DELTAH
      DO JX=-MAXN,MAXN
       DO JY=-MAXN,MAXN
        DO JZ=-MAXN,MAXN
         X2=DFLOAT(JX)*DELTAH
         Y2=DFLOAT(JY)*DELTAH
         Z2=DFLOAT(JZ)*DELTAH
         R=(X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2
         IF (TOGGLE) THEN
          IF (R > 0.0D0) R=1.0D0/SQRT(R)
          X(JX,JY,JZ,IX,IY,IZ)=R*(ORB(JX,JY,JZ,I)*ORB(IX,IY,IZ,J)+ORB(JX,JY,JZ,J)*ORB(IX,IY,IZ,I))
         ELSE
!         IF (R > 0.0D0) R=1.0D0/SQRT(R)
          IF (R > 0.0D0) R=1.0D0
          X(JX,JY,JZ,IX,IY,IZ)=R*(ORB(JX,JY,JZ,A)*ORB(IX,IY,IZ,A)+ORB(JX,JY,JZ,B)*ORB(IX,IY,IZ,A)) &
                              /(EPSILON(I,0)+EPSILON(J,0)-EPSILON(A,0)-EPSILON(B,0))
! debug
X(JX,JY,JZ,IX,IY,IZ)=ORB(JX,JY,JZ,2)
         ENDIF
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO

   RETURN
END SUBROUTINE



SUBROUTINE PROJECTOR(X,TOGGLE)
! ACT A PROJECTOR THAT ELIMINATES OCCUPIED ORBITAL COMPONENTS

   USE STRUCTURE
   USE BASISSET
   USE CARTESIANGRID

   IMPLICIT NONE
   DOUBLE PRECISION :: X(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   LOGICAL :: TOGGLE
   INTEGER :: IX,IY,IZ,JX,JY,JZ
   INTEGER :: K
   DOUBLE PRECISION :: A,SGN

   IF (TOGGLE) THEN
    SGN=1.0D0
   ELSE
    SGN=-1.0D0
   ENDIF

   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      DO K=1,IOCC
       A=0.0D0
       DO JX=-MAXN,MAXN
        DO JY=-MAXN,MAXN
         DO JZ=-MAXN,MAXN
          A=A+X(JX,JY,JZ,IX,IY,IZ)*ORB(JX,JY,JZ,K)
         ENDDO
        ENDDO
       ENDDO
       A=A*DELTAH**3
       DO JX=-MAXN,MAXN
        DO JY=-MAXN,MAXN
         DO JZ=-MAXN,MAXN
          X(JX,JY,JZ,IX,IY,IZ)=X(JX,JY,JZ,IX,IY,IZ)-A*ORB(JX,JY,JZ,K)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO

   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      DO JX=-MAXN,MAXN
       DO JY=-MAXN,MAXN
        DO JZ=-MAXN,MAXN
         X(IX,IY,IZ,JX,JY,JZ)=SGN*X(JX,JY,JZ,IX,IY,IZ)
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO


!  DO IX=-MAXN,MAXN
!   DO IY=-MAXN,MAXN
!    DO IZ=-MAXN,MAXN
!     DO K=1,IOCC
!      A=0.0D0
!      DO JX=-MAXN,MAXN
!       DO JY=-MAXN,MAXN
!        DO JZ=-MAXN,MAXN
!         A=A+X(IX,IY,IZ,JX,JY,JZ)*ORB(JX,JY,JZ,K)
!        ENDDO
!       ENDDO
!      ENDDO
!      A=A*DELTAH**3
!      DO JX=-MAXN,MAXN
!       DO JY=-MAXN,MAXN
!        DO JZ=-MAXN,MAXN
!         X(IX,IY,IZ,JX,JY,JZ)=X(IX,IY,IZ,JX,JY,JZ)-A*ORB(JX,JY,JZ,K)
!        ENDDO
!       ENDDO
!      ENDDO
!     ENDDO
!    ENDDO
!   ENDDO
!  ENDDO
   
   RETURN
END SUBROUTINE



SUBROUTINE CART_ORBITALS
! CALCULATE THE AMPLITUDE OF ORBITALS AT CARTESIAN GRID POINTS.

   USE STRUCTURE
   USE BASISSET
   USE INTEGRAL
   USE CARTESIANGRID

   IMPLICIT NONE
   INTEGER :: N,M,L
   INTEGER :: C2P(NPGS)
   INTEGER :: IX,IY,IZ,I,J,K
   DOUBLE PRECISION :: X1,Y1,Z1,X2,Y2,Z2,R,P
   DOUBLE PRECISION :: CGS(NCGS)

   ORB=0.0D0
   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      X1=DFLOAT(IX)*DELTAH
      Y1=DFLOAT(IY)*DELTAH
      Z1=DFLOAT(IZ)*DELTAH
      CGS=0.0D0
      DO I=1,NCGS
       K=0
       DO J=1,NPGS
        IF (CC(I,J) /= 0.0D0) THEN
         K=K+1
         C2P(K)=J
        ENDIF
       ENDDO
       X2=X1-CGSX(I)
       Y2=Y1-CGSY(I)
       Z2=Z1-CGSZ(I)
       R=X2**2+Y2**2+Z2**2
       P=1.0D0
       IF (PAX(C2P(1)) == 1) THEN
        P=P*X2
       ELSE IF (PAX(C2P(1)) == 2) THEN
        P=P*X2*X2
       ENDIF
       IF (PAY(C2P(1)) == 1) THEN
        P=P*Y2
       ELSE IF (PAY(C2P(1)) == 2) THEN
        P=P*Y2*Y2
       ENDIF
       IF (PAZ(C2P(1)) == 1) THEN
        P=P*Z2
       ELSE IF (PAZ(C2P(1)) == 2) THEN
        P=P*Z2*Z2
       ENDIF
       DO J=1,K
        CGS(I)=CGS(I)+DEXP(-ZT(C2P(J))*R)*CC(I,C2P(J))*P
       ENDDO
      ENDDO
      DO I=1,IALL(0)-IVIRTCORE
       DO J=1,NCGS
        ORB(IX,IY,IZ,I)=ORB(IX,IY,IZ,I)+DREAL(CO(J,I,0))*CGS(J)
       ENDDO
      ENDDO
!if (iz==0) then
!do i=1,ncgs
!write(50+i,*) ix,iy,cgs(i)
!enddo
!do i=1,iall(0)-ivirtcore
!write(40+i,*) ix,iy,orb(ix,iy,iz,i)
!enddo
!endif
     ENDDO
    ENDDO
   ENDDO

   RETURN
END SUBROUTINE



SUBROUTINE CART_FOCK(X,Y,I,J,TOGGLE)
! ACT THE FOCK OPERATOR (MINUS KINETIC) ON X AND ADD IT TO Y

   USE CONSTANTS
   USE STRUCTURE
   USE BASISSET
   USE INTEGRAL
   USE CARTESIANGRID

   IMPLICIT NONE
   INTEGER :: I,J
   DOUBLE PRECISION :: X(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   DOUBLE PRECISION :: Y(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   DOUBLE PRECISION,ALLOCATABLE :: Z(:,:,:)
   INTEGER :: IX,IY,IZ,JX,JY,JZ,KX,KY,KZ
   LOGICAL :: TOGGLE
   DOUBLE PRECISION :: SGN
   INTEGER :: IA
   DOUBLE PRECISION :: X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,R,C
   INTEGER :: INFO

   ALLOCATE(Z(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN))

   IF (TOGGLE) THEN
    SGN=1.0D0
   ELSE
    SGN=-1.0D0
   ENDIF

   DO IX=-MAXN,MAXN
write(*,*) ix
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
!write(*,*) ix,iy,iz

      ! (EPSILON_I + EPSILON_J) X(R1,R2)
      DO JX=-MAXN,MAXN
       DO JY=-MAXN,MAXN
        DO JZ=-MAXN,MAXN
         Y(JX,JY,JZ,IX,IY,IZ)=Y(JX,JY,JZ,IX,IY,IZ)+(EPSILON(I,0)+EPSILON(J,0))*X(JX,JY,JZ,IX,IY,IZ)
! debug
Y(JX,JY,JZ,IX,IY,IZ)=0.0d0
        ENDDO
       ENDDO
      ENDDO

      DO JX=-MAXN,MAXN
       DO JY=-MAXN,MAXN
        DO JZ=-MAXN,MAXN
         X2=DFLOAT(JX)*DELTAH
         Y2=DFLOAT(JY)*DELTAH
         Z2=DFLOAT(JZ)*DELTAH
         ! NUCLEAR ATTRACTION (FOR GRID POINT RIGHT ON TOP OF NUCLEUS, IT IS ZERO AS IT 
         ! CANCELS WITH KINETIC CONTRIBUTION)
         DO IA=1,NATOM
          X1=ATOMX(IA)
          Y1=ATOMY(IA)
          Z1=ATOMZ(IA)
          R=(X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2
          IF (R > 0.0D0) R=1.0D0/SQRT(R)
          Y(JX,JY,JZ,IX,IY,IZ)=Y(JX,JY,JZ,IX,IY,IZ)+DFLOAT(IATOM(IA))*R*X(JX,JY,JZ,IX,IY,IZ)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      ! COULOMB
      C=0.0D0
      DO JX=-MAXN,MAXN
       DO JY=-MAXN,MAXN
        DO JZ=-MAXN,MAXN
         DO IA=1,IOCC
          C=C+2.0D0*ORB(JX,JY,JZ,IA)**2*DELTAH**3
         ENDDO
        ENDDO
       ENDDO
      ENDDO
!     IF (DABS(C-DFLOAT(2*IOCC)) > 0.1D0*DFLOAT(2*IOCC)) THEN
!      CALL WARNING("ERROR IN COULOMB CHARGE")
!      WRITE(6,'(2F20.10)') C,DFLOAT(2*IOCC)
!     ENDIF
      C=DFLOAT(2*IOCC)

      Z=0.0D0
      DO JX=-MAXN,MAXN
       DO JY=-MAXN,MAXN
        DO JZ=-MAXN,MAXN
         DO IA=1,IOCC
          Z(JX,JY,JZ)=Z(JX,JY,JZ)-8.0D0*PI*ORB(JX,JY,JZ,IA)**2
         ENDDO
         X1=DFLOAT(JX)*DELTAH
         Y1=DFLOAT(JY)*DELTAH
         Z1=DFLOAT(JZ)*DELTAH
         R=X1**2+Y1**2+Z1**2
         IF (R > 0.0D0) R=C/DSQRT(R)
         IF (JX == -MAXN+1) Z(JX,JY,JZ)=Z(JX,JY,JZ) +1.0D0/12.0D0/DELTAH**2*R
         IF (JX ==  MAXN-1) Z(JX,JY,JZ)=Z(JX,JY,JZ) +1.0D0/12.0D0/DELTAH**2*R
         IF (JX == -MAXN  ) Z(JX,JY,JZ)=Z(JX,JY,JZ)-11.0D0/12.0D0/DELTAH**2*R
         IF (JX ==  MAXN  ) Z(JX,JY,JZ)=Z(JX,JY,JZ)-11.0D0/12.0D0/DELTAH**2*R
         IF (JY == -MAXN+1) Z(JX,JY,JZ)=Z(JX,JY,JZ) +1.0D0/12.0D0/DELTAH**2*R
         IF (JY ==  MAXN-1) Z(JX,JY,JZ)=Z(JX,JY,JZ) +1.0D0/12.0D0/DELTAH**2*R
         IF (JY == -MAXN  ) Z(JX,JY,JZ)=Z(JX,JY,JZ)-11.0D0/12.0D0/DELTAH**2*R
         IF (JY ==  MAXN  ) Z(JX,JY,JZ)=Z(JX,JY,JZ)-11.0D0/12.0D0/DELTAH**2*R
         IF (JZ == -MAXN+1) Z(JX,JY,JZ)=Z(JX,JY,JZ) +1.0D0/12.0D0/DELTAH**2*R
         IF (JZ ==  MAXN-1) Z(JX,JY,JZ)=Z(JX,JY,JZ) +1.0D0/12.0D0/DELTAH**2*R
         IF (JZ == -MAXN  ) Z(JX,JY,JZ)=Z(JX,JY,JZ)-11.0D0/12.0D0/DELTAH**2*R
         IF (JZ ==  MAXN  ) Z(JX,JY,JZ)=Z(JX,JY,JZ)-11.0D0/12.0D0/DELTAH**2*R
        ENDDO
       ENDDO
      ENDDO
!if ((ix==-maxn).and.(iy==-maxn).and.(iz==-maxn)) then
! do jx=-maxn,maxn
!  do jy=-maxn,maxn
!   write(60,*) jx,jy,z(jx,jy,0)
!  enddo
! enddo
!endif
      CALL DGETRS('N',MAXN3,1,DEL2,MAXN3,IPIV,Z,MAXN3,INFO)
!if ((ix==-maxn).and.(iy==-maxn).and.(iz==-maxn)) then
! write(*,*) "info=",info
! do jx=-maxn,maxn
!  do jy=-maxn,maxn
!   write(61,*) jx,jy,z(jx,jy,0)
!  enddo
! enddo
! z=0.0d0
! do jx=-maxn,maxn
!  do jy=-maxn,maxn
!   x1=dfloat(jx)*deltah
!   y1=dfloat(jy)*deltah
!   z1=0.0d0
!   do kx=-maxn,maxn
!    do ky=-maxn,maxn
!     do kz=-maxn,maxn
!      x2=dfloat(kx)*deltah
!      y2=dfloat(ky)*deltah
!      z2=dfloat(kz)*deltah
!      r=(x1-x2)**2+(y1-y2)**2+(z1-z2)**2
!      if (r > 0.0d0) r=1.0d0/dsqrt(r)
!      z(jx,jy,0)=z(jx,jy,0)+2.0d0*r*orb(kx,ky,kz,1)**2*deltah**3
!     enddo
!    enddo
!   enddo
!   write(62,*) jx,jy,z(jx,jy,0)
!  enddo
! enddo 
! stop
!endif
      DO JX=-MAXN,MAXN
       DO JY=-MAXN,MAXN
        DO JZ=-MAXN,MAXN
         Y(JX,JY,JZ,IX,IY,IZ)=Y(JX,JY,JZ,IX,IY,IZ)-Z(JX,JY,JZ)*X(JX,JY,JZ,IX,IY,IZ)
        ENDDO
       ENDDO
      ENDDO

      ! EXCHANGE
      DO IA=1,IOCC

       C=0.0D0
       DO JX=-MAXN,MAXN
        DO JY=-MAXN,MAXN
         DO JZ=-MAXN,MAXN
          C=C+ORB(JX,JY,JZ,IA)*X(JX,JY,JZ,IX,IY,IZ)*DELTAH**3
         ENDDO
        ENDDO
       ENDDO
!      IF (DABS(C) > 0.1D0) THEN
!       CALL WARNING("ERROR IN EXCHANGE CHARGE")
!       WRITE(6,'(2F20.10)') C,0.0D0
!      ENDIF
! debug
c=1.0d0
!      C=0.0D0

       Z=0.0D0
       DO JX=-MAXN,MAXN
        DO JY=-MAXN,MAXN
         DO JZ=-MAXN,MAXN
          Z(JX,JY,JZ)=Z(JX,JY,JZ)-4.0D0*PI*ORB(JX,JY,JZ,IA)*X(JX,JY,JZ,IX,IY,IZ)
          X1=DFLOAT(JX)*DELTAH
          Y1=DFLOAT(JY)*DELTAH
          Z1=DFLOAT(JZ)*DELTAH
          R=X1**2+Y1**2+Z1**2
          IF (R > 0.0D0) R=C/DSQRT(R)
          IF (JX == -MAXN+1) Z(JX,JY,JZ)=Z(JX,JY,JZ) +1.0D0/12.0D0/DELTAH**2*R
          IF (JX ==  MAXN-1) Z(JX,JY,JZ)=Z(JX,JY,JZ) +1.0D0/12.0D0/DELTAH**2*R
          IF (JX == -MAXN  ) Z(JX,JY,JZ)=Z(JX,JY,JZ)-11.0D0/12.0D0/DELTAH**2*R
          IF (JX ==  MAXN  ) Z(JX,JY,JZ)=Z(JX,JY,JZ)-11.0D0/12.0D0/DELTAH**2*R
          IF (JY == -MAXN+1) Z(JX,JY,JZ)=Z(JX,JY,JZ) +1.0D0/12.0D0/DELTAH**2*R
          IF (JY ==  MAXN-1) Z(JX,JY,JZ)=Z(JX,JY,JZ) +1.0D0/12.0D0/DELTAH**2*R
          IF (JY == -MAXN  ) Z(JX,JY,JZ)=Z(JX,JY,JZ)-11.0D0/12.0D0/DELTAH**2*R
          IF (JY ==  MAXN  ) Z(JX,JY,JZ)=Z(JX,JY,JZ)-11.0D0/12.0D0/DELTAH**2*R
          IF (JZ == -MAXN+1) Z(JX,JY,JZ)=Z(JX,JY,JZ) +1.0D0/12.0D0/DELTAH**2*R
          IF (JZ ==  MAXN-1) Z(JX,JY,JZ)=Z(JX,JY,JZ) +1.0D0/12.0D0/DELTAH**2*R
          IF (JZ == -MAXN  ) Z(JX,JY,JZ)=Z(JX,JY,JZ)-11.0D0/12.0D0/DELTAH**2*R
          IF (JZ ==  MAXN  ) Z(JX,JY,JZ)=Z(JX,JY,JZ)-11.0D0/12.0D0/DELTAH**2*R
         ENDDO
        ENDDO
       ENDDO
       CALL DGETRS('N',MAXN3,1,DEL2,MAXN3,IPIV,Z,MAXN3,INFO)
       DO JX=-MAXN,MAXN
        DO JY=-MAXN,MAXN
         DO JZ=-MAXN,MAXN
          Y(JX,JY,JZ,IX,IY,IZ)=Y(JX,JY,JZ,IX,IY,IZ)+Z(JX,JY,JZ)*X(JX,JY,JZ,IX,IY,IZ)
         ENDDO
        ENDDO
       ENDDO

      ENDDO

     ENDDO
    ENDDO
   ENDDO

   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      DO JX=-MAXN,MAXN
       DO JY=-MAXN,MAXN
        DO JZ=-MAXN,MAXN
! debug
!        Y(IX,IY,IZ,JX,JY,JZ)=SGN*Y(JX,JY,JZ,IX,IY,IZ)
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO

   DEALLOCATE(Z)

   RETURN
END SUBROUTINE



SUBROUTINE CART_DEL2
! CONSTRUCT DEL SQUARE MATRIX

   USE STRUCTURE
   USE BASISSET
   USE INTEGRAL
   USE CARTESIANGRID

   IMPLICIT NONE
   INTEGER :: IX,IY,IZ,JX,JY,JZ

   DEL2=0.0D0
   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      IF ((IX >= -MAXN+2).AND.(IX <= MAXN-2)) THEN
       DEL2(IX,IY,IZ,IX-2,IY,IZ)=DEL2(IX,IY,IZ,IX-2,IY,IZ) -1.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX-1,IY,IZ)=DEL2(IX,IY,IZ,IX-1,IY,IZ)+16.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,  IY,IZ)=DEL2(IX,IY,IZ,IX,  IY,IZ)-30.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX+1,IY,IZ)=DEL2(IX,IY,IZ,IX+1,IY,IZ)+16.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX+2,IY,IZ)=DEL2(IX,IY,IZ,IX+2,IY,IZ) -1.0D0/12.0D0/DELTAH**2
      ELSE IF (IX == -MAXN+1) THEN
       DEL2(IX,IY,IZ,IX-1,IY,IZ)=DEL2(IX,IY,IZ,IX-1,IY,IZ)+16.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,  IY,IZ)=DEL2(IX,IY,IZ,IX,  IY,IZ)-30.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX+1,IY,IZ)=DEL2(IX,IY,IZ,IX+1,IY,IZ)+16.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX+2,IY,IZ)=DEL2(IX,IY,IZ,IX+2,IY,IZ) -1.0D0/12.0D0/DELTAH**2
      ELSE IF (IX == MAXN-1) THEN
       DEL2(IX,IY,IZ,IX-2,IY,IZ)=DEL2(IX,IY,IZ,IX-2,IY,IZ) -1.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX-1,IY,IZ)=DEL2(IX,IY,IZ,IX-1,IY,IZ)+16.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,  IY,IZ)=DEL2(IX,IY,IZ,IX,  IY,IZ)-30.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX+1,IY,IZ)=DEL2(IX,IY,IZ,IX+1,IY,IZ)+16.0D0/12.0D0/DELTAH**2
      ELSE IF (IX == -MAXN) THEN
       DEL2(IX,IY,IZ,IX,  IY,IZ)=DEL2(IX,IY,IZ,IX,  IY,IZ)-20.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX+1,IY,IZ)=DEL2(IX,IY,IZ,IX+1,IY,IZ) +6.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX+2,IY,IZ)=DEL2(IX,IY,IZ,IX+2,IY,IZ) +4.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX+3,IY,IZ)=DEL2(IX,IY,IZ,IX+3,IY,IZ) -1.0D0/12.0D0/DELTAH**2
      ELSE IF (IX == MAXN) THEN
       DEL2(IX,IY,IZ,IX-3,IY,IZ)=DEL2(IX,IY,IZ,IX-3,IY,IZ) -1.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX-2,IY,IZ)=DEL2(IX,IY,IZ,IX-2,IY,IZ) +4.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX-1,IY,IZ)=DEL2(IX,IY,IZ,IX-1,IY,IZ) +6.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,  IY,IZ)=DEL2(IX,IY,IZ,IX,  IY,IZ)-20.0D0/12.0D0/DELTAH**2
      ELSE
       CALL PABORT('A BUG')
      ENDIF
      IF ((IY >= -MAXN+2).AND.(IY <= MAXN-2)) THEN
       DEL2(IX,IY,IZ,IX,IY-2,IZ)=DEL2(IX,IY,IZ,IX,IY-2,IZ) -1.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY-1,IZ)=DEL2(IX,IY,IZ,IX,IY-1,IZ)+16.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,  IZ)=DEL2(IX,IY,IZ,IX,IY,  IZ)-30.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY+1,IZ)=DEL2(IX,IY,IZ,IX,IY+1,IZ)+16.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY+2,IZ)=DEL2(IX,IY,IZ,IX,IY+2,IZ) -1.0D0/12.0D0/DELTAH**2
      ELSE IF (IY == -MAXN+1) THEN
       DEL2(IX,IY,IZ,IX,IY-1,IZ)=DEL2(IX,IY,IZ,IX,IY-1,IZ)+16.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,  IZ)=DEL2(IX,IY,IZ,IX,IY,  IZ)-30.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY+1,IZ)=DEL2(IX,IY,IZ,IX,IY+1,IZ)+16.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY+2,IZ)=DEL2(IX,IY,IZ,IX,IY+2,IZ) -1.0D0/12.0D0/DELTAH**2
      ELSE IF (IY == MAXN-1) THEN
       DEL2(IX,IY,IZ,IX,IY-2,IZ)=DEL2(IX,IY,IZ,IX,IY-2,IZ) -1.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY-1,IZ)=DEL2(IX,IY,IZ,IX,IY-1,IZ)+16.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,  IZ)=DEL2(IX,IY,IZ,IX,IY,  IZ)-30.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY+1,IZ)=DEL2(IX,IY,IZ,IX,IY+1,IZ)+16.0D0/12.0D0/DELTAH**2
      ELSE IF (IY == -MAXN) THEN
       DEL2(IX,IY,IZ,IX,IY,  IZ)=DEL2(IX,IY,IZ,IX,IY,  IZ)-20.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY+1,IZ)=DEL2(IX,IY,IZ,IX,IY+1,IZ) +6.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY+2,IZ)=DEL2(IX,IY,IZ,IX,IY+2,IZ) +4.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY+3,IZ)=DEL2(IX,IY,IZ,IX,IY+3,IZ) -1.0D0/12.0D0/DELTAH**2
      ELSE IF (IY == MAXN) THEN
       DEL2(IX,IY,IZ,IX,IY-3,IZ)=DEL2(IX,IY,IZ,IX,IY-3,IZ) -1.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY-2,IZ)=DEL2(IX,IY,IZ,IX,IY-2,IZ) +4.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY-1,IZ)=DEL2(IX,IY,IZ,IX,IY-1,IZ) +6.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,  IZ)=DEL2(IX,IY,IZ,IX,IY,  IZ)-20.0D0/12.0D0/DELTAH**2
      ELSE
       CALL PABORT('A BUG')
      ENDIF
      IF ((IZ >= -MAXN+2).AND.(IZ <= MAXN-2)) THEN
       DEL2(IX,IY,IZ,IX,IY,IZ-2)=DEL2(IX,IY,IZ,IX,IY,IZ-2) -1.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,IZ-1)=DEL2(IX,IY,IZ,IX,IY,IZ-1)+16.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,IZ  )=DEL2(IX,IY,IZ,IX,IY,IZ  )-30.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,IZ+1)=DEL2(IX,IY,IZ,IX,IY,IZ+1)+16.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,IZ+2)=DEL2(IX,IY,IZ,IX,IY,IZ+2) -1.0D0/12.0D0/DELTAH**2
      ELSE IF (IZ == -MAXN+1) THEN
       DEL2(IX,IY,IZ,IX,IY,IZ-1)=DEL2(IX,IY,IZ,IX,IY,IZ-1)+16.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,IZ  )=DEL2(IX,IY,IZ,IX,IY,IZ  )-30.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,IZ+1)=DEL2(IX,IY,IZ,IX,IY,IZ+1)+16.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,IZ+2)=DEL2(IX,IY,IZ,IX,IY,IZ+2) -1.0D0/12.0D0/DELTAH**2
      ELSE IF (IZ == MAXN-1) THEN
       DEL2(IX,IY,IZ,IX,IY,IZ-2)=DEL2(IX,IY,IZ,IX,IY,IZ-2) -1.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,IZ-1)=DEL2(IX,IY,IZ,IX,IY,IZ-1)+16.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,IZ  )=DEL2(IX,IY,IZ,IX,IY,IZ  )-30.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,IZ+1)=DEL2(IX,IY,IZ,IX,IY,IZ+1)+16.0D0/12.0D0/DELTAH**2
      ELSE IF (IZ == -MAXN) THEN
       DEL2(IX,IY,IZ,IX,IY,IZ  )=DEL2(IX,IY,IZ,IX,IY,IZ  )-20.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,IZ+1)=DEL2(IX,IY,IZ,IX,IY,IZ+1) +6.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,IZ+2)=DEL2(IX,IY,IZ,IX,IY,IZ+2) +4.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,IZ+3)=DEL2(IX,IY,IZ,IX,IY,IZ+3) -1.0D0/12.0D0/DELTAH**2
      ELSE IF (IZ == MAXN) THEN
       DEL2(IX,IY,IZ,IX,IY,IZ-3)=DEL2(IX,IY,IZ,IX,IY,IZ-3) -1.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,IZ-2)=DEL2(IX,IY,IZ,IX,IY,IZ-2) +4.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,IZ-1)=DEL2(IX,IY,IZ,IX,IY,IZ-1) +6.0D0/12.0D0/DELTAH**2
       DEL2(IX,IY,IZ,IX,IY,IZ  )=DEL2(IX,IY,IZ,IX,IY,IZ  )-20.0D0/12.0D0/DELTAH**2
      ELSE
       CALL PABORT('A BUG')
      ENDIF
     ENDDO
    ENDDO
   ENDDO

   RETURN
END SUBROUTINE



SUBROUTINE CART_KINETIC(X,Y,I,J,TOGGLE)
! ACT THE KINETIC-ENERGY OPERATOR (TIMES -1) ON X AND ADD IT TO Y

   USE CONSTANTS
   USE STRUCTURE
   USE BASISSET
   USE INTEGRAL
   USE CARTESIANGRID

   IMPLICIT NONE
   INTEGER :: I,J
   DOUBLE PRECISION :: X(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   DOUBLE PRECISION :: Y(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN)
   INTEGER :: IX,IY,IZ,JX,JY,JZ
   LOGICAL :: TOGGLE
   DOUBLE PRECISION :: X1,Y1,Z1,X2,Y2,Z2,R
   INTEGER :: IA
   LOGICAL :: LCYCLE
   DOUBLE PRECISION :: SGN

   IF (TOGGLE) THEN
    SGN=1.0D0
   ELSE
    SGN=-1.0D0
   ENDIF

   DO JX=-MAXN,MAXN
    DO JY=-MAXN,MAXN
     DO JZ=-MAXN,MAXN
      DO IX=-MAXN,MAXN
       DO IY=-MAXN,MAXN
        DO IZ=-MAXN,MAXN
         X1=DFLOAT(IX)*DELTAH
         Y1=DFLOAT(IY)*DELTAH
         Z1=DFLOAT(IZ)*DELTAH
         LCYCLE=.FALSE.
         DO IA=1,NATOM
          X2=ATOMX(IA)
          Y2=ATOMY(IA)
          Z2=ATOMZ(IA)
          R=(X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2
          IF (R==0.0D0) LCYCLE=.TRUE.
         ENDDO
         IF (LCYCLE) CYCLE
         IF ((IX >= -MAXN+2).AND.(IX <= MAXN-2)) THEN
          Y(IX,IY,IZ,JX,JY,JZ)=Y(IX,IY,IZ,JX,JY,JZ)+(0.5D0/12.0D0/DELTAH**2)*( &
                               -1.0D0*X(IX-2,IY,IZ,JX,JY,JZ) &
                              +16.0D0*X(IX-1,IY,IZ,JX,JY,JZ) &
                              -30.0D0*X(IX,  IY,IZ,JX,JY,JZ) &
                              +16.0D0*X(IX+1,IY,IZ,JX,JY,JZ) &
                               -1.0D0*X(IX+2,IY,IZ,JX,JY,JZ))
         ELSE IF (IX == -MAXN+1) THEN
          Y(IX,IY,IZ,JX,JY,JZ)=Y(IX,IY,IZ,JX,JY,JZ)+(0.5D0/12.0D0/DELTAH**2)*( &
                              +11.0D0*X(IX-1,IY,IZ,JX,JY,JZ) &
                              -20.0D0*X(IX,  IY,IZ,JX,JY,JZ) &
                               +6.0D0*X(IX+1,IY,IZ,JX,JY,JZ) &
                               +4.0D0*X(IX+2,IY,IZ,JX,JY,JZ) &
                               -1.0D0*X(IX+3,IY,IZ,JX,JY,JZ))
         ELSE IF (IX == MAXN-1) THEN
          Y(IX,IY,IZ,JX,JY,JZ)=Y(IX,IY,IZ,JX,JY,JZ)+(0.5D0/12.0D0/DELTAH**2)*( &
                               -1.0D0*X(IX-3,IY,IZ,JX,JY,JZ) &
                               +4.0D0*X(IX-2,IY,IZ,JX,JY,JZ) &
                               +6.0D0*X(IX-1,IY,IZ,JX,JY,JZ) &
                              -20.0D0*X(IX,  IY,IZ,JX,JY,JZ) &
                              +11.0D0*X(IX+1,IY,IZ,JX,JY,JZ))
         ELSE IF (IX == -MAXN) THEN
          Y(IX,IY,IZ,JX,JY,JZ)=Y(IX,IY,IZ,JX,JY,JZ)+(0.5D0/12.0D0/DELTAH**2)*( &
                              +35.0D0*X(IX,  IY,IZ,JX,JY,JZ) &
                             -104.0D0*X(IX+1,IY,IZ,JX,JY,JZ) &
                             +114.0D0*X(IX+2,IY,IZ,JX,JY,JZ) &
                              -56.0D0*X(IX+3,IY,IZ,JX,JY,JZ) &
                              +11.0D0*X(IX+4,IY,IZ,JX,JY,JZ))
         ELSE IF (IX == MAXN) THEN
          Y(IX,IY,IZ,JX,JY,JZ)=Y(IX,IY,IZ,JX,JY,JZ)+(0.5D0/12.0D0/DELTAH**2)*( &
                              +11.0D0*X(IX-4,IY,IZ,JX,JY,JZ) &
                              -56.0D0*X(IX-3,IY,IZ,JX,JY,JZ) &
                             +114.0D0*X(IX-2,IY,IZ,JX,JY,JZ) &
                             -104.0D0*X(IX-1,IY,IZ,JX,JY,JZ) &
                              +35.0D0*X(IX,  IY,IZ,JX,JY,JZ))
         ENDIF
         IF ((IY >= -MAXN+2).AND.(IY <= MAXN-2)) THEN
          Y(IX,IY,IZ,JX,JY,JZ)=Y(IX,IY,IZ,JX,JY,JZ)+(0.5D0/12.0D0/DELTAH**2)*( &
                               -1.0D0*X(IX,IY-2,IZ,JX,JY,JZ) &
                              +16.0D0*X(IX,IY-1,IZ,JX,JY,JZ) &
                              -30.0D0*X(IX,IY,  IZ,JX,JY,JZ) &
                              +16.0D0*X(IX,IY+1,IZ,JX,JY,JZ) &
                               -1.0D0*X(IX,IY+2,IZ,JX,JY,JZ))
         ELSE IF (IY == -MAXN+1) THEN
          Y(IX,IY,IZ,JX,JY,JZ)=Y(IX,IY,IZ,JX,JY,JZ)+(0.5D0/12.0D0/DELTAH**2)*( &
                              +11.0D0*X(IX,IY-1,IZ,JX,JY,JZ) &
                              -20.0D0*X(IX,IY,  IZ,JX,JY,JZ) &
                               +6.0D0*X(IX,IY+1,IZ,JX,JY,JZ) &
                               +4.0D0*X(IX,IY+2,IZ,JX,JY,JZ) &
                               -1.0D0*X(IX,IY+3,IZ,JX,JY,JZ))
         ELSE IF (IY == MAXN-1) THEN
          Y(IX,IY,IZ,JX,JY,JZ)=Y(IX,IY,IZ,JX,JY,JZ)+(0.5D0/12.0D0/DELTAH**2)*( &
                               -1.0D0*X(IX,IY-3,IZ,JX,JY,JZ) &
                               +4.0D0*X(IX,IY-2,IZ,JX,JY,JZ) &
                               +6.0D0*X(IX,IY-1,IZ,JX,JY,JZ) &
                              -20.0D0*X(IX,IY,  IZ,JX,JY,JZ) &
                              +11.0D0*X(IX,IY+1,IZ,JX,JY,JZ))
         ELSE IF (IY == -MAXN) THEN
          Y(IX,IY,IZ,JX,JY,JZ)=Y(IX,IY,IZ,JX,JY,JZ)+(0.5D0/12.0D0/DELTAH**2)*( &
                              +35.0D0*X(IX,IY,  IZ,JX,JY,JZ) &
                             -104.0D0*X(IX,IY+1,IZ,JX,JY,JZ) &
                             +114.0D0*X(IX,IY+2,IZ,JX,JY,JZ) &
                              -56.0D0*X(IX,IY+3,IZ,JX,JY,JZ) &
                              +11.0D0*X(IX,IY+4,IZ,JX,JY,JZ))
         ELSE IF (IY == MAXN) THEN
          Y(IX,IY,IZ,JX,JY,JZ)=Y(IX,IY,IZ,JX,JY,JZ)+(0.5D0/12.0D0/DELTAH**2)*( &
                              +11.0D0*X(IX,IY-4,IZ,JX,JY,JZ) &
                              -56.0D0*X(IX,IY-3,IZ,JX,JY,JZ) &
                             +114.0D0*X(IX,IY-2,IZ,JX,JY,JZ) &
                             -104.0D0*X(IX,IY-1,IZ,JX,JY,JZ) &
                              +35.0D0*X(IX,IY,  IZ,JX,JY,JZ))
         ENDIF
         IF ((IZ >= -MAXN+2).AND.(IZ <= MAXN-2)) THEN
          Y(IX,IY,IZ,JX,JY,JZ)=Y(IX,IY,IZ,JX,JY,JZ)+(0.5D0/12.0D0/DELTAH**2)*( &
                               -1.0D0*X(IX,IY,IZ-2,JX,JY,JZ) &
                              +16.0D0*X(IX,IY,IZ-1,JX,JY,JZ) &
                              -30.0D0*X(IX,IY,IZ,  JX,JY,JZ) &
                              +16.0D0*X(IX,IY,IZ+1,JX,JY,JZ) &
                               -1.0D0*X(IX,IY,IZ+2,JX,JY,JZ))
         ELSE IF (IZ == -MAXN+1) THEN
          Y(IX,IY,IZ,JX,JY,JZ)=Y(IX,IY,IZ,JX,JY,JZ)+(0.5D0/12.0D0/DELTAH**2)*( &
                              +11.0D0*X(IX,IY,IZ-1,JX,JY,JZ) &
                              -20.0D0*X(IX,IY,IZ,  JX,JY,JZ) &
                               +6.0D0*X(IX,IY,IZ+1,JX,JY,JZ) &
                               +4.0D0*X(IX,IY,IZ+2,JX,JY,JZ) &
                               -1.0D0*X(IX,IY,IZ+3,JX,JY,JZ))
         ELSE IF (IZ == MAXN-1) THEN
          Y(IX,IY,IZ,JX,JY,JZ)=Y(IX,IY,IZ,JX,JY,JZ)+(0.5D0/12.0D0/DELTAH**2)*( &
                               -1.0D0*X(IX,IY,IZ-3,JX,JY,JZ) &
                               +4.0D0*X(IX,IY,IZ-2,JX,JY,JZ) &
                               +6.0D0*X(IX,IY,IZ-1,JX,JY,JZ) &
                              -20.0D0*X(IX,IY,IZ,  JX,JY,JZ) &
                              +11.0D0*X(IX,IY,IZ+1,JX,JY,JZ))
         ELSE IF (IZ == -MAXN) THEN
          Y(IX,IY,IZ,JX,JY,JZ)=Y(IX,IY,IZ,JX,JY,JZ)+(0.5D0/12.0D0/DELTAH**2)*( &
                              +35.0D0*X(IX,IY,IZ,  JX,JY,JZ) &
                             -104.0D0*X(IX,IY,IZ+1,JX,JY,JZ) &
                             +114.0D0*X(IX,IY,IZ+2,JX,JY,JZ) &
                              -56.0D0*X(IX,IY,IZ+3,JX,JY,JZ) &
                              +11.0D0*X(IX,IY,IZ+4,JX,JY,JZ))
         ELSE IF (IZ == MAXN) THEN
          Y(IX,IY,IZ,JX,JY,JZ)=Y(IX,IY,IZ,JX,JY,JZ)+(0.5D0/12.0D0/DELTAH**2)*( &
                              +11.0D0*X(IX,IY,IZ-4,JX,JY,JZ) &
                              -56.0D0*X(IX,IY,IZ-3,JX,JY,JZ) &
                             +114.0D0*X(IX,IY,IZ-2,JX,JY,JZ) &
                             -104.0D0*X(IX,IY,IZ-1,JX,JY,JZ) &
                              +35.0D0*X(IX,IY,IZ,  JX,JY,JZ))
         ENDIF
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO

   DO JX=-MAXN,MAXN
    DO JY=-MAXN,MAXN
     DO JZ=-MAXN,MAXN
      DO IX=-MAXN,MAXN
       DO IY=-MAXN,MAXN
        DO IZ=-MAXN,MAXN
! debug
!        Y(JX,JY,JZ,IX,IY,IZ)=SGN*Y(IX,IY,IZ,JX,JY,JZ)
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO

   RETURN
END SUBROUTINE



SUBROUTINE CART_RESIDUAL(I,J,X,Y,D)
! CALCULATE THE RESIDUAL

   USE STRUCTURE
   USE BASISSET
   USE INTEGRAL
   USE CARTESIANGRID

   IMPLICIT NONE
   INTEGER :: I,J ! OCCUPIED ORBITAL PAIR
   DOUBLE PRECISION :: X(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN) ! RHS
   DOUBLE PRECISION :: Y(-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN,-MAXN:MAXN) ! LHS
   DOUBLE PRECISION :: D ! RESIDUAL
   INTEGER :: IX,IY,IZ,JX,JY,JZ
   DOUBLE PRECISION :: X1,Y1,Z1,X2,Y2,Z2,R

   D=0.0D0
   DO IX=-MAXN,MAXN
    DO IY=-MAXN,MAXN
     DO IZ=-MAXN,MAXN
      DO JX=-MAXN,MAXN
       DO JY=-MAXN,MAXN
        DO JZ=-MAXN,MAXN
         D=D+(X(JX,JY,JZ,IX,IY,IZ)-Y(JX,JY,JZ,IX,IY,IZ))**2
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO

   RETURN
END SUBROUTINE
