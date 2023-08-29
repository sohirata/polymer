SUBROUTINE FMT_PRELIMINARY
! PRELIMINARY FOR FMT EVALUATIONS.

   USE CONTROL
   USE FMT
   USE BASISSET
   USE CONSTANTS

   IMPLICIT NONE
   INTEGER :: I,J,K
   DOUBLE PRECISION :: A,B,C

   IF (BASIS_LMAX*4+1 > 21) CALL PABORT('EXTEND FMT')
   IGAMMA=0.0D0
   DO I=0,1700
    A=DFLOAT(I)/20.0D0
    B=DEXP(-A)
    C=0.0D0
    K=50+I/3.
    DO J=K,0,-1
     C=(2.0D0*A*C+B)/DFLOAT(2*J+1)
     IF (J <= 27) IGAMMA(I,J)=C
    ENDDO
   ENDDO

   RETURN
END SUBROUTINE



SUBROUTINE FORM_C2S
! CONSTRUCT CARTERSIAN TO SPHERICAL CONVERSION MATRIX.

   USE MPI
   USE CONTROL
   USE STRUCTURE
   USE BASISSET

   IMPLICIT NONE
   INTEGER :: I,J,K
   INTEGER,ALLOCATABLE :: P2C(:)

   IF (.NOT.LOPTN(69)) THEN
    C2S=0.0D0
    DO I=1,NCGS
     C2S(I,I)=1.0D0
    ENDDO
   ELSE
    ALLOCATE(P2C(NPGS))
    DO I=1,NPGS
     K=0
     DO J=1,NCGS
      IF (CC(J,I) /= 0.0D0) THEN
       K=K+1
       P2C(I)=J
      ENDIF
     ENDDO
     IF (K == 0) THEN
      P2C(I)=1
     ELSE IF (K > 1) THEN
      CALL PABORT('GENERAL CONTRACTION NOT YET IMPLEMENTED')
     ENDIF
    ENDDO
    C2S=0.0D0
    DO I=1,NCGS
     C2S(I,I)=1.0D0
    ENDDO
    ! D ORBITALS: 6 -> 5
    DO I=1,NPGS
     IF ((PAX(I) == 2).AND.(PAY(I) == 0).AND.(PAZ(I) == 0)) THEN
      C2S(P2C(I),P2C(I))=0.0D0
     ELSE IF ((PAX(I) == 0).AND.(PAY(I) == 2).AND.(PAZ(I) == 0)) THEN
      C2S(P2C(I),P2C(I))=-DSQRT(3.0D0)/2.0D0
      C2S(P2C(I),P2C(I-1))=DSQRT(3.0D0)/2.0D0
     ELSE IF ((PAX(I) == 0).AND.(PAY(I) == 0).AND.(PAZ(I) == 2)) THEN
      C2S(P2C(I),P2C(I-2))=-1.0D0/2.0D0
      C2S(P2C(I),P2C(I-1))=-1.0D0/2.0D0
     ENDIF
    ENDDO
    ! F ORBITALS: 10 -> 7
    DO I=1,NPGS
     IF ((PAX(I) == 3).AND.(PAY(I) == 0).AND.(PAZ(I) == 0)) THEN
      C2S(P2C(I),P2C(I+8))=-3.0D0*DSQRT(5.0D0)/10.0D0
      C2S(P2C(I),P2C(I+4))=-3.0D0*DSQRT(5.0D0)/10.0D0
      C2S(P2C(I),P2C(I+2))=1.0D0
      C2S(P2C(I),P2C(I))=0.0D0
     ELSE IF ((PAX(I) == 0).AND.(PAY(I) == 3).AND.(PAZ(I) == 0)) THEN
      C2S(P2C(I),P2C(I-1))=-DSQRT(6.0D0)/4.0D0
      C2S(P2C(I),P2C(I+5))=-DSQRT(30.0D0)/20.0D0
      C2S(P2C(I),P2C(I+4))=DSQRT(30.0D0)/5.0D0
      C2S(P2C(I),P2C(I))=0.0D0
     ELSE IF ((PAX(I) == 0).AND.(PAY(I) == 0).AND.(PAZ(I) == 3)) THEN
      C2S(P2C(I),P2C(I+1))=-DSQRT(30.0D0)/20.0D0
      C2S(P2C(I),P2C(I-1))=-DSQRT(6.0D0)/4.0D0
      C2S(P2C(I),P2C(I+5))=DSQRT(30.0D0)/5.0D0
      C2S(P2C(I),P2C(I))=0.0D0
     ELSE IF ((PAX(I) == 2).AND.(PAY(I) == 1).AND.(PAZ(I) == 0)) THEN
      C2S(P2C(I),P2C(I+5))=DSQRT(3.0D0)/2.0D0
      C2S(P2C(I),P2C(I+1))=-DSQRT(3.0D0)/2.0D0
      C2S(P2C(I),P2C(I))=0.0D0
     ELSE IF ((PAX(I) == 0).AND.(PAY(I) == 2).AND.(PAZ(I) == 1)) THEN
      C2S(P2C(I),P2C(I+5))=1.0D0
      C2S(P2C(I),P2C(I))=0.0D0
     ELSE IF ((PAX(I) == 1).AND.(PAY(I) == 0).AND.(PAZ(I) == 2)) THEN
      C2S(P2C(I),P2C(I-5))=DSQRT(10.0D0)/4.0D0
      C2S(P2C(I),P2C(I+1))=-3.0D0/4.0D0*DSQRT(2.0D0)
      C2S(P2C(I),P2C(I))=0.0D0
     ELSE IF ((PAX(I) == 1).AND.(PAY(I) == 2).AND.(PAZ(I) == 0)) THEN
      C2S(P2C(I),P2C(I-3))=3.0D0/4.0D0*DSQRT(2.0D0)
      C2S(P2C(I),P2C(I-5))=-DSQRT(10.0D0)/4.0D0
      C2S(P2C(I),P2C(I))=0.0D0
     ELSE IF ((PAX(I) == 0).AND.(PAY(I) == 1).AND.(PAZ(I) == 2)) THEN
      C2S(P2C(I),P2C(I))=0.0D0
     ELSE IF ((PAX(I) == 2).AND.(PAY(I) == 0).AND.(PAZ(I) == 1)) THEN
      C2S(P2C(I),P2C(I))=0.0D0
     ELSE IF ((PAX(I) == 1).AND.(PAY(I) == 1).AND.(PAZ(I) == 1)) THEN
      C2S(P2C(I),P2C(I))=0.0D0
     ENDIF
    ENDDO
    IF ((IOPTN(9) >= 2).AND.(MYID == 0)) THEN
     WRITE(6,'(A)') 'CONVERSION MATRIX FOR CARTESIAN TO SPHERICAL GAUSSIANS'
     CALL DUMP5(C2S,NCGS)
    ENDIF
    DEALLOCATE(P2C)
   ENDIF

   RETURN
END SUBROUTINE



SUBROUTINE FORM_C2H
! CONSTRUCT CARTERSIAN TO HELICAL/ZIGZAG CONVERSION MATRIX.

   USE MPI
   USE CONTROL
   USE STRUCTURE
   USE BASISSET
   USE INTEGRAL

   IMPLICIT NONE
   INTEGER :: I,J,K
   INTEGER,ALLOCATABLE :: P2C(:)
   DOUBLE PRECISION :: C,S
   DOUBLE PRECISION :: ANGLE
!double precision,allocatable :: tmp(:,:)
!allocate(tmp(ncgs,ncgs))

   IF (HELIX == 0.0D0) THEN
    C2H=0.0D0
    P2H=0.0D0
    C2H_THETA=0.0D0
    P2H_THETA=0.0D0
    DO J=-CEL2X-CEL1X,CEL2X+CEL1X
     DO I=1,NCGS
      C2H(I,I,J)=1.0D0
     ENDDO
     DO I=1,NPGS
      P2H(I,I,J)=1.0D0
     ENDDO
    ENDDO
   ELSE
    ALLOCATE(P2C(NPGS))
    DO I=1,NPGS
     K=0
     DO J=1,NCGS
      IF (CC(J,I) /= 0.0D0) THEN
       K=K+1
       P2C(I)=J
      ENDIF
     ENDDO
     IF (K == 0) THEN
      P2C(I)=-1
     ELSE IF (K > 1) THEN
      CALL PABORT('INTERNAL PROGRAM ERROR')
     ENDIF
    ENDDO
    C2H=0.0D0
    P2H=0.0D0
    C2H_THETA=0.0D0
    P2H_THETA=0.0D0
    DO J=-CEL2X-CEL1X,CEL2X+CEL1X
     DO I=1,NCGS
      C2H(I,I,J)=1.0D0
     ENDDO
     DO I=1,NPGS
      P2H(I,I,J)=1.0D0
     ENDDO
     ANGLE=HELIX*DFLOAT(J)
     C=DCOS(ANGLE)
     S=DSIN(ANGLE) 
     ! P ORBITALS
     DO I=1,NPGS
      IF ((PAX(I) == 1).AND.(PAY(I) == 0).AND.(PAZ(I) == 0)) THEN
       IF (P2C(I) == -1) CYCLE
       C2H(P2C(I  ),P2C(I  ),J)=1.0D0
       C2H(P2C(I+1),P2C(I+1),J)=C
       C2H(P2C(I+2),P2C(I+2),J)=C
       C2H(P2C(I+1),P2C(I+2),J)=S
       C2H(P2C(I+2),P2C(I+1),J)=-S
       P2H((I  ),(I  ),J)=1.0D0
       P2H((I+1),(I+1),J)=C
       P2H((I+2),(I+2),J)=C
       P2H((I+1),(I+2),J)=S
       P2H((I+2),(I+1),J)=-S
       C2H_THETA(P2C(I  ),P2C(I  ),J)=0.0D0
       C2H_THETA(P2C(I+1),P2C(I+1),J)=-S*DFLOAT(J)
       C2H_THETA(P2C(I+2),P2C(I+2),J)=-S*DFLOAT(J)
       C2H_THETA(P2C(I+1),P2C(I+2),J)=C*DFLOAT(J)
       C2H_THETA(P2C(I+2),P2C(I+1),J)=-C*DFLOAT(J)
       P2H_THETA((I  ),(I  ),J)=0.0D0
       P2H_THETA((I+1),(I+1),J)=-S*DFLOAT(J)
       P2H_THETA((I+2),(I+2),J)=-S*DFLOAT(J)
       P2H_THETA((I+1),(I+2),J)=C*DFLOAT(J)
       P2H_THETA((I+2),(I+1),J)=-C*DFLOAT(J)
      ENDIF
     ENDDO
     ! D ORBITALS
     DO I=1,NPGS
      IF ((PAX(I) == 2).AND.(PAY(I) == 0).AND.(PAZ(I) == 0)) THEN
       IF (P2C(I) == -1) CYCLE
       C2H(P2C(I  ),P2C(I  ),J)=1.0D0
       C2H(P2C(I+1),P2C(I+1),J)=C*C
       C2H(P2C(I+1),P2C(I+2),J)=S*S
       C2H(P2C(I+1),P2C(I+4),J)=2.0D0*C*S/DSQRT(3.0D0)
       C2H(P2C(I+2),P2C(I+1),J)=S*S
       C2H(P2C(I+2),P2C(I+2),J)=C*C
       C2H(P2C(I+2),P2C(I+4),J)=-2.0D0*C*S/DSQRT(3.0D0)
       C2H(P2C(I+3),P2C(I+3),J)=C
       C2H(P2C(I+3),P2C(I+5),J)=S
       C2H(P2C(I+4),P2C(I+1),J)=-C*S*DSQRT(3.0D0)
       C2H(P2C(I+4),P2C(I+2),J)=C*S*DSQRT(3.0D0)
       C2H(P2C(I+4),P2C(I+4),J)=(C*C-S*S)
       C2H(P2C(I+5),P2C(I+3),J)=-S
       C2H(P2C(I+5),P2C(I+5),J)=C
       P2H((I  ),(I  ),J)=1.0D0
       P2H((I+1),(I+1),J)=C*C
       P2H((I+1),(I+2),J)=S*S
       P2H((I+1),(I+4),J)=2.0D0*C*S
       P2H((I+2),(I+1),J)=S*S
       P2H((I+2),(I+2),J)=C*C
       P2H((I+2),(I+4),J)=-2.0D0*C*S
       P2H((I+3),(I+3),J)=C
       P2H((I+3),(I+5),J)=S
       P2H((I+4),(I+1),J)=-C*S
       P2H((I+4),(I+2),J)=C*S
       P2H((I+4),(I+4),J)=(C*C-S*S)
       P2H((I+5),(I+3),J)=-S
       P2H((I+5),(I+5),J)=C
       C2H_THETA(P2C(I  ),P2C(I  ),J)=0.0D0
       C2H_THETA(P2C(I+1),P2C(I+1),J)=-2.0D0*C*S*DFLOAT(J)
       C2H_THETA(P2C(I+1),P2C(I+2),J)=2.0D0*S*C*DFLOAT(J)
       C2H_THETA(P2C(I+1),P2C(I+4),J)=2.0D0*(C*C-S*S)/DSQRT(3.0D0)*DFLOAT(J)
       C2H_THETA(P2C(I+2),P2C(I+1),J)=2.0D0*C*S*DFLOAT(J)
       C2H_THETA(P2C(I+2),P2C(I+2),J)=-2.0D0*S*C*DFLOAT(J)
       C2H_THETA(P2C(I+2),P2C(I+4),J)=-2.0D0*(C*C-S*S)/DSQRT(3.0D0)*DFLOAT(J)
       C2H_THETA(P2C(I+3),P2C(I+3),J)=-S*DFLOAT(J)
       C2H_THETA(P2C(I+3),P2C(I+5),J)=C*DFLOAT(J)
       C2H_THETA(P2C(I+4),P2C(I+1),J)=-(C*C-S*S)*DSQRT(3.0D0)*DFLOAT(J)
       C2H_THETA(P2C(I+4),P2C(I+2),J)=(C*C-S*S)*DSQRT(3.0D0)*DFLOAT(J)
       C2H_THETA(P2C(I+4),P2C(I+4),J)=-4.0D0*C*S*DFLOAT(J)
       C2H_THETA(P2C(I+5),P2C(I+3),J)=-C*DFLOAT(J)
       C2H_THETA(P2C(I+5),P2C(I+5),J)=-S*DFLOAT(J)
       P2H_THETA((I  ),(I  ),J)=0.0D0
       P2H_THETA((I+1),(I+1),J)=-2.0D0*C*S*DFLOAT(J)
       P2H_THETA((I+1),(I+2),J)=2.0D0*S*C*DFLOAT(J)
       P2H_THETA((I+1),(I+4),J)=2.0D0*(C*C-S*S)*DFLOAT(J)
       P2H_THETA((I+2),(I+1),J)=2.0D0*C*S*DFLOAT(J)
       P2H_THETA((I+2),(I+2),J)=-2.0D0*S*C*DFLOAT(J)
       P2H_THETA((I+2),(I+4),J)=-2.0D0*(C*C-S*S)*DFLOAT(J)
       P2H_THETA((I+3),(I+3),J)=-S*DFLOAT(J)
       P2H_THETA((I+3),(I+5),J)=C*DFLOAT(J)
       P2H_THETA((I+4),(I+1),J)=-(C*C-S*S)*DFLOAT(J)
       P2H_THETA((I+4),(I+2),J)=(C*C-S*S)*DFLOAT(J)
       P2H_THETA((I+4),(I+4),J)=-4.0D0*C*S*DFLOAT(J)
       P2H_THETA((I+5),(I+3),J)=-C*DFLOAT(J)
       P2H_THETA((I+5),(I+5),J)=-S*DFLOAT(J)
      ENDIF
     ENDDO
     ! F ORBITALS
     DO I=1,NPGS
      IF ((PAX(I) == 3).AND.(PAY(I) == 0).AND.(PAZ(I) == 0)) THEN
       IF (P2C(I) == -1) CYCLE
       C2H(P2C(I  ),P2C(I  ),J)=1.0D0
       C2H(P2C(I+1),P2C(I+1),J)=C*C*C
       C2H(P2C(I+1),P2C(I+2),J)=S*S*S
       C2H(P2C(I+1),P2C(I+4),J)=3.0D0*C*C*S/DSQRT(5.0D0)
       C2H(P2C(I+1),P2C(I+7),J)=3.0D0*C*S*S/DSQRT(5.0D0)
       C2H(P2C(I+2),P2C(I+1),J)=-S*S*S
       C2H(P2C(I+2),P2C(I+2),J)=C*C*C
       C2H(P2C(I+2),P2C(I+4),J)=3.0D0*C*S*S/DSQRT(5.0D0)
       C2H(P2C(I+2),P2C(I+7),J)=-3.0D0*C*C*S/DSQRT(5.0D0)
       C2H(P2C(I+3),P2C(I+3),J)=C
       C2H(P2C(I+3),P2C(I+8),J)=S
       C2H(P2C(I+4),P2C(I+1),J)=-C*C*S*DSQRT(5.0D0)
       C2H(P2C(I+4),P2C(I+2),J)=C*S*S*DSQRT(5.0D0)
       C2H(P2C(I+4),P2C(I+4),J)=C*C*C-2.0D0*C*S*S
       C2H(P2C(I+4),P2C(I+7),J)=2.0D0*C*C*S-S*S*S
       C2H(P2C(I+5),P2C(I+5),J)=C*C
       C2H(P2C(I+5),P2C(I+6),J)=S*S
       C2H(P2C(I+5),P2C(I+9),J)=-2.0D0*C*S/DSQRT(3.0D0)
       C2H(P2C(I+6),P2C(I+5),J)=S*S
       C2H(P2C(I+6),P2C(I+6),J)=C*C
       C2H(P2C(I+6),P2C(I+9),J)=2.0D0*C*S/DSQRT(3.0D0)
       C2H(P2C(I+7),P2C(I+1),J)=C*S*S*DSQRT(5.0D0)
       C2H(P2C(I+7),P2C(I+2),J)=C*C*S*DSQRT(5.0D0)
       C2H(P2C(I+7),P2C(I+4),J)=S*S*S-2.0D0*C*C*S
       C2H(P2C(I+7),P2C(I+7),J)=C*C*C-2.0D0*C*S*S
       C2H(P2C(I+8),P2C(I+3),J)=-S
       C2H(P2C(I+8),P2C(I+8),J)=C
       C2H(P2C(I+9),P2C(I+5),J)=C*S*DSQRT(3.0D0)
       C2H(P2C(I+9),P2C(I+6),J)=-C*S*DSQRT(3.0D0)
       C2H(P2C(I+9),P2C(I+9),J)=C*C-S*S
       P2H((I  ),(I  ),J)=1.0D0
       P2H((I+1),(I+1),J)=C*C*C
       P2H((I+1),(I+2),J)=S*S*S
       P2H((I+1),(I+4),J)=3.0D0*C*C*S
       P2H((I+1),(I+7),J)=3.0D0*C*S*S
       P2H((I+2),(I+1),J)=-S*S*S
       P2H((I+2),(I+2),J)=C*C*C
       P2H((I+2),(I+4),J)=3.0D0*C*S*S
       P2H((I+2),(I+7),J)=-3.0D0*C*C*S
       P2H((I+3),(I+3),J)=C
       P2H((I+3),(I+8),J)=S
       P2H((I+4),(I+1),J)=-C*C*S
       P2H((I+4),(I+2),J)=C*S*S
       P2H((I+4),(I+4),J)=C*C*C-2.0D0*C*S*S
       P2H((I+4),(I+7),J)=2.0D0*C*C*S-S*S*S
       P2H((I+5),(I+5),J)=C*C
       P2H((I+5),(I+6),J)=S*S
       P2H((I+5),(I+9),J)=-2.0D0*C*S
       P2H((I+6),(I+5),J)=S*S
       P2H((I+6),(I+6),J)=C*C
       P2H((I+6),(I+9),J)=2.0D0*C*S
       P2H((I+7),(I+1),J)=C*S*S
       P2H((I+7),(I+2),J)=C*C*S
       P2H((I+7),(I+4),J)=S*S*S-2.0D0*C*C*S
       P2H((I+7),(I+7),J)=C*C*C-2.0D0*C*S*S
       P2H((I+8),(I+3),J)=-S
       P2H((I+8),(I+8),J)=C
       P2H((I+9),(I+5),J)=C*S
       P2H((I+9),(I+6),J)=-C*S
       P2H((I+9),(I+9),J)=C*C-S*S
       C2H_THETA(P2C(I  ),P2C(I  ),J)=0.0D0
       C2H_THETA(P2C(I+1),P2C(I+1),J)=-3.0D0*C*C*S*DFLOAT(J)
       C2H_THETA(P2C(I+1),P2C(I+2),J)=3.0D0*S*S*C*DFLOAT(J)
       C2H_THETA(P2C(I+1),P2C(I+4),J)=3.0D0*(-2.0D0*C*S*S+C*C*C)/DSQRT(5.0D0)*DFLOAT(J)
       C2H_THETA(P2C(I+1),P2C(I+7),J)=3.0D0*(-S*S*S+2.0D0*C*C*S)/DSQRT(5.0D0)*DFLOAT(J)
       C2H_THETA(P2C(I+2),P2C(I+1),J)=-3.0D0*S*S*C*DFLOAT(J)
       C2H_THETA(P2C(I+2),P2C(I+2),J)=-3.0D0*C*C*S*DFLOAT(J)
       C2H_THETA(P2C(I+2),P2C(I+4),J)=3.0D0*(2.0D0*C*C*S-S*S*S)/DSQRT(5.0D0)*DFLOAT(J)
       C2H_THETA(P2C(I+2),P2C(I+7),J)=-3.0D0*(-2.0D0*S*C*S+C*C*C)/DSQRT(5.0D0)*DFLOAT(J)
       C2H_THETA(P2C(I+3),P2C(I+3),J)=-S*DFLOAT(J)
       C2H_THETA(P2C(I+3),P2C(I+8),J)=C*DFLOAT(J)
       C2H_THETA(P2C(I+4),P2C(I+1),J)=-(-2.0D0*S*C*S+C*C*C)*DSQRT(5.0D0)*DFLOAT(J)
       C2H_THETA(P2C(I+4),P2C(I+2),J)=(-S*S*S+2.0D0*C*C*S)*DSQRT(5.0D0)*DFLOAT(J)
       C2H_THETA(P2C(I+4),P2C(I+4),J)=(-3.0D0*C*C*S-2.0D0*(-S*S*S+2.0D0*C*C*S))*DFLOAT(J)
       C2H_THETA(P2C(I+4),P2C(I+7),J)=(2.0D0*(-2.0D0*S*C*S+C*C*C)-3.0D0*C*S*S)*DFLOAT(J)
       C2H_THETA(P2C(I+5),P2C(I+5),J)=-2.0D0*C*S*DFLOAT(J)
       C2H_THETA(P2C(I+5),P2C(I+6),J)=2.0D0*S*C*DFLOAT(J)
       C2H_THETA(P2C(I+5),P2C(I+9),J)=-2.0D0*(C*C-S*S)/DSQRT(3.0D0)*DFLOAT(J)
       C2H_THETA(P2C(I+6),P2C(I+5),J)=2.0D0*S*C*DFLOAT(J)
       C2H_THETA(P2C(I+6),P2C(I+6),J)=-2.0D0*S*C*DFLOAT(J)
       C2H_THETA(P2C(I+6),P2C(I+9),J)=2.0D0*(C*C-S*S)/DSQRT(3.0D0)*DFLOAT(J)
       C2H_THETA(P2C(I+7),P2C(I+1),J)=(-S*S*S+2.0D0*C*C*S)*DSQRT(5.0D0)*DFLOAT(J)
       C2H_THETA(P2C(I+7),P2C(I+2),J)=(-2.0D0*S*C*S+C*C*C)*DSQRT(5.0D0)*DFLOAT(J)
       C2H_THETA(P2C(I+7),P2C(I+4),J)=(3.0D0*C*S*S-2.0D0*(-2.0D0*S*C*S+C*C*C))*DFLOAT(J)
       C2H_THETA(P2C(I+7),P2C(I+7),J)=(-3.0D0*S*C*C-2.0D0*(-S*S*S+2.0D0*C*C*S))*DFLOAT(J)
       C2H_THETA(P2C(I+8),P2C(I+3),J)=-C*DFLOAT(J)
       C2H_THETA(P2C(I+8),P2C(I+8),J)=-S*DFLOAT(J)
       C2H_THETA(P2C(I+9),P2C(I+5),J)=(C*C-S*S)*DSQRT(3.0D0)*DFLOAT(J)
       C2H_THETA(P2C(I+9),P2C(I+6),J)=-(C*C-S*S)*DSQRT(3.0D0)*DFLOAT(J)
       C2H_THETA(P2C(I+9),P2C(I+9),J)=-4.0D0*S*C*DFLOAT(J)
       P2H_THETA((I  ),(I  ),J)=0.0D0
       P2H_THETA((I+1),(I+1),J)=-3.0D0*C*C*S*DFLOAT(J)
       P2H_THETA((I+1),(I+2),J)=3.0D0*S*S*C*DFLOAT(J)
       P2H_THETA((I+1),(I+4),J)=3.0D0*(-2.0D0*C*S*S+C*C*C)*DFLOAT(J)
       P2H_THETA((I+1),(I+7),J)=3.0D0*(-S*S*S+2.0D0*C*C*S)*DFLOAT(J)
       P2H_THETA((I+2),(I+1),J)=-3.0D0*S*S*C*DFLOAT(J)
       P2H_THETA((I+2),(I+2),J)=-3.0D0*C*C*S*DFLOAT(J)
       P2H_THETA((I+2),(I+4),J)=3.0D0*(2.0D0*C*C*S-S*S*S)*DFLOAT(J)
       P2H_THETA((I+2),(I+7),J)=-3.0D0*(-2.0D0*S*C*S+C*C*C)*DFLOAT(J)
       P2H_THETA((I+3),(I+3),J)=-S*DFLOAT(J)
       P2H_THETA((I+3),(I+8),J)=C*DFLOAT(J)
       P2H_THETA((I+4),(I+1),J)=-(-2.0D0*S*C*S+C*C*C)*DFLOAT(J)
       P2H_THETA((I+4),(I+2),J)=(-S*S*S+2.0D0*C*C*S)*DFLOAT(J)
       P2H_THETA((I+4),(I+4),J)=(-3.0D0*C*C*S-2.0D0*(-S*S*S+2.0D0*C*C*S))*DFLOAT(J)
       P2H_THETA((I+4),(I+7),J)=(2.0D0*(-2.0D0*S*C*S+C*C*C)-3.0D0*C*S*S)*DFLOAT(J)
       P2H_THETA((I+5),(I+5),J)=-2.0D0*C*S*DFLOAT(J)
       P2H_THETA((I+5),(I+6),J)=2.0D0*S*C*DFLOAT(J)
       P2H_THETA((I+5),(I+9),J)=-2.0D0*(C*C-S*S)*DFLOAT(J)
       P2H_THETA((I+6),(I+5),J)=2.0D0*S*C*DFLOAT(J)
       P2H_THETA((I+6),(I+6),J)=-2.0D0*S*C*DFLOAT(J)
       P2H_THETA((I+6),(I+9),J)=2.0D0*(C*C-S*S)*DFLOAT(J)
       P2H_THETA((I+7),(I+1),J)=(-S*S*S+2.0D0*C*C*S)*DFLOAT(J)
       P2H_THETA((I+7),(I+2),J)=(-2.0D0*S*C*S+C*C*C)*DFLOAT(J)
       P2H_THETA((I+7),(I+4),J)=(3.0D0*C*S*S-2.0D0*(-2.0D0*S*C*S+C*C*C))*DFLOAT(J)
       P2H_THETA((I+7),(I+7),J)=(-3.0D0*S*C*C-2.0D0*(-S*S*S+2.0D0*C*C*S))*DFLOAT(J)
       P2H_THETA((I+8),(I+3),J)=-C*DFLOAT(J)
       P2H_THETA((I+8),(I+8),J)=-S*DFLOAT(J)
       P2H_THETA((I+9),(I+5),J)=(C*C-S*S)*DFLOAT(J)
       P2H_THETA((I+9),(I+6),J)=-(C*C-S*S)*DFLOAT(J)
       P2H_THETA((I+9),(I+9),J)=-4.0D0*S*C*DFLOAT(J)
      ENDIF
     ENDDO
    ENDDO
    IF ((IOPTN(9) >= 3).AND.(MYID == 0)) THEN
     WRITE(6,'(A)') 'CONVERSION MATRIX FOR CARTESIAN PRIMITIVES TO HELICAL PRIMITIVES'
     DO J=0,CEL1X
      WRITE(6,'(A,I3)') 'CELL: ',J
      CALL DUMP5(P2H(:,:,J),NPGS)
     ENDDO
     WRITE(6,'(A)') 'CONVERSION MATRIX FOR CARTESIAN CONTRACTED TO HELICAL CONTRACTED'
     DO J=0,CEL1X
      WRITE(6,'(A,I3)') 'CELL: ',J
      CALL DUMP5(C2H(:,:,J),NCGS)
     ENDDO
    ENDIF
    DEALLOCATE(P2C)
   ENDIF

   RETURN
END SUBROUTINE
