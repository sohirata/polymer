SUBROUTINE RESTORE_DENSITYMATRIX
! RESTORE A DENSITY MATRIX FROM CHECKPOINT FILE (30).
! CHEKPOINT FILE IS COPTN(1)//'.chk'.

   USE MPI_F08
   USE CONSTANTS
   USE CONTROL
   USE BASISSET
   USE INTEGRAL

   IMPLICIT NONE
!  INCLUDE "mpif.h"
   INTEGER :: QX,QY,QZ,I,JX,JY,JZ,K,L,EOF
   INTEGER :: ICASE
   DOUBLE PRECISION :: A,B

! --- MYID==0 ONLY
   IF (MYID == 0) THEN
    P_C=0.0D0
    L=0
    OPEN(30,FILE=TRIM(COPTN(1))//'.chk',FORM='FORMATTED')
    REWIND(30)
    DO
     READ(30,*,IOSTAT=EOF) I,JX,JY,JZ
     IF ((EOF /= 0).AND.(L == 0)) THEN
      CLOSE(30)
      WRITE(6,'(A)') 'NO DENSITY MATRIX FOUND IN .chk FILE'
      ICASE=1
      EXIT
     ELSE IF ((EOF /= 0).AND.(L > 0)) THEN
      CLOSE(30)
      WRITE(6,'(A)') '****** DENSITY MATRIX RESTORED FROM .chk FILE ******'
      IF (IOPTN(9) == 3) THEN
       WRITE(6,'(A)') 'INITIAL DENSITY MATRIX'
       CALL DUMP1(P_C,NCGS,CEL1X,CEL1Y,CEL1Z)
      ENDIF
      ICASE=2
      EXIT
     ELSE IF (I /= NCGS) THEN
      CLOSE(30)
      WRITE(6,'(A)') 'DENSITY MATRIX IN INCOMPATIBLE FORMAT'
      ICASE=1
      EXIT
     ELSE
      DO QX=-JX,JX
      DO QY=-JY,JY
      DO QZ=-JZ,JZ
       DO I=1,NCGS
        DO K=1,NCGS
         IF ((QX >= -CEL1X).AND.(QX <= CEL1X).AND. &
             (QY >= -CEL1Y).AND.(QY <= CEL1Y).AND. &
             (QZ >= -CEL1Z).AND.(QZ <= CEL1Z)) THEN
          READ(30,*) P_C(K,I,QX,QY,QZ),B
         ELSE
          READ(30,*) A,B
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      ENDDO
      ENDDO
      L=L+1
     ENDIF
    ENDDO
   ENDIF
! --- MYID==0 ONLY END
   CALL MPI_BCAST(ICASE,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  write(*,*) 'MYID=',MYID,'ICASE=',ICASE
   IF (ICASE == 1) THEN
    CALL GUESS_DENSITYMATRIX
   ELSE IF (ICASE == 2) THEN
    CALL MPI_BCAST(P_C,NCGS**2*(CEL1X*2+1)*(CEL1Y*2+1)*(CEL1Z*2+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   ENDIF
   IF ((IOPTN(9) >= 2).AND.(MYID == 0)) THEN
    WRITE(6,'(A)') 'INITIAL DENSITY MATRIX READ FROM .chk FILE'
    CALL DUMP1(P_C,NCGS,CEL1X,CEL1Y,CEL1Z)
   ENDIF
   RETURN
END SUBROUTINE



SUBROUTINE GUESS_DENSITYMATRIX
! CONSTRUCT INITIAL DENSITY MATRIX GUESS BY SIMPLY
! SCATTERING THE CHARGE OF EACH NUCLEUS TO THE BASIS
! FUNCTIONS WHICH BELONG TO IT.

   USE MPI_F08
   USE CONSTANTS
   USE CONTROL
   USE STRUCTURE
   USE BASISSET
   USE INTEGRAL

   IMPLICIT NONE
!  INCLUDE "mpif.h"
   INTEGER :: I,J,K,L,M
   DOUBLE PRECISION :: A
   LOGICAL :: LDUMMY

   P_C=0.0D0
   DO I=1,NATOM
    J=IATOM(I)
    L=0
    DO K=1,NCGS
     LDUMMY=.TRUE.
     DO M=1,NPGS
      IF (CC(K,M) /= 0.0D0) LDUMMY=.FALSE.
     ENDDO
     IF ((.NOT.LDUMMY).AND.(CGSX(K) == ATOMX(I)).AND.(CGSY(K) == ATOMY(I)).AND.(CGSZ(K) == ATOMZ(I))) L=L+1
    ENDDO
    IF (L == 0) CALL PABORT('AN INTERNAL PROGRAM ERROR IS DETECTED IN GUESS_DENSITYMATRIX')
    A=DFLOAT(J)/DFLOAT(L)
    DO K=1,NCGS
     LDUMMY=.TRUE.
     DO M=1,NPGS
      IF (CC(K,M) /= 0.0D0) LDUMMY=.FALSE.
     ENDDO
     IF ((.NOT.LDUMMY).AND.(CGSX(K) == ATOMX(I)).AND.(CGSY(K) == ATOMY(I)).AND.(CGSZ(K) == ATOMZ(I))) P_C(K,K,0,0,0)=A
    ENDDO
   ENDDO
   IF ((IOPTN(9) == 2).AND.(MYID == 0)) THEN
    WRITE(6,'(A)') 'INITIAL DENSITY MATRIX'
    CALL DUMP1(P_C,NCGS,CEL1X,CEL1Y,CEL1Z)
   ENDIF
   RETURN
END SUBROUTINE



SUBROUTINE STORE_DENSITYMATRIX(ISCF,DCHG)
! STORE THE LAST DENSITY MATRIX AND DENSITY MATRIX CHANGE.
! IN THE CHECKPOINT FILE.  THE MEAN SQUARE CHANGE IS RETURNED AS AN ARGUMENT.

   USE MPI_F08
   USE CONSTANTS
   USE CONTROL
   USE BASISSET
   USE INTEGRAL

   IMPLICIT NONE
!  INCLUDE "mpif.h"
   INTEGER :: ISCF,QX,QY,QZ,I,J,K
   DOUBLE PRECISION :: A,B,DCHG

   DCHG=0.0D0
   DO QX=-CEL1X,CEL1X
   DO QY=-CEL1Y,CEL1Y
   DO QZ=-CEL1Z,CEL1Z
    DO I=1,NCGS
     DO J=1,NCGS
      DCHG=DCHG+(P_C_OUT(J,I,QX,QY,QZ)-P_C(J,I,QX,QY,QZ))**2
     ENDDO
    ENDDO
   ENDDO
   ENDDO
   ENDDO
   PULAY(ISCF,ISCF)=DCHG
   DCHG=DSQRT(DCHG/DFLOAT((2*CEL1X+1)*(2*CEL1Y+1)*(2*CEL1Z+1)*NCGS**2))

   IF (MYID == 0) THEN
    OPEN(30,FILE=TRIM(COPTN(1))//'.chk',FORM='FORMATTED')
    REWIND(30)
    IF (ISCF >= 2) THEN
     DO I=1,ISCF-1
      PULAY(I,ISCF)=0.0D0
      READ(30,*) J,K
      DO QX=-CEL1X,CEL1X
      DO QY=-CEL1Y,CEL1Y
      DO QZ=-CEL1Z,CEL1Z
       DO J=1,NCGS
        DO K=1,NCGS
         READ(30,*) A,B
         PULAY(I,ISCF)=PULAY(I,ISCF)+B*(P_C_OUT(K,J,QX,QY,QZ)-P_C(K,J,QX,QY,QZ))
        ENDDO
       ENDDO
      ENDDO
      ENDDO
      ENDDO
      PULAY(ISCF,I)=PULAY(I,ISCF)
     ENDDO
    ENDIF
!   DUMP PULAY MATRIX
    IF (IOPTN(9) >= 3)THEN
     WRITE(6,'(A)') 'PULAY MATRIX'
     DO I=1,ISCF
      WRITE(6,'(50F10.5:)') (PULAY(J,I),J=1,ISCF)
     ENDDO
    ENDIF
    WRITE(30,*) NCGS,CEL1X,CEL1Y,CEL1Z
    DO QX=-CEL1X,CEL1X
    DO QY=-CEL1Y,CEL1Y
    DO QZ=-CEL1Z,CEL1Z
     DO I=1,NCGS
      DO J=1,NCGS
       WRITE(30,*) P_C(J,I,QX,QY,QZ),P_C_OUT(J,I,QX,QY,QZ)-P_C(J,I,QX,QY,QZ)
      ENDDO
     ENDDO
    ENDDO
    ENDDO
    ENDDO
    CLOSE(30)
   ENDIF

   CALL MPI_BCAST(PULAY,IOPTN(13)**2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

   RETURN
END SUBROUTINE



SUBROUTINE RELAX_DENSITYMATRIX(COEF)
! CONSTRUCT A NEW DENSITY MATRIX BY SIMPLE RELAXATION.
! P_C_OUT MATRIX CONTAINS THE DIFFERECE DENSITY MATRIX.

   USE MPI_F08
   USE CONSTANTS
   USE CONTROL
   USE BASISSET
   USE INTEGRAL

   IMPLICIT NONE
!  INCLUDE "mpif.h"
   INTEGER :: QX,QY,QZ,I,J
   DOUBLE PRECISION :: COEF
   DOUBLE PRECISION :: A,B,C

   DO QX=-CEL1X,CEL1X
   DO QY=-CEL1Y,CEL1Y
   DO QZ=-CEL1Z,CEL1Z
    DO I=1,NCGS
     DO J=1,NCGS
      A=P_C(J,I,QX,QY,QZ)
      B=P_C_OUT(J,I,QX,QY,QZ)
      C=COEF*B+(1.0D0-COEF)*A
      P_C(J,I,QX,QY,QZ)=C
      P_C_OUT(J,I,QX,QY,QZ)=C-A
     ENDDO
    ENDDO
   ENDDO
   ENDDO
   ENDDO
!  DUMP NEXT DENSITY MATRIX
   IF ((IOPTN(9) == 2).AND.(MYID == 0)) THEN
    WRITE(6,'(A)') 'NEXT DENSITY MATRIX'
    CALL DUMP1(P_C,NCGS,CEL1X,CEL1Y,CEL1Z)
   ENDIF
   RETURN
END SUBROUTINE



SUBROUTINE DIIS(ISCF,ITERM)
! CONSTRUCT A NEW DENSITY MATRIX BY DIIS EXTRAPOLATION.

   USE MPI_F08
   USE CONSTANTS
   USE CONTROL
   USE BASISSET
   USE INTEGRAL

   IMPLICIT NONE
!  INCLUDE "mpif.h"
   INTEGER :: ISCF,ITERM,QX,QY,QZ,I,J,K
   INTEGER,ALLOCATABLE :: INDX(:)
   DOUBLE PRECISION :: A,D,E
   DOUBLE PRECISION,ALLOCATABLE :: B(:,:),BS(:,:),C(:),CS(:)
   
   ALLOCATE(INDX(IOPTN(13)),B(IOPTN(13),IOPTN(13)),BS(IOPTN(13),IOPTN(13)),C(IOPTN(13)),CS(IOPTN(13)))

   DO I=1,ITERM
    DO J=1,ITERM
     B(J,I)=PULAY(ISCF-J+1,ISCF-I+1)
    ENDDO
    B(I,ITERM+1)=-1.0D0
    B(ITERM+1,I)=-1.0D0
    C(I)=0.0D0
   ENDDO
   B(ITERM+1,ITERM+1)=0.0D0
   C(ITERM+1)=-1.0D0
   DO I=1,ITERM+1
    DO J=1,ITERM+1
     BS(J,I)=B(J,I)
    ENDDO
    CS(I)=C(I)
   ENDDO

   CALL LUDCMP(B,ITERM+1,IOPTN(13),INDX,D)
   CALL LUBKSB(B,ITERM+1,IOPTN(13),INDX,C)
   CALL MPROVE(BS,B,ITERM+1,IOPTN(13),INDX,CS,C)

   P_C_OUT=0.0D0

   IF (MYID == 0) OPEN(30,FILE=TRIM(COPTN(1))//'.chk',FORM='FORMATTED')
   IF (MYID == 0) REWIND(30)
   DO I=1,ISCF
    IF (MYID == 0) READ(30,*) J,K
    CALL MPI_BCAST(J,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    CALL MPI_BCAST(K,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    DO QX=-CEL1X,CEL1X
    DO QY=-CEL1Y,CEL1Y
    DO QZ=-CEL1Z,CEL1Z
     DO J=1,NCGS
      DO K=1,NCGS
       IF (MYID == 0) READ(30,*) A,D
       CALL MPI_BCAST(A,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(D,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
       IF (ISCF-I+1 <= ITERM) P_C_OUT(K,J,QX,QY,QZ)=P_C_OUT(K,J,QX,QY,QZ)+C(ISCF-I+1)*A
      ENDDO
     ENDDO
    ENDDO
    ENDDO
    ENDDO
   ENDDO

   DO QX=-CEL1X,CEL1X
   DO QY=-CEL1Y,CEL1Y
   DO QZ=-CEL1Z,CEL1Z
    DO I=1,NCGS
     DO J=1,NCGS
      E=P_C(J,I,QX,QY,QZ)
      P_C(J,I,QX,QY,QZ)=P_C_OUT(J,I,QX,QY,QZ)
      P_C_OUT(J,I,QX,QY,QZ)=P_C(J,I,QX,QY,QZ)-E
     ENDDO
    ENDDO
   ENDDO
   ENDDO
   ENDDO

!  DUMP NEXT DENSITY MATRIX
   IF ((IOPTN(9) == 2).AND.(MYID == 0)) THEN
    WRITE(6,'(A)') 'NEXT DENSITY MATRIX'
    CALL DUMP1(P_C,NCGS,CEL1X,CEL1Y,CEL1Z)
   ENDIF
   DEALLOCATE(INDX,B,BS,C,CS)
   IF (MYID == 0) CLOSE(30)

   RETURN
END SUBROUTINE



SUBROUTINE DECONTRACT_DENSITYMATRIX
! CONSTRUCT DENSITY MATRIX AND ENERGY-WEIGHTED DENSITY MATRIX OVER PRIMITIVE GAUSSIANS.

   USE MPI_F08
   USE CONSTANTS
   USE CONTROL
   USE BASISSET
   USE INTEGRAL
   
   IMPLICIT NONE
!  INCLUDE "mpif.h"
   INTEGER :: QX,QY,QZ
   DOUBLE PRECISION,ALLOCATABLE :: W1(:,:),W2(:,:),W3(:,:)

   ALLOCATE(W1(NCGS,NCGS),W2(NCGS,NPGS),W3(NPGS,NPGS))
   DO QX=-CEL1X,CEL1X
   DO QY=-CEL1Y,CEL1Y
   DO QZ=-CEL1Z,CEL1Z
    W1=P_C(:,:,QX,QY,QZ)
    W2=MATMUL(W1,CC)
    W3=MATMUL(TRANSPOSE(CC),W2)
    P_P(:,:,QX,QY,QZ)=W3
    W1=W_C(:,:,QX,QY,QZ)
    W2=MATMUL(W1,CC)
    W3=MATMUL(TRANSPOSE(CC),W2)
    W_P(:,:,QX,QY,QZ)=W3
   ENDDO
   ENDDO
   ENDDO
   DEALLOCATE(W1,W2,W3)
!  DUMP DENSITY MATRIX OVER PRIMITIVES
   IF ((IOPTN(9) == 3).AND.(MYID == 0)) THEN
    WRITE(6,'(A)') 'DENSITY MATRIX OVER PRIMITIVE GAUSSIANS'
    CALL DUMP1(P_P,NPGS,CEL1X,CEL1Y,CEL1Z)
    WRITE(6,'(A)') 'ENERGY-WEIGHTED DENSITY MATRIX OVER PRIMITIVE GAUSSIANS'
    CALL DUMP1(W_P,NPGS,CEL1X,CEL1Y,CEL1Z)
   ENDIF
   RETURN
END SUBROUTINE



SUBROUTINE LUDCMP(A,N,NP,INDX,D)
! LU DECOMPOSITION FOR REAL MATRICES.
! SEE NUMERICAL RECIPES IN FORTRAN 2ND ED. P.38.

   IMPLICIT NONE
   INTEGER,PARAMETER :: NMAX = 1000
   REAL,PARAMETER :: TINY = 1.0E-20
   INTEGER :: N,NP,INDX(N)
   INTEGER :: I,IMAX,J,K
   DOUBLE PRECISION :: D,A(NP,NP)
   DOUBLE PRECISION :: AAMAX,DUM,SUM,VV(NMAX)

   D=1.0D0
   DO 10 I=1,N
    AAMAX=0.0D0
    DO 11 J=1,N
     IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
11  CONTINUE
    IF (AAMAX.EQ.0.0D0) CALL PABORT('SINGULAR MATRIX IN LU DECOMPOSITION')
    VV(I)=1.0D0/AAMAX
10 CONTINUE
   DO 12 J=1,N
    DO 13 I=1,J-1
     SUM=A(I,J)
     DO 14 K=1,I-1
      SUM=SUM-A(I,K)*A(K,J)
14   CONTINUE
     A(I,J)=SUM
13  CONTINUE
    AAMAX=0.0D0
    DO 15 I=J,N
     SUM=A(I,J)
     DO 16 K=1,J-1
      SUM=SUM-A(I,K)*A(K,J)
16   CONTINUE
     A(I,J)=SUM
     DUM=VV(I)*DABS(SUM)
     IF (DUM.GE.AAMAX) THEN
      IMAX=I
      AAMAX=DUM
     END IF
15  CONTINUE
    IF (J.NE.IMAX) THEN
     DO 17 K=1,N
      DUM=A(IMAX,K)
      A(IMAX,K)=A(J,K)
      A(J,K)=DUM
17   CONTINUE
     D=-D
     VV(IMAX)=VV(J)
    END IF
    INDX(J)=IMAX
    IF (A(J,J).EQ.0.0D0) A(J,J)=TINY
    IF (J.NE.N) THEN
     DUM=1.0D0/A(J,J)
     DO 18 I=J+1,N
      A(I,J)=A(I,J)*DUM
18   CONTINUE
    END IF
12 CONTINUE

   RETURN
END SUBROUTINE



SUBROUTINE LUBKSB(A,N,NP,INDX,B)
! FORWARD AND BACKSUBSTITUTION FOR LINEAR EQUATIONS AX=B.
! SEE NUMERICAL RECIPES IN FORTRAN 2ND ED. P.39.

   IMPLICIT NONE
   INTEGER :: N,NP,INDX(N)
   INTEGER :: I,II,J,LL
   DOUBLE PRECISION :: A(NP,NP),B(N)
   DOUBLE PRECISION :: SUM

   II=0
   DO 10 I=1,N
    LL=INDX(I)
    SUM=B(LL)
    B(LL)=B(I)
    IF (II.NE.0) THEN
     DO 11 J=II,I-1
      SUM=SUM-A(I,J)*B(J)
11   CONTINUE
    ELSE IF (SUM.NE.0.0D0) THEN
     II=I
    END IF
    B(I)=SUM
10 CONTINUE
   DO 12 I=N,1,-1
    SUM=B(I)
    DO 13 J=I+1,N
     SUM=SUM-A(I,J)*B(J)
13  CONTINUE
    B(I)=SUM/A(I,I)
12 CONTINUE

   RETURN
END SUBROUTINE



SUBROUTINE MPROVE(A,ALUD,N,NP,INDX,B,X)
! ITERATIVE IMPROVEMENT OF SOLUTION VECTORS FOR LINEAR EQUATIONS AX=B.
! SEE NUMERICAL RECIPES IN FORTRAN 2ND ED. P.48.

   IMPLICIT NONE
   INTEGER,PARAMETER :: NMAX = 1000
   INTEGER :: N,NP,INDX(N)
   INTEGER :: I,J
   DOUBLE PRECISION :: A(NP,NP),ALUD(NP,NP),B(N),X(N)
   DOUBLE PRECISION :: R(NMAX),SDP

   DO 10 I=1,N
    SDP=-B(I)
    DO 11 J=1,N
     SDP=SDP+A(I,J)*X(J)
11  CONTINUE
    R(I)=SDP
10 CONTINUE
   CALL LUBKSB(ALUD,N,NP,INDX,R)
   DO 12 I=1,N
    X(I)=X(I)-R(I)
12 CONTINUE

   RETURN
END SUBROUTINE



SUBROUTINE ROTATE_DENSITYMATRIX

   USE MPI_F08
   USE CONSTANTS
   USE CONTROL
   USE STRUCTURE
   USE BASISSET
   USE INTEGRAL
   USE GRADIENT

   INTEGER :: Q1X,Q2X
   INTEGER :: I,J,II,JJ
   DOUBLE PRECISION,ALLOCATABLE :: TMP1(:,:),TMP2(:,:),TMP3(:,:),TMP4(:,:)

   P_P_HELIX=0.0D0
   W_P_HELIX=0.0D0
   P_P_HELIX_THETA=0.0D0
   W_P_HELIX_THETA=0.0D0
   ALLOCATE(TMP1(NPGS,NPGS),TMP2(NPGS,NPGS),TMP3(NPGS,NPGS),TMP4(NPGS,NPGS))
   DO Q2X=-CEL2X,CEL2X
    DO Q1X=-CEL2X-CEL1X,CEL2X+CEL1X
     IF ((Q1X-Q2X >= -CEL1X).AND.(Q1X-Q2X <= CEL1X)) THEN
      DO I=1,NPGS
       DO J=1,NPGS
        TMP1(I,J)=0.0D0
        TMP2(I,J)=0.0D0
        TMP3(I,J)=0.0D0
        TMP4(I,J)=0.0D0
        DO II=1,NPGS
         TMP1(I,J)=TMP1(I,J)+P2H(II,I,Q2X)*P_P(II,J,Q1X-Q2X,0,0)
         TMP2(I,J)=TMP2(I,J)+P2H(II,I,Q2X)*W_P(II,J,Q1X-Q2X,0,0)
         TMP3(I,J)=TMP3(I,J)+P2H_THETA(II,I,Q2X)*P_P(II,J,Q1X-Q2X,0,0)
         TMP4(I,J)=TMP4(I,J)+P2H_THETA(II,I,Q2X)*W_P(II,J,Q1X-Q2X,0,0)
        ENDDO
       ENDDO
      ENDDO
      DO I=1,NPGS
       DO J=1,NPGS
        P_P_HELIX(I,J,Q2X,Q1X)=0.0D0
        W_P_HELIX(I,J,Q2X,Q1X)=0.0D0
        P_P_HELIX_THETA(I,J,Q2X,Q1X)=0.0D0
        W_P_HELIX_THETA(I,J,Q2X,Q1X)=0.0D0
        DO JJ=1,NPGS
         P_P_HELIX(I,J,Q2X,Q1X)=P_P_HELIX(I,J,Q2X,Q1X)+P2H(JJ,J,Q1X)*TMP1(I,JJ)
         W_P_HELIX(I,J,Q2X,Q1X)=W_P_HELIX(I,J,Q2X,Q1X)+P2H(JJ,J,Q1X)*TMP2(I,JJ)
         P_P_HELIX_THETA(I,J,Q2X,Q1X)=P_P_HELIX_THETA(I,J,Q2X,Q1X)+P2H_THETA(JJ,J,Q1X)*TMP1(I,JJ)+P2H(JJ,J,Q1X)*TMP3(I,JJ)
         W_P_HELIX_THETA(I,J,Q2X,Q1X)=W_P_HELIX_THETA(I,J,Q2X,Q1X)+P2H_THETA(JJ,J,Q1X)*TMP2(I,JJ)+P2H(JJ,J,Q1X)*TMP4(I,JJ)
        ENDDO
       ENDDO
      ENDDO
     ENDIF
    ENDDO
   ENDDO
   DEALLOCATE(TMP1,TMP2,TMP3,TMP4)

   P_C_HELIX=0.0D0
   W_C_HELIX=0.0D0
   P_C_HELIX_THETA=0.0D0
   W_C_HELIX_THETA=0.0D0
   ALLOCATE(TMP1(NCGS,NCGS),TMP2(NCGS,NCGS),TMP3(NCGS,NCGS),TMP4(NCGS,NCGS))
   DO Q2X=-CEL2X,CEL2X
    DO Q1X=-CEL2X-CEL1X,CEL2X+CEL1X
     IF ((Q1X-Q2X >= -CEL1X).AND.(Q1X-Q2X <= CEL1X)) THEN
      DO I=1,NCGS
       DO J=1,NCGS
        TMP1(I,J)=0.0D0
        TMP2(I,J)=0.0D0
        TMP3(I,J)=0.0D0
        TMP4(I,J)=0.0D0
        DO II=1,NCGS
         TMP1(I,J)=TMP1(I,J)+C2H(II,I,Q2X)*P_C(II,J,Q1X-Q2X,0,0)
         TMP2(I,J)=TMP2(I,J)+C2H(II,I,Q2X)*W_C(II,J,Q1X-Q2X,0,0)
         TMP3(I,J)=TMP3(I,J)+C2H_THETA(II,I,Q2X)*P_C(II,J,Q1X-Q2X,0,0)
         TMP4(I,J)=TMP4(I,J)+C2H_THETA(II,I,Q2X)*W_C(II,J,Q1X-Q2X,0,0)
        ENDDO
       ENDDO
      ENDDO
      DO I=1,NCGS
       DO J=1,NCGS
        P_C_HELIX(I,J,Q2X,Q1X)=0.0D0
        W_C_HELIX(I,J,Q2X,Q1X)=0.0D0
        P_C_HELIX_THETA(I,J,Q2X,Q1X)=0.0D0
        W_C_HELIX_THETA(I,J,Q2X,Q1X)=0.0D0
        DO JJ=1,NCGS
         P_C_HELIX(I,J,Q2X,Q1X)=P_C_HELIX(I,J,Q2X,Q1X)+C2H(JJ,J,Q1X)*TMP1(I,JJ)
         W_C_HELIX(I,J,Q2X,Q1X)=W_C_HELIX(I,J,Q2X,Q1X)+C2H(JJ,J,Q1X)*TMP2(I,JJ)
         P_C_HELIX_THETA(I,J,Q2X,Q1X)=P_C_HELIX_THETA(I,J,Q2X,Q1X)+C2H_THETA(JJ,J,Q1X)*TMP1(I,JJ)+C2H(JJ,J,Q1X)*TMP3(I,JJ)
         W_C_HELIX_THETA(I,J,Q2X,Q1X)=W_C_HELIX_THETA(I,J,Q2X,Q1X)+C2H_THETA(JJ,J,Q1X)*TMP2(I,JJ)+C2H(JJ,J,Q1X)*TMP4(I,JJ)
        ENDDO
       ENDDO
      ENDDO
     ENDIF
    ENDDO
   ENDDO
   DEALLOCATE(TMP1,TMP2,TMP3,TMP4)

   RETURN
END SUBROUTINE
