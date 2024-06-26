SUBROUTINE FOURINDEX_ERI_GRADIENT
! CALCULATE FOUR-CENTER TWO-ELECTRON ELECTRON REPULSION INTEGRAL DERIVATIVES
! USING THE OBARA-SAIKA RECURSION METHOD.

   USE MPI_F08
   USE CONSTANTS
   USE CONTROL
   USE STRUCTURE
   USE BASISSET
   USE INTEGRAL
   USE GRADIENT
   USE FMT
   
   IMPLICIT NONE
!  INCLUDE "mpif.h"
   DOUBLE PRECISION :: ICPUS,ICPUE
   INTEGER :: Q1X,Q1Y,Q1Z,Q2X,Q2Y,Q2Z,Q3X,Q3Y,Q3Z
   INTEGER :: I,J,K,L,II,JJ,KK,LL,IJ,KL,M
   INTEGER :: I1,I2,I3,J1,J2,J3,IA,JA
   INTEGER :: IJT,KLT
   INTEGER :: NPSHELL_SQ
   INTEGER :: IJANG,SIJ,SKL,NIJ,NIJ_G,IJKLANG,TS
   INTEGER :: DOWN1(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12,(BASIS_LMAX+1)**2,3)
   INTEGER :: DOWN2(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12,(BASIS_LMAX+1)**2,3)
   INTEGER :: DOWN3(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12,(BASIS_LMAX+1)**2)
   INTEGER :: DOWN4(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12,(BASIS_LMAX+1)**2)
   INTEGER :: DOWN5(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12,(BASIS_LMAX+1)**2)
   INTEGER :: DOWNXYZ(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12,(BASIS_LMAX+1)**2)
   INTEGER :: UP1(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/36,(BASIS_LMAX+1)**2,3)
   INTEGER :: UP2(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/36,(BASIS_LMAX+1)**2,3)
   INTEGER :: NKL((BASIS_LMAX+1)**2),NKL_G((BASIS_LMAX+1)**2),NKL2
   INTEGER :: PGSI(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12)
   INTEGER :: PGSJ(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12)
   INTEGER :: IJMAP(0:BASIS_LMAX+1,0:BASIS_LMAX+1,0:BASIS_LMAX+1,0:BASIS_LMAX+1,0:BASIS_LMAX+1,0:BASIS_LMAX+1)
   INTEGER :: KLANGINDX(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12,(BASIS_LMAX+1)**2)
   INTEGER :: IATM,JATM
   INTEGER,ALLOCATABLE :: KLT_S(:),KLANG(:),PGSK(:,:),PGSL(:,:),KATM(:),LATM(:)
   DOUBLE PRECISION :: ERI(0:(BASIS_LMAX+1)**2,((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12, &
                                               ((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12)
   DOUBLE PRECISION :: PI_AI(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12),IJIJ
   DOUBLE PRECISION :: PISUB,F1(0:(BASIS_LMAX+1)*4),F2(0:(BASIS_LMAX+1)*4),WX,WY,WZ,WI_PI(3),WI_QI(3)
   DOUBLE PRECISION :: DOWN1_C(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12,(BASIS_LMAX+1)**2,3)
   DOUBLE PRECISION :: DOWN2_C(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12,(BASIS_LMAX+1)**2,3)
   DOUBLE PRECISION :: DOWN4_C(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12,(BASIS_LMAX+1)**2)
   DOUBLE PRECISION :: DOWN5_C(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12,(BASIS_LMAX+1)**2)
   DOUBLE PRECISION :: P1X,P1Y,P1Z,P2X,P2Y,P2Z,P1C,P1S,P2C,P2S
   DOUBLE PRECISION :: AX,AY,AZ,BX,BY,BZ,CX,CY,CZ,DX,DY,DZ,PX,PY,PZ,QX,QY,QZ
   DOUBLE PRECISION :: COMZ,COMZ2,C1,C2,LARGECOMZ,RHO_COMZ1,RHO_COMZ2,RHO,T,H
   DOUBLE PRECISION,ALLOCATABLE :: COMZ_S(:),PX_S(:),PY_S(:),PZ_S(:),OV(:),QI_CI(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: ZTK(:),ZTL(:),KLKL(:)
   DOUBLE PRECISION,ALLOCATABLE :: P3X(:),P3Y(:),P3Z(:),P3C(:),P3S(:)
   DOUBLE PRECISION,ALLOCATABLE :: MAXC(:)
   DOUBLE PRECISION :: A,FN,GN1,GN2,ANGLE
   DOUBLE PRECISION :: XI1,YI1,ZI1,XJ1,YJ1,ZJ1,XK1,YK1,ZK1,XL1,YL1,ZL1
   LOGICAL,ALLOCATABLE :: LCEL1X(:),LCEL1Y(:),LCEL1Z(:)
   DOUBLE PRECISION,ALLOCATABLE :: FORCEX_TMP(:),FORCEY_TMP(:),FORCEZ_TMP(:)
   DOUBLE PRECISION :: FORCETX_TMP,FORCETY_TMP,FORCETZ_TMP,FORCETA_TMP
   INTEGER :: TICKET

   NPSHELL_SQ=NPSHELL**2
   ALLOCATE(P3X(-CEL2X-CEL1X:CEL2X+CEL1X),P3Y(-CEL2Y-CEL1Y:CEL2Y+CEL1Y),P3Z(-CEL2Z-CEL1Z:CEL2Z+CEL1Z))
   ALLOCATE(P3C(-CEL2X-CEL1X:CEL2X+CEL1X),P3S(-CEL2X-CEL1X:CEL2X+CEL1X))
   ALLOCATE(COMZ_S(NPSHELL_SQ),PX_S(NPSHELL_SQ),PY_S(NPSHELL_SQ),PZ_S(NPSHELL_SQ),OV(NPSHELL_SQ))
   ALLOCATE(QI_CI(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12,NPSHELL_SQ))
   ALLOCATE(KLT_S(NPSHELL_SQ),KLANG(NPSHELL_SQ),PGSK(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12,NPSHELL_SQ))
   ALLOCATE(PGSL(((BASIS_LMAX+1)*(BASIS_LMAX+2)*(BASIS_LMAX+3))**2/12,NPSHELL_SQ),KLKL(NPSHELL_SQ))
   ALLOCATE(KATM(NPSHELL_SQ),LATM(NPSHELL_SQ),LCEL1X(-CEL2X:CEL2X),LCEL1Y(-CEL2Y:CEL2Y),LCEL1Z(-CEL2Z:CEL2Z))
   ALLOCATE(MAXC(NPGS),ZTK(NPSHELL_SQ),ZTL(NPSHELL_SQ))
   ALLOCATE(FORCEX_TMP(NATOM),FORCEY_TMP(NATOM),FORCEZ_TMP(NATOM))

   IF (MYID == 0) WRITE(6,'(A)') 'CALCULATE TWO-ELECTRON INTEGRAL DERIVATIVES'

   FORCEX_TMP=0.0D0
   FORCEY_TMP=0.0D0
   FORCEZ_TMP=0.0D0
   FORCETX_TMP=0.0D0
   FORCETY_TMP=0.0D0
   FORCETZ_TMP=0.0D0
   FORCETA_TMP=0.0D0
  
   ! PRELIMINARY & ALLOCATION
   PISUB=2.0D0/DSQRT(PI)
   F1(0)=1.0D0/PISUB
   F2(0)=-0.5D0
   DO I=1,(BASIS_LMAX+1)*4
    F1(I)=F1(I-1)*DFLOAT(2*I-1)/2.0D0
    F2(I)=F2(I-1)-1.0D0
   ENDDO
   DO Q3X=-CEL2X-CEL1X,CEL2X+CEL1X
    ANGLE=DFLOAT(Q3X)*HELIX
    P3X(Q3X)=DFLOAT(Q3X)*PERIODX
    P3C(Q3X)=DCOS(ANGLE)
    P3S(Q3X)=DSIN(ANGLE)
   ENDDO
   DO Q3Y=-CEL2Y-CEL1Y,CEL2Y+CEL1Y
    P3Y(Q3Y)=DFLOAT(Q3Y)*PERIODY
   ENDDO
   DO Q3Z=-CEL2Z-CEL1Z,CEL2Z+CEL1Z
    P3Z(Q3Z)=DFLOAT(Q3Z)*PERIODZ
   ENDDO
   MAXC=0.0D0
   DO I=1,NPGS
    DO J=1,NCGS
     IF (DABS(CC(J,I)) > MAXC(I)) MAXC(I)=DABS(CC(J,I))
    ENDDO
   ENDDO

   ! PRELIMINARY FOR THE OBARA-SAIKA RECURSION FORMULA
   DOWN1=0
   DOWN2=0
   DOWN3=0
   DOWN4=0
   DOWN5=0
   DO IA=0,BASIS_LMAX
    DO JA=0,BASIS_LMAX
     KLT=IA*(BASIS_LMAX+1)+JA+1
     NKL(KLT)=0
     DO I1=0,IA
      DO I2=0,IA-I1
       DO I3=0,IA-I1-I2
        DO J1=0,JA
         DO J2=0,JA-J1
          DO J3=0,JA-J1-J2
           NKL(KLT)=NKL(KLT)+1
           KLANGINDX(NKL(KLT),KLT)=I1+I2+I3+J1+J2+J3
           IJMAP(I1,I2,I3,J1,J2,J3)=NKL(KLT)
           IF ((I1 >= 1).AND.(J1 >= 1)) THEN
            DOWN1(NKL(KLT),KLT,1)=IJMAP(I1-1,I2,I3,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,1)=DFLOAT(I1)/2.0D0
            DOWN2(NKL(KLT),KLT,1)=IJMAP(I1,I2,I3,J1-1,J2,J3)
            DOWN2_C(NKL(KLT),KLT,1)=DFLOAT(J1)/2.0D0
           ELSE IF ((I1 >= 1).AND.(J1 == 0)) THEN
            DOWN1(NKL(KLT),KLT,1)=IJMAP(I1-1,I2,I3,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,1)=DFLOAT(I1)/2.0D0
           ELSE IF ((J1 >= 1).AND.(I1 == 0)) THEN
            DOWN2(NKL(KLT),KLT,1)=IJMAP(I1,I2,I3,J1-1,J2,J3)
            DOWN2_C(NKL(KLT),KLT,1)=DFLOAT(J1)/2.0D0
           ENDIF
           IF ((I2 >= 1).AND.(J2 >= 1)) THEN
            DOWN1(NKL(KLT),KLT,2)=IJMAP(I1,I2-1,I3,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,2)=DFLOAT(I2)/2.0D0
            DOWN2(NKL(KLT),KLT,2)=IJMAP(I1,I2,I3,J1,J2-1,J3)
            DOWN2_C(NKL(KLT),KLT,2)=DFLOAT(J2)/2.0D0
           ELSE IF ((I2 >= 1).AND.(J2 == 0)) THEN
            DOWN1(NKL(KLT),KLT,2)=IJMAP(I1,I2-1,I3,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,2)=DFLOAT(I2)/2.0D0
            ELSE IF ((J2 >= 1).AND.(I2 == 0)) THEN
            DOWN2(NKL(KLT),KLT,2)=IJMAP(I1,I2,I3,J1,J2-1,J3)
            DOWN2_C(NKL(KLT),KLT,2)=DFLOAT(J2)/2.0D0
           ENDIF
           IF ((I3 >= 1).AND.(J3 >= 1)) THEN
            DOWN1(NKL(KLT),KLT,3)=IJMAP(I1,I2,I3-1,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,3)=DFLOAT(I3)/2.0D0
            DOWN2(NKL(KLT),KLT,3)=IJMAP(I1,I2,I3,J1,J2,J3-1)
            DOWN2_C(NKL(KLT),KLT,3)=DFLOAT(J3)/2.0D0
           ELSE IF ((I3 >= 1).AND.(J3 == 0)) THEN
            DOWN1(NKL(KLT),KLT,3)=IJMAP(I1,I2,I3-1,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,3)=DFLOAT(I3)/2.0D0
           ELSE IF ((J3 >= 1).AND.(I3 == 0)) THEN
            DOWN2(NKL(KLT),KLT,3)=IJMAP(I1,I2,I3,J1,J2,J3-1)
            DOWN2_C(NKL(KLT),KLT,3)=DFLOAT(J3)/2.0D0
           ENDIF

           IF (I1 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=1
            DOWN3(NKL(KLT),KLT)=IJMAP(I1-1,I2,I3,J1,J2,J3)
            IF (I1 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1-2,I2,I3,J1,J2,J3)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(I1-1)/2.0D0
            ENDIF
            IF (J1 >= 1) THEN
             DOWN5(NKL(KLT),KLT)=IJMAP(I1-1,I2,I3,J1-1,J2,J3)
             DOWN5_C(NKL(KLT),KLT)=DFLOAT(J1)/2.0D0
            ENDIF
           ELSE IF (I2 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=2
            DOWN3(NKL(KLT),KLT)=IJMAP(I1,I2-1,I3,J1,J2,J3)
            IF (I2 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1,I2-2,I3,J1,J2,J3)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(I2-1)/2.0D0
            ENDIF
            IF (J2 >= 1) THEN
             DOWN5(NKL(KLT),KLT)=IJMAP(I1,I2-1,I3,J1,J2-1,J3)
             DOWN5_C(NKL(KLT),KLT)=DFLOAT(J2)/2.0D0
            ENDIF
           ELSE IF (I3 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=3
            DOWN3(NKL(KLT),KLT)=IJMAP(I1,I2,I3-1,J1,J2,J3)
            IF (I3 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1,I2,I3-2,J1,J2,J3)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(I3-1)/2.0D0
            ENDIF
            IF (J3 >= 1) THEN
             DOWN5(NKL(KLT),KLT)=IJMAP(I1,I2,I3-1,J1,J2,J3-1)
             DOWN5_C(NKL(KLT),KLT)=DFLOAT(J3)/2.0D0
            ENDIF
           ELSE IF (J1 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=1
            DOWN3(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1-1,J2,J3)
            IF (J1 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1-2,J2,J3)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(J1-1)/2.0D0
            ENDIF
           ELSE IF (J2 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=2
            DOWN3(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1,J2-1,J3)
            IF (J2 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1,J2-2,J3)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(J2-1)/2.0D0
            ENDIF
           ELSE IF (J3 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=3
            DOWN3(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1,J2,J3-1)
            IF (J3 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1,J2,J3-2)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(J3-1)/2.0D0
            ENDIF
           ELSE
            DOWNXYZ(NKL(KLT),KLT)=0
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
     DO I1=0,IA+1
      DO I2=0,IA+1-I1
       DO I3=IA+1-I1-I2,IA+1-I1-I2
        DO J1=0,JA
         DO J2=0,JA-J1
          DO J3=0,JA-J1-J2
           NKL(KLT)=NKL(KLT)+1
           KLANGINDX(NKL(KLT),KLT)=I1+I2+I3+J1+J2+J3
           IJMAP(I1,I2,I3,J1,J2,J3)=NKL(KLT)
           IF ((I1 >= 1).AND.(J1 >= 1)) THEN
            DOWN1(NKL(KLT),KLT,1)=IJMAP(I1-1,I2,I3,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,1)=DFLOAT(I1)/2.0D0
            DOWN2(NKL(KLT),KLT,1)=IJMAP(I1,I2,I3,J1-1,J2,J3)
            DOWN2_C(NKL(KLT),KLT,1)=DFLOAT(J1)/2.0D0
           ELSE IF ((I1 >= 1).AND.(J1 == 0)) THEN
            DOWN1(NKL(KLT),KLT,1)=IJMAP(I1-1,I2,I3,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,1)=DFLOAT(I1)/2.0D0
           ELSE IF ((J1 >= 1).AND.(I1 == 0)) THEN
            DOWN2(NKL(KLT),KLT,1)=IJMAP(I1,I2,I3,J1-1,J2,J3)
            DOWN2_C(NKL(KLT),KLT,1)=DFLOAT(J1)/2.0D0
           ENDIF
           IF ((I2 >= 1).AND.(J2 >= 1)) THEN
            DOWN1(NKL(KLT),KLT,2)=IJMAP(I1,I2-1,I3,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,2)=DFLOAT(I2)/2.0D0
            DOWN2(NKL(KLT),KLT,2)=IJMAP(I1,I2,I3,J1,J2-1,J3)
            DOWN2_C(NKL(KLT),KLT,2)=DFLOAT(J2)/2.0D0
           ELSE IF ((I2 >= 1).AND.(J2 == 0)) THEN
            DOWN1(NKL(KLT),KLT,2)=IJMAP(I1,I2-1,I3,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,2)=DFLOAT(I2)/2.0D0
           ELSE IF ((J2 >= 1).AND.(I2 == 0)) THEN
            DOWN2(NKL(KLT),KLT,2)=IJMAP(I1,I2,I3,J1,J2-1,J3)
            DOWN2_C(NKL(KLT),KLT,2)=DFLOAT(J2)/2.0D0
           ENDIF
           IF ((I3 >= 1).AND.(J3 >= 1)) THEN
            DOWN1(NKL(KLT),KLT,3)=IJMAP(I1,I2,I3-1,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,3)=DFLOAT(I3)/2.0D0
            DOWN2(NKL(KLT),KLT,3)=IJMAP(I1,I2,I3,J1,J2,J3-1)
            DOWN2_C(NKL(KLT),KLT,3)=DFLOAT(J3)/2.0D0
           ELSE IF ((I3 >= 1).AND.(J3 == 0)) THEN
            DOWN1(NKL(KLT),KLT,3)=IJMAP(I1,I2,I3-1,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,3)=DFLOAT(I3)/2.0D0
           ELSE IF ((J3 >= 1).AND.(I3 == 0)) THEN
            DOWN2(NKL(KLT),KLT,3)=IJMAP(I1,I2,I3,J1,J2,J3-1)
            DOWN2_C(NKL(KLT),KLT,3)=DFLOAT(J3)/2.0D0
           ENDIF

           IF (I1 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=1
            DOWN3(NKL(KLT),KLT)=IJMAP(I1-1,I2,I3,J1,J2,J3)
            IF (I1 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1-2,I2,I3,J1,J2,J3)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(I1-1)/2.0D0
            ENDIF
            IF (J1 >= 1) THEN
             DOWN5(NKL(KLT),KLT)=IJMAP(I1-1,I2,I3,J1-1,J2,J3)
             DOWN5_C(NKL(KLT),KLT)=DFLOAT(J1)/2.0D0
            ENDIF
           ELSE IF (I2 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=2
            DOWN3(NKL(KLT),KLT)=IJMAP(I1,I2-1,I3,J1,J2,J3)
            IF (I2 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1,I2-2,I3,J1,J2,J3)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(I2-1)/2.0D0
            ENDIF
            IF (J2 >= 1) THEN
             DOWN5(NKL(KLT),KLT)=IJMAP(I1,I2-1,I3,J1,J2-1,J3)
             DOWN5_C(NKL(KLT),KLT)=DFLOAT(J2)/2.0D0
            ENDIF
           ELSE IF (I3 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=3
            DOWN3(NKL(KLT),KLT)=IJMAP(I1,I2,I3-1,J1,J2,J3)
            IF (I3 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1,I2,I3-2,J1,J2,J3)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(I3-1)/2.0D0
            ENDIF
            IF (J3 >= 1) THEN
             DOWN5(NKL(KLT),KLT)=IJMAP(I1,I2,I3-1,J1,J2,J3-1)
             DOWN5_C(NKL(KLT),KLT)=DFLOAT(J3)/2.0D0
            ENDIF
           ELSE IF (J1 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=1
            DOWN3(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1-1,J2,J3)
            IF (J1 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1-2,J2,J3)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(J1-1)/2.0D0
            ENDIF
           ELSE IF (J2 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=2
            DOWN3(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1,J2-1,J3)
            IF (J2 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1,J2-2,J3)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(J2-1)/2.0D0
            ENDIF
           ELSE IF (J3 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=3
            DOWN3(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1,J2,J3-1)
            IF (J3 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1,J2,J3-2)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(J3-1)/2.0D0
            ENDIF
           ELSE
            DOWNXYZ(NKL(KLT),KLT)=0
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
     DO I1=0,IA
      DO I2=0,IA-I1
       DO I3=0,IA-I1-I2
        DO J1=0,JA+1
         DO J2=0,JA+1-J1
          DO J3=JA+1-J1-J2,JA+1-J1-J2
           NKL(KLT)=NKL(KLT)+1
           KLANGINDX(NKL(KLT),KLT)=I1+I2+I3+J1+J2+J3
           IJMAP(I1,I2,I3,J1,J2,J3)=NKL(KLT)
           IF ((I1 >= 1).AND.(J1 >= 1)) THEN
            DOWN1(NKL(KLT),KLT,1)=IJMAP(I1-1,I2,I3,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,1)=DFLOAT(I1)/2.0D0
            DOWN2(NKL(KLT),KLT,1)=IJMAP(I1,I2,I3,J1-1,J2,J3)
            DOWN2_C(NKL(KLT),KLT,1)=DFLOAT(J1)/2.0D0
           ELSE IF ((I1 >= 1).AND.(J1 == 0)) THEN
            DOWN1(NKL(KLT),KLT,1)=IJMAP(I1-1,I2,I3,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,1)=DFLOAT(I1)/2.0D0
           ELSE IF ((J1 >= 1).AND.(I1 == 0)) THEN
            DOWN2(NKL(KLT),KLT,1)=IJMAP(I1,I2,I3,J1-1,J2,J3)
            DOWN2_C(NKL(KLT),KLT,1)=DFLOAT(J1)/2.0D0
           ENDIF
           IF ((I2 >= 1).AND.(J2 >= 1)) THEN
            DOWN1(NKL(KLT),KLT,2)=IJMAP(I1,I2-1,I3,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,2)=DFLOAT(I2)/2.0D0
            DOWN2(NKL(KLT),KLT,2)=IJMAP(I1,I2,I3,J1,J2-1,J3)
            DOWN2_C(NKL(KLT),KLT,2)=DFLOAT(J2)/2.0D0
           ELSE IF ((I2 >= 1).AND.(J2 == 0)) THEN
            DOWN1(NKL(KLT),KLT,2)=IJMAP(I1,I2-1,I3,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,2)=DFLOAT(I2)/2.0D0
           ELSE IF ((J2 >= 1).AND.(I2 == 0)) THEN
            DOWN2(NKL(KLT),KLT,2)=IJMAP(I1,I2,I3,J1,J2-1,J3)
            DOWN2_C(NKL(KLT),KLT,2)=DFLOAT(J2)/2.0D0
           ENDIF
           IF ((I3 >= 1).AND.(J3 >= 1)) THEN
            DOWN1(NKL(KLT),KLT,3)=IJMAP(I1,I2,I3-1,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,3)=DFLOAT(I3)/2.0D0
            DOWN2(NKL(KLT),KLT,3)=IJMAP(I1,I2,I3,J1,J2,J3-1)
            DOWN2_C(NKL(KLT),KLT,3)=DFLOAT(J3)/2.0D0
           ELSE IF ((I3 >= 1).AND.(J3 == 0)) THEN
            DOWN1(NKL(KLT),KLT,3)=IJMAP(I1,I2,I3-1,J1,J2,J3)
            DOWN1_C(NKL(KLT),KLT,3)=DFLOAT(I3)/2.0D0
           ELSE IF ((J3 >= 1).AND.(I3 == 0)) THEN
            DOWN2(NKL(KLT),KLT,3)=IJMAP(I1,I2,I3,J1,J2,J3-1)
            DOWN2_C(NKL(KLT),KLT,3)=DFLOAT(J3)/2.0D0
           ENDIF

           IF (I1 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=1
            DOWN3(NKL(KLT),KLT)=IJMAP(I1-1,I2,I3,J1,J2,J3)
            IF (I1 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1-2,I2,I3,J1,J2,J3)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(I1-1)/2.0D0
            ENDIF
            IF (J1 >= 1) THEN
             DOWN5(NKL(KLT),KLT)=IJMAP(I1-1,I2,I3,J1-1,J2,J3)
             DOWN5_C(NKL(KLT),KLT)=DFLOAT(J1)/2.0D0
            ENDIF
           ELSE IF (I2 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=2
            DOWN3(NKL(KLT),KLT)=IJMAP(I1,I2-1,I3,J1,J2,J3)
            IF (I2 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1,I2-2,I3,J1,J2,J3)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(I2-1)/2.0D0
            ENDIF
            IF (J2 >= 1) THEN
             DOWN5(NKL(KLT),KLT)=IJMAP(I1,I2-1,I3,J1,J2-1,J3)
             DOWN5_C(NKL(KLT),KLT)=DFLOAT(J2)/2.0D0
            ENDIF
           ELSE IF (I3 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=3
            DOWN3(NKL(KLT),KLT)=IJMAP(I1,I2,I3-1,J1,J2,J3)
            IF (I3 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1,I2,I3-2,J1,J2,J3)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(I3-1)/2.0D0
            ENDIF
            IF (J3 >= 1) THEN
             DOWN5(NKL(KLT),KLT)=IJMAP(I1,I2,I3-1,J1,J2,J3-1)
             DOWN5_C(NKL(KLT),KLT)=DFLOAT(J3)/2.0D0
            ENDIF
           ELSE IF (J1 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=1
            DOWN3(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1-1,J2,J3)
            IF (J1 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1-2,J2,J3)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(J1-1)/2.0D0
            ENDIF
           ELSE IF (J2 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=2
            DOWN3(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1,J2-1,J3)
            IF (J2 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1,J2-2,J3)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(J2-1)/2.0D0
            ENDIF
           ELSE IF (J3 >= 1) THEN
            DOWNXYZ(NKL(KLT),KLT)=3
            DOWN3(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1,J2,J3-1)
            IF (J3 >= 2) THEN
             DOWN4(NKL(KLT),KLT)=IJMAP(I1,I2,I3,J1,J2,J3-2)
             DOWN4_C(NKL(KLT),KLT)=DFLOAT(J3-1)/2.0D0
            ENDIF
           ELSE
            DOWNXYZ(NKL(KLT),KLT)=0
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
     NKL_G(KLT)=0
     DO I1=0,IA
      DO I2=0,IA-I1
       DO I3=0,IA-I1-I2
        DO J1=0,JA
         DO J2=0,JA-J1
          DO J3=0,JA-J1-J2
           NKL_G(KLT)=NKL_G(KLT)+1
           UP1(NKL_G(KLT),KLT,1)=IJMAP(I1+1,I2,I3,J1,J2,J3)
           UP1(NKL_G(KLT),KLT,2)=IJMAP(I1,I2+1,I3,J1,J2,J3)
           UP1(NKL_G(KLT),KLT,3)=IJMAP(I1,I2,I3+1,J1,J2,J3)
           UP2(NKL_G(KLT),KLT,1)=IJMAP(I1,I2,I3,J1+1,J2,J3)
           UP2(NKL_G(KLT),KLT,2)=IJMAP(I1,I2,I3,J1,J2+1,J3)
           UP2(NKL_G(KLT),KLT,3)=IJMAP(I1,I2,I3,J1,J2,J3+1)
           IF ((IJMAP(I1+1,I2,I3,J1,J2,J3) == 0).OR.(IJMAP(I1,I2+1,I3,J1,J2,J3) == 0) &
           .OR.(IJMAP(I1,I2,I3+1,J1,J2,J3) == 0).OR.(IJMAP(I1,I2,I3,J1+1,J2,J3) == 0) &
           .OR.(IJMAP(I1,I2,I3,J1,J2+1,J3) == 0).OR.(IJMAP(I1,I2,I3,J1,J2,J3+1) == 0)) &
           CALL PABORT('AN INTERNAL PROGRAM ERROR IS DETECTED')
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
   
   ! LOOP OVER UNIT CELL INDEX 1
   CALL PCPU_TIME(ICPUS)
   TICKET=0
   DO Q1X=-REDUCED_CEL1X,REDUCED_CEL1X
   DO Q1Y=-REDUCED_CEL1Y,REDUCED_CEL1Y
   DO Q1Z=-REDUCED_CEL1Z,REDUCED_CEL1Z
    P1X=DFLOAT(Q1X)*PERIODX
    P1Y=DFLOAT(Q1Y)*PERIODY
    P1Z=DFLOAT(Q1Z)*PERIODZ
    P1C=P3C(Q1X)
    P1S=P3S(Q1X)
    ! LOOP OVER UNIT CELL INDEX 2
    DO Q2X=-REDUCED_CEL1X,REDUCED_CEL1X
    DO Q2Y=-REDUCED_CEL1Y,REDUCED_CEL1Y
    DO Q2Z=-REDUCED_CEL1Z,REDUCED_CEL1Z
     TICKET=TICKET+1
     IF (TICKET == NUMPROCS) TICKET=TICKET-NUMPROCS
!--- PARALLEL LOOP
     IF (MYID == TICKET) THEN
!--- PARALLEL LOOP
     P2X=DFLOAT(Q2X)*PERIODX
     P2Y=DFLOAT(Q2Y)*PERIODY
     P2Z=DFLOAT(Q2Z)*PERIODZ
     P2C=P3C(Q2X)
     P2S=P3S(Q2X)
     ! SET CELL SYMMETRY INDECES
     LCEL1X=.FALSE.
     LCEL1Y=.FALSE.
     LCEL1Z=.FALSE.
     DO Q3X=-CEL2X,CEL2X
      IF ((Q3X >= -CEL1X).AND.(Q3X <= CEL1X).AND.(Q2X+Q3X-Q1X >= -CEL1X).AND.(Q2X+Q3X-Q1X <= CEL1X)) LCEL1X(Q3X)=.TRUE.
     ENDDO
     DO Q3Y=-CEL2Y,CEL2Y
      IF ((Q3Y >= -CEL1Y).AND.(Q3Y <= CEL1Y).AND.(Q2Y+Q3Y-Q1Y >= -CEL1Y).AND.(Q2Y+Q3Y-Q1Y <= CEL1Y)) LCEL1Y(Q3Y)=.TRUE.
     ENDDO
     DO Q3Z=-CEL2Z,CEL2Z
      IF ((Q3Z >= -CEL1Z).AND.(Q3Z <= CEL1Z).AND.(Q2Z+Q3Z-Q1Z >= -CEL1Z).AND.(Q2Z+Q3Z-Q1Z <= CEL1Z)) LCEL1Z(Q3Z)=.TRUE.
     ENDDO
     ! LOOP OVER KL SHELL
     DO K=1,NPSHELL
      KK=P_SHELL(K,0,0,0)
      CX=PGSX(KK)
      CY=PGSY(KK)
      CZ=PGSZ(KK)
      DO L=1,NPSHELL
       LL=P_SHELL(L,0,0,0)
       DX=PGSX(LL)+P2X
       DY=PGSY(LL)*P2C-PGSZ(LL)*P2S+P2Y
       DZ=PGSY(LL)*P2S+PGSZ(LL)*P2C+P2Z
       COMZ=ZT(KK)+ZT(LL)
       PX=(CX*ZT(KK)+DX*ZT(LL))/COMZ
       PY=(CY*ZT(KK)+DY*ZT(LL))/COMZ
       PZ=(CZ*ZT(KK)+DZ*ZT(LL))/COMZ
       SKL=(K-1)*NPSHELL+L   ! SHELL INDEX
       COMZ_S(SKL)=COMZ      ! STORE ZETA+ZETA
       ZTK(SKL)=ZT(KK)       ! STORE ZETA
       ZTL(SKL)=ZT(LL)       ! STORE ZETA
       PX_S(SKL)=PX          ! STORE PX, PY, & PZ
       PY_S(SKL)=PY
       PZ_S(SKL)=PZ
       OV(SKL)=S_P(KK,LL,Q2X,Q2Y,Q2Z)   ! STORE UNNORMALIZED OVERLAP INTEGRALS
       KLANG(SKL)=P_SHELL_ANG(K)+P_SHELL_ANG(L)     ! STORE SUM OF ANGULAR MOMENTA
       KLT_S(SKL)=P_SHELL_ANG(K)*(BASIS_LMAX+1)+P_SHELL_ANG(L)+1 ! STORE SHELL TYPE
       KATM(SKL)=PGS_ATOM(KK)
       LATM(SKL)=PGS_ATOM(LL)
       NKL2=0
       ! ELECTRON 2
       DO I1=0,P_SHELL_ANG(K)
        DO I2=0,P_SHELL_ANG(K)-I1
         DO I3=0,P_SHELL_ANG(K)-I1-I2
          DO J1=0,P_SHELL_ANG(L)
           DO J2=0,P_SHELL_ANG(L)-J1
            DO J3=0,P_SHELL_ANG(L)-J1-J2
             NKL2=NKL2+1
             PGSK(NKL2,SKL)=P_SHELL(K,I1,I2,I3)
             PGSL(NKL2,SKL)=P_SHELL(L,J1,J2,J3)
             IF (I1 >= 1) THEN
              QI_CI(NKL2,SKL)=PX-CX
             ELSE IF (I2 >= 1) THEN
              QI_CI(NKL2,SKL)=PY-CY
             ELSE IF (I3 >= 1) THEN
              QI_CI(NKL2,SKL)=PZ-CZ
             ELSE IF (J1 >= 1) THEN
              QI_CI(NKL2,SKL)=PX-DX
             ELSE IF (J2 >= 1) THEN
              QI_CI(NKL2,SKL)=PY-DY
             ELSE IF (J3 >= 1) THEN
              QI_CI(NKL2,SKL)=PZ-DZ
             ENDIF
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       DO I1=0,P_SHELL_ANG(K)+1
        DO I2=0,P_SHELL_ANG(K)+1-I1
         DO I3=P_SHELL_ANG(K)+1-I1-I2,P_SHELL_ANG(K)+1-I1-I2
          DO J1=0,P_SHELL_ANG(L)
           DO J2=0,P_SHELL_ANG(L)-J1
            DO J3=0,P_SHELL_ANG(L)-J1-J2
             NKL2=NKL2+1
             IF (I1 >= 1) THEN
              QI_CI(NKL2,SKL)=PX-CX
             ELSE IF (I2 >= 1) THEN
              QI_CI(NKL2,SKL)=PY-CY
             ELSE IF (I3 >= 1) THEN
              QI_CI(NKL2,SKL)=PZ-CZ
             ELSE IF (J1 >= 1) THEN
              QI_CI(NKL2,SKL)=PX-DX
             ELSE IF (J2 >= 1) THEN
              QI_CI(NKL2,SKL)=PY-DY
             ELSE IF (J3 >= 1) THEN
              QI_CI(NKL2,SKL)=PZ-DZ
             ENDIF
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       DO I1=0,P_SHELL_ANG(K)
        DO I2=0,P_SHELL_ANG(K)-I1
         DO I3=0,P_SHELL_ANG(K)-I1-I2
          DO J1=0,P_SHELL_ANG(L)+1
           DO J2=0,P_SHELL_ANG(L)+1-J1
            DO J3=P_SHELL_ANG(L)+1-J1-J2,P_SHELL_ANG(L)+1-J1-J2
             NKL2=NKL2+1
             IF (I1 >= 1) THEN
              QI_CI(NKL2,SKL)=PX-CX
             ELSE IF (I2 >= 1) THEN
              QI_CI(NKL2,SKL)=PY-CY
             ELSE IF (I3 >= 1) THEN
              QI_CI(NKL2,SKL)=PZ-CZ
             ELSE IF (J1 >= 1) THEN
              QI_CI(NKL2,SKL)=PX-DX
             ELSE IF (J2 >= 1) THEN
              QI_CI(NKL2,SKL)=PY-DY
             ELSE IF (J3 >= 1) THEN
              QI_CI(NKL2,SKL)=PZ-DZ
             ENDIF
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       IF (NKL2 /= NKL(KLT_S(SKL))) CALL PABORT('AN INTERNAL PROGRAM ERROR IS DETECTED')
      ENDDO
     ENDDO

     ! FIND THE MAXIMUM (KL|KL) VALUES FOR ALL SHELLS
     KLKL=0.0D0
     DO SKL=1,NPSHELL_SQ
      KLT=KLT_S(SKL)
      LARGECOMZ=2.0D0*COMZ_S(SKL)
      RHO=0.5D0*COMZ_S(SKL)
      IJKLANG=2*KLANG(SKL)
      C2=2.0D0*DSQRT(RHO/PI)*OV(SKL)**2
      IF (IJKLANG > 0) THEN
       DO M=0,IJKLANG
        ERI(M,1,1)=C2/DFLOAT(2*M+1)
       ENDDO
       DO IJ=1,NKL_G(KLT)
        IF (IJ == 1) THEN
         DO KL=2,NKL_G(KLT)
          DO M=0,IJKLANG-KLANGINDX(KL,KLT)
           ERI(M,KL,1)=QI_CI(KL,SKL)*ERI(M,DOWN3(KL,KLT),1)
           IF (DOWN4(KL,KLT) /= 0) ERI(M,KL,1)=ERI(M,KL,1) &
           +DOWN4_C(KL,KLT)/COMZ_S(SKL)*(ERI(M,DOWN4(KL,KLT),1)-0.5D0*ERI(M+1,DOWN4(KL,KLT),1))
           IF (DOWN5(KL,KLT) /= 0) ERI(M,KL,1)=ERI(M,KL,1) &
           +DOWN5_C(KL,KLT)/COMZ_S(SKL)*(ERI(M,DOWN5(KL,KLT),1)-0.5D0*ERI(M+1,DOWN5(KL,KLT),1))
          ENDDO
         ENDDO
        ELSE
         DO KL=1,NKL_G(KLT)
          DO M=0,IJKLANG-KLANGINDX(IJ,KLT)-KLANGINDX(KL,KLT)
           ERI(M,KL,IJ)=QI_CI(IJ,SKL)*ERI(M,KL,DOWN3(IJ,KLT))
           IF (DOWN4(IJ,KLT) /= 0) ERI(M,KL,IJ)=ERI(M,KL,IJ) &
           +DOWN4_C(IJ,KLT)/COMZ_S(SKL)*(ERI(M,KL,DOWN4(IJ,KLT))-0.5D0*ERI(M+1,KL,DOWN4(IJ,KLT)))
           IF (DOWN5(IJ,KLT) /= 0) ERI(M,KL,IJ)=ERI(M,KL,IJ) &
           +DOWN5_C(IJ,KLT)/COMZ_S(SKL)*(ERI(M,KL,DOWN5(IJ,KLT))-0.5D0*ERI(M+1,KL,DOWN5(IJ,KLT)))
           IF (DOWN1(KL,KLT,DOWNXYZ(IJ,KLT)) /= 0) ERI(M,KL,IJ)=ERI(M,KL,IJ)+DOWN1_C(KL,KLT,DOWNXYZ(IJ,KLT))/LARGECOMZ &
                                                   *ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,KLT)),DOWN3(IJ,KLT))
           IF (DOWN2(KL,KLT,DOWNXYZ(IJ,KLT)) /= 0) ERI(M,KL,IJ)=ERI(M,KL,IJ)+DOWN2_C(KL,KLT,DOWNXYZ(IJ,KLT))/LARGECOMZ &
                                                   *ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,KLT)),DOWN3(IJ,KLT))
          ENDDO
         ENDDO
        ENDIF
       ENDDO
       DO KL=2,NKL_G(KLT)
        A=DSQRT(ERI(0,KL,KL)*MAXC(PGSK(KL,SKL))*MAXC(PGSL(KL,SKL)))
        IF (A > KLKL(SKL)) KLKL(SKL)=A
       ENDDO
      ELSE
       A=DSQRT(C2*MAXC(PGSK(1,SKL))*MAXC(PGSL(1,SKL)))
       IF (A > KLKL(SKL)) KLKL(SKL)=A
      ENDIF
     ENDDO

     ! MAIN LOOP STARTS
     DO I=1,NPSHELL
      II=P_SHELL(I,0,0,0)
      IATM=PGS_ATOM(II)
      AX=PGSX(II)
      AY=PGSY(II)
      AZ=PGSZ(II)
      DO J=1,NPSHELL
       JJ=P_SHELL(J,0,0,0)
       JATM=PGS_ATOM(JJ)
       BX=PGSX(JJ)+P1X
       BY=PGSY(JJ)*P1C-PGSZ(JJ)*P1S+P1Y
       BZ=PGSY(JJ)*P1S+PGSZ(JJ)*P1C+P1Z
       COMZ=ZT(II)+ZT(JJ)
       PX=(AX*ZT(II)+BX*ZT(JJ))/COMZ
       PY=(AY*ZT(II)+BY*ZT(JJ))/COMZ
       PZ=(AZ*ZT(II)+BZ*ZT(JJ))/COMZ
       C1=S_P(II,JJ,Q1X,Q1Y,Q1Z)*2.0D0/DSQRT(PI)
       SIJ=(I-1)*NPSHELL+J   ! SHELL INDEX
       IJT=P_SHELL_ANG(I)*(BASIS_LMAX+1)+P_SHELL_ANG(J)+1
       IJANG=P_SHELL_ANG(I)+P_SHELL_ANG(J)
       
       ! ELECTRON 1
       NIJ=0
       NIJ_G=0
       DO I1=0,P_SHELL_ANG(I)
        DO I2=0,P_SHELL_ANG(I)-I1
         DO I3=0,P_SHELL_ANG(I)-I1-I2
          DO J1=0,P_SHELL_ANG(J)
           DO J2=0,P_SHELL_ANG(J)-J1
            DO J3=0,P_SHELL_ANG(J)-J1-J2
             NIJ=NIJ+1
             NIJ_G=NIJ_G+1
             PGSI(NIJ)=P_SHELL(I,I1,I2,I3)
             PGSJ(NIJ)=P_SHELL(J,J1,J2,J3)
             IF (I1 >= 1) THEN
              PI_AI(NIJ)=PX-AX
             ELSE IF (I2 >= 1) THEN
              PI_AI(NIJ)=PY-AY
             ELSE IF (I3 >= 1) THEN
              PI_AI(NIJ)=PZ-AZ
             ELSE IF (J1 >= 1) THEN
              PI_AI(NIJ)=PX-BX
             ELSE IF (J2 >= 1) THEN
              PI_AI(NIJ)=PY-BY
             ELSE IF (J3 >= 1) THEN
              PI_AI(NIJ)=PZ-BZ
             ENDIF
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       DO I1=0,P_SHELL_ANG(I)+1
        DO I2=0,P_SHELL_ANG(I)+1-I1
         DO I3=P_SHELL_ANG(I)+1-I1-I2,P_SHELL_ANG(I)+1-I1-I2
          DO J1=0,P_SHELL_ANG(J)
           DO J2=0,P_SHELL_ANG(J)-J1
            DO J3=0,P_SHELL_ANG(J)-J1-J2
             NIJ=NIJ+1
             IF (I1 >= 1) THEN
              PI_AI(NIJ)=PX-AX
             ELSE IF (I2 >= 1) THEN
              PI_AI(NIJ)=PY-AY
             ELSE IF (I3 >= 1) THEN
              PI_AI(NIJ)=PZ-AZ
             ELSE IF (J1 >= 1) THEN
              PI_AI(NIJ)=PX-BX
             ELSE IF (J2 >= 1) THEN
              PI_AI(NIJ)=PY-BY
             ELSE IF (J3 >= 1) THEN
              PI_AI(NIJ)=PZ-BZ
             ENDIF
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       IF (NIJ == 1) CALL PABORT('AN INTERNAL PROGRAM ERROR IS DETECTED')

       ! FIND THE MAXIMUM (IJ|IJ) VALUE FOR THE CURRENT SHELL
       IJIJ=0.0D0
       LARGECOMZ=2.0D0*COMZ
       RHO=0.5D0*COMZ
       IJKLANG=2*IJANG
       C2=2.0D0*DSQRT(RHO/PI)*S_P(II,JJ,Q1X,Q1Y,Q1Z)**2
       IF (IJKLANG > 0) THEN
        DO M=0,IJKLANG
         ERI(M,1,1)=C2/DFLOAT(2*M+1)
        ENDDO
        DO IJ=1,NIJ_G
         IF (IJ == 1) THEN
          DO KL=2,NIJ_G
           DO M=0,IJKLANG-KLANGINDX(KL,IJT)
            ERI(M,KL,1)=QI_CI(KL,SIJ)*ERI(M,DOWN3(KL,IJT),1)
            IF (DOWN4(KL,IJT) /= 0) ERI(M,KL,1)=ERI(M,KL,1) &
           +DOWN4_C(KL,IJT)/COMZ*(ERI(M,DOWN4(KL,IJT),1)-0.5D0*ERI(M+1,DOWN4(KL,IJT),1))
            IF (DOWN5(KL,IJT) /= 0) ERI(M,KL,1)=ERI(M,KL,1) &
           +DOWN5_C(KL,IJT)/COMZ*(ERI(M,DOWN5(KL,IJT),1)-0.5D0*ERI(M+1,DOWN5(KL,IJT),1))
           ENDDO
          ENDDO
         ELSE
          DO KL=1,NIJ_G
           DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,IJT)
            ERI(M,KL,IJ)=QI_CI(IJ,SIJ)*ERI(M,KL,DOWN3(IJ,IJT))
            IF (DOWN4(IJ,IJT) /= 0) ERI(M,KL,IJ)=ERI(M,KL,IJ) &
           +DOWN4_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN4(IJ,IJT))-0.5D0*ERI(M+1,KL,DOWN4(IJ,IJT)))
            IF (DOWN5(IJ,IJT) /= 0) ERI(M,KL,IJ)=ERI(M,KL,IJ) &
           +DOWN5_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN5(IJ,IJT))-0.5D0*ERI(M+1,KL,DOWN5(IJ,IJT)))
            IF (DOWN1(KL,IJT,DOWNXYZ(IJ,IJT)) /= 0) ERI(M,KL,IJ)=ERI(M,KL,IJ)+DOWN1_C(KL,IJT,DOWNXYZ(IJ,IJT))/LARGECOMZ &
                                                    *ERI(M+1,DOWN1(KL,IJT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
            IF (DOWN2(KL,IJT,DOWNXYZ(IJ,IJT)) /= 0) ERI(M,KL,IJ)=ERI(M,KL,IJ)+DOWN2_C(KL,IJT,DOWNXYZ(IJ,IJT))/LARGECOMZ &
                                                    *ERI(M+1,DOWN2(KL,IJT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
           ENDDO
          ENDDO
         ENDIF
        ENDDO
        DO IJ=2,NIJ_G
         A=DSQRT(ERI(0,IJ,IJ)*MAXC(PGSI(IJ))*MAXC(PGSJ(IJ)))
         IF (A > IJIJ) IJIJ=A
        ENDDO
       ELSE
        A=DSQRT(C2*MAXC(PGSI(1))*MAXC(PGSJ(1)))
        IF (A > IJIJ) IJIJ=A
       ENDIF
       
       ! LOOP OVER UNIT CELL INDEX 3
       DO Q3X=-CEL2X,CEL2X
       DO Q3Y=-CEL2Y,CEL2Y
       DO Q3Z=-CEL2Z,CEL2Z
        ! LOOP OVER K SHELL
        DO K=1,NPSHELL
         KK=P_SHELL(K,0,0,0)
         CX=PGSX(KK)+P3X(Q3X)
         CY=PGSY(KK)*P3C(Q3X)-PGSZ(KK)*P3S(Q3X)+P3Y(Q3Y)
         CZ=PGSY(KK)*P3S(Q3X)+PGSZ(KK)*P3C(Q3X)+P3Z(Q3Z)
         ! LOOP OVER L SHELL
         DO L=1,NPSHELL
          LL=P_SHELL(L,0,0,0)
          DX=PGSX(LL)+P3X(Q3X+Q2X)
          DY=PGSY(LL)*P3C(Q3X+Q2X)-PGSZ(LL)*P3S(Q3X+Q2X)+P3Y(Q3Y+Q2Y)
          DZ=PGSY(LL)*P3S(Q3X+Q2X)+PGSZ(LL)*P3C(Q3X+Q2X)+P3Z(Q3Z+Q2Z)
          COMZ2=ZT(KK)+ZT(LL)
          QX=(CX*ZT(KK)+DX*ZT(LL))/COMZ2
          QY=(CY*ZT(KK)+DY*ZT(LL))/COMZ2
          QZ=(CZ*ZT(KK)+DZ*ZT(LL))/COMZ2
          SKL=(K-1)*NPSHELL+L   ! SHELL INDEX
          COMZ_S(SKL)=COMZ2     ! STORE ZETA+ZETA
          ZTK(SKL)=ZT(KK)       ! STORE ZETA
          ZTL(SKL)=ZT(LL)       ! STORE ZETA
          PX_S(SKL)=QX          ! STORE PX, PY, & PZ
          PY_S(SKL)=QY
          PZ_S(SKL)=QZ
          OV(SKL)=S_P(KK,LL,Q2X,Q2Y,Q2Z)   ! STORE UNNORMALIZED OVERLAP INTEGRALS
          KLANG(SKL)=P_SHELL_ANG(K)+P_SHELL_ANG(L)     ! STORE SUM OF ANGULAR MOMENTA
          KLT_S(SKL)=P_SHELL_ANG(K)*(BASIS_LMAX+1)+P_SHELL_ANG(L)+1 ! STORE SHELL TYPE
          KATM(SKL)=PGS_ATOM(KK)
          LATM(SKL)=PGS_ATOM(LL)

          NKL2=0
          ! ELECTRON 2
          DO I1=0,P_SHELL_ANG(K)
           DO I2=0,P_SHELL_ANG(K)-I1
            DO I3=0,P_SHELL_ANG(K)-I1-I2
             DO J1=0,P_SHELL_ANG(L)
              DO J2=0,P_SHELL_ANG(L)-J1
               DO J3=0,P_SHELL_ANG(L)-J1-J2
                NKL2=NKL2+1
                PGSK(NKL2,SKL)=P_SHELL(K,I1,I2,I3)
                PGSL(NKL2,SKL)=P_SHELL(L,J1,J2,J3)
                IF (I1 >= 1) THEN
                 QI_CI(NKL2,SKL)=QX-CX
                ELSE IF (I2 >= 1) THEN
                 QI_CI(NKL2,SKL)=QY-CY
                ELSE IF (I3 >= 1) THEN
                 QI_CI(NKL2,SKL)=QZ-CZ
                ELSE IF (J1 >= 1) THEN
                 QI_CI(NKL2,SKL)=QX-DX
                ELSE IF (J2 >= 1) THEN
                 QI_CI(NKL2,SKL)=QY-DY
                ELSE IF (J3 >= 1) THEN
                 QI_CI(NKL2,SKL)=QZ-DZ
                ENDIF
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
          DO I1=0,P_SHELL_ANG(K)+1
           DO I2=0,P_SHELL_ANG(K)+1-I1
            DO I3=P_SHELL_ANG(K)+1-I1-I2,P_SHELL_ANG(K)+1-I1-I2
             DO J1=0,P_SHELL_ANG(L)
              DO J2=0,P_SHELL_ANG(L)-J1
               DO J3=0,P_SHELL_ANG(L)-J1-J2
                NKL2=NKL2+1
                IF (I1 >= 1) THEN
                 QI_CI(NKL2,SKL)=QX-CX
                ELSE IF (I2 >= 1) THEN
                 QI_CI(NKL2,SKL)=QY-CY
                ELSE IF (I3 >= 1) THEN
                 QI_CI(NKL2,SKL)=QZ-CZ
                ELSE IF (J1 >= 1) THEN
                 QI_CI(NKL2,SKL)=QX-DX
                ELSE IF (J2 >= 1) THEN
                 QI_CI(NKL2,SKL)=QY-DY
                ELSE IF (J3 >= 1) THEN
                 QI_CI(NKL2,SKL)=QZ-DZ
                ENDIF
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
          DO I1=0,P_SHELL_ANG(K)
           DO I2=0,P_SHELL_ANG(K)-I1
            DO I3=0,P_SHELL_ANG(K)-I1-I2
             DO J1=0,P_SHELL_ANG(L)+1
              DO J2=0,P_SHELL_ANG(L)+1-J1
               DO J3=P_SHELL_ANG(L)+1-J1-J2,P_SHELL_ANG(L)+1-J1-J2
                NKL2=NKL2+1
                IF (I1 >= 1) THEN
                 QI_CI(NKL2,SKL)=QX-CX
                ELSE IF (I2 >= 1) THEN
                 QI_CI(NKL2,SKL)=QY-CY
                ELSE IF (I3 >= 1) THEN
                 QI_CI(NKL2,SKL)=QZ-CZ
                ELSE IF (J1 >= 1) THEN
                 QI_CI(NKL2,SKL)=QX-DX
                ELSE IF (J2 >= 1) THEN
                 QI_CI(NKL2,SKL)=QY-DY
                ELSE IF (J3 >= 1) THEN
                 QI_CI(NKL2,SKL)=QZ-DZ
                ENDIF
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
          IF (NKL2 /= NKL(KLT_S(SKL))) CALL PABORT('AN INTERNAL PROGRAM ERROR IS DETECTED')
         ENDDO
        ENDDO
        
        DO SKL=1,NPSHELL_SQ
         IF (IJIJ*KLKL(SKL) < DOPTN(12)) CYCLE
         LARGECOMZ=COMZ+COMZ_S(SKL)
         RHO_COMZ1=COMZ_S(SKL)/LARGECOMZ
         RHO_COMZ2=COMZ/LARGECOMZ
         RHO=COMZ*RHO_COMZ1
         IJKLANG=IJANG+KLANG(SKL)+1
         C2=C1*DSQRT(RHO)*OV(SKL)
         KLT=KLT_S(SKL)
         QX=PX_S(SKL)
         QY=PY_S(SKL)
         QZ=PZ_S(SKL)
         T=RHO*((QX-PX)**2+(QY-PY)**2+(QZ-PZ)**2)
         DO M=0,IJKLANG
          IF (T < TF(M)) THEN
           TS=NINT(T*20.0D0)
           H=0.05D0*DFLOAT(TS)-T
           ERI(M,1,1)=((((((IGAMMA(TS,M+6)*H*0.166666666666667D0+IGAMMA(TS,M+5))*H*0.2D0+IGAMMA(TS,M+4))*H*0.25D0+ &
           IGAMMA(TS,M+3))*H*0.333333333333333D0+IGAMMA(TS,M+2))*H*0.5D0+IGAMMA(TS,M+1))*H+IGAMMA(TS,M))*C2
          ELSE
           ERI(M,1,1)=C2*F1(M)*DEXP(F2(M)*DLOG(T))
          ENDIF
         ENDDO
         WX=(PX*COMZ+QX*COMZ_S(SKL))/LARGECOMZ
         WY=(PY*COMZ+QY*COMZ_S(SKL))/LARGECOMZ
         WZ=(PZ*COMZ+QZ*COMZ_S(SKL))/LARGECOMZ
         WI_PI(1)=WX-PX
         WI_PI(2)=WY-PY
         WI_PI(3)=WZ-PZ
         WI_QI(1)=WX-QX
         WI_QI(2)=WY-QY
         WI_QI(3)=WZ-QZ
         ! CORE INTEGRAL EVALUATION
         IF (NKL_G(KLT) /= 1) THEN
          DO KL=2,NKL_G(KLT)
           DO M=0,IJKLANG-KLANGINDX(KL,KLT)
            ERI(M,KL,1)=QI_CI(KL,SKL)*ERI(M,DOWN3(KL,KLT),1)+WI_QI(DOWNXYZ(KL,KLT))*ERI(M+1,DOWN3(KL,KLT),1)
            IF (DOWN4(KL,KLT) /= 0) ERI(M,KL,1)=ERI(M,KL,1) &
           +DOWN4_C(KL,KLT)/COMZ_S(SKL)*(ERI(M,DOWN4(KL,KLT),1)-RHO_COMZ2*ERI(M+1,DOWN4(KL,KLT),1))
            IF (DOWN5(KL,KLT) /= 0) ERI(M,KL,1)=ERI(M,KL,1) &
           +DOWN5_C(KL,KLT)/COMZ_S(SKL)*(ERI(M,DOWN5(KL,KLT),1)-RHO_COMZ2*ERI(M+1,DOWN5(KL,KLT),1))
           ENDDO
          ENDDO
         ENDIF
         IF (NIJ_G /= 1) THEN
          DO IJ=2,NIJ_G
           IF ((DOWN4(IJ,IJT) == 0).AND.(DOWN5(IJ,IJT) == 0)) THEN
            DO KL=1,NKL_G(KLT)
             IF ((DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0).AND.(DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0)) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT))
              ENDDO
             ELSE IF (DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
              ENDDO
             ELSE IF (DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
              ENDDO
             ELSE
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               +(DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)) &
               + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)))/LARGECOMZ
              ENDDO
             ENDIF
            ENDDO
           ELSE IF (DOWN5(IJ,IJT) == 0) THEN
            DO KL=1,NKL_G(KLT)
             IF ((DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0).AND.(DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0)) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN4_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT)))
              ENDDO
             ELSE IF (DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN4_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
               + DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
              ENDDO
             ELSE IF (DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN4_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
               + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
              ENDDO
             ELSE
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN4_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
               +(DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)) &
               + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)))/LARGECOMZ
              ENDDO
             ENDIF
            ENDDO
           ELSE IF (DOWN4(IJ,IJT) == 0) THEN
            DO KL=1,NKL_G(KLT)
             IF ((DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0).AND.(DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0)) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN5_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT)))
              ENDDO
             ELSE IF (DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN5_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))) &
               + DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
              ENDDO
             ELSE IF (DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN5_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))) &
               + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
              ENDDO
             ELSE
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN5_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))) &
               +(DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)) &
               + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)))/LARGECOMZ
              ENDDO
             ENDIF
            ENDDO
           ELSE
            DO KL=1,NKL_G(KLT)
             IF ((DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0).AND.(DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0)) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               +(DOWN4_C(IJ,IJT)*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
               + DOWN5_C(IJ,IJT)*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))))/COMZ
              ENDDO
             ELSE IF (DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               +(DOWN4_C(IJ,IJT)*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
               + DOWN5_C(IJ,IJT)*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))))/COMZ &
               + DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
              ENDDO
             ELSE IF (DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               +(DOWN4_C(IJ,IJT)*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
               + DOWN5_C(IJ,IJT)*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))))/COMZ &
               + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
              ENDDO
             ELSE
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               +(DOWN4_C(IJ,IJT)*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
               + DOWN5_C(IJ,IJT)*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))))/COMZ &
               +(DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)) &
               + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)))/LARGECOMZ
              ENDDO
             ENDIF
            ENDDO
           ENDIF
          ENDDO
         ENDIF
         ! INTEGRAL DERIVATIVE EVALUATION
         DO KL=NKL_G(KLT)+1,NKL(KLT)
          DO M=0,IJKLANG-KLANGINDX(KL,KLT)
           ERI(M,KL,1)=QI_CI(KL,SKL)*ERI(M,DOWN3(KL,KLT),1)+WI_QI(DOWNXYZ(KL,KLT))*ERI(M+1,DOWN3(KL,KLT),1)
           IF (DOWN4(KL,KLT) /= 0) ERI(M,KL,1)=ERI(M,KL,1) &
           +DOWN4_C(KL,KLT)/COMZ_S(SKL)*(ERI(M,DOWN4(KL,KLT),1)-RHO_COMZ2*ERI(M+1,DOWN4(KL,KLT),1))
           IF (DOWN5(KL,KLT) /= 0) ERI(M,KL,1)=ERI(M,KL,1) &
           +DOWN5_C(KL,KLT)/COMZ_S(SKL)*(ERI(M,DOWN5(KL,KLT),1)-RHO_COMZ2*ERI(M+1,DOWN5(KL,KLT),1))
          ENDDO
         ENDDO
         IF (NIJ_G /= 1) THEN
          DO IJ=2,NIJ_G
           IF ((DOWN4(IJ,IJT) == 0).AND.(DOWN5(IJ,IJT) == 0)) THEN
            DO KL=NKL_G(KLT)+1,NKL(KLT)
             IF ((DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0).AND.(DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0)) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT))
              ENDDO
             ELSE IF (DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
              ENDDO
             ELSE IF (DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
              ENDDO
             ELSE
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               +(DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)) &
               + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)))/LARGECOMZ
              ENDDO
             ENDIF
            ENDDO
           ELSE IF (DOWN5(IJ,IJT) == 0) THEN
            DO KL=NKL_G(KLT)+1,NKL(KLT)
             IF ((DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0).AND.(DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0)) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN4_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT)))
              ENDDO
             ELSE IF (DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN4_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
               + DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
              ENDDO
             ELSE IF (DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN4_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
               + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
              ENDDO
             ELSE
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN4_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
               +(DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)) &
               + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)))/LARGECOMZ
              ENDDO
             ENDIF
            ENDDO
           ELSE IF (DOWN4(IJ,IJT) == 0) THEN
            DO KL=NKL_G(KLT)+1,NKL(KLT)
             IF ((DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0).AND.(DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0)) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN5_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT)))
              ENDDO
             ELSE IF (DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN5_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))) &
               + DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
              ENDDO
             ELSE IF (DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN5_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))) &
               + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
              ENDDO
             ELSE
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               + DOWN5_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))) &
               +(DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)) &
               + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)))/LARGECOMZ
              ENDDO
             ENDIF
            ENDDO
           ELSE
            DO KL=NKL_G(KLT)+1,NKL(KLT)
             IF ((DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0).AND.(DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0)) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               +(DOWN4_C(IJ,IJT)*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
               + DOWN5_C(IJ,IJT)*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))))/COMZ
              ENDDO
             ELSE IF (DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               +(DOWN4_C(IJ,IJT)*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
               + DOWN5_C(IJ,IJT)*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))))/COMZ &
               + DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
              ENDDO
             ELSE IF (DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               +(DOWN4_C(IJ,IJT)*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
               + DOWN5_C(IJ,IJT)*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))))/COMZ &
               + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
              ENDDO
             ELSE
              DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
               ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
               +(DOWN4_C(IJ,IJT)*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
               + DOWN5_C(IJ,IJT)*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))))/COMZ &
               +(DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)) &
               + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)))/LARGECOMZ
              ENDDO
             ENDIF
            ENDDO
           ENDIF
          ENDDO
         ENDIF
         ! INTEGRAL DERIVATIVE EVALUATION
         IF (NKL_G(KLT) /= 1) THEN
          DO KL=2,NKL_G(KLT)
           DO M=0,IJKLANG-KLANGINDX(KL,KLT)
            ERI(M,KL,1)=QI_CI(KL,SKL)*ERI(M,DOWN3(KL,KLT),1)+WI_QI(DOWNXYZ(KL,KLT))*ERI(M+1,DOWN3(KL,KLT),1)
            IF (DOWN4(KL,KLT) /= 0) ERI(M,KL,1)=ERI(M,KL,1) &
           +DOWN4_C(KL,KLT)/COMZ_S(SKL)*(ERI(M,DOWN4(KL,KLT),1)-RHO_COMZ2*ERI(M+1,DOWN4(KL,KLT),1))
            IF (DOWN5(KL,KLT) /= 0) ERI(M,KL,1)=ERI(M,KL,1) &
           +DOWN5_C(KL,KLT)/COMZ_S(SKL)*(ERI(M,DOWN5(KL,KLT),1)-RHO_COMZ2*ERI(M+1,DOWN5(KL,KLT),1))
           ENDDO
          ENDDO
         ENDIF
         DO IJ=NIJ_G+1,NIJ
          IF ((DOWN4(IJ,IJT) == 0).AND.(DOWN5(IJ,IJT) == 0)) THEN
           DO KL=1,NKL_G(KLT)
            IF ((DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0).AND.(DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0)) THEN
             DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
              ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT))
             ENDDO
            ELSE IF (DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
             DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
              ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
              + DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
             ENDDO
            ELSE IF (DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
             DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
              ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
              + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
             ENDDO
            ELSE
             DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
              ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
              +(DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)) &
              + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)))/LARGECOMZ
             ENDDO
            ENDIF
           ENDDO
          ELSE IF (DOWN5(IJ,IJT) == 0) THEN
           DO KL=1,NKL_G(KLT)
            IF ((DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0).AND.(DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0)) THEN
             DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
              ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
              + DOWN4_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT)))
             ENDDO
            ELSE IF (DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
             DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
              ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
              + DOWN4_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
              + DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
             ENDDO
            ELSE IF (DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
             DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
              ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
              + DOWN4_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
              + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
             ENDDO
            ELSE
             DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
              ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
              + DOWN4_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
              +(DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)) &
              + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)))/LARGECOMZ
             ENDDO
            ENDIF
           ENDDO
          ELSE IF (DOWN4(IJ,IJT) == 0) THEN
           DO KL=1,NKL_G(KLT)
            IF ((DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0).AND.(DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0)) THEN
             DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
              ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
              + DOWN5_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT)))
             ENDDO
            ELSE IF (DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
             DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
              ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
              + DOWN5_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))) &
              + DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
             ENDDO
            ELSE IF (DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
             DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
              ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
              + DOWN5_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))) &
              + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
             ENDDO
            ELSE
             DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
              ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
              + DOWN5_C(IJ,IJT)/COMZ*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))) &
              +(DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)) &
              + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)))/LARGECOMZ
             ENDDO
            ENDIF
           ENDDO
          ELSE
           DO KL=1,NKL_G(KLT)
            IF ((DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0).AND.(DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0)) THEN
             DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
              ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
              +(DOWN4_C(IJ,IJT)*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
              + DOWN5_C(IJ,IJT)*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))))/COMZ
             ENDDO
            ELSE IF (DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
             DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
              ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
              +(DOWN4_C(IJ,IJT)*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
              + DOWN5_C(IJ,IJT)*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))))/COMZ &
              + DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
             ENDDO
            ELSE IF (DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)) == 0) THEN
             DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
              ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
              +(DOWN4_C(IJ,IJT)*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
              + DOWN5_C(IJ,IJT)*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))))/COMZ &
              + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))/LARGECOMZ*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT))
             ENDDO
            ELSE
             DO M=0,IJKLANG-KLANGINDX(IJ,IJT)-KLANGINDX(KL,KLT)
              ERI(M,KL,IJ)=PI_AI(IJ)*ERI(M,KL,DOWN3(IJ,IJT))+WI_PI(DOWNXYZ(IJ,IJT))*ERI(M+1,KL,DOWN3(IJ,IJT)) &
              +(DOWN4_C(IJ,IJT)*(ERI(M,KL,DOWN4(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN4(IJ,IJT))) &
              + DOWN5_C(IJ,IJT)*(ERI(M,KL,DOWN5(IJ,IJT))-RHO_COMZ1*ERI(M+1,KL,DOWN5(IJ,IJT))))/COMZ &
              +(DOWN1_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN1(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)) &
              + DOWN2_C(KL,KLT,DOWNXYZ(IJ,IJT))*ERI(M+1,DOWN2(KL,KLT,DOWNXYZ(IJ,IJT)),DOWN3(IJ,IJT)))/LARGECOMZ
             ENDDO
            ENDIF
           ENDDO
          ENDIF
         ENDDO
         ! FORCE ACCUMULATION
         DO IJ=1,NIJ_G
          DO KL=1,NKL_G(KLT)
           IF (HELIX == 0.0D0) THEN
            GN1=0.5D0*P_P(PGSI(IJ),PGSJ(IJ),Q1X,Q1Y,Q1Z)*P_P(PGSK(KL,SKL),PGSL(KL,SKL),Q2X,Q2Y,Q2Z)
            IF (LCEL1X(Q3X).AND.LCEL1Y(Q3Y).AND.LCEL1Z(Q3Z)) GN1=GN1-0.25D0*P_P(PGSI(IJ),PGSK(KL,SKL),Q3X,Q3Y,Q3Z) &
                                         *P_P(PGSJ(IJ),PGSL(KL,SKL),Q2X+Q3X-Q1X,Q2Y+Q3Y-Q1Y,Q2Z+Q3Z-Q1Z)*DOPTN(19)
           ELSE ! --- HELIX
            GN1=0.5D0*P_P_HELIX(PGSI(IJ),PGSJ(IJ),0,Q1X)*P_P_HELIX(PGSK(KL,SKL),PGSL(KL,SKL),Q3X,Q3X+Q2X)
            IF (LCEL1X(Q3X).AND.LCEL1Y(Q3Y).AND.LCEL1Z(Q3Z)) GN1=GN1-0.25D0*P_P_HELIX(PGSI(IJ),PGSK(KL,SKL),0,Q3X) &
                                         *P_P_HELIX(PGSJ(IJ),PGSL(KL,SKL),Q1X,Q3X+Q2X)*DOPTN(19)
            GN2=0.5D0*P_P_HELIX_THETA(PGSI(IJ),PGSJ(IJ),0,Q1X)*P_P_HELIX(PGSK(KL,SKL),PGSL(KL,SKL),Q3X,Q3X+Q2X) &
               +0.5D0*P_P_HELIX(PGSI(IJ),PGSJ(IJ),0,Q1X)*P_P_HELIX_THETA(PGSK(KL,SKL),PGSL(KL,SKL),Q3X,Q3X+Q2X)
            IF (LCEL1X(Q3X).AND.LCEL1Y(Q3Y).AND.LCEL1Z(Q3Z)) GN2=GN2-0.25D0*P_P_HELIX_THETA(PGSI(IJ),PGSK(KL,SKL),0,Q3X) &
                                                                    *P_P_HELIX(PGSJ(IJ),PGSL(KL,SKL),Q1X,Q3X+Q2X)*DOPTN(19) &
                                                                    -0.25D0*P_P_HELIX(PGSI(IJ),PGSK(KL,SKL),0,Q3X) &
                                                                    *P_P_HELIX_THETA(PGSJ(IJ),PGSL(KL,SKL),Q1X,Q3X+Q2X)*DOPTN(19)
           ENDIF ! --- HELIX END

           FN=2.0D0*ZT(II)*ERI(0,KL,UP1(IJ,IJT,1))
           IF (DOWN1(IJ,IJT,1) /= 0) FN=FN-2.0D0*DOWN1_C(IJ,IJT,1)*ERI(0,KL,DOWN1(IJ,IJT,1))
           XI1=FN*GN1
           FN=2.0D0*ZT(II)*ERI(0,KL,UP1(IJ,IJT,2))
           IF (DOWN1(IJ,IJT,2) /= 0) FN=FN-2.0D0*DOWN1_C(IJ,IJT,2)*ERI(0,KL,DOWN1(IJ,IJT,2))
           YI1=FN*GN1
           FN=2.0D0*ZT(II)*ERI(0,KL,UP1(IJ,IJT,3))
           IF (DOWN1(IJ,IJT,3) /= 0) FN=FN-2.0D0*DOWN1_C(IJ,IJT,3)*ERI(0,KL,DOWN1(IJ,IJT,3))
           ZI1=FN*GN1

           FN=2.0D0*ZTK(SKL)*ERI(0,UP1(KL,KLT,1),IJ)
           IF (DOWN1(KL,KLT,1) /= 0) FN=FN-2.0D0*DOWN1_C(KL,KLT,1)*ERI(0,DOWN1(KL,KLT,1),IJ)
           XK1=FN*GN1
           FN=2.0D0*ZTK(SKL)*ERI(0,UP1(KL,KLT,2),IJ)
           IF (DOWN1(KL,KLT,2) /= 0) FN=FN-2.0D0*DOWN1_C(KL,KLT,2)*ERI(0,DOWN1(KL,KLT,2),IJ)
           YK1=FN*GN1
           FN=2.0D0*ZTK(SKL)*ERI(0,UP1(KL,KLT,3),IJ)
           IF (DOWN1(KL,KLT,3) /= 0) FN=FN-2.0D0*DOWN1_C(KL,KLT,3)*ERI(0,DOWN1(KL,KLT,3),IJ)
           ZK1=FN*GN1

           FN=2.0D0*ZTL(SKL)*ERI(0,UP2(KL,KLT,1),IJ)
           IF (DOWN2(KL,KLT,1) /= 0) FN=FN-2.0D0*DOWN2_C(KL,KLT,1)*ERI(0,DOWN2(KL,KLT,1),IJ)
           XL1=FN*GN1
           FN=2.0D0*ZTL(SKL)*ERI(0,UP2(KL,KLT,2),IJ)
           IF (DOWN2(KL,KLT,2) /= 0) FN=FN-2.0D0*DOWN2_C(KL,KLT,2)*ERI(0,DOWN2(KL,KLT,2),IJ)
           YL1=FN*GN1
           FN=2.0D0*ZTL(SKL)*ERI(0,UP2(KL,KLT,3),IJ)
           IF (DOWN2(KL,KLT,3) /= 0) FN=FN-2.0D0*DOWN2_C(KL,KLT,3)*ERI(0,DOWN2(KL,KLT,3),IJ)
           ZL1=FN*GN1

           IF (HELIX == 0.0D0) THEN
            FORCEX_TMP(KATM(SKL))=FORCEX_TMP(KATM(SKL))+XK1
            FORCEY_TMP(KATM(SKL))=FORCEY_TMP(KATM(SKL))+YK1
            FORCEZ_TMP(KATM(SKL))=FORCEZ_TMP(KATM(SKL))+ZK1
            FORCETX_TMP=FORCETX_TMP+XK1*DFLOAT(Q3X)
            FORCETY_TMP=FORCETY_TMP+YK1*DFLOAT(Q3Y)
            FORCETZ_TMP=FORCETZ_TMP+ZK1*DFLOAT(Q3Z)
            FORCEX_TMP(LATM(SKL))=FORCEX_TMP(LATM(SKL))+XL1
            FORCEY_TMP(LATM(SKL))=FORCEY_TMP(LATM(SKL))+YL1
            FORCEZ_TMP(LATM(SKL))=FORCEZ_TMP(LATM(SKL))+ZL1
            FORCETX_TMP=FORCETX_TMP+XL1*DFLOAT(Q3X+Q2X)
            FORCETY_TMP=FORCETY_TMP+YL1*DFLOAT(Q3Y+Q2Y)
            FORCETZ_TMP=FORCETZ_TMP+ZL1*DFLOAT(Q3Z+Q2Z)
            FORCEX_TMP(IATM)=FORCEX_TMP(IATM)+XI1
            FORCEY_TMP(IATM)=FORCEY_TMP(IATM)+YI1
            FORCEZ_TMP(IATM)=FORCEZ_TMP(IATM)+ZI1
            FORCEX_TMP(JATM)=FORCEX_TMP(JATM)-(XI1+XK1+XL1)
            FORCEY_TMP(JATM)=FORCEY_TMP(JATM)-(YI1+YK1+YL1)
            FORCEZ_TMP(JATM)=FORCEZ_TMP(JATM)-(ZI1+ZK1+ZL1)
            FORCETX_TMP=FORCETX_TMP-((XI1+XK1+XL1)*DFLOAT(Q1X))
            FORCETY_TMP=FORCETY_TMP-((YI1+YK1+YL1)*DFLOAT(Q1Y))
            FORCETZ_TMP=FORCETZ_TMP-((ZI1+ZK1+ZL1)*DFLOAT(Q1Z))
           ELSE ! --- HELIX
            FORCEX_TMP(KATM(SKL))=FORCEX_TMP(KATM(SKL))+XK1
            FORCEY_TMP(KATM(SKL))=FORCEY_TMP(KATM(SKL))+YK1*P3C(Q3X)+ZK1*P3S(Q3X)
            FORCEZ_TMP(KATM(SKL))=FORCEZ_TMP(KATM(SKL))-YK1*P3S(Q3X)+ZK1*P3C(Q3X)
            FORCETX_TMP=FORCETX_TMP+XK1*DFLOAT(Q3X)
            FORCETA_TMP=FORCETA_TMP+DFLOAT(Q3X)*(-ATOMY(KATM(SKL))*P3S(Q3X)-ATOMZ(KATM(SKL))*P3C(Q3X))*YK1
            FORCETA_TMP=FORCETA_TMP+DFLOAT(Q3X)*( ATOMY(KATM(SKL))*P3C(Q3X)-ATOMZ(KATM(SKL))*P3S(Q3X))*ZK1
            FORCEX_TMP(LATM(SKL))=FORCEX_TMP(LATM(SKL))+XL1
            FORCEY_TMP(LATM(SKL))=FORCEY_TMP(LATM(SKL))+YL1*P3C(Q3X+Q2X)+ZL1*P3S(Q3X+Q2X)
            FORCEZ_TMP(LATM(SKL))=FORCEZ_TMP(LATM(SKL))-YL1*P3S(Q3X+Q2X)+ZL1*P3C(Q3X+Q2X)
            FORCETX_TMP=FORCETX_TMP+XL1*DFLOAT(Q3X+Q2X)
            FORCETA_TMP=FORCETA_TMP+DFLOAT(Q3X+Q2X)*(-ATOMY(LATM(SKL))*P3S(Q3X+Q2X)-ATOMZ(LATM(SKL))*P3C(Q3X+Q2X))*YL1
            FORCETA_TMP=FORCETA_TMP+DFLOAT(Q3X+Q2X)*( ATOMY(LATM(SKL))*P3C(Q3X+Q2X)-ATOMZ(LATM(SKL))*P3S(Q3X+Q2X))*ZL1
            FORCEX_TMP(IATM)=FORCEX_TMP(IATM)+XI1
            FORCEY_TMP(IATM)=FORCEY_TMP(IATM)+YI1
            FORCEZ_TMP(IATM)=FORCEZ_TMP(IATM)+ZI1
            FORCEX_TMP(JATM)=FORCEX_TMP(JATM)-(XI1+XK1+XL1)
            FORCEY_TMP(JATM)=FORCEY_TMP(JATM)-(YI1+YK1+YL1)*P3C(Q1X)-(ZI1+ZK1+ZL1)*P3S(Q1X)
            FORCEZ_TMP(JATM)=FORCEZ_TMP(JATM)+(YI1+YK1+YL1)*P3S(Q1X)-(ZI1+ZK1+ZL1)*P3C(Q1X)
            FORCETX_TMP=FORCETX_TMP-((XI1+XK1+XL1)*DFLOAT(Q1X))
            FORCETA_TMP=FORCETA_TMP+DFLOAT(Q1X)*(-ATOMY(JATM)*P3S(Q1X)-ATOMZ(JATM)*P3C(Q1X))*(-(YI1+YK1+YL1))
            FORCETA_TMP=FORCETA_TMP+DFLOAT(Q1X)*( ATOMY(JATM)*P3C(Q1X)-ATOMZ(JATM)*P3S(Q1X))*(-(ZI1+ZK1+ZL1))
            ! BASIS ROTATION FORCE
            FORCETA_TMP=FORCETA_TMP+GN2*ERI(0,KL,IJ)
           ENDIF  ! --- HELIX END

          ENDDO
         ENDDO
        ENDDO
       ENDDO
       ENDDO
       ENDDO
      ENDDO
     ENDDO
!--- PARALLEL LOOP
     ENDIF
!--- PARALLEL LOOP
    ENDDO
    ENDDO
    ENDDO
    CALL PCPU_TIME(ICPUE)
    IF (MYID == 0) WRITE(6,'(A,F6.2,A,F10.1,A)') ' ... ', &
      DFLOAT(((Q1X+REDUCED_CEL1X)*(2*REDUCED_CEL1Y+1)+(Q1Y+REDUCED_CEL1Y))*(2*REDUCED_CEL1Z+1)+(Q1Z+REDUCED_CEL1Z+1))/ &
      DFLOAT((2*REDUCED_CEL1X+1)*(2*REDUCED_CEL1Y+1)*(2*REDUCED_CEL1Z+1))*100.0D0, &
      ' % OF ERI GRADIENTS GENERATION (CPU / SEC = ',ICPUE-ICPUS,')'
    CALL PCPU_TIME(ICPUS)
    CALL PFLUSH(6)
   ENDDO
   ENDDO
   ENDDO

   CALL MPI_ALLREDUCE(MPI_IN_PLACE,FORCEX_TMP,NATOM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,FORCEY_TMP,NATOM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,FORCEZ_TMP,NATOM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,FORCETX_TMP,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,FORCETY_TMP,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,FORCETZ_TMP,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,FORCETA_TMP,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   FORCEX=FORCEX+FORCEX_TMP
   FORCEY=FORCEY+FORCEY_TMP
   FORCEZ=FORCEZ+FORCEZ_TMP
   FORCETX=FORCETX+FORCETX_TMP
   FORCETY=FORCETY+FORCETY_TMP
   FORCETZ=FORCETZ+FORCETZ_TMP
   FORCETA=FORCETA+FORCETA_TMP

   DEALLOCATE(FORCEX_TMP,FORCEY_TMP,FORCEZ_TMP)
   DEALLOCATE(P3X,P3Y,P3Z,P3C,P3S,COMZ_S,PX_S,PY_S,PZ_S,OV,QI_CI)
   DEALLOCATE(KLT_S,KLANG,PGSK,PGSL,KLKL,ZTK,ZTL)
   DEALLOCATE(LCEL1X,LCEL1Y,LCEL1Z,MAXC)
   IF ((IOPTN(9) == 3).AND.(MYID == 0)) CALL DUMP2

   RETURN
END SUBROUTINE
