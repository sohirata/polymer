SUBROUTINE DIRECT_CONTRACT_FOURINDEX_ERI(LTDR,A_X,P_AX,B_Y,P_BY,B_X,P_BX,A_Y,P_AY,N,S,IFLAG)
! CALCULATE FOUR-CENTER TWO-ELECTRON ELECTRON REPULSION INTEGRALS
! USING THE OBARA-SAIKA RECURSION METHOD AND CONTRACT THEM WITH
! TRIAL DENSITY MATRICES PROVIDED AS ARGUMENTS.

   USE CONTROL
   USE STRUCTURE
   USE BASISSET
   USE INTEGRAL
   USE GRADIENT
   USE FMT
   USE CONSTANTS
   
   IMPLICIT NONE
   LOGICAL :: LTDR
   INTEGER :: IFLAG,N,S
   INTEGER :: Q1,Q2,Q3
   INTEGER :: I,J,K,L,II,JJ,KK,LL,IJ,KL,M
   INTEGER :: I1,I2,I3,J1,J2,J3,IA,JA
   INTEGER :: IJT,KLT
   INTEGER :: NPSHELL_SQ
   INTEGER :: IJANG,SIJ,SKL,NIJ,IJKLANG,TS
   INTEGER :: DOWN1(300,9,3),DOWN2(300,9,3)
   INTEGER :: DOWN3(300,9),DOWN4(300,9),DOWN5(300,9),DOWNXYZ(300,9)
   INTEGER :: NKL(9),PGSI(300),PGSJ(300)
   INTEGER :: IJMAP(0:3,0:3,0:3,0:3,0:3,0:3)
   INTEGER :: KLANGINDX(300,9)
   INTEGER,ALLOCATABLE :: KLT_S(:),KLANG(:),PGSK(:,:),PGSL(:,:)
   DOUBLE COMPLEX :: A_X(N,N,-S:S),P_AX(N,N,-S:S),A_Y(N,N,-S:S),P_AY(N,N,-S:S)
   DOUBLE COMPLEX :: B_X(N,N,-S:S),P_BX(N,N,-S:S),B_Y(N,N,-S:S),P_BY(N,N,-S:S)
   DOUBLE COMPLEX,ALLOCATABLE :: A_XP(:,:,:),P_AXP(:,:,:),A_YP(:,:,:),P_AYP(:,:,:)
   DOUBLE COMPLEX,ALLOCATABLE :: B_XP(:,:,:),P_BXP(:,:,:),B_YP(:,:,:),P_BYP(:,:,:)
   DOUBLE PRECISION :: ERI(0:9,300,300)
   DOUBLE PRECISION :: Q1X,Q1C,Q1S,Q2X,Q2C,Q2S
   DOUBLE PRECISION :: AX,AY,AZ,BX,BY,BZ,CX,CY,CZ,DX,DY,DZ,PX,PY,PZ,QX,QY,QZ
   DOUBLE PRECISION :: COMZ,C1,C2,LARGECOMZ,RHO_COMZ1,RHO_COMZ2,RHO,T,H
   DOUBLE PRECISION :: PI_AI(300),IJIJ,PMAX
   DOUBLE PRECISION :: PISUB,F1(0:9),F2(0:9),WX,WY,WZ,WI_PI(3),WI_QI(3)
   DOUBLE PRECISION :: DOWN1_C(300,9,3),DOWN2_C(300,9,3),DOWN4_C(300,9),DOWN5_C(300,9)
   DOUBLE PRECISION,ALLOCATABLE :: COMZ_S(:),PX_S(:),PY_S(:),PZ_S(:),OV(:),QI_CI(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: KLKL(:)
   DOUBLE PRECISION,ALLOCATABLE :: Q3X(:),Q3C(:),Q3S(:)
   DOUBLE COMPLEX,ALLOCATABLE :: W1(:,:),W2(:,:),W3(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: MAXC(:)
   DOUBLE PRECISION :: A
   LOGICAL,ALLOCATABLE :: LCEL1X(:),LCEL2X(:)

   IF (HELIX /= 0.0D0) CALL PABORT('AN INTERNAL PROGRAM ERROR IS DETECTED')
   IF ((N /= NCGS).OR.(S /= CEL1)) CALL PABORT('INCORRECT ARGUMENTS IN DIRECT_CONTRACT_FOURINDEX_ERI')

   NPSHELL_SQ=NPSHELL**2
   ALLOCATE(Q3X(0:CEL2),Q3C(0:CEL2),Q3S(0:CEL2))
   ALLOCATE(COMZ_S(NPSHELL_SQ),PX_S(NPSHELL_SQ),PY_S(NPSHELL_SQ),PZ_S(NPSHELL_SQ),OV(NPSHELL_SQ),QI_CI(300,NPSHELL_SQ))
   ALLOCATE(KLT_S(NPSHELL_SQ),KLANG(NPSHELL_SQ),PGSK(300,NPSHELL_SQ),PGSL(300,NPSHELL_SQ),KLKL(NPSHELL_SQ))
   ALLOCATE(LCEL1X(0:CEL2),LCEL2X(0:CEL2),MAXC(NPGS))
   ALLOCATE(A_XP(NPGS,NPGS,-CEL1:CEL1),P_AXP(NPGS,NPGS,-CEL1:CEL1))
   A_XP=DCMPLX(0.0D0,0.0D0)
   IF (LTDR) THEN
    ALLOCATE(A_YP(NPGS,NPGS,-CEL1:CEL1),P_AYP(NPGS,NPGS,-CEL1:CEL1))
    ALLOCATE(B_XP(NPGS,NPGS,-CEL1:CEL1),P_BXP(NPGS,NPGS,-CEL1:CEL1))
    ALLOCATE(B_YP(NPGS,NPGS,-CEL1:CEL1),P_BYP(NPGS,NPGS,-CEL1:CEL1))
    A_YP=DCMPLX(0.0D0,0.0D0)
    B_XP=DCMPLX(0.0D0,0.0D0)
    B_YP=DCMPLX(0.0D0,0.0D0)
   ENDIF

   ! DECONTRACT DENSITY MATRICES
   ALLOCATE(W1(NCGS,NCGS),W2(NCGS,NPGS),W3(NPGS,NPGS))
   DO Q1=-CEL1,CEL1
    W1=P_AX(:,:,Q1)
    W2=MATMUL(W1,CC)
    W3=MATMUL(TRANSPOSE(CC),W2)
    P_AXP(:,:,Q1)=W3
   ENDDO
   IF (LTDR) THEN
    DO Q1=-CEL1,CEL1
     W1=P_AY(:,:,Q1)
     W2=MATMUL(W1,CC)
     W3=MATMUL(TRANSPOSE(CC),W2)
     P_AYP(:,:,Q1)=W3
    ENDDO
    DO Q1=-CEL1,CEL1
     W1=P_BX(:,:,Q1)
     W2=MATMUL(W1,CC)
     W3=MATMUL(TRANSPOSE(CC),W2)
     P_BXP(:,:,Q1)=W3
    ENDDO
    DO Q1=-CEL1,CEL1
     W1=P_BY(:,:,Q1)
     W2=MATMUL(W1,CC)
     W3=MATMUL(TRANSPOSE(CC),W2)
     P_BYP(:,:,Q1)=W3
    ENDDO
   ENDIF
   DEALLOCATE(W1,W2,W3)
  
   ! PRELIMINARY
   PISUB=2.0D0/DSQRT(PI)
   F1(0)=1.0D0/PISUB
   F1(1)=1.0D0/2.0D0/PISUB
   F1(2)=3.0D0/4.0D0/PISUB
   F1(3)=15.0D0/8.0D0/PISUB
   F1(4)=105.0D0/16.0D0/PISUB
   F1(5)=945.0D0/32.0D0/PISUB
   F1(6)=10395.0D0/64.0D0/PISUB
   F1(7)=135135.0D0/128.0D0/PISUB
   F1(8)=2027025.D0/256.0D0/PISUB
   F1(9)=34459425.0D0/512.0D0/PISUB
   F2(0)=-0.5D0
   F2(1)=-1.5D0
   F2(2)=-2.5D0
   F2(3)=-3.5D0
   F2(4)=-4.5D0
   F2(5)=-5.5D0
   F2(6)=-6.5D0
   F2(7)=-7.5D0
   F2(8)=-8.5D0
   F2(9)=-9.5D0
   DO Q3=0,CEL2
    Q3X(Q3)=DFLOAT(Q3)*PERIOD
    Q3C(Q3)=DCOS(DFLOAT(Q3)*HELIX)
    Q3S(Q3)=DSIN(DFLOAT(Q3)*HELIX)
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
   DO IA=0,2
    DO JA=0,2
     KLT=IA*3+JA+1
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
    ENDDO
   ENDDO
   
   ! LOOP OVER UNIT CELL INDEX 1
   DO Q1=-REDUCED_CEL1,REDUCED_CEL1
    Q1X=DFLOAT(Q1)*PERIOD
    Q1C=DCOS(DFLOAT(Q1)*HELIX)
    Q1S=DSIN(DFLOAT(Q1)*HELIX)
    ! LOOP OVER UNIT CELL INDEX 2
    DO Q2=-REDUCED_CEL1,REDUCED_CEL1
     Q2X=DFLOAT(Q2)*PERIOD
     Q2C=DCOS(DFLOAT(Q2)*HELIX)
     Q2S=DSIN(DFLOAT(Q2)*HELIX)
     ! SET CELL SYMMETRY INDECES
     LCEL1X=.FALSE.
     LCEL2X=.FALSE.
     DO Q3=0,CEL2
      IF ((Q3 <= CEL1).AND.(-Q2-Q3+Q1 >= -CEL1).AND.(-Q2-Q3+Q1 <= CEL1)) LCEL1X(Q3)=.TRUE.
      IF ((Q3 /= 0).AND.(Q3 <= CEL1).AND.(Q2+Q3-Q1 >= -CEL1).AND.(Q2+Q3-Q1 <= CEL1)) LCEL2X(Q3)=.TRUE.
     ENDDO
     ! LOOP OVER KL SHELL
     DO K=1,NPSHELL
      KK=P_SHELL(K,0,0,0)
      CX=PGSX(KK)
      CY=PGSY(KK)
      CZ=PGSZ(KK)
      DO L=1,NPSHELL
       LL=P_SHELL(L,0,0,0)
       DX=PGSX(LL)+Q2X
       DY=PGSY(LL)*Q2C-PGSZ(LL)*Q2S
       DZ=PGSY(LL)*Q2S+PGSZ(LL)*Q2C
       COMZ=ZT(KK)+ZT(LL)
       PX=(CX*ZT(KK)+DX*ZT(LL))/COMZ
       PY=(CY*ZT(KK)+DY*ZT(LL))/COMZ
       PZ=(CZ*ZT(KK)+DZ*ZT(LL))/COMZ
       SKL=(K-1)*NPSHELL+L   ! SHELL INDEX
       COMZ_S(SKL)=COMZ      ! STORE ZETA+ZETA
       PX_S(SKL)=PX          ! STORE PX, PY, & PZ
       PY_S(SKL)=PY
       PZ_S(SKL)=PZ
       OV(SKL)=S_P(KK,LL,Q2)   ! STORE UNNORMALIZED OVERLAP INTEGRALS
       KLANG(SKL)=P_SHELL_ANG(K)+P_SHELL_ANG(L)   ! STORE SUM OF ANGULAR MOMENTA
       KLT_S(SKL)=P_SHELL_ANG(K)*3+P_SHELL_ANG(L)+1 ! STORE SHELL TYPE
       NIJ=0
       ! ELECTRON 2
       DO I1=0,P_SHELL_ANG(K)
        DO I2=0,P_SHELL_ANG(K)-I1
         DO I3=0,P_SHELL_ANG(K)-I1-I2
          DO J1=0,P_SHELL_ANG(L)
           DO J2=0,P_SHELL_ANG(L)-J1
            DO J3=0,P_SHELL_ANG(L)-J1-J2
             NIJ=NIJ+1
             PGSK(NIJ,SKL)=P_SHELL(K,I1,I2,I3)
             PGSL(NIJ,SKL)=P_SHELL(L,J1,J2,J3)
             IF (I1 >= 1) THEN
              QI_CI(NIJ,SKL)=PX-CX
             ELSE IF (I2 >= 1) THEN
              QI_CI(NIJ,SKL)=PY-CY
             ELSE IF (I3 >= 1) THEN
              QI_CI(NIJ,SKL)=PZ-CZ
             ELSE IF (J1 >= 1) THEN
              QI_CI(NIJ,SKL)=PX-DX
             ELSE IF (J2 >= 1) THEN
              QI_CI(NIJ,SKL)=PY-DY
             ELSE IF (J3 >= 1) THEN
              QI_CI(NIJ,SKL)=PZ-DZ
             ENDIF
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       IF (NIJ /= NKL(KLT_S(SKL))) CALL PABORT('AN INTERNAL PROGRAM ERROR IS DETECTED')
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
       DO IJ=1,NKL(KLT)
        IF (IJ == 1) THEN
         DO KL=2,NKL(KLT)
          DO M=0,IJKLANG-KLANGINDX(KL,KLT)
           ERI(M,KL,1)=QI_CI(KL,SKL)*ERI(M,DOWN3(KL,KLT),1)
           IF (DOWN4(KL,KLT) /= 0) ERI(M,KL,1)=ERI(M,KL,1) &
           +DOWN4_C(KL,KLT)/COMZ_S(SKL)*(ERI(M,DOWN4(KL,KLT),1)-0.5D0*ERI(M+1,DOWN4(KL,KLT),1))
           IF (DOWN5(KL,KLT) /= 0) ERI(M,KL,1)=ERI(M,KL,1) &
           +DOWN5_C(KL,KLT)/COMZ_S(SKL)*(ERI(M,DOWN5(KL,KLT),1)-0.5D0*ERI(M+1,DOWN5(KL,KLT),1))
          ENDDO
         ENDDO
        ELSE
         DO KL=1,NKL(KLT)
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
       DO KL=2,NKL(KLT)
        A=DSQRT(ERI(0,KL,KL)*MAXC(PGSK(KL,SKL))*MAXC(PGSL(KL,SKL)))
        IF (A > KLKL(SKL)) KLKL(SKL)=A
       ENDDO
      ELSE
       A=DSQRT(C2*MAXC(PGSK(1,SKL))*MAXC(PGSL(1,SKL)))
       IF (A > KLKL(SKL)) KLKL(SKL)=A
      ENDIF
     ENDDO

     ! MAIN LOOP STARTED
     DO I=1,NPSHELL
      II=P_SHELL(I,0,0,0)
      AX=PGSX(II)
      AY=PGSY(II)
      AZ=PGSZ(II)
      DO J=1,NPSHELL
       JJ=P_SHELL(J,0,0,0)
       BX=PGSX(JJ)+Q1X
       BY=PGSY(JJ)*Q1C-PGSZ(JJ)*Q1S
       BZ=PGSY(JJ)*Q1S+PGSZ(JJ)*Q1C
       COMZ=ZT(II)+ZT(JJ)
       PX=(AX*ZT(II)+BX*ZT(JJ))/COMZ
       PY=(AY*ZT(II)+BY*ZT(JJ))/COMZ
       PZ=(AZ*ZT(II)+BZ*ZT(JJ))/COMZ
       C1=S_P(II,JJ,Q1)*2.0D0/DSQRT(PI)
       SIJ=(I-1)*NPSHELL+J   ! SHELL INDEX
       IJT=P_SHELL_ANG(I)*3+P_SHELL_ANG(J)+1
       IJANG=P_SHELL_ANG(I)+P_SHELL_ANG(J)
       
       ! ELECTRON 1
       NIJ=0
       DO I1=0,P_SHELL_ANG(I)
        DO I2=0,P_SHELL_ANG(I)-I1
         DO I3=0,P_SHELL_ANG(I)-I1-I2
          DO J1=0,P_SHELL_ANG(J)
           DO J2=0,P_SHELL_ANG(J)-J1
            DO J3=0,P_SHELL_ANG(J)-J1-J2
             NIJ=NIJ+1
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
       IF (NIJ /= NKL(KLT_S(SIJ))) CALL PABORT('AN INTERNAL PROGRAM ERROR IS DETECTED')
    
       ! FIND THE MAXIMUM (IJ|IJ) VALUE FOR THE CURRENT SHELL
       IJIJ=0.0D0
       LARGECOMZ=2.0D0*COMZ
       RHO=0.5D0*COMZ
       IJKLANG=2*IJANG
       C2=2.0D0*DSQRT(RHO/PI)*S_P(II,JJ,Q1)**2
       IF (IJKLANG > 0) THEN
        DO M=0,IJKLANG
         ERI(M,1,1)=C2/DFLOAT(2*M+1)
        ENDDO
        DO IJ=1,NIJ
         IF (IJ == 1) THEN
          DO KL=2,NIJ
           DO M=0,IJKLANG-KLANGINDX(KL,IJT)
            ERI(M,KL,1)=QI_CI(KL,SIJ)*ERI(M,DOWN3(KL,IJT),1)
            IF (DOWN4(KL,IJT) /= 0) ERI(M,KL,1)=ERI(M,KL,1) &
           +DOWN4_C(KL,IJT)/COMZ*(ERI(M,DOWN4(KL,IJT),1)-0.5D0*ERI(M+1,DOWN4(KL,IJT),1))
            IF (DOWN5(KL,IJT) /= 0) ERI(M,KL,1)=ERI(M,KL,1) &
           +DOWN5_C(KL,IJT)/COMZ*(ERI(M,DOWN5(KL,IJT),1)-0.5D0*ERI(M+1,DOWN5(KL,IJT),1))
           ENDDO
          ENDDO
         ELSE
          DO KL=1,NIJ
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
        DO IJ=2,NIJ
         A=DSQRT(ERI(0,IJ,IJ)*MAXC(PGSI(IJ))*MAXC(PGSJ(IJ)))
         IF (A > IJIJ) IJIJ=A
        ENDDO
       ELSE
        A=DSQRT(C2*MAXC(PGSI(1))*MAXC(PGSJ(1)))
        IF (A > IJIJ) IJIJ=A
       ENDIF
       
       ! LOOP OVER KL SHELL
       DO SKL=1,NPSHELL_SQ
        IF (IJIJ*KLKL(SKL) < DOPTN(12)) CYCLE
        LARGECOMZ=COMZ+COMZ_S(SKL)
        RHO_COMZ1=COMZ_S(SKL)/LARGECOMZ
        RHO_COMZ2=COMZ/LARGECOMZ
        RHO=COMZ*RHO_COMZ1
        IJKLANG=IJANG+KLANG(SKL)
        C2=C1*DSQRT(RHO)*OV(SKL)
        KLT=KLT_S(SKL)
        ! LOOP OVER UNIT CELL INDEX 3
        IF (IJKLANG > 0) THEN
         DO Q3=0,CEL2
          QX=PX_S(SKL)+Q3X(Q3)
          QY=PY_S(SKL)*Q3C(Q3)-PZ_S(SKL)*Q3S(Q3)
          QZ=PY_S(SKL)*Q3S(Q3)+PZ_S(SKL)*Q3C(Q3)
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
          IF (NKL(KLT) /= 1) THEN
           DO KL=2,NKL(KLT)
            DO M=0,IJKLANG-KLANGINDX(KL,KLT)
             ERI(M,KL,1)=QI_CI(KL,SKL)*ERI(M,DOWN3(KL,KLT),1)+WI_QI(DOWNXYZ(KL,KLT))*ERI(M+1,DOWN3(KL,KLT),1)
             IF (DOWN4(KL,KLT) /= 0) ERI(M,KL,1)=ERI(M,KL,1) &
           +DOWN4_C(KL,KLT)/COMZ_S(SKL)*(ERI(M,DOWN4(KL,KLT),1)-RHO_COMZ2*ERI(M+1,DOWN4(KL,KLT),1))
             IF (DOWN5(KL,KLT) /= 0) ERI(M,KL,1)=ERI(M,KL,1) &
           +DOWN5_C(KL,KLT)/COMZ_S(SKL)*(ERI(M,DOWN5(KL,KLT),1)-RHO_COMZ2*ERI(M+1,DOWN5(KL,KLT),1))
            ENDDO
           ENDDO
          ENDIF
          IF (NIJ /= 1) THEN
           DO IJ=2,NIJ
            IF ((DOWN4(IJ,IJT) == 0).AND.(DOWN5(IJ,IJT) == 0)) THEN
             DO KL=1,NKL(KLT)
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
             DO KL=1,NKL(KLT)
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
             DO KL=1,NKL(KLT)
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
             DO KL=1,NKL(KLT)
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
          IF (IFLAG == 1) THEN
           DO IJ=1,NIJ
            DO KL=1,NKL(KLT)
             A_XP(PGSI(IJ),PGSJ(IJ),Q1)=A_XP(PGSI(IJ),PGSJ(IJ),Q1)+2.0D0*P_AXP(PGSK(KL,SKL),PGSL(KL,SKL),Q2)*ERI(0,KL,IJ)
            ENDDO
           ENDDO
           IF (Q3 /= 0) THEN
            DO IJ=1,NIJ
             DO KL=1,NKL(KLT)
              A_XP(PGSK(KL,SKL),PGSL(KL,SKL),Q2)= &
              A_XP(PGSK(KL,SKL),PGSL(KL,SKL),Q2)+2.0D0*P_AXP(PGSI(IJ),PGSJ(IJ),Q1)*ERI(0,KL,IJ)
             ENDDO
            ENDDO
           ENDIF
          ENDIF
          IF (LCEL1X(Q3)) THEN
           DO IJ=1,NIJ
            DO KL=1,NKL(KLT)
             A_XP(PGSI(IJ),PGSK(KL,SKL),Q3)= &
              A_XP(PGSI(IJ),PGSK(KL,SKL),Q3)-P_AXP(PGSL(KL,SKL),PGSJ(IJ),-Q2-Q3+Q1)*ERI(0,KL,IJ)*DOPTN(76)
            ENDDO
           ENDDO
          ENDIF
          IF (LCEL2X(Q3)) THEN
           DO IJ=1,NIJ
            DO KL=1,NKL(KLT)
             A_XP(PGSK(KL,SKL),PGSI(IJ),-Q3)= &
              A_XP(PGSK(KL,SKL),PGSI(IJ),-Q3)-P_AXP(PGSJ(IJ),PGSL(KL,SKL),Q2+Q3-Q1)*ERI(0,KL,IJ)*DOPTN(76)
            ENDDO
           ENDDO
          ENDIF
          IF (LTDR) THEN
           IF (IFLAG == 1) THEN
            DO IJ=1,NIJ
             DO KL=1,NKL(KLT)
              A_YP(PGSI(IJ),PGSJ(IJ),Q1)=A_YP(PGSI(IJ),PGSJ(IJ),Q1)+2.0D0*P_AYP(PGSK(KL,SKL),PGSL(KL,SKL),Q2)*ERI(0,KL,IJ)
              B_XP(PGSI(IJ),PGSJ(IJ),Q1)=B_XP(PGSI(IJ),PGSJ(IJ),Q1)+2.0D0*P_BXP(PGSK(KL,SKL),PGSL(KL,SKL),Q2)*ERI(0,KL,IJ)
              B_YP(PGSI(IJ),PGSJ(IJ),Q1)=B_YP(PGSI(IJ),PGSJ(IJ),Q1)+2.0D0*P_BYP(PGSK(KL,SKL),PGSL(KL,SKL),Q2)*ERI(0,KL,IJ)
             ENDDO
            ENDDO
            IF (Q3 /= 0) THEN
             DO IJ=1,NIJ
              DO KL=1,NKL(KLT)
               A_YP(PGSK(KL,SKL),PGSL(KL,SKL),Q2)=A_YP(PGSK(KL,SKL),PGSL(KL,SKL),Q2)+2.0D0*P_AYP(PGSI(IJ),PGSJ(IJ),Q1)*ERI(0,KL,IJ)
               B_XP(PGSK(KL,SKL),PGSL(KL,SKL),Q2)=B_XP(PGSK(KL,SKL),PGSL(KL,SKL),Q2)+2.0D0*P_BXP(PGSI(IJ),PGSJ(IJ),Q1)*ERI(0,KL,IJ)
               B_YP(PGSK(KL,SKL),PGSL(KL,SKL),Q2)=B_YP(PGSK(KL,SKL),PGSL(KL,SKL),Q2)+2.0D0*P_BYP(PGSI(IJ),PGSJ(IJ),Q1)*ERI(0,KL,IJ)
              ENDDO
             ENDDO
            ENDIF
           ENDIF
           IF (LCEL1X(Q3)) THEN
            DO IJ=1,NIJ
             DO KL=1,NKL(KLT)
              A_YP(PGSI(IJ),PGSK(KL,SKL),Q3)= &
              A_YP(PGSI(IJ),PGSK(KL,SKL),Q3)-P_AYP(PGSL(KL,SKL),PGSJ(IJ),-Q2-Q3+Q1)*ERI(0,KL,IJ)*DOPTN(76)
              B_XP(PGSI(IJ),PGSK(KL,SKL),Q3)= &
              B_XP(PGSI(IJ),PGSK(KL,SKL),Q3)-P_BXP(PGSL(KL,SKL),PGSJ(IJ),-Q2-Q3+Q1)*ERI(0,KL,IJ)*DOPTN(76)
              B_YP(PGSI(IJ),PGSK(KL,SKL),Q3)= &
              B_YP(PGSI(IJ),PGSK(KL,SKL),Q3)-P_BYP(PGSL(KL,SKL),PGSJ(IJ),-Q2-Q3+Q1)*ERI(0,KL,IJ)*DOPTN(76)
             ENDDO
            ENDDO
           ENDIF
           IF (LCEL2X(Q3)) THEN
            DO IJ=1,NIJ
             DO KL=1,NKL(KLT)
              A_YP(PGSK(KL,SKL),PGSI(IJ),-Q3)= &
              A_YP(PGSK(KL,SKL),PGSI(IJ),-Q3)-P_AYP(PGSJ(IJ),PGSL(KL,SKL),Q2+Q3-Q1)*ERI(0,KL,IJ)*DOPTN(76)
              B_XP(PGSK(KL,SKL),PGSI(IJ),-Q3)= &
              B_XP(PGSK(KL,SKL),PGSI(IJ),-Q3)-P_BXP(PGSJ(IJ),PGSL(KL,SKL),Q2+Q3-Q1)*ERI(0,KL,IJ)*DOPTN(76)
              B_YP(PGSK(KL,SKL),PGSI(IJ),-Q3)= &
              B_YP(PGSK(KL,SKL),PGSI(IJ),-Q3)-P_BYP(PGSJ(IJ),PGSL(KL,SKL),Q2+Q3-Q1)*ERI(0,KL,IJ)*DOPTN(76)
             ENDDO
            ENDDO
           ENDIF
          ENDIF
         ENDDO
        ELSE
         DO Q3=0,CEL2
          QX=PX_S(SKL)+Q3X(Q3)
          QY=PY_S(SKL)*Q3C(Q3)-PZ_S(SKL)*Q3S(Q3)
          QZ=PY_S(SKL)*Q3S(Q3)+PZ_S(SKL)*Q3C(Q3)
          T=RHO*((QX-PX)**2+(QY-PY)**2+(QZ-PZ)**2)
          IF (T < TF(0)) THEN
           TS=NINT(T*20.0D0)
           H=0.05D0*DFLOAT(TS)-T
           ERI(0,1,1)=((((((IGAMMA(TS,6)*H*0.166666666666667D0+IGAMMA(TS,5))*H*0.2D0+IGAMMA(TS,4))*H*0.25D0+ &
           IGAMMA(TS,3))*H*0.333333333333333D0+IGAMMA(TS,2))*H*0.5D0+IGAMMA(TS,1))*H+IGAMMA(TS,0))*C2
          ELSE
           ERI(0,1,1)=C2*F1(0)*DEXP(-0.5D0*DLOG(T))
          ENDIF
          IF (IFLAG == 1) THEN
           A_XP(PGSI(1),PGSJ(1),Q1)=A_XP(PGSI(1),PGSJ(1),Q1)+2.0D0*P_AXP(PGSK(1,SKL),PGSL(1,SKL),Q2)*ERI(0,1,1)
           IF (Q3 /= 0) A_XP(PGSK(1,SKL),PGSL(1,SKL),Q2)= &
              A_XP(PGSK(1,SKL),PGSL(1,SKL),Q2)+2.0D0*P_AXP(PGSI(1),PGSJ(1),Q1)*ERI(0,1,1)
          ENDIF
          ! NOTE THE ORDER OF CGS ARGUMENTS IN DENSITY MATRIX
          IF (LCEL1X(Q3)) A_XP(PGSI(1),PGSK(1,SKL),Q3)= &
              A_XP(PGSI(1),PGSK(1,SKL),Q3)-P_AXP(PGSL(1,SKL),PGSJ(1),-Q2-Q3+Q1)*ERI(0,1,1)*DOPTN(76)
          IF (LCEL2X(Q3)) A_XP(PGSK(1,SKL),PGSI(1),-Q3)= &
              A_XP(PGSK(1,SKL),PGSI(1),-Q3)-P_AXP(PGSJ(1),PGSL(1,SKL),Q2+Q3-Q1)*ERI(0,1,1)*DOPTN(76)
          IF (LTDR) THEN
           IF (IFLAG == 1) THEN
            A_YP(PGSI(1),PGSJ(1),Q1)=A_YP(PGSI(1),PGSJ(1),Q1)+2.0D0*P_AYP(PGSK(1,SKL),PGSL(1,SKL),Q2)*ERI(0,1,1)
            B_XP(PGSI(1),PGSJ(1),Q1)=B_XP(PGSI(1),PGSJ(1),Q1)+2.0D0*P_BXP(PGSK(1,SKL),PGSL(1,SKL),Q2)*ERI(0,1,1)
            B_YP(PGSI(1),PGSJ(1),Q1)=B_YP(PGSI(1),PGSJ(1),Q1)+2.0D0*P_BYP(PGSK(1,SKL),PGSL(1,SKL),Q2)*ERI(0,1,1)
            IF (Q3 /= 0) THEN
             A_YP(PGSK(1,SKL),PGSL(1,SKL),Q2)=A_YP(PGSK(1,SKL),PGSL(1,SKL),Q2)+2.0D0*P_AYP(PGSI(1),PGSJ(1),Q1)*ERI(0,1,1)
             B_XP(PGSK(1,SKL),PGSL(1,SKL),Q2)=B_XP(PGSK(1,SKL),PGSL(1,SKL),Q2)+2.0D0*P_BXP(PGSI(1),PGSJ(1),Q1)*ERI(0,1,1)
             B_YP(PGSK(1,SKL),PGSL(1,SKL),Q2)=B_YP(PGSK(1,SKL),PGSL(1,SKL),Q2)+2.0D0*P_BYP(PGSI(1),PGSJ(1),Q1)*ERI(0,1,1)
            ENDIF
           ENDIF
           ! NOTE THE ORDER OF CGS ARGUMENTS IN DENSITY MATRIX
           IF (LCEL1X(Q3)) THEN
            A_YP(PGSI(1),PGSK(1,SKL),Q3)=A_YP(PGSI(1),PGSK(1,SKL),Q3)-P_AYP(PGSL(1,SKL),PGSJ(1),-Q2-Q3+Q1)*ERI(0,1,1)*DOPTN(76)
            B_XP(PGSI(1),PGSK(1,SKL),Q3)=B_XP(PGSI(1),PGSK(1,SKL),Q3)-P_BXP(PGSL(1,SKL),PGSJ(1),-Q2-Q3+Q1)*ERI(0,1,1)*DOPTN(76)
            B_YP(PGSI(1),PGSK(1,SKL),Q3)=B_YP(PGSI(1),PGSK(1,SKL),Q3)-P_BYP(PGSL(1,SKL),PGSJ(1),-Q2-Q3+Q1)*ERI(0,1,1)*DOPTN(76)
           ENDIF
           IF (LCEL2X(Q3)) THEN
            A_YP(PGSK(1,SKL),PGSI(1),-Q3)=A_YP(PGSK(1,SKL),PGSI(1),-Q3)-P_AYP(PGSJ(1),PGSL(1,SKL),Q2+Q3-Q1)*ERI(0,1,1)*DOPTN(76)
            B_XP(PGSK(1,SKL),PGSI(1),-Q3)=B_XP(PGSK(1,SKL),PGSI(1),-Q3)-P_BXP(PGSJ(1),PGSL(1,SKL),Q2+Q3-Q1)*ERI(0,1,1)*DOPTN(76)
            B_YP(PGSK(1,SKL),PGSI(1),-Q3)=B_YP(PGSK(1,SKL),PGSI(1),-Q3)-P_BYP(PGSJ(1),PGSL(1,SKL),Q2+Q3-Q1)*ERI(0,1,1)*DOPTN(76)
           ENDIF
          ENDIF
         ENDDO
        ENDIF
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO

   DEALLOCATE(Q3X,Q3C,Q3S,COMZ_S,PX_S,PY_S,PZ_S,OV,QI_CI)
   DEALLOCATE(KLT_S,KLANG,PGSK,PGSL,KLKL)
   DEALLOCATE(LCEL1X,LCEL2X,MAXC)

!  CONTRACT COULOMB & EXCHANGE
   ALLOCATE(W1(NPGS,NPGS),W2(NCGS,NPGS),W3(NCGS,NCGS))
   DO Q1=-CEL1,CEL1
    W1=A_XP(:,:,Q1)
    W2=MATMUL(CC,W1)
    W3=MATMUL(W2,TRANSPOSE(CC))
    A_X(:,:,Q1)=W3
    IF (LTDR) THEN
     W1=A_YP(:,:,Q1)
     W2=MATMUL(CC,W1)
     W3=MATMUL(W2,TRANSPOSE(CC))
     A_Y(:,:,Q1)=W3
     W1=B_XP(:,:,Q1)
     W2=MATMUL(CC,W1)
     W3=MATMUL(W2,TRANSPOSE(CC))
     B_X(:,:,Q1)=W3
     W1=B_YP(:,:,Q1)
     W2=MATMUL(CC,W1)
     W3=MATMUL(W2,TRANSPOSE(CC))
     B_Y(:,:,Q1)=W3
    ENDIF
   ENDDO
   DEALLOCATE(W1,W2,W3)
   DEALLOCATE(A_XP,P_AXP)
   IF (LTDR) THEN
    DEALLOCATE(A_YP,P_AYP)
    DEALLOCATE(B_XP,P_BXP)
    DEALLOCATE(B_YP,P_BYP)
   ENDIF

   RETURN
END SUBROUTINE



SUBROUTINE CONTRACT_FOURINDEX_ERI(X,T,N,S,IFLAG)
! RESTORE FOUR-CENTER TWO-ELECTRON ELECTRON REPULSION INTEGRALS FROM THE INTEGRAL FILE.

   USE CONTROL
   USE BASISSET
   USE INTEGRAL
   
   IMPLICIT NONE
   INTEGER,PARAMETER :: CASHESIZE = 1000 ! CASHE SIZE MUST BE EXACTLY THE SAME IN THE STORE & RESTORE SUBROUTINES
   INTEGER :: IFLAG
   INTEGER :: N,S
   INTEGER :: EOF
   INTEGER :: Q1,Q2,Q3
   INTEGER :: I,J,K,L,P,Q,R
   INTEGER :: ICASHECOUNT,ICOUNT
   INTEGER(4),ALLOCATABLE :: ICASHE1(:),ICASHE2(:)
   DOUBLE PRECISION,ALLOCATABLE :: DCASHE(:)
   DOUBLE COMPLEX :: X(N,N,-S:S),T(N,N,-S:S)

   ALLOCATE(ICASHE1(CASHESIZE),ICASHE2(CASHESIZE),DCASHE(CASHESIZE))

   REWIND(31)
   ICOUNT=0

   ! DO NOT ZERO SCRATCH K_C
   DO
    READ(31,IOSTAT=EOF) ICASHE1,ICASHE2,DCASHE
    IF (EOF /= 0) EXIT
    DO ICASHECOUNT=1,CASHESIZE
     IF (ICASHE1(ICASHECOUNT) == -1) EXIT
     ICOUNT=ICOUNT+1
     P=0
     Q=0
     R=0
     I=0
     J=0
     K=0
     L=0
     CALL MVBITS(ICASHE1(ICASHECOUNT),16,8,P,0)
     CALL MVBITS(ICASHE1(ICASHECOUNT), 8,8,Q,0)
     CALL MVBITS(ICASHE1(ICASHECOUNT), 0,8,R,0)
     CALL MVBITS(ICASHE2(ICASHECOUNT),24,8,I,0)
     CALL MVBITS(ICASHE2(ICASHECOUNT),16,8,J,0)
     CALL MVBITS(ICASHE2(ICASHECOUNT), 8,8,K,0)
     CALL MVBITS(ICASHE2(ICASHECOUNT), 0,8,L,0)
     Q1=P-CEL1
     Q2=Q-CEL1
     Q3=R-CEL2
     IF ((Q1 < -REDUCED_CEL1).OR.(Q1 > REDUCED_CEL1).OR.(Q2 < -REDUCED_CEL1).OR.(Q2 > REDUCED_CEL1).OR.(Q3 < 0).OR.(Q3 > CEL2)) &
     CALL PABORT('INTEGRAL FILE HAS BEEN DEGRADED')
     IF (IFLAG == 1) THEN
      X(I,J,Q1)=X(I,J,Q1)+2.0D0*T(K,L,Q2)*DCASHE(ICASHECOUNT)
      IF (Q3 /= 0) X(K,L,Q2)=X(K,L,Q2)+2.0D0*T(I,J,Q1)*DCASHE(ICASHECOUNT)
     ENDIF
     ! NOTE THE ORDER OF CGS ARGUMENTS IN DENSITY MATRIX (PT_C)
     IF ((Q3 <= CEL1).AND.(-Q2-Q3+Q1 <= CEL1).AND.(-Q2-Q3+Q1 >= -CEL1)) &
               X(I,K,Q3)=X(I,K,Q3)-T(L,J,-Q2-Q3+Q1)*DCASHE(ICASHECOUNT)*DOPTN(76)
     IF ((Q3 /= 0).AND.(Q3 <= CEL1).AND.(Q2+Q3-Q1 <= CEL1).AND.(Q2+Q3-Q1 >= -CEL1)) &
               X(K,I,-Q3)=X(K,I,-Q3)-T(J,L,Q2+Q3-Q1)*DCASHE(ICASHECOUNT)*DOPTN(76)
!    THE FOLLOWING TWO LINES MAY SUBSTITUTE THE PRECEDING TWO LINES WITH INCREASED ARGUMENT RANGE OF CEL1 FOR T AND X
!    X(I,L,Q3+Q2)=X(I,L,Q3+Q2)-T(K,J,Q1-Q3)*DCASHE(ICASHECOUNT)*DOPTN(76)
!    IF (Q3 /= 0) X(K,J,Q1-Q3)=X(K,J,Q1-Q3)-T(I,L,Q2+Q3)*DCASHE(ICASHECOUNT)*DOPTN(76)
    ENDDO
   ENDDO
   IF (IOPTN(9) >= 2) WRITE(6,'(I9,A,F7.3,A)') ICOUNT, &
   ' ERIS (',DFLOAT(ICOUNT)/DFLOAT((CEL1*2+1)**2*(CEL2+1)*NCGS**4)*100.0,'% OF TOTAL ERIS) HAVE BEEN RESTORED'

   DEALLOCATE(ICASHE1,ICASHE2,DCASHE)
   
   RETURN
END SUBROUTINE
