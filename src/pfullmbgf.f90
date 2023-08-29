SUBROUTINE FULL_MBGF(OMEGA,EN0,SELFENERGY,ISELF)
! CONSTRUCT ONE-PARTICLE EXACT MBGF FOR A GIVEN BASIS SET.
! SEE EQ.(13) OF PICKUP & GOSCINSKI, MOL.PHYS. 26, 1013 (1973).

   USE CONSTANTS
   USE CONTROL
   USE GRADIENT
   USE STRUCTURE
   USE INTEGRAL
   USE BASISSET
   USE FULLCI

   IMPLICIT NONE
   INTEGER :: ISELF ! 1: SELF-ENERGY ; 2: GREEN'S FUNCTION ; 3: IP ; 4: EA ; 5 ENERGY
   DOUBLE PRECISION :: OMEGA
   REAL :: MEM,ICPUS,ICPUE
   INTEGER :: MOP,MOQ,MOR
   INTEGER :: IA,IB,IC,ID
   INTEGER(4) :: CFONE
   DOUBLE PRECISION :: SELFENERGY(IALL(0,0,0)-IVIRTCORE-ICORE,IALL(0,0,0)-IVIRTCORE-ICORE)
   DOUBLE PRECISION,ALLOCATABLE :: VEC0(:,:),VEC1(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: EIP(:),WIP(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: EEA(:),WEA(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: VIP1(:,:),VIP2(:,:),VIP3(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: VEA1(:,:),VEA2(:,:),VEA3(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: GREEN(:,:)
   DOUBLE PRECISION :: EN0,X,SGN
   DOUBLE PRECISION,ALLOCATABLE :: VL(:,:),VR(:,:),ER(:),EI(:),WK(:)
   INTEGER :: INFO
   DOUBLE PRECISION,ALLOCATABLE :: A(:,:),B(:,:),AS(:,:),BS(:,:)
   INTEGER,ALLOCATABLE :: INDX(:)
   DOUBLE PRECISION :: C

   CALL PCPU_TIME(ICPUS)
!  WRITE(6,'(A,F20.15,A)') 'ONE-PARTICLE EXACT MBGF WILL BE CONSTRUCTED AT OMEGA=',OMEGA,' HARTREE'

!  MEM=0.0
!  IF (MEM > 1000000.0) THEN
!   WRITE(6,'(A,F7.1,A)') 'ESTIMATED MEMORY USAGE WILL BE ',MEM/1000000.0,' MB'
!  ELSE IF (MEM > 1000.0) THEN
!   WRITE(6,'(A,F7.1,A)') 'ESTIMATED MEMORY USAGE WILL BE ',MEM/1000.0,' KB'
!  ELSE
!   WRITE(6,'(A,F7.1,A)') 'ESTIMATED MEMORY USAGE WILL BE ',MEM,' B'
!  ENDIF
!  IF (MEM > DOPTN(28)*1000000.0) CALL PABORT('OUT OF MEMORY') 

   ALLOCATE(VEC0(NCF,NCF),VEC1(NCF,NCF))
   ALLOCATE(EIP(NCF*IP_NCF),WIP(NCF*IP_NCF,NCF*IP_NCF))
   ALLOCATE(EEA(NCF*EA_NCF),WEA(NCF*EA_NCF,NCF*EA_NCF))
   ALLOCATE(VIP1(IP_NCF,NCF),VIP2(IP_NCF,NCF),VIP3(IP_NCF,NCF))
   ALLOCATE(VEA1(EA_NCF,NCF),VEA2(EA_NCF,NCF),VEA3(EA_NCF,NCF))
   ALLOCATE(GREEN(IALL(0,0,0)-IVIRTCORE-ICORE,IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(VL(1,IALL(0,0,0)-IVIRTCORE-ICORE),ER(IALL(0,0,0)-IVIRTCORE-ICORE),EI(IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(VR(IALL(0,0,0)-IVIRTCORE-ICORE,IALL(0,0,0)-IVIRTCORE-ICORE),WK(4*IALL(0,0,0)-IVIRTCORE-ICORE))

   ! EVALUATE E^(N)_0  
   OPEN(50,FILE=TRIM(COPTN(1))//'.gf0',FORM='UNFORMATTED')
!  WRITE(6,'(A,A)') '(N) WAVE FUNCTION FOR MBGF READ FROM ',TRIM(COPTN(1))//'.gf0'
   REWIND(50)
   READ(50) VEC0
!  FULL NORMALIZATION FROM HERE ...
   X=0.0D0
   DO IA=1,NCF
    DO IB=1,NCF
     X=X+VEC0(IB,IA)**2
    ENDDO
   ENDDO
   VEC0=VEC0/DSQRT(X)
   REWIND(50)
   WRITE(50) VEC0
!  ... TO HERE
   IF (EN0==0.0D0) THEN
    OPEN(51,FILE=TRIM(COPTN(1))//'.gfa',FORM='UNFORMATTED')
    CALL HAMILTONIAN_PRODUCT(50,51,0,2*MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE))
    REWIND(51)
    READ(51) VEC1
!!!INTERMEDIATE NORMALIZATION FROM HERE ...
!!!EN0=VEC1(1,1)
!!!... TO HERE; OR FULL NORMALIZATION FROM HERE ...
    EN0=0.0D0
    DO IA=1,NCF
     DO IB=1,NCF
      EN0=EN0+VEC0(IB,IA)*VEC1(IB,IA)
     ENDDO
    ENDDO
   ENDIF
!  WRITE(6,'(A,F20.15,A)') '(N) ENERGY = ',EN0,' HARTREE'
!  ... TO HERE

   ! READ IONIZED ENERGIES AND WAVE FUNCTIONS
   OPEN(52,FILE=TRIM(COPTN(1))//'.gf1',FORM='UNFORMATTED')
!  WRITE(6,'(A,A)') '(N-1) WAVE FUNCTIONS FOR MBGF READ FROM ',TRIM(COPTN(1))//'.gf1'
   REWIND(52) 
   READ(52) EIP
   READ(52) WIP
   CLOSE(52)
   ! ORTHONORMALIZE WAVE FUNCTIONS
   DO IC=1,NCF*IP_NCF
    DO IA=1,NCF
     DO IB=1,IP_NCF
      VIP1(IB,IA)=WIP(IB+(IA-1)*IP_NCF,IC)
     ENDDO
    ENDDO
    IF (IC > 1) THEN
     DO ID=1,IC-1
      DO IA=1,NCF
       DO IB=1,IP_NCF
        VIP2(IB,IA)=WIP(IB+(IA-1)*IP_NCF,ID)
       ENDDO
      ENDDO
      X=0.0D0
      DO IA=1,NCF
       DO IB=1,IP_NCF
        X=X+VIP1(IB,IA)*VIP2(IB,IA)
       ENDDO
      ENDDO
      DO IA=1,NCF
       DO IB=1,IP_NCF
        VIP1(IB,IA)=VIP1(IB,IA)-X*VIP2(IB,IA)
       ENDDO
      ENDDO
     ENDDO
    ENDIF
    X=0.0D0
    DO IA=1,NCF
     DO IB=1,IP_NCF
      X=X+VIP1(IB,IA)*VIP1(IB,IA)
     ENDDO
    ENDDO
    VIP1=VIP1/DSQRT(X)
    DO IA=1,NCF
     DO IB=1,IP_NCF
      WIP(IB+(IA-1)*IP_NCF,IC)=VIP1(IB,IA)
     ENDDO
    ENDDO
!   WRITE(6,'(A,I3,A,F20.10)') 'IP ',IC,' = ',EN0-EIP(IC)
   ENDDO

   ! READ ELECTRON-ATTACHED ENERGIES AND WAVE FUNCTIONS
   OPEN(52,FILE=TRIM(COPTN(1))//'.gf2',FORM='UNFORMATTED')
!  WRITE(6,'(A,A)') '(N+1) WAVE FUNCTIONS FOR MBGF READ FROM ',TRIM(COPTN(1))//'.gf2'
   REWIND(52) 
   READ(52) EEA
   READ(52) WEA
   CLOSE(52)
   ! ORTHONORMALIZE WAVE FUNCTIONS
   DO IC=1,NCF*EA_NCF
    DO IA=1,NCF
     DO IB=1,EA_NCF
      VEA1(IB,IA)=WEA(IB+(IA-1)*EA_NCF,IC)
     ENDDO
    ENDDO
    IF (IC > 1) THEN
     DO ID=1,IC-1
      DO IA=1,NCF
       DO IB=1,EA_NCF
        VEA2(IB,IA)=WEA(IB+(IA-1)*EA_NCF,ID)
       ENDDO
      ENDDO
      X=0.0D0
      DO IA=1,NCF
       DO IB=1,EA_NCF
        X=X+VEA1(IB,IA)*VEA2(IB,IA)
       ENDDO
      ENDDO
      DO IA=1,NCF
       DO IB=1,EA_NCF
        VEA1(IB,IA)=VEA1(IB,IA)-X*VEA2(IB,IA)
       ENDDO
      ENDDO
     ENDDO
    ENDIF
    X=0.0D0
    DO IA=1,NCF
     DO IB=1,EA_NCF
      X=X+VEA1(IB,IA)*VEA1(IB,IA)
     ENDDO
    ENDDO
    VEA1=VEA1/DSQRT(X)
    DO IA=1,NCF
     DO IB=1,EA_NCF
      WEA(IB+(IA-1)*EA_NCF,IC)=VEA1(IB,IA)
     ENDDO
    ENDDO
!   WRITE(6,'(A,I3,A,F20.10)') 'EA ',IC,' = ',EEA(IC)-EN0
   ENDDO

   REWIND(50)
   READ(50) VEC0

   ! CONSTRUCT G(PQ)
   GREEN=0.0D0
   DO MOQ=ICORE+1, IALL(0,0,0)-IVIRTCORE

    ! ANNIHILATE AN ELECTRON IN MOQ
    VIP1=0.0D0
    DO IB=1,NCF
     IF (BTEST(CFHALF(IB),MOQ-1)) THEN
      CFONE=IBCLR(CFHALF(IB),MOQ-1)
      IC=IP_ADDRSS(CFONE)
      SGN=1.0D0
      IF (MOQ > ICORE+1) THEN
       DO MOR=ICORE+1,MOQ-1
        IF (BTEST(CFHALF(IB),MOR-1)) SGN=-SGN
       ENDDO
      ENDIF
      DO IA=1,NCF
       VIP1(IC,IA)=VIP1(IC,IA)+SGN*VEC0(IB,IA)
      ENDDO
     ENDIF
    ENDDO

    ! OVERLAP WITH IONIZED STATES
    VIP2=0.0D0
    DO IC=1,NCF*IP_NCF
     DO IA=1,NCF
      DO IB=1,IP_NCF
       VIP3(IB,IA)=WIP(IB+(IA-1)*IP_NCF,IC)
      ENDDO
     ENDDO
     X=0.0D0
     DO IA=1,NCF
      DO IB=1,IP_NCF
       X=X+VIP1(IB,IA)*VIP3(IB,IA)
      ENDDO
     ENDDO
     VIP2=VIP2+VIP3*X/(OMEGA-EN0+EIP(IC))
!    VIP2=VIP2+VIP3*X
    ENDDO
!write(*,*) '********** No energy denom!'

    DO MOP=ICORE+1, IALL(0,0,0)-IVIRTCORE
     DO IB=1,IP_NCF
      IF (.NOT.BTEST(IP_CFHALF(IB),MOP-1)) THEN
       CFONE=IBSET(IP_CFHALF(IB),MOP-1)
       IC=ADDRSS(CFONE)
       SGN=1.0D0
       IF (MOP > ICORE+1) THEN
        DO MOR=ICORE+1,MOP-1
         IF (BTEST(IP_CFHALF(IB),MOR-1)) SGN=-SGN
        ENDDO
       ENDIF
       DO IA=1,NCF
        GREEN(MOP-ICORE,MOQ-ICORE)=GREEN(MOP-ICORE,MOQ-ICORE)+SGN*VIP2(IB,IA)*VEC0(IC,IA) ! FULL NORMALIZATION
!       IF ((IC == 1).AND.(IA == 1)) GREEN(MOP,MOQ)=GREEN(MOP,MOQ)+SGN*VIP2(IB,IA) ! INTERMEDIATE NORMALIZATION
       ENDDO
      ENDIF
     ENDDO
    ENDDO
!write(*,*) '********* EA ONLY!'

    ! CREATE AN ELECTRON IN MOQ
    VEA1=0.0D0
    DO IB=1,NCF
     IF (.NOT.BTEST(CFHALF(IB),MOQ-1)) THEN
      CFONE=IBSET(CFHALF(IB),MOQ-1)
      IC=EA_ADDRSS(CFONE)
      SGN=1.0D0
      IF (MOQ > ICORE+1) THEN
       DO MOR=ICORE+1,MOQ-1
        IF (BTEST(CFHALF(IB),MOR-1)) SGN=-SGN
       ENDDO
      ENDIF
      DO IA=1,NCF
       VEA1(IC,IA)=VEA1(IC,IA)+SGN*VEC0(IB,IA)
      ENDDO
     ENDIF
    ENDDO

    ! OVERLAP WITH ELECTRON-ATTACHED STATES
    VEA2=0.0D0
    DO IC=1,NCF*EA_NCF
     DO IA=1,NCF
      DO IB=1,EA_NCF
       VEA3(IB,IA)=WEA(IB+(IA-1)*EA_NCF,IC)
      ENDDO
     ENDDO
     X=0.0D0
     DO IA=1,NCF
      DO IB=1,EA_NCF
       X=X+VEA1(IB,IA)*VEA3(IB,IA)
      ENDDO
     ENDDO
     VEA2=VEA2+VEA3*X/(OMEGA+EN0-EEA(IC))
!    VEA2=VEA2+VEA3*X
    ENDDO
!write(*,*) '********** No energy denom!'

    DO MOP=ICORE+1, IALL(0,0,0)-IVIRTCORE
     DO IB=1,EA_NCF
      IF (BTEST(EA_CFHALF(IB),MOP-1)) THEN
       CFONE=IBCLR(EA_CFHALF(IB),MOP-1)
       IC=ADDRSS(CFONE)
       SGN=1.0D0
       IF (MOP > ICORE+1) THEN
        DO MOR=ICORE+1,MOP-1
         IF (BTEST(EA_CFHALF(IB),MOR-1)) SGN=-SGN
        ENDDO
       ENDIF
       DO IA=1,NCF
        GREEN(MOQ-ICORE,MOP-ICORE)=GREEN(MOQ-ICORE,MOP-ICORE)+SGN*VEA2(IB,IA)*VEC0(IC,IA) ! FULL NORMALIZATION
!       IF ((IC == 1).AND.(IA == 1)) GREEN(MOQ,MOP)=GREEN(MOQ,MOP)+SGN*VEA2(IB,IA)*VEC0(IC,IA) ! INTERMEDIATE NORMALIZATION
       ENDDO
      ENDIF
     ENDDO
    ENDDO
!write(*,*) '********* IP ONLY!'

   ENDDO

!  WRITE(7,'(20F15.8)') OMEGA,(GREEN(MOP-ICORE,MOP-ICORE),MOP=ICORE+1,IALL(0,0,0)-IVIRTCORE)
!  WRITE(6,'(A)') "EXACT GREEN'S FUNCTION"
!  CALL DUMP5(GREEN,IALL(0,0,0)-IVIRTCORE-ICORE)

   IF (ISELF==1) THEN

    IC=IALL(0,0,0)-IVIRTCORE-ICORE
    ALLOCATE(A(IC,IC),B(IC,IC),AS(IC,IC),BS(IC,IC),INDX(IC))
    B=0.0D0
    DO IA=1,IC
     B(IA,IA)=1.0D0
    ENDDO
    A=GREEN
    AS=A
    BS=B
    CALL LUDCMP(A,IC,IC,INDX,C)
    DO IA=1,IC
     CALL LUBKSB(A,IC,IC,INDX,B(:,IA))
     CALL MPROVE(AS,A,IC,IC,INDX,BS(:,IA),B(:,IA))
    ENDDO
!   WRITE(6,'(A)') "INVERSE OF EXACT GREEN'S FUNCTION"
!   CALL DUMP5(B,IALL(0,0,0)-IVIRTCORE-ICORE)
    BS=0.0D0
    DO IA=1,IC
     BS(IA,IA)=OMEGA-EPSILON(IA+ICORE,0,0,0)
     DO IB=1,IC
      BS(IA,IB)=BS(IA,IB)-B(IA,IB)
     ENDDO
    ENDDO   
    SELFENERGY=BS
!   WRITE(6,'(A)') "EXACT SIGMA"
!   CALL DUMP5(BS,IALL(0,0,0)-IVIRTCORE-ICORE)
    DEALLOCATE(A,B,AS,BS,INDX)
!   WRITE(6,'(20F15.8)') OMEGA,(EPSILON(MOP,0,0,0)+OMEGA-EPSILON(MOP,0,0,0)-1.0D0/GREEN(MOP-ICORE,MOP-ICORE), &
!                       MOP=ICORE+1,IALL(0,0,0)-IVIRTCORE)
!   WRITE(9,'(20F15.8)') OMEGA,(EPSILON(MOP,0,0,0)+OMEGA-EPSILON(MOP,0,0,0)-1.0D0/GREEN(MOP-ICORE,MOP-ICORE), &
!                       MOP=ICORE+1,IALL(0,0,0)-IVIRTCORE)
!   WRITE(9,'(7F15.8)') OMEGA,(1.0D0/GREEN(MOP,MOP),MOP=1,5)
 
!   CALL DGEEV('N','V',IALL(0,0,0)-IVIRTCORE-ICORE,GREEN,IALL(0,0,0)-IVIRTCORE-ICORE,ER,EI, &
!     VL,1,VR,IALL(0,0,0)-IVIRTCORE-ICORE,WK,4*IALL(0,0,0)-IVIRTCORE-ICORE,INFO)
!   IF (INFO /= 0) CALL PABORT('DGEEV FAILED TO DIAGONALIZE A MATRIX')
!   WRITE(10,'(20F15.8)') OMEGA,(ER(MOP-ICORE),MOP=ICORE+1,IALL(0,0,0)-IVIRTCORE)
!   WRITE(11,'(20F15.8)') OMEGA,(EPSILON(MOP,0,0,0)+OMEGA-EPSILON(MOP,0,0,0)-1.0D0/ER(MOP-ICORE), &
!                       MOP=ICORE+1,IALL(0,0,0)-IVIRTCORE)
   
   ELSE IF (ISELF==2) THEN

    SELFENERGY=GREEN

   ELSE IF (ISELF==3) THEN

    DO MOP=1,IALL(0,0,0)-IVIRTCORE-ICORE
     DO MOQ=1,IALL(0,0,0)-IVIRTCORE-ICORE
      SELFENERGY(MOQ,MOP)=EIP((MOP-1)*(IALL(0,0,0)-IVIRTCORE-ICORE)+MOQ)
     ENDDO
    ENDDO

   ELSE IF (ISELF==4) THEN

    DO MOP=1,IALL(0,0,0)-IVIRTCORE-ICORE
     DO MOQ=1,IALL(0,0,0)-IVIRTCORE-ICORE
      SELFENERGY(MOQ,MOP)=EEA((MOP-1)*(IALL(0,0,0)-IVIRTCORE-ICORE)+MOQ)
     ENDDO
    ENDDO

   ELSE IF (ISELF==5) THEN

    DO MOP=1,IALL(0,0,0)-IVIRTCORE-ICORE
     DO MOQ=1,IALL(0,0,0)-IVIRTCORE-ICORE
      SELFENERGY(MOQ,MOP)=EN0
     ENDDO
    ENDDO


   ENDIF

   DEALLOCATE(EIP,WIP,EEA,WEA,VEC0,VEC1,VIP1,VIP2,VIP3,VEA1,VEA2,VEA3,GREEN)
   DEALLOCATE(VL,VR,ER,EI,WK)

   RETURN

END SUBROUTINE
