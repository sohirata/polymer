SUBROUTINE TCE_DRIVER
! TCE-GENERATED ELECTRON-CORRELATION MODULE DRIVER

   USE CONTROL
   USE STRUCTURE
   USE INTEGRAL
   USE BASISSET
   USE GRADIENT
   USE THR_FULLCI

   IMPLICIT NONE
   LOGICAL,PARAMETER :: LPERMUTATION=.TRUE.
!  LOGICAL,PARAMETER :: LPERMUTATION=.FALSE.
   LOGICAL,PARAMETER :: LDISCONNECTED=.TRUE.
!  LOGICAL,PARAMETER :: LDISCONNECTED=.FALSE.
   INTEGER,PARAMETER :: MAXITER=100
   INTEGER :: ITER
   INTEGER :: NA,NO
   INTEGER :: I,J,K,L,M,N
   LOGICAL :: LSINGLES,LDOUBLES,LTRIPLES
   LOGICAL :: LDIIS
   INTEGER :: IDIIS,IORDER
   INTEGER,ALLOCATABLE :: MAPA(:),MAPB(:)
   DOUBLE PRECISION :: EMP2,ECORR,R1,R2,R3
   DOUBLE PRECISION,ALLOCATABLE :: HCORE1(:,:),EPS1(:)
   DOUBLE PRECISION,ALLOCATABLE :: F1(:,:),F1E(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: V2(:,:,:,:),V2E(:,:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: T1(:),T1E(:),T2(:),T2E(:),T3(:),T3E(:)
   DOUBLE PRECISION,ALLOCATABLE :: I0_0(:),I0E_0(:),I0_1(:),I0E_1(:),I0_2(:),I0E_2(:),I0_3(:),I0E_3(:)
   DOUBLE PRECISION,ALLOCATABLE :: S1(:),S1E(:),S2(:),S2E(:),S3(:),S3E(:)
   DOUBLE PRECISION,ALLOCATABLE :: I1(:),I1E(:),I2(:),I2E(:),I3(:),I3E(:)
   CHARACTER(3)  :: APPEND
   DOUBLE PRECISION,ALLOCATABLE :: B(:,:),BS(:,:),C(:),CS(:)
   INTEGER,ALLOCATABLE :: INDX(:)
   DOUBLE PRECISION :: D
   INTEGER :: IA,IB
   DOUBLE PRECISION,ALLOCATABLE :: E0(:),EHF
   INTEGER,ALLOCATABLE :: ITARGETS(:)
   DOUBLE PRECISION :: PREVIOUSMIN,CURRENTMIN
   INTEGER :: ISTATE,JSTATE,NDEGEN
   INTEGER,ALLOCATABLE :: MAP(:,:),SPINMAP(:,:),BACKMAP(:,:)
   CHARACTER(40) :: CBITSA,CBITSB
   DOUBLE PRECISION,ALLOCATABLE :: D_F1(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: D_T1(:,:),D_T1E(:,:),D_T2(:,:),D_T2E(:,:),D_T3(:,:),D_T3E(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: D_I0_0(:,:),D_I0E_0(:,:),D_I0_1(:,:),D_I0E_1(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: D_I0_2(:,:),D_I0E_2(:,:),D_I0_3(:,:),D_I0E_3(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: D_S1(:,:),D_S1E(:,:),D_S2(:,:),D_S2E(:,:),D_S3(:,:),D_S3E(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: HMAT(:,:),EMAT(:,:),SMAT(:,:),TMP(:,:)
   INTEGER :: H1,H2,P1,P2
   INTEGER,ALLOCATABLE :: MAP2(:)
   LOGICAL :: LEXIST
   INTEGER,ALLOCATABLE :: INTERNAL1(:,:),INTERNAL2(:,:),INTERNAL3(:),INTERNAL4(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: ER(:),EI(:),VR(:,:),VL(:,:),WK(:)
   DOUBLE PRECISION,ALLOCATABLE :: DEVSQ1(:),DEVSQ2(:),DEVSQ3(:)
   INTEGER :: INFO
   CHARACTER(2),PARAMETER :: SPINNAME(2) = (/'A','B'/)
   INTEGER :: IPRTY,JPRTY,KPRTY

   WRITE(6,'(A)') '*************************'
   WRITE(6,'(A)') 'TENSOR CONTRACTION ENGINE'
   WRITE(6,'(A)') '*************************'
   IF (LPERMUTATION) THEN
    WRITE(6,'(A)') 'INDEX PERMUTATIONS ARE ON'
   ELSE
    WRITE(6,'(A)') 'INDEX PERMUTATIONS ARE OFF'
   ENDIF
   IF (LDISCONNECTED) THEN
    WRITE(6,'(A)') 'DISCONNECTED DIAGRAMS ARE ON'
   ELSE
    WRITE(6,'(A)') 'DISCONNECTED DIAGRAMS ARE OFF'
   ENDIF

   IF (IOPTN(27) /= 0) CALL WARNING('FROZEN CORE APPROXIMATION WILL BE USED')
   IF (IOPTN(55) /= 0) CALL WARNING('FROZEN VIRTUAL APPROXIMATION WILL BE USED')
   IF (IOPTN(54) > 0) THEN
    LDIIS=.TRUE.
    IDIIS=IOPTN(54)
    WRITE(6,'(A,I3,A)') 'DIIS EXTRAPOLATION WILL BE PERFORMED ONCE IN',IDIIS,' ITERATIONS'
   ELSE
    LDIIS=.FALSE.
    IDIIS=1
   ENDIF

   NO=2*(IOCC-ICORE)
   NA=2*(IALL(0,0,0)-IVIRTCORE-ICORE)

   IF ((KVCX /= 0).OR.(KVCY /= 0).OR.(KVCZ /= 0)) CALL PABORT('FOR MOLECULES ONLY')
   IF (LOPTN(25)) CALL PABORT('USE NON-DIRECT ALGORITHM')
   IF (NO <= 0) CALL PABORT('NO OCCUPIED ORBITAL')
   IF (NA-NO <= 0) CALL PABORT('NO VIRTUAL ORBITAL')

   ALLOCATE(H(IALLMAX-IVIRTCORE,IALLMAX-IVIRTCORE))
   ALLOCATE(G(IALLMAX-IVIRTCORE,IALLMAX-IVIRTCORE,IALLMAX-IVIRTCORE,IALLMAX-IVIRTCORE))
   ALLOCATE(MAPA(IALLMAX-IVIRTCORE),MAPB(IALLMAX-IVIRTCORE))
   ALLOCATE(HCORE1(NA,NA),EPS1(NA))
   ALLOCATE(F1(NA,NA),V2(NA,NA,NA,NA))
   ALLOCATE(F1E(NA,NA),V2E(NA,NA,NA,NA))
   MAPA=-1
   MAPB=-1
   DO I=ICORE+1,IOCC
    MAPA(I)=I-ICORE
    MAPB(I)=NO/2+I-ICORE
   ENDDO
   DO I=IOCC+1,IALLMAX-IVIRTCORE
    MAPA(I)=NO+I-IOCC
    MAPB(I)=NO+(NA-NO)/2+I-IOCC
   ENDDO
!  DO I=1,IALLMAX-IVIRTCORE
!   WRITE(*,*) I,MAPA(I),MAPB(I)
!  ENDDO
!  CALL ONE_ELECTRON_INTEGRALS
!  CALL TWO_ELECTRON_INTEGRALS
   CALL THERMAL_ONE_ELECTRON_INTEGRALS
   CALL THERMAL_TWO_ELECTRON_INTEGRALS
   EPS1=0.0D0
   F1=0.0D0
   DO I=ICORE+1,IALLMAX-IVIRTCORE
    EPS1(MAPA(I))=EPSILON(I,0,0,0)
    EPS1(MAPB(I))=EPSILON(I,0,0,0)
    F1(MAPA(I),MAPA(I))=EPSILON(I,0,0,0) ! THIS IS CORRECT ONLY FOR GROUND-STATE REF; OVERWRITTEN FOR DEGENERATE REF
    F1(MAPB(I),MAPB(I))=EPSILON(I,0,0,0) ! THIS IS CORRECT ONLY FOR GROUND-STATE REF; OVERWRITTEN FOR DEGENERATE REF
   ENDDO
   HCORE1=0.0D0
   DO I=ICORE+1,IALLMAX-IVIRTCORE
    DO J=ICORE+1,IALLMAX-IVIRTCORE
     HCORE1(MAPA(I),MAPA(J))=H(I,J)
     HCORE1(MAPB(I),MAPB(J))=H(I,J)
    ENDDO
   ENDDO
   V2=0.0D0
   DO I=ICORE+1,IALLMAX-IVIRTCORE
    DO J=ICORE+1,IALLMAX-IVIRTCORE
     DO K=ICORE+1,IALLMAX-IVIRTCORE
      DO L=ICORE+1,IALLMAX-IVIRTCORE
       V2(MAPA(I),MAPA(J),MAPA(K),MAPA(L))=G(I,K,J,L)-G(I,L,J,K) ! <alpha,alpha||alpha,alpha>
       V2(MAPA(I),MAPB(J),MAPA(K),MAPB(L))=G(I,K,J,L)            ! <alpha,beta ||alpha,beta >
       V2(MAPA(I),MAPB(J),MAPB(K),MAPA(L))=          -G(I,L,J,K) ! <alpha,beta ||beta, alpha>
       V2(MAPB(I),MAPA(J),MAPB(K),MAPA(L))=G(I,K,J,L)            ! <beta, alpha||beta, alpha>
       V2(MAPB(I),MAPA(J),MAPA(K),MAPB(L))=          -G(I,L,J,K) ! <beta, alpha||alpha,beta >
       V2(MAPB(I),MAPB(J),MAPB(K),MAPB(L))=G(I,K,J,L)-G(I,L,J,K) ! <beta, beta ||beta, beta >
      ENDDO
     ENDDO
    ENDDO
   ENDDO

   EMP2=0.0D0
   DO I=1,NO
    DO J=1,NO
     IF (I<=J) CYCLE
     DO K=NO+1,NA
      DO L=NO+1,NA
       IF (K<=L) CYCLE
       EMP2=EMP2+V2(I,J,K,L)*V2(K,L,I,J)/(EPS1(I)+EPS1(J)-EPS1(K)-EPS1(L))
      ENDDO
     ENDDO
    ENDDO
   ENDDO
   WRITE(6,'(A,F20.15,A)') 'MP2 ENERGY = ',EMP2,' HARTREE'

   IF (.NOT.LOPTN(109)) THEN

!  -----------------
!  NON-DEGENERATE CC
!  -----------------

   DO N=1,1
    WRITE(APPEND,'(I3)') N
    OPEN(100+N,FILE=TRIM(COPTN(1))//'.T.'//APPEND,FORM='UNFORMATTED')
    OPEN(400+N,FILE=TRIM(COPTN(1))//'.delta.'//APPEND,FORM='UNFORMATTED')
   ENDDO
   IF (COPTN(110)=='CCS') THEN
    WRITE(6,'(A)') '-----------------------------------------------------'
    WRITE(6,'(A)') 'COUPLED-CLUSTER SINGLES (CCS) AKA HARTREE-FOCK BY TCE'
    WRITE(6,'(A)') '-----------------------------------------------------'
    LSINGLES=.TRUE.
    LDOUBLES=.FALSE.
    LTRIPLES=.FALSE.
   ELSE IF (COPTN(110)=='CCD') THEN
    WRITE(6,'(A)') '------------------------------------'
    WRITE(6,'(A)') 'COUPLED-CLUSTER DOUBLES (CCD) BY TCE'
    WRITE(6,'(A)') '------------------------------------'
    LSINGLES=.FALSE.
    LDOUBLES=.TRUE.
    LTRIPLES=.FALSE.
   ELSE IF (COPTN(110)=='CCSD') THEN
    WRITE(6,'(A)') '-----------------------------------------------'
    WRITE(6,'(A)') 'COUPLED-CLUSTER SINGLES & DOUBLES (CCSD) BY TCE'
    WRITE(6,'(A)') '-----------------------------------------------'
    LSINGLES=.TRUE.
    LDOUBLES=.TRUE.
    LTRIPLES=.FALSE.
   ELSE IF (COPTN(110)=='CCSDT') THEN
    WRITE(6,'(A)') '---------------------------------------------------------'
    WRITE(6,'(A)') 'COUPLED-CLUSTER SINGLES, DOUBLES & TRIPLES (CCSDT) BY TCE'
    WRITE(6,'(A)') '---------------------------------------------------------'
    LSINGLES=.TRUE.
    LDOUBLES=.TRUE.
    LTRIPLES=.TRUE.
   ELSE
    WRITE(6,'(A)') 'ILLEGAL TCE MODULE NAME ; EXITING'
    RETURN
   ENDIF
   WRITE(6,'(/,A)') '  ITER       RESIDUAL    CORRELATION        TOTAL ENERGY'
   WRITE(6,'(A)')   '--------------------------------------------------------'

   ALLOCATE(I0_0(1),I0E_0(1))
   IF (LSINGLES) ALLOCATE(I0_1(NA**2),T1(NA**2))
   IF (LDOUBLES) ALLOCATE(I0_2(NA**4),T2(NA**4))
   IF (LTRIPLES) ALLOCATE(I0_3(NA**6),T3(NA**6))
   IF (LSINGLES) ALLOCATE(I0E_1(NA**2),T1E(NA**2))
   IF (LDOUBLES) ALLOCATE(I0E_2(NA**4),T2E(NA**4))
   IF (LTRIPLES) ALLOCATE(I0E_3(NA**6),T3E(NA**6))
   IF (LDISCONNECTED.AND.LSINGLES) ALLOCATE(S1(NA**2),S1E(NA**2))
   IF (LDISCONNECTED.AND.LSINGLES.AND.LDOUBLES) ALLOCATE(S2(NA**4),S2E(NA**4))
   IF (LDISCONNECTED.AND.LSINGLES.AND.LDOUBLES.AND.LTRIPLES) ALLOCATE(S3(NA**6),S3E(NA**6))

   IF (LSINGLES) T1=0.0D0
   IF (LDOUBLES) T2=0.0D0
   IF (LTRIPLES) T3=0.0D0

   DO ITER=1,MAXITER

    REWIND(100+1)
    IF (ITER > 1) THEN
     DO I=1,ITER-1
      IF (LSINGLES) READ(100+1) t1e
      IF (LDOUBLES) READ(100+1) t2e
      IF (LTRIPLES) READ(100+1) t3e
     ENDDO
    ENDIF
    IF (LSINGLES) WRITE(100+1) t1
    IF (LDOUBLES) WRITE(100+1) t2
    IF (LTRIPLES) WRITE(100+1) t3

    IF (COPTN(110)=='CCS') THEN
     ALLOCATE(I1(NA**2),I1E(NA**2),I2(NA**2),I2E(NA**2))
     IF (LDISCONNECTED) THEN
      IF (LPERMUTATION) THEN
       CALL p_ccs_e(NA,f1,f1e,i0_0,i0e_0,i1,i1e,NO,t1,t1e,v2,v2e)
       CALL p_ccs_t1_disconnected(NA,f1,f1e,i0_1,i0e_1,i1,i1e,i2,i2e,NO,t1,t1e,v2,v2e)
      ELSE
       CALL np_ccs_e(NA,f1,i0_0,i1,NO,t1,v2)
       CALL np_ccs_t1_disconnected(NA,f1,i0_1,i1,i2,NO,t1,v2)
      ENDIF
      I0_1=I0_1-I0_0(1)*T1
     ELSE
      IF (LPERMUTATION) THEN
       CALL p_ccs_e(NA,f1,f1e,i0_0,i0e_0,i1,i1e,NO,t1,t1e,v2,v2e)
       CALL p_ccs_t1(NA,f1,f1e,i0_1,i0e_1,i1,i1e,i2,i2e,NO,t1,t1e,v2,v2e)
      ELSE
       CALL np_ccs_e(NA,f1,i0_0,i1,NO,t1,v2)
       CALL np_ccs_t1(NA,f1,i0_1,i1,i2,NO,t1,v2)
      ENDIF
     ENDIF
     DEALLOCATE(I1,I1E,I2,I2E)

    ELSE IF (COPTN(110)=='CCD') THEN
     ALLOCATE(I1(NA**4),I1E(NA**4))
     IF (LDISCONNECTED) THEN
      IF (LPERMUTATION) THEN
       CALL p_ccd_e(NA,i0_0,i0e_0,NO,t2,t2e,v2,v2e)
       CALL p_ccd_t2_disconnected(NA,f1,f1e,i0_2,i0e_2,i1,i1e,NO,t2,t2e,v2,v2e)
      ELSE
       CALL np_ccd_e(NA,i0_0,NO,t2,v2)
       CALL np_ccd_t2_disconnected(NA,f1,i0_2,i1,NO,t2,v2)
      ENDIF
      I0_2=I0_2-I0_0(1)*T2
     ELSE
      IF (LPERMUTATION) THEN
       CALL p_ccd_e(NA,i0_0,i0e_0,NO,t2,t2e,v2,v2e)
       CALL p_ccd_t2(NA,f1,f1e,i0_2,i0e_2,i1,i1e,NO,t2,t2e,v2,v2e)
      ELSE 
       CALL np_ccd_e(NA,i0_0,NO,t2,v2)
       CALL np_ccd_t2(NA,f1,i0_2,i1,NO,t2,v2)
      ENDIF
     ENDIF
     DEALLOCATE(I1,I1E)

    ELSE IF (COPTN(110)=='CCSD') THEN
     ALLOCATE(I1(NA**4),I1E(NA**4),I2(NA**4),I2E(NA**4),I3(NA**4),I3E(NA**4))
     IF (LDISCONNECTED) THEN
      IF (LPERMUTATION) THEN
       CALL p_ccsd_e(NA,f1,f1e,i0_0,i0e_0,i1,i1e,NO,t1,t1e,t2,t2e,v2,v2e)
       CALL p_ccsd_t1_disconnected(NA,f1,f1e,i0_1,i0e_1,i1,i1e,i2,i2e,NO,t1,t1e,t2,t2e,v2,v2e)
       CALL p_ccsd_t2_disconnected(NA,f1,f1e,i0_2,i0e_2,i1,i1e,i2,i2e,i3,i3e,NO,t1,t1e,t2,t2e,v2,v2e)
       CALL p_ccsd_t2_overlap(NA,s2,s2e,NO,t1,t1e,t2,t2e)
       S1=T1
      ELSE
       CALL np_ccsd_e(NA,f1,i0_0,i1,NO,t1,t2,v2)
       CALL np_ccsd_t1_disconnected(NA,f1,i0_1,i1,i2,NO,t1,t2,v2)
       CALL np_ccsd_t2_disconnected(NA,f1,i0_2,i1,i2,i3,NO,t1,t2,v2)
       CALL np_ccsd_t2_overlap(NA,s2,NO,t1,t2)
       S1=T1
      ENDIF
      I0_1=I0_1-I0_0(1)*S1
      I0_2=I0_2-I0_0(1)*S2
     ELSE
      IF (LPERMUTATION) THEN
       CALL p_ccsd_e(NA,f1,f1e,i0_0,i0e_0,i1,i1e,NO,t1,t1e,t2,t2e,v2,v2e)
       CALL p_ccsd_t1(NA,f1,f1e,i0_1,i0e_1,i1,i1e,i2,i2e,NO,t1,t1e,t2,t2e,v2,v2e)
       CALL p_ccsd_t2(NA,f1,f1e,i0_2,i0e_2,i1,i1e,i2,i2e,i3,i3e,NO,t1,t1e,t2,t2e,v2,v2e)
      ELSE
       CALL np_ccsd_e(NA,f1,i0_0,i1,NO,t1,t2,v2)
       CALL np_ccsd_t1(NA,f1,i0_1,i1,i2,NO,t1,t2,v2)
       CALL np_ccsd_t2(NA,f1,i0_2,i1,i2,i3,NO,t1,t2,v2)
      ENDIF
     ENDIF
     DEALLOCATE(I1,I1E,I2,I2E,I3,I3E)

    ELSE IF (COPTN(110)=='CCSDT') THEN
     ALLOCATE(I1(NA**6),I1E(NA**6),I2(NA**6),I2E(NA**6),I3(NA**6),I3E(NA**6))
     IF (LDISCONNECTED) THEN
      IF (LPERMUTATION) THEN
       CALL p_ccsdt_e(NA,f1,f1e,i0_0,i0e_0,i1,i1e,NO,t1,t1e,t2,t2e,v2,v2e)
       CALL p_ccsdt_t1_disconnected(NA,f1,f1e,i0_1,i0e_1,i1,i1e,i2,i2e,NO,t1,t1e,t2,t2e,t3,t3e,v2,v2e) 
       CALL p_ccsdt_t2_disconnected(NA,f1,f1e,i0_2,i0e_2,i1,i1e,i2,i2e,i3,i3e,NO,t1,t1e,t2,t2e,t3,t3e,v2,v2e)
       CALL p_ccsdt_t3_disconnected(NA,f1,f1e,i0_3,i0e_3,i1,i1e,i2,i2e,i3,i3e,NO,t1,t1e,t2,t2e,t3,t3e,v2,v2e)
       CALL p_ccsdt_t3_overlap(NA,s3,s3e,i1,i1e,NO,t1,t1e,t2,t2e,t3,t3e)
       CALL p_ccsd_t2_overlap(NA,s2,s2e,NO,t1,t1e,t2,t2e)
       S1=T1
      ELSE
       CALL np_ccsdt_e(NA,f1,i0_0,i1,NO,t1,t2,v2)
       CALL np_ccsdt_t1_disconnected(NA,f1,i0_1,i1,i2,NO,t1,t2,t3,v2) 
       CALL np_ccsdt_t2_disconnected(NA,f1,i0_2,i1,i2,i3,NO,t1,t2,t3,v2)
       CALL np_ccsdt_t3_disconnected(NA,f1,i0_3,i1,i2,i3,NO,t1,t2,t3,v2)
       CALL np_ccsdt_t3_overlap(NA,s3,i1,NO,t1,t2,t3)
       CALL np_ccsd_t2_overlap(NA,s2,NO,t1,t2)
       S1=T1
      ENDIF
      I0_1=I0_1-I0_0(1)*S1
      I0_2=I0_2-I0_0(1)*S2
      I0_3=I0_3-I0_0(1)*S3
     ELSE
      IF (LPERMUTATION) THEN
       CALL p_ccsdt_e(NA,f1,f1e,i0_0,i0e_0,i1,i1e,NO,t1,t1e,t2,t2e,v2,v2e)
       CALL p_ccsdt_t1(NA,f1,f1e,i0_1,i0e_1,i1,i1e,i2,i2e,NO,t1,t1e,t2,t2e,t3,t3e,v2,v2e) 
       CALL p_ccsdt_t2(NA,f1,f1e,i0_2,i0e_2,i1,i1e,i2,i2e,i3,i3e,NO,t1,t1e,t2,t2e,t3,t3e,v2,v2e)
       CALL p_ccsdt_t3(NA,f1,f1e,i0_3,i0e_3,i1,i1e,i2,i2e,i3,i3e,NO,t1,t1e,t2,t2e,t3,t3e,v2,v2e)
      ELSE
       CALL np_ccsdt_e(NA,f1,i0_0,i1,NO,t1,t2,v2)
       CALL np_ccsdt_t1(NA,f1,i0_1,i1,i2,NO,t1,t2,t3,v2) 
       CALL np_ccsdt_t2(NA,f1,i0_2,i1,i2,i3,NO,t1,t2,t3,v2)
       CALL np_ccsdt_t3(NA,f1,i0_3,i1,i2,i3,NO,t1,t2,t3,v2)
      ENDIF
     ENDIF
     DEALLOCATE(I1,I1E,I2,I2E,I3,I3E)
    ENDIF

    REWIND(400+1)
    IF (ITER > 1) THEN
     DO I=1,ITER-1
      IF (LSINGLES) READ(400+1) i0e_1
      IF (LDOUBLES) READ(400+1) i0e_2
      IF (LTRIPLES) READ(400+1) i0e_3
     ENDDO
    ENDIF
    IF (LSINGLES) WRITE(400+1) i0_1
    IF (LDOUBLES) WRITE(400+1) i0_2
    IF (LTRIPLES) WRITE(400+1) i0_3

    R1=0.0D0
    R2=0.0D0
    R3=0.0D0
    ! SIMPLE RELAXATION
    IF (LSINGLES) CALL T_RELAX(1,NA,NO,I0_1,T1,F1,R1,LPERMUTATION)
    IF (LDOUBLES) CALL T_RELAX(2,NA,NO,I0_2,T2,F1,R2,LPERMUTATION)
    IF (LTRIPLES) CALL T_RELAX(3,NA,NO,I0_3,T3,F1,R3,LPERMUTATION)
    IF (LSINGLES) WRITE(100+1) t1
    IF (LDOUBLES) WRITE(100+1) t2
    IF (LTRIPLES) WRITE(100+1) t3

    ECORR=I0_0(1)
    WRITE(6,'(I6,2F15.10,F20.10)') ITER,DSQRT(R1+R2+R3),ECORR,EHFKS+ECORR
    WRITE(6,'(A)')   '........................................................'
    IF (DSQRT(R1+R2+R3) < DOPTN(62)) EXIT

    ! DIIS
    IF ((LDIIS).AND.(MOD(ITER,IDIIS) == 0)) THEN
     ALLOCATE(INDX(21),B(21,21),BS(21,21),C(21),CS(21))
     IORDER=MIN(ITER,20)
     B=0.0D0
     DO I=1,IORDER
      REWIND(400+1)
      IF (ITER > I) THEN
       DO K=1,ITER-I
        IF (LSINGLES) READ(400+1) i0e_1
        IF (LDOUBLES) READ(400+1) i0e_2
        IF (LTRIPLES) READ(400+1) i0e_3
       ENDDO
      ENDIF
      IF (LSINGLES) READ(400+1) i0_1
      IF (LDOUBLES) READ(400+1) i0_2
      IF (LTRIPLES) READ(400+1) i0_3
      DO J=1,IORDER
       REWIND(400+1)
       IF (ITER > J) THEN
        DO K=1,ITER-J
         IF (LSINGLES) READ(400+1) i0e_1
         IF (LDOUBLES) READ(400+1) i0e_2
         IF (LTRIPLES) READ(400+1) i0e_3
        ENDDO
       ENDIF
       R1=0.0D0
       R2=0.0D0
       R3=0.0D0
       IF (LSINGLES) READ(400+1) i0e_1
       IF (LDOUBLES) READ(400+1) i0e_2
       IF (LTRIPLES) READ(400+1) i0e_3
       IF (LSINGLES) CALL T_OVERLAP(1,NA,NO,i0_1,i0e_1,R1,LPERMUTATION)
       IF (LDOUBLES) CALL T_OVERLAP(2,NA,NO,i0_2,i0e_2,R2,LPERMUTATION)
       IF (LTRIPLES) CALL T_OVERLAP(3,NA,NO,i0_3,i0e_3,R3,LPERMUTATION)
       B(J,I)=R1+R2+R3
      ENDDO
      B(I,IORDER+1)=-1.0D0
      B(IORDER+1,I)=-1.0D0
      C(I)=0.0D0
     ENDDO
     B(IORDER+1,IORDER+1)=0.0D0
     C(IORDER+1)=-1.0D0
     BS=B
     CS=C
!    WRITE(6,'(A)') 'PULAY MATRIX'
!    CALL DUMP5(B,21)
     CALL LUDCMP(B,IORDER+1,21,INDX,D)
     CALL LUBKSB(B,IORDER+1,21,INDX,C)
     CALL MPROVE(BS,B,IORDER+1,21,INDX,CS,C)
     WRITE(6,'(A,100F7.4:)') 'DIIS VECTOR',(C(I),I=1,IORDER+1)
     IF (LSINGLES) t1=0.0D0
     IF (LDOUBLES) t2=0.0D0
     IF (LTRIPLES) t3=0.0D0
     REWIND(100+1)
     IF (ITER > IORDER) THEN
      DO I=1,ITER-IORDER
       IF (LSINGLES) READ(100+1) t1e
       IF (LDOUBLES) READ(100+1) t2e
       IF (LTRIPLES) READ(100+1) t3e
      ENDDO
     ENDIF
     DO I=1,IORDER
      IF (LSINGLES) READ(100+1) t1e
      IF (LDOUBLES) READ(100+1) t2e
      IF (LTRIPLES) READ(100+1) t3e
      IF (LSINGLES) t1=t1+t1e*C(IORDER-I+1)
      IF (LDOUBLES) t2=t2+t2e*C(IORDER-I+1)
      IF (LTRIPLES) t3=t3+t3e*C(IORDER-I+1)
     ENDDO
     REWIND(100+1)
     DO I=1,ITER
      IF (LSINGLES) READ(100+1) t1e
      IF (LDOUBLES) READ(100+1) t2e
      IF (LTRIPLES) READ(100+1) t3e
     ENDDO
     IF (LSINGLES) WRITE(100+1) t1
     IF (LDOUBLES) WRITE(100+1) t2
     IF (LTRIPLES) WRITE(100+1) t3
     DEALLOCATE(INDX,B,BS,C,CS)
    ENDIF
   ENDDO
   IF (ITER>=MAXITER) CALL PABORT('MAXITER REACHED')
   WRITE(6,'(A)') '--------------------------------------------------------'
   WRITE(6,'(A)') 'CC ITERATION CONVERGED'          
 
   DEALLOCATE(I0_0,I0E_0)
   IF (LSINGLES) DEALLOCATE(I0_1,I0E_1,T1,T1E)
   IF (LDOUBLES) DEALLOCATE(I0_2,I0E_2,T2,T2E)
   IF (LTRIPLES) DEALLOCATE(I0_3,I0E_3,T3,T3E)
   IF (LDISCONNECTED.AND.LSINGLES) DEALLOCATE(S1,S1E)
   IF (LDISCONNECTED.AND.LSINGLES.AND.LDOUBLES) DEALLOCATE(S2,S2E)
   IF (LDISCONNECTED.AND.LSINGLES.AND.LDOUBLES.AND.LTRIPLES) DEALLOCATE(S3,S3E)

   DO N=1,1
    CLOSE(100+N)
    CLOSE(400+N)
   ENDDO

   ELSE

!  -------------
!  DEGENERATE CC
!  -------------

   IF (.NOT.LPERMUTATION)  CALL PABORT('DEGENERATE CC MUST USE PERMUTATION SYMMETRY')
   ! Because T(ab,ij) and T(ba,ji) are always degenerate, confusing the algorithms as they appear internal excitations.
   IF (.NOT.LDISCONNECTED) CALL PABORT('DEGENERATE CC EQUATIONS ARE DISCONNECTED')
   ! Diagonal disconnected diagrams can be deleted, but off-diagonal ones remain (they are considered connected),
   ! and hence the "connected" formulation does not bring much benefit in this case.

   IF (COPTN(110)=='CCS') THEN
    WRITE(6,'(A)') '-----------------------------------------------------------------'
    WRITE(6,'(A)') 'DEGENERATE COUPLED-CLUSTER SINGLES (CCS) AKA DEGENERATE HF BY TCE'
    WRITE(6,'(A)') '-----------------------------------------------------------------'
    LSINGLES=.TRUE.
    LDOUBLES=.FALSE.
    LTRIPLES=.FALSE.
   ELSE IF (COPTN(110)=='CCSD') THEN
    WRITE(6,'(A)') '----------------------------------------------------------'
    WRITE(6,'(A)') 'DEGENERATE COUPLED-CLUSTER SINGLES & DOUBLES (CCSD) BY TCE'
    WRITE(6,'(A)') '----------------------------------------------------------'
    LSINGLES=.TRUE.
    LDOUBLES=.TRUE.
    LTRIPLES=.FALSE.
   ELSE IF (COPTN(110)=='CCSDT') THEN
    CALL PABORT('DEGENERATE CCSDT NOT IMPLEMENTED')
    WRITE(6,'(A)') '--------------------------------------------------------------------'
    WRITE(6,'(A)') 'DEGENERATE COUPLED-CLUSTER SINGLES, DOUBLES & TRIPLES (CCSDT) BY TCE'
    WRITE(6,'(A)') '--------------------------------------------------------------------'
    LSINGLES=.TRUE.
    LDOUBLES=.TRUE.
    LTRIPLES=.TRUE.
   ELSE
    WRITE(6,'(A)') 'ILLEGAL TCE MODULE NAME ; EXITING'
    RETURN
   ENDIF

   IF (LOPTN(60)) THEN
    NELEA=IOCC-1
    NELEB=IOCC
    WRITE(6,'(A)') 'IP ROOTS ARE SOUGHT'
    WRITE(6,'(A,2I3)') 'NELEA,NELEB = ',NELEA,NELEB
    WRITE(6,'(A,I3)') 'NUMBER OF ROOTS = ',IOPTN(59)
   ELSE IF (LOPTN(61)) THEN
    NELEA=IOCC+1
    NELEB=IOCC
    WRITE(6,'(A)') 'EA ROOTS ARE SOUGHT'
    WRITE(6,'(A,2I3)') 'NELEA,NELEB = ',NELEA,NELEB
    WRITE(6,'(A,I3)') 'NUMBER OF ROOTS = ',IOPTN(59)
   ELSE
    NELEA=IOCC
    NELEB=IOCC
    WRITE(6,'(A)') 'EXCITED ROOTS ARE SOUGHT'
    WRITE(6,'(A,2I3)') 'NELEA,NELEB = ',NELEA,NELEB
    WRITE(6,'(A,I3)') 'NUMBER OF ROOTS = ',IOPTN(59)
   ENDIF

   NO=NELEA+NELEB-2*ICORE
   NA=2*(IALL(0,0,0)-IVIRTCORE-ICORE)
   WRITE(6,'(A,2I5)') 'NUMBERS OF ALL AND OCCUPIED SPIN-ORBITALS = ',NA,NO
   !
   ! ZEROTH-ORDER ENERGIES
   !
   CALL THERMAL_GENERATE_CONFIGURATIONSA(NELEA,MIN(NELEA,IALL(0,0,0)-IVIRTCORE-NELEA))
   CALL THERMAL_GENERATE_CONFIGURATIONSB(NELEB,MIN(NELEB,IALL(0,0,0)-IVIRTCORE-NELEB))
   NALL=NCFA*NCFB
   WRITE(6,'(A,I10)') 'DIMENSION OF E0 ARRAY = ',NALL
   ALLOCATE(E0(NALL))
   DO IA=1,NCFA
    DO IB=1,NCFB
     E0((IA-1)*NCFB+IB)=NUCLEAR_REPULSION
     DO I=1,IALL(0,0,0)-IVIRTCORE
      IF (BTEST(CFHALFA(IA),I-1)) &
       E0((IA-1)*NCFB+IB)=E0((IA-1)*NCFB+IB)+EPSILON(I,0,0,0)
      IF (BTEST(CFHALFB(IB),I-1)) &
       E0((IA-1)*NCFB+IB)=E0((IA-1)*NCFB+IB)+EPSILON(I,0,0,0)
     ENDDO
    ENDDO
   ENDDO
   !
   ! POPULATE TARGET STATES
   ! 
   ALLOCATE(ITARGETS(IOPTN(59)))
   PREVIOUSMIN=-1.0D99
   DO I=1,IOPTN(59)
    CURRENTMIN=1.0D99
    DO J=1,NALL
     IF ((E0(J)-PREVIOUSMIN > 1.0D-9).AND.(E0(J) < CURRENTMIN)) THEN
      CURRENTMIN=E0(J)
      ITARGETS(I)=J
     ENDIF
    ENDDO
    PREVIOUSMIN=CURRENTMIN
   ENDDO
   DO I=1,IOPTN(59)
    K=0
    DO J=1,NALL
     IF (DABS(E0(J)-E0(ITARGETS(I))) < 1.0D-9) K=K+1
    ENDDO
    WRITE(6,'(A,I3,A,I10,A,F15.10,I5,A)') 'TARGET',I,' = STATE',ITARGETS(I),' ZEROTH-ORDER ENERGY = ',E0(ITARGETS(I)),&
                                           K,'-FOLD DEGENERATE'
   ENDDO
   ! **********************
   ! LOOP OVER TARGET STATE
   ! **********************
   DO ISTATE=1,IOPTN(59)

    WRITE(6,'(/,A )') '/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\'
    WRITE(6,'(A,I3,A,I6)') '***** ',ISTATE,'-TH TARGET STATE:',ITARGETS(ISTATE)
    WRITE(6,'(A   )') '\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/'

    ! Here, we enumerate all degenerate determinants with the target state
    ! by using a string-based logic rather than looping over excitations.
    ! This is because we do not know the highest excitation level of degenerate
    ! references a priori. The parity is computed using the string-based logic 
    ! relative to the ground-state determinant with subroutine PRTY.
    NDEGEN=0
    DO I=1,NALL
     IF (DABS(E0(I)-E0(ITARGETS(ISTATE))) < 1.0D-9) NDEGEN=NDEGEN+1
    ENDDO
    WRITE(6,'(A,I3)') 'DEG.DEGENERACY = ',NDEGEN
    ALLOCATE(MAP(NA,NDEGEN),SPINMAP(NA,NDEGEN),BACKMAP(NA,NDEGEN))
    ALLOCATE(INTERNAL1(NDEGEN,NDEGEN),INTERNAL2(NDEGEN,NDEGEN),INTERNAL3(NDEGEN),INTERNAL4(NDEGEN,NDEGEN))
    JSTATE=0
    DO IA=1,NCFA
     DO IB=1,NCFB
      IF (DABS(E0((IA-1)*NCFB+IB)-E0(ITARGETS(ISTATE))) < 1.0D-9) THEN
       JSTATE=JSTATE+1

       CBITSA=""
       CBITSB=""
       DO K=1,IALL(0,0,0)-IVIRTCORE
        IF (BTEST(CFHALFA(IA),K-1)) THEN
         CBITSA="1"//CBITSA
        ELSE
         CBITSA="0"//CBITSA
        ENDIF
        IF (BTEST(CFHALFB(IB),K-1)) THEN
         CBITSB="1"//CBITSB
        ELSE
         CBITSB="0"//CBITSB
        ENDIF
       ENDDO
       WRITE(6,'(/,I3,A,2(A,A,A))') JSTATE,'-TH DEGENERATE STATE:',' (',TRIM(CBITSA),')',' (',TRIM(CBITSB),')'

       K=0
       DO I=ICORE+1,IALL(0,0,0)-IVIRTCORE
        IF (BTEST(CFHALFA(IA),I-1)) THEN
         K=K+1
         MAP(K,JSTATE)=MAPA(I)
         SPINMAP(K,JSTATE)=1
         BACKMAP(K,JSTATE)=I
        ENDIF
       ENDDO
       DO I=ICORE+1,IALL(0,0,0)-IVIRTCORE
        IF (BTEST(CFHALFB(IB),I-1)) THEN
         K=K+1
         MAP(K,JSTATE)=MAPB(I)
         SPINMAP(K,JSTATE)=2
         BACKMAP(K,JSTATE)=I
        ENDIF
       ENDDO
       IF (K /= NO) CALL PABORT('ERROR IN DEGENERATE CC 1')
       DO I=ICORE+1,IALL(0,0,0)-IVIRTCORE
        IF (.NOT.BTEST(CFHALFA(IA),I-1)) THEN
         K=K+1
         MAP(K,JSTATE)=MAPA(I)
         SPINMAP(K,JSTATE)=1
         BACKMAP(K,JSTATE)=I
        ENDIF
       ENDDO
       DO I=ICORE+1,IALL(0,0,0)-IVIRTCORE
        IF (.NOT.BTEST(CFHALFB(IB),I-1)) THEN
         K=K+1
         MAP(K,JSTATE)=MAPB(I)
         SPINMAP(K,JSTATE)=2
         BACKMAP(K,JSTATE)=I
        ENDIF
       ENDDO
       IF (K /= 2*(IALL(0,0,0)-IVIRTCORE-ICORE)) CALL PABORT('ERROR IN DEGENERATE CC 2')

       WRITE(*,'(A,100I3:)') 'OCCUPIED ORBITALS:',(MAP(K,JSTATE),K=1,NO)
       WRITE(*,'(A,100A3:)') '                   ',(SPINNAME(SPINMAP(K,JSTATE)),K=1,NO)
       WRITE(*,'(A,100I3:)') '         ORIGINAL:',(BACKMAP(K,JSTATE),K=1,NO)
       WRITE(*,'(A,100I3:)') 'VIRTUAL  ORBITALS:',(MAP(K,JSTATE),K=NO+1,NA)
       WRITE(*,'(A,100A3:)') '                   ',(SPINNAME(SPINMAP(K,JSTATE)),K=NO+1,NA)
       WRITE(*,'(A,100I3:)') '         ORIGINAL:',(BACKMAP(K,JSTATE),K=NO+1,NA)
!      SORT
       CALL ISORT3(MAP(:,JSTATE),SPINMAP(:,JSTATE),BACKMAP(:,JSTATE),NO)
       CALL ISORT3(MAP(NO+1,JSTATE),SPINMAP(NO+1,JSTATE),BACKMAP(NO+1,JSTATE),NA-NO)
       CALL PRTY(NO,NA,MAP(:,JSTATE),INTERNAL3(JSTATE))
       WRITE(*,'(A,100I3:)') 'SORTED OCCUPIED ORBITALS:',(MAP(K,JSTATE),K=1,NO)
       WRITE(*,'(A,100A3:)') '                          ',(SPINNAME(SPINMAP(K,JSTATE)),K=1,NO)
       WRITE(*,'(A,100I3:)') '                ORIGINAL:',(BACKMAP(K,JSTATE),K=1,NO)
       WRITE(*,'(A,100I3:)') 'SORTED VIRTUAL  ORBITALS:',(MAP(K,JSTATE),K=NO+1,NA)
       WRITE(*,'(A,100A3:)') '                          ',(SPINNAME(SPINMAP(K,JSTATE)),K=NO+1,NA)
       WRITE(*,'(A,100I3:)') '                ORIGINAL:',(BACKMAP(K,JSTATE),K=NO+1,NA)
       WRITE(*,'(A,I3)')     'PARITY RE GROUND STATE  :',INTERNAL3(JSTATE)

      ENDIF
     ENDDO
    ENDDO

    ! IDENTIFY INTERNAL T AMPLITUDES
    DO N=1,NDEGEN
     WRITE(6,'(/,A,I3,A)') ' === ',N,'-TH DEGENERATE STATE'

     INTERNAL1(:,N)=-999 ! ERROR TRAP
     INTERNAL2(:,N)=-999 ! ERROR TRAP
     INTERNAL4(:,N)=1    ! DIAGONAL CASES

     ! ZERO-ELECTRON DIFFERENCE
     INTERNAL1(N,N)=0
     INTERNAL2(N,N)=1  ! FIRST DEGENERATE STATE IS ALWAYS ITSELF
     WRITE(6,'(/,A,I3,A,I6,A,I3,A)') '  --- STATE FOUND:',N,'-TH DEGEN.STATE =  E(',INTERNAL2(N,N),',',N,')'

     ALLOCATE(MAP2(NA))

     ! ONE-ELECTRON DIFFERENCE
     IF (LSINGLES) THEN
      DO P1=NO+1,NA
       DO H1=1,NO
        IF (DABS(EPS1(MAP(P1,N))-EPS1(MAP(H1,N))) < 1.0D-9) THEN
!        WRITE(6,'(/,A,I3,A,I3)') 'INTERNAL T1:',MAP(H1,N),' ->',MAP(P1,N)
         MAP2=MAP(:,N)
         CALL PRTY1(MAP(P1,N),MAP(H1,N),NO,NA,MAP2,KPRTY)
         IF (KPRTY==0) CYCLE
         DO JSTATE=1,NDEGEN
          LEXIST=.TRUE.
          DO I=1,NA
           IF (MAP(I,JSTATE) /= MAP2(I)) LEXIST=.FALSE.
          ENDDO
          IF (LEXIST) THEN
           INTERNAL1(JSTATE,N)=1
           INTERNAL2(JSTATE,N)=(MAP(P1,N)-1)*NA+MAP(H1,N)
           INTERNAL4(JSTATE,N)=KPRTY
           WRITE(6,'(/,A,I3,A,I6,A,I3,A,A,I3,A,I3,A,I3)') '  --- STATE FOUND:',JSTATE,&
            '-TH DEGEN.STATE = T1(',INTERNAL2(JSTATE,N),',',N,')',&
            ' INTERNAL T1:',MAP(H1,N),' ->',MAP(P1,N),&
            ' WITH EXCITATION PARITY ',KPRTY
          ENDIF
         ENDDO
        ENDIF
       ENDDO
      ENDDO
     ENDIF

     ! TWO-ELECTRON DIFFERENCE
     IF (LDOUBLES) THEN
      DO P1=NO+1,NA
       DO P2=NO+1,NA
        IF (LPERMUTATION.AND.(P1>=P2)) CYCLE
        DO H1=1,NO
         DO H2=1,NO
          IF (LPERMUTATION.AND.(H1>=H2)) CYCLE
          IF (DABS(EPS1(MAP(P1,N))+EPS1(MAP(P2,N))-EPS1(MAP(H1,N))-EPS1(MAP(H2,N))) < 1.0D-9) THEN
!          WRITE(6,'(/,A,2I3,A,2I3)') 'INTERNAL T2:',MAP(H1,N),MAP(H2,N),' ->',MAP(P1,N),MAP(P2,N)
           MAP2=MAP(:,N)
           CALL PRTY2(MAP(P1,N),MAP(P2,N),MAP(H1,N),MAP(H2,N),NO,NA,MAP2,KPRTY)
           IF (KPRTY==0) CYCLE
           DO JSTATE=1,NDEGEN
            LEXIST=.TRUE.
            DO I=1,NA
             IF (MAP(I,JSTATE) /= MAP2(I)) LEXIST=.FALSE.
            ENDDO
            IF (LEXIST) THEN
             INTERNAL1(JSTATE,N)=2
             INTERNAL2(JSTATE,N)=(((MAP(P1,N)-1)*NA+MAP(P2,N)-1)*NA+MAP(H1,N)-1)*NA+MAP(H2,N)
             INTERNAL4(JSTATE,N)=KPRTY
             WRITE(6,'(/,A,I3,A,I6,A,I3,A,A,2I3,A,2I3,A,I3)') '  --- STATE FOUND:',JSTATE,&
              '-TH DEGEN.STATE = T2(',INTERNAL2(JSTATE,N),',',N,')',&
              ' INTERNAL T2:',MAP(H1,N),MAP(H2,N),' ->',MAP(P1,N),MAP(P2,N),&
              ' WITH EXCITATION PARITY ',KPRTY
            ENDIF
           ENDDO
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDIF

     DEALLOCATE(MAP2)
    ENDDO
    WRITE(6,'(/,A)') 'EXCITATION MATRIX'
    CALL DUMP17(INTERNAL1,NDEGEN)
    WRITE(6,'(/,A)') 'INTERNAL T AMPLITUDES'
    CALL DUMP17(INTERNAL2,NDEGEN)
    WRITE(6,'(/,A)') 'INTERNAL T AMPLITUDES PARITY'
    CALL DUMP17(INTERNAL4,NDEGEN)
    WRITE(6,'(/,A,I5,100I6:)') 'STATE PARITY',(INTERNAL3(N),N=1,NDEGEN)

    ! CONSTRUCT NORMAL-ORDERED HAMILTONIAN
    ALLOCATE(D_F1(NA,NA,NDEGEN))
    DO N=1,NDEGEN
     D_F1(:,:,N)=0.0D0
     DO P1=1,NA
      DO P2=1,NA
       D_F1(MAP(P1,N),MAP(P2,N),N)=HCORE1(MAP(P1,N),MAP(P2,N))
       DO H1=1,NO
        D_F1(MAP(P1,N),MAP(P2,N),N)=D_F1(MAP(P1,N),MAP(P2,N),N)+V2(MAP(P1,N),MAP(H1,N),MAP(P2,N),MAP(H1,N))
       ENDDO
      ENDDO
     ENDDO
     IF (ICORE > 0) THEN
      DO P1=1,NA
!write(*,*) P1,' map,spinmap,backmap=',map(P1,N),spinmap(P1,N),backmap(P1,N)
       DO P2=1,NA
        DO I=1,ICORE
         IF (SPINMAP(P1,N)==SPINMAP(P2,N)) D_F1(MAP(P1,N),MAP(P2,N),N)=D_F1(MAP(P1,N),MAP(P2,N),N) &
          +2.0D0*G(BACKMAP(P1,N),BACKMAP(P2,N),I,I)-G(BACKMAP(P1,N),I,I,BACKMAP(P2,N))
        ENDDO
       ENDDO
      ENDDO
     ENDIF
!    WRITE(6,'(/,A,I3)') 'FOCK MATRIX FOR DEGENERATE REF ',N
!    CALL DUMP5(D_F1(:,:,N),NA)
     EHF=NUCLEAR_REPULSION
     DO H1=1,NO
!write(*,'(2i3,2f20.10)') h1,map(h1,n),F1(MAP(H1,N),MAP(H1,N)),F1(h1,h1)
      EHF=EHF+HCORE1(MAP(H1,N),MAP(H1,N))
      DO H2=1,NO
!write(*,'(4i3,2f20.10)') h1,h2,map(h1,n),map(h2,n),V2(MAP(H1,N),MAP(H2,N),MAP(H1,N),MAP(H2,N)),v2(h1,h2,h1,h2)
       EHF=EHF+0.5D0*V2(MAP(H1,N),MAP(H2,N),MAP(H1,N),MAP(H2,N))
      ENDDO
     ENDDO
     IF (ICORE > 0) THEN
      DO I=1,ICORE
       EHF=EHF+2.0D0*H(I,I)
       DO J=1,ICORE
        EHF=EHF+2.0D0*G(I,I,J,J)-G(I,J,J,I)
       ENDDO
       DO P1=1,NO
        EHF=EHF+2.0D0*G(I,I,BACKMAP(P1,N),BACKMAP(P1,N))-G(I,BACKMAP(P1,N),BACKMAP(P1,N),I)
       ENDDO
      ENDDO
     ENDIF
     WRITE(6,'(/,A,I3,A,F20.10,A)') 'HF ENERGY FOR DEGENERATE REF ',N,' = ',EHF,' HARTREE'
    ENDDO

!   DO N=1,NDEGEN
!    IF (INTERNAL1(N)==0) THEN
!     WRITE(6,'(I3,A,I9,A)') N,'-TH DEG.REF  E(',INTERNAL2(N),')'
!    ELSE IF (INTERNAL1(N)==1) THEN
!     WRITE(6,'(I3,A,I9,A)') N,'-TH DEG.REF T1(',INTERNAL2(N),')'
!    ELSE IF (INTERNAL1(N)==2) THEN
!     WRITE(6,'(I3,A,I9,A)') N,'-TH DEG.REF T2(',INTERNAL2(N),')'
!    ELSE IF (INTERNAL1(N)==3) THEN
!     WRITE(6,'(I3,A,I9,A)') N,'-TH DEG.REF T3(',INTERNAL2(N),')'
!    ELSE
!     CALL PABORT('INTERNAL T-AMPLITUDES NOT IDENTIFIED')
!    ENDIF
!   ENDDO

    DO N=1,NDEGEN
     WRITE(APPEND,'(I3)') N
     OPEN(100+N,FILE=TRIM(COPTN(1))//'.T.'//APPEND,FORM='UNFORMATTED')
     OPEN(400+N,FILE=TRIM(COPTN(1))//'.delta.'//APPEND,FORM='UNFORMATTED')
    ENDDO

    ALLOCATE(D_I0_0(1,NDEGEN),D_I0E_0(1,NDEGEN))
    IF (LSINGLES) ALLOCATE(D_I0_1(NA**2,NDEGEN),D_T1(NA**2,NDEGEN))
    IF (LDOUBLES) ALLOCATE(D_I0_2(NA**4,NDEGEN),D_T2(NA**4,NDEGEN))
    IF (LTRIPLES) ALLOCATE(D_I0_3(NA**6,NDEGEN),D_T3(NA**6,NDEGEN))
    IF (LSINGLES) ALLOCATE(D_I0E_1(NA**2,NDEGEN),D_T1E(NA**2,NDEGEN))
    IF (LDOUBLES) ALLOCATE(D_I0E_2(NA**4,NDEGEN),D_T2E(NA**4,NDEGEN))
    IF (LTRIPLES) ALLOCATE(D_I0E_3(NA**6,NDEGEN),D_T3E(NA**6,NDEGEN))
    IF (LDISCONNECTED.AND.LSINGLES) ALLOCATE(D_S1(NA**2,NDEGEN),D_S1E(NA**2,NDEGEN))
    IF (LDISCONNECTED.AND.LSINGLES.AND.LDOUBLES) ALLOCATE(D_S2(NA**4,NDEGEN),D_S2E(NA**4,NDEGEN))
    IF (LDISCONNECTED.AND.LSINGLES.AND.LDOUBLES.AND.LTRIPLES) ALLOCATE(D_S3(NA**6,NDEGEN),D_S3E(NA**6,NDEGEN))
 
    IF (LSINGLES) D_T1=0.0D0
    IF (LDOUBLES) D_T2=0.0D0
    IF (LTRIPLES) D_T3=0.0D0
    DO N=1,NDEGEN
     REWIND(100+N)
     IF (LSINGLES) WRITE(100+N) d_t1(:,N)
     IF (LDOUBLES) WRITE(100+N) d_t2(:,N)
     IF (LTRIPLES) WRITE(100+N) d_t3(:,N)
    ENDDO
  
    ALLOCATE(HMAT(NDEGEN,NDEGEN),EMAT(NDEGEN,NDEGEN),SMAT(NDEGEN,NDEGEN),TMP(NDEGEN,NDEGEN))
    ALLOCATE(ER(NDEGEN),EI(NDEGEN),VR(NDEGEN,NDEGEN),VL(1,NDEGEN),WK(4*NDEGEN))
    ALLOCATE(DEVSQ1(NDEGEN),DEVSQ2(NDEGEN),DEVSQ3(NDEGEN))

    WRITE(6,'(/,A)') '  ITER       RESIDUAL    CORRELATION        TOTAL ENERGY'
    WRITE(6,'(A)')   '--------------------------------------------------------'
 
    DO ITER=1,MAXITER
!write(*,*) 'iter=',iter

     SMAT=0.0D0
     HMAT=0.0D0
 
     DO N=1,NDEGEN
!write(*,*) 'n=',n,'/',ndegen

      REWIND(100+N)
      DO I=1,ITER
       IF (LSINGLES) READ(100+N) D_T1(:,N)
       IF (LDOUBLES) READ(100+N) D_T2(:,N)
       IF (LTRIPLES) READ(100+N) D_T3(:,N)
      ENDDO

      IF (COPTN(110)=='CCS') THEN

       ALLOCATE(I1(NA**2),I1E(NA**2),I2(NA**2),I2E(NA**2))
       IF (LDISCONNECTED) THEN
        CALL p_m_ccs_e(NA,d_f1(:,:,N),f1e,d_i0_0(:,N),d_i0e_0,i1,i1e,MAP(:,N),NO,d_t1(:,N),d_t1e(:,N),v2,v2e)
        CALL p_m_ccs_t1_disconnected(NA,d_f1(:,:,N),f1e,d_i0_1(:,N),d_i0e_1(:,N),i1,i1e,i2,i2e,MAP(:,N),NO,&
                                     d_t1(:,N),d_t1e(:,N),v2,v2e)
        D_S1(:,N)=D_T1(:,N)
        DO M=1,NDEGEN
         IF (INTERNAL1(M,N)==0) THEN
          HMAT(M,N)=D_I0_0(1,N)
          SMAT(M,N)=1.0D0
         ELSE IF (IABS(INTERNAL1(M,N))==1) THEN
          HMAT(M,N)=D_I0_1(INTERNAL2(M,N),N)*DFLOAT(INTERNAL3(M)*INTERNAL3(N)*INTERNAL4(M,N))
          SMAT(M,N)=D_T1(INTERNAL2(M,N),N)*DFLOAT(INTERNAL3(M)*INTERNAL3(N)*INTERNAL4(M,N))
         ENDIF
        ENDDO
       ELSE
        CALL PABORT('DEGENERATE CC EQUATIONS ARE DISCONNECTED')
!       IF (LPERMUTATION) THEN
!       ELSE
!       ENDIF
       ENDIF
       DEALLOCATE(I1,I1E,I2,I2E)

      ELSE IF (COPTN(110)=='CCSD') THEN

       ALLOCATE(I1(NA**4),I1E(NA**4),I2(NA**4),I2E(NA**4),I3(NA**4),I3E(NA**4))
       IF (LDISCONNECTED) THEN
        CALL p_m_ccsd_e(NA,d_f1(:,:,N),f1e,d_i0_0(:,N),d_i0e_0(:,N),i1,i1e,MAP(:,N),NO,&
                        d_t1(:,N),d_t1e(:,N),d_t2(:,N),d_t2e(:,N),v2,v2e)
        CALL p_m_ccsd_t1_disconnected(NA,d_f1(:,:,N),f1e,d_i0_1(:,N),d_i0e_1(:,N),i1,i1e,i2,i2e,MAP(:,N),NO,&
                                      d_t1(:,N),d_t1e(:,N),d_t2(:,N),d_t2e(:,N),v2,v2e)
        CALL p_m_ccsd_t2_disconnected(NA,d_f1(:,:,N),f1e,d_i0_2(:,N),d_i0e_2(:,N),i1,i1e,i2,i2e,i3,i3e,&
                                      MAP(:,N),NO,d_t1(:,N),d_t1e(:,N),d_t2(:,N),d_t2e(:,N),v2,v2e)
        CALL p_m_ccsd_t2_overlap(NA,d_s2(:,N),d_s2e(:,N),MAP(:,N),NO,d_t1(:,N),d_t1e(:,N),d_t2(:,N),d_t2e(:,N))
        D_S1(:,N)=D_T1(:,N)
        DO M=1,NDEGEN
         IF (INTERNAL1(M,N)==0) THEN
          HMAT(M,N)=D_I0_0(1,N)
          SMAT(M,N)=1.0D0
         ELSE IF (IABS(INTERNAL1(M,N))==1) THEN
          HMAT(M,N)=D_I0_1(INTERNAL2(M,N),N)*DFLOAT(INTERNAL3(M)*INTERNAL3(N)*INTERNAL4(M,N))
          SMAT(M,N)=D_T1(INTERNAL2(M,N),N)*DFLOAT(INTERNAL3(M)*INTERNAL3(N)*INTERNAL4(M,N))
         ELSE IF (IABS(INTERNAL1(M,N))==2) THEN
          HMAT(M,N)=D_I0_2(INTERNAL2(M,N),N)*DFLOAT(INTERNAL3(M)*INTERNAL3(N)*INTERNAL4(M,N))
          SMAT(M,N)=D_S2(INTERNAL2(M,N),N)*DFLOAT(INTERNAL3(M)*INTERNAL3(N)*INTERNAL4(M,N))
         ENDIF
        ENDDO
!do i=no+1,na
!do k=1,no
! if (dabs(d_t1((map(i,n)-1)*na+map(k,n),n))>1.0d-9) then
!  write(*,'(i5,2i3,f20.10)') n,map(i,n),map(k,n),d_t1((map(i,n)-1)*na+map(k,n),n)
! endif
!enddo
!enddo
 do i=no+1,na
 do j=no+1,na
 if (i>j) cycle
 do k=1,no
 do l=1,no
!if (k>l) cycle
! if (dabs(d_t2((((map(i,n)-1)*na+map(j,n)-1)*na+map(k,n)-1)*na+map(l,n),n))>1.0d-9) then
 if ((n==1).and.(map(i,n)==11).and.(map(j,n)==14).and.(map(k,n)==4).and.(map(l,n)==8)) &
  write(*,'(i5,4i3,2f20.10)') n,map(i,n),map(j,n),map(k,n),map(l,n),&
  d_t2((((map(i,n)-1)*na+map(j,n)-1)*na+map(k,n)-1)*na+map(l,n),n),&
  d_i0_2((((map(i,n)-1)*na+map(j,n)-1)*na+map(k,n)-1)*na+map(l,n),n)
! endif
 enddo
 enddo
 enddo
 enddo
       ELSE
        CALL PABORT('DEGENERATE CC EQUATIONS ARE DISCONNECTED')
!       IF (LPERMUTATION) THEN
!        CALL p_m_ccsd_e(NA,f1,f1e,i0_0,i0e_0,i1,i1e,NO,t1,t1e,t2,t2e,v2,v2e)
!        CALL p_m_ccsd_t1(NA,f1,f1e,i0_1,i0e_1,i1,i1e,i2,i2e,NO,t1,t1e,t2,t2e,v2,v2e)
!        CALL p_m_ccsd_t2(NA,f1,f1e,i0_2,i0e_2,i1,i1e,i2,i2e,i3,i3e,NO,t1,t1e,t2,t2e,v2,v2e)
!       ELSE
!        CALL np_m_ccsd_e(NA,f1,i0_0,i1,NO,t1,t2,v2)
!        CALL np_m_ccsd_t1(NA,f1,i0_1,i1,i2,NO,t1,t2,v2)
!        CALL np_m_ccsd_t2(NA,f1,i0_2,i1,i2,i3,NO,t1,t2,v2)
!       ENDIF
       ENDIF
       DEALLOCATE(I1,I1E,I2,I2E,I3,I3E)
      ENDIF

     ENDDO

     IF (IOPTN(9) > 1) THEN
      WRITE(6,'(A)') 'OVERLAP MATRIX'
      CALL DUMP5LONG(SMAT,NDEGEN)
      TMP=HMAT
      DO M=1,NDEGEN
       TMP(M,M)=TMP(M,M)+EHF
      ENDDO
      WRITE(6,'(A)') 'HAMILTONIAN MATRIX'
      CALL DUMP5LONG(TMP,NDEGEN)
!     WRITE(6,'(A)') 'HAMILTONIAN MATRIX (RAW)'
!     CALL DUMP5(HMAT,NDEGEN)
     ENDIF

     ! MULTIPLY SMAT-INVERSE WITH HMAT
     ALLOCATE(INDX(NDEGEN),B(NDEGEN,NDEGEN),BS(NDEGEN,NDEGEN),C(NDEGEN),CS(NDEGEN))
     B=SMAT
     BS=B
     DO N=1,NDEGEN
      C=HMAT(:,N)
      CS=C
      CALL LUDCMP(B,NDEGEN,NDEGEN,INDX,D)
      CALL LUBKSB(B,NDEGEN,NDEGEN,INDX,C)
      CALL MPROVE(BS,B,NDEGEN,NDEGEN,INDX,CS,C)
      EMAT(:,N)=C
     ENDDO
     DEALLOCATE(INDX,B,BS,C,CS)
     TMP=EMAT
     DO M=1,NDEGEN
      TMP(M,M)=TMP(M,M)+EHF
     ENDDO
     IF (IOPTN(9) > 1) THEN
      WRITE(6,'(A)') 'ENERGY MATRIX'
      CALL DUMP5(TMP,NDEGEN)
      WRITE(6,'(A)') 'ENERGY MATRIX (RAW)'
      CALL DUMP5(EMAT,NDEGEN)
     ENDIF

     ! UPDATE T-AMPLITUDES
     DEVSQ1=0.0D0
     DEVSQ2=0.0D0
     DEVSQ3=0.0D0
     DO N=1,NDEGEN 
       
      ! SIMPLE RELAXATION
      IF (COPTN(110)=='CCSD') THEN
       CALL T_M_RELAX2(1,LSINGLES,LDOUBLES,N,NA,NO,MAP,NDEGEN,SMAT,EMAT,D_I0_1(:,N),NA**2,D_T1,D_S1,D_S2,EPS1,&
                       INTERNAL1,INTERNAL2,INTERNAL3,DEVSQ1(N),LPERMUTATION)
       CALL T_M_RELAX2(2,LSINGLES,LDOUBLES,N,NA,NO,MAP,NDEGEN,SMAT,EMAT,D_I0_2(:,N),NA**4,D_T2,D_S1,D_S2,EPS1,&
                       INTERNAL1,INTERNAL2,INTERNAL3,DEVSQ2(N),LPERMUTATION)
      ELSE IF (COPTN(110)=='CCS') THEN
       CALL T_M_RELAX1(1,LSINGLES,LDOUBLES,N,NA,NO,MAP,NDEGEN,SMAT,EMAT,D_I0_1(:,N),NA**2,D_T1,D_S1,EPS1,&
                       INTERNAL1,INTERNAL2,INTERNAL3,DEVSQ1(N),LPERMUTATION)
      ENDIF

      REWIND(400+N)
      IF (ITER > 1) THEN
       DO I=1,ITER-1
        IF (LSINGLES) READ(400+N) D_I0E_1(:,N)
        IF (LDOUBLES) READ(400+N) D_I0E_2(:,N)
        IF (LTRIPLES) READ(400+N) D_I0E_3(:,N)
       ENDDO
      ENDIF
      IF (LSINGLES) WRITE(400+N) D_I0_1(:,N)
      IF (LDOUBLES) WRITE(400+N) D_I0_2(:,N)
      IF (LTRIPLES) WRITE(400+N) D_I0_3(:,N)

      REWIND(100+N)
      DO I=1,ITER
       IF (LSINGLES) READ(100+N) D_T1E(:,N)
       IF (LDOUBLES) READ(100+N) D_T2E(:,N)
       IF (LTRIPLES) READ(100+N) D_T3E(:,N)
      ENDDO
      IF (LSINGLES) WRITE(100+N) D_T1(:,N)
      IF (LDOUBLES) WRITE(100+N) D_T2(:,N)
      IF (LTRIPLES) WRITE(100+N) D_T3(:,N)

     ENDDO

     ! DIAGONALIZE EFFECTIVE ENERGY MATRIX
     CALL DGEEV('N','V',NDEGEN,EMAT,NDEGEN,ER,EI,VL,1,VR,NDEGEN,WK,4*NDEGEN,INFO)
     IF (INFO /= 0) CALL PABORT('DGEEV FAILED TO DIAGONALIZE A MATRIX')
     CALL PIKSRT2(NDEGEN,NDEGEN,ER,EI,VR,WK)

     LEXIST=.FALSE.
     DO N=1,NDEGEN 
      WRITE(6,'(I6,2F15.10,F20.10)') ITER,DSQRT(DEVSQ1(N)+DEVSQ2(N)+DEVSQ3(N)),ER(N),EHF+ER(N)
      IF (DSQRT(DEVSQ1(N)+DEVSQ2(N)+DEVSQ3(N)) > DOPTN(62)) LEXIST=.TRUE.
     ENDDO
     WRITE(6,'(A)')   '........................................................'
     IF (.NOT.LEXIST) THEN
      EXIT
     ENDIF
       
     ! DIIS
     IF ((LDIIS).AND.(MOD(ITER,IDIIS) == 0)) THEN

      ALLOCATE(INDX(21),B(21,21),BS(21,21),C(21),CS(21))
      IORDER=MIN(ITER,20)

      DO N=1,NDEGEN

      B=0.0D0
      DO I=1,IORDER
       REWIND(400+N)
       IF (ITER > I) THEN
        DO K=1,ITER-I
         IF (LSINGLES) READ(400+N) d_i0e_1(:,N)
         IF (LDOUBLES) READ(400+N) d_i0e_2(:,N)
         IF (LTRIPLES) READ(400+N) d_i0e_3(:,N)
        ENDDO
       ENDIF
       IF (LSINGLES) READ(400+N) d_i0_1(:,N)
       IF (LDOUBLES) READ(400+N) d_i0_2(:,N)
       IF (LTRIPLES) READ(400+N) d_i0_3(:,N)
       DO J=1,IORDER
        REWIND(400+N)
        IF (ITER > J) THEN
         DO K=1,ITER-J
          IF (LSINGLES) READ(400+N) d_i0e_1(:,N)
          IF (LDOUBLES) READ(400+N) d_i0e_2(:,N)
          IF (LTRIPLES) READ(400+N) d_i0e_3(:,N)
         ENDDO
        ENDIF
        IF (LSINGLES) READ(400+N) d_i0e_1(:,N)
        IF (LDOUBLES) READ(400+N) d_i0e_2(:,N)
        IF (LTRIPLES) READ(400+N) d_i0e_3(:,N)
        R1=0.0D0
        R2=0.0D0
        R3=0.0D0
        IF (LSINGLES) CALL T_M_OVERLAP(1,NA,NO,MAP(:,N),d_i0_1(:,N),d_i0e_1(:,N),R1,LPERMUTATION)
        IF (LDOUBLES) CALL T_M_OVERLAP(2,NA,NO,MAP(:,N),d_i0_2(:,N),d_i0e_2(:,N),R2,LPERMUTATION)
        IF (LTRIPLES) CALL T_M_OVERLAP(3,NA,NO,MAP(:,N),d_i0_3(:,N),d_i0e_3(:,N),R3,LPERMUTATION)
        B(J,I)=R1+R2+R3
       ENDDO
       B(I,IORDER+1)=-1.0D0
       B(IORDER+1,I)=-1.0D0
       C(I)=0.0D0
      ENDDO
      B(IORDER+1,IORDER+1)=0.0D0
      C(IORDER+1)=-1.0D0
      BS=B
      CS=C
!     WRITE(6,'(A)') 'PULAY MATRIX'
!     CALL DUMP5(B,21)
      CALL LUDCMP(B,IORDER+1,21,INDX,D)
      CALL LUBKSB(B,IORDER+1,21,INDX,C)
      CALL MPROVE(BS,B,IORDER+1,21,INDX,CS,C)
      WRITE(6,'(A,100F7.4:)') 'DIIS VECTOR',(C(I),I=1,IORDER+1)
      IF (LSINGLES) d_t1(:,N)=0.0D0
      IF (LDOUBLES) d_t2(:,N)=0.0D0
      IF (LTRIPLES) d_t3(:,N)=0.0D0
      REWIND(100+N)
      IF (ITER > IORDER) THEN
       DO I=1,ITER-IORDER
        IF (LSINGLES) READ(100+N) d_t1e(:,N)
        IF (LDOUBLES) READ(100+N) d_t2e(:,N)       
        IF (LTRIPLES) READ(100+N) d_t3e(:,N)
       ENDDO
      ENDIF
      DO I=1,IORDER 
       IF (LSINGLES) READ(100+N) d_t1e(:,N)
       IF (LDOUBLES) READ(100+N) d_t2e(:,N)
       IF (LTRIPLES) READ(100+N) d_t3e(:,N)
       IF (LSINGLES) d_t1(:,N)=d_t1(:,N)+d_t1e(:,N)*C(IORDER-I+1)
       IF (LDOUBLES) d_t2(:,N)=d_t2(:,N)+d_t2e(:,N)*C(IORDER-I+1)
       IF (LTRIPLES) d_t3(:,N)=d_t3(:,N)+d_t3e(:,N)*C(IORDER-I+1)
      ENDDO
      REWIND(100+N) 
      DO I=1,ITER
       IF (LSINGLES) READ(100+N) d_t1e(:,N)
       IF (LDOUBLES) READ(100+N) d_t2e(:,N)
       IF (LTRIPLES) READ(100+N) d_t3e(:,N)
      ENDDO
      IF (LSINGLES) WRITE(100+N) d_t1(:,N)
      IF (LDOUBLES) WRITE(100+N) d_t2(:,N)
      IF (LTRIPLES) WRITE(100+N) d_t3(:,N)

      ENDDO ! END LOOP OVER NDEGEN

      DEALLOCATE(INDX,B,BS,C,CS)
     ENDIF ! END OF DIIS

    ENDDO ! END LOOP OVER ITERATIONS
    IF (ITER>=MAXITER) CALL PABORT('MAXITER REACHED')
    WRITE(6,'(A)') '------------------------- TCE --------------------------'
    WRITE(6,'(A)') 'CC ITERATION CONVERGED'

    DEALLOCATE(ER,EI,VR,VL,WK)
    DEALLOCATE(DEVSQ1,DEVSQ2,DEVSQ3)
    DEALLOCATE(HMAT,EMAT,SMAT,TMP)
    DEALLOCATE(D_I0_0,D_I0E_0)
    IF (LSINGLES) DEALLOCATE(D_I0_1,D_I0E_1,D_T1,D_T1E)
    IF (LDOUBLES) DEALLOCATE(D_I0_2,D_I0E_2,D_T2,D_T2E)
    IF (LTRIPLES) DEALLOCATE(D_I0_3,D_I0E_3,D_T3,D_T3E)
    IF (LDISCONNECTED.AND.LSINGLES) DEALLOCATE(D_S1,D_S1E)
    IF (LDISCONNECTED.AND.LSINGLES.AND.LDOUBLES) DEALLOCATE(D_S2,D_S2E)
    IF (LDISCONNECTED.AND.LSINGLES.AND.LDOUBLES.AND.LTRIPLES) DEALLOCATE(D_S3,D_S3E)
 
    DO N=1,NDEGEN
     CLOSE(100+N)
     CLOSE(400+N)
    ENDDO
 
    DEALLOCATE(D_F1)
    DEALLOCATE(INTERNAL1,INTERNAL2,INTERNAL3,INTERNAL4)
    DEALLOCATE(MAP,SPINMAP,BACKMAP)

   ENDDO ! END LOOP OVER TARGET STATE
 
   DEALLOCATE(CFHALFA,ADDRSSA,NORDERA)
   DEALLOCATE(CFHALFB,ADDRSSB,NORDERB)

   ENDIF

   DEALLOCATE(H,G,MAPA,MAPB)
   DEALLOCATE(EPS1,HCORE1,F1,V2)

   RETURN
END SUBROUTINE



SUBROUTINE T_RELAX(ORDER,N,NOCC,I0,T,F1,R,LPERMUTATION)

   IMPLICIT NONE
   INTEGER :: ORDER
   LOGICAL :: LPERMUTATION
   INTEGER :: N
   INTEGER :: nocc
   INTEGER :: p4
   INTEGER :: p7
   INTEGER :: h1
   INTEGER :: p8
   INTEGER :: h2 
   INTEGER :: h9
   INTEGER :: p1
   INTEGER :: p2
   INTEGER :: h3
   INTEGER :: h4
   REAL*8 :: i0(*)
   REAL*8 :: f1(*)
   REAL*8 :: t(*)
   REAL*8 :: R
   
   R=0.0D0
   IF (ORDER==3) THEN
   DO p4=nocc+1,N
   DO p7=nocc+1,N
   DO p8=nocc+1,N
   IF (LPERMUTATION.AND.(p4>p7)) CYCLE
   IF (LPERMUTATION.AND.(p7>p8)) CYCLE
   DO h1=1,nocc
   DO h2=1,nocc
   DO h9=1,nocc
   IF (LPERMUTATION.AND.(h1>h2)) CYCLE
   IF (LPERMUTATION.AND.(h2>h9)) CYCLE
   t(((((((p4-1)*N+p7-1)*N+p8-1)*N+h1-1)*N+h2-1)*N+h9))=&
   t(((((((p4-1)*N+p7-1)*N+p8-1)*N+h1-1)*N+h2-1)*N+h9))-&
   i0(((((((p4-1)*N+p7-1)*N+p8-1)*N+h1-1)*N+h2-1)*N+h9))/&
   (f1((p4-1)*N+p4)+f1((p7-1)*N+p7)+f1((p8-1)*N+p8)-f1((h1-1)*N+h1)-f1((h2-1)*N+h2)-f1((h9-1)*N+h9))
   R=R+i0(((((((p4-1)*N+p7-1)*N+p8-1)*N+h1-1)*N+h2-1)*N+h9))**2
   ENDDO
   ENDDO
   ENDDO
   ENDDO
   ENDDO
   ENDDO
   ELSE IF (ORDER==2) THEN
   DO p1=nocc+1,N
   DO p2=nocc+1,N
   IF (LPERMUTATION.AND.(p1>p2)) CYCLE
   DO h3=1,nocc
   DO h4=1,nocc
   IF (LPERMUTATION.AND.(h3>h4)) CYCLE
   t(((((p1-1)*N+p2-1)*N+h3-1)*N+h4))=t(((((p1-1)*N+p2-1)*N+h3-1)*N+h4))-i0(((((p1-1)*N+p2-1)*N+h3-1)*N+h4)) &
   /(F1((p1-1)*N+p1)+F1((p2-1)*N+p2)-F1((h3-1)*N+h3)-F1((h4-1)*N+h4))
   R=R+i0(((((p1-1)*N+p2-1)*N+h3-1)*N+h4))**2
!write(*,'(4i4,f20.10,f12.9)') p1,p2,h3,h4,i0(((((p1-1)*N+p2-1)*N+h3-1)*N+h4)),r
   ENDDO
   ENDDO
   ENDDO
   ENDDO
   ELSE IF (ORDER==1) THEN
   DO p1=nocc+1,N
   DO h2=1,nocc
   t(((p1-1)*N+h2))=t(((p1-1)*N+h2))-i0(((p1-1)*N+h2)) &
   /(F1((p1-1)*N+p1)-F1((h2-1)*N+h2))
   R=R+i0(((p1-1)*N+h2))**2
!write(*,'(2i4,f20.10,f12.9)') p1,h2,i0((((p1)-1)*N+(h2))),r
   ENDDO
   ENDDO
   ENDIF
   RETURN
END SUBROUTINE



SUBROUTINE T_OVERLAP(ORDER,N,NOCC,I0,T,R,LPERMUTATION)

   IMPLICIT NONE
   INTEGER :: ORDER
   LOGICAL :: LPERMUTATION
   INTEGER :: N
   INTEGER :: nocc
   INTEGER :: p4
   INTEGER :: p7
   INTEGER :: h1
   INTEGER :: p8
   INTEGER :: h2
   INTEGER :: h9
   INTEGER :: p1
   INTEGER :: p2
   INTEGER :: h3
   INTEGER :: h4
   REAL*8 :: i0(*)
   REAL*8 :: t(*)
   REAL*8 :: R
  
   R=0.0D0
   IF (ORDER==3) THEN
   DO p4=nocc+1,N
   DO p7=nocc+1,N
   DO p8=nocc+1,N
   IF (LPERMUTATION.AND.(p4>p7)) CYCLE
   IF (LPERMUTATION.AND.(p7>p8)) CYCLE
   DO h1=1,nocc
   DO h2=1,nocc
   DO h9=1,nocc
   IF (LPERMUTATION.AND.(h1>h2)) CYCLE
   IF (LPERMUTATION.AND.(h2>h9)) CYCLE
   R=R+i0(((((((p4-1)*N+p7-1)*N+p8-1)*N+h1-1)*N+h2-1)*N+h9))*t(((((((p4-1)*N+p7-1)*N+p8-1)*N+h1-1)*N+h2-1)*N+h9))
   ENDDO
   ENDDO
   ENDDO
   ENDDO
   ENDDO
   ENDDO
   ELSE IF (ORDER==2) THEN 
   DO p1=nocc+1,N
   DO p2=nocc+1,N
   IF (LPERMUTATION.AND.(p1>p2)) CYCLE
   DO h3=1,nocc
   DO h4=1,nocc
   IF (LPERMUTATION.AND.(h3>h4)) CYCLE
   R=R+i0(((((p1-1)*N+p2-1)*N+h3-1)*N+h4))*t(((((p1-1)*N+p2-1)*N+h3-1)*N+h4))
   ENDDO
   ENDDO
   ENDDO
   ENDDO
   ELSE IF (ORDER==1) THEN
   DO p1=nocc+1,N
   DO h2=1,nocc
   R=R+i0(((p1-1)*N+h2))*t(((p1-1)*N+h2))
   ENDDO
   ENDDO
   ENDIF
   RETURN
END SUBROUTINE



SUBROUTINE T_M_OVERLAP(ORDER,N,NOCC,M,I0,T,R,LPERMUTATION)

   IMPLICIT NONE
   INTEGER :: ORDER
   LOGICAL :: LPERMUTATION
   INTEGER :: N
   INTEGER :: nocc
   INTEGER :: p4
   INTEGER :: p7
   INTEGER :: h1
   INTEGER :: p8
   INTEGER :: h2
   INTEGER :: h9
   INTEGER :: p1
   INTEGER :: p2
   INTEGER :: h3
   INTEGER :: h4
   INTEGER :: m(*)
   REAL*8 :: i0(*)
   REAL*8 :: t(*)
   REAL*8 :: R

   R=0.0D0
   IF (ORDER==3) THEN
   DO p4=nocc+1,N
   DO p7=nocc+1,N
   DO p8=nocc+1,N
   IF (LPERMUTATION.AND.(p4>=p7)) CYCLE
   IF (LPERMUTATION.AND.(p7>=p8)) CYCLE
   DO h1=1,nocc
   DO h2=1,nocc
   DO h9=1,nocc
   IF (LPERMUTATION.AND.(h1>=h2)) CYCLE
   IF (LPERMUTATION.AND.(h2>=h9)) CYCLE
   R=R+i0(((((((m(p4)-1)*N+m(p7)-1)*N+m(p8)-1)*N+m(h1)-1)*N+m(h2)-1)*N+m(h9)))&
   *t(((((((m(p4)-1)*N+m(p7)-1)*N+m(p8)-1)*N+m(h1)-1)*N+m(h2)-1)*N+m(h9)))
   ENDDO
   ENDDO
   ENDDO
   ENDDO
   ENDDO
   ENDDO
   ELSE IF (ORDER==2) THEN
   DO p1=nocc+1,N
   DO p2=nocc+1,N
   IF (LPERMUTATION.AND.(p1>=p2)) CYCLE
   DO h3=1,nocc
   DO h4=1,nocc
   IF (LPERMUTATION.AND.(h3>=h4)) CYCLE
   R=R+i0(((((m(p1)-1)*N+m(p2)-1)*N+m(h3)-1)*N+m(h4)))*t(((((m(p1)-1)*N+m(p2)-1)*N+m(h3)-1)*N+m(h4)))
   ENDDO
   ENDDO
   ENDDO
   ENDDO
   ELSE IF (ORDER==1) THEN
   DO p1=nocc+1,N
   DO h2=1,nocc
   R=R+i0(((m(p1)-1)*N+m(h2)))*t(((m(p1)-1)*N+m(h2)))
   ENDDO
   ENDDO
   ENDIF
   RETURN
END SUBROUTINE



SUBROUTINE T_M_RELAX1(ORDER,LSINGLES,LDOUBLES,N,NA,NO,M,NDEGEN,S,E,I0,NT,T,S1,EPS,&
                      INTERNAL1,INTERNAL2,INTERNAL3,R,LPERMUTATION)
! FOR CCS

   IMPLICIT NONE
   INTEGER :: ORDER
   LOGICAL :: LSINGLES
   LOGICAL :: LDOUBLES
   LOGICAL :: LPERMUTATION
   INTEGER :: N
   INTEGER :: NDEGEN
   INTEGER :: NT
   INTEGER :: NA,NO
   INTEGER :: H1
   INTEGER :: H2 
   INTEGER :: P1
   INTEGER :: P2
   INTEGER :: H3
   INTEGER :: H4
   INTEGER :: M(NA,NDEGEN)          ! MAP
   INTEGER :: INTERNAL1(NDEGEN,NDEGEN)
   INTEGER :: INTERNAL2(NDEGEN,NDEGEN)
   INTEGER :: INTERNAL3(NDEGEN)
   INTEGER :: I,J,K,L
   REAL*8 :: I0(NT)
   REAL*8 :: EPS(NA)
   REAL*8 :: T(NT,NDEGEN)
   REAL*8 :: S1(NA**2,NDEGEN)
   REAL*8 :: R
   REAL*8 :: S(NDEGEN,NDEGEN)       ! SMAT
   REAL*8 :: E(NDEGEN,NDEGEN)       ! EMAT
   INTEGER :: MAP2(NA),MAP3(NA)
   INTEGER :: PA,PB,HA,HB
   LOGICAL :: LEXIST
   INTEGER :: NEXIST1,NEXIST2
   INTEGER :: IPRTY,JPRTY

   R=0.0D0

   IF (ORDER==1) THEN

   DO p1=NO+1,NA
   DO h2=1,NO
   IF (DABS(EPS(m(p1,N))-EPS(m(h2,N))) < 1.0D-9) THEN
    J=-1
    DO I=1,NDEGEN
     IF ((INTERNAL1(I,N)==1).AND.(INTERNAL2(I,N)==(m(P1,N)-1)*NA+m(h2,N))) J=I
    ENDDO
    IF (J /= -1) THEN
     i0(((m(p1,N)-1)*NA+m(h2,N)))=S(J,N)
     t(((m(p1,N)-1)*NA+m(h2,N)),N)=t(((m(p1,N)-1)*NA+m(h2,N)),N)-i0(((m(p1,N)-1)*NA+m(h2,N))) 
     R=R+i0(((m(p1,N)-1)*NA+m(h2,N)))**2
    ENDIF
   ELSE
    MAP2(:)=m(:,N)
    CALL PRTY1(m(P1,N),m(H2,N),NO,NA,MAP2,IPRTY)
    DO I=1,NDEGEN
     IF (LSINGLES) THEN
      NEXIST1=0
      DO pa=NO+1,NA
       DO ha=1,NO
        MAP3(:)=m(:,I)
        CALL PRTY1(m(PA,I),m(HA,I),NO,NA,MAP3,JPRTY)
        LEXIST=.TRUE.
        DO J=1,NA
         IF (MAP2(J)/=MAP3(J)) LEXIST=.FALSE.
        ENDDO
        IF (LEXIST) NEXIST1=NEXIST1+1
        IF (NEXIST1 > 1) EXIT
        IF (LEXIST) THEN
         i0(((m(p1,N)-1)*NA+m(h2,N)))=i0(((m(p1,N)-1)*NA+m(h2,N)))&
          -s1(((m(pa,I)-1)*NA+m(ha,I)),I)*E(I,N)&
          *DFLOAT(INTERNAL3(N)*IPRTY*INTERNAL3(I)*JPRTY)
        ENDIF
       ENDDO
      ENDDO
     ENDIF
     IF (LPERMUTATION.AND.(NEXIST1 > 1)) THEN
      WRITE(6,'(A)') 'ERROR: AN EXTERNAL DETERMINANT COULD NOT BE MAPPED ON A UNIQUE T-AMPLITUDE'
      STOP
     ENDIF
    ENDDO
    t(((m(p1,N)-1)*NA+m(h2,N)),N)=t(((m(p1,N)-1)*NA+m(h2,N)),N)-i0(((m(p1,N)-1)*NA+m(h2,N)))/(EPS(m(p1,N))-EPS(m(h2,N)))
    R=R+i0(((m(p1,N)-1)*NA+m(h2,N)))**2
   ENDIF
   ENDDO
   ENDDO

   ELSE
    WRITE(6,'(A)') 'ERROR IN T_M_RELAX1'
    STOP
   ENDIF

   RETURN
END SUBROUTINE



SUBROUTINE T_M_RELAX2(ORDER,LSINGLES,LDOUBLES,N,NA,NO,M,NDEGEN,S,E,I0,NT,T,S1,S2,EPS,&
                      INTERNAL1,INTERNAL2,INTERNAL3,R,LPERMUTATION)
! FOR CCSD

   IMPLICIT NONE
   INTEGER :: ORDER
   LOGICAL :: LSINGLES
   LOGICAL :: LDOUBLES
   LOGICAL :: LPERMUTATION
   INTEGER :: N
   INTEGER :: NDEGEN
   INTEGER :: NT
   INTEGER :: NA,NO
   INTEGER :: H1
   INTEGER :: H2 
   INTEGER :: P1
   INTEGER :: P2
   INTEGER :: H3
   INTEGER :: H4
   INTEGER :: M(NA,NDEGEN)          ! MAP
   INTEGER :: INTERNAL1(NDEGEN,NDEGEN)
   INTEGER :: INTERNAL2(NDEGEN,NDEGEN)
   INTEGER :: INTERNAL3(NDEGEN)
   INTEGER :: I,J,K,L
   REAL*8 :: I0(NT)
   REAL*8 :: EPS(NA)
   REAL*8 :: T(NT,NDEGEN)
   REAL*8 :: S1(NA**2,NDEGEN)
   REAL*8 :: S2(NA**4,NDEGEN)
   REAL*8 :: R
   REAL*8 :: S(NDEGEN,NDEGEN)       ! SMAT
   REAL*8 :: E(NDEGEN,NDEGEN)       ! EMAT
   INTEGER :: MAP2(NA),MAP3(NA)
   INTEGER :: PA,PB,HA,HB
   LOGICAL :: LEXIST
   INTEGER :: NEXIST1,NEXIST2
   INTEGER :: IPRTY,JPRTY

   R=0.0D0

   IF (ORDER==2) THEN

   DO p1=NO+1,NA
   DO p2=NO+1,NA
   IF (LPERMUTATION.AND.(p1>=p2)) CYCLE
   IF (p1==p2) CYCLE
   DO h3=1,NO
   DO h4=1,NO
   IF (LPERMUTATION.AND.(h3>=h4)) CYCLE
   IF (h3==h4) CYCLE
   IF (DABS(EPS(m(p1,N))+EPS(m(p2,N))-EPS(m(h3,N))-EPS(m(h4,N))) < 1.0D-9) THEN
    J=-1
    DO I=1,NDEGEN
     IF ((INTERNAL1(I,N)==2).AND.(INTERNAL2(I,N)==(((m(P1,N)-1)*NA+m(P2,N)-1)*NA+m(H3,N)-1)*NA+m(H4,N))) J=I
    ENDDO
    IF (J /= -1) THEN
     i0(((((m(p1,N)-1)*NA+m(p2,N)-1)*NA+m(h3,N)-1)*NA+m(h4,N)))=S(J,N)
     t(((((m(p1,N)-1)*NA+m(p2,N)-1)*NA+m(h3,N)-1)*NA+m(h4,N)),N)=t(((((m(p1,N)-1)*NA+m(p2,N)-1)*NA+m(h3,N)-1)*NA+m(h4,N)),N)&
     -i0(((((m(p1,N)-1)*NA+m(p2,N)-1)*NA+m(h3,N)-1)*NA+m(h4,N)))
     R=R+i0(((((m(p1,N)-1)*NA+m(p2,N)-1)*NA+m(h3,N)-1)*NA+m(h4,N)))**2
    ENDIF
   ELSE
    MAP2(:)=m(:,N)
    CALL PRTY2(m(P1,N),m(P2,N),m(H3,N),m(H4,N),NO,NA,MAP2,IPRTY)
    DO I=1,NDEGEN
     IF (LSINGLES) THEN
      NEXIST1=0
      DO pa=NO+1,NA
       DO ha=1,NO
        MAP3(:)=m(:,I)
        CALL PRTY1(m(PA,I),m(HA,I),NO,NA,MAP3,JPRTY)
        LEXIST=.TRUE.
        DO J=1,NA
         IF (MAP2(J)/=MAP3(J)) LEXIST=.FALSE.
        ENDDO
        IF (LEXIST) NEXIST1=NEXIST1+1
        IF (NEXIST1 > 1) EXIT
        IF (LEXIST) THEN
         i0(((((m(p1,N)-1)*NA+m(p2,N)-1)*NA+m(h3,N)-1)*NA+m(h4,N)))=i0(((((m(p1,N)-1)*NA+m(p2,N)-1)*NA+m(h3,N)-1)*NA+m(h4,N)))&
          -s1(((m(pa,I)-1)*NA+m(ha,I)),I)*E(I,N)&
          *DFLOAT(INTERNAL3(N)*IPRTY*INTERNAL3(I)*JPRTY)
        ENDIF
       ENDDO
      ENDDO
     ENDIF
     IF (LDOUBLES) THEN
      NEXIST2=0
      DO pa=NO+1,NA
       DO pb=NO+1,NA
        IF (LPERMUTATION.AND.(pa>=pb)) CYCLE
        IF (pa==pb) CYCLE
        DO ha=1,NO
         DO hb=1,NO
          IF (LPERMUTATION.AND.(ha>=hb)) CYCLE
          IF (ha==hb) CYCLE
          MAP3(:)=m(:,I)
          CALL PRTY2(m(PA,I),m(PB,I),m(HA,I),m(HB,I),NO,NA,MAP3,JPRTY)
          LEXIST=.TRUE.
          DO J=1,NA
           IF (MAP2(J)/=MAP3(J)) LEXIST=.FALSE.
          ENDDO
          IF (LEXIST) NEXIST2=NEXIST2+1
          IF (NEXIST2 > 1) EXIT
          IF (LEXIST) THEN
           i0(((((m(p1,N)-1)*NA+m(p2,N)-1)*NA+m(h3,N)-1)*NA+m(h4,N)))=i0(((((m(p1,N)-1)*NA+m(p2,N)-1)*NA+m(h3,N)-1)*NA+m(h4,N)))&
            -s2(((((m(pa,I)-1)*NA+m(pb,I)-1)*NA+m(ha,I)-1)*NA+m(hb,I)),I)*E(I,N)&
            *DFLOAT(INTERNAL3(N)*IPRTY*INTERNAL3(I)*JPRTY)
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDIF
     IF (LPERMUTATION.AND.(NEXIST1+NEXIST2 > 1)) THEN
      WRITE(6,'(A)') 'ERROR: AN EXTERNAL DETERMINANT COULD NOT BE MAPPED ON A UNIQUE T-AMPLITUDE 2'
      STOP
     ENDIF
    ENDDO
    t(((((m(p1,N)-1)*NA+m(p2,N)-1)*NA+m(h3,N)-1)*NA+m(h4,N)),N)=t(((((m(p1,N)-1)*NA+m(p2,N)-1)*NA+m(h3,N)-1)*NA+m(h4,N)),N)&
     -i0(((((m(p1,N)-1)*NA+m(p2,N)-1)*NA+m(h3,N)-1)*NA+m(h4,N)))/(EPS(m(p1,N))+EPS(m(p2,N))-EPS(m(h3,N))-EPS(m(h4,N)))
    R=R+i0(((((m(p1,N)-1)*NA+m(p2,N)-1)*NA+m(h3,N)-1)*NA+m(h4,N)))**2
   ENDIF
   ENDDO
   ENDDO
   ENDDO
   ENDDO

   ELSE IF (ORDER==1) THEN

   DO p1=NO+1,NA
   DO h2=1,NO
   IF (DABS(EPS(m(p1,N))-EPS(m(h2,N))) < 1.0D-9) THEN
    J=-1
    DO I=1,NDEGEN
     IF ((INTERNAL1(I,N)==1).AND.(INTERNAL2(I,N)==(m(P1,N)-1)*NA+m(h2,N))) J=I
    ENDDO
    IF (J /= -1) THEN
     i0(((m(p1,N)-1)*NA+m(h2,N)))=S(J,N)
     t(((m(p1,N)-1)*NA+m(h2,N)),N)=t(((m(p1,N)-1)*NA+m(h2,N)),N)-i0(((m(p1,N)-1)*NA+m(h2,N))) 
     R=R+i0(((m(p1,N)-1)*NA+m(h2,N)))**2
    ENDIF
   ELSE
    MAP2(:)=m(:,N)
    CALL PRTY1(m(P1,N),m(H2,N),NO,NA,MAP2,IPRTY)
    DO I=1,NDEGEN
     IF (LSINGLES) THEN
      NEXIST1=0
      DO pa=NO+1,NA
       DO ha=1,NO
        MAP3(:)=m(:,I)
        CALL PRTY1(m(PA,I),m(HA,I),NO,NA,MAP3,JPRTY)
        LEXIST=.TRUE.
        DO J=1,NA
         IF (MAP2(J)/=MAP3(J)) LEXIST=.FALSE.
        ENDDO
        IF (LEXIST) NEXIST1=NEXIST1+1
        IF (NEXIST1 > 1) EXIT
        IF (LEXIST) THEN
         i0(((m(p1,N)-1)*NA+m(h2,N)))=i0(((m(p1,N)-1)*NA+m(h2,N)))&
          -s1(((m(pa,I)-1)*NA+m(ha,I)),I)*E(I,N)&
          *DFLOAT(INTERNAL3(N)*IPRTY*INTERNAL3(I)*JPRTY)
        ENDIF
       ENDDO
      ENDDO
     ENDIF
     IF (LDOUBLES) THEN
      NEXIST2=0
      DO pa=NO+1,NA
       DO pb=NO+1,NA
        IF (LPERMUTATION.AND.(pa>=pb)) CYCLE
        IF (pa==pb) CYCLE
        DO ha=1,NO
         DO hb=1,NO
          IF (LPERMUTATION.AND.(ha>=hb)) CYCLE
          IF (ha==hb) CYCLE
          MAP3(:)=m(:,I)
          CALL PRTY2(m(PA,I),m(PB,I),m(HA,I),m(HB,I),NO,NA,MAP3,JPRTY)
          LEXIST=.TRUE.
          DO J=1,NA
           IF (MAP2(J)/=MAP3(J)) LEXIST=.FALSE.
          ENDDO
          IF (LEXIST) NEXIST2=NEXIST2+1
          IF (NEXIST2 > 1) EXIT
          IF (LEXIST) THEN
           i0(((m(p1,N)-1)*NA+m(h2,N)))=i0(((m(p1,N)-1)*NA+m(h2,N)))&
            -s2(((((m(pa,I)-1)*NA+m(pb,I)-1)*NA+m(ha,I)-1)*NA+m(hb,I)),I)*E(I,N)&
            *DFLOAT(INTERNAL3(N)*IPRTY*INTERNAL3(I)*JPRTY)
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO 
     ENDIF
     IF (LPERMUTATION.AND.(NEXIST1+NEXIST2 > 1)) THEN
      WRITE(6,'(A)') 'ERROR: AN EXTERNAL DETERMINANT COULD NOT BE MAPPED ON A UNIQUE T-AMPLITUDE 3'
      STOP
     ENDIF
    ENDDO
    t(((m(p1,N)-1)*NA+m(h2,N)),N)=t(((m(p1,N)-1)*NA+m(h2,N)),N)-i0(((m(p1,N)-1)*NA+m(h2,N)))/(EPS(m(p1,N))-EPS(m(h2,N)))
    R=R+i0(((m(p1,N)-1)*NA+m(h2,N)))**2
   ENDIF
   ENDDO
   ENDDO
   ENDIF
   RETURN
END SUBROUTINE



SUBROUTINE ISORT(IARRAY,N)

   IMPLICIT NONE
   INTEGER :: N
   INTEGER :: IARRAY(*)
   INTEGER,ALLOCATABLE :: ISUB(:)
   INTEGER :: I,J,K
   INTEGER :: LARGE

   ALLOCATE(ISUB(N))
   DO I=1,N
    LARGE=10**8
    DO J=1,N
     IF (IARRAY(J) < LARGE) THEN
      K=J
      LARGE=IARRAY(K)
     ENDIF
    ENDDO
    ISUB(I)=IARRAY(K)
    IARRAY(K)=2*10**8
   ENDDO
   IARRAY(1:N)=ISUB(1:N)
   DEALLOCATE(ISUB)
   RETURN
END SUBROUTINE
        
         
         
SUBROUTINE ISORT2(IARRAY,JARRAY,N)
         
   IMPLICIT NONE
   INTEGER :: N
   INTEGER :: IARRAY(*)
   INTEGER :: JARRAY(*)
   INTEGER,ALLOCATABLE :: ISUB(:),JSUB(:)
   INTEGER :: I,J,K
   INTEGER :: LARGE
        
   ALLOCATE(ISUB(N),JSUB(N))
   DO I=1,N
    LARGE=10**8
    DO J=1,N
     IF (IARRAY(J) < LARGE) THEN
      K=J
      LARGE=IARRAY(K)
     ENDIF
    ENDDO 
    ISUB(I)=IARRAY(K)
    JSUB(I)=JARRAY(K)
    IARRAY(K)=2*10**8
   ENDDO 
   IARRAY(1:N)=ISUB(1:N)
   JARRAY(1:N)=JSUB(1:N)
   DEALLOCATE(ISUB,JSUB)
   RETURN 
END SUBROUTINE

        

SUBROUTINE ISORT3(IARRAY,JARRAY,KARRAY,N)
         
   IMPLICIT NONE
   INTEGER :: N
   INTEGER :: IARRAY(*)
   INTEGER :: JARRAY(*)
   INTEGER :: KARRAY(*)
   INTEGER,ALLOCATABLE :: ISUB(:),JSUB(:),KSUB(:)
   INTEGER :: I,J,K
   INTEGER :: LARGE
         
   ALLOCATE(ISUB(N),JSUB(N),KSUB(N))
   DO I=1,N
    LARGE=10**8 
    DO J=1,N
     IF (IARRAY(J) < LARGE) THEN
      K=J
      LARGE=IARRAY(K)
     ENDIF
    ENDDO 
    ISUB(I)=IARRAY(K)
    JSUB(I)=JARRAY(K)
    KSUB(I)=KARRAY(K)
    IARRAY(K)=2*10**8
   ENDDO 
   IARRAY(1:N)=ISUB(1:N)
   JARRAY(1:N)=JSUB(1:N)
   KARRAY(1:N)=KSUB(1:N)
   DEALLOCATE(ISUB,JSUB,KSUB)
   RETURN 
END SUBROUTINE
 


SUBROUTINE PRTY(NO,NA,MAP,KPRTY)

   INTEGER :: NO
   INTEGER :: NA
   INTEGER :: MAP(NA)
   INTEGER :: IPRTY,JPRTY
   INTEGER :: KPRTY
   INTEGER :: I,J,K

   IPRTY=0
   DO I=1,NO
    K=0
    DO J=1,NO
     IF (MAP(J) < I) THEN
      K=K+1
     ELSE IF (MAP(J)==I) THEN
      EXIT
     ELSE IF (MAP(J) > I) THEN
      IPRTY=IPRTY+K
!write(*,*) 'occ ',i,'is added with ',k,' preceding'
      EXIT
     ENDIF
    ENDDO
   ENDDO
   JPRTY=0
   DO I=NA,NO+1,-1
    K=0
    DO J=1,NO
     IF (MAP(J) <= NO) THEN
      K=K+1
     ELSE IF (MAP(J)==I) THEN
      JPRTY=JPRTY+K
!write(*,*) 'vir ',i,'is added with ',k,' preceding'
      EXIT
     ELSE IF (MAP(J) > I) THEN
      EXIT
     ENDIF
    ENDDO
   ENDDO
!write(*,*) 'occ,vir parities = ',IPRTY,JPRTY
   KPRTY=(-1)**(IPRTY+JPRTY)

   RETURN
END SUBROUTINE



SUBROUTINE PRTY1(P,H,NO,NA,MAP,KPRTY)

   INTEGER :: P,H
   INTEGER :: NO
   INTEGER :: NA
   INTEGER :: MAP(NA)
   INTEGER :: IPRTY,JPRTY
   INTEGER :: KPRTY
   INTEGER :: PA,HA
   INTEGER :: I,J,K

!write(*,'(A,I3,A,I3)') 'Excitation ',H,' -> ',P
!write(*,'(A,100I3:)') 'MAP =',(MAP(I),I=1,NA)
   HA=-1
   K=0
   DO I=1,NO
    IF (MAP(I) < H) THEN
     K=K+1
    ELSE IF (MAP(I)==H) THEN
     HA=I
     IPRTY=(-1)**K
!write(*,*) H,' ANNIHILATED WITH PARITY ',K
     EXIT
    ELSE IF (MAP(I) > H) THEN
     KPRTY=0
     RETURN
    ENDIF
   ENDDO
   IF (HA==-1) THEN
    WRITE(6,'(A)') 'ERROR IN PRTY1: HA NOT DETERMINED'
    STOP
   ENDIF
   K=0
   DO I=1,NO
    IF (MAP(I) > P) THEN
     JPRTY=(-1)**K
!write(*,*) P,' CREATED WITH PARITY ',K
     EXIT
    ELSE IF (MAP(I) < P) THEN
     IF (I /= HA) K=K+1
     IF (I==NO) THEN
      JPRTY=(-1)**K
!write(*,*) P,' CREATED WITH PARITY ',K
      EXIT
     ENDIF
    ELSE IF (MAP(I)==P) THEN
     KPRTY=0
     RETURN
    ENDIF
   ENDDO
   PA=-1
   DO I=NO+1,NA
    IF (MAP(I)==P) PA=I
   ENDDO
   IF (PA==-1) THEN
    WRITE(6,'(A)') 'ERROR IN PRTY1: PA NOT DETERMINED'
    STOP
   ENDIF
   MAP(HA)=P
   MAP(PA)=H
   CALL ISORT(MAP,NO)
   CALL ISORT(MAP(NO+1),NA-NO)
!write(*,'(A,100I3:)') 'SORTED MAP =',(MAP(I),I=1,NA)
   KPRTY=IPRTY*JPRTY
!write(*,'(A,I3)') 'PARITY OF THIS 1-E EXCITATION = ',KPRTY
   RETURN
END SUBROUTINE



SUBROUTINE PRTY2(P1,P2,H1,H2,NO,NA,MAP,KPRTY)

   INTEGER :: P1,P2,H1,H2
   INTEGER :: NO
   INTEGER :: NA
   INTEGER :: MAP(NA)
   INTEGER :: IPRTY,JPRTY
   INTEGER :: KPRTY
   INTEGER :: PA,PB,HA,HB
   INTEGER :: I,J,K
   LOGICAL :: LCHECK

!write(*,*) 'Excitation ',H1,' <',H2,' -> ',P1,' <',P2
!write(*,'(A,100I3:)') 'MAP =',(MAP(I),I=1,NA)
   HA=-1
   K=0
   DO I=1,NO
    IF (MAP(I) < H1) THEN
     K=K+1
    ELSE IF (MAP(I)==H1) THEN
     HA=I
     IPRTY=(-1)**K
!write(*,*) H1,' ANNIHILATED WITH PARITY ',K
     EXIT
    ELSE IF (MAP(I) > H1) THEN
     KPRTY=0
     RETURN
    ENDIF
   ENDDO
   IF (HA==-1) THEN
    WRITE(6,'(A)') 'ERROR IN PRTY2: HA NOT DETERMINED'
    STOP
   ENDIF
   HB=-1
   K=0
   DO I=1,NO
    IF (MAP(I) < H2) THEN
     IF (I /= HA) K=K+1
    ELSE IF (MAP(I)==H2) THEN
     HB=I
     IPRTY=IPRTY*(-1)**K
!write(*,*) H2,' ANNIHILATED WITH PARITY ',K
     LCHECK=.TRUE.
     EXIT
    ELSE IF (MAP(I) > H2) THEN
     KPRTY=0
     RETURN
    ENDIF
   ENDDO
   IF (HB==-1) THEN
    WRITE(6,'(A)') 'ERROR IN PRTY2: HB NOT DETERMINED'
    STOP
   ENDIF
   LCHECK=.FALSE.
   K=0
   DO I=1,NO
    IF (MAP(I) > P2) THEN
     JPRTY=(-1)**K
!write(*,*) P2,' CREATED WITH PARITY ',K
     LCHECK=.TRUE.
     EXIT
    ELSE IF (MAP(I) < P2) THEN
     IF ((I /= HA).AND.(I /= HB)) K=K+1
     IF (I==NO) THEN
      JPRTY=(-1)**K
!write(*,*) P2,' CREATED WITH PARITY ',K
      LCHECK=.TRUE.
      EXIT
     ENDIF
    ELSE IF (MAP(I)==P2) THEN
     KPRTY=0
     RETURN
    ENDIF
   ENDDO
   IF (.NOT.LCHECK) THEN
    WRITE(6,'(I3,A)') P2,' NOT CREATED: ERROR'
    STOP
   ENDIF
   PB=-1
   DO I=NO+1,NA
    IF (MAP(I)==P2) PB=I
   ENDDO
   IF (PB==-1) THEN
    WRITE(6,'(A)') 'ERROR IN PRTY2: PB NOT DETERMINED'
    STOP
   ENDIF
   LCHECK=.FALSE.
   K=0
   DO I=1,NO
    IF (MAP(I) > P1) THEN
     JPRTY=JPRTY*(-1)**K
!write(*,*) P1,' CREATED WITH PARITY ',K
     LCHECK=.TRUE.
     EXIT
    ELSE IF (MAP(I) < P1) THEN
     IF ((I /= HA).AND.(I /= HB)) K=K+1
     IF (I==NO) THEN
      JPRTY=JPRTY*(-1)**K
!write(*,*) P1,' CREATED WITH PARITY ',K
      LCHECK=.TRUE.
      EXIT
     ENDIF
    ELSE IF (MAP(I)==P1) THEN
     KPRTY=0
     RETURN
    ENDIF
   ENDDO
   IF (.NOT.LCHECK) THEN
    WRITE(6,'(I3,A)') P1,' NOT CREATED: ERROR'
    STOP
   ENDIF
   PA=-1
   DO I=NO+1,NA
    IF (MAP(I)==P1) PA=I
   ENDDO
   IF (PA==-1) THEN
    WRITE(6,'(A)') 'ERROR IN PRTY2: PA NOT DETERMINED'
    STOP
   ENDIF
   MAP(HA)=P1
   MAP(HB)=P2
   MAP(PA)=H1
   MAP(PB)=H2
   CALL ISORT(MAP,NO)
   CALL ISORT(MAP(NO+1),NA-NO)
!write(*,'(A,100I3:)') 'SORTED MAP =',(MAP(I),I=1,NA)
   KPRTY=IPRTY*JPRTY
!write(*,'(A,I3)') 'PARITY OF THIS 2-E EXCITATION = ',KPRTY
   RETURN
END SUBROUTINE
