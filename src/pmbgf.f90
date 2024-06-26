SUBROUTINE HIGHORDER_GF(OMEGA,ORDER)
! CONSTRUCT ONE-PARTICLE MBGF FOR A GIVEN DETERMINANTAL WAVE FUNCTION IN FILE 
! SEE EQ.(13) OF PICKUP & GOSCINSKI, MOL.PHYS. 26, 1013 (1973) -- Dec, 2012 This paper looks wrong.
! Dec, 2012, SEE JORGENSEN & SIMONS EQ. (6.67)

   USE CONSTANTS
   USE CONTROL
   USE GRADIENT
   USE STRUCTURE
   USE INTEGRAL
   USE BASISSET
   USE FULLCI

   IMPLICIT NONE
   DOUBLE PRECISION :: OMEGA
   INTEGER :: ORDER,BRAORDER,KETORDER
   REAL :: MEM,ICPUS,ICPUE
   INTEGER :: I,PA,PB,QA,QB,RA,RB
   INTEGER :: P,Q
   INTEGER,ALLOCATABLE :: IPEA(:)
   DOUBLE PRECISION :: HA,HB
   DOUBLE PRECISION,ALLOCATABLE :: V_IPEA(:,:,:,:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: W_IPEA(:,:,:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: S_IPEA(:,:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: T_IPEA(:,:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: SIGMA(:,:,:),TIGMA(:,:),DENOM(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: VEC0(:,:),VEC1(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: BRA(:,:),KET(:,:)
   INTEGER :: IORDER,JORDER
   INTEGER,ALLOCATABLE :: RANK(:,:)
   DOUBLE PRECISION :: NONHER
   INTEGER :: IRANK,JRANK,MAXRANK
! debug
!integer :: moi,moj,moa,mob
!integer :: hn,pn,h1(10),p1(10),hab(10),pab(10)
!double precision :: gg,gg1,gg2,gg3,gg4,sig
!integer :: moia,moja,moaa,moba
! end debug

!  CALL PCPU_TIME(ICPUS)
!  IF (ORDER < 2) CALL PABORT('ORDER IS TOO LOW')
   IF (ICORE /= 0) CALL WARNING('FROZEN CORE IS DANGEROUS FOR MBGF')

   WRITE(6,'(A,F20.15,A)') 'SELF-ENERGY MATRIX WILL BE COMPUTED AT OMEGA=',OMEGA,' HARTREE'

   MEM=0.0
   IF (MEM > 1000000.0) THEN
    WRITE(6,'(A,F7.1,A)') 'ESTIMATED MEMORY USAGE WILL BE ',MEM/1000000.0,' MB'
   ELSE IF (MEM > 1000.0) THEN
    WRITE(6,'(A,F7.1,A)') 'ESTIMATED MEMORY USAGE WILL BE ',MEM/1000.0,' KB'
   ELSE
    WRITE(6,'(A,F7.1,A)') 'ESTIMATED MEMORY USAGE WILL BE ',MEM,' B'
   ENDIF
   IF (MEM > DOPTN(28)*1000000.0) CALL PABORT('OUT OF MEMORY') 

   ALLOCATE(V_IPEA(NCF,IP_NCF+EA_NCF,NCF,IP_NCF+EA_NCF,0:ORDER-1,0:ORDER-1))
   ALLOCATE(W_IPEA(NCF,IP_NCF+EA_NCF,NCF,IP_NCF+EA_NCF,1:ORDER))
   ALLOCATE(S_IPEA(NCF,IP_NCF+EA_NCF,NCF,IP_NCF+EA_NCF))
   ALLOCATE(T_IPEA(NCF,IP_NCF+EA_NCF,NCF,IP_NCF+EA_NCF))
   ALLOCATE(SIGMA(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,0:ORDER))
   ALLOCATE(TIGMA(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE))
   ALLOCATE(DENOM(NCF,IP_NCF+EA_NCF))
   ALLOCATE(RANK(NCF,IP_NCF+EA_NCF))
   ALLOCATE(IPEA(1:IALL(0,0,0)-IVIRTCORE))

   ! DISK AND MEMORY SPACES
   OPEN(50,FILE=TRIM(COPTN(1))//'.gfa',FORM='UNFORMATTED')
   OPEN(51,FILE=TRIM(COPTN(1))//'.gfb',FORM='UNFORMATTED')
   OPEN(52,FILE=TRIM(COPTN(1))//'.gfc',FORM='UNFORMATTED')
   ALLOCATE(VEC0(NCF,NCF),VEC1(NCF,NCF))
   ALLOCATE(BRA(NCF,NCF),KET(NCF,NCF))

   ! DENOMINATOR AND RANK
   MAXRANK=0
   DO PA=1,NCF
    HA=0.0D0
    DO I=1,IOCC
     IF (.NOT.BTEST(CFHALF(PA),I-1)) HA=HA+EPSILON(I,0,0,0)
    ENDDO
    DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
     IF (BTEST(CFHALF(PA),I-1)) HA=HA-EPSILON(I,0,0,0)
    ENDDO
    DO PB=1,IP_NCF
     HB=0.0D0
     DO I=1,IOCC
      IF (.NOT.BTEST(IP_CFHALF(PB),I-1)) HB=HB+EPSILON(I,0,0,0)
     ENDDO
     DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
      IF (BTEST(IP_CFHALF(PB),I-1)) HB=HB-EPSILON(I,0,0,0)
     ENDDO
     DENOM(PA,PB)=HA+HB
     RANK(PA,PB)=NORDER(PA)+IP_NORDER(PB)
     IF (RANK(PA,PB) > MAXRANK) MAXRANK=RANK(PA,PB)
    ENDDO
   ENDDO
   DO PA=1,NCF
    HA=0.0D0
    DO I=1,IOCC
     IF (.NOT.BTEST(CFHALF(PA),I-1)) HA=HA-EPSILON(I,0,0,0)
    ENDDO
    DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
     IF (BTEST(CFHALF(PA),I-1)) HA=HA+EPSILON(I,0,0,0)
    ENDDO
    DO PB=1,EA_NCF
     HB=0.0D0
     DO I=1,IOCC
      IF (.NOT.BTEST(EA_CFHALF(PB),I-1)) HB=HB-EPSILON(I,0,0,0)
     ENDDO
     DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
      IF (BTEST(EA_CFHALF(PB),I-1)) HB=HB+EPSILON(I,0,0,0)
     ENDDO
     DENOM(PA,IP_NCF+PB)=HA+HB
     RANK(PA,IP_NCF+PB)=NORDER(PA)+EA_NORDER(PB)
     IF (RANK(PA,IP_NCF+PB) > MAXRANK) MAXRANK=RANK(PA,IP_NCF+PB)
    ENDDO
   ENDDO

   ! MAP ORBITALS TO DETERMINANTS
   DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
    IF (P <= IOCC) THEN
     DO PA=1,IP_NCF
      IF ((IP_NORDER(PA)==1).AND.(.NOT.BTEST(IP_CFHALF(PA),P-1))) IPEA(P)=PA
     ENDDO
    ELSE
     DO PA=1,EA_NCF
      IF ((EA_NORDER(PA)==1).AND.(BTEST(EA_CFHALF(PA),P-1))) IPEA(P)=PA+IP_NCF
     ENDDO
    ENDIF
   ENDDO

   ! CONSTRUCT SIGMA(:,:,0)
   SIGMA=0.0D0
   DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
    SIGMA(P,P,0)=EPSILON(P,0,0,0)
   ENDDO
   WRITE(6,'(I3,A)') 0,'TH-ORDER SIGMA'
   CALL DUMP5(SIGMA(:,:,0),IALL(0,0,0)-IVIRTCORE)

   DO IORDER=ORDER,ORDER
    DO BRAORDER=0,IORDER-1
     DO KETORDER=0,IORDER-1-BRAORDER

      ! CONSTRUCT (BRA|V KET) SUPEROPERATOR MATRIX
      WRITE(6,'(/,A,I1,A,I1,A)') 'CONSTRUCTING (',BRAORDER,'|V',KETORDER,') SUPEROPERATOR MATRIX'

      ! Note added on April 27, 2016. Superoperator matrix of H
      ! is not Hermitian. In particular (5|H1) = 0 but (1|H5) <> 0.
      ! It is said that improvement of wave function to MP1 makes
      ! (5|H1) and (1|H5) the same (if it's zero is not clear).

      OPEN(97,FILE=TRIM(COPTN(1))//'.gf0',FORM='UNFORMATTED')
      REWIND(97)
      DO I=0,BRAORDER
       READ(97) BRA
      ENDDO
      WRITE(6,'(/,I3,A,A)') BRAORDER,'TH-ORDER BRA REFERENCE WAVE FUNCTION READ FROM ',TRIM(COPTN(1))//'.gf0'
      REWIND(97)
      DO I=0,KETORDER
       READ(97) KET
      ENDDO
      WRITE(6,'(I3,A,A)') KETORDER,'TH-ORDER KET REFERENCE WAVE FUNCTION READ FROM ',TRIM(COPTN(1))//'.gf0'
      CLOSE(97)

      REWIND(50)
      WRITE(50) KET
!write(*,*) 'Ket'
!call dump5(ket,ncf)
!write(*,*) 'Bra'
!call dump5(bra,ncf)
      V_IPEA(:,:,:,:,BRAORDER,KETORDER)=0.0D0

      DO QA=1,NCF
       DO QB=1,IP_NCF+EA_NCF
        DO PA=1,NCF
         DO PB=1,IP_NCF+EA_NCF
   
          VEC1=0.0D0
   
          ! p+ q H
          CALL HAMILTONIAN_PRODUCT(50,51,0,2*MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE))
          CALL IONIZE(51,52,QA,QB,.FALSE.,0)
          CALL IONIZE(52,51,PA,PB,.TRUE. ,1)
          REWIND(51)
          READ(51) VEC0
          VEC1=VEC1+VEC0

!if ((braorder==0).and.(ketorder==1).and.(rank(pa,pb)==1).and.(rank(qa,qb)==1)) then
!write(*,*) "++++++++++++++++++"
!write(*,*) pa,pb,qa,qb
!write(*,*) "++++++++++++++++++"
!call dump5(vec0,ncf)
!endif
   
          ! - p+ H q
          CALL IONIZE(50,51,QA,QB,.FALSE.,0)
          CALL IP_HAMILTONIAN_PRODUCT(51,52,0,MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE) &
                                             +MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE+1))
          CALL IONIZE(52,51,PA,PB,.TRUE. ,1)
          REWIND(51)
          READ(51) VEC0
          VEC1=VEC1-VEC0

!if ((braorder==0).and.(ketorder==1).and.(rank(pa,pb)==1).and.(rank(qa,qb)==1)) then
!call dump5(vec0,ncf)
!endif
   
          ! q H p+
          CALL IONIZE(50,51,PA,PB,.TRUE. ,0)
          CALL EA_HAMILTONIAN_PRODUCT(51,52,0,MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE) &
                                             +MIN(IOCC-ICORE+1,IALL(0,0,0)-IOCC-IVIRTCORE))
          CALL IONIZE(52,51,QA,QB,.FALSE.,-1)
          REWIND(51)
          READ(51) VEC0
          VEC1=VEC1+VEC0

!if ((braorder==0).and.(ketorder==1).and.(rank(pa,pb)==1).and.(rank(qa,qb)==1)) then
!call dump5(vec0,ncf)
!endif
   
          ! - H q p+
          CALL IONIZE(50,51,PA,PB,.TRUE. ,0)
          CALL IONIZE(51,52,QA,QB,.FALSE.,-1)
          CALL HAMILTONIAN_PRODUCT(52,51,0,2*MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE))
          REWIND(51)
          READ(51) VEC0
          VEC1=VEC1-VEC0

!if ((braorder==0).and.(ketorder==1).and.(rank(pa,pb)==1).and.(rank(qa,qb)==1)) then
!call dump5(vec0,ncf)
!write(*,*) "sum so far"
!call dump5(vec1,ncf)
!endif
   
          ! - p+ q F
          CALL H0_PRODUCT(50,51,0)
          CALL IONIZE(51,52,QA,QB,.FALSE.,0)
          CALL IONIZE(52,51,PA,PB,.TRUE. ,1)
          REWIND(51)
          READ(51) VEC0
          VEC1=VEC1-VEC0
   
          ! p+ F q
          CALL IONIZE(50,51,QA,QB,.FALSE.,0)
          CALL H0_PRODUCT(51,52,1)
          CALL IONIZE(52,51,PA,PB,.TRUE. ,1)
          REWIND(51)
          READ(51) VEC0
          VEC1=VEC1+VEC0
   
          ! - q F p+
          CALL IONIZE(50,51,PA,PB,.TRUE. ,0)
          CALL H0_PRODUCT(51,52,-1)
          CALL IONIZE(52,51,QA,QB,.FALSE.,-1)
          REWIND(51)
          READ(51) VEC0
          VEC1=VEC1-VEC0
   
          ! F q p+
          CALL IONIZE(50,51,PA,PB,.TRUE. ,0)
          CALL IONIZE(51,52,QA,QB,.FALSE.,-1)
          CALL H0_PRODUCT(52,51,0)
          REWIND(51)
          READ(51) VEC0
          VEC1=VEC1+VEC0
   
          DO RA=1,NCF
           DO RB=1,NCF
            V_IPEA(PA,PB,QA,QB,BRAORDER,KETORDER)= &
            V_IPEA(PA,PB,QA,QB,BRAORDER,KETORDER)+BRA(RA,RB)*VEC1(RA,RB)
           ENDDO
          ENDDO

!if ((braorder==0).and.(ketorder==1).and.(rank(pa,pb)==1).and.(rank(qa,qb)==1)) &
!write(*,'(4i3,f20.10,A)') pa,pb,qa,qb,v_ipea(pa,pb,qa,qb,braorder,ketorder),'XXXXXXXXXXX'
   
         ENDDO
        ENDDO
       ENDDO
      ENDDO

     ENDDO
    ENDDO
   ENDDO

   W_IPEA=0.0D0
   DO IORDER=1,ORDER
    DO BRAORDER=0,IORDER-1
     W_IPEA(:,:,:,:,IORDER)=W_IPEA(:,:,:,:,IORDER)+V_IPEA(:,:,:,:,BRAORDER,IORDER-BRAORDER-1)
    ENDDO
   ENDDO

!write(*,*) ' +++++++++++ <0|(f1|V f3)|0>, <0|(f3|V f1)|0> '
!DO QA=1,NCF
! DO QB=1,IP_NCF+EA_NCF
!  DO PA=1,NCF
!   DO PB=1,IP_NCF+EA_NCF
!    if ((rank(qa,qb)==2).and.(rank(pa,pb)==1)) &
!write(*,'(4i3,2f20.10)') pa,pb,qa,qb,&
!v_ipea(pa,pb,qa,qb,0,0),v_ipea(qa,qb,pa,pb,0,0)
!   ENDDO
!  ENDDO
! ENDDO
!ENDDO
!write(*,*) ' +++++++++++ <0|(f1|V f3)|1>, <1|(f3|V f1)|0> '
!DO QA=1,NCF
! DO QB=1,IP_NCF+EA_NCF
!  DO PA=1,NCF
!   DO PB=1,IP_NCF+EA_NCF
!    if ((rank(qa,qb)==2).and.(rank(pa,pb)==1)) &
!write(*,'(4i3,2f20.10)') pa,pb,qa,qb,&
!v_ipea(pa,pb,qa,qb,0,1)+v_ipea(pa,pb,qa,qb,1,0), &
!v_ipea(qa,qb,pa,pb,0,1)+v_ipea(qa,qb,pa,pb,1,0)
!   ENDDO
!  ENDDO
! ENDDO
!ENDDO
!  WRITE(6,'(/,A)') 'MAXIMUM NON-HERMITIVITY'
!  DO IORDER=ORDER,ORDER
!   DO BRAORDER=0,IORDER-1
!    DO KETORDER=0,IORDER-1-BRAORDER
!     NONHER=0.0D0
!     DO QA=1,NCF
!      DO QB=1,IP_NCF+EA_NCF
!       DO PA=1,NCF
!        DO PB=1,IP_NCF+EA_NCF
!    IF (DABS(V_IPEA(PA,PB,QA,QB,BRAORDER,KETORDER)-V_IPEA(QA,QB,PA,PB,BRAORDER,KETORDER)) > NONHER) &
!    NONHER=DABS(V_IPEA(PA,PB,QA,QB,BRAORDER,KETORDER)-V_IPEA(QA,QB,PA,PB,BRAORDER,KETORDER))
!        ENDDO
!       ENDDO
!      ENDDO
!     ENDDO
!     WRITE(6,'(A,I1,A,I1,A,I1,A,I1,A,F20.10)') '  (',BRAORDER,'|V',KETORDER,')-(',KETORDER,'|V',BRAORDER,')T =',NONHER
!    ENDDO
!   ENDDO
!  ENDDO

!  WRITE(*,*) "FORCING (f1|V f5)=0 !!!!"
!  DO QA=1,NCF
!   DO QB=1,IP_NCF+EA_NCF
!    DO PA=1,NCF
!     DO PB=1,IP_NCF+EA_NCF
!      IF (RANK(QA,QB)-RANK(PA,PB)==2) THEN
!       V_IPEA(PA,PB,QA,QB,0,0)=0.0D0
!      ENDIF
!     ENDDO
!    ENDDO
!   ENDDO
!  ENDDO

!  DO IORDER=ORDER,ORDER
!   DO BRAORDER=0,IORDER-1
!    DO KETORDER=0,IORDER-1-BRAORDER
!     WRITE(6,'(/,A,I1,A,I1,A)') 'NON-HERMITICITY OF (',BRAORDER,'|V',KETORDER,') SUPEROPERATOR MATRIX'
!  DO IRANK=1,MAXRANK
!   DO JRANK=1,MAXRANK
!    NONHER=-1.0D0
!    DO QA=1,NCF
!     DO QB=1,IP_NCF+EA_NCF
!      IF (RANK(QA,QB)==IRANK) THEN
!      DO PA=1,NCF
!       DO PB=1,IP_NCF+EA_NCF
!        IF (RANK(PA,PB)==JRANK) THEN
!         IF (DABS(V_IPEA(PA,PB,QA,QB,BRAORDER,KETORDER)-V_IPEA(QA,QB,PA,PB,KETORDER,BRAORDER)) > NONHER) &
!         NONHER=DABS(V_IPEA(PA,PB,QA,QB,BRAORDER,KETORDER)-V_IPEA(QA,QB,PA,PB,KETORDER,BRAORDER))
!        ENDIF
!       ENDDO
!      ENDDO
!      ENDIF
!     ENDDO
!    ENDDO
!    WRITE(6,'(A,I1,A,I1,A,I1,A,I1,A,F20.10)') &
!    '  (f',2*IRANK-1,'|V f',2*JRANK-1,')-(f',2*JRANK-1,'|V f',2*IRANK-1,') =',NONHER
!   ENDDO
!  ENDDO
!    ENDDO
!   ENDDO
!  ENDDO

! debug
!write(*,*) '************ H13^(1)'
!DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
! QA=1
! QB=IPEA(Q)
! write(*,*) '================ ',q
! sig=0.0d0
! DO PA=1,NCF
!  DO PB=1,IP_NCF+EA_NCF
!   IF (RANK(PA,PB)==2) THEN
!    hn=0
!    pn=0
!    DO I=1,IOCC
!     IF (.NOT.BTEST(CFHALF(PA),I-1)) then
!      hn=hn+1
!      h1(hn)=i
!      hab(hn)=0
!     endif
!    ENDDO
!    DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
!     IF (BTEST(CFHALF(PA),I-1)) then
!      pn=pn+1
!      p1(pn)=i
!      pab(pn)=0
!     endif
!    ENDDO
!    if (pb <= ip_ncf) then
!     DO I=1,IOCC
!      IF (.NOT.BTEST(IP_CFHALF(PB),I-1)) then
!       hn=hn+1
!       h1(hn)=i
!       hab(hn)=1
!      endif
!     ENDDO
!     DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
!      IF (BTEST(IP_CFHALF(PB),I-1)) then
!       pn=pn+1
!       p1(pn)=i
!       pab(pn)=1
!      endif
!     ENDDO
!    else
!     DO I=1,IOCC
!      IF (.NOT.BTEST(EA_CFHALF(PB-ip_ncf),I-1)) then
!       hn=hn+1
!       h1(hn)=i
!       hab(hn)=1
!      endif
!     ENDDO
!     DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
!      IF (BTEST(EA_CFHALF(PB-ip_ncf),I-1)) then
!       pn=pn+1
!       p1(pn)=i
!       pab(pn)=1
!      endif
!     ENDDO
!    endif
!    if (hn+pn /= 3) call pabort('something wrong')
!    if (hn==1) then
!     gg=0.0d0
!     if ((pab(1)==1).and.(hab(1)==pab(2))) gg=gg+g(q,p1(1),h1(1),p1(2))
!     if ((pab(2)==1).and.(hab(1)==pab(1))) gg=gg-g(q,p1(2),h1(1),p1(1))
!     write(*,'(4i3,3f20.10)') q,h1(1),p1(1),p1(2),gg,V_IPEA(PA,PB,QA,QB,0,0)
!     sig=sig+gg*gg/(omega-denom(pa,pb))
!    else if (hn==2) then
!     gg=0.0d0
!     if ((hab(1)==1).and.(pab(1)==hab(2))) gg=gg+g(q,h1(1),p1(1),h1(2))
!     if ((hab(2)==1).and.(pab(1)==hab(1))) gg=gg-g(q,h1(2),p1(1),h1(1))
!     write(*,'(4i3,3f20.10)') q,h1(1),h1(2),p1(1),gg,V_IPEA(PA,PB,QA,QB,0,0)
!     sig=sig+gg*gg/(omega-denom(pa,pb))
!    else
!     call pabort('bug')
!    endif
!   ENDIF
!  ENDDO
! ENDDO
! write(*,'(A,F20.10)') 'sig(2)=',sig
!ENDDO
!write(*,*) '************ H13^(2)'
!DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
! QA=1
! QB=IPEA(Q)
! write(*,*) '================ ',q
! sig=0.0d0
! DO PA=1,NCF
!  DO PB=1,IP_NCF+EA_NCF
!   IF (RANK(PA,PB)==2) THEN
!    hn=0
!    pn=0
!    DO I=1,IOCC
!     IF (.NOT.BTEST(CFHALF(PA),I-1)) then
!      hn=hn+1
!      h1(hn)=i
!      hab(hn)=0
!     endif
!    ENDDO
!    DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
!     IF (BTEST(CFHALF(PA),I-1)) then
!      pn=pn+1
!      p1(pn)=i
!      pab(pn)=0
!     endif
!    ENDDO
!    if (pb <= ip_ncf) then
!     DO I=1,IOCC
!      IF (.NOT.BTEST(IP_CFHALF(PB),I-1)) then
!       hn=hn+1
!       h1(hn)=i
!       hab(hn)=1
!      endif
!     ENDDO
!     DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
!      IF (BTEST(IP_CFHALF(PB),I-1)) then
!       pn=pn+1
!       p1(pn)=i
!       pab(pn)=1
!      endif
!     ENDDO
!    else
!     DO I=1,IOCC
!      IF (.NOT.BTEST(EA_CFHALF(PB-ip_ncf),I-1)) then
!       hn=hn+1
!       h1(hn)=i
!       hab(hn)=1
!      endif
!     ENDDO
!     DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
!      IF (BTEST(EA_CFHALF(PB-ip_ncf),I-1)) then
!       pn=pn+1
!       p1(pn)=i
!       pab(pn)=1
!      endif
!     ENDDO
!    endif
!    if (hn+pn /= 3) call pabort('something wrong')
!    if (hn==1) then
!     gg=0.0d0
!     do moia=0,1
!     do moi=icore+1,iocc
!     do moja=0,1
!     do moj=icore+1,iocc
!     gg1=0.0d0
!     if ((hab(1)==moia).and.(moja==1)) gg1=gg1+g(h1(1),moi,q,moj)
!     if ((hab(1)==moja).and.(moia==1)) gg1=gg1-g(h1(1),moj,q,moi)
!     gg2=0.0d0
!     if ((pab(1)==moia).and.(pab(2)==moja)) gg2=gg2+g(moi,p1(1),moj,p1(2))
!     if ((pab(1)==moja).and.(pab(2)==moia)) gg2=gg2-g(moi,p1(2),moj,p1(1))
!     gg=gg+0.5d0*gg1*gg2/(epsilon(moi,0)+epsilon(moj,0)-epsilon(p1(1),0)-epsilon(p1(2),0))
!     enddo
!     enddo
!     enddo
!     enddo
!     do moia=0,1
!     do moi=icore+1,iocc
!     do moaa=0,1
!     do moa=iocc+1,iall(0,0,0)-ivirtcore
!     gg1=0.0d0
!     if ((moia==1).and.(moaa==pab(1))) gg1=gg1+g(q,moi,moa,p1(1))
!     if ((moia==pab(1)).and.(moia==moaa)) gg1=gg1-g(q,p1(1),moa,moi)
!     gg2=0.0d0
!     if ((hab(1)==pab(2)).and.(moia==moaa)) gg2=gg2+g(h1(1),p1(2),moi,moa)
!     if ((hab(1)==moaa).and.(moia==pab(2))) gg2=gg2-g(h1(1),moa,moi,p1(2))
!     gg=gg+gg1*gg2/(epsilon(h1(1),0)+epsilon(moi,0)-epsilon(p1(2),0)-epsilon(moa,0))
!     gg1=0.0d0
!     if ((moia==1).and.(moaa==pab(2))) gg1=gg1+g(q,moi,moa,p1(2))
!     if ((moia==pab(2)).and.(moia==moaa)) gg1=gg1-g(q,p1(2),moa,moi)
!     gg2=0.0d0
!     if ((hab(1)==pab(1)).and.(moia==moaa)) gg2=gg2+g(h1(1),p1(1),moi,moa)
!     if ((hab(1)==moaa).and.(moia==pab(1))) gg2=gg2-g(h1(1),moa,moi,p1(1))
!     gg=gg-gg1*gg2/(epsilon(h1(1),0)+epsilon(moi,0)-epsilon(p1(1),0)-epsilon(moa,0))
!     enddo
!     enddo
!     enddo
!     enddo
!     gg1=0.0d0
!     if ((pab(1)==1).and.(hab(1)==pab(2))) gg1=gg1+g(q,p1(1),h1(1),p1(2))
!     if ((pab(2)==1).and.(hab(1)==pab(1))) gg1=gg1-g(q,p1(2),h1(1),p1(1))
!     sig=sig+gg*gg1/(omega-denom(pa,pb))
!     write(*,'(4i3,4f20.10)') q,h1(1),p1(1),p1(2),gg1,gg,V_IPEA(PA,PB,QA,QB,0,1),V_IPEA(PA,PB,QA,QB,1,0)
!    else if (hn==2) then
!     gg=0.0d0
!     do moaa=0,1
!     do moa=iocc+1,iall(0,0,0)-ivirtcore
!     do moba=0,1
!     do mob=iocc+1,iall(0,0,0)-ivirtcore
!     gg1=0.0d0
!     if ((pab(1)==moaa).and.(moba==1)) gg1=gg1+g(p1(1),moa,q,mob)
!     if ((pab(1)==moba).and.(moaa==1)) gg1=gg1-g(p1(1),mob,q,moa)
!     gg2=0.0d0
!     if ((hab(1)==moaa).and.(hab(2)==moba)) gg2=gg2+g(h1(1),moa,h1(2),mob)
!     if ((hab(1)==moba).and.(hab(2)==moaa)) gg2=gg2-g(h1(1),mob,h1(2),moa)
!     gg=gg+0.5d0*gg1*gg2/(epsilon(h1(1),0)+epsilon(h1(2),0)-epsilon(moa,0)-epsilon(mob,0))
!     enddo
!     enddo
!     enddo
!     enddo
!     do moia=0,1
!     do moi=icore+1,iocc
!     do moaa=0,1
!     do moa=iocc+1,iall(0,0,0)-ivirtcore
!     gg1=0.0d0
!     if ((moaa==1).and.(moia==hab(1))) gg1=gg1+g(q,moa,moi,h1(1))
!     if ((hab(1)==1).and.(moia==moaa)) gg1=gg1-g(q,h1(1),moi,moa)
!     gg2=0.0d0
!     if ((hab(2)==pab(1)).and.(moia==moaa)) gg2=gg2+g(h1(2),p1(1),moi,moa)
!     if ((hab(2)==moaa).and.(moia==pab(1))) gg2=gg2-g(h1(2),moa,moi,p1(1))
!     gg=gg+gg1*gg2/(epsilon(h1(2),0)+epsilon(moi,0)-epsilon(p1(1),0)-epsilon(moa,0))
!     gg1=0.0d0
!     if ((moaa==1).and.(moia==hab(2))) gg1=gg1+g(q,moa,moi,h1(2))
!     if ((hab(2)==1).and.(moia==moaa)) gg1=gg1-g(q,h1(2),moi,moa)
!     gg2=0.0d0
!     if ((hab(1)==pab(1)).and.(moia==moaa)) gg2=gg2+g(h1(1),p1(1),moi,moa)
!     if ((hab(1)==moaa).and.(moia==pab(1))) gg2=gg2-g(h1(1),moa,moi,p1(1))
!     gg=gg-gg1*gg2/(epsilon(h1(1),0)+epsilon(moi,0)-epsilon(p1(1),0)-epsilon(moa,0))
!     enddo
!     enddo
!     enddo
!     enddo
!     gg1=0.0d0
!     if ((hab(1)==1).and.(pab(1)==hab(2))) gg1=gg1+g(q,h1(1),p1(1),h1(2))
!     if ((hab(2)==1).and.(pab(1)==hab(1))) gg1=gg1-g(q,h1(2),p1(1),h1(1))
!     sig=sig+gg*gg1/(omega-denom(pa,pb))
!     write(*,'(4i3,4f20.10)') q,h1(1),h1(2),p1(1),gg1,gg,V_IPEA(PA,PB,QA,QB,0,1),V_IPEA(PA,PB,QA,QB,1,0)
!    else
!     call pabort('bug')
!    endif
!   ENDIF
!  ENDDO
! ENDDO
! write(*,'(A,F20.10)') 'sig(3)=',sig
!ENDDO
SIGMA=0.0D0
TIGMA=0.0D0
write(*,*) '============== 2nd order'
write(*,*) '************ <0(f1|V f35)(f35|G0 f35)(f35|V f1)|0> '
DO QA=1,NCF
 DO QB=1,IP_NCF+EA_NCF
  DO PA=1,NCF
   DO PB=1,IP_NCF+EA_NCF
    S_IPEA(PA,PB,QA,QB)=0.0D0
    DO RA=1,NCF
     DO RB=1,IP_NCF+EA_NCF
      S_IPEA(PA,PB,QA,QB)=S_IPEA(PA,PB,QA,QB)+ &
      V_IPEA(PA,PB,RA,RB,0,0)*V_IPEA(RA,RB,QA,QB,0,0)/(OMEGA-DENOM(RA,RB))
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
 DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
  TIGMA(P,Q)=S_IPEA(1,IPEA(P),1,IPEA(Q))
 ENDDO
ENDDO
SIGMA(:,:,2)=SIGMA(:,:,2)+TIGMA
CALL DUMP5(TIGMA(:,:),IALL(0,0,0)-IVIRTCORE)
write(*,*) '************ <0(f1|V f3)(f3|G0 f3)(f3|V f1)|0> '
DO QA=1,NCF
 DO QB=1,IP_NCF+EA_NCF
  DO PA=1,NCF
   DO PB=1,IP_NCF+EA_NCF
    S_IPEA(PA,PB,QA,QB)=0.0D0
    DO RA=1,NCF
     DO RB=1,IP_NCF+EA_NCF
      IF (RANK(RA,RB)==2) THEN
       S_IPEA(PA,PB,QA,QB)=S_IPEA(PA,PB,QA,QB)+ &
       V_IPEA(PA,PB,RA,RB,0,0)*V_IPEA(RA,RB,QA,QB,0,0)/(OMEGA-DENOM(RA,RB))
      ENDIF
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
 DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
  TIGMA(P,Q)=S_IPEA(1,IPEA(P),1,IPEA(Q))
 ENDDO
ENDDO
CALL DUMP5(TIGMA(:,:),IALL(0,0,0)-IVIRTCORE)
write(*,*) '************ <0|(f1|V f1)|1> '
DO QA=1,NCF
 DO QB=1,IP_NCF+EA_NCF
  DO PA=1,NCF
   DO PB=1,IP_NCF+EA_NCF
    S_IPEA(PA,PB,QA,QB)=V_IPEA(PA,PB,QA,QB,0,1)
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
 DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
  TIGMA(P,Q)=S_IPEA(1,IPEA(P),1,IPEA(Q))
 ENDDO
ENDDO
SIGMA(:,:,2)=SIGMA(:,:,2)+TIGMA
CALL DUMP5(TIGMA(:,:),IALL(0,0,0)-IVIRTCORE)
write(*,*) '************ <1|(f1|V f1)|0> '
DO QA=1,NCF
 DO QB=1,IP_NCF+EA_NCF
  DO PA=1,NCF
   DO PB=1,IP_NCF+EA_NCF
    S_IPEA(PA,PB,QA,QB)=V_IPEA(PA,PB,QA,QB,1,0)
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
 DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
  TIGMA(P,Q)=S_IPEA(1,IPEA(P),1,IPEA(Q))
 ENDDO
ENDDO
SIGMA(:,:,2)=SIGMA(:,:,2)+TIGMA
CALL DUMP5(TIGMA(:,:),IALL(0,0,0)-IVIRTCORE)
write(*,*) "=========== 2ND-ORDER SIGMA"
CALL DUMP5(SIGMA(:,:,2),IALL(0,0,0)-IVIRTCORE)
write(*,*) '============== 3rd order'
write(*,*) '************ (1,f3|V|f3,1) '
DO QA=1,NCF
 DO QB=1,IP_NCF+EA_NCF
  DO PA=1,NCF
   DO PB=1,IP_NCF+EA_NCF
    S_IPEA(PA,PB,QA,QB)=V_IPEA(PA,PB,QA,QB,1,1)
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
 DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
  TIGMA(P,Q)=S_IPEA(1,IPEA(P),1,IPEA(Q))
 ENDDO
ENDDO
CALL DUMP5(TIGMA(:,:),IALL(0,0,0)-IVIRTCORE)
SIGMA(:,:,3)=SIGMA(:,:,3)+TIGMA
write(*,*) '************ (2,f3|V|f3,0) '
DO QA=1,NCF
 DO QB=1,IP_NCF+EA_NCF
  DO PA=1,NCF
   DO PB=1,IP_NCF+EA_NCF
    S_IPEA(PA,PB,QA,QB)=V_IPEA(PA,PB,QA,QB,2,0)
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
 DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
  TIGMA(P,Q)=S_IPEA(1,IPEA(P),1,IPEA(Q))
 ENDDO
ENDDO
CALL DUMP5(TIGMA(:,:),IALL(0,0,0)-IVIRTCORE)
SIGMA(:,:,3)=SIGMA(:,:,3)+TIGMA
write(*,*) '************ (0,f3|V|f3,2) '
DO QA=1,NCF
 DO QB=1,IP_NCF+EA_NCF
  DO PA=1,NCF
   DO PB=1,IP_NCF+EA_NCF
    S_IPEA(PA,PB,QA,QB)=V_IPEA(PA,PB,QA,QB,0,2)
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
 DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
  TIGMA(P,Q)=S_IPEA(1,IPEA(P),1,IPEA(Q))
 ENDDO
ENDDO
CALL DUMP5(TIGMA(:,:),IALL(0,0,0)-IVIRTCORE)
SIGMA(:,:,3)=SIGMA(:,:,3)+TIGMA
write(*,*) '************ (0,f3|V G0 V G0 V|f3,0) '
DO QA=1,NCF
 DO QB=1,IP_NCF+EA_NCF
  DO PA=1,NCF
   DO PB=1,IP_NCF+EA_NCF
    S_IPEA(PA,PB,QA,QB)=0.0D0
    DO RA=1,NCF
     DO RB=1,IP_NCF+EA_NCF
      IF (RANK(RA,RB)==2) THEN
       S_IPEA(PA,PB,QA,QB)=S_IPEA(PA,PB,QA,QB)+ &
       V_IPEA(PA,PB,RA,RB,0,0)*V_IPEA(RA,RB,QA,QB,0,0)/(OMEGA-DENOM(RA,RB))
      ENDIF
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO QA=1,NCF
 DO QB=1,IP_NCF+EA_NCF
  DO PA=1,NCF
   DO PB=1,IP_NCF+EA_NCF
    T_IPEA(PA,PB,QA,QB)=0.0D0
    DO RA=1,NCF
     DO RB=1,IP_NCF+EA_NCF
      IF (RANK(RA,RB)==2) THEN
       T_IPEA(PA,PB,QA,QB)=T_IPEA(PA,PB,QA,QB)+ &
       V_IPEA(PA,PB,RA,RB,0,0)*S_IPEA(RA,RB,QA,QB)/(OMEGA-DENOM(RA,RB))
      ENDIF
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
 DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
  TIGMA(P,Q)=T_IPEA(1,IPEA(P),1,IPEA(Q))
 ENDDO
ENDDO
CALL DUMP5(TIGMA(:,:),IALL(0,0,0)-IVIRTCORE)
SIGMA(:,:,3)=SIGMA(:,:,3)+TIGMA
write(*,*) '************ (0,f3| V G0 V |f3,1) '
DO QA=1,NCF
 DO QB=1,IP_NCF+EA_NCF
  DO PA=1,NCF
   DO PB=1,IP_NCF+EA_NCF
    S_IPEA(PA,PB,QA,QB)=0.0D0
    DO RA=1,NCF
     DO RB=1,IP_NCF+EA_NCF
      IF (RANK(RA,RB)==2) THEN
       S_IPEA(PA,PB,QA,QB)=S_IPEA(PA,PB,QA,QB)+ &
       V_IPEA(PA,PB,RA,RB,0,0)*(V_IPEA(RA,RB,QA,QB,0,1)+V_IPEA(RA,RB,QA,QB,1,0)) &
!      V_IPEA(PA,PB,RA,RB,0,0)*V_IPEA(RA,RB,QA,QB,1,0) &
       /(OMEGA-DENOM(RA,RB))
      ENDIF
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
 DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
  TIGMA(P,Q)=S_IPEA(1,IPEA(P),1,IPEA(Q))
 ENDDO
ENDDO
CALL DUMP5(TIGMA(:,:),IALL(0,0,0)-IVIRTCORE)
SIGMA(:,:,3)=SIGMA(:,:,3)+TIGMA
SIGMA(:,:,3)=SIGMA(:,:,3)+TRANSPOSE(TIGMA)
write(*,*) "=========== 3RD-ORDER SIGMA (TRANSPOSE OF THE ABOVE ADDED)"
CALL DUMP5(SIGMA(:,:,3),IALL(0,0,0)-IVIRTCORE)
SIGMA(:,:,3)=SIGMA(:,:,3)-TRANSPOSE(TIGMA)
write(*,*) '************ (1,f3| V G0 V |f3,0) '
DO QA=1,NCF
 DO QB=1,IP_NCF+EA_NCF
  DO PA=1,NCF
   DO PB=1,IP_NCF+EA_NCF
    S_IPEA(PA,PB,QA,QB)=0.0D0
    DO RA=1,NCF
     DO RB=1,IP_NCF+EA_NCF
      IF (RANK(RA,RB)==2) THEN
       S_IPEA(PA,PB,QA,QB)=S_IPEA(PA,PB,QA,QB)+ &
       (V_IPEA(PA,PB,RA,RB,1,0)+V_IPEA(PA,PB,RA,RB,0,1))*V_IPEA(RA,RB,QA,QB,0,0) &
!      V_IPEA(PA,PB,RA,RB,0,1)*V_IPEA(RA,RB,QA,QB,0,0) &
       /(OMEGA-DENOM(RA,RB))
      ENDIF
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
 DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
  TIGMA(P,Q)=S_IPEA(1,IPEA(P),1,IPEA(Q))
 ENDDO
ENDDO
CALL DUMP5(TIGMA(:,:),IALL(0,0,0)-IVIRTCORE)
SIGMA(:,:,3)=SIGMA(:,:,3)+TIGMA
write(*,*) "=========== 3RD-ORDER SIGMA (ABOVE USED INSTEAD)"
CALL DUMP5(SIGMA(:,:,3),IALL(0,0,0)-IVIRTCORE)
stop
!write(*,*) '************ (0,f3|V G0 V|f3,1) '
!DO QA=1,NCF
! DO QB=1,IP_NCF+EA_NCF
!  DO PA=1,NCF
!   DO PB=1,IP_NCF+EA_NCF
!    S_IPEA(PA,PB,QA,QB)=0.0D0
!    DO RA=1,NCF
!     DO RB=1,IP_NCF+EA_NCF
!      IF (RANK(RA,RB)<=2) THEN
!       S_IPEA(PA,PB,QA,QB)=S_IPEA(PA,PB,QA,QB)+V_IPEA(PA,PB,RA,RB,0,0)*V_IPEA(RA,RB,QA,QB,0,1)/(OMEGA-DENOM(RA,RB))
!      ENDIF
!     ENDDO
!    ENDDO
!   ENDDO
!  ENDDO
! ENDDO
!ENDDO
!DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
! DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
!  TIGMA(P,Q)=S_IPEA(1,IPEA(P),1,IPEA(Q))
! ENDDO
!ENDDO
!CALL DUMP5(TIGMA(:,:),IALL(0,0,0)-IVIRTCORE)
!write(*,*) '************ (1,f3|V G0 V|f3,0) '
!DO QA=1,NCF
! DO QB=1,IP_NCF+EA_NCF
!  DO PA=1,NCF
!   DO PB=1,IP_NCF+EA_NCF
!    S_IPEA(PA,PB,QA,QB)=0.0D0
!    DO RA=1,NCF
!     DO RB=1,IP_NCF+EA_NCF
!      IF (RANK(RA,RB)<=2) THEN
!       S_IPEA(PA,PB,QA,QB)=S_IPEA(PA,PB,QA,QB)+V_IPEA(PA,PB,RA,RB,1,0)*V_IPEA(RA,RB,QA,QB,0,0)/(OMEGA-DENOM(RA,RB))
!      ENDIF
!     ENDDO
!    ENDDO
!   ENDDO
!  ENDDO
! ENDDO
!ENDDO
!DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
! DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
!  TIGMA(P,Q)=S_IPEA(1,IPEA(P),1,IPEA(Q))
! ENDDO
!ENDDO
!CALL DUMP5(TIGMA(:,:),IALL(0,0,0)-IVIRTCORE)
! end debug

   ! CONSTRUCT SIGMA
   DO IORDER=1,ORDER
    WRITE(6,'(/,A,I2,A)') 'CONSTRUCTING',IORDER,'TH SELF-ENERGY'
    DO BRAORDER=0,IORDER-1
     DO KETORDER=0,IORDER-1-BRAORDER
      TIGMA=0.0D0
      WRITE(6,'(A,I3)') 'NUMBER OF ACTIONS OF V =',IORDER-BRAORDER-KETORDER
      IF (IORDER-BRAORDER-KETORDER==1) THEN
       DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
        DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
         TIGMA(P,Q)=V_IPEA(1,IPEA(P),1,IPEA(Q),BRAORDER,KETORDER)
        ENDDO
       ENDDO
      ELSE
       S_IPEA(:,:,:,:)=V_IPEA(:,:,:,:,0,KETORDER)
       DO JORDER=1,IORDER-BRAORDER-KETORDER
        ! G0**(-1)
        DO QA=1,NCF
         DO QB=1,IP_NCF+EA_NCF
          DO PA=1,NCF
           DO PB=1,IP_NCF+EA_NCF
            IF (RANK(PA,PB)==1) THEN
             S_IPEA(PA,PB,QA,QB)=0.0D0
            ELSE IF (RANK(PA,PB)>=3) THEN
             S_IPEA(PA,PB,QA,QB)=0.0D0
            ELSE
             S_IPEA(PA,PB,QA,QB)=S_IPEA(PA,PB,QA,QB)/(OMEGA-DENOM(PA,PB))
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        ! V
        DO QA=1,NCF
         DO QB=1,IP_NCF+EA_NCF
          DO PA=1,NCF
           DO PB=1,IP_NCF+EA_NCF
            T_IPEA(PA,PB,QA,QB)=0.0D0
            IF (JORDER==IORDER-BRAORDER-KETORDER) THEN
             DO RA=1,NCF
              DO RB=1,IP_NCF+EA_NCF
               T_IPEA(PA,PB,QA,QB)=T_IPEA(PA,PB,QA,QB)+V_IPEA(PA,PB,RA,RB,BRAORDER,0)*S_IPEA(RA,RB,QA,QB)
              ENDDO
             ENDDO
            ELSE
             DO RA=1,NCF
              DO RB=1,IP_NCF+EA_NCF
               T_IPEA(PA,PB,QA,QB)=T_IPEA(PA,PB,QA,QB)+V_IPEA(PA,PB,RA,RB,0,0)*S_IPEA(RA,RB,QA,QB)
              ENDDO
             ENDDO
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        S_IPEA=T_IPEA
       ENDDO
       DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
        DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
         TIGMA(P,Q)=S_IPEA(1,IPEA(P),1,IPEA(Q))
        ENDDO
       ENDDO
      ENDIF
      WRITE(6,'(A,I2,A,I1,A,I2,A)') '(',BRAORDER,'|(VG)^',IORDER-BRAORDER-KETORDER-1,'V',KETORDER,')'
      CALL DUMP5(TIGMA(:,:),IALL(0,0,0)-IVIRTCORE)
      SIGMA(:,:,IORDER)=SIGMA(:,:,IORDER)+TIGMA(:,:)
     ENDDO
    ENDDO
    WRITE(6,'(I3,A)') IORDER,'TH-ORDER SELF-ENERGY'
    CALL DUMP5(SIGMA(:,:,IORDER),IALL(0,0,0)-IVIRTCORE)
   ENDDO

!  ! CONSTRUCT SIGMA
!  DO IORDER=1,ORDER
!   WRITE(6,'(/,A,I2,A)') 'CONSTRUCTING',IORDER,'TH SELF-ENERGY'
!   DO BRAORDER=0,IORDER-1
!    DO KETORDER=0,IORDER-1-BRAORDER
!     TIGMA=0.0D0
!     WRITE(6,'(A,I3)') 'NUMBER OF ACTIONS OF V =',IORDER-BRAORDER-KETORDER
!     IF (IORDER-BRAORDER-KETORDER==1) THEN
!      DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
!       DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
!        TIGMA(P,Q)=W_IPEA(1,IPEA(P),1,IPEA(Q),IORDER)
!       ENDDO
!      ENDDO
!     ELSE
!      S_IPEA(:,:,:,:)=W_IPEA(:,:,:,:,KETORDER+1)
!      DO JORDER=1,IORDER-BRAORDER-KETORDER-1
!       ! G0**(-1)
!       DO QA=1,NCF
!        DO QB=1,IP_NCF+EA_NCF
!         DO PA=1,NCF
!          DO PB=1,IP_NCF+EA_NCF
!           IF (RANK(PA,PB)==1) THEN
!            S_IPEA(PA,PB,QA,QB)=0.0D0
!           ELSE IF (RANK(PA,PB)>=IORDER) THEN
!            S_IPEA(PA,PB,QA,QB)=0.0D0
!           ELSE
!            S_IPEA(PA,PB,QA,QB)=-S_IPEA(PA,PB,QA,QB)/(OMEGA-DENOM(PA,PB))
!           ENDIF
!          ENDDO
!         ENDDO
!        ENDDO
!       ENDDO
!       ! V
!       DO QA=1,NCF
!        DO QB=1,IP_NCF+EA_NCF
!         DO PA=1,NCF
!          DO PB=1,IP_NCF+EA_NCF
!           T_IPEA(PA,PB,QA,QB)=0.0D0
!           IF (JORDER==IORDER-BRAORDER-KETORDER) THEN
!            DO RA=1,NCF
!             DO RB=1,IP_NCF+EA_NCF
!              T_IPEA(PA,PB,QA,QB)=T_IPEA(PA,PB,QA,QB)+W_IPEA(PA,PB,RA,RB,BRAORDER+1)*S_IPEA(RA,RB,QA,QB)
!             ENDDO
!            ENDDO
!           ELSE
!            DO RA=1,NCF
!             DO RB=1,IP_NCF+EA_NCF
!              T_IPEA(PA,PB,QA,QB)=T_IPEA(PA,PB,QA,QB)+V_IPEA(PA,PB,RA,RB,0,0)*S_IPEA(RA,RB,QA,QB)
!             ENDDO
!            ENDDO
!           ENDIF
!          ENDDO
!         ENDDO
!        ENDDO
!       ENDDO
!       S_IPEA=T_IPEA
!      ENDDO
!      DO P=ICORE+1,IALL(0,0,0)-IVIRTCORE
!       DO Q=ICORE+1,IALL(0,0,0)-IVIRTCORE
!        TIGMA(P,Q)=S_IPEA(1,IPEA(P),1,IPEA(Q))
!       ENDDO
!      ENDDO
!     ENDIF
!     IF (IORDER-BRAORDER-KETORDER==1) THEN
!      WRITE(6,'(A,I1,A,I1,A)') '(',BRAORDER,'|V',KETORDER,')'
!     ELSE
!      WRITE(6,'(A,I1,A,I1,A,I1,A)') '(',BRAORDER,'|(VG)^',IORDER-BRAORDER-KETORDER-1,'V',KETORDER,')'
!     ENDIF
!     CALL DUMP5(TIGMA(:,:),IALL(0,0,0)-IVIRTCORE)
!     SIGMA(:,:,IORDER)=SIGMA(:,:,IORDER)+TIGMA(:,:)
!    ENDDO
!   ENDDO
!   WRITE(6,'(I3,A)') IORDER,'TH-ORDER SELF-ENERGY'
!   CALL DUMP5(SIGMA(:,:,IORDER),IALL(0,0,0)-IVIRTCORE)
!  ENDDO

   DEALLOCATE(VEC0,VEC1,BRA,KET,T_IPEA,S_IPEA,V_IPEA,W_IPEA,SIGMA,TIGMA,DENOM,RANK,IPEA)
   CLOSE(50)
   CLOSE(51)
   CLOSE(52)

!  CALL PCPU_TIME(ICPUE)
   CALL PFLUSH(6)

END SUBROUTINE



SUBROUTINE IONIZE(INFILE,OUTFILE,PA,PB,DAGGER,STATE)
! DAGGER=.FALSE.: IONIZE
! DAGGER=.TRUE.:  ELECTRON ATTACH
! STATE=0:  NEUTRAL REFERENCE
! STATE=+1: CATION REFERENCE
! STATE=-1: ANION REFERENCE

   USE CONTROL
   USE STRUCTURE
   USE INTEGRAL
   USE BASISSET
   USE FULLCI

   IMPLICIT NONE
   INTEGER :: INFILE,OUTFILE,STATE
   INTEGER :: PA,PB
   LOGICAL :: DAGGER,EA
   DOUBLE PRECISION,ALLOCATABLE :: TRL(:,:),PRD(:,:)
   INTEGER :: I,J,IC,ID
   INTEGER(4) :: CFONE,CFTWO
   INTEGER :: NHA,NPA,NHB,NPB,JSIGN,KSIGN
   INTEGER :: HOLEA(IALL(0,0,0)),HOLEB(IALL(0,0,0)),PARTA(IALL(0,0,0)),PARTB(IALL(0,0,0))
   LOGICAL :: LCYCLE

   IF (PB <= IP_NCF) THEN
    EA=.FALSE.
   ELSE
    EA=.TRUE.
!   DAGGER=.NOT.DAGGER
   ENDIF

   IF (STATE==0) THEN
    ALLOCATE(TRL(NCF,NCF))
    IF (.NOT.DAGGER) THEN
     ALLOCATE(PRD(IP_NCF,NCF))
    ELSE
     ALLOCATE(PRD(EA_NCF,NCF))
    ENDIF
   ELSE IF (STATE==1) THEN
    ALLOCATE(TRL(IP_NCF,NCF))
    IF (.NOT.DAGGER) THEN
     CALL PABORT('REACHED N-2 STATES')
    ELSE
     ALLOCATE(PRD(NCF,NCF))
    ENDIF
   ELSE IF (STATE==-1) THEN
    ALLOCATE(TRL(EA_NCF,NCF))
    IF (.NOT.DAGGER) THEN
     ALLOCATE(PRD(NCF,NCF))
    ELSE
     CALL PABORT('REACHED N+2 STATES')
    ENDIF
   ELSE
    CALL PABORT('ILLEGAL STATE ARGUMENT')
   ENDIF

   REWIND(INFILE)
   READ(INFILE) TRL
   PRD=0.0D0

   IF (EA) THEN

    IF (DAGGER) THEN
 
     ! =================
     ! ELECTRON ATTACHED 
     ! =================
 
     NHA=0
     NPA=0
     DO I=IOCC,1,-1 ! DAGGER
      IF (.NOT.BTEST(CFHALF(PA),I-1)) THEN
       NHA=NHA+1
       HOLEA(NHA)=I
      ENDIF
     ENDDO
     DO I=IALL(0,0,0)-IVIRTCORE,IOCC+1,-1 ! DAGGER
      IF (BTEST(CFHALF(PA),I-1)) THEN
       NPA=NPA+1
       PARTA(NPA)=I
      ENDIF
     ENDDO
 
     NHB=0
     NPB=0
     DO I=IOCC,1,-1 ! DAGGER
      IF (.NOT.BTEST(EA_CFHALF(PB-IP_NCF),I-1)) THEN
       NHB=NHB+1
       HOLEB(NHB)=I
      ENDIF
     ENDDO
     DO I=IALL(0,0,0)-IVIRTCORE,IOCC+1,-1 ! DAGGER
      IF (BTEST(EA_CFHALF(PB-IP_NCF),I-1)) THEN
       NPB=NPB+1
       PARTB(NPB)=I
      ENDIF
     ENDDO
 
    ELSE
 
     ! =======
     ! IONIZED
     ! =======
 
     NHA=0
     NPA=0
     DO I=1,IOCC
      IF (.NOT.BTEST(CFHALF(PA),I-1)) THEN
       NPA=NPA+1
       PARTA(NPA)=I
      ENDIF
     ENDDO
     DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
      IF (BTEST(CFHALF(PA),I-1)) THEN
       NHA=NHA+1
       HOLEA(NHA)=I
      ENDIF
     ENDDO
 
     NHB=0
     NPB=0
     DO I=1,IOCC
      IF (.NOT.BTEST(EA_CFHALF(PB-IP_NCF),I-1)) THEN
       NPB=NPB+1
       PARTB(NPB)=I
      ENDIF
     ENDDO
     DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
      IF (BTEST(EA_CFHALF(PB-IP_NCF),I-1)) THEN
       NHB=NHB+1
       HOLEB(NHB)=I
      ENDIF
     ENDDO
 
    ENDIF

   ELSE

    IF (.NOT.DAGGER) THEN
 
     ! =======
     ! IONIZED
     ! =======
 
     NHA=0
     NPA=0
     DO I=1,IOCC
      IF (.NOT.BTEST(CFHALF(PA),I-1)) THEN
       NHA=NHA+1
       HOLEA(NHA)=I
      ENDIF
     ENDDO
     DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
      IF (BTEST(CFHALF(PA),I-1)) THEN
       NPA=NPA+1
       PARTA(NPA)=I
      ENDIF
     ENDDO
 
     NHB=0
     NPB=0
     DO I=1,IOCC
      IF (.NOT.BTEST(IP_CFHALF(PB),I-1)) THEN
       NHB=NHB+1
       HOLEB(NHB)=I
      ENDIF
     ENDDO
     DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
      IF (BTEST(IP_CFHALF(PB),I-1)) THEN
       NPB=NPB+1
       PARTB(NPB)=I
      ENDIF
     ENDDO
 
    ELSE
 
     ! =================
     ! ELECTRON ATTACHED
     ! =================
 
     NHA=0
     NPA=0
     DO I=IOCC,1,-1 ! DAGGER
      IF (.NOT.BTEST(CFHALF(PA),I-1)) THEN
       NPA=NPA+1
       PARTA(NPA)=I
      ENDIF
     ENDDO
     DO I=IALL(0,0,0)-IVIRTCORE,IOCC+1,-1 ! DAGGER
      IF (BTEST(CFHALF(PA),I-1)) THEN
       NHA=NHA+1
       HOLEA(NHA)=I
      ENDIF
     ENDDO
 
     NHB=0
     NPB=0
     DO I=IOCC,1,-1 ! DAGGER
      IF (.NOT.BTEST(IP_CFHALF(PB),I-1)) THEN
       NPB=NPB+1
       PARTB(NPB)=I
      ENDIF
     ENDDO
     DO I=IALL(0,0,0)-IVIRTCORE,IOCC+1,-1 ! DAGGER
      IF (BTEST(IP_CFHALF(PB),I-1)) THEN
       NHB=NHB+1
       HOLEB(NHB)=I
      ENDIF
     ENDDO
 
    ENDIF

   ENDIF

   IF (STATE==0) THEN

    ! -----------------
    ! NEUTRAL REFERENCE
    ! -----------------

    DO IC=1,NCF
     JSIGN=1
     CFONE=CFHALF(IC)
     LCYCLE=.FALSE.
!    IF (.NOT.DAGGER) THEN
!     IF (NPA >= 1) THEN
!      DO I=1,NPA
!       IF (.NOT.BTEST(CFONE,PARTA(I)-1)) THEN
!        IF (PARTA(I) > 1) THEN
!         DO J=1,PARTA(I)-1
!          IF (BTEST(CFONE,J-1)) JSIGN=-JSIGN
!         ENDDO
!        ENDIF
!        CFONE=IBSET(CFONE,PARTA(I)-1)
!       ELSE
!        LCYCLE=.TRUE.
!       ENDIF
!      ENDDO
!      IF (LCYCLE) CYCLE
!     ENDIF
!     IF (NHA >= 1) THEN
!      DO I=1,NHA
!       IF (BTEST(CFONE,HOLEA(I)-1)) THEN
!        CFONE=IBCLR(CFONE,HOLEA(I)-1)
!        IF (HOLEA(I) > 1) THEN
!         DO J=1,HOLEA(I)-1
!          IF (BTEST(CFONE,J-1)) JSIGN=-JSIGN
!         ENDDO
!        ENDIF
!       ELSE
!        LCYCLE=.TRUE.
!       ENDIF
!      ENDDO
!      IF (LCYCLE) CYCLE
!     ENDIF
!    ELSE
      IF (NHA >= 1) THEN
       DO I=1,NHA
        IF (BTEST(CFONE,HOLEA(I)-1)) THEN
         CFONE=IBCLR(CFONE,HOLEA(I)-1)
         IF (HOLEA(I) > 1) THEN
          DO J=1,HOLEA(I)-1
           IF (BTEST(CFONE,J-1)) JSIGN=-JSIGN
          ENDDO
         ENDIF
        ELSE
         LCYCLE=.TRUE.
        ENDIF
       ENDDO
       IF (LCYCLE) CYCLE
      ENDIF
      IF (NPA >= 1) THEN
       DO I=1,NPA
        IF (.NOT.BTEST(CFONE,PARTA(I)-1)) THEN
         IF (PARTA(I) > 1) THEN
          DO J=1,PARTA(I)-1
           IF (BTEST(CFONE,J-1)) JSIGN=-JSIGN
          ENDDO
         ENDIF
         CFONE=IBSET(CFONE,PARTA(I)-1)
        ELSE
         LCYCLE=.TRUE.
        ENDIF
       ENDDO
       IF (LCYCLE) CYCLE
      ENDIF
!    ENDIF

     DO ID=1,NCF
      KSIGN=1
      CFTWO=CFHALF(ID)
      LCYCLE=.FALSE.
!     IF (.NOT.DAGGER) THEN
!      IF (NPB >= 1) THEN
!       DO I=1,NPB
!        IF (.NOT.BTEST(CFTWO,PARTB(I)-1)) THEN
!         IF (PARTB(I) > 1) THEN
!          DO J=1,PARTB(I)-1
!           IF (BTEST(CFTWO,J-1)) KSIGN=-KSIGN
!          ENDDO
!         ENDIF
!         CFTWO=IBSET(CFTWO,PARTB(I)-1)
!        ELSE
!         LCYCLE=.TRUE.
!        ENDIF
!       ENDDO
!       IF (LCYCLE) CYCLE
!      ENDIF
!      IF (NHB >= 1) THEN
!       DO I=1,NHB
!        IF (BTEST(CFTWO,HOLEB(I)-1)) THEN
!         CFTWO=IBCLR(CFTWO,HOLEB(I)-1)
!         IF (HOLEB(I) > 1) THEN
!          DO J=1,HOLEB(I)-1
!           IF (BTEST(CFTWO,J-1)) KSIGN=-KSIGN
!          ENDDO
!         ENDIF
!        ELSE
!         LCYCLE=.TRUE.
!        ENDIF
!       ENDDO
!       IF (LCYCLE) CYCLE
!      ENDIF
!     ELSE
       IF (NHB >= 1) THEN
        DO I=1,NHB
         IF (BTEST(CFTWO,HOLEB(I)-1)) THEN
          CFTWO=IBCLR(CFTWO,HOLEB(I)-1)
          IF (HOLEB(I) > 1) THEN
           DO J=1,HOLEB(I)-1
            IF (BTEST(CFTWO,J-1)) KSIGN=-KSIGN
           ENDDO
          ENDIF
         ELSE
          LCYCLE=.TRUE.
         ENDIF
        ENDDO
        IF (LCYCLE) CYCLE
       ENDIF
       IF (NPB >= 1) THEN
        DO I=1,NPB
         IF (.NOT.BTEST(CFTWO,PARTB(I)-1)) THEN
          IF (PARTB(I) > 1) THEN
           DO J=1,PARTB(I)-1
            IF (BTEST(CFTWO,J-1)) KSIGN=-KSIGN
           ENDDO
          ENDIF
          CFTWO=IBSET(CFTWO,PARTB(I)-1)
         ELSE
          LCYCLE=.TRUE.
         ENDIF
        ENDDO
        IF (LCYCLE) CYCLE
       ENDIF
!     ENDIF
 
      IF (.NOT.DAGGER) THEN
       PRD(IP_ADDRSS(CFTWO),ADDRSS(CFONE))=PRD(IP_ADDRSS(CFTWO),ADDRSS(CFONE))+JSIGN*KSIGN*TRL(ID,IC)
      ELSE
       PRD(EA_ADDRSS(CFTWO),ADDRSS(CFONE))=PRD(EA_ADDRSS(CFTWO),ADDRSS(CFONE))+JSIGN*KSIGN*TRL(ID,IC)
      ENDIF

     ENDDO
    ENDDO

   ELSE IF (STATE==1) THEN

    ! ----------------
    ! CATION REFERENCE
    ! ----------------

    DO IC=1,NCF
     JSIGN=1
     CFONE=CFHALF(IC)
     LCYCLE=.FALSE.
!    IF (.NOT.DAGGER) THEN
!     IF (NPA >= 1) THEN
!      DO I=1,NPA
!       IF (.NOT.BTEST(CFONE,PARTA(I)-1)) THEN
!        IF (PARTA(I) > 1) THEN
!         DO J=1,PARTA(I)-1
!          IF (BTEST(CFONE,J-1)) JSIGN=-JSIGN
!         ENDDO
!        ENDIF
!        CFONE=IBSET(CFONE,PARTA(I)-1)
!       ELSE
!        LCYCLE=.TRUE.
!       ENDIF
!      ENDDO
!      IF (LCYCLE) CYCLE
!     ENDIF
!     IF (NHA >= 1) THEN
!      DO I=1,NHA
!       IF (BTEST(CFONE,HOLEA(I)-1)) THEN
!        CFONE=IBCLR(CFONE,HOLEA(I)-1)
!        IF (HOLEA(I) > 1) THEN
!         DO J=1,HOLEA(I)-1
!          IF (BTEST(CFONE,J-1)) JSIGN=-JSIGN
!         ENDDO
!        ENDIF
!       ELSE
!        LCYCLE=.TRUE.
!       ENDIF
!      ENDDO
!      IF (LCYCLE) CYCLE
!     ENDIF
!    ELSE
      IF (NHA >= 1) THEN
       DO I=1,NHA
        IF (BTEST(CFONE,HOLEA(I)-1)) THEN
         CFONE=IBCLR(CFONE,HOLEA(I)-1)
         IF (HOLEA(I) > 1) THEN
          DO J=1,HOLEA(I)-1
           IF (BTEST(CFONE,J-1)) JSIGN=-JSIGN
          ENDDO
         ENDIF
        ELSE
         LCYCLE=.TRUE.
        ENDIF
       ENDDO
       IF (LCYCLE) CYCLE
      ENDIF
      IF (NPA >= 1) THEN
       DO I=1,NPA
        IF (.NOT.BTEST(CFONE,PARTA(I)-1)) THEN
         IF (PARTA(I) > 1) THEN
          DO J=1,PARTA(I)-1
           IF (BTEST(CFONE,J-1)) JSIGN=-JSIGN
          ENDDO
         ENDIF
         CFONE=IBSET(CFONE,PARTA(I)-1)
        ELSE
         LCYCLE=.TRUE.
        ENDIF
       ENDDO
       IF (LCYCLE) CYCLE
      ENDIF
!    ENDIF
 
     DO ID=1,IP_NCF
      KSIGN=1
      CFTWO=IP_CFHALF(ID)
      LCYCLE=.FALSE.
!     IF (.NOT.DAGGER) THEN
!      IF (NPB >= 1) THEN
!       DO I=1,NPB
!        IF (.NOT.BTEST(CFTWO,PARTB(I)-1)) THEN
!         IF (PARTB(I) > 1) THEN
!          DO J=1,PARTB(I)-1
!           IF (BTEST(CFTWO,J-1)) KSIGN=-KSIGN
!          ENDDO
!         ENDIF
!         CFTWO=IBSET(CFTWO,PARTB(I)-1)
!        ELSE
!         LCYCLE=.TRUE.
!        ENDIF
!       ENDDO
!       IF (LCYCLE) CYCLE
!      ENDIF
!      IF (NHB >= 1) THEN
!       DO I=1,NHB
!        IF (BTEST(CFTWO,HOLEB(I)-1)) THEN
!         CFTWO=IBCLR(CFTWO,HOLEB(I)-1)
!         IF (HOLEB(I) > 1) THEN
!          DO J=1,HOLEB(I)-1
!           IF (BTEST(CFTWO,J-1)) KSIGN=-KSIGN
!          ENDDO
!         ENDIF
!        ELSE
!         LCYCLE=.TRUE.
!        ENDIF
!       ENDDO
!       IF (LCYCLE) CYCLE
!      ENDIF
!     ELSE
       IF (NHB >= 1) THEN
        DO I=1,NHB
         IF (BTEST(CFTWO,HOLEB(I)-1)) THEN
          CFTWO=IBCLR(CFTWO,HOLEB(I)-1)
          IF (HOLEB(I) > 1) THEN
           DO J=1,HOLEB(I)-1
            IF (BTEST(CFTWO,J-1)) KSIGN=-KSIGN
           ENDDO
          ENDIF
         ELSE
          LCYCLE=.TRUE.
         ENDIF
        ENDDO
        IF (LCYCLE) CYCLE
       ENDIF
       IF (NPB >= 1) THEN
        DO I=1,NPB
         IF (.NOT.BTEST(CFTWO,PARTB(I)-1)) THEN
          IF (PARTB(I) > 1) THEN
           DO J=1,PARTB(I)-1
            IF (BTEST(CFTWO,J-1)) KSIGN=-KSIGN
           ENDDO
          ENDIF
          CFTWO=IBSET(CFTWO,PARTB(I)-1)
         ELSE
          LCYCLE=.TRUE.
         ENDIF
        ENDDO
        IF (LCYCLE) CYCLE
       ENDIF
!     ENDIF

      PRD(ADDRSS(CFTWO),ADDRSS(CFONE))=PRD(ADDRSS(CFTWO),ADDRSS(CFONE))+JSIGN*KSIGN*TRL(ID,IC)

     ENDDO
    ENDDO

   ELSE IF (STATE==-1) THEN

    ! ---------------
    ! ANION REFERENCE
    ! ---------------

    DO IC=1,NCF
     JSIGN=1
     CFONE=CFHALF(IC)
     LCYCLE=.FALSE.
!    IF (.NOT.DAGGER) THEN
!     IF (NPA >= 1) THEN
!      DO I=1,NPA
!       IF (.NOT.BTEST(CFONE,PARTA(I)-1)) THEN
!        IF (PARTA(I) > 1) THEN
!         DO J=1,PARTA(I)-1
!          IF (BTEST(CFONE,J-1)) JSIGN=-JSIGN
!         ENDDO
!        ENDIF
!        CFONE=IBSET(CFONE,PARTA(I)-1)
!       ELSE
!        LCYCLE=.TRUE.
!       ENDIF
!      ENDDO
!      IF (LCYCLE) CYCLE
!     ENDIF
!     IF (NHA >= 1) THEN
!      DO I=1,NHA
!       IF (BTEST(CFONE,HOLEA(I)-1)) THEN
!        CFONE=IBCLR(CFONE,HOLEA(I)-1)
!        IF (HOLEA(I) > 1) THEN
!         DO J=1,HOLEA(I)-1
!          IF (BTEST(CFONE,J-1)) JSIGN=-JSIGN
!         ENDDO
!        ENDIF
!       ELSE
!        LCYCLE=.TRUE.
!       ENDIF
!      ENDDO
!      IF (LCYCLE) CYCLE
!     ENDIF
!    ELSE
      IF (NHA >= 1) THEN
       DO I=1,NHA
        IF (BTEST(CFONE,HOLEA(I)-1)) THEN
         CFONE=IBCLR(CFONE,HOLEA(I)-1)
         IF (HOLEA(I) > 1) THEN
          DO J=1,HOLEA(I)-1
           IF (BTEST(CFONE,J-1)) JSIGN=-JSIGN
          ENDDO
         ENDIF
        ELSE
         LCYCLE=.TRUE.
        ENDIF
       ENDDO
       IF (LCYCLE) CYCLE
      ENDIF
      IF (NPA >= 1) THEN
       DO I=1,NPA
        IF (.NOT.BTEST(CFONE,PARTA(I)-1)) THEN
         IF (PARTA(I) > 1) THEN
          DO J=1,PARTA(I)-1
           IF (BTEST(CFONE,J-1)) JSIGN=-JSIGN
          ENDDO
         ENDIF
         CFONE=IBSET(CFONE,PARTA(I)-1)
        ELSE
         LCYCLE=.TRUE.
        ENDIF
       ENDDO
       IF (LCYCLE) CYCLE
      ENDIF
!    ENDIF

     DO ID=1,EA_NCF
      KSIGN=1
      CFTWO=EA_CFHALF(ID)
      LCYCLE=.FALSE.
!     IF (.NOT.DAGGER) THEN
!      IF (NPB >= 1) THEN
!       DO I=1,NPB
!        IF (.NOT.BTEST(CFTWO,PARTB(I)-1)) THEN
!         IF (PARTB(I) > 1) THEN
!          DO J=1,PARTB(I)-1
!           IF (BTEST(CFTWO,J-1)) KSIGN=-KSIGN
!          ENDDO
!         ENDIF
!         CFTWO=IBSET(CFTWO,PARTB(I)-1)
!        ELSE
!         LCYCLE=.TRUE.
!        ENDIF
!       ENDDO
!       IF (LCYCLE) CYCLE
!      ENDIF
!      IF (NHB >= 1) THEN
!       DO I=1,NHB
!        IF (BTEST(CFTWO,HOLEB(I)-1)) THEN
!         CFTWO=IBCLR(CFTWO,HOLEB(I)-1)
!         IF (HOLEB(I) > 1) THEN
!          DO J=1,HOLEB(I)-1
!           IF (BTEST(CFTWO,J-1)) KSIGN=-KSIGN
!          ENDDO
!         ENDIF
!        ELSE
!         LCYCLE=.TRUE.
!        ENDIF
!       ENDDO
!       IF (LCYCLE) CYCLE
!      ENDIF
!     ELSE
       IF (NHB >= 1) THEN
        DO I=1,NHB
         IF (BTEST(CFTWO,HOLEB(I)-1)) THEN
          CFTWO=IBCLR(CFTWO,HOLEB(I)-1)
          IF (HOLEB(I) > 1) THEN
           DO J=1,HOLEB(I)-1
            IF (BTEST(CFTWO,J-1)) KSIGN=-KSIGN
           ENDDO
          ENDIF
         ELSE
          LCYCLE=.TRUE.
         ENDIF
        ENDDO
        IF (LCYCLE) CYCLE
       ENDIF
       IF (NPB >= 1) THEN
        DO I=1,NPB
         IF (.NOT.BTEST(CFTWO,PARTB(I)-1)) THEN
          IF (PARTB(I) > 1) THEN
           DO J=1,PARTB(I)-1
            IF (BTEST(CFTWO,J-1)) KSIGN=-KSIGN
           ENDDO
          ENDIF
          CFTWO=IBSET(CFTWO,PARTB(I)-1)
         ELSE
          LCYCLE=.TRUE.
         ENDIF
        ENDDO
        IF (LCYCLE) CYCLE
       ENDIF
!     ENDIF
 
      PRD(ADDRSS(CFTWO),ADDRSS(CFONE))=PRD(ADDRSS(CFTWO),ADDRSS(CFONE))+JSIGN*KSIGN*TRL(ID,IC)

     ENDDO
    ENDDO

   ENDIF

   ! STORE PRODUCT VECTOR
   REWIND(OUTFILE)
   WRITE(OUTFILE) PRD

   DEALLOCATE(TRL,PRD)

   RETURN
END SUBROUTINE



SUBROUTINE H0_PRODUCT(REFFILE,INFILE,STATE)
! ACT H0 ON A REFERENCE WAVE FUNCTION IN REFFILE AND WRITE THE RESULT IN INFILE
! STATE=0:  NEUTRAL REFERENCE
! STATE=+1: CATION REFERENCE
! STATE=-1: ANION REFERENCE

   USE CONTROL
   USE STRUCTURE
   USE INTEGRAL
   USE BASISSET
   USE FULLCI

   IMPLICIT NONE
   INTEGER :: REFFILE,INFILE,STATE
   DOUBLE PRECISION,ALLOCATABLE :: TRL(:,:),PRD(:,:)
   DOUBLE PRECISION :: HA,HB
   INTEGER :: IA,IB,I

   IF (STATE==0) THEN
    ! NEUTRAL
    ALLOCATE(TRL(NCF,NCF),PRD(NCF,NCF))
   ELSE IF (STATE==1) THEN
    ! IONIZED
    ALLOCATE(TRL(IP_NCF,NCF),PRD(IP_NCF,NCF))
   ELSE IF (STATE==-1) THEN
    ! ELECTRON ATTACHED
    ALLOCATE(TRL(EA_NCF,NCF),PRD(EA_NCF,NCF))
   ELSE
    CALL PABORT('ILLEGAL STATE ARGUMENT')
   ENDIF

   ! READ TRIAL VECTOR FROM FILE INFILE
   REWIND(REFFILE)
   READ(REFFILE) TRL

   PRD=0.0D0

   IF (STATE==0) THEN
    DO IB=1,NCF
     HB=0.0D0
     DO I=1,IOCC
      IF (BTEST(CFHALF(IB),I-1)) HB=HB+EPSILON(I,0,0,0)
     ENDDO
     DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
      IF (BTEST(CFHALF(IB),I-1)) HB=HB+EPSILON(I,0,0,0)
     ENDDO
     DO IA=1,NCF
      HA=0.0D0
      DO I=1,IOCC
       IF (BTEST(CFHALF(IA),I-1)) HA=HA+EPSILON(I,0,0,0)
      ENDDO
      DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
       IF (BTEST(CFHALF(IA),I-1)) HA=HA+EPSILON(I,0,0,0)
      ENDDO
      PRD(IB,IA)=PRD(IB,IA)+(HA+HB)*TRL(IB,IA)
     ENDDO
    ENDDO
   ELSE IF (STATE==1) THEN
    DO IB=1,IP_NCF
     HB=0.0D0
     DO I=1,IOCC
      IF (BTEST(IP_CFHALF(IB),I-1)) HB=HB+EPSILON(I,0,0,0)
     ENDDO
     DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
      IF (BTEST(IP_CFHALF(IB),I-1)) HB=HB+EPSILON(I,0,0,0)
     ENDDO
     DO IA=1,NCF
      HA=0.0D0
      DO I=1,IOCC
       IF (BTEST(CFHALF(IA),I-1)) HA=HA+EPSILON(I,0,0,0)
      ENDDO
      DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
       IF (BTEST(CFHALF(IA),I-1)) HA=HA+EPSILON(I,0,0,0)
      ENDDO
      PRD(IB,IA)=PRD(IB,IA)+(HA+HB)*TRL(IB,IA)
     ENDDO
    ENDDO
   ELSE IF (STATE==-1) THEN
    DO IB=1,EA_NCF
     HB=0.0D0
     DO I=1,IOCC
      IF (BTEST(EA_CFHALF(IB),I-1)) HB=HB+EPSILON(I,0,0,0)
     ENDDO
     DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
      IF (BTEST(EA_CFHALF(IB),I-1)) HB=HB+EPSILON(I,0,0,0)
     ENDDO
     DO IA=1,NCF
      HA=0.0D0
      DO I=1,IOCC
       IF (BTEST(CFHALF(IA),I-1)) HA=HA+EPSILON(I,0,0,0)
      ENDDO
      DO I=IOCC+1,IALL(0,0,0)-IVIRTCORE
       IF (BTEST(CFHALF(IA),I-1)) HA=HA+EPSILON(I,0,0,0)
      ENDDO
      PRD(IB,IA)=PRD(IB,IA)+(HA+HB)*TRL(IB,IA)
     ENDDO
    ENDDO
   ENDIF

   ! STORE PRODUCT VECTOR
   REWIND(INFILE)
   WRITE(INFILE) PRD

   DEALLOCATE(TRL,PRD)

   RETURN
END SUBROUTINE



SUBROUTINE HIGHORDER_MP_REDUX(MOP,MOQ,OMEGA,ORDER)
! PERFORM HIGH-ORDER MOELLER-PLESSET PERTURBATION CALCULATIONS IN A RECURSIVE ALGORITHM.

   USE CONTROL
   USE GRADIENT
   USE STRUCTURE
   USE INTEGRAL
   USE BASISSET
   USE FULLCI

   IMPLICIT NONE
   INTEGER,PARAMETER :: MAXFILE = 100
   INTEGER :: MOP,MOQ
   DOUBLE PRECISION :: OMEGA
   INTEGER :: ORDER
   INTEGER :: MOI,MOJ
   INTEGER :: IA,IB
   INTEGER :: I,J,K,IFILE
   REAL :: MEM,ICPUS,ICPUE
   DOUBLE PRECISION :: HA,HB
   DOUBLE PRECISION :: EMP,DMP(0:ORDER)
   DOUBLE PRECISION,ALLOCATABLE :: VEC1(:,:),VEC2(:,:)

!  CALL PCPU_TIME(ICPUS)
   WRITE(6,'(A,2I3)') 'RECURSIVE HIGH-ORDER MOELLER-PLESSET PERTURBATION CALCULATIONS WILL BE PERFORMED',MOP,MOQ
   WRITE(6,'(A,I2,A,F20.15,A)') 'ORBITAL ENERGY OF',MOP,' REPLACED BY OMEGA=',OMEGA,' HARTREE'
   MP_STORED=0.0D0
   IF (MOP/=MOQ) THEN
    WRITE(6,'(A)') '*********************************************'
    WRITE(6,'(A)') '* TOTAL & CORRELATION ENERGIES ARE ALL NULL *'
    WRITE(6,'(A)') '*********************************************'
    RETURN
   ENDIF
   MEM=16.0*2.0*NCF**2
   IF (MEM > 1000000.0) THEN
    WRITE(6,'(A,F7.1,A)') 'ESTIMATED MEMORY USAGE WILL BE ',MEM/1000000.0,' MB'
   ELSE IF (MEM > 1000.0) THEN
    WRITE(6,'(A,F7.1,A)') 'ESTIMATED MEMORY USAGE WILL BE ',MEM/1000.0,' KB'
   ELSE
    WRITE(6,'(A,F7.1,A)') 'ESTIMATED MEMORY USAGE WILL BE ',MEM,' B'
   ENDIF
   IF (MEM > DOPTN(28)*1000000.0) CALL PABORT('OUT OF MEMORY')
   WRITE(6,'(A)') '------------------------------------------------------'
   WRITE(6,'(A)') 'ORDER      CORRELATION        TOTAL ENERGY   CPU / SEC'
   ALLOCATE(VEC1(NCF,NCF),VEC2(NCF,NCF))
   OPEN(50,FILE=TRIM(COPTN(1))//'.fi0',FORM='UNFORMATTED')
   OPEN(51,FILE=TRIM(COPTN(1))//'.fo0',FORM='UNFORMATTED')

   ! INITIALIZE WAVEFUNCTION TO THE HARTREE-FOCK DETERMINANT
   EMP=NUCLEAR_REPULSION
   VEC1=0.0D0
   VEC1(1,1)=1.0D0
   REWIND(50)
   WRITE(50) VEC1

   ! DMP(0)
   DMP(0)=0.0D0
   DO MOI=1,IOCC
!!! For f5, comment in the next 6 comments
    IF (MOI==MOP) THEN
     DMP(0)=DMP(0)+EPSILON(MOI,0,0,0)+OMEGA
    ELSE IF (MOI==MOQ) THEN
     DMP(0)=DMP(0)+EPSILON(MOI,0,0,0)+OMEGA
    ELSE
     DMP(0)=DMP(0)+2.0D0*EPSILON(MOI,0,0,0)
    ENDIF
   ENDDO
!  The following line is never used in the code correct for GF2
!  DMP(0)=DMP(0)-OMEGA
   EMP=EMP+DMP(0)
   WRITE(6,'(I2,F20.10,F20.10,F12.1)') 0,DMP(0),EMP,0.0
   MP_STORED(0)=DMP(0)

   ! DMP(1) ONWARD
   DO IFILE=0,ORDER-1
    IF (IOPTN(9) >= 3) THEN
     REWIND(50)
     DO I=0,IFILE
      READ(50) VEC1
     ENDDO
     WRITE(6,'(I3,A)') IFILE,'TH ORDER WAVE FUNCTION'
     CALL DUMP5(VEC1,NCF)
    ENDIF
    DEALLOCATE(VEC1,VEC2)
    CALL HAMILTONIAN_PRODUCT(50,51,IFILE,2*MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE))
    ALLOCATE(VEC1(NCF,NCF),VEC2(NCF,NCF))
    REWIND(51)
    DO I=0,IFILE
     READ(51) VEC1
    ENDDO
    IF (IOPTN(9) >= 3) THEN
     WRITE(6,'(A,I1,A)') 'H|',IFILE,'>'
     CALL DUMP5(VEC1,NCF)
    ENDIF
!   IF (IFILE == 0) THEN
!    DMP(1)=-DMP(0)
!    DO MOI=1,IOCC
!     DMP(1)=DMP(1)+2.0D0*H(MOI,MOI)
!    ENDDO
!    DO MOI=1,IOCC
!     DO MOJ=1,IOCC
!      DMP(1)=DMP(1)+(2.0D0*G(MOI,MOI,MOJ,MOJ)-G(MOI,MOJ,MOJ,MOI))
!     ENDDO
!    ENDDO
!   ELSE
!    DMP(IFILE+1)=VEC1(1,1)
!   ENDIF
    DMP(IFILE+1)=VEC1(1,1)
    IF (IFILE == 0) DMP(1)=DMP(1)-DMP(0)
    EMP=EMP+DMP(IFILE+1)
!   CALL PCPU_TIME(ICPUE)
    WRITE(6,'(I2,F20.10,F20.10,F12.1)') IFILE+1,DMP(IFILE+1),EMP,ICPUE-ICPUS
    MP_STORED(IFILE+1)=DMP(IFILE+1)
!   CALL PCPU_TIME(ICPUS)
    CALL PFLUSH(6)
    IF (DABS(DMP(IFILE+1)) < 1.0D-15) EXIT
    VEC1=-VEC1
    REWIND(50)
    DO K=0,IFILE
     READ(50) VEC2
     DO IA=1,NCF
      HA=0.0D0
      IF (K == IFILE) THEN
       DO I=1,IALL(0,0,0)-IVIRTCORE
        IF (BTEST(CFHALF(IA),I-1)) HA=HA+EPSILON(I,0,0,0)
       ENDDO
      ENDIF
      DO IB=1,NCF
       HB=0.0D0
       IF (K == IFILE) THEN
        DO I=1,IALL(0,0,0)-IVIRTCORE
         IF (BTEST(CFHALF(IB),I-1)) THEN
!!! For f5, comment in the next 6 comments
          IF (I==MOP) THEN
           HB=HB+OMEGA
          ELSE IF (I==MOQ) THEN
           HB=HB+OMEGA
          ELSE
           HB=HB+EPSILON(I,0,0,0)
          ENDIF
         ENDIF
        ENDDO
       ENDIF
       IF (K == IFILE) THEN
        VEC1(IB,IA)=VEC1(IB,IA)+(DMP(IFILE-K+1)+HA+HB)*VEC2(IB,IA)
       ELSE
        VEC1(IB,IA)=VEC1(IB,IA)+DMP(IFILE-K+1)*VEC2(IB,IA)
       ENDIF
      ENDDO
     ENDDO
    ENDDO
    DO IA=1,NCF
     HA=0.0D0
     DO I=1,IALL(0,0,0)-IVIRTCORE
      IF (BTEST(CFHALF(IA),I-1)) HA=HA+EPSILON(I,0,0,0)
     ENDDO
     DO IB=1,NCF
      HB=0.0D0
      DO I=1,IALL(0,0,0)-IVIRTCORE
       IF (BTEST(CFHALF(IB),I-1)) THEN
!!! For f5, comment in the next 6 comments
        IF (I==MOP) THEN
         HB=HB+OMEGA
        ELSE IF (I==MOQ) THEN
         HB=HB+OMEGA
        ELSE
         HB=HB+EPSILON(I,0,0,0)
        ENDIF
       ENDIF
      ENDDO
      IF ((IA == 1).AND.(IB == 1)) THEN
       VEC2(IB,IA)=0.0D0
      ELSE
       VEC2(IB,IA)=VEC1(IB,IA)/(HA+HB-DMP(0))
      ENDIF
     ENDDO
    ENDDO
    IF (IFILE == MAXFILE) EXIT
    WRITE(50) VEC2
   ENDDO

   WRITE(6,'(A)') '------------------------------------------------------'
!  IF (IOPTN(97) > 0) THEN
!   OPEN(97,FILE=TRIM(COPTN(1))//'.gf0',FORM='UNFORMATTED')
!   VEC2=0.0D0
!   REWIND(50)
!   DO I=0,ORDER-1
!    READ(50) VEC1
!    VEC2=VEC2+VEC1
!   ENDDO
!   REWIND(97)
!   WRITE(97) VEC2
!   CLOSE(97)
!   WRITE(6,'(A,A)') 'REFERENCE WAVE FUNCTION FOR MBGF STORED IN ',TRIM(COPTN(1))//'.gf0'
!  ENDIF

   DEALLOCATE(VEC1,VEC2)
   CLOSE(50)
   CLOSE(51)
   RETURN
END SUBROUTINE



SUBROUTINE HIGHORDER_GF_REDUX_OBSOLETE(OMEGA,ORDER)
! PERFORM HIGH-ORDER MBGF CALCULATIONS

   USE CONTROL
   USE GRADIENT
   USE STRUCTURE
   USE INTEGRAL
   USE BASISSET
   USE FULLCI

   IMPLICIT NONE
   INTEGER,PARAMETER :: MAXFILE = 100
   DOUBLE PRECISION :: OMEGA
   INTEGER :: ORDER
   INTEGER :: I,J,K,L,M,N1,N2,N3
   INTEGER :: PA,PB,QA,QB,RX,RY,II
   INTEGER :: MOX,MOY,MOZ
   DOUBLE PRECISION,ALLOCATABLE :: IPVEC1(:,:),IPVEC2(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: EAVEC1(:,:),EAVEC2(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: WFN(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: IP_WFN(:,:,:,:,:),IP_E(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: EA_WFN(:,:,:,:,:),EA_E(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: GN(:,:,:),SN(:,:,:)
   INTEGER,ALLOCATABLE :: TEMP1(:),TEMP2(:,:)
   INTEGER :: N,NP
   INTEGER,ALLOCATABLE :: NCOMB(:),NRANKS(:,:,:),DCOMB(:),DRANK(:,:),DRANKS(:,:,:)
   DOUBLE PRECISION :: DDD

   WRITE(6,'(A,F20.15,A)') 'SELF-ENERGY MATRIX WILL BE COMPUTED AT OMEGA=',OMEGA,' HARTREE'

   ALLOCATE(IPVEC1(IP_NCF,NCF),IPVEC2(IP_NCF,NCF))
   ALLOCATE(EAVEC1(EA_NCF,NCF),EAVEC2(EA_NCF,NCF))

   ! COMBINATORICS
   NP=0
   DO I=1,ORDER
    NP=NP+ORDER**I
   ENDDO
   ALLOCATE(TEMP1(NP),TEMP2(ORDER,NP))
   ALLOCATE(DCOMB(0:ORDER),DRANK(NP,ORDER),DRANKS(ORDER,NP,ORDER))
   ALLOCATE(NCOMB(0:ORDER),NRANKS(2,NP,ORDER))
   DRANK=0
   DRANKS=0
   NRANKS=0
   DCOMB=0
   NCOMB=0
   DO I=1,ORDER
    DO J=1,I
     CALL COMBINATORICS(1,J,TEMP1,TEMP2,N,NP,ORDER)
    ENDDO
!   WRITE(*,*) 'ORDER=',I,' N=',N
!   DO J=1,N
!    WRITE(*,*) J,TEMP1(J),(TEMP2(K,J),K=1,TEMP1(J))
!   ENDDO
    M=0
    DO J=1,N
     L=0
     DO K=1,TEMP1(J)
      L=L+TEMP2(K,J)
     ENDDO
     IF (L==I) THEN
      M=M+1
      DRANK(M,I)=TEMP1(J)
      DO K=1,TEMP1(J)
       DRANKS(K,M,I)=TEMP2(K,J)
      ENDDO
     ENDIF
    ENDDO
    DCOMB(I)=M
    WRITE(*,*) 'ORDER=',I,' DCOMB=',M
!   DO J=1,M
!    WRITE(*,*) J,'RANK=',DRANK(J,I),'ORDERS=',(DRANKS(K,J,I),K=1,DRANK(J,I))
!   ENDDO
   ENDDO
   DO I=0,ORDER
    M=0
    DO J=0,I
     DO K=0,I
      IF (J+K==I) THEN
       M=M+1
       NRANKS(1,M,I)=J
       NRANKS(2,M,I)=K
      ENDIF
     ENDDO
    ENDDO
    NCOMB(I)=M
    WRITE(*,*) 'ORDER=',I,' NCOMB=',NCOMB(I)
!   DO J=1,NCOMB(I)
!    WRITE(*,*) J,'ORDERS=',(NRANKS(K,J,I),K=1,2)
!   ENDDO
   ENDDO
   DO I=1,ORDER
    N3=0
    DO J=0,I
     K=I-J
     DO N1=1,NCOMB(J)
      DO N2=1,NCOMB(K)
       N3=N3+1
      ENDDO
     ENDDO
    ENDDO
    WRITE(*,*) 'ORDER=',I,' NCOMB*NCOMB=',N3
   ENDDO
   DO I=1,ORDER
!   WRITE(6,'(I3,A)') I,"TH-ORDER GREEN'S FUNCTION IS ..."
    DO J=0,I
     DO K=0,I-J
      L=I-J-K
      DO N1=1,NCOMB(J)
       DO N2=1,NCOMB(K)
        IF (DCOMB(L)==0) THEN
!        WRITE(6,'(A,I1,A,I1,A,I1,A,I1,A)') &
!          '    <',NRANKS(1,N1,J),'|',NRANKS(2,N1,J),'><',NRANKS(2,N2,K),'|',NRANKS(1,N2,K),'>'
        ELSE
         DO M=1,DCOMB(L)
!         WRITE(6,'(A,I1,A,I1,A,I1,A,I1,A,10(A,I1))') &
!          '    <',NRANKS(1,N1,J),'|',NRANKS(2,N1,J),'><',NRANKS(2,N2,K),'|',NRANKS(1,N2,K),'>' &
!           ,(' D',DRANKS(N3,M,L),N3=1,DRANK(M,L))
         ENDDO
        ENDIF
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO

   ! READ ALL N-ELECTRON WAVE FUNCTIONS
   ALLOCATE(WFN(NCF,NCF,0:ORDER-1))
   OPEN(97,FILE=TRIM(COPTN(1))//'.gf0',FORM='UNFORMATTED')
   REWIND(97)
   DO I=0,ORDER-1
    READ(97) WFN(:,:,I)
   ENDDO
   CLOSE(97)
   WRITE(6,'(A,I3,A,A)') '0 THRU',ORDER-1,'TH-ORDER N-ELECTRON WAVE FUNCTIONS READ FROM ',TRIM(COPTN(1))//'.gf0'

   ! READ ALL (N+/-1)-ELECTRON WAVE FUNCTIONS
   ALLOCATE(IP_WFN(IP_NCF,NCF,IP_NCF,NCF,0:ORDER-1))
   OPEN(60,FILE=TRIM(COPTN(1))//'.ip1',FORM='UNFORMATTED')
   REWIND(60)
   DO PA=1,NCF
    DO PB=1,IP_NCF
     DO I=0,ORDER-1
      READ(60) IP_WFN(:,:,PB,PA,I)
     ENDDO
    ENDDO
   ENDDO
   CLOSE(60)
   WRITE(6,'(A,I3,A,A)') '0 THRU',ORDER-1,'TH-ORDER (N-1)-ELECTRON WAVE FUNCTIONS READ FROM ',TRIM(COPTN(1))//'.ip1'
   ALLOCATE(EA_WFN(EA_NCF,NCF,EA_NCF,NCF,0:ORDER-1))
   OPEN(60,FILE=TRIM(COPTN(1))//'.ea1',FORM='UNFORMATTED')
   REWIND(60)
   DO PA=1,NCF
    DO PB=1,EA_NCF
     DO I=0,ORDER-1
      READ(60) EA_WFN(:,:,PB,PA,I)
     ENDDO
    ENDDO
   ENDDO
   CLOSE(60)
   WRITE(6,'(A,I3,A,A)') '0 THRU',ORDER-1,'TH-ORDER (N+1)-ELECTRON WAVE FUNCTIONS READ FROM ',TRIM(COPTN(1))//'.ea1'

   ! READ ALL STATIC SIGMA
   ALLOCATE(IP_E(IP_NCF,NCF,0:ORDER-1))
   OPEN(60,FILE=TRIM(COPTN(1))//'.ip0',FORM='UNFORMATTED')
   REWIND(60)
   DO PA=1,NCF
    DO PB=1,IP_NCF
     DO I=0,ORDER-1
      READ(60) IP_E(PB,PA,I)
     ENDDO
    ENDDO
   ENDDO
   CLOSE(60)
   WRITE(6,'(A,A)') '(N-1)-ELECTRON STATIC SIGMA READ FROM ',TRIM(COPTN(1))//'.ip0'
   ALLOCATE(EA_E(EA_NCF,NCF,0:ORDER-1))
   OPEN(60,FILE=TRIM(COPTN(1))//'.ea0',FORM='UNFORMATTED')
   REWIND(60)
   DO PA=1,NCF
    DO PB=1,EA_NCF
     DO I=0,ORDER-1
      READ(60) EA_E(PB,PA,I)
     ENDDO
    ENDDO
   ENDDO
   CLOSE(60)
   WRITE(6,'(A,A)') '(N+1)-ELECTRON STATIC SIGMA READ FROM ',TRIM(COPTN(1))//'.ea0'

!--------------------------- check orthogonality
goto 1
do pa=1,ncf
do pb=1,ip_ncf
write(*,*) 'pa,pb=',pa,pb

!goto 2

REWIND(50)
WRITE(50) WFN(:,:,0)
write(*,*) '0-order N-e wfn'
call dump5(wfn(:,:,0),ncf)
CALL IONIZE(50,51,pa,pb,.FALSE.,0)
REWIND(51)
READ(51) IPVEC1
write(*,*) '0-order (N-1)-e wfn'
call dump16(wfn(:,:,0),ip_ncf,ncf)

write(*,*) 'overlap with 0-order (N-1)-e all wfns'
do qa=1,ncf
 do qb=1,ip_ncf
  ddd=0.0d0
  do rx=1,ncf
   do ry=1,ip_ncf
    ddd=ddd+ipvec1(ry,rx)*ip_wfn(ry,rx,qb,qa,0)
   enddo
  enddo
  if (dabs(ddd) > 1.0d-10) write(*,*) qa,qb,ddd
 enddo
enddo

write(*,*) 'overlap between 0-order GS and 1-order (N-1)-e all wfns'
do qa=1,ncf
 do qb=1,ip_ncf
  ddd=0.0d0
  do rx=1,ncf
   do ry=1,ip_ncf
    ddd=ddd+ipvec1(ry,rx)*ip_wfn(ry,rx,qb,qa,1)
   enddo
  enddo
  if (dabs(ddd) > 1.0d-10) write(*,*) qa,qb,ddd
 enddo
enddo

write(*,*) 'overlap between 1-order (N-1)-e GS and all wfns'
do qa=1,ncf
 do qb=1,ip_ncf
  ddd=0.0d0
  do rx=1,ncf
   do ry=1,ip_ncf
    ddd=ddd+ip_wfn(ry,rx,pb,pa,1)*ip_wfn(ry,rx,qb,qa,1)
   enddo
  enddo
  if (dabs(ddd) > 1.0d-10) write(*,*) qa,qb,ddd
 enddo
enddo

write(*,*) 'overlap between 0-order GS and all-order GS wfns'
do ii=0,order-1
  ddd=0.0d0
  do rx=1,ncf
   do ry=1,ip_ncf
    ddd=ddd+ip_wfn(ry,rx,pb,pa,0)*ip_wfn(ry,rx,pb,pa,ii)
   enddo
  enddo
  write(*,*) ii,ddd
enddo

write(*,*) 'overlap between 1-order GS and all-order GS wfns'
do ii=0,order-1
  ddd=0.0d0
  do rx=1,ncf
   do ry=1,ip_ncf
    ddd=ddd+ip_wfn(ry,rx,pb,pa,1)*ip_wfn(ry,rx,pb,pa,ii)
   enddo
  enddo
  write(*,*) ii,ddd
enddo

2 continue 

write(*,*) 'calculate all-order energies'
OPEN(50,FILE=TRIM(COPTN(1))//'.fi0',FORM='UNFORMATTED')
OPEN(51,FILE=TRIM(COPTN(1))//'.fo0',FORM='UNFORMATTED')
do ii=0,order-1
REWIND(50)
WRITE(50) ip_wfn(:,:,pb,pa,ii)
CALL IP_HAMILTONIAN_PRODUCT(50,51,0,MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE)+MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE+1))
REWIND(51)
READ(51) ipvec1
write(*,*) ii,ipvec1(pb,pa),ip_e(pb,pa,ii)
enddo
close(50)
close(51)

enddo
enddo
1 continue
!---------------------------

   ALLOCATE(GN(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,ORDER))
   GN=0.0D0
   OPEN(50,FILE=TRIM(COPTN(1))//'.fi0',FORM='UNFORMATTED')
   OPEN(51,FILE=TRIM(COPTN(1))//'.fo0',FORM='UNFORMATTED')

   DO I=1,ORDER-1

    WRITE(6,'(A,I3)') "GREEN'S FUNCTION ORDER=",I

    DO J=0,I
     DO K=0,I-J
      L=I-J-K
      DO N1=1,NCOMB(J)
       DO N2=1,NCOMB(K)

        IF (DCOMB(L)==0) THEN

         WRITE(6,'(A,I1,A,I1,A,I1,A,I1,A)') &
           '    <',NRANKS(1,N1,J),'|',NRANKS(2,N1,J),'><',NRANKS(2,N2,K),'|',NRANKS(1,N2,K),'>'

         DO RX=1,IP_NCF+EA_NCF
          DO RY=1,IP_NCF+EA_NCF

           IF ((RX <= IP_NCF).AND.(IP_NORDER(RX)==1)) THEN
            DO II=1,IOCC
             IF (.NOT.BTEST(IP_CFHALF(RX),II-1)) MOX=II
            ENDDO
           ELSE IF (EA_NORDER(RX-IP_NCF)==1) THEN
            DO II=IOCC+1,IALL(0,0,0)-IVIRTCORE
             IF (BTEST(EA_CFHALF(RX-IP_NCF),II-1)) MOX=II
            ENDDO
           ELSE
            CYCLE
           ENDIF
            
           IF ((RY <= IP_NCF).AND.(IP_NORDER(RY)==1)) THEN
            DO II=1,IOCC
             IF (.NOT.BTEST(IP_CFHALF(RY),II-1)) MOY=II
            ENDDO
           ELSE IF (EA_NORDER(RY-IP_NCF)==1) THEN
            DO II=IOCC+1,IALL(0,0,0)-IVIRTCORE
             IF (BTEST(EA_CFHALF(RY-IP_NCF),II-1)) MOY=II
            ENDDO
           ELSE
            CYCLE
           ENDIF

           ! IP SECTOR
           REWIND(50)
           WRITE(50) WFN(:,:,NRANKS(1,N1,J))
           CALL IONIZE(50,51,1,RX,.FALSE.,0)
           REWIND(51)
           READ(51) IPVEC1
           REWIND(50)
           WRITE(50) WFN(:,:,NRANKS(1,N2,K))
           CALL IONIZE(50,51,1,RY,.FALSE.,0)
           REWIND(51)
           READ(51) IPVEC2

           DO PA=1,NCF
            DO PB=1,IP_NCF
             DO QA=1,NCF
              DO QB=1,IP_NCF
               GN(MOX,MOY,I)=GN(MOX,MOY,I) &
                +IPVEC1(PB,PA)*IP_WFN(PB,PA,QB,QA,NRANKS(2,N1,J)) &
                *IPVEC2(PB,PA)*IP_WFN(PB,PA,QB,QA,NRANKS(2,N2,K)) 
!               /(OMEGA-IP_E(QB,QA,0))
              ENDDO
             ENDDO
            ENDDO
           ENDDO
   
           ! EA SECTOR
           REWIND(50)
           WRITE(50) WFN(:,:,NRANKS(1,N1,J))
           CALL IONIZE(50,51,1,RX,.TRUE.,0)
           REWIND(51)
           READ(51) EAVEC1
           REWIND(50)
           WRITE(50) WFN(:,:,NRANKS(1,N2,K))
           CALL IONIZE(50,51,1,RY,.TRUE.,0)
           REWIND(51)
           READ(51) EAVEC2
   
           DO PA=1,NCF
            DO PB=1,EA_NCF
             DO QA=1,NCF
              DO QB=1,EA_NCF
               GN(MOX,MOY,I)=GN(MOX,MOY,I) &
                +EAVEC1(PB,PA)*EA_WFN(PB,PA,QB,QA,NRANKS(2,N1,J)) &
                *EAVEC2(PB,PA)*EA_WFN(PB,PA,QB,QA,NRANKS(2,N2,K)) 
!               /(OMEGA-EA_E(QB,QA,0))
              ENDDO
             ENDDO
            ENDDO
           ENDDO
 
          ENDDO
         ENDDO
!call dump5(GN(:,:,I),iall(0,0,0)-ivirtcore)

        ELSE

         DO M=1,DCOMB(L)
cycle
          WRITE(6,'(A,I1,A,I1,A,I1,A,I1,A,10(A,I1))') &
           '    <',NRANKS(1,N1,J),'|',NRANKS(2,N1,J), &
           '><',NRANKS(2,N2,K),'|',NRANKS(1,N2,K),'>' &
           ,(' D',DRANKS(N3,M,L),N3=1,DRANK(M,L))
 
          DO RX=1,IP_NCF+EA_NCF
           DO RY=1,IP_NCF+EA_NCF

            IF ((RX <= IP_NCF).AND.(IP_NORDER(RX)==1)) THEN
             DO II=1,IOCC
              IF (.NOT.BTEST(IP_CFHALF(RX),II-1)) MOX=II
             ENDDO
            ELSE IF (EA_NORDER(RX-IP_NCF)==1) THEN
             DO II=IOCC+1,IALL(0,0,0)-IVIRTCORE
              IF (BTEST(EA_CFHALF(RX-IP_NCF),II-1)) MOX=II
             ENDDO
            ELSE
             CYCLE
            ENDIF

            IF ((RY <= IP_NCF).AND.(IP_NORDER(RY)==1)) THEN
             DO II=1,IOCC
              IF (.NOT.BTEST(IP_CFHALF(RY),II-1)) MOY=II
             ENDDO
            ELSE IF (EA_NORDER(RY-IP_NCF)==1) THEN
             DO II=IOCC+1,IALL(0,0,0)-IVIRTCORE
              IF (BTEST(EA_CFHALF(RY-IP_NCF),II-1)) MOY=II
             ENDDO
            ELSE
             CYCLE
            ENDIF
   
            ! IP SECTOR
            REWIND(50)
            WRITE(50) WFN(:,:,NRANKS(1,N1,J))
            CALL IONIZE(50,51,1,RX,.FALSE.,0)
            REWIND(51)
            READ(51) IPVEC1
            REWIND(50)
            WRITE(50) WFN(:,:,NRANKS(1,N2,K))
            CALL IONIZE(50,51,1,RY,.FALSE.,0)
            REWIND(51)
            READ(51) IPVEC2
   
            DO PA=1,NCF
             DO PB=1,IP_NCF
              DO QA=1,NCF
               DO QB=1,IP_NCF
                DDD=IPVEC1(PB,PA)*IP_WFN(PB,PA,QB,QA,NRANKS(2,N1,J)) &
                   *IPVEC2(PB,PA)*IP_WFN(PB,PA,QB,QA,NRANKS(2,N2,K)) &
                 /(OMEGA-IP_E(QB,QA,0))
                DO N3=1,DRANK(M,L)
                 DDD=DDD*IP_E(QB,QA,DRANKS(N3,M,L))/(OMEGA-IP_E(QB,QA,0))
                ENDDO
                GN(MOX,MOY,I)=GN(MOX,MOY,I)+DDD
               ENDDO
              ENDDO
             ENDDO
            ENDDO
   
            ! EA SECTOR
            REWIND(50)
            WRITE(50) WFN(:,:,NRANKS(1,N1,J))
            CALL IONIZE(50,51,1,RX,.TRUE.,0)
            REWIND(51)
            READ(51) EAVEC1
            REWIND(50)
            WRITE(50) WFN(:,:,NRANKS(1,N2,K))
            CALL IONIZE(50,51,1,RY,.TRUE.,0)
            REWIND(51)
            READ(51) EAVEC2
    
            DO PA=1,NCF
             DO PB=1,EA_NCF
              DO QA=1,NCF
               DO QB=1,EA_NCF
                DDD=EAVEC1(PB,PA)*EA_WFN(PB,PA,QB,QA,NRANKS(2,N1,J)) &
                   *EAVEC2(PB,PA)*EA_WFN(PB,PA,QB,QA,NRANKS(2,N2,K)) &
                 /(OMEGA-EA_E(QB,QA,0))
                DO N3=1,DRANK(M,L)
                 DDD=DDD*EA_E(QB,QA,DRANKS(N3,M,L))/(OMEGA-EA_E(QB,QA,0))
                ENDDO
                GN(MOX,MOY,I)=GN(MOX,MOY,I)+DDD
               ENDDO
              ENDDO
             ENDDO
            ENDDO
    
           ENDDO
          ENDDO
         ENDDO

        ENDIF

       ENDDO
      ENDDO
     ENDDO
    ENDDO

    WRITE(6,'(I3,A)') I,"TH-ORDER GREEN'S FUNCTION"
    CALL DUMP5(GN(:,:,I),IALL(0,0,0)-IVIRTCORE)

   ENDDO

   CLOSE(50)
   CLOSE(51)

   ALLOCATE(SN(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,ORDER))

   DO I=1,ORDER-1

    WRITE(6,'(A,I3)') "SELF-ENERGY ORDER=",I

    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
     DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
      SN(MOX,MOY,I)=(OMEGA-EPSILON(MOX,0,0,0))*GN(MOX,MOY,I)*(OMEGA-EPSILON(MOY,0,0,0))
     ENDDO
    ENDDO

    IF (I > 1) THEN
     DO J=1,I-1
      DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
       DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
        DO MOZ=ICORE+1,IALL(0,0,0)-IVIRTCORE
         SN(MOX,MOY,I)=SN(MOX,MOY,I)-SN(MOX,MOZ,J)*GN(MOZ,MOY,I-J)*(OMEGA-EPSILON(MOY,0,0,0))
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDIF

!   WRITE(6,'(I3,A)') I,"TH-ORDER SELF ENERGY"
!   CALL DUMP5(SN(:,:,I),IALL(0,0,0)-IVIRTCORE)

   ENDDO

   DEALLOCATE(IPVEC1,IPVEC2)
   DEALLOCATE(EAVEC1,EAVEC2)
   DEALLOCATE(GN,SN)

   RETURN

END SUBROUTINE



RECURSIVE SUBROUTINE COMBINATORICS(IORDER,JORDER,A,B,N,NP,ORDER)

   IMPLICIT NONE
   INTEGER :: IORDER,JORDER,N,NP,ORDER
   INTEGER :: A(NP),B(ORDER,NP)
   INTEGER :: I,J,K

   K=-1
   DO I=0,JORDER-1
    K=K+ORDER**I
   ENDDO
   IF ((ORDER**JORDER)+K > NP) CALL PABORT('A BUG')
   DO I=1,ORDER**JORDER
    J=(I-1)/(ORDER**(IORDER-1))-ORDER*((I-1)/ORDER**IORDER)+1
    A(K+I)=JORDER
    B(IORDER,K+I)=J
   ENDDO
   N=K+ORDER**JORDER
   IF (IORDER < JORDER) CALL COMBINATORICS(IORDER+1,JORDER,A,B,N,NP,ORDER)

   RETURN
END SUBROUTINE



SUBROUTINE HIGHORDER_GF_REDUX_OBSOLETE2(OMEGA,ORDER)
! PERFORM HIGH-ORDER MBGF CALCULATIONS

   USE CONTROL
   USE GRADIENT
   USE STRUCTURE
   USE INTEGRAL
   USE BASISSET
   USE FULLCI

   IMPLICIT NONE
   DOUBLE PRECISION,PARAMETER :: TOL = 1.0D-10
   INTEGER,PARAMETER :: MAXFILE = 100
   DOUBLE PRECISION :: OMEGA
   INTEGER :: ORDER
   INTEGER :: I,J,K,N1,N2,N3,N4
   INTEGER :: PA,PB,QA,QB,RA,RB,XA,XB,YA,YB
   INTEGER :: MOX,MOY,MOZ
   DOUBLE PRECISION,ALLOCATABLE :: GN(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: SN(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: IPVEC1(:,:),IPVEC2(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: EAVEC1(:,:),EAVEC2(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: WFN(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: IP_WFN(:,:,:,:,:),IP_E(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: EA_WFN(:,:,:,:,:),EA_E(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: NIP(:,:,:,:,:),NEA(:,:,:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: DIP(:,:,:),DEA(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: GIP(:,:,:,:,:),GEA(:,:,:,:,:)
   INTEGER :: N,NP
   DOUBLE PRECISION :: T1,T2
   INTEGER :: SGN
   DOUBLE PRECISION,ALLOCATABLE :: VEC1(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: HMAT(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: EIP(:),EEA(:),VIP(:,:),VEA(:,:)
double precision :: t3,t4,t5

   WRITE(6,'(A,F20.15,A)') 'SELF-ENERGY MATRIX WILL BE COMPUTED AT OMEGA=',OMEGA,' HARTREE'

   ALLOCATE(IPVEC1(IP_NCF,NCF),IPVEC2(IP_NCF,NCF))
   ALLOCATE(EAVEC1(EA_NCF,NCF),EAVEC2(EA_NCF,NCF))

   ! FORM VIP,VEA
   ALLOCATE(EIP(NCF*IP_NCF))
   ALLOCATE(EEA(NCF*EA_NCF))
   ALLOCATE(VIP(NCF*IP_NCF,NCF*IP_NCF))
   ALLOCATE(VEA(NCF*EA_NCF,NCF*EA_NCF))

   ! IP SECTOR
   DO PA=1,NCF
    DO PB=1,IP_NCF
     IPVEC1=0.0D0
     IPVEC1(PB,PA)=1.0D0
     REWIND(50)
     WRITE(50) IPVEC1
     CALL IP_HAMILTONIAN_PRODUCT(50,51,0,MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE)+MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE+1))
     REWIND(51)
     READ(51) IPVEC1
     CALL H0_PRODUCT(50,51,1)
     REWIND(51)
     READ(51) IPVEC2
     EIP((PA-1)*IP_NCF+PB)=IPVEC2(PB,PA)
     IPVEC1=IPVEC1-IPVEC2 ! V=H-H0
     DO QA=1,NCF
      DO QB=1,IP_NCF
       VIP((QA-1)*IP_NCF+QB,(PA-1)*IP_NCF+PB)=IPVEC1(QB,QA)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
!  WRITE(6,*) 'EIP VECTOR'
!  CALL DUMP16(EIP,NCF*IP_NCF,1)
!  WRITE(6,*) 'VIP MATRIX'
!  CALL DUMP5(VIP,NCF*IP_NCF)

   ! EA SECTOR
   DO PA=1,NCF
    DO PB=1,EA_NCF
     EAVEC1=0.0D0
     EAVEC1(PB,PA)=1.0D0
     REWIND(50)
     WRITE(50) EAVEC1
     CALL EA_HAMILTONIAN_PRODUCT(50,51,0,MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE)+MIN(IOCC-ICORE+1,IALL(0,0,0)-IOCC-IVIRTCORE))
     REWIND(51)
     READ(51) EAVEC1
     CALL H0_PRODUCT(50,51,-1)
     REWIND(51)
     READ(51) EAVEC2
     EEA((PA-1)*EA_NCF+PB)=EAVEC2(PB,PA)
     EAVEC1=EAVEC1-EAVEC2 ! V=H-H0
     DO QA=1,NCF
      DO QB=1,EA_NCF
       VEA((QA-1)*EA_NCF+QB,(PA-1)*EA_NCF+PB)=EAVEC1(QB,QA)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
!  WRITE(6,*) 'EEA VECTOR'
!  CALL DUMP16(EEA,NCF*EA_NCF,1)
!  WRITE(6,*) 'VEA MATRIX'
!  CALL DUMP5(VEA,NCF*EA_NCF)

   ! READ ALL STATIC SIGMA
   ALLOCATE(IP_E(IP_NCF,NCF,0:ORDER-1))
   OPEN(60,FILE=TRIM(COPTN(1))//'.ip0',FORM='UNFORMATTED')
   REWIND(60)
   DO PA=1,NCF
    DO PB=1,IP_NCF
     DO I=0,ORDER-1
      READ(60) IP_E(PB,PA,I)
     ENDDO
!write(*,*) (IP_E(PB,PA,I),i=0,order-1)
    ENDDO
   ENDDO
   CLOSE(60)
   WRITE(6,'(A,A)') '(N-1)-ELECTRON STATIC SIGMA READ FROM ',TRIM(COPTN(1))//'.ip0'
   ALLOCATE(EA_E(EA_NCF,NCF,0:ORDER-1))
   OPEN(60,FILE=TRIM(COPTN(1))//'.ea0',FORM='UNFORMATTED')
   REWIND(60)
   DO PA=1,NCF
    DO PB=1,EA_NCF
     DO I=0,ORDER-1
      READ(60) EA_E(PB,PA,I)
     ENDDO
!write(*,*) (ea_E(PB,PA,I),i=0,order-1)
    ENDDO
   ENDDO
   CLOSE(60)
   WRITE(6,'(A,A)') '(N+1)-ELECTRON STATIC SIGMA READ FROM ',TRIM(COPTN(1))//'.ea0'

   ! READ ALL N-ELECTRON WAVE FUNCTIONS
   ALLOCATE(WFN(NCF,NCF,0:ORDER-1))
   OPEN(97,FILE=TRIM(COPTN(1))//'.gf0',FORM='UNFORMATTED')
   REWIND(97)
   DO I=0,ORDER-1
    READ(97) WFN(:,:,I)
   ENDDO
   CLOSE(97)
   WRITE(6,'(A,I3,A,A)') '0 THRU',ORDER-1,'TH-ORDER N-ELECTRON WAVE FUNCTIONS READ FROM ',TRIM(COPTN(1))//'.gf0'
   ! INTERMEDIATE NORMALIZATION SEE EQ(2.28) OF SHAVITT & BARTLETT
   DO I=1,ORDER-1
    T1=0.0D0
    DO PA=1,NCF
     DO PB=1,NCF
      T1=T1+WFN(PB,PA,I)*WFN(PB,PA,0)
     ENDDO
    ENDDO
    IF (DABS(T1) > TOL) THEN
     WRITE(6,*) I,T1
     CALL WARNING('NONORTHOGONAL N')
    ENDIF
   ENDDO
   ! TRIPLE CHECK
   OPEN(50,FILE=TRIM(COPTN(1))//'.fi0',FORM='UNFORMATTED')
   OPEN(51,FILE=TRIM(COPTN(1))//'.fo0',FORM='UNFORMATTED')
   ALLOCATE(VEC1(NCF,NCF))
   DO I=0,ORDER-1
    VEC1=0.0D0
    DO J=0,I
     VEC1(:,:)=VEC1(:,:)+WFN(:,:,J)
    ENDDO
    REWIND(50)
    WRITE(50) VEC1
    CALL HAMILTONIAN_PRODUCT(50,51,0,2*MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE))
    REWIND(51)
    READ(51) VEC1
    T1=0.0D0
    DO PA=1,NCF
     DO PB=1,NCF
      T1=T1+VEC1(PB,PA)*WFN(PB,PA,0)
     ENDDO
    ENDDO
    IF (DABS(T1+NUCLEAR_REPULSION-MP_STORED(I+1)) > TOL) THEN
     WRITE(6,'(I3,2F20.15)') I, T1+NUCLEAR_REPULSION,MP_STORED(I+1)
     CALL PABORT('ENERGY INCONSISTENCY N')
    ENDIF
   ENDDO
   DEALLOCATE(VEC1)
   CLOSE(50)
   CLOSE(51)

   ! READ ALL (N-1)-ELECTRON WAVE FUNCTIONS
   ALLOCATE(IP_WFN(IP_NCF,NCF,IP_NCF,NCF,0:ORDER-1))
   OPEN(60,FILE=TRIM(COPTN(1))//'.ip1',FORM='UNFORMATTED')
   REWIND(60)
   DO PA=1,NCF
    DO PB=1,IP_NCF
     DO I=0,ORDER-1
      READ(60) IP_WFN(:,:,PB,PA,I)
     ENDDO
!write(*,*) pb,pa
!call dump16(ip_wfn(:,:,pb,pa,0),ip_ncf,ncf)
    ENDDO
   ENDDO
   CLOSE(60)
   WRITE(6,'(A,I3,A,A)') '0 THRU',ORDER-1,'TH-ORDER (N-1)-ELECTRON WAVE FUNCTIONS READ FROM ',TRIM(COPTN(1))//'.ip1'
   ! INTERMEDIATE NORMALIZATION SEE EQ(2.28) OF SHAVITT & BARTLETT
!  ALLOCATE(HMAT(NCF*IP_NCF,NCF*IP_NCF))
!  DO I=0,ORDER-1
!   DO RA=1,NCF
!    DO RB=1,IP_NCF
!   DO PA=1,NCF
!    DO PB=1,IP_NCF
!     T1=0.0D0
!     DO XA=1,NCF
!      DO XB=1,IP_NCF
!       T1=T1+IP_WFN(XB,XA,RB,RA,I)*IP_WFN(XB,XA,PB,PA,0)
!     DO YA=1,NCF
!      DO YB=1,IP_NCF
!       T1=T1+IP_WFN(YB,YA,RB,RA,I)*(VIP((YA-1)*IP_NCF+YB,(XA-1)*IP_NCF+XB))*IP_WFN(XB,XA,PB,PA,0)
!      ENDDO
!     ENDDO
!      ENDDO
!     ENDDO
!     HMAT((RA-1)*IP_NCF+RB,(PA-1)*IP_NCF+PB)=T1
!!    IF (DABS(T1) > TOL) THEN
!!     WRITE(6,*) I,T1
!!     CALL WARNING('NONORTHOGONAL N-1')
!!    ENDIF
!    ENDDO
!   ENDDO
!    ENDDO
!   ENDDO
!write(*,*) 'Order=',I
!call dump5(hmat,ncf*ip_ncf)
!  ENDDO
!  DEALLOCATE(HMAT)

   ! READ ALL (N+1)-ELECTRON WAVE FUNCTIONS
   ALLOCATE(EA_WFN(EA_NCF,NCF,EA_NCF,NCF,0:ORDER-1))
   OPEN(60,FILE=TRIM(COPTN(1))//'.ea1',FORM='UNFORMATTED')
   REWIND(60)
   DO PA=1,NCF
    DO PB=1,EA_NCF
     DO I=0,ORDER-1
      READ(60) EA_WFN(:,:,PB,PA,I)
     ENDDO
!write(*,*) pb,pa
!call dump16(ea_wfn(:,:,pb,pa,0),ea_ncf,ncf)
    ENDDO
   ENDDO
   CLOSE(60)
   WRITE(6,'(A,I3,A,A)') '0 THRU',ORDER-1,'TH-ORDER (N+1)-ELECTRON WAVE FUNCTIONS READ FROM ',TRIM(COPTN(1))//'.ea1'
   ! INTERMEDIATE NORMALIZATION SEE EQ(2.28) OF SHAVITT & BARTLETT
   DO I=1,ORDER-1
    DO PA=1,NCF
     DO PB=1,EA_NCF
      T1=0.0D0
      DO QA=1,NCF
       DO QB=1,EA_NCF
        T1=T1+EA_WFN(QB,QA,PB,PA,I)*EA_WFN(QB,QA,PB,PA,0)
       ENDDO
      ENDDO
      IF (DABS(T1) > TOL) THEN
       WRITE(6,*) I,T1
       CALL WARNING('NONORTHOGONAL N+1')
      ENDIF
     ENDDO
    ENDDO
   ENDDO

   ! CONSTRUCT NIP/NEA AND GIP/GEA
   ALLOCATE(NIP(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,IP_NCF,NCF,0:ORDER-1))
   ALLOCATE(NEA(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,EA_NCF,NCF,0:ORDER-1))
   ALLOCATE(DIP(IP_NCF,NCF,0:ORDER-1))
   ALLOCATE(DEA(EA_NCF,NCF,0:ORDER-1))
   
!write(*,*) "================= G WITHOUT DENOM, G0=1 and G1=G2=...=0"
   NIP=0.0D0
   NEA=0.0D0
   DIP=0.0D0
   DEA=0.0D0

   DO I=0,ORDER-1
!   WRITE(6,'(A,I3)') 'ORDER=',I

!t5=0.0d0
    DO N1=0,I
     DO N2=0,I
      DO N3=0,I
       DO N4=0,I
        IF (N1+N2+N3+N4 /= I) CYCLE
!IF (.NOT.(N3==0)) CYCLE
!       IF (.NOT.((N2==0).AND.(N3==0))) CYCLE
!       IF (.NOT.((N1==0).AND.(N4==0))) CYCLE
!       WRITE(6,'(A,I1,A,I1,A,I1,A,I1,A)') '    <',N1,'|',N2,'><',N3,'|',N4,'>'

!if (i==0) write(*,*) n1,n2,n3,n4
        DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
         DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
!if (i==0) write(*,*) 'mox,moy=',mox,moy

          ! IP SECTOR
! - Old ionizer
!         REWIND(50)
!         WRITE(50) WFN(:,:,N1)
!         XA=1
!         IF (MOX <= IOCC) THEN
!          XB=IP_ADDRSS(IBCLR(CFHALF(1),MOX-1))
!         ELSE
!          XB=IP_NCF+EA_ADDRSS(IBSET(CFHALF(1),MOX-1))
!         ENDIF
!         CALL IONIZE(50,51,XA,XB,.FALSE.,0)
!         REWIND(51)
!         READ(51) IPVEC1
!         REWIND(50)
!         WRITE(50) WFN(:,:,N4)
!         YA=1
!         IF (MOY <= IOCC) THEN
!          YB=IP_ADDRSS(IBCLR(CFHALF(1),MOY-1))
!         ELSE
!          YB=IP_NCF+EA_ADDRSS(IBSET(CFHALF(1),MOY-1))
!         ENDIF
!         CALL IONIZE(50,51,YA,YB,.FALSE.,0)
!         REWIND(51)
!         READ(51) IPVEC2
! - New ionizer
          IPVEC1=0.0D0
          DO PA=1,NCF
           DO PB=1,NCF
            IF (.NOT.BTEST(CFHALF(PB),MOX-1)) CYCLE
            QA=PA
            QB=IP_ADDRSS(IBCLR(CFHALF(PB),MOX-1))
            SGN=1
            IF (MOX /= 1) THEN
             DO J=1,MOX-1
              IF (BTEST(CFHALF(PB),J-1)) SGN=-SGN
             ENDDO
            ENDIF
            IPVEC1(QB,QA)=IPVEC1(QB,QA)+WFN(PB,PA,N1)*DFLOAT(SGN)
           ENDDO
          ENDDO
!if ((moy==icore+1).and.(n2+n3+n4==0)) then
!write(*,*) 'mox=',mox,' n1=',n1
!call dump16(ipvec1,ip_ncf,ncf)
!endif
          IPVEC2=0.0D0
          DO PA=1,NCF
           DO PB=1,NCF
            IF (.NOT.BTEST(CFHALF(PB),MOY-1)) CYCLE
            QA=PA
            QB=IP_ADDRSS(IBCLR(CFHALF(PB),MOY-1))
            SGN=1
            IF (MOY /= 1) THEN
             DO J=1,MOY-1
              IF (BTEST(CFHALF(PB),J-1)) SGN=-SGN
             ENDDO
            ENDIF
            IPVEC2(QB,QA)=IPVEC2(QB,QA)+WFN(PB,PA,N4)*DFLOAT(SGN)
           ENDDO
          ENDDO
! - end 
!if (i==0) call dump16(ipvec1,ip_ncf,ncf)
!if (i==0) call dump16(ipvec2,ip_ncf,ncf)

!t3=0.0d0
!t4=0.0d0
!if ((mox==2).and.(moy==icore+1).and.(n3+n4==0)) then
!write(*,*) 'mox=',mox,' n1=',n1,' n2=',n2
!write(*,*) 'ipvec1'
!DO QA=1,NCF
!DO QB=1,IP_NCF
!write(*,'(f15.10)') ipvec1(qb,qa)
!enddo
!enddo
!write(*,*) 'ipvec1 x ip_wfn'
!endif
          DO PA=1,NCF
           DO PB=1,IP_NCF
!if ((n2+n3==0).and.(norder(pA)+ip_norder(pb) == 1)) cycle
            T1=0.0D0
            DO QA=1,NCF
             DO QB=1,IP_NCF
              T1=T1+IPVEC1(QB,QA)*IP_WFN(QB,QA,PB,PA,N2)
             ENDDO
            ENDDO
!if ((mox==2).and.(moy==icore+1).and.(n3+n4==0)) then
!write(*,'(F15.10)') t1
!endif
            T2=0.0D0
            DO QA=1,NCF
             DO QB=1,IP_NCF
              T2=T2+IPVEC2(QB,QA)*IP_WFN(QB,QA,PB,PA,N3)
             ENDDO
            ENDDO
            NIP(MOX,MOY,PB,PA,I)=NIP(MOX,MOY,PB,PA,I)+T1*T2/(OMEGA-IP_E(PB,PA,0))
!           NIP(MOX,MOY,PB,PA,I)=NIP(MOX,MOY,PB,PA,I)+T1*T2
!if ((i<=9).and.(mox==2).and.(moy==3)) then
!t3=t3+t1*t2
!write(*,'(a,2i3,2f20.15)') 'IP pa,pb,t1*t2 =',pa,pb,t1*t2,nip(mox,moy,pb,pa,i)
!write(*,*) 't1,t2=',t1,t2
!call dump5(nip(:,:,pb,pa,i),iall(0,0,0)-ivirtcore)
!endif
!if (moy==4) stop
           ENDDO
          ENDDO
!if ((i<=9).and.(mox==2).and.(moy==3)) then
!write(*,'(a,f20.15)') 'IP sum t1*t2 =',t3
!t5=t5+t3
!endif

          ! EA SECTOR
! - Old electron attacher
!         REWIND(50)
!         WRITE(50) WFN(:,:,N4)
!         XA=1
!         IF (MOX <= IOCC) THEN
!          XB=IP_ADDRSS(IBCLR(CFHALF(1),MOX-1))
!         ELSE
!          XB=IP_NCF+EA_ADDRSS(IBSET(CFHALF(1),MOX-1))
!         ENDIF
!         CALL IONIZE(50,51,XA,XB,.TRUE.,0)
!         REWIND(51)
!         READ(51) EAVEC1
!         REWIND(50)
!         WRITE(50) WFN(:,:,N1)
!         YA=1
!         IF (MOY <= IOCC) THEN
!          YB=IP_ADDRSS(IBCLR(CFHALF(1),MOY-1))
!         ELSE
!          YB=IP_NCF+EA_ADDRSS(IBSET(CFHALF(1),MOY-1))
!         ENDIF
!         CALL IONIZE(50,51,YA,YB,.TRUE.,0)
!         REWIND(51)
!         READ(51) EAVEC2
! - New electron attacher
          EAVEC1=0.0D0
          DO PA=1,NCF
           DO PB=1,NCF
            IF (BTEST(CFHALF(PB),MOX-1)) CYCLE
            QA=PA
            QB=EA_ADDRSS(IBSET(CFHALF(PB),MOX-1))
            SGN=1
            IF (MOX /= 1) THEN
             DO J=1,MOX-1
              IF (BTEST(CFHALF(PB),J-1)) SGN=-SGN
             ENDDO
            ENDIF
            EAVEC1(QB,QA)=EAVEC1(QB,QA)+WFN(PB,PA,N4)*DFLOAT(SGN)
           ENDDO
          ENDDO
          EAVEC2=0.0D0
          DO PA=1,NCF
           DO PB=1,NCF
            IF (BTEST(CFHALF(PB),MOY-1)) CYCLE
            QA=PA
            QB=EA_ADDRSS(IBSET(CFHALF(PB),MOY-1))
            SGN=1
            IF (MOY /= 1) THEN
             DO J=1,MOY-1
              IF (BTEST(CFHALF(PB),J-1)) SGN=-SGN
             ENDDO
            ENDIF
            EAVEC2(QB,QA)=EAVEC2(QB,QA)+WFN(PB,PA,N1)*DFLOAT(SGN)
           ENDDO
          ENDDO
! - end

          DO PA=1,NCF
           DO PB=1,EA_NCF
!if ((n2+n3==0).and.(norder(pA)+ea_norder(pb) == 1)) cycle
            T1=0.0D0
            DO QA=1,NCF
             DO QB=1,EA_NCF
              T1=T1+EAVEC2(QB,QA)*EA_WFN(QB,QA,PB,PA,N2)
             ENDDO
            ENDDO
            T2=0.0D0
            DO QA=1,NCF
             DO QB=1,EA_NCF
              T2=T2+EAVEC1(QB,QA)*EA_WFN(QB,QA,PB,PA,N3)
             ENDDO
            ENDDO
            NEA(MOX,MOY,PB,PA,I)=NEA(MOX,MOY,PB,PA,I)+T1*T2/(OMEGA-EA_E(PB,PA,0))
!           NEA(MOX,MOY,PB,PA,I)=NEA(MOX,MOY,PB,PA,I)+T1*T2
!if ((i<=9).and.(mox==2).and.(moy==3)) then
!t4=t4+t1*t2
!write(*,'(a,2i3,2f20.15)') 'EA pa,pb,t1*t2 =',pa,pb,t1*t2,nea(mox,moy,pb,pa,i)
!endif
           ENDDO
          ENDDO
!if ((i<=9).and.(mox==2).and.(moy==3)) then
!write(*,'(a,f20.15)') 'EA sum t1*t2 =',t4
!t5=t5+t4
!endif

         ENDDO
        ENDDO

        DO PA=1,NCF
         DO PB=1,IP_NCF
          T1=0.0D0
          DO QA=1,NCF
           DO QB=1,NCF
            T1=T1+WFN(QB,QA,N1)*WFN(QB,QA,N4)
           ENDDO
          ENDDO
          T2=0.0D0
          DO QA=1,NCF
           DO QB=1,IP_NCF
            T2=T2+IP_WFN(QB,QA,PB,PA,N2)*IP_WFN(QB,QA,PB,PA,N3)
           ENDDO
          ENDDO
          DIP(PB,PA,I)=DIP(PB,PA,I)+T1*T2
         ENDDO
        ENDDO

        DO PA=1,NCF
         DO PB=1,EA_NCF
          T1=0.0D0
          DO QA=1,NCF
           DO QB=1,NCF
            T1=T1+WFN(QB,QA,N1)*WFN(QB,QA,N4)
           ENDDO
          ENDDO
          T2=0.0D0
          DO QA=1,NCF
           DO QB=1,EA_NCF
            T2=T2+EA_WFN(QB,QA,PB,PA,N2)*EA_WFN(QB,QA,PB,PA,N3)
           ENDDO
          ENDDO
          DEA(PB,PA,I)=DEA(PB,PA,I)+T1*T2
         ENDDO
        ENDDO

!if (i==0) then
!write(*,*) '++++++++ IP'
!do pa=1,ncf
!do pb=1,ip_ncf
! write(*,*) 'pa,pb=',pa,pb
! call dump5(nip(:,:,pb,pa,0),iall(0,0,0)-ivirtcore)
!enddo
!enddo
!write(*,*) '++++++++ EA'
!do pa=1,ncf
!do pb=1,ea_ncf
! write(*,*) 'pa,pb=',pa,pb
! call dump5(nea(:,:,pb,pa,0),iall(0,0,0)-ivirtcore)
!enddo
!enddo
!endif

!write(*,*) 'i=',i
!call dump16(dip(:,:,i),ip_ncf,ncf)
!call dump16(dea(:,:,i),ea_ncf,ncf)

       ENDDO
      ENDDO
     ENDDO
    ENDDO
!write(*,'(a,f20.15)') 'sum all t1*t2 =',t5

   ENDDO
!  CLOSE(50)
!  CLOSE(51)

   ! CONSTRUCT GIP AND GEA

   ALLOCATE(GIP(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,IP_NCF,NCF,0:ORDER-1))
   ALLOCATE(GEA(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,EA_NCF,NCF,0:ORDER-1))

   GIP=1.0D99
   GEA=1.0D99

!goto 3
   DO I=0,ORDER-1

    GIP(:,:,:,:,I)=NIP(:,:,:,:,I)
    GEA(:,:,:,:,I)=NEA(:,:,:,:,I)
    
    IF (I > 0) THEN
     DO J=0,I-1
      K=I-J
      DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
       DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
 
        ! IP SECTOR
        DO PA=1,NCF
         DO PB=1,IP_NCF
!if ((j==2).and.(k==2).and.(norder(pa)+ip_norder(pb)==1)) cycle
          GIP(MOX,MOY,PB,PA,I)=GIP(MOX,MOY,PB,PA,I) &
           +GIP(MOX,MOY,PB,PA,J)*(IP_E(PB,PA,K)/(OMEGA-IP_E(PB,PA,0))-DIP(PB,PA,K))
!          +GIP(MOX,MOY,PB,PA,J)*(IP_E(PB,PA,K)/(OMEGA-IP_E(PB,PA,0)))
!          +GIP(MOX,MOY,PB,PA,J)*(-DIP(PB,PA,K))
         ENDDO
        ENDDO
 
        ! EA SECTOR
        DO PA=1,NCF
         DO PB=1,EA_NCF
!if ((j==2).and.(k==2).and.(norder(pa)+ea_norder(pb)==1)) cycle
          GEA(MOX,MOY,PB,PA,I)=GEA(MOX,MOY,PB,PA,I) &
           +GEA(MOX,MOY,PB,PA,J)*(EA_E(PB,PA,K)/(OMEGA-EA_E(PB,PA,0))-DEA(PB,PA,K))
!          +GEA(MOX,MOY,PB,PA,J)*(EA_E(PB,PA,K)/(OMEGA-EA_E(PB,PA,0)))
!          +GEA(MOX,MOY,PB,PA,J)*(-DEA(PB,PA,K))
         ENDDO
        ENDDO
 
       ENDDO
      ENDDO
     ENDDO
    ENDIF

   ENDDO

! --- manual 2- and 3-order GF without denom
3 continue
goto 4
write(*,*) 'DEBUG CODE IN OPERATION ~~~~~~~~~~~'
DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
DO PA=1,NCF
DO PB=1,IP_NCF
 GIP(MOX,MOY,PB,PA,2)= &
 +NIP(MOX,MOY,PB,PA,2) &
!-NIP(MOX,MOY,PB,PA,1)*DIP(PB,PA,1) &
!-NIP(MOX,MOY,PB,PA,0)*DIP(PB,PA,2) &
!+NIP(MOX,MOY,PB,PA,0)*DIP(PB,PA,1)*DIP(PB,PA,1) &
+0.0D0
ENDDO
ENDDO
DO PA=1,NCF
DO PB=1,EA_NCF
 GEA(MOX,MOY,PB,PA,2)= &
 +NEA(MOX,MOY,PB,PA,2) &
!-NEA(MOX,MOY,PB,PA,1)*DEA(PB,PA,1) &
!-NEA(MOX,MOY,PB,PA,0)*DEA(PB,PA,2) &
!+NEA(MOX,MOY,PB,PA,0)*DEA(PB,PA,1)*DEA(PB,PA,1) &
+0.0D0
ENDDO
ENDDO
ENDDO
ENDDO

DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
DO PA=1,NCF
DO PB=1,IP_NCF
 GIP(MOX,MOY,PB,PA,3)= &
 +NIP(MOX,MOY,PB,PA,3) &
!-NIP(MOX,MOY,PB,PA,2)*DIP(PB,PA,1) &
!-NIP(MOX,MOY,PB,PA,1)*DIP(PB,PA,2) &
!-NIP(MOX,MOY,PB,PA,0)*DIP(PB,PA,3) &
!+NIP(MOX,MOY,PB,PA,0)*DIP(PB,PA,1)*DIP(PB,PA,2) &
!+NIP(MOX,MOY,PB,PA,0)*DIP(PB,PA,2)*DIP(PB,PA,1) &
!-NIP(MOX,MOY,PB,PA,0)*DIP(PB,PA,1)*DIP(PB,PA,1)*DIP(PB,PA,1) &
+0.0D0
ENDDO
ENDDO
DO PA=1,NCF
DO PB=1,EA_NCF
 GEA(MOX,MOY,PB,PA,3)= &
 +NEA(MOX,MOY,PB,PA,3) &
!-NEA(MOX,MOY,PB,PA,2)*DEA(PB,PA,1) &
!-NEA(MOX,MOY,PB,PA,1)*DEA(PB,PA,2) &
!-NEA(MOX,MOY,PB,PA,0)*DEA(PB,PA,3) &
!+NEA(MOX,MOY,PB,PA,0)*DEA(PB,PA,1)*DEA(PB,PA,2) &
!+NEA(MOX,MOY,PB,PA,0)*DEA(PB,PA,2)*DEA(PB,PA,1) &
!-NEA(MOX,MOY,PB,PA,0)*DEA(PB,PA,1)*DEA(PB,PA,1)*DEA(PB,PA,1) &
+0.0D0
ENDDO
ENDDO
ENDDO
ENDDO
4 continue
! --- manual 2- and 3-order GF without denom end

   ! CONSTRUCT GN
   ALLOCATE(GN(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   GN=0.0D0

   DO I=0,ORDER-1

    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
     DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE

      ! IP SECTOR
!write(*,*) 'EA sector only!'
      DO PA=1,NCF
       DO PB=1,IP_NCF
        GN(MOX,MOY,I)=GN(MOX,MOY,I)+GIP(MOX,MOY,PB,PA,I)
       ENDDO
      ENDDO

      ! EA SECTOR
!write(*,*) 'IP sector only!'
      DO PA=1,NCF
       DO PB=1,EA_NCF
        GN(MOX,MOY,I)=GN(MOX,MOY,I)+GEA(MOX,MOY,PB,PA,I)
       ENDDO
      ENDDO

     ENDDO
    ENDDO

    WRITE(6,'(I3,A)') I,"TH-ORDER GREEN'S FUNCTION"
    CALL DUMP5(GN(:,:,I),IALL(0,0,0)-IVIRTCORE)

   ENDDO
   DEALLOCATE(NIP,NEA,DIP,DEA)
   DEALLOCATE(GIP,GEA)

   ALLOCATE(SN(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,ORDER))
   SN=0.0D0

   DO I=1,ORDER-1

!   WRITE(6,'(A,I3)') "SELF-ENERGY ORDER=",I

    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
     DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
      SN(MOX,MOY,I)=(OMEGA-EPSILON(MOX,0,0,0))*GN(MOX,MOY,I)*(OMEGA-EPSILON(MOY,0,0,0))
     ENDDO
    ENDDO

    IF (I > 1) THEN
     DO J=1,I-1
      DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
       DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
        DO MOZ=ICORE+1,IALL(0,0,0)-IVIRTCORE
         SN(MOX,MOY,I)=SN(MOX,MOY,I)-SN(MOX,MOZ,J)*GN(MOZ,MOY,I-J)*(OMEGA-EPSILON(MOY,0,0,0))
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDIF

!   WRITE(6,'(I3,A)') I,"TH-ORDER SELF ENERGY"
!   CALL DUMP5(SN(:,:,I),IALL(0,0,0)-IVIRTCORE)

   ENDDO

   DEALLOCATE(EIP,EEA,VIP,VEA)
   DEALLOCATE(IPVEC1,IPVEC2)
   DEALLOCATE(EAVEC1,EAVEC2)
   DEALLOCATE(GN,SN)

   RETURN

END SUBROUTINE



SUBROUTINE HIGHORDER_GF_REDUX1(OMEGA,ORDER,DIAG,EIGN,GREN,GDYS,NORB)
! PERFORM HIGH-ORDER MBGF CALCULATIONS WITHOUT IP/EA PERTURBATION EXPANSION

   USE CONTROL
   USE GRADIENT
   USE STRUCTURE
   USE INTEGRAL
   USE BASISSET
   USE FULLCI

   IMPLICIT NONE
   DOUBLE PRECISION,PARAMETER :: TOL = 1.0D-10
   INTEGER,PARAMETER :: MAXFILE = 100
   DOUBLE PRECISION :: OMEGA
   INTEGER :: ORDER
   INTEGER :: NORB
   DOUBLE PRECISION,ALLOCATABLE :: GN(:,:,:),GN2(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: SN1(:,:,:),SN2(:,:,:),SN3(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: GD(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: GIP(:,:,:),GEA(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: SIP(:,:,:),SEA(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: D(:)
   DOUBLE PRECISION,ALLOCATABLE :: VIP(:,:),VEA(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: EIP(:),EEA(:),E0(:)
   DOUBLE PRECISION,ALLOCATABLE :: WFN(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: XIP(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: XEA(:,:,:)
   INTEGER :: I,J,K,L,M,N
   INTEGER :: PA,PB,QA,QB
   INTEGER :: XA,XB
   INTEGER :: MOX,MOY,MOZ,MOW
   INTEGER :: SGN
   DOUBLE PRECISION :: T1
   DOUBLE PRECISION,ALLOCATABLE :: VEC1(:,:),VEC2(:,:)
   LOGICAL,ALLOCATABLE :: LKOOPIP(:),LKOOPEA(:)
   INTEGER,ALLOCATABLE :: IPMAP(:),EAMAP(:)
   INTEGER :: INFO
   DOUBLE PRECISION,ALLOCATABLE :: AMAT(:,:),VL(:,:),VR(:,:),ER(:),EI(:),WK(:)
   DOUBLE PRECISION,ALLOCATABLE :: B(:,:),BS(:,:),C(:),CS(:)
   INTEGER,ALLOCATABLE :: INDX(:)
   DOUBLE PRECISION :: X
   DOUBLE PRECISION :: DIAG(NORB,0:ORDER-1),EIGN(NORB,0:ORDER-1),GREN(NORB,0:ORDER-1),GDYS(NORB,0:ORDER-1)
!double precision,allocatable :: largevip(:,:,:),largevea(:,:,:)

   WRITE(6,'(/,A,F20.15,A)') "GENERAL-ORDER MANY-BODY GREEN'S FUNCTION METHOD USING SINGLE-DETERMINANT N+1/N-1 WAVE FUNCTIONS"
   WRITE(6,'(A,F20.15,A)') 'SELF-ENERGY MATRIX WILL BE COMPUTED AT OMEGA=',OMEGA,' HARTREE'
   OPEN(50,FILE=TRIM(COPTN(1))//'.fi0',FORM='UNFORMATTED')
   OPEN(51,FILE=TRIM(COPTN(1))//'.fo0',FORM='UNFORMATTED')

   ! READ ALL N-ELECTRON WAVE FUNCTIONS
   ALLOCATE(WFN(NCF,NCF,0:ORDER-1))
   OPEN(97,FILE=TRIM(COPTN(1))//'.gf0',FORM='UNFORMATTED')
   REWIND(97)
   DO I=0,ORDER-1
    READ(97) WFN(:,:,I)
   ENDDO
   CLOSE(97)
   WRITE(6,'(A,I3,A,A)') '0 THRU',ORDER-1,'TH-ORDER N-ELECTRON WAVE FUNCTIONS READ FROM ',TRIM(COPTN(1))//'.gf0'
   ! INTERMEDIATE NORMALIZATION SEE EQ(2.28) OF SHAVITT & BARTLETT
   DO I=1,ORDER-1
    T1=0.0D0
    DO PA=1,NCF
     DO PB=1,NCF
      T1=T1+WFN(PB,PA,I)*WFN(PB,PA,0)
     ENDDO
    ENDDO
    IF (DABS(T1) > TOL) THEN
     WRITE(6,*) I,T1
     CALL WARNING('NONORTHOGONAL N')
    ENDIF
   ENDDO
   ! TRIPLE CHECK
   ALLOCATE(VEC1(NCF,NCF))
   DO I=0,ORDER-1
    VEC1=0.0D0
    DO J=0,I
     VEC1(:,:)=VEC1(:,:)+WFN(:,:,J)
    ENDDO
    REWIND(50)
    WRITE(50) VEC1
    CALL HAMILTONIAN_PRODUCT(50,51,0,2*MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE))
    REWIND(51)
    READ(51) VEC1
    T1=0.0D0
    DO PA=1,NCF
     DO PB=1,NCF
      T1=T1+VEC1(PB,PA)*WFN(PB,PA,0)
     ENDDO
    ENDDO
    IF (DABS(T1+NUCLEAR_REPULSION-MP_STORED(I+1)) > TOL) THEN
     WRITE(6,'(I3,2F20.15)') I, T1+NUCLEAR_REPULSION,MP_STORED(I+1)
     CALL PABORT('ENERGY INCONSISTENCY N')
    ENDIF
   ENDDO
   DEALLOCATE(VEC1)

   ALLOCATE(E0(0:ORDER-1))
   E0(0)=MP_STORED(0)-NUCLEAR_REPULSION
   DO I=1,ORDER-1
    E0(I)=MP_STORED(I)-MP_STORED(I-1)
   ENDDO
!  WRITE(6,*) 'E0 VECTOR (ORDER IS SHIFTED BY 1)'
!  CALL DUMP16(E0,ORDER,1)

   ! FORM XIP,XEA
   ALLOCATE(XIP(NCF*IP_NCF,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(XEA(NCF*EA_NCF,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   XIP=1.0D99 ! safety
   XEA=1.0D99 ! safety
   DO I=0,ORDER-1

    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE

     ! IP SECTOR
     ALLOCATE(VEC1(IP_NCF,NCF))
     VEC1=0.0D0
     DO PA=1,NCF
      DO PB=1,NCF
       IF (.NOT.BTEST(CFHALF(PB),MOX-1)) CYCLE
       QA=PA
       QB=IP_ADDRSS(IBCLR(CFHALF(PB),MOX-1))
       SGN=1
       IF (MOX /= 1) THEN
        DO J=1,MOX-1
         IF (BTEST(CFHALF(PB),J-1)) SGN=-SGN
        ENDDO
       ENDIF
       VEC1(QB,QA)=VEC1(QB,QA)+WFN(PB,PA,I)*DFLOAT(SGN)
      ENDDO
     ENDDO
!    REWIND(50)
!    WRITE(50) WFN(:,:,I)
!    XA=1
!    IF (MOX <= IOCC) THEN
!     XB=IP_ADDRSS(IBCLR(CFHALF(1),MOX-1))
!    ELSE
!     XB=IP_NCF+EA_ADDRSS(IBSET(CFHALF(1),MOX-1))
!    ENDIF
!    CALL IONIZE(50,51,XA,XB,.FALSE.,0)
!    REWIND(51)
!    READ(51) VEC1
     DO PA=1,NCF
      DO PB=1,IP_NCF
       XIP((PA-1)*IP_NCF+PB,MOX,I)=VEC1(PB,PA)
      ENDDO
     ENDDO
     DEALLOCATE(VEC1)

     ! EA SECTOR
     ALLOCATE(VEC1(EA_NCF,NCF))
     VEC1=0.0D0
     DO PA=1,NCF
      DO PB=1,NCF
       IF (BTEST(CFHALF(PB),MOX-1)) CYCLE
       QA=PA
       QB=EA_ADDRSS(IBSET(CFHALF(PB),MOX-1))
       SGN=1
       IF (MOX /= 1) THEN
        DO J=1,MOX-1
         IF (BTEST(CFHALF(PB),J-1)) SGN=-SGN
        ENDDO
       ENDIF
       VEC1(QB,QA)=VEC1(QB,QA)+WFN(PB,PA,I)*DFLOAT(SGN)
      ENDDO
     ENDDO
!    REWIND(50)
!    WRITE(50) WFN(:,:,I)
!    XA=1
!    IF (MOX <= IOCC) THEN
!     XB=IP_ADDRSS(IBCLR(CFHALF(1),MOX-1))
!    ELSE
!     XB=IP_NCF+EA_ADDRSS(IBSET(CFHALF(1),MOX-1))
!    ENDIF
!    CALL IONIZE(50,51,XA,XB,.TRUE.,0)
!    REWIND(51)
!    READ(51) VEC1
     DO PA=1,NCF
      DO PB=1,EA_NCF
       XEA((PA-1)*EA_NCF+PB,MOX,I)=VEC1(PB,PA)
      ENDDO
     ENDDO
     DEALLOCATE(VEC1)

    ENDDO

!   WRITE(6,*) "XIP VECTOR AT ORDER=",I
!   CALL DUMP16(XIP(:,:,I),NCF*IP_NCF,IALL(0,0,0)-IVIRTCORE)
!   WRITE(6,*) "XEA VECTOR AT ORDER=",I
!   CALL DUMP16(XEA(:,:,I),NCF*EA_NCF,IALL(0,0,0)-IVIRTCORE)

   ENDDO

   ! FORM VIP,VEA
   ALLOCATE(EIP(NCF*IP_NCF))
   ALLOCATE(EEA(NCF*EA_NCF))
   ALLOCATE(VIP(NCF*IP_NCF,NCF*IP_NCF))
   ALLOCATE(VEA(NCF*EA_NCF,NCF*EA_NCF))

   ! IP SECTOR
   ALLOCATE(VEC1(IP_NCF,NCF),VEC2(IP_NCF,NCF))
   DO PA=1,NCF
    DO PB=1,IP_NCF
     VEC1=0.0D0
     VEC1(PB,PA)=1.0D0
     REWIND(50)
     WRITE(50) VEC1
     CALL IP_HAMILTONIAN_PRODUCT(50,51,0,MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE)+MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE+1))
     REWIND(51)
     READ(51) VEC1
     CALL H0_PRODUCT(50,51,1)
     REWIND(51)
     READ(51) VEC2
     EIP((PA-1)*IP_NCF+PB)=VEC2(PB,PA)
     VEC1=VEC1-VEC2 ! V=H-H0
     DO QA=1,NCF
      DO QB=1,IP_NCF
       VIP((QA-1)*IP_NCF+QB,(PA-1)*IP_NCF+PB)=VEC1(QB,QA)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
   DEALLOCATE(VEC1,VEC2)
!  WRITE(6,*) 'EIP VECTOR'
!  CALL DUMP16(EIP,NCF*IP_NCF,1)
!  WRITE(6,*) 'VIP MATRIX'
!  CALL DUMP5(VIP,NCF*IP_NCF)

   ! EA SECTOR
   ALLOCATE(VEC1(EA_NCF,NCF),VEC2(EA_NCF,NCF))
   DO PA=1,NCF
    DO PB=1,EA_NCF
     VEC1=0.0D0
     VEC1(PB,PA)=1.0D0
     REWIND(50)
     WRITE(50) VEC1
     CALL EA_HAMILTONIAN_PRODUCT(50,51,0,MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE)+MIN(IOCC-ICORE+1,IALL(0,0,0)-IOCC-IVIRTCORE))
     REWIND(51)
     READ(51) VEC1
     CALL H0_PRODUCT(50,51,-1)
     REWIND(51)
     READ(51) VEC2
     EEA((PA-1)*EA_NCF+PB)=VEC2(PB,PA)
     VEC1=VEC1-VEC2 ! V=H-H0
     DO QA=1,NCF
      DO QB=1,EA_NCF
       VEA((QA-1)*EA_NCF+QB,(PA-1)*EA_NCF+PB)=VEC1(QB,QA)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
   DEALLOCATE(VEC1,VEC2)
!  WRITE(6,*) 'EEA VECTOR'
!  CALL DUMP16(EEA,NCF*EA_NCF,1)
!  WRITE(6,*) 'VEA MATRIX'
!  CALL DUMP5(VEA,NCF*EA_NCF)

   ! FORM D
   ALLOCATE(D(0:ORDER-1))
   DO I=0,ORDER-1
    D(I)=0.0D0
    DO J=0,I
     DO PA=1,NCF
      DO PB=1,NCF
       D(I)=D(I)+WFN(PB,PA,J)*WFN(PB,PA,I-J)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
!  WRITE(6,*) 'D VECTOR'
!  CALL DUMP16(D,ORDER,1)

   ! FORM LKOOPIP AND LKOOPEA
   ALLOCATE(LKOOPIP(NCF*IP_NCF))
   ALLOCATE(LKOOPEA(NCF*EA_NCF))
   LKOOPIP=.FALSE.
   DO PA=1,NCF
    DO PB=1,IP_NCF
     IF (NORDER(PA)+IP_NORDER(PB) < 2) LKOOPIP((PA-1)*IP_NCF+PB)=.TRUE.
    ENDDO
   ENDDO
   LKOOPEA=.FALSE.
   DO PA=1,NCF
    DO PB=1,EA_NCF
     IF (NORDER(PA)+EA_NORDER(PB) < 2) LKOOPEA((PA-1)*EA_NCF+PB)=.TRUE.
    ENDDO
   ENDDO
   ! CHECK
   J=0
   DO PA=1,NCF*IP_NCF
    IF (LKOOPIP(PA)) J=J+1
   ENDDO
!  WRITE(*,*) J
   J=0
   DO PA=1,NCF*EA_NCF
    IF (LKOOPEA(PA)) J=J+1
   ENDDO
!  WRITE(*,*) J

   ! MAP ORBITALS TO DETERMINANTS
   ALLOCATE(IPMAP(IALL(0,0,0)-IVIRTCORE),EAMAP(IALL(0,0,0)-IVIRTCORE))
   IPMAP=99999999
   EAMAP=99999999
!  WRITE(6,'(A)') 'ORBITAL -> KOOPMANS DETERMINANT'
   DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
    IF (MOX <= IOCC) THEN
     DO PA=1,IP_NCF
      IF ((IP_NORDER(PA)==1).AND.(.NOT.BTEST(IP_CFHALF(PA),MOX-1))) IPMAP(MOX)=PA
     ENDDO
    ELSE
     DO PA=1,EA_NCF
      IF ((EA_NORDER(PA)==1).AND.(BTEST(EA_CFHALF(PA),MOX-1))) EAMAP(MOX)=PA
     ENDDO
    ENDIF
!   WRITE(6,'(3I5)') MOX,IPMAP(MOX),EAMAP(MOX)
   ENDDO

   ALLOCATE(GN(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(GN2(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(SN1(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(SN2(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
!  ALLOCATE(SN3(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(GD(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(GIP(NCF*IP_NCF,NCF*IP_NCF,0:ORDER-1))
   ALLOCATE(GEA(NCF*EA_NCF,NCF*EA_NCF,0:ORDER-1))
!  ALLOCATE(SIP(NCF*IP_NCF,NCF*IP_NCF,0:ORDER-1))
!  ALLOCATE(SEA(NCF*EA_NCF,NCF*EA_NCF,0:ORDER-1))
   GIP=1.0D99 ! safety
   GEA=1.0D99 ! safety
!  SIP=1.0D99 ! safety
!  SEA=1.0D99 ! safety
   SN1=1.0D99
   SN2=1.0D99
   GN=1.0D99

   ! GRAND LOOP
   DO I=0,ORDER-1
   
    ! FORM GIP AND GEA 
    IF (I==0) THEN

     ! IP SECTOR
     GIP(:,:,0)=0.0D0
     DO PA=1,NCF*IP_NCF
!     GIP(PA,PA,0)=0.0D0
      GIP(PA,PA,0)=1.0D0/(OMEGA-E0(0)+EIP(PA))
     ENDDO
     ! EA SECTOR
     GEA(:,:,0)=0.0D0
     DO PA=1,NCF*EA_NCF
!     GEA(PA,PA,0)=0.0D0
      GEA(PA,PA,0)=1.0D0/(OMEGA-EEA(PA)+E0(0))
     ENDDO
!    WRITE(6,*) 'ZEROTH-ORDER GIP'
!    CALL DUMP5(GIP(:,:,0),NCF*IP_NCF)
!    WRITE(6,*) 'ZEROTH-ORDER GEA'
!    CALL DUMP5(GEA(:,:,0),NCF*EA_NCF)

    ELSE

     ! IP SECTOR
     DO N=1,NCF*IP_NCF
      DO M=1,NCF*IP_NCF
       GIP(N,M,I)=0.0D0
       DO L=1,NCF*IP_NCF
        GIP(N,M,I)=GIP(N,M,I)-GIP(N,L,I-1)*VIP(L,M)*GIP(M,M,0)
       ENDDO
       DO J=1,I
        GIP(N,M,I)=GIP(N,M,I)+GIP(N,M,I-J)*E0(J)*GIP(M,M,0)
!       GIP(N,M,I)=GIP(N,M,I)-D(J)*GIP(N,M,I-J) This gives wrong results 4th and higher orders
!                                               because D is the derivative of x & y vectors
!                                               which are not even specified at this stage
       ENDDO
      ENDDO
     ENDDO

     ! EA SECTOR
     DO N=1,NCF*EA_NCF
      DO M=1,NCF*EA_NCF
       GEA(N,M,I)=0.0D0
       DO L=1,NCF*EA_NCF
        GEA(N,M,I)=GEA(N,M,I)+GEA(N,L,I-1)*VEA(L,M)*GEA(M,M,0)
       ENDDO
       DO J=1,I
        GEA(N,M,I)=GEA(N,M,I)-GEA(N,M,I-J)*E0(J)*GEA(M,M,0)
!       GEA(N,M,I)=GEA(N,M,I)-D(J)*GEA(N,M,I-J) This gives wrong results 4th and higher orders
!                                               because D is the derivative of x & y vectors
!                                               which are not even specified at this stage
       ENDDO
      ENDDO
     ENDDO

    ENDIF

    ! FORM GN
    GN(:,:,I)=0.0D0
    GN2(:,:,I)=0.0D0
    DO J=0,I
     DO K=0,I
      DO L=0,I
       IF (J+K+L==I) THEN
        DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
         DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
          ! IP SECTOR
          DO PA=1,NCF*IP_NCF
           DO PB=1,NCF*IP_NCF
            GN(MOX,MOY,I)=GN(MOX,MOY,I)+XIP(PA,MOX,J)*GIP(PA,PB,K)*XIP(PB,MOY,L) 
           ENDDO
          ENDDO
          ! EA SECTOR
          DO PA=1,NCF*EA_NCF
           DO PB=1,NCF*EA_NCF
            GN(MOX,MOY,I)=GN(MOX,MOY,I)+XEA(PA,MOY,J)*GEA(PA,PB,K)*XEA(PB,MOX,L)
           ENDDO
          ENDDO
         ENDDO 
        ENDDO
       ENDIF
      ENDDO
     ENDDO
    ENDDO
    IF (I > 0) THEN
     DO J=1,I
      DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
       DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
        GN(MOX,MOY,I)=GN(MOX,MOY,I)-D(J)*GN(MOX,MOY,I-J)
       ENDDO
      ENDDO
     ENDDO
    ENDIF
    DO J=0,I
     DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
      DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
       GN2(MOX,MOY,I)=GN2(MOX,MOY,I)+GN(MOX,MOY,J)
      ENDDO
     ENDDO
    ENDDO

    ! FORM SN WITH NUMERICAL REDUCIBLE DELETION
    IF (I==0) THEN
     DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
      DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
       SN2(MOX,MOY,I)=0.0D0
      ENDDO
      SN2(MOX,MOX,I)=EPSILON(MOX,0,0,0)
     ENDDO
!    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
!     DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
!      SN1(MOX,MOY,I)=0.0D0
!     ENDDO
!     SN1(MOX,MOX,I)=OMEGA-EPSILON(MOX,0,0,0)
!    ENDDO
    ELSE IF (I==1) THEN
     DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
      DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
!      SN1(MOX,MOY,I)=0.0d0
       SN1(MOX,MOY,I)=GN(MOX,MOY,1)/GN(MOX,MOX,0)/GN(MOY,MOY,0)
       SN2(MOX,MOY,I)=SN2(MOX,MOY,I-1)+SN1(MOX,MOY,I)
      ENDDO
     ENDDO
    ELSE
     DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
      DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
       SN1(MOX,MOY,I)=GN(MOX,MOY,I)
       DO J=1,I-1
        DO MOZ=ICORE+1,IALL(0,0,0)-IVIRTCORE
         DO MOW=ICORE+1,IALL(0,0,0)-IVIRTCORE
! Note: the renormalization term below makes difference only at fourth and higher orders
          SN1(MOX,MOY,I)=SN1(MOX,MOY,I)-GN(MOX,MOZ,0)*SN1(MOZ,MOW,J)*GN(MOW,MOY,I-J)
         ENDDO
        ENDDO
       ENDDO
       SN1(MOX,MOY,I)=SN1(MOX,MOY,I)/GN(MOX,MOX,0)/GN(MOY,MOY,0)
       SN2(MOX,MOY,I)=SN2(MOX,MOY,I-1)+SN1(MOX,MOY,I)
      ENDDO
     ENDDO
    ENDIF
!   IF (I > 0) THEN
!    WRITE(6,'(I3,A)') I,"TH-ORDER SELF-ENERGY WITH NUMERICAL REDUCIBLE CANCELLATION"
!    CALL DUMP5(SN(:,:,I),IALL(0,0,0)-IVIRTCORE)
!   ENDIF

!   ! FORM SIP AND SEA 
!   IF (I==0) THEN

!    SIP(:,:,0)=1.0D99
!    SEA(:,:,0)=1.0D99
!    DO PA=1,NCF*IP_NCF
!     SIP(PA,PA,0)=(OMEGA-E0(0)+EIP(PA))
!    ENDDO
!    DO PA=1,NCF*EA_NCF
!     SEA(PA,PA,0)=(OMEGA-EEA(PA)+E0(0))
!    ENDDO

!   ELSE IF (I==1) THEN

!    ! IP SECTOR
!    DO PA=1,NCF*IP_NCF
!     DO PB=1,NCF*IP_NCF
!      SIP(PA,PB,1)=-VIP(PA,PB)
!     ENDDO
!     SIP(PA,PA,1)=SIP(PA,PA,1)+E0(1)
!    ENDDO
!    ! EA SECTOR
!    DO PA=1,NCF*EA_NCF
!     DO PB=1,NCF*EA_NCF
!      SEA(PA,PB,1)=VEA(PA,PB)
!     ENDDO
!     SEA(PA,PA,1)=SEA(PA,PA,1)-E0(1)
!    ENDDO

!   ELSE

!    ! IP SECTOR
!    DO N=1,NCF*IP_NCF
!     DO M=1,NCF*IP_NCF
!      SIP(N,M,I)=0.0D0
!      DO L=1,NCF*IP_NCF
!       IF (.NOT.LKOOPIP(L)) &
!if (dabs(gip(l,l,0)) < 1.0d6) &
!        SIP(N,M,I)=SIP(N,M,I)-SIP(N,L,I-1)*GIP(L,L,0)*VIP(L,M)
!      ENDDO
!      DO J=1,I-1
!       IF (.NOT.LKOOPIP(M)) &
!if (dabs(gip(m,m,0)) < 1.0d6) &
!        SIP(N,M,I)=SIP(N,M,I)+SIP(N,M,I-J)*GIP(M,M,0)*E0(J)
!      ENDDO
!     ENDDO
!     SIP(N,N,I)=SIP(N,N,I)+E0(I)
!    ENDDO

!    ! EA SECTOR
!    DO N=1,NCF*EA_NCF
!     DO M=1,NCF*EA_NCF
!      SEA(N,M,I)=0.0D0
!      DO L=1,NCF*EA_NCF
!       IF (.NOT.LKOOPEA(L)) &
!if (dabs(gea(l,l,0)) < 1.0d6) &
!        SEA(N,M,I)=SEA(N,M,I)+SEA(N,L,I-1)*GEA(L,L,0)*VEA(L,M)
!      ENDDO
!      DO J=1,I-1
!       IF (.NOT.LKOOPEA(M)) &
!if (dabs(gea(m,m,0)) < 1.0d6) &
!        SEA(N,M,I)=SEA(N,M,I)-SEA(N,M,I-J)*GEA(M,M,0)*E0(J)
!      ENDDO
!     ENDDO
!     SEA(N,N,I)=SEA(N,N,I)-E0(I)
!    ENDDO

!   ENDIF

    ! FORM SN
!   SN2(:,:,I)=0.0D0
!   SN3(:,:,I)=0.0D0
!   DO J=0,I
!    DO L=0,I
!     DO K=1,I
!      IF (J+K+L==I) THEN
!       DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
!        DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
!         ! IP SECTOR
!         DO PA=1,NCF*IP_NCF
!          DO PB=1,NCF*IP_NCF
!           T1=XIP(PA,MOX,J)*SIP(PA,PB,K)*XIP(PB,MOY,L)
!           IF (IPMAP(MOX)/=PA) &
!            T1=T1*GIP(PA,PA,0)/GN(MOX,MOX,0)
!           IF (IPMAP(MOY)/=PB) &
!            T1=T1*GIP(PB,PB,0)/GN(MOY,MOY,0)
!           SN2(MOX,MOY,I)=SN2(MOX,MOY,I)+T1
!           IF ((.NOT.LKOOPIP(PA)).OR.(.NOT.LKOOPIP(PB))) SN3(MOX,MOY,I)=SN3(MOX,MOY,I)+T1
!          ENDDO
!         ENDDO
!         ! EA SECTOR
!         DO PA=1,NCF*EA_NCF
!          DO PB=1,NCF*EA_NCF
!           T1=XEA(PA,MOY,J)*SEA(PA,PB,K)*XEA(PB,MOX,L)
!           IF (EAMAP(MOX)/=PA) &
!            T1=T1*GEA(PA,PA,0)/GN(MOX,MOX,0)
!           IF (EAMAP(MOY)/=PB) &
!            T1=T1*GEA(PB,PB,0)/GN(MOY,MOY,0)
!           SN2(MOX,MOY,I)=SN2(MOX,MOY,I)+T1
!           IF ((.NOT.LKOOPEA(PA)).OR.(.NOT.LKOOPEA(PB))) SN3(MOX,MOY,I)=SN3(MOX,MOY,I)+T1
!          ENDDO
!         ENDDO
!        ENDDO 
!       ENDDO
!      ENDIF
!     ENDDO
!     IF (J+L==I) THEN
!      DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
!       DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
!        ! IP SECTOR
!        DO PA=1,NCF*IP_NCF
!         T1=XIP(PA,MOX,J)/GIP(PA,PA,0)*XIP(PA,MOY,L)
!         IF (IPMAP(MOX)/=PA) &
!          T1=T1*GIP(PA,PA,0)/GN(MOX,MOX,0)
!         IF (IPMAP(MOY)/=PA) &
!          T1=T1*GIP(PA,PA,0)/GN(MOY,MOY,0)
!         SN2(MOX,MOY,I)=SN2(MOX,MOY,I)+T1
!         IF (.NOT.LKOOPIP(PA)) SN3(MOX,MOY,I)=SN3(MOX,MOY,I)+T1
!        ENDDO
!        ! EA SECTOR
!        DO PA=1,NCF*EA_NCF
!         T1=XEA(PA,MOY,J)/GEA(PA,PA,0)*XEA(PA,MOX,L)
!         IF (EAMAP(MOX)/=PA) &
!          T1=T1*GEA(PA,PA,0)/GN(MOX,MOX,0)
!         IF (EAMAP(MOY)/=PA) &
!          T1=T1*GEA(PA,PA,0)/GN(MOY,MOY,0)
!         SN2(MOX,MOY,I)=SN2(MOX,MOY,I)+T1
!         IF (.NOT.LKOOPEA(PA)) SN3(MOX,MOY,I)=SN3(MOX,MOY,I)+T1
!        ENDDO
!       ENDDO
!      ENDDO
!     ENDIF
!    ENDDO
!   ENDDO
!   DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
!    DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
!     SN2(MOX,MOY,I)=SN2(MOX,MOY,I)/GN(MOX,MOX,0)/GN(MOY,MOY,0)
!    ENDDO
!   ENDDO
!   IF (I==0) THEN
!    CYCLE
!   ELSE 
!    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
!     SN2(MOX,MOX,I)=SN2(MOX,MOX,I)-D(I)*(OMEGA-EPSILON(MOX,0,0,0))
!     SN3(MOX,MOX,I)=SN3(MOX,MOX,I)-D(I)*(OMEGA-EPSILON(MOX,0,0,0))
!    ENDDO
!    IF (I>=2) THEN
!     DO J=1,I-1
!      DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
!       DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
!        SN2(MOX,MOY,I)=SN2(MOX,MOY,I)-D(J)*SN2(MOX,MOY,I-J)
!        SN3(MOX,MOY,I)=SN3(MOX,MOY,I)-D(J)*SN3(MOX,MOY,I-J)
!       ENDDO
!      ENDDO
!     ENDDO
!    ENDIF
!   ENDIF
!   IF (I>=2) THEN
!    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
!     DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
!      DO J=1,I-1
!       DO MOZ=ICORE+1,IALL(0,0,0)-IVIRTCORE
!        DO MOW=ICORE+1,IALL(0,0,0)-IVIRTCORE
!         SN2(MOX,MOY,I)=SN2(MOX,MOY,I)-SN2(MOX,MOW,J)*GN(MOW,MOY,I-J)/GN(MOY,MOY,0)
!         SN2(MOX,MOY,I)=SN2(MOX,MOY,I)-SN3(MOX,MOW,J)*GN(MOW,MOW,0)*SN3(MOW,MOY,I-J)
!if ( ((mox > iocc).and.(mow <= iocc)).or.((mox <= iocc).and.(mow > iocc)).or. &
!    ((mow > iocc).and.(moy <= iocc)).or.((mow <= iocc).and.(moy > iocc)) ) &
!if (dabs(gn(mow,mow,0)) < 1.0d6) &
!         SN2(MOX,MOY,I)=SN2(MOX,MOY,I)-SN2(MOX,MOW,J)*GN(MOW,MOW,0)*SN2(MOW,MOY,I-J)
!        ENDDO
!       ENDDO
!      ENDDO
!     ENDDO
!    ENDDO
!   ENDIF


! brute force GN(2) 
goto 999
gn=1.0D99
GN(:,:,0)=0.0d0
GN(:,:,1)=0.0d0
GN(:,:,2)=0.0d0
do mox=icore+1,iall(0,0,0)-ivirtcore
do moy=icore+1,iall(0,0,0)-ivirtcore
! first 3 terms
do pa=1,ncf*ip_ncf
 gn(mox,moy,0)=gn(mox,moy,0)+xip(pa,mox,0)*gip(pa,pa,0)*xip(pa,moy,0)

 gn(mox,moy,1)=gn(mox,moy,1)+xip(pa,mox,0)*gip(pa,pa,0)*xip(pa,moy,1)
 gn(mox,moy,1)=gn(mox,moy,1)+xip(pa,mox,1)*gip(pa,pa,0)*xip(pa,moy,0)

 gn(mox,moy,2)=gn(mox,moy,2)+xip(pa,mox,0)*gip(pa,pa,0)*xip(pa,moy,2)
 gn(mox,moy,2)=gn(mox,moy,2)+xip(pa,mox,1)*gip(pa,pa,0)*xip(pa,moy,1)
 gn(mox,moy,2)=gn(mox,moy,2)+xip(pa,mox,2)*gip(pa,pa,0)*xip(pa,moy,0)
!gn(mox,moy,2)=gn(mox,moy,2)+xip(pa,mox,0)*xip(pa,moy,2)
!gn(mox,moy,2)=gn(mox,moy,2)+xip(pa,mox,1)*xip(pa,moy,1)
!gn(mox,moy,2)=gn(mox,moy,2)+xip(pa,mox,2)*xip(pa,moy,0)
enddo
do pa=1,ncf*ea_ncf
 gn(mox,moy,0)=gn(mox,moy,0)+xea(pa,moy,0)*gea(pa,pa,0)*xea(pa,mox,0)

 gn(mox,moy,1)=gn(mox,moy,1)+xea(pa,moy,0)*gea(pa,pa,0)*xea(pa,mox,1)
 gn(mox,moy,1)=gn(mox,moy,1)+xea(pa,moy,1)*gea(pa,pa,0)*xea(pa,mox,0)

 gn(mox,moy,2)=gn(mox,moy,2)+xea(pa,moy,0)*gea(pa,pa,0)*xea(pa,mox,2)
 gn(mox,moy,2)=gn(mox,moy,2)+xea(pa,moy,1)*gea(pa,pa,0)*xea(pa,mox,1)
 gn(mox,moy,2)=gn(mox,moy,2)+xea(pa,moy,2)*gea(pa,pa,0)*xea(pa,mox,0)
!gn(mox,moy,2)=gn(mox,moy,2)+xea(pa,moy,0)*xea(pa,mox,2)
!gn(mox,moy,2)=gn(mox,moy,2)+xea(pa,moy,1)*xea(pa,mox,1)
!gn(mox,moy,2)=gn(mox,moy,2)+xea(pa,moy,2)*xea(pa,mox,0)
enddo
!goto 111
! next 4 terms
do pa=1,ncf*ip_ncf
do pb=1,ncf*ip_ncf
gn(mox,moy,1)=gn(mox,moy,1)-xip(pa,mox,0)*gip(pa,pa,0)*vip(pa,pb)*gip(pb,pb,0)*xip(pb,moy,0)

gn(mox,moy,2)=gn(mox,moy,2)-xip(pa,mox,0)*gip(pa,pa,0)*vip(pa,pb)*gip(pb,pb,0)*xip(pb,moy,1)
gn(mox,moy,2)=gn(mox,moy,2)-xip(pa,mox,1)*gip(pa,pa,0)*vip(pa,pb)*gip(pb,pb,0)*xip(pb,moy,0)
enddo
gn(mox,moy,1)=gn(mox,moy,1)+xip(pa,mox,0)*gip(pa,pa,0)*e0(1)*gip(pa,pa,0)*xip(pa,moy,0)

gn(mox,moy,2)=gn(mox,moy,2)+xip(pa,mox,0)*gip(pa,pa,0)*e0(1)*gip(pa,pa,0)*xip(pa,moy,1)
gn(mox,moy,2)=gn(mox,moy,2)+xip(pa,mox,1)*gip(pa,pa,0)*e0(1)*gip(pa,pa,0)*xip(pa,moy,0)
enddo
do pa=1,ncf*ea_ncf
do pb=1,ncf*ea_ncf
gn(mox,moy,1)=gn(mox,moy,1)+xea(pa,moy,0)*gea(pa,pa,0)*vea(pa,pb)*gea(pb,pb,0)*xea(pb,mox,0)

gn(mox,moy,2)=gn(mox,moy,2)+xea(pa,moy,0)*gea(pa,pa,0)*vea(pa,pb)*gea(pb,pb,0)*xea(pb,mox,1)
gn(mox,moy,2)=gn(mox,moy,2)+xea(pa,moy,1)*gea(pa,pa,0)*vea(pa,pb)*gea(pb,pb,0)*xea(pb,mox,0)
enddo
gn(mox,moy,1)=gn(mox,moy,1)-xea(pa,moy,0)*gea(pa,pa,0)*e0(1)*gea(pa,pa,0)*xea(pa,mox,0)

gn(mox,moy,2)=gn(mox,moy,2)-xea(pa,moy,0)*gea(pa,pa,0)*e0(1)*gea(pa,pa,0)*xea(pa,mox,1)
gn(mox,moy,2)=gn(mox,moy,2)-xea(pa,moy,1)*gea(pa,pa,0)*e0(1)*gea(pa,pa,0)*xea(pa,mox,0)
enddo
!goto 112
! next 4 terms
do pa=1,ncf*ip_ncf
do pb=1,ncf*ip_ncf
do qa=1,ncf*ip_ncf
!   if (lkoopip(pa).or.lkoopip(pb).or.lkoopip(qa)) cycle
gn(mox,moy,2)=gn(mox,moy,2)+xip(pa,mox,0)*gip(pa,pa,0)*vip(pa,pb)*gip(pb,pb,0)*vip(pb,qa)*gip(qa,qa,0)*xip(qa,moy,0)
enddo
!   if (lkoopip(pa).or.lkoopip(pb)) cycle
gn(mox,moy,2)=gn(mox,moy,2)-xip(pa,mox,0)*gip(pa,pa,0)*vip(pa,pb)*gip(pb,pb,0)*e0(1)*gip(pb,pb,0)*xip(pb,moy,0)
!   if (lkoopip(pa).or.lkoopip(pb)) cycle
gn(mox,moy,2)=gn(mox,moy,2)-xip(pa,mox,0)*gip(pa,pa,0)*e0(1)*gip(pa,pa,0)*vip(pa,pb)*gip(pb,pb,0)*xip(pb,moy,0)
enddo
!   if (lkoopip(pa)) cycle
gn(mox,moy,2)=gn(mox,moy,2)+xip(pa,mox,0)*gip(pa,pa,0)*e0(1)*gip(pa,pa,0)*e0(1)*gip(pa,pa,0)*xip(pa,moy,0)
enddo
do pa=1,ncf*ea_ncf
do pb=1,ncf*ea_ncf
do qa=1,ncf*ea_ncf
!   if (lkoopea(pa).or.lkoopea(pb).or.lkoopea(qa)) cycle
gn(mox,moy,2)=gn(mox,moy,2)+xea(pa,moy,0)*gea(pa,pa,0)*vea(pa,pb)*gea(pb,pb,0)*vea(pb,qa)*gea(qa,qa,0)*xea(qa,mox,0)
enddo
!   if (lkoopea(pa).or.lkoopea(pb)) cycle
gn(mox,moy,2)=gn(mox,moy,2)-xea(pa,moy,0)*gea(pa,pa,0)*vea(pa,pb)*gea(pb,pb,0)*e0(1)*gea(pb,pb,0)*xea(pb,mox,0)
!   if (lkoopea(pa).or.lkoopea(pb)) cycle
gn(mox,moy,2)=gn(mox,moy,2)-xea(pa,moy,0)*gea(pa,pa,0)*e0(1)*gea(pa,pa,0)*vea(pa,pb)*gea(pb,pb,0)*xea(pb,mox,0)
enddo
!   if (lkoopea(pa)) cycle
gn(mox,moy,2)=gn(mox,moy,2)+xea(pa,moy,0)*gea(pa,pa,0)*e0(1)*gea(pa,pa,0)*e0(1)*gea(pa,pa,0)*xea(pa,mox,0)
enddo
112 continue
! next term
do pa=1,ncf*ip_ncf
gn(mox,moy,2)=gn(mox,moy,2)+xip(pa,mox,0)*gip(pa,pa,0)*e0(2)*gip(pa,pa,0)*xip(pa,moy,0)
enddo
do pa=1,ncf*ea_ncf
gn(mox,moy,2)=gn(mox,moy,2)-xea(pa,moy,0)*gea(pa,pa,0)*e0(2)*gea(pa,pa,0)*xea(pa,mox,0)
enddo
! D's
111 continue
do pa=1,ncf*ip_ncf
gn(mox,moy,1)=gn(mox,moy,1)-d(1)*xip(pa,mox,0)*gip(pa,pa,0)*xip(pa,moy,0)

gn(mox,moy,2)=gn(mox,moy,2)-(d(2)-d(1)*d(1))*xip(pa,mox,0)*gip(pa,pa,0)*xip(pa,moy,0)
!gn(mox,moy,2)=gn(mox,moy,2)-(d(2)-d(1)*d(1))*xip(pa,mox,0)*xip(pa,moy,0)
enddo
do pa=1,ncf*ea_ncf
gn(mox,moy,1)=gn(mox,moy,1)-d(1)*xea(pa,moy,0)*gea(pa,pa,0)*xea(pa,mox,0)

gn(mox,moy,2)=gn(mox,moy,2)-(d(2)-d(1)*d(1))*xea(pa,moy,0)*gea(pa,pa,0)*xea(pa,mox,0)
!gn(mox,moy,2)=gn(mox,moy,2)-(d(2)-d(1)*d(1))*xea(pa,moy,0)*xea(pa,mox,0)
enddo

enddo
enddo
999 continue

   ENDDO

   WRITE(6,'(/,A,F20.15,A)') "METHOD 1: SINGLE-DETERMINANT N+1/N-1 WAVE FUNCTIONS"
   DO I=0,ORDER-1
    WRITE(6,'(I3,A)') I,"TH-ORDER GREEN'S FUNCTION"
    CALL DUMP5(GN(:,:,I),IALL(0,0,0)-IVIRTCORE)
   ENDDO

   WRITE(6,'(/,A,F20.15,A)') "METHOD 1: SINGLE-DETERMINANT N+1/N-1 WAVE FUNCTIONS"
   DO I=0,ORDER-1
    WRITE(6,'(I3,A)') I,"TH-ORDER GREEN'S FUNCTION (CUMULATIVE)"
    CALL DUMP5(GN2(:,:,I),IALL(0,0,0)-IVIRTCORE)
    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
     GREN(MOX,I)=GN2(MOX,MOX,I)
    ENDDO
   ENDDO

   WRITE(6,'(/,A,F20.15,A)') "METHOD 1: SINGLE-DETERMINANT N+1/N-1 WAVE FUNCTIONS"
   DO I=1,ORDER-1
    WRITE(6,'(I3,A)') I,"TH-ORDER SELF-ENERGY WITH NUMERICAL REDUCIBLE CANCELLATION"
    CALL DUMP5(SN1(:,:,I),IALL(0,0,0)-IVIRTCORE)
   ENDDO

   WRITE(6,'(/,A,F20.15,A)') "METHOD 1: SINGLE-DETERMINANT N+1/N-1 WAVE FUNCTIONS" 
   DO I=1,ORDER-1 
    WRITE(6,'(I3,A)') I,"TH-ORDER SELF-ENERGY (CUMULATIVE)" 
    CALL DUMP5(SN2(:,:,I),IALL(0,0,0)-IVIRTCORE)
!   WRITE(100,'(I3,F10.5,20F15.8:)') I,OMEGA,(SN2(MOX,MOX,I),MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE)
    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
     DIAG(MOX,I)=SN2(MOX,MOX,I)
    ENDDO
   ENDDO
!  WRITE(6,'(A)') 'fort.100 contains cumulative diagonal self-energies'

   ALLOCATE(AMAT(IALL(0,0,0)-IVIRTCORE-ICORE,IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(ER(IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(EI(IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(VR(IALL(0,0,0)-IVIRTCORE-ICORE,IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(VL(1,IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(WK(4*(IALL(0,0,0)-IVIRTCORE-ICORE)))
   WRITE(6,'(/,A,F20.15,A)') "METHOD 1: SINGLE-DETERMINANT N+1/N-1 WAVE FUNCTIONS" 
   DO I=1,ORDER-1 
    WRITE(6,'(I3,A)') I,"TH-ORDER SELF-ENERGY (EIGENVALUES)" 
    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
     DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
      AMAT(MOX-ICORE,MOY-ICORE)=SN2(MOX,MOY,I)
     ENDDO
    ENDDO
    CALL DGEEV('N','V',IALL(0,0,0)-IVIRTCORE-ICORE,AMAT,IALL(0,0,0)-IVIRTCORE-ICORE,ER,EI,VL,1,VR,&
     IALL(0,0,0)-IVIRTCORE-ICORE,WK,4*(IALL(0,0,0)-IVIRTCORE-ICORE),INFO)
    IF (INFO /= 0) CALL PABORT('DGEEV FAILED TO DIAGONALIZE A MATRIX')
    CALL PIKSRT(IALL(0,0,0)-IVIRTCORE-ICORE,IALL(0,0,0)-IVIRTCORE-ICORE,ER,VR,EI)
    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
     WRITE(6,'(F20.12)') ER(MOX-ICORE)
     EIGN(MOX,I)=ER(MOX-ICORE)
    ENDDO
!   WRITE(101,'(I3,F10.5,20F15.8:)') I,OMEGA,(ER(MOX-ICORE),MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE)
   ENDDO
   DEALLOCATE(AMAT,ER,EI,VR,VL,WK)
!  WRITE(6,'(A)') 'fort.101 contains cumulative eigen self-energies'

   ALLOCATE(B(IALL(0,0,0)-IVIRTCORE-ICORE,IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(BS(IALL(0,0,0)-IVIRTCORE-ICORE,IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(C(IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(CS(IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(INDX(IALL(0,0,0)-IVIRTCORE-ICORE))
   WRITE(6,'(/,A,F20.15,A)') "METHOD 1: SINGLE-DETERMINANT N+1/N-1 WAVE FUNCTIONS"
   DO I=0,ORDER-1
    B=0.0D0
    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
     B(MOX-ICORE,MOX-ICORE)=OMEGA
     DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
      B(MOX-ICORE,MOY-ICORE)=B(MOX-ICORE,MOY-ICORE)-SN2(MOX,MOY,I)
     ENDDO
    ENDDO
    BS=B
    CALL LUDCMP(B,IALL(0,0,0)-IVIRTCORE-ICORE,IALL(0,0,0)-IVIRTCORE-ICORE,INDX,X)
    GD=0.0D0
    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
     C=0.0D0
     C(MOX-ICORE)=1.0D0
     CS=C
     CALL LUBKSB(B,IALL(0,0,0)-IVIRTCORE-ICORE,IALL(0,0,0)-IVIRTCORE-ICORE,INDX,C)
     CALL MPROVE(BS,B,IALL(0,0,0)-IVIRTCORE-ICORE,IALL(0,0,0)-IVIRTCORE-ICORE,INDX,CS,C)
     DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
      GD(MOY,MOX,I)=C(MOY-ICORE)
     ENDDO
     GDYS(MOX,I)=GD(MOX,MOX,I)
    ENDDO
    WRITE(6,'(I3,A)') I,"TH-ORDER DYSON-GREEN'S FUNCTION (INFINITE PARTIAL SUM)"
    CALL DUMP5(GD(:,:,I),IALL(0,0,0)-IVIRTCORE)
!-- check ...
!   B=0.0d0
!   do mox=ICORE+1,IALL(0,0,0)-IVIRTCORE
!   do moy=ICORE+1,IALL(0,0,0)-IVIRTCORE
!   do pa=ICORE+1,IALL(0,0,0)-IVIRTCORE
!    B(mox-icore,moy-icore)=B(mox-icore,moy-icore)+gd(mox,pa,i)*bs(pa-icore,moy-icore)
!   enddo
!   enddo
!   enddo
!   write(6,*) 'Should be a unit matrix:'
!   call dump5(b,iall(0,0,0)-ivirtcore-icore)
!-- ... check
   ENDDO
   DEALLOCATE(B,BS,C,CS,INDX)

   DEALLOCATE(LKOOPIP,LKOOPEA,IPMAP,EAMAP)
   DEALLOCATE(GN,GN2,SN1,SN2,GD)
!  DEALLOCATE(SN3)
!  DEALLOCATE(SIP,SEA)
   DEALLOCATE(GIP,GEA)
   DEALLOCATE(EIP,EEA)
   DEALLOCATE(XIP,XEA)
   DEALLOCATE(VIP,VEA)
   DEALLOCATE(WFN,D,E0)
   CLOSE(50)
   CLOSE(51)

   RETURN

END SUBROUTINE



SUBROUTINE HIGHORDER_GF_REDUX2(OMEGA,ORDER)
! PERFORM HIGH-ORDER MBGF CALCULATIONS WITHOUT IP/EA PERTURBATION EXPANSION

   USE CONTROL
   USE GRADIENT
   USE STRUCTURE
   USE INTEGRAL
   USE BASISSET
   USE FULLCI

   IMPLICIT NONE
   DOUBLE PRECISION,PARAMETER :: TOL = 1.0D-10
   INTEGER,PARAMETER :: MAXFILE = 100
   DOUBLE PRECISION :: OMEGA
   INTEGER :: ORDER
   DOUBLE PRECISION,ALLOCATABLE :: IP_WFN(:,:,:,:),IP_E(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: EA_WFN(:,:,:,:),EA_E(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: E0(:)
   DOUBLE PRECISION,ALLOCATABLE :: GN(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: SN(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: SN1(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: SN2(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: GIP(:,:,:),GEA(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: SIP(:,:,:),SEA(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: D(:)
   DOUBLE PRECISION,ALLOCATABLE :: VIP(:,:,:),VEA(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: WFN(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: XIP(:,:,:),YIP(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: XEA(:,:,:),YEA(:,:,:)
   INTEGER :: I,J,K,L,M,N
   INTEGER :: PA,PB,QA,QB,RA,RB
   INTEGER :: XA,XB
   INTEGER :: MOX,MOY,MOW,MOZ
   INTEGER :: SGN
   DOUBLE PRECISION :: T1
   DOUBLE PRECISION,ALLOCATABLE :: VEC1(:,:),VEC2(:,:)
   LOGICAL,ALLOCATABLE :: LKOOPIP(:),LKOOPEA(:)
   INTEGER,ALLOCATABLE :: IPMAP(:),EAMAP(:)
   INTEGER :: INFO
   DOUBLE PRECISION,ALLOCATABLE :: AMAT(:,:),VL(:,:),VR(:,:),ER(:),EI(:),WK(:)

   WRITE(6,'(/,A,F20.15,A)') "GENERAL-ORDER MANY-BODY GREEN'S FUNCTION METHOD USING DEGENERATE RSPT N+1/N-1 WAVE FUNCTIONS"
   WRITE(6,'(A,F20.15,A)') 'SELF-ENERGY MATRIX WILL BE COMPUTED AT OMEGA=',OMEGA,' HARTREE'
   OPEN(50,FILE=TRIM(COPTN(1))//'.fi0',FORM='UNFORMATTED')
   OPEN(51,FILE=TRIM(COPTN(1))//'.fo0',FORM='UNFORMATTED')

   ! READ ALL N-ELECTRON WAVE FUNCTIONS
   ALLOCATE(WFN(NCF,NCF,0:ORDER-1))
   OPEN(97,FILE=TRIM(COPTN(1))//'.gf0',FORM='UNFORMATTED')
   REWIND(97)
   DO I=0,ORDER-1
    READ(97) WFN(:,:,I)
   ENDDO
   CLOSE(97)
   WRITE(6,'(A,I3,A,A)') '0 THRU',ORDER-1,'TH-ORDER N-ELECTRON WAVE FUNCTIONS READ FROM ',TRIM(COPTN(1))//'.gf0'
   ALLOCATE(E0(0:ORDER-1))
   E0(0)=MP_STORED(0)-NUCLEAR_REPULSION
   DO I=1,ORDER-1
    E0(I)=MP_STORED(I)-MP_STORED(I-1)
   ENDDO
!  WRITE(6,*) 'E0 VECTOR (ORDER IS SHIFTED BY 1)'
!  CALL DUMP16(E0,ORDER,1)

   ! INTERMEDIATE NORMALIZATION SEE EQ(2.28) OF SHAVITT & BARTLETT
   DO I=1,ORDER-1
    T1=0.0D0
    DO PA=1,NCF
     DO PB=1,NCF
      T1=T1+WFN(PB,PA,I)*WFN(PB,PA,0)
     ENDDO
    ENDDO
    IF (DABS(T1) > TOL) THEN
     WRITE(6,*) I,T1
     CALL WARNING('NONORTHOGONAL N')
    ENDIF
   ENDDO
   ! TRIPLE CHECK
   ALLOCATE(VEC1(NCF,NCF))
   DO I=0,ORDER-1
    VEC1=0.0D0
    DO J=0,I
     VEC1(:,:)=VEC1(:,:)+WFN(:,:,J)
    ENDDO
    REWIND(50)
    WRITE(50) VEC1
    CALL HAMILTONIAN_PRODUCT(50,51,0,2*MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE))
    REWIND(51)
    READ(51) VEC1
    T1=0.0D0
    DO PA=1,NCF
     DO PB=1,NCF
      T1=T1+VEC1(PB,PA)*WFN(PB,PA,0)
     ENDDO
    ENDDO
    IF (DABS(T1+NUCLEAR_REPULSION-MP_STORED(I+1)) > TOL) THEN
     WRITE(6,'(I3,2F20.15)') I, T1+NUCLEAR_REPULSION,MP_STORED(I+1)
     CALL PABORT('ENERGY INCONSISTENCY N')
    ENDIF
   ENDDO
   DEALLOCATE(VEC1)

!  ALLOCATE(E0(0:ORDER-1))
!  E0(0)=MP_STORED(0)-NUCLEAR_REPULSION
!  DO I=1,ORDER-1
!   E0(I)=MP_STORED(I)-MP_STORED(I-1)
!  ENDDO
!  WRITE(6,*) 'E0 VECTOR (ORDER IS SHIFTED BY 1)'
!  CALL DUMP16(E0,ORDER,1)

   ! READ ALL SMALL SIGMA
   ALLOCATE(IP_E(IP_NCF*NCF,0:ORDER-1))
   OPEN(60,FILE=TRIM(COPTN(1))//'.ip0',FORM='UNFORMATTED')
   REWIND(60)
   DO PA=1,NCF*IP_NCF
    DO I=0,ORDER-1
     READ(60) IP_E(PA,I)
    ENDDO
!write(6,'(I3,20F12.5:)') PA,(IP_E(PA,I),I=0,ORDER-1)
   ENDDO
   CLOSE(60)
   WRITE(6,'(A,A)') '(N-1)-ELECTRON STATIC SIGMA READ FROM ',TRIM(COPTN(1))//'.ip0'
   ALLOCATE(EA_E(EA_NCF*NCF,0:ORDER-1))
   OPEN(60,FILE=TRIM(COPTN(1))//'.ea0',FORM='UNFORMATTED')
   REWIND(60)
   DO PA=1,NCF*EA_NCF
    DO I=0,ORDER-1
     READ(60) EA_E(PA,I)
    ENDDO
!write(6,'(I3,20F12.5:)') PA,(EA_E(PA,I),I=0,ORDER-1)
   ENDDO
   CLOSE(60)
   WRITE(6,'(A,A)') '(N+1)-ELECTRON STATIC SIGMA READ FROM ',TRIM(COPTN(1))//'.ea0'

   ! READ ALL (N-1)-ELECTRON WAVE FUNCTIONS
   ALLOCATE(IP_WFN(IP_NCF,NCF,IP_NCF*NCF,0:ORDER-1))
   OPEN(60,FILE=TRIM(COPTN(1))//'.ip1',FORM='UNFORMATTED')
   REWIND(60)
   DO PA=1,NCF*IP_NCF
    DO I=0,ORDER-1
     READ(60) IP_WFN(:,:,PA,I)
    ENDDO
   ENDDO
   CLOSE(60)
   WRITE(6,'(A,I3,A,A)') '0 THRU',ORDER-1,'TH-ORDER (N-1)-ELECTRON WAVE FUNCTIONS READ FROM ',TRIM(COPTN(1))//'.ip1'

   ! READ ALL (N+1)-ELECTRON WAVE FUNCTIONS
   ALLOCATE(EA_WFN(EA_NCF,NCF,EA_NCF*NCF,0:ORDER-1))
   OPEN(60,FILE=TRIM(COPTN(1))//'.ea1',FORM='UNFORMATTED')
   REWIND(60)
   DO PA=1,NCF*EA_NCF
    DO I=0,ORDER-1
     READ(60) EA_WFN(:,:,PA,I)
    ENDDO
   ENDDO
   CLOSE(60)
   WRITE(6,'(A,I3,A,A)') '0 THRU',ORDER-1,'TH-ORDER (N+1)-ELECTRON WAVE FUNCTIONS READ FROM ',TRIM(COPTN(1))//'.ea1'

   ! FORM XIP,XEA
   ALLOCATE(XIP(NCF*IP_NCF,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(YIP(NCF*IP_NCF,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(XEA(NCF*EA_NCF,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(YEA(NCF*EA_NCF,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   XIP=0.0D0
   YIP=0.0D0
   XEA=0.0D0
   YEA=0.0D0
   DO I=0,ORDER-1
    DO K=0,I

     DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE

      ! IP SECTOR
      ALLOCATE(VEC1(IP_NCF,NCF))
      VEC1=0.0D0
      DO PA=1,NCF
       DO PB=1,NCF
        IF (.NOT.BTEST(CFHALF(PB),MOX-1)) CYCLE
        QA=PA
        QB=IP_ADDRSS(IBCLR(CFHALF(PB),MOX-1))
        SGN=1
        IF (MOX /= 1) THEN
         DO J=1,MOX-1
          IF (BTEST(CFHALF(PB),J-1)) SGN=-SGN
         ENDDO
        ENDIF
        VEC1(QB,QA)=VEC1(QB,QA)+WFN(PB,PA,K)*DFLOAT(SGN)
       ENDDO
      ENDDO
      DO PA=1,NCF*IP_NCF
       DO QA=1,NCF
        DO QB=1,IP_NCF
         XIP(PA,MOX,I)=XIP(PA,MOX,I)+VEC1(QB,QA)*IP_WFN(QB,QA,PA,I-K)
        ENDDO
       ENDDO
      ENDDO
      DO PA=1,NCF*IP_NCF
       DO QA=1,NCF
        DO QB=1,IP_NCF
         YIP(PA,MOX,I)=YIP(PA,MOX,I)+VEC1(QB,QA)*IP_WFN(QB,QA,PA,I-K)
        ENDDO
       ENDDO
      ENDDO
      DEALLOCATE(VEC1)

      ! EA SECTOR
      ALLOCATE(VEC1(EA_NCF,NCF))
      VEC1=0.0D0
      DO PA=1,NCF
       DO PB=1,NCF
        IF (BTEST(CFHALF(PB),MOX-1)) CYCLE
        QA=PA
        QB=EA_ADDRSS(IBSET(CFHALF(PB),MOX-1))
        SGN=1
        IF (MOX /= 1) THEN
         DO J=1,MOX-1
          IF (BTEST(CFHALF(PB),J-1)) SGN=-SGN
         ENDDO
        ENDIF
        VEC1(QB,QA)=VEC1(QB,QA)+WFN(PB,PA,K)*DFLOAT(SGN)
       ENDDO
      ENDDO
      DO PA=1,NCF*EA_NCF
       DO QA=1,NCF
        DO QB=1,EA_NCF
         YEA(PA,MOX,I)=YEA(PA,MOX,I)+VEC1(QB,QA)*EA_WFN(QB,QA,PA,I-K)
        ENDDO
       ENDDO
      ENDDO
      DO PA=1,NCF*EA_NCF
       DO QA=1,NCF
        DO QB=1,EA_NCF
         XEA(PA,MOX,I)=XEA(PA,MOX,I)+VEC1(QB,QA)*EA_WFN(QB,QA,PA,I-K)
        ENDDO
       ENDDO
      ENDDO
      DEALLOCATE(VEC1)

     ENDDO
    ENDDO

!   WRITE(6,*) "XIP VECTOR AT ORDER=",I
!   CALL DUMP16(XIP(:,:,I),NCF*IP_NCF,IALL(0,0,0)-IVIRTCORE)
!   WRITE(6,*) "XEA VECTOR AT ORDER=",I
!   CALL DUMP16(XEA(:,:,I),NCF*EA_NCF,IALL(0,0,0)-IVIRTCORE)

   ENDDO

   ! FORM VIP,VEA
   ALLOCATE(VIP(NCF*IP_NCF,NCF*IP_NCF,0:ORDER-1))
   ALLOCATE(VEA(NCF*EA_NCF,NCF*EA_NCF,0:ORDER-1))

   ! IP SECTOR
   ALLOCATE(VEC1(IP_NCF,NCF),VEC2(IP_NCF,NCF))
   VIP=0.0D0
   DO I=0,ORDER-1
    DO J=0,I
     DO PA=1,NCF*IP_NCF
      REWIND(50)
      WRITE(50) IP_WFN(:,:,PA,J)
      CALL IP_HAMILTONIAN_PRODUCT(50,51,0,MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE)+MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE+1))
      REWIND(51)
      READ(51) VEC1
      CALL H0_PRODUCT(50,51,1)
      REWIND(51)
      READ(51) VEC2
      DO K=0,I
       DO L=0,I
        IF (J+K+L==I) THEN
         IF (K==0) THEN
          ! omega + H0
          DO QA=1,NCF*IP_NCF
           DO RA=1,NCF
            DO RB=1,IP_NCF
             VIP(QA,PA,I)=VIP(QA,PA,I)+IP_WFN(RB,RA,QA,L)*OMEGA*IP_WFN(RB,RA,PA,J)
             VIP(QA,PA,I)=VIP(QA,PA,I)+IP_WFN(RB,RA,QA,L)*VEC2(RB,RA)
            ENDDO
           ENDDO
          ENDDO
         ELSE IF (K==1) THEN
          ! +V
          DO QA=1,NCF*IP_NCF
           DO RA=1,NCF
            DO RB=1,IP_NCF
             VIP(QA,PA,I)=VIP(QA,PA,I)+IP_WFN(RB,RA,QA,L)*(VEC1(RB,RA)-VEC2(RB,RA))
            ENDDO
           ENDDO
          ENDDO
         ENDIF
         ! -EN(K)
         DO QA=1,NCF*IP_NCF
          DO RA=1,NCF
           DO RB=1,IP_NCF
            VIP(QA,PA,I)=VIP(QA,PA,I)-IP_WFN(RB,RA,QA,L)*E0(K)*IP_WFN(RB,RA,PA,J)
           ENDDO
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO
     ENDDO
    ENDDO
!   WRITE(6,*) 'VIP MATRIX, ORDER=',I
!   CALL DUMP5(VIP(:,:,I),NCF*IP_NCF)
   ENDDO
   ! CHECK DIAGONALITY OF VIP(:,:,0)
   DO PA=1,NCF*IP_NCF
    DO QA=1,NCF*IP_NCF
     IF ((PA/=QA).AND.(DABS(VIP(PA,QA,0)) > TOL)) CALL PABORT('G0 NONDIAGONAL IN IP')
!    DO I=0,ORDER-1
!     IF (DABS(VIP(PA,QA,I)-VIP(QA,PA,I)) > TOL) THEN
!      WRITE(*,*) VIP(PA,QA,I),VIP(QA,PA,I)
!      CALL PABORT('G0 NONSYMMETRIC IN IP')
!     ENDIF
!    ENDDO
    ENDDO
   ENDDO
   DEALLOCATE(VEC1,VEC2)
   ! COMPARE DIAGONAL VIP VS SMALL SIGMA
!  WRITE(6,'(A,20I10:)') '   ',(I,I,I=0,ORDER-1)
!  DO PA=1,NCF*IP_NCF
!   WRITE(6,'(I3,20F10.6:)') PA,(IP_E(PA,I),VIP(PA,PA,I),I=0,ORDER-1)
!  ENDDO

   ! EA SECTOR
   ALLOCATE(VEC1(EA_NCF,NCF),VEC2(EA_NCF,NCF))
   VEA=0.0D0
   DO I=0,ORDER-1
    DO J=0,I
     DO PA=1,NCF*EA_NCF
      REWIND(50)
      WRITE(50) EA_WFN(:,:,PA,J)
      CALL EA_HAMILTONIAN_PRODUCT(50,51,0,MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE)+MIN(IOCC-ICORE+1,IALL(0,0,0)-IOCC-IVIRTCORE))
      REWIND(51)
      READ(51) VEC1
      CALL H0_PRODUCT(50,51,-1)
      REWIND(51)
      READ(51) VEC2
      DO K=0,I
       DO L=0,I
        IF (J+K+L==I) THEN
         IF (K==0) THEN
          ! omega - H0
          DO QA=1,NCF*EA_NCF
           DO RA=1,NCF
            DO RB=1,EA_NCF
             VEA(QA,PA,I)=VEA(QA,PA,I)+EA_WFN(RB,RA,QA,L)*OMEGA*EA_WFN(RB,RA,PA,J)
             VEA(QA,PA,I)=VEA(QA,PA,I)-EA_WFN(RB,RA,QA,L)*VEC2(RB,RA)
            ENDDO
           ENDDO
          ENDDO
         ELSE IF (K==1) THEN
          ! -V
          DO QA=1,NCF*EA_NCF
           DO RA=1,NCF
            DO RB=1,EA_NCF
             VEA(QA,PA,I)=VEA(QA,PA,I)-EA_WFN(RB,RA,QA,L)*(VEC1(RB,RA)-VEC2(RB,RA))
            ENDDO
           ENDDO
          ENDDO
         ENDIF
         ! +EN(K)
         DO QA=1,NCF*EA_NCF
          DO RA=1,NCF
           DO RB=1,EA_NCF
            VEA(QA,PA,I)=VEA(QA,PA,I)+EA_WFN(RB,RA,QA,L)*E0(K)*EA_WFN(RB,RA,PA,J)
           ENDDO
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO
     ENDDO
    ENDDO
!   WRITE(6,*) 'VEA MATRIX, ORDER=',I
!   CALL DUMP5(VEA(:,:,I),NCF*EA_NCF)
   ENDDO
   ! CHECK DIAGONALITY OF VEA(:,:,0)
   DO PA=1,NCF*EA_NCF
    DO QA=1,NCF*EA_NCF
     IF ((PA/=QA).AND.(DABS(VEA(PA,QA,0)) > TOL)) CALL PABORT('G0 NONDIAGONAL IN EA')
!    DO I=0,ORDER-1
!     IF (DABS(VEA(PA,QA,I)-VEA(QA,PA,I)) > TOL) THEN
!      WRITE(*,*) VEA(PA,QA,I),VEA(QA,PA,I)
!      CALL PABORT('G0 NONSYMMETRIC IN EA')
!     ENDIF
!    ENDDO
    ENDDO
   ENDDO
   DEALLOCATE(VEC1,VEC2)
   ! COMPARE DIAGONAL VEA VS SMALL SIGMA
!  WRITE(6,'(A,20I10:)') '   ',(I,I,I=0,ORDER-1)
!  DO PA=1,NCF*EA_NCF
!   WRITE(6,'(I3,20F10.6:)') PA,(EA_E(PA,I),VEA(PA,PA,I),I=0,ORDER-1)
!  ENDDO

   ! FORM D
   ALLOCATE(D(0:ORDER-1))
   DO I=0,ORDER-1
    D(I)=0.0D0
    DO J=0,I
     DO PA=1,NCF
      DO PB=1,NCF
       D(I)=D(I)+WFN(PB,PA,J)*WFN(PB,PA,I-J)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
!  WRITE(6,*) 'D VECTOR'
!  CALL DUMP16(D,ORDER,1)

   ! FORM LKOOPIP AND LKOOPEA
   ALLOCATE(LKOOPIP(NCF*IP_NCF))
   ALLOCATE(LKOOPEA(NCF*EA_NCF))
   LKOOPIP=.FALSE.
   DO PA=1,NCF
    DO PB=1,IP_NCF
     IF (NORDER(PA)+IP_NORDER(PB) < 2) LKOOPIP((PA-1)*IP_NCF+PB)=.TRUE.
    ENDDO
   ENDDO
   LKOOPEA=.FALSE.
   DO PA=1,NCF
    DO PB=1,EA_NCF
     IF (NORDER(PA)+EA_NORDER(PB) < 2) LKOOPEA((PA-1)*EA_NCF+PB)=.TRUE.
    ENDDO
   ENDDO

   ! MAP ORBITALS TO DETERMINANTS
   ALLOCATE(IPMAP(IALL(0,0,0)-IVIRTCORE),EAMAP(IALL(0,0,0)-IVIRTCORE))
   IPMAP=99999999
   EAMAP=99999999
!  WRITE(6,'(A)') 'ORBITAL -> KOOPMANS DETERMINANT'
   DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
    IF (MOX <= IOCC) THEN
     DO PA=1,IP_NCF
      IF ((IP_NORDER(PA)==1).AND.(.NOT.BTEST(IP_CFHALF(PA),MOX-1))) IPMAP(MOX)=PA
     ENDDO
    ELSE
     DO PA=1,EA_NCF
      IF ((EA_NORDER(PA)==1).AND.(BTEST(EA_CFHALF(PA),MOX-1))) EAMAP(MOX)=PA
     ENDDO
    ENDIF
!   WRITE(6,'(3I5)') MOX,IPMAP(MOX),EAMAP(MOX)
   ENDDO

   ALLOCATE(GN(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(SN(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(SN1(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(SN2(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(GIP(NCF*IP_NCF,NCF*IP_NCF,0:ORDER-1))
   ALLOCATE(GEA(NCF*EA_NCF,NCF*EA_NCF,0:ORDER-1))
   ALLOCATE(SIP(NCF*IP_NCF,NCF*IP_NCF,0:ORDER-1))
   ALLOCATE(SEA(NCF*EA_NCF,NCF*EA_NCF,0:ORDER-1))
   GIP=1.0D99 ! safety
   GEA=1.0D99 ! safety
   SIP=1.0D99 ! safety
   SEA=1.0D99 ! safety
   SN1=1.0D99
   SN2=1.0D99

   ! GRAND LOOP
   DO I=0,ORDER-1
   
    ! FORM GIP AND GEA 
    IF (I==0) THEN

     ! IP SECTOR
     GIP(:,:,0)=0.0D0
     DO PA=1,NCF*IP_NCF
!     GIP(PA,PA,0)=1.0D0
      GIP(PA,PA,0)=1.0D0/(OMEGA-IP_E(PA,0))
     ENDDO
     ! EA SECTOR
     GEA(:,:,0)=0.0D0
     DO PA=1,NCF*EA_NCF
!     GEA(PA,PA,0)=1.0D0
      GEA(PA,PA,0)=1.0D0/(OMEGA-EA_E(PA,0))
     ENDDO
!    WRITE(6,*) 'ZEROTH-ORDER GIP'
!    CALL DUMP5(GIP(:,:,0),NCF*IP_NCF)
!    WRITE(6,*) 'ZEROTH-ORDER GEA'
!    CALL DUMP5(GEA(:,:,0),NCF*EA_NCF)

    ELSE

     ! IP SECTOR
     DO N=1,NCF*IP_NCF
      DO M=1,NCF*IP_NCF
       GIP(N,M,I)=0.0D0
       DO J=1,I
        DO L=1,NCF*IP_NCF
         GIP(N,M,I)=GIP(N,M,I)-GIP(N,L,I-J)*VIP(L,M,J)*GIP(M,M,0)
!        GIP(N,M,I)=GIP(N,M,I)-GIP(N,L,I-J)*VIP(M,L,J)
        ENDDO
!       GIP(N,M,I)=GIP(N,M,I)+GIP(N,M,I-J)*IP_E(M,J)*GIP(M,M,0)
       ENDDO
      ENDDO
     ENDDO

     ! EA SECTOR
     DO N=1,NCF*EA_NCF
      DO M=1,NCF*EA_NCF
       GEA(N,M,I)=0.0D0
       DO J=1,I
        DO L=1,NCF*EA_NCF
         GEA(N,M,I)=GEA(N,M,I)-GEA(N,L,I-J)*VEA(L,M,J)*GEA(M,M,0)
!        GEA(N,M,I)=GEA(N,M,I)+GEA(N,L,I-J)*VEA(L,M,J)
        ENDDO
!       GEA(N,M,I)=GEA(N,M,I)+GEA(N,M,I-J)*EA_E(M,J)*GEA(M,M,0)
       ENDDO
      ENDDO
     ENDDO

    ENDIF

    ! FORM GN
    GN(:,:,I)=0.0D0
    DO J=0,I
     DO K=0,I
      DO L=0,I
       IF (J+K+L==I) THEN
        DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
         DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
          ! IP SECTOR
          DO PA=1,NCF*IP_NCF
           DO PB=1,NCF*IP_NCF
            GN(MOX,MOY,I)=GN(MOX,MOY,I)+XIP(PA,MOX,J)*GIP(PA,PB,K)*YIP(PB,MOY,L) 
           ENDDO
          ENDDO
          ! EA SECTOR
          DO PA=1,NCF*EA_NCF
           DO PB=1,NCF*EA_NCF
            GN(MOX,MOY,I)=GN(MOX,MOY,I)+YEA(PA,MOY,J)*GEA(PA,PB,K)*XEA(PB,MOX,L)
           ENDDO
          ENDDO
         ENDDO 
        ENDDO
       ENDIF
      ENDDO
     ENDDO
    ENDDO
    IF (I > 0) THEN
     DO J=1,I
      DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
       DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
        GN(MOX,MOY,I)=GN(MOX,MOY,I)-D(J)*GN(MOX,MOY,I-J)
       ENDDO
      ENDDO
     ENDDO
    ENDIF

    ! FORM SN WITH NUMERICAL REDUCIBLE DELETION
    IF (I==0) THEN
     DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
      DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
       SN2(MOX,MOY,I)=0.0D0
      ENDDO
      SN2(MOX,MOX,I)=EPSILON(MOX,0,0,0)
     ENDDO
!    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
!     DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
!      SN1(MOX,MOY,I)=0.0D0
!     ENDDO
!     SN1(MOX,MOX,I)=OMEGA-EPSILON(MOX,0,0,0)
!    ENDDO
    ELSE IF (I==1) THEN
     DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
      DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
!      SN1(MOX,MOY,I)=0.0d0
       SN1(MOX,MOY,I)=GN(MOX,MOY,1)/GN(MOX,MOX,0)/GN(MOY,MOY,0)
       SN2(MOX,MOY,I)=SN2(MOX,MOY,I-1)+SN1(MOX,MOY,I)
      ENDDO
     ENDDO
    ELSE
     DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
      DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
       SN1(MOX,MOY,I)=GN(MOX,MOY,I)
       DO J=1,I-1
        DO MOZ=ICORE+1,IALL(0,0,0)-IVIRTCORE
         DO MOW=ICORE+1,IALL(0,0,0)-IVIRTCORE
! Note: the renormalization term below makes difference only at fourth and higher orders
          SN1(MOX,MOY,I)=SN1(MOX,MOY,I)-GN(MOX,MOZ,0)*SN1(MOZ,MOW,J)*GN(MOW,MOY,I-J)
         ENDDO
        ENDDO
       ENDDO
       SN1(MOX,MOY,I)=SN1(MOX,MOY,I)/GN(MOX,MOX,0)/GN(MOY,MOY,0)
       SN2(MOX,MOY,I)=SN2(MOX,MOY,I-1)+SN1(MOX,MOY,I)
      ENDDO
     ENDDO
    ENDIF

   ENDDO

   WRITE(6,'(/,A,F20.15,A)') "METHOD 2: DEGENERATE RSPT N+1/N-1 WAVE FUNCTIONS"
   DO I=0,ORDER-1
    WRITE(6,'(I3,A)') I,"TH-ORDER GREEN'S FUNCTION"
    CALL DUMP5(GN(:,:,I),IALL(0,0,0)-IVIRTCORE)
   ENDDO

   WRITE(6,'(/,A,F20.15,A)') "METHOD 2: DEGENERATE RSPT N+1/N-1 WAVE FUNCTIONS"
   DO I=1,ORDER-1
    WRITE(6,'(I3,A)') I,"TH-ORDER SELF-ENERGY WITH NUMERICAL REDUCIBLE CANCELLATION"
    CALL DUMP5(SN1(:,:,I),IALL(0,0,0)-IVIRTCORE)
   ENDDO

   WRITE(6,'(/,A,F20.15,A)') "METHOD 2: DEGENERATE RSPT N+1/N-1 WAVE FUNCTIONS"
   DO I=1,ORDER-1
    WRITE(6,'(I3,A)') I,"TH-ORDER SELF-ENERGY (CUMULATIVE)"
    CALL DUMP5(SN2(:,:,I),IALL(0,0,0)-IVIRTCORE)
   ENDDO

   ALLOCATE(AMAT(IALL(0,0,0)-IVIRTCORE-ICORE,IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(ER(IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(EI(IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(VR(IALL(0,0,0)-IVIRTCORE-ICORE,IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(VL(1,IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(WK(4*(IALL(0,0,0)-IVIRTCORE-ICORE)))
   WRITE(6,'(/,A,F20.15,A)') "METHOD 2: DEGENERATE RSPT N+1/N-1 WAVE FUNCTIONS"
   DO I=1,ORDER-1
    WRITE(6,'(I3,A)') I,"TH-ORDER SELF-ENERGY (EIGENVALUES)"
    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
     DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
      AMAT(MOX-ICORE,MOY-ICORE)=SN2(MOX,MOY,I)
     ENDDO
    ENDDO
    CALL DGEEV('N','V',IALL(0,0,0)-IVIRTCORE-ICORE,AMAT,IALL(0,0,0)-IVIRTCORE-ICORE,ER,EI,VL,1,VR,&
     IALL(0,0,0)-IVIRTCORE-ICORE,WK,4*(IALL(0,0,0)-IVIRTCORE-ICORE),INFO)
    IF (INFO /= 0) CALL PABORT('DGEEV FAILED TO DIAGONALIZE A MATRIX')
    CALL PIKSRT(IALL(0,0,0)-IVIRTCORE-ICORE,IALL(0,0,0)-IVIRTCORE-ICORE,ER,VR,EI)
    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
     WRITE(6,'(F20.12)') ER(MOX-ICORE)
    ENDDO
!   WRITE(101,'(I3,F10.5,20F15.8:)') I,OMEGA,(ER(MOX-ICORE),MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE)
   ENDDO
   DEALLOCATE(AMAT,ER,EI,VR,VL,WK)
!  WRITE(6,'(A)') 'fort.101 contains cumulative eigen self-energies'

!  WRITE(6,'(/,A,F20.15,A)') "METHOD 2: DEGENERATE RSPT N+1/N-1 WAVE FUNCTIONS"
!  WRITE(6,'(/,A,F20.15,A)') " =========== W R O N G !!! =========== "
!  DO I=0,ORDER-1
!   WRITE(6,'(I3,A)') I,"TH-ORDER SELF-ENERGY WITH NUMERICAL REDUCIBLE CANCELLATION"
!   CALL DUMP5(SN(:,:,I),IALL(0,0,0)-IVIRTCORE)
!  ENDDO

   DEALLOCATE(IP_WFN,EA_WFN,IP_E,EA_E)
   DEALLOCATE(LKOOPIP,LKOOPEA,IPMAP,EAMAP)
   DEALLOCATE(SN,SN1,SN2)
   DEALLOCATE(SIP,SEA)
   DEALLOCATE(GN)
   DEALLOCATE(GIP,GEA)
   DEALLOCATE(XIP,XEA)
   DEALLOCATE(VIP,VEA)
   DEALLOCATE(WFN,D)
   CLOSE(50)
   CLOSE(51)

   RETURN

END SUBROUTINE



SUBROUTINE HIGHORDER_GF_REDUX3(OMEGA,ORDER)
! PERFORM HIGH-ORDER MBGF CALCULATIONS WITHOUT IP/EA PERTURBATION EXPANSION

   USE CONTROL
   USE GRADIENT
   USE STRUCTURE
   USE INTEGRAL
   USE BASISSET
   USE FULLCI

   IMPLICIT NONE
   DOUBLE PRECISION,PARAMETER :: TOL = 1.0D-10
   INTEGER,PARAMETER :: MAXFILE = 100
   DOUBLE PRECISION :: OMEGA
   INTEGER :: ORDER
   DOUBLE PRECISION,ALLOCATABLE :: IP_WFN(:,:,:,:),IP_E(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: EA_WFN(:,:,:,:),EA_E(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: E0(:)
   DOUBLE PRECISION,ALLOCATABLE :: GN(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: SN(:,:,:),SN_CUM(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: GIP(:,:,:),GEA(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: D(:)
   DOUBLE PRECISION,ALLOCATABLE :: VIP(:,:,:),VEA(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: WFN(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: XIP(:,:,:),YIP(:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: XEA(:,:,:),YEA(:,:,:)
   INTEGER :: I,J,K,L,M,N
   INTEGER :: PA,PB,QA,QB,RA,RB
   INTEGER :: XA,XB
   INTEGER :: MOX,MOY,MOW,MOZ
   INTEGER :: SGN
   DOUBLE PRECISION :: T1
   DOUBLE PRECISION,ALLOCATABLE :: VEC1(:,:),VEC2(:,:)
   LOGICAL,ALLOCATABLE :: LKOOPIP(:),LKOOPEA(:)
   INTEGER :: INFO
   DOUBLE PRECISION,ALLOCATABLE :: AMAT(:,:),VL(:,:),VR(:,:),ER(:),EI(:),WK(:)
   double precision :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7
   double precision,allocatable :: ta1(:,:),ta2(:,:)

   WRITE(6,'(/,A,F20.15,A)') "GENERAL-ORDER MANY-BODY GREEN'S FUNCTION METHOD USING BIORTHOGONAL WAVE FUNCTIONS"
   WRITE(6,'(A,F20.15,A)') 'SELF-ENERGY MATRIX WILL BE COMPUTED AT OMEGA=',OMEGA,' HARTREE'
   OPEN(50,FILE=TRIM(COPTN(1))//'.fi0',FORM='UNFORMATTED')
   OPEN(51,FILE=TRIM(COPTN(1))//'.fo0',FORM='UNFORMATTED')

   ! READ ALL N-ELECTRON WAVE FUNCTIONS
   ALLOCATE(WFN(NCF,NCF,0:ORDER-1))
   OPEN(97,FILE=TRIM(COPTN(1))//'.gf0',FORM='UNFORMATTED')
   REWIND(97)
   DO I=0,ORDER-1
    READ(97) WFN(:,:,I)
   ENDDO
   CLOSE(97)
   WRITE(6,'(A,I3,A,A)') '0 THRU',ORDER-1,'TH-ORDER N-ELECTRON WAVE FUNCTIONS READ FROM ',TRIM(COPTN(1))//'.gf0'
   ALLOCATE(E0(0:ORDER-1))
   E0(0)=MP_STORED(0)-NUCLEAR_REPULSION
   DO I=1,ORDER-1
    E0(I)=MP_STORED(I)-MP_STORED(I-1)
   ENDDO
!  WRITE(6,*) 'E0 VECTOR (ORDER IS SHIFTED BY 1)'
!  CALL DUMP16(E0,ORDER,1)

   ! INTERMEDIATE NORMALIZATION SEE EQ(2.28) OF SHAVITT & BARTLETT
   DO I=1,ORDER-1
    T1=0.0D0
    DO PA=1,NCF
     DO PB=1,NCF
      T1=T1+WFN(PB,PA,I)*WFN(PB,PA,0)
     ENDDO
    ENDDO
    IF (DABS(T1) > TOL) THEN
     WRITE(6,*) I,T1
     CALL WARNING('NONORTHOGONAL N')
    ENDIF
   ENDDO
   ! TRIPLE CHECK
   ALLOCATE(VEC1(NCF,NCF))
   DO I=0,ORDER-1
    VEC1=0.0D0
    DO J=0,I
     VEC1(:,:)=VEC1(:,:)+WFN(:,:,J)
    ENDDO
    REWIND(50)
    WRITE(50) VEC1
    CALL HAMILTONIAN_PRODUCT(50,51,0,2*MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE))
    REWIND(51)
    READ(51) VEC1
    T1=0.0D0
    DO PA=1,NCF
     DO PB=1,NCF
      T1=T1+VEC1(PB,PA)*WFN(PB,PA,0)
     ENDDO
    ENDDO
    IF (DABS(T1+NUCLEAR_REPULSION-MP_STORED(I+1)) > TOL) THEN
     WRITE(6,'(I3,2F20.15)') I, T1+NUCLEAR_REPULSION,MP_STORED(I+1)
     CALL PABORT('ENERGY INCONSISTENCY N')
    ENDIF
   ENDDO
   DEALLOCATE(VEC1)

!  ALLOCATE(E0(0:ORDER-1))
!  E0(0)=MP_STORED(0)-NUCLEAR_REPULSION
!  DO I=1,ORDER-1
!   E0(I)=MP_STORED(I)-MP_STORED(I-1)
!  ENDDO
!  WRITE(6,*) 'E0 VECTOR (ORDER IS SHIFTED BY 1)'
!  CALL DUMP16(E0,ORDER,1)

   ! READ ALL SMALL SIGMA
   ALLOCATE(IP_E(IP_NCF*NCF,0:ORDER-1))
   OPEN(60,FILE=TRIM(COPTN(1))//'.ip0',FORM='UNFORMATTED')
   REWIND(60)
   DO PA=1,NCF*IP_NCF
    DO I=0,ORDER-1
     READ(60) IP_E(PA,I)
    ENDDO
!write(6,'(I3,20F12.5:)') PA,(IP_E(PA,I),I=0,ORDER-1)
   ENDDO
   CLOSE(60)
   WRITE(6,'(A,A)') '(N-1)-ELECTRON STATIC SIGMA READ FROM ',TRIM(COPTN(1))//'.ip0'
   ALLOCATE(EA_E(EA_NCF*NCF,0:ORDER-1))
   OPEN(60,FILE=TRIM(COPTN(1))//'.ea0',FORM='UNFORMATTED')
   REWIND(60)
   DO PA=1,NCF*EA_NCF
    DO I=0,ORDER-1
     READ(60) EA_E(PA,I)
    ENDDO
!write(6,'(I3,20F12.5:)') PA,(EA_E(PA,I),I=0,ORDER-1)
   ENDDO
   CLOSE(60)
   WRITE(6,'(A,A)') '(N+1)-ELECTRON STATIC SIGMA READ FROM ',TRIM(COPTN(1))//'.ea0'

   ! READ ALL (N-1)-ELECTRON WAVE FUNCTIONS
   ALLOCATE(IP_WFN(IP_NCF,NCF,IP_NCF*NCF,0:ORDER-1))
   OPEN(60,FILE=TRIM(COPTN(1))//'.ip1',FORM='UNFORMATTED')
   REWIND(60)
   DO PA=1,NCF*IP_NCF
    DO I=0,ORDER-1
     READ(60) IP_WFN(:,:,PA,I)
    ENDDO
   ENDDO
   CLOSE(60)
   WRITE(6,'(A,I3,A,A)') '0 THRU',ORDER-1,'TH-ORDER (N-1)-ELECTRON WAVE FUNCTIONS READ FROM ',TRIM(COPTN(1))//'.ip1'

   ! READ ALL (N+1)-ELECTRON WAVE FUNCTIONS
   ALLOCATE(EA_WFN(EA_NCF,NCF,EA_NCF*NCF,0:ORDER-1))
   OPEN(60,FILE=TRIM(COPTN(1))//'.ea1',FORM='UNFORMATTED')
   REWIND(60)
   DO PA=1,NCF*EA_NCF
    DO I=0,ORDER-1
     READ(60) EA_WFN(:,:,PA,I)
    ENDDO
   ENDDO
   CLOSE(60)
   WRITE(6,'(A,I3,A,A)') '0 THRU',ORDER-1,'TH-ORDER (N+1)-ELECTRON WAVE FUNCTIONS READ FROM ',TRIM(COPTN(1))//'.ea1'

   ! FORM XIP,XEA
   ALLOCATE(XIP(NCF*IP_NCF,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(YIP(NCF*IP_NCF,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(XEA(NCF*EA_NCF,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(YEA(NCF*EA_NCF,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   XIP=0.0D0
   YIP=0.0D0
   XEA=0.0D0
   YEA=0.0D0
   DO I=0,ORDER-1
    DO K=0,I

     DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE

      ! IP SECTOR
      ALLOCATE(VEC1(IP_NCF,NCF))
      VEC1=0.0D0
      DO PA=1,NCF
       DO PB=1,NCF
        IF (.NOT.BTEST(CFHALF(PB),MOX-1)) CYCLE
        QA=PA
        QB=IP_ADDRSS(IBCLR(CFHALF(PB),MOX-1))
        SGN=1
        IF (MOX /= 1) THEN
         DO J=1,MOX-1
          IF (BTEST(CFHALF(PB),J-1)) SGN=-SGN
         ENDDO
        ENDIF
        VEC1(QB,QA)=VEC1(QB,QA)+WFN(PB,PA,K)*DFLOAT(SGN)
       ENDDO
      ENDDO
      DO PA=1,NCF*IP_NCF
       DO QA=1,NCF
        DO QB=1,IP_NCF
         XIP(PA,MOX,I)=XIP(PA,MOX,I)+VEC1(QB,QA)*IP_WFN(QB,QA,PA,I-K)
        ENDDO
       ENDDO
      ENDDO
! aymmetric
if (k==i) then
      DO PA=1,NCF*IP_NCF
       DO QA=1,NCF
        DO QB=1,IP_NCF
         YIP(PA,MOX,I)=YIP(PA,MOX,I)+VEC1(QB,QA)*IP_WFN(QB,QA,PA,I-K)
        ENDDO
       ENDDO
      ENDDO
! aymmetric
endif
      DEALLOCATE(VEC1)

      ! EA SECTOR
      ALLOCATE(VEC1(EA_NCF,NCF))
      VEC1=0.0D0
      DO PA=1,NCF
       DO PB=1,NCF
        IF (BTEST(CFHALF(PB),MOX-1)) CYCLE
        QA=PA
        QB=EA_ADDRSS(IBSET(CFHALF(PB),MOX-1))
        SGN=1
        IF (MOX /= 1) THEN
         DO J=1,MOX-1
          IF (BTEST(CFHALF(PB),J-1)) SGN=-SGN
         ENDDO
        ENDIF
        VEC1(QB,QA)=VEC1(QB,QA)+WFN(PB,PA,K)*DFLOAT(SGN)
       ENDDO
      ENDDO
      DO PA=1,NCF*EA_NCF
       DO QA=1,NCF
        DO QB=1,EA_NCF
         YEA(PA,MOX,I)=YEA(PA,MOX,I)+VEC1(QB,QA)*EA_WFN(QB,QA,PA,I-K)
        ENDDO
       ENDDO
      ENDDO
! aymmetric
if (k==i) then
      DO PA=1,NCF*EA_NCF
       DO QA=1,NCF
        DO QB=1,EA_NCF
         XEA(PA,MOX,I)=XEA(PA,MOX,I)+VEC1(QB,QA)*EA_WFN(QB,QA,PA,I-K)
        ENDDO
       ENDDO
      ENDDO
! aymmetric
endif
      DEALLOCATE(VEC1)

     ENDDO
    ENDDO

!   WRITE(6,*) "XIP VECTOR AT ORDER=",I
!   CALL DUMP16(XIP(:,:,I),NCF*IP_NCF,IALL(0,0,0)-IVIRTCORE)
!   WRITE(6,*) "XEA VECTOR AT ORDER=",I
!   CALL DUMP16(XEA(:,:,I),NCF*EA_NCF,IALL(0,0,0)-IVIRTCORE)
!   WRITE(6,*) "YIP VECTOR AT ORDER=",I
!   CALL DUMP16(YIP(:,:,I),NCF*IP_NCF,IALL(0,0,0)-IVIRTCORE)
!   WRITE(6,*) "YEA VECTOR AT ORDER=",I
!   CALL DUMP16(YEA(:,:,I),NCF*EA_NCF,IALL(0,0,0)-IVIRTCORE)

   ENDDO

   ! FORM VIP,VEA
   ALLOCATE(VIP(NCF*IP_NCF,NCF*IP_NCF,0:ORDER-1))
   ALLOCATE(VEA(NCF*EA_NCF,NCF*EA_NCF,0:ORDER-1))

   ! IP SECTOR
   ALLOCATE(VEC1(IP_NCF,NCF),VEC2(IP_NCF,NCF))
   VIP=0.0D0
   DO I=0,ORDER-1
    DO J=0,I
     DO PA=1,NCF*IP_NCF
      REWIND(50)
      WRITE(50) IP_WFN(:,:,PA,J)
      CALL IP_HAMILTONIAN_PRODUCT(50,51,0,MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE)+MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE+1))
      REWIND(51)
      READ(51) VEC1
      CALL H0_PRODUCT(50,51,1)
      REWIND(51)
      READ(51) VEC2
      DO K=0,I
       DO L=0,I
        IF (J+K+L==I) THEN
! asymmetric
if (l==0) then
         IF (K==0) THEN
          ! omega + H0
          DO QA=1,NCF*IP_NCF
           DO RA=1,NCF
            DO RB=1,IP_NCF
             VIP(QA,PA,I)=VIP(QA,PA,I)+IP_WFN(RB,RA,QA,L)*OMEGA*IP_WFN(RB,RA,PA,J)
!            VIP(QA,PA,I)=VIP(QA,PA,I)+IP_WFN(RB,RA,QA,L)*VEC2(RB,RA)
            ENDDO
           ENDDO
          ENDDO
         ELSE IF (K==1) THEN
          ! +V
          DO QA=1,NCF*IP_NCF
           DO RA=1,NCF
            DO RB=1,IP_NCF
!            VIP(QA,PA,I)=VIP(QA,PA,I)+IP_WFN(RB,RA,QA,L)*(VEC1(RB,RA)-VEC2(RB,RA))
            ENDDO
           ENDDO
          ENDDO
         ENDIF
         ! -EN(K)
         DO QA=1,NCF*IP_NCF
          DO RA=1,NCF
           DO RB=1,IP_NCF
!           VIP(QA,PA,I)=VIP(QA,PA,I)-IP_WFN(RB,RA,QA,L)*E0(K)*IP_WFN(RB,RA,PA,J)
            VIP(QA,PA,I)=VIP(QA,PA,I)-IP_WFN(RB,RA,QA,L)*IP_E(PA,K)*IP_WFN(RB,RA,PA,J)
           ENDDO
          ENDDO
         ENDDO
        ENDIF
! asymmmetric
endif
       ENDDO
      ENDDO
     ENDDO
    ENDDO
!   WRITE(6,*) 'VIP MATRIX, ORDER=',I
!   CALL DUMP5(VIP(:,:,I),NCF*IP_NCF)
   ENDDO
   ! CHECK DIAGONALITY OF VIP(:,:,0)
   DO PA=1,NCF*IP_NCF
    DO QA=1,NCF*IP_NCF
     IF ((PA/=QA).AND.(DABS(VIP(PA,QA,0)) > TOL)) CALL PABORT('G0 NONDIAGONAL IN IP')
!    DO I=0,ORDER-1
!     IF (DABS(VIP(PA,QA,I)-VIP(QA,PA,I)) > TOL) THEN
!      WRITE(*,*) VIP(PA,QA,I),VIP(QA,PA,I)
!      CALL PABORT('G0 NONSYMMETRIC IN IP')
!     ENDIF
!    ENDDO
    ENDDO
   ENDDO
   DEALLOCATE(VEC1,VEC2)
   ! COMPARE DIAGONAL VIP VS SMALL SIGMA
!  WRITE(6,'(A,20I10:)') '   ',(I,I,I=0,ORDER-1)
!  DO PA=1,NCF*IP_NCF
!   WRITE(6,'(I3,20F10.6:)') PA,(IP_E(PA,I),VIP(PA,PA,I),I=0,ORDER-1)
!  ENDDO

   ! EA SECTOR
   ALLOCATE(VEC1(EA_NCF,NCF),VEC2(EA_NCF,NCF))
   VEA=0.0D0
   DO I=0,ORDER-1
    DO J=0,I
     DO PA=1,NCF*EA_NCF
      REWIND(50)
      WRITE(50) EA_WFN(:,:,PA,J)
      CALL EA_HAMILTONIAN_PRODUCT(50,51,0,MIN(IOCC-ICORE,IALL(0,0,0)-IOCC-IVIRTCORE)+MIN(IOCC-ICORE+1,IALL(0,0,0)-IOCC-IVIRTCORE))
      REWIND(51)
      READ(51) VEC1
      CALL H0_PRODUCT(50,51,-1)
      REWIND(51)
      READ(51) VEC2
      DO K=0,I
       DO L=0,I
        IF (J+K+L==I) THEN
! asymmetric
if (l==0) then
         IF (K==0) THEN
          ! omega - H0
          DO QA=1,NCF*EA_NCF
           DO RA=1,NCF
            DO RB=1,EA_NCF
             VEA(QA,PA,I)=VEA(QA,PA,I)+EA_WFN(RB,RA,QA,L)*OMEGA*EA_WFN(RB,RA,PA,J)
!            VEA(QA,PA,I)=VEA(QA,PA,I)-EA_WFN(RB,RA,QA,L)*VEC2(RB,RA)
            ENDDO
           ENDDO
          ENDDO
         ELSE IF (K==1) THEN
          ! -V
          DO QA=1,NCF*EA_NCF
           DO RA=1,NCF
            DO RB=1,EA_NCF
!            VEA(QA,PA,I)=VEA(QA,PA,I)-EA_WFN(RB,RA,QA,L)*(VEC1(RB,RA)-VEC2(RB,RA))
            ENDDO
           ENDDO
          ENDDO
         ENDIF
         ! +EN(K)
         DO QA=1,NCF*EA_NCF
          DO RA=1,NCF
           DO RB=1,EA_NCF
!           VEA(QA,PA,I)=VEA(QA,PA,I)+EA_WFN(RB,RA,QA,L)*E0(K)*EA_WFN(RB,RA,PA,J)
            VEA(QA,PA,I)=VEA(QA,PA,I)-EA_WFN(RB,RA,QA,L)*EA_E(PA,K)*EA_WFN(RB,RA,PA,J)
           ENDDO
          ENDDO
         ENDDO
        ENDIF
! asymmetric
endif
       ENDDO
      ENDDO
     ENDDO
    ENDDO
!   WRITE(6,*) 'VEA MATRIX, ORDER=',I
!   CALL DUMP5(VEA(:,:,I),NCF*EA_NCF)
   ENDDO
   ! CHECK DIAGONALITY OF VEA(:,:,0)
   DO PA=1,NCF*EA_NCF
    DO QA=1,NCF*EA_NCF
     IF ((PA/=QA).AND.(DABS(VEA(PA,QA,0)) > TOL)) CALL PABORT('G0 NONDIAGONAL IN EA')
!    DO I=0,ORDER-1
!     IF (DABS(VEA(PA,QA,I)-VEA(QA,PA,I)) > TOL) THEN
!      WRITE(*,*) VEA(PA,QA,I),VEA(QA,PA,I)
!      CALL PABORT('G0 NONSYMMETRIC IN EA')
!     ENDIF
!    ENDDO
    ENDDO
   ENDDO
   DEALLOCATE(VEC1,VEC2)
   ! COMPARE DIAGONAL VEA VS SMALL SIGMA
!  WRITE(6,'(A,20I10:)') '   ',(I,I,I=0,ORDER-1)
!  DO PA=1,NCF*EA_NCF
!   WRITE(6,'(I3,20F10.6:)') PA,(EA_E(PA,I),VEA(PA,PA,I),I=0,ORDER-1)
!  ENDDO

   ! FORM D
   ALLOCATE(D(0:ORDER-1))
   DO I=0,ORDER-1
    D(I)=0.0D0
    DO J=0,I
     DO PA=1,NCF
      DO PB=1,NCF
       D(I)=D(I)+WFN(PB,PA,J)*WFN(PB,PA,I-J)
      ENDDO
     ENDDO
    ENDDO
!   WRITE(6,*) 'D(',I,')=',D(I)
   ENDDO

   ! FORM LKOOPIP AND LKOOPEA
   ALLOCATE(LKOOPIP(NCF*IP_NCF))
   ALLOCATE(LKOOPEA(NCF*EA_NCF))
   LKOOPIP=.FALSE.
   DO PA=1,NCF
    DO PB=1,IP_NCF
     IF (NORDER(PA)+IP_NORDER(PB) < 2) LKOOPIP((PA-1)*IP_NCF+PB)=.TRUE.
    ENDDO
   ENDDO
   LKOOPEA=.FALSE.
   DO PA=1,NCF
    DO PB=1,EA_NCF
     IF (NORDER(PA)+EA_NORDER(PB) < 2) LKOOPEA((PA-1)*EA_NCF+PB)=.TRUE.
    ENDDO
   ENDDO

   ALLOCATE(GN(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(GIP(NCF*IP_NCF,NCF*IP_NCF,0:ORDER-1))
   ALLOCATE(GEA(NCF*EA_NCF,NCF*EA_NCF,0:ORDER-1))
   GIP=1.0D99 ! safety
   GEA=1.0D99 ! safety

   ! GRAND LOOP
   DO I=0,ORDER-1
   
    ! FORM GIP AND GEA 
    IF (I==0) THEN

     ! IP SECTOR
     GIP(:,:,0)=0.0D0
     DO PA=1,NCF*IP_NCF
!     GIP(PA,PA,0)=1.0D0
      GIP(PA,PA,0)=1.0D0/(OMEGA-IP_E(PA,0))
     ENDDO
     ! EA SECTOR
     GEA(:,:,0)=0.0D0
     DO PA=1,NCF*EA_NCF
!     GEA(PA,PA,0)=1.0D0
      GEA(PA,PA,0)=1.0D0/(OMEGA-EA_E(PA,0))
     ENDDO
!    WRITE(6,*) 'ZEROTH-ORDER GIP'
!    CALL DUMP5(GIP(:,:,0),NCF*IP_NCF)
!    WRITE(6,*) 'ZEROTH-ORDER GEA'
!    CALL DUMP5(GEA(:,:,0),NCF*EA_NCF)

    ELSE

     ! IP SECTOR
     DO N=1,NCF*IP_NCF
      DO M=1,NCF*IP_NCF
       GIP(N,M,I)=0.0D0
       DO J=1,I
        DO L=1,NCF*IP_NCF
         GIP(N,M,I)=GIP(N,M,I)-GIP(N,L,I-J)*VIP(L,M,J)*GIP(M,M,0)
!        GIP(N,M,I)=GIP(N,M,I)-GIP(N,L,I-J)*VIP(M,L,J)
        ENDDO
!       GIP(N,M,I)=GIP(N,M,I)+GIP(N,M,I-J)*IP_E(M,J)*GIP(M,M,0)
       ENDDO
      ENDDO
     ENDDO

     ! EA SECTOR
     DO N=1,NCF*EA_NCF
      DO M=1,NCF*EA_NCF
       GEA(N,M,I)=0.0D0
       DO J=1,I
        DO L=1,NCF*EA_NCF
         GEA(N,M,I)=GEA(N,M,I)-GEA(N,L,I-J)*VEA(L,M,J)*GEA(M,M,0)
!        GEA(N,M,I)=GEA(N,M,I)+GEA(N,L,I-J)*VEA(L,M,J)
        ENDDO
!       GEA(N,M,I)=GEA(N,M,I)+GEA(N,M,I-J)*EA_E(M,J)*GEA(M,M,0)
       ENDDO
      ENDDO
     ENDDO

    ENDIF

!   WRITE(6,*) "GIP AT ORDER=",I
!   CALL DUMP16(GIP(:,:,I),NCF*IP_NCF,NCF*IP_NCF)
!   WRITE(6,*) "GEA AT ORDER=",I
!   CALL DUMP16(GEA(:,:,I),NCF*EA_NCF,NCF*EA_NCF)

    ! FORM GN
    GN(:,:,I)=0.0D0
    DO J=0,I
     DO K=0,I
      DO L=0,I
       IF (J+K+L==I) THEN
        DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
         DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
          ! IP SECTOR
          DO PA=1,NCF*IP_NCF
           DO PB=1,NCF*IP_NCF
            GN(MOX,MOY,I)=GN(MOX,MOY,I)+XIP(PA,MOX,J)*GIP(PA,PB,K)*YIP(PB,MOY,L) 
           ENDDO
          ENDDO
          ! EA SECTOR
          DO PA=1,NCF*EA_NCF
           DO PB=1,NCF*EA_NCF
            GN(MOX,MOY,I)=GN(MOX,MOY,I)+YEA(PA,MOY,J)*GEA(PA,PB,K)*XEA(PB,MOX,L)
           ENDDO
          ENDDO
         ENDDO 
        ENDDO
       ENDIF
      ENDDO
     ENDDO
    ENDDO
    IF (I > 0) THEN
     DO J=1,I
      DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
       DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
        GN(MOX,MOY,I)=GN(MOX,MOY,I)-D(J)*GN(MOX,MOY,I-J)
       ENDDO
      ENDDO
     ENDDO
    ENDIF

    WRITE(6,'(/,A,F20.15,A)') "METHOD 3: TWO-SET EXPANSION"
    WRITE(6,'(I3,A)') I,"TH-ORDER GREEN'S FUNCTION"
    CALL DUMP5(GN(:,:,I),IALL(0,0,0)-IVIRTCORE)

   ENDDO

! ---------------------------------------------------------
!  MBGF versus DeltaMP analysis ... comment out from here
! ---------------------------------------------------------
goto 111
   DO I=0,ORDER-1

    ! FORM GIP AND GEA 
    IF (I==0) THEN

     ! IP SECTOR
     GIP(:,:,0)=0.0D0
     DO PA=1,NCF*IP_NCF
      IF (DABS(OMEGA-IP_E(PA,0)) < 1.0D-10) THEN
       GIP(PA,PA,0)=0.0D0
      ELSE
       GIP(PA,PA,0)=1.0D0/(OMEGA-IP_E(PA,0))
      ENDIF
     ENDDO
     ! EA SECTOR
     GEA(:,:,0)=0.0D0
     DO PA=1,NCF*EA_NCF
      IF (DABS(OMEGA-EA_E(PA,0)) < 1.0D-10) THEN
       GEA(PA,PA,0)=0.0D0
      ELSE
       GEA(PA,PA,0)=1.0D0/(OMEGA-EA_E(PA,0))
      ENDIF
     ENDDO
!    WRITE(6,*) 'ZEROTH-ORDER GIP'
!    CALL DUMP5(GIP(:,:,0),NCF*IP_NCF)
!    WRITE(6,*) 'ZEROTH-ORDER GEA'
!    CALL DUMP5(GEA(:,:,0),NCF*EA_NCF)

    ELSE

     ! IP SECTOR
     DO N=1,NCF*IP_NCF
      DO M=1,NCF*IP_NCF
       GIP(N,M,I)=0.0D0
       DO J=1,I
        DO L=1,NCF*IP_NCF
         GIP(N,M,I)=GIP(N,M,I)-GIP(N,L,I-J)*VIP(L,M,J)*GIP(M,M,0)
!        GIP(N,M,I)=GIP(N,M,I)-GIP(N,L,I-J)*VIP(M,L,J)
        ENDDO
!       GIP(N,M,I)=GIP(N,M,I)+GIP(N,M,I-J)*IP_E(M,J)*GIP(M,M,0)
       ENDDO
      ENDDO
     ENDDO

     ! EA SECTOR
     DO N=1,NCF*EA_NCF
      DO M=1,NCF*EA_NCF
       GEA(N,M,I)=0.0D0
       DO J=1,I
        DO L=1,NCF*EA_NCF
         GEA(N,M,I)=GEA(N,M,I)-GEA(N,L,I-J)*VEA(L,M,J)*GEA(M,M,0)
!        GEA(N,M,I)=GEA(N,M,I)+GEA(N,L,I-J)*VEA(L,M,J)
        ENDDO
!       GEA(N,M,I)=GEA(N,M,I)+GEA(N,M,I-J)*EA_E(M,J)*GEA(M,M,0)
       ENDDO
      ENDDO
     ENDDO

    ENDIF

!   WRITE(6,*) "GIP AT ORDER=",I
!   CALL DUMP16(GIP(:,:,I),NCF*IP_NCF,NCF*IP_NCF)
!   WRITE(6,*) "GEA AT ORDER=",I
!   CALL DUMP16(GEA(:,:,I),NCF*EA_NCF,NCF*EA_NCF)

    ! FORM GN
    GN(:,:,I)=0.0D0
    DO J=0,I
     DO K=0,I
      DO L=0,I
       IF (J+K+L==I) THEN
        DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
         DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
          ! IP SECTOR
          DO PA=1,NCF*IP_NCF
           DO PB=1,NCF*IP_NCF
            GN(MOX,MOY,I)=GN(MOX,MOY,I)+XIP(PA,MOX,J)*GIP(PA,PB,K)*YIP(PB,MOY,L)
           ENDDO
          ENDDO
          ! EA SECTOR
          DO PA=1,NCF*EA_NCF
           DO PB=1,NCF*EA_NCF
            GN(MOX,MOY,I)=GN(MOX,MOY,I)+YEA(PA,MOY,J)*GEA(PA,PB,K)*XEA(PB,MOX,L)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDIF
      ENDDO
     ENDDO
    ENDDO

! overwrite GF(3) !!!
goto 114
write(*,*) "GF3 overwritten!"
if (i==3) then
DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
 DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
  GN(MOX,MOY,3)=0.0D0
  ! X0 G0 V3 G0 Y0
  DO PA=1,NCF*IP_NCF
   DO PB=1,NCF*IP_NCF
    GN(MOX,MOY,3)=GN(MOX,MOY,3)-XIP(PA,MOX,0)*GIP(PA,PA,0)*VIP(PA,PB,3)*GIP(PB,PB,0)*YIP(PB,MOY,0)
   ENDDO
  ENDDO
 ENDDO
ENDDO
endif
114 continue

! overwrite GF(4) !!!
write(*,*) "GF4 overwritten!"
if (i==4) then
DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
 DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
  GN(MOX,MOY,4)=0.0D0
  ! X0 G0 V4 G0 Y0
  DO PA=1,NCF*IP_NCF
!  DO PB=1,NCF*IP_NCF
!   GN(MOX,MOY,4)=GN(MOX,MOY,4)-XIP(PA,MOX,0)*GIP(PA,PA,0)*VIP(PA,PB,4)*GIP(PB,PB,0)*YIP(PB,MOY,0)
    GN(MOX,MOY,4)=GN(MOX,MOY,4)+XIP(PA,MOX,0)*GIP(PA,PA,0)*IP_E(PA,4)*GIP(PA,PA,0)*YIP(PA,MOY,0)
!  ENDDO
  ENDDO
  ! X0 G0 V2 G0 V2 G0 Y0
  DO PA=1,NCF*IP_NCF
!  DO PB=1,NCF*IP_NCF
!   DO QA=1,NCF*IP_NCF
!    GN(MOX,MOY,4)=GN(MOX,MOY,4)+XIP(PA,MOX,0)*GIP(PA,PA,0)*VIP(PA,QA,2)*GIP(QA,QA,0)*VIP(QA,PB,2)*GIP(PB,PB,0)*YIP(PB,MOY,0)
     GN(MOX,MOY,4)=GN(MOX,MOY,4)+XIP(PA,MOX,0)*GIP(PA,PA,0)*IP_E(PA,2)*GIP(PA,PA,0)*IP_E(PA,2)*GIP(PA,PA,0)*YIP(PA,MOY,0)
!   ENDDO
!  ENDDO
  ENDDO
!goto 113
  ! X0 G0 V3 G0 V1 G0 Y0 = zero
  DO PA=1,NCF*IP_NCF
   DO PB=1,NCF*IP_NCF
    DO QA=1,NCF*IP_NCF
     GN(MOX,MOY,4)=GN(MOX,MOY,4)+XIP(PA,MOX,0)*GIP(PA,PA,0)*VIP(PA,QA,3)*GIP(QA,QA,0)*VIP(QA,PB,1)*GIP(PB,PB,0)*YIP(PB,MOY,0)
!    GN(MOX,MOY,4)=GN(MOX,MOY,4)+XIP(PA,MOX,0)*GIP(PA,PA,0)*IP_E(PA,2)*GIP(PA,PA,0)*IP_E(PA,2)*GIP(PA,PA,0)*YIP(PA,MOY,0)
    ENDDO
   ENDDO
  ENDDO
  ! X0 G0 V1 G0 V3 G0 Y0
  DO PA=1,NCF*IP_NCF
   DO PB=1,NCF*IP_NCF
    DO QA=1,NCF*IP_NCF
     GN(MOX,MOY,4)=GN(MOX,MOY,4)+XIP(PA,MOX,0)*GIP(PA,PA,0)*VIP(PA,QA,1)*GIP(QA,QA,0)*VIP(QA,PB,3)*GIP(PB,PB,0)*YIP(PB,MOY,0)
!    GN(MOX,MOY,4)=GN(MOX,MOY,4)+XIP(PA,MOX,0)*GIP(PA,PA,0)*IP_E(PA,2)*GIP(PA,PA,0)*IP_E(PA,2)*GIP(PA,PA,0)*YIP(PA,MOY,0)
    ENDDO
   ENDDO
  ENDDO
  ! X0 G0 V2 G0 V1 G0 V1 G0 Y0
  DO PA=1,NCF*IP_NCF
   DO PB=1,NCF*IP_NCF
    DO QA=1,NCF*IP_NCF
     DO QB=1,NCF*IP_NCF
     GN(MOX,MOY,4)=GN(MOX,MOY,4)-XIP(PA,MOX,0)*GIP(PA,PA,0)*VIP(PA,QA,2)*GIP(QA,QA,0)*VIP(QA,QB,1)*GIP(QB,QB,0) &
                  *VIP(QB,PB,1)*GIP(PB,PB,0)*YIP(PB,MOY,0)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  ! X0 G0 V1 G0 V1 G0 V2 G0 Y0
  DO PA=1,NCF*IP_NCF
   DO PB=1,NCF*IP_NCF
    DO QA=1,NCF*IP_NCF
     DO QB=1,NCF*IP_NCF
     GN(MOX,MOY,4)=GN(MOX,MOY,4)-XIP(PA,MOX,0)*GIP(PA,PA,0)*VIP(PA,QA,1)*GIP(QA,QA,0)*VIP(QA,QB,1)*GIP(QB,QB,0) &
                  *VIP(QB,PB,2)*GIP(PB,PB,0)*YIP(PB,MOY,0)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  ! X0 G0 V1 G0 V2 G0 V1 G0 Y0 = zero
  DO PA=1,NCF*IP_NCF
   DO PB=1,NCF*IP_NCF
    DO QA=1,NCF*IP_NCF
     DO QB=1,NCF*IP_NCF
     GN(MOX,MOY,4)=GN(MOX,MOY,4)-XIP(PA,MOX,0)*GIP(PA,PA,0)*VIP(PA,QA,1)*GIP(QA,QA,0)*VIP(QA,QB,2)*GIP(QB,QB,0) &
                  *VIP(QB,PB,1)*GIP(PB,PB,0)*YIP(PB,MOY,0)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  ! X0 G0 V2 G0 V1 G0 Y1 
  DO PA=1,NCF*IP_NCF
   DO PB=1,NCF*IP_NCF
    DO QA=1,NCF*IP_NCF
     GN(MOX,MOY,4)=GN(MOX,MOY,4)+XIP(PA,MOX,0)*GIP(PA,PA,0)*VIP(PA,QA,2)*GIP(QA,QA,0)*VIP(QA,PB,1)*GIP(PB,PB,0) &
                  *YIP(PB,MOY,1)
    ENDDO
   ENDDO
  ENDDO
  ! X2 G0 V2 G0 Y0 
  DO PA=1,NCF*IP_NCF
   DO PB=1,NCF*IP_NCF
    GN(MOX,MOY,4)=GN(MOX,MOY,4)-XIP(PA,MOX,2)*GIP(PA,PA,0)*VIP(PA,PB,2)*GIP(PB,PB,0)*YIP(PB,MOY,0)
   ENDDO
  ENDDO
113 continue
 ENDDO
ENDDO
endif

    IF (I > 0) THEN
     DO J=1,I
      DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
       DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
        GN(MOX,MOY,I)=GN(MOX,MOY,I)-D(J)*GN(MOX,MOY,I-J)
       ENDDO
      ENDDO
     ENDDO
    ENDIF

    WRITE(6,'(I3,A)') I,"TH-ORDER GREEN'S FUNCTION (STABLE BUT WRONG AT OMEGA=EPSILON !!!)"
    CALL DUMP5(GN(:,:,I),IALL(0,0,0)-IVIRTCORE)

   ENDDO
111 continue

   ALLOCATE(SN(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   ALLOCATE(SN_CUM(IALL(0,0,0)-IVIRTCORE,IALL(0,0,0)-IVIRTCORE,0:ORDER-1))
   SN=1.0D99
   SN_CUM=1.0D99
   DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
    DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
     SN_CUM(MOX,MOY,0)=0.0D0
    ENDDO
    SN_CUM(MOX,MOX,0)=EPSILON(MOX,0,0,0)
   ENDDO
   DO I=0,ORDER-1
    ! FORM SN WITH NUMERICAL REDUCIBLE DELETION
    IF (I==0) THEN
     SN(:,:,I)=1.0D99
    ELSE IF (I==1) THEN
     DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
      DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
!      IF ((DABS(GN(MOX,MOX,0)) > 1.0D-10).AND. &
!          (DABS(GN(MOY,MOY,0)) > 1.0D-10)) THEN
        SN(MOX,MOY,I)=GN(MOX,MOY,1)/GN(MOX,MOX,0)/GN(MOY,MOY,0)
!      ELSE
!       SN(MOX,MOY,I)=0.0D0
!      ENDIF
      ENDDO
     ENDDO
    ELSE
     DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
      DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
       SN(MOX,MOY,I)=GN(MOX,MOY,I)/GN(MOX,MOX,0)/GN(MOY,MOY,0)
!if ((mox==3).and.(moy==3)) then
!write(*,*) SN(3,3,2)*GN(3,3,0)*SN(3,3,2),IP_E(1,2)*GIP(1,1,0)*IP_E(1,2),'    Sigma * G0 * Sigma'
!write(*,*) SN(3,3,2),IP_E(1,2),'    Sigma , Sigmabar'
!write(*,*) GN(3,3,0),GIP(1,1,0),'    G , Fancy G'
!endif
       DO J=1,I-1
!       IF ((I==4).AND.(J==2).AND.(MOX==MOY)) &
!        SN(MOX,MOX,I)=SN(MOX,MOX,I)-SN(MOX,MOX,2)*GN(MOX,MOX,0)*SN(MOX,MOX,2)
        DO MOZ=ICORE+1,IALL(0,0,0)-IVIRTCORE
!       IF ((I==4).AND.(J==2).AND.(MOX==MOY).AND.(MOZ/=MOX)) &
!        SN(MOX,MOX,I)=SN(MOX,MOX,I)-SN(MOX,MOZ,2)*GN(MOZ,MOZ,0)*SN(MOZ,MOX,2)
         DO MOW=ICORE+1,IALL(0,0,0)-IVIRTCORE
! Note: the renormalization term below makes difference only at fourth and higher orders
!         IF ((MOX==3).AND.(MOY==3).AND.(MOZ==3).AND.(MOW==3).AND.(I==4).AND.(J==2)) THEN
!          SN(MOX,MOY,I)=SN(MOX,MOY,I)-IP_E(1,2)*GIP(1,1,0)*IP_E(1,2)
!         IF (MOX==MOZ) THEN
!          SN(MOX,MOY,I)=SN(MOX,MOY,I)-SN(MOZ,MOW,J)*GN(MOW,MOY,I-J)/GN(MOY,MOY,0)
!if ((mox==3).and.(moy==3).and.(i==4).and.(j==2)) write(*,*) sn(moz,mow,j),GN(MOW,MOY,I-J)/GN(MOY,MOY,0)/gn(mow,mow,0),&
!SN(MOZ,MOW,J)*GN(MOW,MOY,I-J)/GN(MOY,MOY,0)
!          ELSE
           SN(MOX,MOY,I)=SN(MOX,MOY,I)-GN(MOX,MOZ,0)*SN(MOZ,MOW,J)*GN(MOW,MOY,I-J)/GN(MOX,MOX,0)/GN(MOY,MOY,0)
!          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDIF
    IF (I > 0) SN_CUM(:,:,I)=SN_CUM(:,:,I-1)+SN(:,:,I)
   ENDDO

   WRITE(6,'(/,A,F20.15,A)') "METHOD 3: TWO-SET EXPANSION"
   DO I=0,ORDER-1
    WRITE(6,'(I3,A)') I,"TH-ORDER SELF-ENERGY"
    CALL DUMP5(SN(:,:,I),IALL(0,0,0)-IVIRTCORE)
   ENDDO

   WRITE(6,'(/,A,F20.15,A)') "METHOD 3: TWO-SET EXPANSION"
   DO I=0,ORDER-1
    WRITE(6,'(I3,A)') I,"TH-ORDER SELF-ENERGY (CUMULATIVE)"
    CALL DUMP5(SN_CUM(:,:,I),IALL(0,0,0)-IVIRTCORE)
   ENDDO

   ALLOCATE(AMAT(IALL(0,0,0)-IVIRTCORE-ICORE,IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(ER(IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(EI(IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(VR(IALL(0,0,0)-IVIRTCORE-ICORE,IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(VL(1,IALL(0,0,0)-IVIRTCORE-ICORE))
   ALLOCATE(WK(4*(IALL(0,0,0)-IVIRTCORE-ICORE)))
   WRITE(6,'(/,A,F20.15,A)') "METHOD 3: TWO-SET EXPANSION"
   DO I=1,ORDER-1
    WRITE(6,'(I3,A)') I,"TH-ORDER SELF-ENERGY (EIGENVALUES)"
    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
     DO MOY=ICORE+1,IALL(0,0,0)-IVIRTCORE
      AMAT(MOX-ICORE,MOY-ICORE)=SN_CUM(MOX,MOY,I)
     ENDDO
    ENDDO
    CALL DGEEV('N','V',IALL(0,0,0)-IVIRTCORE-ICORE,AMAT,IALL(0,0,0)-IVIRTCORE-ICORE,ER,EI,VL,1,VR,&
     IALL(0,0,0)-IVIRTCORE-ICORE,WK,4*(IALL(0,0,0)-IVIRTCORE-ICORE),INFO)
    IF (INFO /= 0) CALL PABORT('DGEEV FAILED TO DIAGONALIZE A MATRIX')
    CALL PIKSRT(IALL(0,0,0)-IVIRTCORE-ICORE,IALL(0,0,0)-IVIRTCORE-ICORE,ER,VR,EI)
    DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
     WRITE(6,'(F20.12)') ER(MOX-ICORE)
    ENDDO
!   WRITE(101,'(I3,F10.5,20F15.8:)') I,OMEGA,(ER(MOX-ICORE),MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE)
   ENDDO
   DEALLOCATE(AMAT,ER,EI,VR,VL,WK)
!  WRITE(6,'(A)') 'fort.101 contains cumulative eigen self-energies'

goto 112
   WRITE(*,*) 'n   MBGF[n](eps)      DeltaMPn'
   WRITE(*,'(I2,2F15.10)') 0,SN_CUM(3,3,0),IP_E(1,0)
   WRITE(*,'(I2,2F15.10)') 1,SN_CUM(3,3,1),IP_E(1,0)+IP_E(1,1)
   WRITE(*,'(I2,2F15.10)') 2,SN_CUM(3,3,2),IP_E(1,0)+IP_E(1,1)+IP_E(1,2)
   WRITE(*,'(I2,2F15.10)') 3,SN_CUM(3,3,3),IP_E(1,0)+IP_E(1,1)+IP_E(1,2)+IP_E(1,3)
   WRITE(*,'(I2,2F15.10)') 4,SN_CUM(3,3,4),IP_E(1,0)+IP_E(1,1)+IP_E(1,2)+IP_E(1,3)+IP_E(1,4)
   tmp1=0.0d0
   DO MOX=ICORE+1,IALL(0,0,0)-IVIRTCORE
    if (mox /= 3) tmp1=tmp1-sn(3,mox,2)*gn(mox,mox,0)*sn(mox,3,2)
   enddo
   WRITE(*,'(I2,F15.10,A)') 4,tmp1,'  Semi-reducible (-Sigma_pq Gqq Sigma_qp)'
   tmp1=0.0d0
   do ra=1,ncf
    do rb=1,ip_ncf
     tmp1=tmp1+ip_wfn(rb,ra,1,0)*ip_wfn(rb,ra,1,2)
    enddo
   enddo
   tmp1=-tmp1*IP_E(1,2)
   WRITE(*,'(I2,F15.10,A)') 4,tmp1,'  - x0 DMP2 <0|2> y0 '
   tmp2=XIP(1,3,2)*IP_E(1,2)
   WRITE(*,'(I2,F15.10,A)') 4,tmp2,'  Linked-disconnected (x2 DMP2 y0)'
   tmp2=0.0d0
   do ra=1,ncf
    do rb=1,ip_ncf
!if (norder(ra)+ip_norder(rb)==3) &
     tmp2=tmp2+ip_wfn(rb,ra,1,1)*ip_wfn(rb,ra,1,1)
    enddo
   enddo
   WRITE(*,'(I2,F15.10,A)') 4,tmp2,'  <N-1,(1)|N-1,(1)> '
   tmp2=tmp2*IP_E(1,2)
   WRITE(*,'(I2,F15.10,A)') 4,tmp2,'  Linked-disconnected (<1|1> DMP2 ; a different type)'
   tmp3=-D(2)*IP_E(1,2)
   WRITE(*,'(I2,F15.10,A)') 4,tmp3,'  - D2 x DMP2'
   tmp4=0.0d0
   DO PA=1,NCF*IP_NCF
    tmp4=tmp4+XIP(PA,3,1)*IP_E(PA,2)*YIP(PA,3,1)
   ENDDO
   WRITE(*,'(I2,F15.10,A)') 4,tmp4,'  x1 DMP2 y1 (IP)'
   tmp5=0.0d0
   DO PA=1,NCF*EA_NCF
!   tmp5=tmp5+YEA(PA,3,1)*EA_E(PA,2)*XEA(PA,3,1)
    tmp5=tmp5+YEA(PA,3,1)*GEA(PA,PA,2)*XEA(PA,3,1)/GN(3,3,0)/GN(3,3,0)
   ENDDO
   WRITE(*,'(I2,F15.10,A)') 4,tmp5,'  x1 DMP2 y1 (EA) should be zero'
   WRITE(*,*) 'XIP:',XIP(1,3,0),XIP(1,3,1),XIP(1,3,2),XIP(1,3,3),' should be 1,0,?,?'
   WRITE(*,*) 'YIP:',YIP(1,3,0),YIP(1,3,1),YIP(1,3,2),YIP(1,3,3),' should be 1,0,0,0'
   allocate(ta1(IP_NCF*NCF,IP_NCF*NCF),ta2(IP_NCF*NCF,IP_NCF*NCF))
   ta1=0.0d0
   do pa=1,ip_ncf*ncf
    do qa=1,ip_ncf*ncf
     do ra=1,ncf
      do rb=1,ip_ncf
       ta1(qa,pa)=ta1(qa,pa)+ip_wfn(rb,ra,qa,0)*ip_wfn(rb,ra,pa,1)
      enddo
     enddo
    enddo
   enddo
   tmp6=0.0d0
   do pa=1,ncf*ip_ncf
    tmp6=tmp6+ta1(1,PA)*ta1(pa,1)
   enddo
   WRITE(*,'(I2,F15.10,A)') 4,tmp6,'  <0|1><0|1>'
   tmp6=tmp6*2.0d0*IP_E(1,2)
   WRITE(*,'(I2,F15.10,A)') 4,tmp6,'  2 <0|1><0|1> DMP2'
   tmp7=0.0d0
   do pa=1,IP_ncf*NCF
   do pb=1,IP_ncf*NCF
    tmp7=tmp7+XIP(PA,3,0)*GIP(PA,PB,3)*YIP(PB,3,1)/GN(3,3,0)/GN(3,3,0)
   enddo
   enddo
   WRITE(*,'(I2,F15.10,A)') 4,tmp7,'  x0 G3 y1'
   tmp7=0.0d0
   do pa=1,ncf*ip_ncf
    tmp7=tmp7-IP_E(1,2)*ta1(1,PA)*YIP(PA,3,1)
   enddo
   WRITE(*,'(I2,F15.10,A)') 4,tmp7,'  -DMP2<0|1>y1'
   WRITE(*,'(I2,F15.10,A)') 2,D(2),'  D(2)'
   WRITE(*,'(I2,F15.10,A)') 2,XIP(1,3,2),'  XIP(gamma,p,2)'
   WRITE(*,'(I2,2F15.10,A)') 1,IP_E(1,1),SN(3,3,1),'  DMP1 SN1'
   WRITE(*,'(I2,2F15.10,A)') 2,IP_E(1,2),SN(3,3,2),'  DMP2 SN2'
   WRITE(*,'(I2,2F15.10,A)') 3,IP_E(1,3),SN(3,3,3),'  DMP3 SN3'
   WRITE(*,'(I2,2F15.10,A)') 4,IP_E(1,4),SN(3,3,4),'  DMP4 SN4'
   WRITE(*,'(I2,2F25.15,A)') 2,IP_E(1,2),SN(3,3,2),'  DMP2 vs SN2'
   WRITE(*,'(I2,2F25.15,A)') 2,IP_E(1,2)*IP_E(1,2)*GN(3,3,0),SN(3,3,2)*SN(3,3,2)*GN(3,3,0),'  DMP2 G DMP2 vs SN2 G SN2'
   deallocate(ta1,ta2)
112 continue

! ---------------------------------------------------------
!  MBGF versus DeltaMP analysis ... up to here
! ---------------------------------------------------------

   DEALLOCATE(IP_WFN,EA_WFN,IP_E,EA_E)
   DEALLOCATE(LKOOPIP,LKOOPEA)
   DEALLOCATE(GN,SN,SN_CUM)
   DEALLOCATE(GIP,GEA)
   DEALLOCATE(XIP,XEA)
   DEALLOCATE(YIP,YEA)
   DEALLOCATE(VIP,VEA)
   DEALLOCATE(WFN,D)
   CLOSE(50)
   CLOSE(51)

   RETURN

END SUBROUTINE
