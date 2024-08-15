SUBROUTINE CLASSIFY_ORBITALS
! CLASSIFY ORBITALS INTO CORE, OCCUPIED AND VIRTUAL ORBITALS.

   USE MPI_F08
   USE CONSTANTS
   USE CONTROL
   USE STRUCTURE
   USE INTEGRAL
   USE BASISSET

   IMPLICIT NONE
!  INCLUDE "mpif.h"
   INTEGER :: I,J,KX,KY,KZ,L

   J=0
   DO I=1,NATOM
    J=J+IATOM(I)
   ENDDO
   IOCC=(J-CHARGE)/2
   ICORE=0
   IF (IOPTN(27) > IOCC) THEN
    CALL PABORT('ILLEGAL NUMBER FOR FROZEN CORE ORBITALS')
   ELSE IF (IOPTN(27) > 0) THEN
    ICORE=IOPTN(27)
   ELSE IF (IOPTN(27) < 0) THEN
    DO I=1,NATOM
     IF ((IATOM(I) >= 3).AND.(IATOM(I) <= 10)) ICORE=ICORE+1
    ENDDO
   ENDIF
   DO KX=-KVCX,MAX(0,KVCX-1)
   DO KY=-KVCY,MAX(0,KVCY-1)
   DO KZ=-KVCZ,MAX(0,KVCZ-1)
    IF (IOPTN(55) > IALL(KX,KY,KZ)-IOCC) CALL PABORT('ILLEGAL NUMBER FOR FROZEN VIRTUAL ORBITALS')
   ENDDO
   ENDDO
   ENDDO
   IF (IOPTN(55) > 0) THEN
    IVIRTCORE=IOPTN(55)
   ELSE IF (IOPTN(55) < 0) THEN
    DO I=1,NATOM
     IF ((IATOM(I) >= 3).AND.(IATOM(I) <= 10)) IVIRTCORE=IVIRTCORE+1
    ENDDO
   ENDIF

   IF (MYID == 0) THEN
   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
   KX=KVCX/11
   IF (KX > 0) THEN
    DO L=0,KX-1
     WRITE(6,'(A8,11I8)') 'KX =    ',(J,J=L*11,L*11+10)
     WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
     WRITE(6,'(A5,3X,11I8)') 'LNDEP',(NCGS-IALL(J,0,0),J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'FCORE',(ICORE,J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'FVIRT',(IVIRTCORE,J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'OCC  ',(IOCC,J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'VIRT ',(IALL(J,0,0)-IOCC,J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'AOCC ',(IOCC-ICORE,J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'AVIRT',(IALL(J,0,0)-IOCC-IVIRTCORE,J=L*11,L*11+10)
     WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
    ENDDO
   ENDIF
   WRITE(6,'(A8,11I8:)') 'KX =    ',(J,J=KX*11,KVCX)
   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
   WRITE(6,'(A5,3X,11I8:)') 'LNDEP',(NCGS-IALL(J,0,0),J=KX*11,KVCX)
   WRITE(6,'(A5,3X,11I8:)') 'FCORE',(ICORE,J=KX*11,KVCX)
   WRITE(6,'(A5,3X,11I8:)') 'FVIRT',(IVIRTCORE,J=KX*11,KVCX)
   WRITE(6,'(A5,3X,11I8:)') 'OCC  ',(IOCC,J=KX*11,KVCX)
   WRITE(6,'(A5,3X,11I8:)') 'VIRT ',(IALL(J,0,0)-IOCC,J=KX*11,KVCX)
   WRITE(6,'(A5,3X,11I8:)') 'AOCC ',(IOCC-ICORE,J=KX*11,KVCX)
   WRITE(6,'(A5,3X,11I8:)') 'AVIRT',(IALL(J,0,0)-IOCC-IVIRTCORE,J=KX*11,KVCX)
   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'

   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
   KY=KVCY/11
   IF (KY > 0) THEN
    DO L=0,KY-1
     WRITE(6,'(A8,11I8)') 'KY =    ',(J,J=L*11,L*11+10)
     WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
     WRITE(6,'(A5,3X,11I8)') 'LNDEP',(NCGS-IALL(0,J,0),J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'FCORE',(ICORE,J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'FVIRT',(IVIRTCORE,J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'OCC  ',(IOCC,J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'VIRT ',(IALL(0,J,0)-IOCC,J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'AOCC ',(IOCC-ICORE,J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'AVIRT',(IALL(0,J,0)-IOCC-IVIRTCORE,J=L*11,L*11+10)
     WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
    ENDDO
   ENDIF
   WRITE(6,'(A8,11I8:)') 'KY =    ',(J,J=KY*11,KVCY)
   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
   WRITE(6,'(A5,3X,11I8:)') 'LNDEP',(NCGS-IALL(0,J,0),J=KY*11,KVCY)
   WRITE(6,'(A5,3X,11I8:)') 'FCORE',(ICORE,J=KY*11,KVCY)
   WRITE(6,'(A5,3X,11I8:)') 'FVIRT',(IVIRTCORE,J=KY*11,KVCY)
   WRITE(6,'(A5,3X,11I8:)') 'OCC  ',(IOCC,J=KY*11,KVCY)
   WRITE(6,'(A5,3X,11I8:)') 'VIRT ',(IALL(0,J,0)-IOCC,J=KY*11,KVCY)
   WRITE(6,'(A5,3X,11I8:)') 'AOCC ',(IOCC-ICORE,J=KY*11,KVCY)
   WRITE(6,'(A5,3X,11I8:)') 'AVIRT',(IALL(0,J,0)-IOCC-IVIRTCORE,J=KY*11,KVCY)
   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'

   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
   KZ=KVCZ/11
   IF (KZ > 0) THEN
    DO L=0,KZ-1
     WRITE(6,'(A8,11I8)') 'KZ =    ',(J,J=L*11,L*11+10)
     WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
     WRITE(6,'(A5,3X,11I8)') 'LNDEP',(NCGS-IALL(0,0,J),J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'FCORE',(ICORE,J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'FVIRT',(IVIRTCORE,J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'OCC  ',(IOCC,J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'VIRT ',(IALL(0,0,J)-IOCC,J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'AOCC ',(IOCC-ICORE,J=L*11,L*11+10)
     WRITE(6,'(A5,3X,11I8)') 'AVIRT',(IALL(0,0,J)-IOCC-IVIRTCORE,J=L*11,L*11+10)
     WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
    ENDDO
   ENDIF
   WRITE(6,'(A8,11I8:)') 'KZ =    ',(J,J=KZ*11,KVCZ)
   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
   WRITE(6,'(A5,3X,11I8:)') 'LNDEP',(NCGS-IALL(0,0,J),J=KZ*11,KVCZ)
   WRITE(6,'(A5,3X,11I8:)') 'FCORE',(ICORE,J=KZ*11,KVCZ)
   WRITE(6,'(A5,3X,11I8:)') 'FVIRT',(IVIRTCORE,J=KZ*11,KVCZ)
   WRITE(6,'(A5,3X,11I8:)') 'OCC  ',(IOCC,J=KZ*11,KVCZ)
   WRITE(6,'(A5,3X,11I8:)') 'VIRT ',(IALL(0,0,J)-IOCC,J=KZ*11,KVCZ)
   WRITE(6,'(A5,3X,11I8:)') 'AOCC ',(IOCC-ICORE,J=KZ*11,KVCZ)
   WRITE(6,'(A5,3X,11I8:)') 'AVIRT',(IALL(0,0,J)-IOCC-IVIRTCORE,J=KZ*11,KVCZ)
   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
   ENDIF

   RETURN
END SUBROUTINE



SUBROUTINE PRINT_ENERGYBANDS(IALG)
! PRINT (QUASI-PARTICLE) ONE-ELECTRON ENERGY LEVELS.

   USE MPI_F08
   USE CONSTANTS
   USE CONTROL
   USE STRUCTURE
   USE INTEGRAL
   USE BASISSET

   IMPLICIT NONE
!  INCLUDE "mpif.h"
   INTEGER :: IALG ! = 1: HF/DFT; =2: MP2; =3: GF2
   INTEGER :: IX,IY,IZ,I,J,KX,KY,KZ,L,M
   DOUBLE PRECISION :: DGAP,IGAP
   DOUBLE PRECISION :: MAXE,MINE,F1,F2,F12,N1,N2,N12
   DOUBLE PRECISION :: TEMP
   INTEGER :: N,NB
   INTEGER :: NMOD
   DOUBLE PRECISION,ALLOCATABLE :: E(:,:,:,:)

   IF ((IALG > 1).AND.(IOPTN(44) <= 0)) RETURN

   ALLOCATE(E(NCGS,-KVCX:KVCX,-KVCY:KVCY,-KVCZ:KVCZ))
   IF (IALG == 1) THEN
    NB=NCGS
    NMOD=1
    E=EPSILON
    IF (DOPTN(98) > 0.0D0) THEN
     TEMP=DOPTN(98)
     IF (MYID == 0) WRITE(6,'(A)') 'HF/DFT ENERGY BANDS AT MP TEMPERATURE'
     IF ((DOPTN(108) > 0.0D0).AND.(MYID == 0)) THEN
      WRITE(6,'(A)') '*********************************************'
      WRITE(6,'(A,E22.15,A)') '* MP TEMPERATURE = ',TEMP,' K *'
      WRITE(6,'(A)') '*********************************************'
     ENDIF
    ELSE
     TEMP=DOPTN(108)
     IF (MYID == 0) WRITE(6,'(A)') 'HF/DFT ENERGY BANDS'
     IF ((DOPTN(108) > 0.0D0).AND.(MYID == 0)) THEN
      WRITE(6,'(A)') '*********************************************'
      WRITE(6,'(A,E22.15,A)') '* HF TEMPERATURE = ',TEMP,' K *'
      WRITE(6,'(A)') '*********************************************'
     ENDIF
    ENDIF
   ELSE IF (IALG == 2) THEN
    NB=IOPTN(44)
    NMOD=IOPTN(99)
    E=QEPSILON
    TEMP=DOPTN(98)
    IF (MYID == 0) WRITE(6,'(A)') 'MP2 QUASI-PARTICLE ENERGY BANDS (DIAGONAL, FREQ-INDEP)'
    IF ((TEMP > 0.0D0).AND.(MYID == 0)) THEN
     WRITE(6,'(A)') '*********************************************'
     WRITE(6,'(A,E22.15,A)') '* MP TEMPERATURE = ',TEMP,' K *'
     WRITE(6,'(A)') '*********************************************'
    ENDIF
   ELSE IF (IALG == 3) THEN
    NB=IOPTN(44)
    NMOD=IOPTN(99)
    E=DEPSILON
    TEMP=DOPTN(98)
    IF ((MYID == 0).AND.(LOPTN(106)))      WRITE(6,'(A)') 'GF2 QUASI-PARTICLE ENERGY BANDS (DIAGONAL, FREQ-DEP)'
    IF ((MYID == 0).AND.(.NOT.LOPTN(106))) WRITE(6,'(A)') 'GF2 QUASI-PARTICLE ENERGY BANDS (NONDIAGONAL, FREQ-DEP)'
    IF ((TEMP > 0.0D0).AND.(MYID == 0)) THEN
     WRITE(6,'(A)') '*********************************************'
     WRITE(6,'(A,E22.15,A)') '* MP TEMPERATURE = ',TEMP,' K *'
     WRITE(6,'(A)') '*********************************************'
    ENDIF
   ENDIF

   IF (TEMP == 0.0D0) THEN
    TEMP=1.0D-4
   ENDIF

!  COUNT THE NUMBER OF ELECTRONS IN A UNIT CELL.  IF THE NUMBER IS ODD, THEN THE PROGRAM PABORTS.
   N=0
   DO I=1,NATOM
    N=N+IATOM(I)
   ENDDO
   N=N-IOPTN(5)
   IF (MOD(N,2) /= 0) CALL PABORT('ODD NUMBER OF ELECTRONS')
   IF (MYID == 0) WRITE(6,'(A,I3)') 'NUMBER OF ELECTRONS = ',N

!  DETERMINE THE FERMI LEVEL
   IF (MYID == 0) WRITE(6,'(A)') 'DETERMINING THE FERMI LEVEL'
   MAXE=-1.0D99
   MINE=1.0D99
   DO IX=-KVCX,MAX(0,KVCX-1)
   DO IY=-KVCY,MAX(0,KVCY-1)
   DO IZ=-KVCZ,MAX(0,KVCZ-1)
    DO J=1,NCGS
     IF (EPSILON(J,IX,IY,IZ) > MAXE) MAXE=EPSILON(J,IX,IY,IZ)
     IF (EPSILON(J,IX,IY,IZ) < MINE) MINE=EPSILON(J,IX,IY,IZ)
    ENDDO
   ENDDO
   ENDDO
   ENDDO
   F1=MINE
   F12=(MINE+MAXE)/2.0D0
   F2=MAXE
   DO WHILE (.TRUE.)
    N1=0.0D0
    DO IX=-KVCX,MAX(0,KVCX-1)
    DO IY=-KVCY,MAX(0,KVCY-1)
    DO IZ=-KVCZ,MAX(0,KVCZ-1)
     DO J=1,NCGS
      N1=N1+2.0D0/(DEXP((EPSILON(J,IX,IY,IZ)-F1)/BOLTZMANN/TEMP)+1.0D0)
     ENDDO
    ENDDO
    ENDDO
    ENDDO
    N1=N1/DFLOAT(MAX(1,2*KVCX)*MAX(1,2*KVCY)*MAX(1,2*KVCZ))-DFLOAT(N)
    N2=0.0D0
    DO IX=-KVCX,MAX(0,KVCX-1)
    DO IY=-KVCY,MAX(0,KVCY-1)
    DO IZ=-KVCZ,MAX(0,KVCZ-1)
     DO J=1,NCGS
      N2=N2+2.0D0/(DEXP((EPSILON(J,IX,IY,IZ)-F2)/BOLTZMANN/TEMP)+1.0D0)
     ENDDO
    ENDDO
    ENDDO
    ENDDO
    N2=N2/DFLOAT(MAX(1,2*KVCX)*MAX(1,2*KVCY)*MAX(1,2*KVCZ))-DFLOAT(N)
    N12=0.0D0
    DO IX=-KVCX,MAX(0,KVCX-1)
    DO IY=-KVCY,MAX(0,KVCY-1)
    DO IZ=-KVCZ,MAX(0,KVCZ-1)
     DO J=1,NCGS
      N12=N12+2.0D0/(DEXP((EPSILON(J,IX,IY,IZ)-F12)/BOLTZMANN/TEMP)+1.0D0)
     ENDDO
    ENDDO
    ENDDO
    ENDDO
    N12=N12/DFLOAT(MAX(1,2*KVCX)*MAX(1,2*KVCY)*MAX(1,2*KVCZ))-DFLOAT(N)

    IF (IOPTN(9) >= 2) WRITE(6,'(3(A,F7.3,A,F10.7))') 'MUMIN=',F1, ': DeltaN=',N1, &
                                                  '    MUMID=',F12,': DeltaN=',N12,&
                                                  '    MUMAX=',F2, ': DeltaN=',N2
    IF (N1 > N12) CALL PABORT('DETERMINATION OF FERMI ENERGY FAILED 2')
    IF (N12 > N2) CALL PABORT('DETERMINATION OF FERMI ENERGY FAILED 3')
    IF (DABS(N12) < 1.0D-10) THEN
     FERMI=F12
     IF (IOPTN(9) >= 2) WRITE(6,'(A,F20.15,A)') 'FERMI ENERGY = ',FERMI
     EXIT
    ELSE IF ((N1 > 0.0D0).AND.(N2 > 0.0D0)) THEN
     F2=F1
     F1=F1-0.1D0
    ELSE IF ((N1 < 0.0D0).AND.(N2 < 0.0D0)) THEN
     F1=F2
     F2=F2+0.1D0
    ELSE IF ((N1 < 0.0D0).AND.(N12 >= 0.0D0)) THEN
     F2=F12
    ELSE IF ((N12 <= 0.0D0).AND.(N2 > 0.0D0)) THEN
     F1=F12
    ELSE
     CALL PABORT('DETERMINATION OF FERMI ENERGY FAILED 4')
    ENDIF
    F12=(F1+F2)/2.0D0

!   IF (F2 > MAXE) CALL PABORT('DETERMINATION OF FERMI ENERGY FAILED')
   ENDDO

   IF (IOPTN(9) == 3) THEN
    DO IX=-KVCX/NMOD,MAX(0,KVCX/NMOD-1)
    DO IY=-KVCY/NMOD,MAX(0,KVCY/NMOD-1)
    DO IZ=-KVCZ/NMOD,MAX(0,KVCZ/NMOD-1)
     IF (MYID == 0) WRITE(6,'(A,100F10.5:)') 'EPSILON = ',(EPSILON(J,IX*NMOD,IY*NMOD,IZ*NMOD),J=1,NCGS)
     IF (MYID == 0) WRITE(6,'(A,100F10.5:)') 'WEIGHT =  ', &
      (2.0D0/(DEXP((EPSILON(J,IX*NMOD,IY*NMOD,IZ*NMOD)-FERMI)/BOLTZMANN/TEMP)+1.0D0),J=1,NCGS)
    ENDDO
    ENDDO
    ENDDO
   ENDIF

   M=0
   DO I=1,NATOM
    M=M+IATOM(I)
   ENDDO
   M=(M-CHARGE)/2
   IF (M >= NCGS) THEN
    IF (MYID == 0) WRITE(6,'(A)') 'UNABLE TO ESTIMATE BAND GAP'
    RETURN
   ENDIF

   IF (MYID == 0) THEN

   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
   KX=KVCX/NMOD/11
   IF (KX > 0) THEN
    DO L=0,KX-1
     WRITE(6,'(A8,I6,10I8)') 'KX =    ',(J*NMOD,J=L*11,L*11+10)
     WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
     DO I=NCGS,1,-1
      IF (I > IOCC+NB) THEN
       CYCLE
      ELSE IF (I < IOCC+1-NB) THEN
       CYCLE
      ELSE IF (I == M) THEN
       IF (TEMP == 0.0D0) THEN
        WRITE(6,'(A5,3X,11F8.5)') ' HOMO',(E(I,-J*NMOD,0,0),J=L*11,L*11+10)
        IF (IALG == 3) &
         WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,-J*NMOD,0,0),J=L*11,L*11+10)
       ELSE
        WRITE(6,'(I5,3X,11F8.5)') I,(E(I,-J*NMOD,0,0),J=L*11,L*11+10)
        WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,-J*NMOD,0,0)-FERMI)/BOLTZMANN/TEMP)+1),J=L*11,L*11+10)
       ENDIF
      ELSE IF (I == M+1) THEN
       IF (TEMP == 0.0D0) THEN
        WRITE(6,'(A5,3X,11F8.5)') ' LUMO',(E(I,-J*NMOD,0,0),J=L*11,L*11+10)
        IF (IALG == 3) &
         WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,-J*NMOD,0,0),J=L*11,L*11+10)
       ELSE
        WRITE(6,'(I5,3X,11F8.5)') I,(E(I,-J*NMOD,0,0),J=L*11,L*11+10)
        WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,-J*NMOD,0,0)-FERMI)/BOLTZMANN/TEMP)+1),J=L*11,L*11+10)
       ENDIF
      ELSE
       IF (TEMP == 0.0D0) THEN
        WRITE(6,'(I5,3X,11F8.5)') I,(E(I,-J*NMOD,0,0),J=L*11,L*11+10)
        IF (IALG == 3) &
         WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,-J*NMOD,0,0),J=L*11,L*11+10)
       ELSE
        WRITE(6,'(I5,3X,11F8.5)') I,(E(I,-J*NMOD,0,0),J=L*11,L*11+10)
        WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,-J*NMOD,0,0)-FERMI)/BOLTZMANN/TEMP)+1),J=L*11,L*11+10)
       ENDIF
      ENDIF
     ENDDO
     WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
    ENDDO
   ENDIF
   WRITE(6,'(A8,I6,10I8:)') 'KX =    ',(J*NMOD,J=KX*11,KVCX/NMOD)
   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
   DO I=NCGS,1,-1
    IF (I > IOCC+NB) THEN
     CYCLE
    ELSE IF (I < IOCC+1-NB) THEN
     CYCLE
    ELSE IF (I == M) THEN
     IF (TEMP == 0.0D0) THEN
      WRITE(6,'(A5,3X,11F8.5)') ' HOMO',(E(I,-J*NMOD,0,0),J=KX*11,KVCX/NMOD)
      IF (IALG == 3) &
       WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,-J*NMOD,0,0),J=KX*11,KVCX/NMOD)
     ELSE
      WRITE(6,'(I5,3X,11F8.5)') I,(E(I,-J*NMOD,0,0),J=KX*11,KVCX/NMOD)
      WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,-J*NMOD,0,0)-FERMI)/BOLTZMANN/TEMP)+1),J=KX*11,KVCX/NMOD)
     ENDIF
    ELSE IF (I == M+1) THEN
     IF (TEMP == 0.0D0) THEN
      WRITE(6,'(A5,3X,11F8.5)') ' LUMO',(E(I,-J*NMOD,0,0),J=KX*11,KVCX/NMOD)
      IF (IALG == 3) &
       WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,-J*NMOD,0,0),J=KX*11,KVCX/NMOD)
     ELSE
      WRITE(6,'(I5,3X,11F8.5)') I,(E(I,-J*NMOD,0,0),J=KX*11,KVCX/NMOD)
      WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,-J*NMOD,0,0)-FERMI)/BOLTZMANN/TEMP)+1),J=KX*11,KVCX/NMOD)
     ENDIF
    ELSE
     IF (TEMP == 0.0D0) THEN
      WRITE(6,'(I5,3X,11F8.5)') I,(E(I,-J*NMOD,0,0),J=KX*11,KVCX/NMOD)
      IF (IALG == 3) &
       WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,-J*NMOD,0,0),J=KX*11,KVCX/NMOD)
     ELSE
      WRITE(6,'(I5,3X,11F8.5)') I,(E(I,-J*NMOD,0,0),J=KX*11,KVCX/NMOD)
      WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,-J*NMOD,0,0)-FERMI)/BOLTZMANN/TEMP)+1),J=KX*11,KVCX/NMOD)
     ENDIF
    ENDIF
   ENDDO
   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'

   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
   KY=KVCY/NMOD/11
   IF (KY > 0) THEN
    DO L=0,KY-1
     WRITE(6,'(A8,I6,10I8)') 'KY =    ',(J*NMOD,J=L*11,L*11+10)
     WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
     DO I=NCGS,1,-1
      IF (I > IOCC+NB) THEN
       CYCLE
      ELSE IF (I < IOCC+1-NB) THEN
       CYCLE
      ELSE IF (I == M) THEN
       IF (TEMP == 0.0D0) THEN
        WRITE(6,'(A5,3X,11F8.5)') ' HOMO',(E(I,0,-J*NMOD,0),J=L*11,L*11+10)
        IF (IALG == 3) &
         WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,0,-J*NMOD,0),J=L*11,L*11+10)
       ELSE
        WRITE(6,'(I5,3X,11F8.5)') I,(E(I,0,-J*NMOD,0),J=L*11,L*11+10)
        WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,0,-J*NMOD,0)-FERMI)/BOLTZMANN/TEMP)+1),J=L*11,L*11+10)
       ENDIF
      ELSE IF (I == M+1) THEN
       IF (TEMP == 0.0D0) THEN
        WRITE(6,'(A5,3X,11F8.5)') ' LUMO',(E(I,0,-J*NMOD,0),J=L*11,L*11+10)
        IF (IALG == 3) &
         WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,0,-J*NMOD,0),J=L*11,L*11+10)
       ELSE
        WRITE(6,'(I5,3X,11F8.5)') I,(E(I,0,-J*NMOD,0),J=L*11,L*11+10)
        WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,0,-J*NMOD,0)-FERMI)/BOLTZMANN/TEMP)+1),J=L*11,L*11+10)
       ENDIF
      ELSE
       IF (TEMP == 0.0D0) THEN
        WRITE(6,'(I5,3X,11F8.5)') I,(E(I,0,-J*NMOD,0),J=L*11,L*11+10)
        IF (IALG == 3) &
         WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,0,-J*NMOD,0),J=L*11,L*11+10)
       ELSE
        WRITE(6,'(I5,3X,11F8.5)') I,(E(I,0,-J*NMOD,0),J=L*11,L*11+10)
        WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,0,-J*NMOD,0)-FERMI)/BOLTZMANN/TEMP)+1),J=L*11,L*11+10)
       ENDIF
      ENDIF
     ENDDO
     WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
    ENDDO
   ENDIF
   WRITE(6,'(A8,I6,10I8:)') 'KY =    ',(J*NMOD,J=KY*11,KVCY/NMOD)
   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
   DO I=NCGS,1,-1
    IF (I > IOCC+NB) THEN
     CYCLE
    ELSE IF (I < IOCC+1-NB) THEN
     CYCLE
    ELSE IF (I == M) THEN
     IF (TEMP == 0.0D0) THEN
      WRITE(6,'(A5,3X,11F8.5)') ' HOMO',(E(I,0,-J*NMOD,0),J=KY*11,KVCY/NMOD)
      IF (IALG == 3) &
       WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,0,-J*NMOD,0),J=KY*11,KVCY/NMOD)
     ELSE
      WRITE(6,'(I5,3X,11F8.5)') I,(E(I,0,-J*NMOD,0),J=KY*11,KVCY/NMOD)
      WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,0,-J*NMOD,0)-FERMI)/BOLTZMANN/TEMP)+1),J=KY*11,KVCY/NMOD)
     ENDIF
    ELSE IF (I == M+1) THEN
     IF (TEMP == 0.0D0) THEN
      WRITE(6,'(A5,3X,11F8.5)') ' LUMO',(E(I,0,-J*NMOD,0),J=KY*11,KVCY/NMOD)
      IF (IALG == 3) &
       WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,0,-J*NMOD,0),J=KY*11,KVCY/NMOD)
     ELSE
      WRITE(6,'(I5,3X,11F8.5)') I,(E(I,0,-J*NMOD,0),J=KY*11,KVCY/NMOD)
      WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,0,-J*NMOD,0)-FERMI)/BOLTZMANN/TEMP)+1),J=KY*11,KVCY/NMOD)
     ENDIF
    ELSE
     IF (TEMP == 0.0D0) THEN
      WRITE(6,'(I5,3X,11F8.5)') I,(E(I,0,-J*NMOD,0),J=KY*11,KVCY/NMOD)
      IF (IALG == 3) &
       WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,0,-J*NMOD,0),J=KY*11,KVCY/NMOD)
     ELSE
      WRITE(6,'(I5,3X,11F8.5)') I,(E(I,0,-J*NMOD,0),J=KY*11,KVCY/NMOD)
      WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,0,-J*NMOD,0)-FERMI)/BOLTZMANN/TEMP)+1),J=KY*11,KVCY/NMOD)
     ENDIF
    ENDIF
   ENDDO
   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'

   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
   KZ=KVCZ/NMOD/11
   IF (KZ > 0) THEN
    DO L=0,KZ-1
     WRITE(6,'(A8,I6,10I8)') 'KZ =    ',(J*NMOD,J=L*11,L*11+10)
     WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
     DO I=NCGS,1,-1
      IF (I > IOCC+NB) THEN
       CYCLE
      ELSE IF (I < IOCC+1-NB) THEN
       CYCLE
      ELSE IF (I == M) THEN
       IF (TEMP == 0.0D0) THEN
        WRITE(6,'(A5,3X,11F8.5)') ' HOMO',(E(I,0,0,-J*NMOD),J=L*11,L*11+10)
        IF (IALG == 3) &
         WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,0,0,-J*NMOD),J=L*11,L*11+10)
       ELSE
        WRITE(6,'(I5,3X,11F8.5)') I,(E(I,0,0,-J*NMOD),J=L*11,L*11+10)
        WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,0,0,-J*NMOD)-FERMI)/BOLTZMANN/TEMP)+1),J=L*11,L*11+10)
       ENDIF
      ELSE IF (I == M+1) THEN
       IF (TEMP == 0.0D0) THEN
        WRITE(6,'(A5,3X,11F8.5)') ' LUMO',(E(I,0,0,-J*NMOD),J=L*11,L*11+10)
        IF (IALG == 3) &
         WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,0,0,-J*NMOD),J=L*11,L*11+10)
       ELSE
        WRITE(6,'(I5,3X,11F8.5)') I,(E(I,0,0,-J*NMOD),J=L*11,L*11+10)
        WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,0,0,-J*NMOD)-FERMI)/BOLTZMANN/TEMP)+1),J=L*11,L*11+10)
       ENDIF
      ELSE
       IF (TEMP == 0.0D0) THEN
        WRITE(6,'(I5,3X,11F8.5)') I,(E(I,0,0,-J*NMOD),J=L*11,L*11+10)
        IF (IALG == 3) &
         WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,0,0,-J*NMOD),J=L*11,L*11+10)
       ELSE
        WRITE(6,'(I5,3X,11F8.5)') I,(E(I,0,0,-J*NMOD),J=L*11,L*11+10)
        WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,0,0,-J*NMOD)-FERMI)/BOLTZMANN/TEMP)+1),J=L*11,L*11+10)
       ENDIF
      ENDIF
     ENDDO
     WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
    ENDDO
   ENDIF
   WRITE(6,'(A8,I6,10I8:)') 'KZ =    ',(J*NMOD,J=KZ*11,KVCZ/NMOD)
   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
   DO I=NCGS,1,-1
    IF (I > IOCC+NB) THEN
     CYCLE
    ELSE IF (I < IOCC+1-NB) THEN
     CYCLE
    ELSE IF (I == M) THEN
     IF (TEMP == 0.0D0) THEN
      WRITE(6,'(A5,3X,11F8.5)') ' HOMO',(E(I,0,0,-J*NMOD),J=KZ*11,KVCZ/NMOD)
      IF (IALG == 3) &
       WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,0,0,-J*NMOD),J=KZ*11,KVCZ/NMOD)
     ELSE
      WRITE(6,'(I5,3X,11F8.5)') I,(E(I,0,0,-J*NMOD),J=KZ*11,KVCZ/NMOD)
      WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,0,0,-J*NMOD)-FERMI)/BOLTZMANN/TEMP)+1),J=KZ*11,KVCZ/NMOD)
     ENDIF
    ELSE IF (I == M+1) THEN
     IF (TEMP == 0.0D0) THEN
      WRITE(6,'(A5,3X,11F8.5)') ' LUMO',(E(I,0,0,-J*NMOD),J=KZ*11,KVCZ/NMOD)
      IF (IALG == 3) &
       WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,0,0,-J*NMOD),J=KZ*11,KVCZ/NMOD)
     ELSE
      WRITE(6,'(I5,3X,11F8.5)') I,(E(I,0,0,-J*NMOD),J=KZ*11,KVCZ/NMOD)
      WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,0,0,-J*NMOD)-FERMI)/BOLTZMANN/TEMP)+1),J=KZ*11,KVCZ/NMOD)
     ENDIF
    ELSE
     IF (TEMP == 0.0D0) THEN
      WRITE(6,'(I5,3X,11F8.5)') I,(E(I,0,0,-J*NMOD),J=KZ*11,KVCZ/NMOD)
      IF (IALG == 3) &
       WRITE(6,'(8X,11(" (",F5.3,")"))') (DPOLE(I,0,0,-J*NMOD),J=KZ*11,KVCZ/NMOD)
     ELSE
      WRITE(6,'(I5,3X,11F8.5)') I,(E(I,0,0,-J*NMOD),J=KZ*11,KVCZ/NMOD)
      WRITE(6,'(8X,11(" (",F5.3,")"))') (2.0D0/(DEXP((E(I,0,0,-J*NMOD)-FERMI)/BOLTZMANN/TEMP)+1),J=KZ*11,KVCZ/NMOD)
     ENDIF
    ENDIF
   ENDDO
   WRITE(6,'(A)') '------------------------------------------------------------------------------------------------'
   ENDIF

   IF (TEMP == 0.0D0) THEN
    DGAP=1.0D99
    DO IX=-KVCX/NMOD,MAX(0,KVCX/NMOD-1)
    DO IY=-KVCY/NMOD,MAX(0,KVCY/NMOD-1)
    DO IZ=-KVCZ/NMOD,MAX(0,KVCZ/NMOD-1)
     IF (E(M+1,IX*NMOD,IY*NMOD,IZ*NMOD)-E(M,IX*NMOD,IY*NMOD,IZ*NMOD) < DGAP) &
      DGAP=E(M+1,IX*NMOD,IY*NMOD,IZ*NMOD)-E(M,IX*NMOD,IY*NMOD,IZ*NMOD)
    ENDDO
    ENDDO
    ENDDO
    IF (MYID == 0) WRITE(6,'(A,F10.7,A,F11.7,A)') 'DIRECT BAND GAP   = ',DGAP,' HARTREE (',DGAP*EV,' EV )'
    IGAP=1.0D99
    DO IX=-KVCX/NMOD,MAX(0,KVCX/NMOD-1)
    DO IY=-KVCY/NMOD,MAX(0,KVCY/NMOD-1)
    DO IZ=-KVCZ/NMOD,MAX(0,KVCZ/NMOD-1)
     DO KX=-KVCX/NMOD,MAX(0,KVCX/NMOD-1)
     DO KY=-KVCY/NMOD,MAX(0,KVCY/NMOD-1)
     DO KZ=-KVCZ/NMOD,MAX(0,KVCZ/NMOD-1)
      IF (E(M+1,IX*NMOD,IY*NMOD,IZ*NMOD)-E(M,KX*NMOD,KY*NMOD,KZ*NMOD) < IGAP) &
       IGAP=E(M+1,IX*NMOD,IY*NMOD,IZ*NMOD)-E(M,KX*NMOD,KY*NMOD,KZ*NMOD)
     ENDDO
     ENDDO
     ENDDO
    ENDDO
    ENDDO
    ENDDO
    IF (MYID == 0) WRITE(6,'(A,F10.7,A,F11.7,A)') 'INDIRECT BAND GAP = ',IGAP,' HARTREE (',IGAP*EV,' EV )'
    IF ((DGAP*EV < 0.1).AND.(MYID == 0)) THEN
     CALL WARNING('********************************')
     CALL WARNING('* CAUTION: SYSTEM MAY BE METAL *')
     CALL WARNING('********************************')
    ELSE IF ((IGAP*EV < 0.1).AND.(MYID == 0)) THEN
     CALL WARNING('*************************************')
     CALL WARNING('* CAUTION: SYSTEM MAY BE SEMI-METAL *')
     CALL WARNING('*************************************')
    ENDIF
   ELSE
    IF (MYID == 0) WRITE(6,'(A,F10.7,A,F11.7,A)') 'FERMI ENERGY      = ',FERMI,' HARTREE (',FERMI*EV,' EV )'
    DGAP=1.0D99
    DO IX=-KVCX/NMOD,MAX(0,KVCX/NMOD-1)
    DO IY=-KVCY/NMOD,MAX(0,KVCY/NMOD-1)
    DO IZ=-KVCZ/NMOD,MAX(0,KVCZ/NMOD-1)
     DO J=1,NCGS
      IF ((E(J,IX*NMOD,IY*NMOD,IZ*NMOD) < FERMI).AND.(FERMI-E(J,IX*NMOD,IY*NMOD,IZ*NMOD) < DGAP)) &
       DGAP=FERMI-E(J,IX*NMOD,IY*NMOD,IZ*NMOD)
     ENDDO
    ENDDO
    ENDDO
    ENDDO
    IGAP=1.0D99
    DO IX=-KVCX/NMOD,MAX(0,KVCX/NMOD-1)
    DO IY=-KVCY/NMOD,MAX(0,KVCY/NMOD-1)
    DO IZ=-KVCZ/NMOD,MAX(0,KVCZ/NMOD-1)
     DO J=1,NCGS
      IF ((E(J,IX*NMOD,IY*NMOD,IZ*NMOD) > FERMI).AND.(E(J,IX*NMOD,IY*NMOD,IZ*NMOD)-FERMI < IGAP)) &
       IGAP=E(J,IX*NMOD,IY*NMOD,IZ*NMOD)-FERMI
     ENDDO
    ENDDO
    ENDDO
    ENDDO
    DGAP=DGAP+IGAP
    IF (MYID == 0) WRITE(6,'(A,F10.7,A,F11.7,A)') 'BAND GAP          = ',DGAP,' HARTREE (',DGAP*EV,' EV )'
    IF (DGAP*EV < 0.1) CALL WARNING('METAL')
   ENDIF

!  DUMP ONE-ELECTRON ENERGIES IN FULL PRECISION
!  COMMENT OUT FOR REGULAR USE FROM HERE ...
!  IF (MYID == 0) THEN
!  DO J=1,NCGS
!   DO IX=-KVCX,MAX(0,KVCX-1)
!   DO IY=-KVCY,MAX(0,KVCY-1)
!   DO IZ=-KVCZ,MAX(0,KVCZ-1)
!    WRITE(77,*) J,IX,IY,IZ,E(J,IX,IY,IZ),E(J,IX,IY,IZ)-EPSILON(J,IX,IY,IZ)
!   ENDDO
!   ENDDO
!   ENDDO
!  ENDDO
!  ENDIF
!  ... TO HERE

   DEALLOCATE(E)
   RETURN
END SUBROUTINE



SUBROUTINE DUMP_CRYSTALORBITALS
! DUMP THE ENTIRE CRYSTALLINE ORBITAL COEFFICIENTS.

   USE MPI_F08
   USE CONSTANTS
   USE CONTROL
   USE STRUCTURE
   USE BASISSET
   USE INTEGRAL
   USE GRADIENT

   IMPLICIT NONE
!  INCLUDE "mpif.h"
   INTEGER :: KX,KY,KZ,MO,I

   IF (MYID == 0) THEN
   DO KX=-KVCX,MAX(0,KVCX-1)
   DO KY=-KVCY,MAX(0,KVCY-1)
   DO KZ=-KVCZ,MAX(0,KVCZ-1)
    DO MO=1,IALL(KX,KY,KZ)
     WRITE(6,'(A,I3,A,3I3,A)') 'CRYSTAL ORBITAL',MO,'[',KX,KY,KZ,']'
     WRITE(6,'(A)') 'CGS       REAL           IMAG'
     DO I=1,NCGS
      WRITE(6,'(I3,2F15.10)') I,CO(I,MO,KX,KY,KZ)
     ENDDO
    ENDDO
   ENDDO
   ENDDO
   ENDDO
   ENDIF

   RETURN
END SUBROUTINE



SUBROUTINE CALC_ENERGYBANDS(IALG)
! INTERPOLATE THE MP2/GF2 ENERGY BANDS ON FINE MESH

   USE MPI_F08
   USE CONSTANTS
   USE CONTROL
   USE STRUCTURE
   USE BASISSET
   USE INTEGRAL
   USE GRADIENT
   USE MULTIPOLE
   USE OEP

   IMPLICIT NONE
!  INCLUDE "mpif.h"
   INTEGER :: IALG ! = 1: HF/DFT; =2: MP2; =3: GF2
   INTEGER :: MESH
   DOUBLE PRECISION,PARAMETER :: MESH2=10000.0D0
   DOUBLE PRECISION,PARAMETER :: THRESH=0.05D0
   DOUBLE PRECISION :: ICPUS,ICPUE
   INTEGER :: QX,QY,QZ,IX,IY,IZ,I,J,K,L,N,JX,JY,JZ,KX,KY,KZ
   INTEGER :: REDUCED_NCGS
   DOUBLE PRECISION :: ASYMMETRY,W,A
   DOUBLE PRECISION,ALLOCATABLE :: E(:,:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: E2(:,:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: E3(:,:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: P(:,:,:,:)
   DOUBLE COMPLEX :: PX,PY,PZ,PQX,PQY,PQZ,AVE
   DOUBLE COMPLEX,ALLOCATABLE :: WS(:,:),WF(:,:),WC(:,:)
   DOUBLE COMPLEX,ALLOCATABLE :: WC_PREVIOUS(:,:)
   DOUBLE COMPLEX,ALLOCATABLE :: WK(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: SIMILARITY(:,:)
   DOUBLE COMPLEX :: OVERLAP
   DOUBLE PRECISION :: MAXIMUM
   LOGICAL,ALLOCATABLE :: LSWAPPED(:)
   INTEGER :: JMAX
   INTEGER :: NSWAPS
   INTEGER :: SWAP(10000,2)
   INTEGER :: NMOD
   DOUBLE PRECISION :: E2IXJYJZ,E2IXKYJZ,E2IXJYKZ,E2IXKYKZ,E2IXIYJZ,E2IXIYKZ
   DOUBLE PRECISION :: E3IXJYJZ,E3IXKYJZ,E3IXJYKZ,E3IXKYKZ,E3IXIYJZ,E3IXIYKZ
   DOUBLE PRECISION :: PIXJYJZ,PIXKYJZ,PIXJYKZ,PIXKYKZ,PIXIYJZ,PIXIYKZ
   DOUBLE PRECISION,ALLOCATABLE :: DOSHF(:,:),DOSMP(:,:),DOSGF(:,:),DOSGFW(:,:)
   INTEGER :: EMIN,EMAX,ES,EE
! -- debug
!  integer :: moa,kax,kay,kaz
!  integer :: mob,kbx,kby,kbz
!  integer :: moi,kix,kiy,kiz
!  integer :: moj,kjx,kjy,kjz
! -- debug end

   MESH=MAX(10,50/MAX(1,MAX(KVCX,MAX(KVCY,KVCZ))))

   EMIN=INT(10.0D0*MESH2)
   EMAX=-INT(10.0D0*MESH2)
   DO IX=-KVCX,MAX(0,KVCX-1)
   DO IY=-KVCY,MAX(0,KVCY-1)
   DO IZ=-KVCZ,MAX(0,KVCZ-1)
    DO I=ICORE+1,IALL(IX,IY,IZ)-IVIRTCORE
     IF (INT(EPSILON(I,IX,IY,IZ)*MESH2) < EMIN) EMIN=INT(EPSILON(I,IX,IY,IZ)*MESH2)
     IF (INT(EPSILON(I,IX,IY,IZ)*MESH2) > EMAX) EMAX=INT(EPSILON(I,IX,IY,IZ)*MESH2)
    ENDDO
   ENDDO
   ENDDO
   ENDDO
   EMIN=EMIN-10
   EMAX=EMAX+10
   IF (EMIN > EMAX) CALL PABORT('AN INTERNAL ERROR')
   ALLOCATE(DOSHF(EMIN:EMAX,NCGS),DOSMP(EMIN:EMAX,NCGS),DOSGF(EMIN:EMAX,NCGS),DOSGFW(EMIN:EMAX,NCGS))

   ALLOCATE(WS(NCGS,NCGS),WF(NCGS,NCGS),WC(NCGS,NCGS),WC_PREVIOUS(NCGS,NCGS),WK(NCGS,NCGS))
   ALLOCATE(E(NCGS,-KVCX*MESH:KVCX*MESH,-KVCY*MESH:KVCY*MESH,-KVCZ*MESH:KVCZ*MESH))
   ALLOCATE(SIMILARITY(NCGS,NCGS),LSWAPPED(NCGS))
   ASYMMETRY=0.0D0

   NMOD=IOPTN(99)
   IF (IALG == 1) THEN
    IF (MYID == 0) WRITE(6,'(A,I0,A)') 'SORTING HF/DFT ENERGY BANDS ON A ',MESH,'x MESH'
    IF (MYID == 0) WRITE(6,'(A)') ' ... KVCX, KVCY, KVCZ, HF ENERGY FOR BAND y IN fort.10y'
    IF (MYID == 0) WRITE(6,'(A)') ' ... ENERGY, HF DOS IN fort.200'
   ELSE IF (IALG == 2) THEN
    IF (MYID == 0) WRITE(6,'(A,I0,A)') 'INTERPOLATING MP2 ENERGY BANDS ON A ',MESH,'x MESH'
    IF (MYID == 0) WRITE(6,'(A)') ' ... KVCX, KVCY, KVCZ, HF ENERGY, MP2 ENERGY FOR BAND y IN fort.10y'
    IF (MYID == 0) WRITE(6,'(A)') ' ... ENERGY, HF DOS IN fort.200'
    IF (MYID == 0) WRITE(6,'(A)') ' ... ENERGY, MP2 DOS IN fort.201'
    ALLOCATE(E2(NCGS,-KVCX*MESH:KVCX*MESH,-KVCY*MESH:KVCY*MESH,-KVCZ*MESH:KVCZ*MESH))
    E2=0.0D0
   ELSE IF (IALG == 3) THEN
    IF (MYID == 0) WRITE(6,'(A,I0,A)') 'INTERPOLATING MP2 & GF2 ENERGY BANDS ON A ',MESH,'x MESH'
    IF (MYID == 0) WRITE(6,'(A)') ' ... KVCX, KVCY, KVCZ, HF ENERGY, MP2 ENERGY, GF2 ENERGY, GF2 POLE STR FOR BAND y IN fort.10y'
    IF (MYID == 0) WRITE(6,'(A)') ' ... ENERGY, HF DOS IN fort.200'
    IF (MYID == 0) WRITE(6,'(A)') ' ... ENERGY, MP2 DOS IN fort.201'
    IF (MYID == 0) WRITE(6,'(A)') ' ... ENERGY, GF2 DOS IN fort.202'
    IF (MYID == 0) WRITE(6,'(A)') ' ... ENERGY, GF2 DOS x POLE STR IN fort.203'
    ALLOCATE(E2(NCGS,-KVCX*MESH:KVCX*MESH,-KVCY*MESH:KVCY*MESH,-KVCZ*MESH:KVCZ*MESH))
    ALLOCATE(E3(NCGS,-KVCX*MESH:KVCX*MESH,-KVCY*MESH:KVCY*MESH,-KVCZ*MESH:KVCZ*MESH))
    ALLOCATE(P(NCGS,-KVCX*MESH:KVCX*MESH,-KVCY*MESH:KVCY*MESH,-KVCZ*MESH:KVCZ*MESH))
    E2=0.0D0
    E3=0.0D0
    P=0.0D0
   ENDIF

!  DUMP FOCK MATRIX
   IF ((IOPTN(9) >= 2).AND.(MYID == 0)) THEN
    WRITE(6,'(A)') 'FOCK MATRIX FOR CONTRACTED GAUSSIANS'
    CALL DUMP1(F_C,NCGS,CEL1X,CEL1Y,CEL1Z)
   ENDIF

!  ZERO SCRATCH DENSITY MATRIX AND ENERGY BAND
   P_C_OUT=0.0D0
   W_C=0.0D0
   E=0.0D0
   SWAP=0
   NSWAPS=0

!  LOOP OVER K VECTORS
   CALL PCPU_TIME(ICPUS)
   DO IX=-KVCX*MESH,KVCX*MESH
   DO IY=-KVCY*MESH,KVCY*MESH
   DO IZ=-KVCZ*MESH,KVCZ*MESH
    IF (IOPTN(9) == 3) WRITE(*,*) 'LOOP 1:',IX,IY,IZ

!   CONSTRUCT DYNAMICAL OVERLAP & FOCK MATRICES
    IF (KVCX == 0) THEN
     PX=DCMPLX(0.0D0,0.0D0)
    ELSE
     PX=DCMPLX(0.0D0,PI*DFLOAT(IX)/DFLOAT(KVCX*MESH))
    ENDIF
    IF (KVCY == 0) THEN
     PY=DCMPLX(0.0D0,0.0D0)
    ELSE
     PY=DCMPLX(0.0D0,PI*DFLOAT(IY)/DFLOAT(KVCY*MESH))
    ENDIF
    IF (KVCZ == 0) THEN
     PZ=DCMPLX(0.0D0,0.0D0)
    ELSE
     PZ=DCMPLX(0.0D0,PI*DFLOAT(IZ)/DFLOAT(KVCZ*MESH))
    ENDIF
    WS=DCMPLX(0.0D0,0.0D0)
    WF=DCMPLX(0.0D0,0.0D0)
    DO J=1,NCGS
     DO K=1,NCGS
      DO QX=-CEL1X,CEL1X
      DO QY=-CEL1Y,CEL1Y
      DO QZ=-CEL1Z,CEL1Z
       PQX=PX*DFLOAT(QX)
       PQY=PY*DFLOAT(QY)
       PQZ=PZ*DFLOAT(QZ)
       WS(K,J)=WS(K,J)+S_C(K,J,QX,QY,QZ)*CDEXP(PQX+PQY+PQZ)
       WF(K,J)=WF(K,J)+F_C(K,J,QX,QY,QZ)*CDEXP(PQX+PQY+PQZ)
      ENDDO
      ENDDO
      ENDDO
     ENDDO
    ENDDO
    DO J=1,NCGS
     DO K=1,NCGS
      IF (CDABS(WF(K,J)-DCONJG(WF(J,K))) > ASYMMETRY) THEN
       ASYMMETRY=CDABS(WF(K,J)-DCONJG(WF(J,K)))
      ENDIF
      AVE=(WF(K,J)+DCONJG(WF(J,K)))/2.0D0
      WF(K,J)=AVE
      WF(J,K)=DCONJG(AVE)
     ENDDO
    ENDDO

!   CONVERT DYNAMICAL FOCK & OVERLAP FROM CARTESIAN TO SPHERICAL GAUSSIANS
    DO J=1,NCGS
     DO K=1,NCGS
      WK(K,J)=DCMPLX(0.0D0,0.0D0)
      DO L=1,NCGS
       WK(K,J)=WK(K,J)+C2S(K,L)*WS(L,J)
      ENDDO
     ENDDO
    ENDDO
    DO J=1,NCGS
     DO K=1,NCGS
      WS(K,J)=DCMPLX(0.0D0,0.0D0)
      DO L=1,NCGS
       WS(K,J)=WS(K,J)+C2S(J,L)*WK(K,L)
      ENDDO
     ENDDO
    ENDDO
    DO J=1,NCGS
     DO K=1,NCGS
      WK(K,J)=DCMPLX(0.0D0,0.0D0)
      DO L=1,NCGS
       WK(K,J)=WK(K,J)+C2S(K,L)*WF(L,J)
      ENDDO
     ENDDO
    ENDDO
    DO J=1,NCGS
     DO K=1,NCGS
      WF(K,J)=DCMPLX(0.0D0,0.0D0)
      DO L=1,NCGS
       WF(K,J)=WF(K,J)+C2S(J,L)*WK(K,L)
      ENDDO
     ENDDO
    ENDDO

!   CANONICAL ORTHOGONALIZATION OF DYNAMICAL OVERLAP MATRIX

    CALL HHBS(NCGS,NCGS,WS,E(:,IX,IY,IZ),WC)

    IF (IOPTN(9) == 3) WRITE(6,'(A,100F10.5:)') 'METRIC  = ',(E(J,IX,IY,IZ),J=1,NCGS)
    REDUCED_NCGS=0
    DO J=1,NCGS
     IF (E(J,IX,IY,IZ) > DOPTN(85)) THEN
      REDUCED_NCGS=REDUCED_NCGS+1
      DO K=1,NCGS
       WS(K,REDUCED_NCGS)=WC(K,J)/DSQRT(E(J,IX,IY,IZ))
      ENDDO
     ENDIF
    ENDDO

!   DIAGONALIZATION OF DYNAMICAL FOCK MATRIX
    DO J=1,REDUCED_NCGS
     DO K=1,NCGS
      WC(K,J)=DCMPLX(0.0D0,0.0D0)
      DO L=1,NCGS
       WC(K,J)=WC(K,J)+WF(K,L)*WS(L,J)
      ENDDO
     ENDDO
    ENDDO
    DO J=1,REDUCED_NCGS
     DO K=1,REDUCED_NCGS
      WF(K,J)=DCMPLX(0.0D0,0.0D0)
      DO L=1,NCGS
       WF(K,J)=WF(K,J)+DCONJG(WS(L,K))*WC(L,J)
      ENDDO
     ENDDO
    ENDDO

    CALL HHBS(NCGS,REDUCED_NCGS,WF,E(:,IX,IY,IZ),WC)

    IF (IALG >= 2) THEN
     IF ((IX/MESH*MESH==IX).AND.(IY/MESH*MESH==IY).AND.(IZ/MESH*MESH==IZ).AND. &
         (IX/MESH/NMOD*NMOD==IX/MESH).AND.(IY/MESH/NMOD*NMOD==IY/MESH).AND.(IZ/MESH/NMOD*NMOD==IZ/MESH)) THEN
      DO J=1,REDUCED_NCGS
       JX=IX/MESH
       IF (JX==KVCX) JX=-KVCX
       JY=IY/MESH
       IF (JY==KVCY) JY=-KVCY
       JZ=IZ/MESH
       IF (JZ==KVCZ) JZ=-KVCZ
       IF (IALG >= 2) E2(J,IX,IY,IZ)=QEPSILON(J,JX,JY,JZ)
       IF (IALG == 3) E3(J,IX,IY,IZ)=DEPSILON(J,JX,JY,JZ)
       IF (IALG == 3) P(J,IX,IY,IZ) =DPOLE(J,JX,JY,JZ)
      ENDDO
     ENDIF
    ENDIF

!   SWAP 
    IF (NSWAPS > 0) THEN
!write(*,*) ix,iy,iz
     DO J=1,NSWAPS
!write(*,*) swap(j,1),swap(j,2),' are swapped'
      MAXIMUM=E(SWAP(J,1),IX,IY,IZ)
      E(SWAP(J,1),IX,IY,IZ)=E(SWAP(J,2),IX,IY,IZ)
      E(SWAP(J,2),IX,IY,IZ)=MAXIMUM
!     IF (LMP2) THEN
!      MAXIMUM=E2(SWAP(J,1),IX,IY,IZ)
!      E2(SWAP(J,1),IX,IY,IZ)=E2(SWAP(J,2),IX,IY,IZ)
!      E2(SWAP(J,2),IX,IY,IZ)=MAXIMUM
!     ENDIF
      DO K=1,NCGS
       OVERLAP=WC(K,SWAP(J,1))
       WC(K,SWAP(J,1))=WC(K,SWAP(J,2))
       WC(K,SWAP(J,2))=OVERLAP
      ENDDO
      IF ((IALG == 1).AND.(IX/MESH*MESH==IX).AND.(IY/MESH*MESH==IY).AND.(IZ/MESH*MESH==IZ)) THEN
       MAXIMUM=EPSILON(SWAP(J,1),IX/MESH,IY/MESH,IZ/MESH)
       EPSILON(SWAP(J,1),IX/MESH,IY/MESH,IZ/MESH)=EPSILON(SWAP(J,2),IX/MESH,IY/MESH,IZ/MESH)
       EPSILON(SWAP(J,2),IX/MESH,IY/MESH,IZ/MESH)=MAXIMUM
       DO K=1,NCGS
        OVERLAP=CO(K,SWAP(J,1),IX/MESH,IY/MESH,IZ/MESH)
        CO(K,SWAP(J,1),IX/MESH,IY/MESH,IZ/MESH)=CO(K,SWAP(J,2),IX/MESH,IY/MESH,IZ/MESH)
        CO(K,SWAP(J,2),IX/MESH,IY/MESH,IZ/MESH)=OVERLAP
       ENDDO
      ENDIF
     ENDDO
    ENDIF

! SORT ENERGY BANDS
    
    IF ((IX /= -KVCX*MESH).OR.(IY /= -KVCY*MESH).OR.(IZ /= -KVCZ*MESH)) THEN
     SIMILARITY=0.0D0
     LSWAPPED=.FALSE.
     DO I=1,REDUCED_NCGS
      DO J=1,REDUCED_NCGS
       OVERLAP=DCMPLX(0.0D0,0.0D0)
       DO K=1,NCGS
!       SIMILARITY(J,I)=SIMILARITY(J,I)+DABS(DREAL(DCONJG(WC_PREVIOUS(K,J))*WC(K,I)))
        OVERLAP=OVERLAP+DABS(DREAL(DCONJG(WC_PREVIOUS(K,J))*WC(K,I)))
!       SIMILARITY(J,I)=SIMILARITY(J,I)+(DCONJG(WC_PREVIOUS(K,J))*WC_PREVIOUS(K,J)*DCONJG(WC(K,I))*WC(K,I))
       ENDDO
       SIMILARITY(J,I)=OVERLAP
      ENDDO
     ENDDO
     DO I=1,REDUCED_NCGS
      MAXIMUM=0.0D0
      DO J=1,REDUCED_NCGS
       IF (SIMILARITY(J,I) > MAXIMUM) THEN
        MAXIMUM=SIMILARITY(J,I)
        JMAX=J
       ENDIF
      ENDDO
!if (myid == 0) write(*,*) I,' maximum:',JMAX
      IF ((I /= JMAX).AND.(DABS(E(I,IX,IY,IZ)-E(JMAX,IX,IY,IZ)) < THRESH).AND.(.NOT.LSWAPPED(I)).AND.(.NOT.LSWAPPED(JMAX))) THEN
       IF (SIMILARITY(JMAX,JMAX) > MAXIMUM) EXIT
       IF (MYID == 0) WRITE(6,'(A,I3,A,I3,A,3I4)') ' ... BANDS',I,' AND ',JMAX,' CROSSED AT K =',IX,IY,IZ
       IF ((IALG >= 2).AND.(I <= IOCC+IOPTN(44)).AND.(JMAX > IOCC+IOPTN(44))) CALL WARNING('INCREASE/DECREASE QPBANDS')
       IF ((IALG >= 2).AND.(JMAX <= IOCC+IOPTN(44)).AND.(I > IOCC+IOPTN(44))) CALL WARNING('INCREASE/DECREASE QPBANDS')
       IF ((IALG >= 2).AND.(I >= IOCC-IOPTN(44)).AND.(JMAX < IOCC-IOPTN(44))) CALL WARNING('INCREASE/DECREASE QPBANDS')
       IF ((IALG >= 2).AND.(JMAX >= IOCC-IOPTN(44)).AND.(I < IOCC-IOPTN(44))) CALL WARNING('INCREASE/DECREASE QPBANDS')
       IF ((I == IOCC).AND.(JMAX == IOCC+1)) CALL PABORT('HOMO AND LUMO CROSSED')
       IF ((JMAX == IOCC).AND.(I == IOCC+1)) CALL PABORT('HOMO AND LUMO CROSSED')
!write(*,*) 'Previous point ',i,jmax,' This point',i,jmax
!do j=1,ncgs
! write(*,'(8F15.8)') wc_previous(j,i),wc_previous(j,jmax),wc(j,i),wc(j,jmax)
!enddo
!overlap=dcmplx(0.0d0,0.0d0)
!do j=1,ncgs
! overlap=overlap+dabs(dreal(dconjg(wc_previous(j,i))*wc(j,i)))
!enddo
!write(*,'(a,2i4,a,2f15.8)') 'overlap',i,i,'=',overlap
!overlap=dcmplx(0.0d0,0.0d0)
!do j=1,ncgs
! overlap=overlap+dabs(dreal(dconjg(wc_previous(j,i))*wc(j,jmax)))
!enddo
!write(*,'(a,2i4,a,2f15.8)') 'overlap',i,jmax,'=',overlap
!overlap=dcmplx(0.0d0,0.0d0)
!do j=1,ncgs
! overlap=overlap+dabs(dreal(dconjg(wc_previous(j,jmax))*wc(j,jmax)))
!enddo
!write(*,'(a,2i4,a,2f15.8)') 'overlap',jmax,jmax,'=',overlap
!if (myid == 0) CALL DUMP5(SIMILARITY,NCGS)
       NSWAPS=NSWAPS+1
       IF (NSWAPS == 10001) CALL PABORT('TOO MANY BAND SWAPS')
       SWAP(NSWAPS,1)=I
       SWAP(NSWAPS,2)=JMAX
       LSWAPPED(I)=.TRUE.
       LSWAPPED(JMAX)=.TRUE.
       MAXIMUM=E(I,IX,IY,IZ)
       E(I,IX,IY,IZ)=E(JMAX,IX,IY,IZ)
       E(JMAX,IX,IY,IZ)=MAXIMUM
!      IF (LMP2) THEN
!       MAXIMUM=E2(I,IX,IY,IZ)
!       E2(I,IX,IY,IZ)=E2(JMAX,IX,IY,IZ)
!       E2(JMAX,IX,IY,IZ)=MAXIMUM
!      ENDIF
       DO K=1,NCGS
        OVERLAP=WC(K,I)
        WC(K,I)=WC(K,JMAX)
        WC(K,JMAX)=OVERLAP
       ENDDO
       IF ((IALG == 1).AND.(IX/MESH*MESH==IX).AND.(IY/MESH*MESH==IY).AND.(IZ/MESH*MESH==IZ)) THEN
        MAXIMUM=EPSILON(I,IX/MESH,IY/MESH,IZ/MESH)
        EPSILON(I,IX/MESH,IY/MESH,IZ/MESH)=EPSILON(JMAX,IX/MESH,IY/MESH,IZ/MESH)
        EPSILON(JMAX,IX/MESH,IY/MESH,IZ/MESH)=MAXIMUM
        DO K=1,NCGS
         OVERLAP=CO(K,I,IX/MESH,IY/MESH,IZ/MESH)
         CO(K,I,IX/MESH,IY/MESH,IZ/MESH)=CO(K,JMAX,IX/MESH,IY/MESH,IZ/MESH)
         CO(K,JMAX,IX/MESH,IY/MESH,IZ/MESH)=OVERLAP
        ENDDO
       ENDIF
      ENDIF
     ENDDO
!if (myid == 0) write(*,*) IX,IY,IZ
!if (myid == 0) CALL DUMP5(SIMILARITY,NCGS)
    ENDIF
    WC_PREVIOUS=WC

! --- PARALLEL LOOP
!   ENDIF
! --- PARALLEL LOOP
   ENDDO
   ENDDO
   IF (IOPTN(9) == 3) WRITE(*,*) 'LOOP 1:END'
   CALL PCPU_TIME(ICPUE)
!  IF (MYID == 0) WRITE(6,'(A,F6.2,A,F10.1,A)') ' ... ', &
!    DFLOAT(IX+KVCX*MESH)/DFLOAT(MAX(1,2*KVCX*MESH+1))*100.0D0, &
!    ' % OF ENERGY BAND CALCULATION (CPU / SEC = ',ICPUE-ICPUS,')'
   CALL PCPU_TIME(ICPUS)
   CALL PFLUSH(6)
   ENDDO

! INTERPOLATE QUASI-PARTICLE ENERGY BANDS
   IF (IALG >= 2) THEN
    DO J=1,REDUCED_NCGS
     DO IX=-KVCX*MESH,KVCX*MESH
     DO IY=-KVCY*MESH,KVCY*MESH
     DO IZ=-KVCZ*MESH,KVCZ*MESH
      IF ((IX/MESH*MESH==IX).AND.(IY/MESH*MESH==IY).AND.(IZ/MESH*MESH==IZ).AND. &
          (IX/MESH/NMOD*NMOD==IX/MESH).AND.(IY/MESH/NMOD*NMOD==IY/MESH).AND.(IZ/MESH/NMOD*NMOD==IZ/MESH)) THEN
       CYCLE
      ELSE
       JX=IABS(IX)/MESH/NMOD*NMOD*ISIGN(1,IX)
       IF ((IX/MESH*MESH==IX).AND.(JX==IX/MESH)) THEN
        KX=JX
       ELSE
        KX=JX+NMOD*ISIGN(1,IX)
       ENDIF
       JY=IABS(IY)/MESH/NMOD*NMOD*ISIGN(1,IY)
       IF ((IY/MESH*MESH==IY).AND.(JY==IY/MESH)) THEN
        KY=JY
       ELSE
        KY=JY+NMOD*ISIGN(1,IY)
       ENDIF
       JZ=IABS(IZ)/MESH/NMOD*NMOD*ISIGN(1,IZ)
       IF ((IZ/MESH*MESH==IZ).AND.(JZ==IZ/MESH)) THEN
        KZ=JZ
       ELSE
        KZ=JZ+NMOD*ISIGN(1,IZ)
       ENDIF
      ENDIF
      JX=JX*MESH
      JY=JY*MESH
      JZ=JZ*MESH
      KX=KX*MESH
      KY=KY*MESH
      KZ=KZ*MESH
! --- INTERPOLATE E2
      IF (KX /= JX) THEN
       E2IXJYJZ=E2(J,JX,JY,JZ)+(E2(J,KX,JY,JZ)-E2(J,JX,JY,JZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX) &
               +E(J,IX,JY,JZ)-E(J,JX,JY,JZ)-(E(J,KX,JY,JZ)-E(J,JX,JY,JZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX)
       E2IXKYJZ=E2(J,JX,KY,JZ)+(E2(J,KX,KY,JZ)-E2(J,JX,KY,JZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX) &
               +E(J,IX,KY,JZ)-E(J,JX,KY,JZ)-(E(J,KX,KY,JZ)-E(J,JX,KY,JZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX)
       E2IXJYKZ=E2(J,JX,JY,KZ)+(E2(J,KX,JY,KZ)-E2(J,JX,JY,KZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX) &
               +E(J,IX,JY,KZ)-E(J,JX,JY,KZ)-(E(J,KX,JY,KZ)-E(J,JX,JY,KZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX)
       E2IXKYKZ=E2(J,JX,KY,KZ)+(E2(J,KX,KY,KZ)-E2(J,JX,KY,KZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX) &
               +E(J,IX,KY,KZ)-E(J,JX,KY,KZ)-(E(J,KX,KY,KZ)-E(J,JX,KY,KZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX)
      ELSE
       E2IXJYJZ=E2(J,JX,JY,JZ)
       E2IXKYJZ=E2(J,JX,KY,JZ)
       E2IXJYKZ=E2(J,JX,JY,KZ)
       E2IXKYKZ=E2(J,JX,KY,KZ)
      ENDIF
      IF (KY /= JY) THEN 
       E2IXIYJZ=E2IXJYJZ+(E2IXKYJZ-E2IXJYJZ)*DFLOAT(IY-JY)/DFLOAT(KY-JY) &
               +E(J,IX,IY,JZ)-E(J,IX,JY,JZ)-(E(J,IX,KY,JZ)-E(J,IX,JY,JZ))*DFLOAT(IY-JY)/DFLOAT(KY-JY)
       E2IXIYKZ=E2IXJYKZ+(E2IXKYKZ-E2IXJYKZ)*DFLOAT(IY-JY)/DFLOAT(KY-JY) &
               +E(J,IX,IY,KZ)-E(J,IX,JY,KZ)-(E(J,IX,KY,KZ)-E(J,IX,JY,KZ))*DFLOAT(IY-JY)/DFLOAT(KY-JY)
      ELSE
       E2IXIYJZ=E2IXJYJZ
       E2IXIYKZ=E2IXJYKZ
      ENDIF
      IF (KZ /= JZ) THEN 
       E2(J,IX,IY,IZ)=E2(J,IX,JY,JZ)+(E2IXIYJZ-E2IXIYKZ)*DFLOAT(IZ-JZ)/DFLOAT(KZ-JZ) &
               +E(J,IX,IY,IZ)-E(J,IX,IY,JZ)-(E(J,IX,IY,KZ)-E(J,IX,IY,JZ))*DFLOAT(IZ-JZ)/DFLOAT(KZ-JZ)
      ELSE
       E2(J,IX,IY,JZ)=E2IXJYJZ
      ENDIF
      IF (IALG == 3) THEN
! ---- INTERPOLATE E3
       IF (KX /= JX) THEN
        E3IXJYJZ=E3(J,JX,JY,JZ)+(E3(J,KX,JY,JZ)-E3(J,JX,JY,JZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX) &
                +E(J,IX,JY,JZ)-E(J,JX,JY,JZ)-(E(J,KX,JY,JZ)-E(J,JX,JY,JZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX)
        E3IXKYJZ=E3(J,JX,KY,JZ)+(E3(J,KX,KY,JZ)-E3(J,JX,KY,JZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX) &
                +E(J,IX,KY,JZ)-E(J,JX,KY,JZ)-(E(J,KX,KY,JZ)-E(J,JX,KY,JZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX)
        E3IXJYKZ=E3(J,JX,JY,KZ)+(E3(J,KX,JY,KZ)-E3(J,JX,JY,KZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX) &
                +E(J,IX,JY,KZ)-E(J,JX,JY,KZ)-(E(J,KX,JY,KZ)-E(J,JX,JY,KZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX)
        E3IXKYKZ=E3(J,JX,KY,KZ)+(E3(J,KX,KY,KZ)-E3(J,JX,KY,KZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX) &
                +E(J,IX,KY,KZ)-E(J,JX,KY,KZ)-(E(J,KX,KY,KZ)-E(J,JX,KY,KZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX)
       ELSE
        E3IXJYJZ=E3(J,JX,JY,JZ)
        E3IXKYJZ=E3(J,JX,KY,JZ)
        E3IXJYKZ=E3(J,JX,JY,KZ)
        E3IXKYKZ=E3(J,JX,KY,KZ)
       ENDIF
       IF (KY /= JY) THEN
        E3IXIYJZ=E3IXJYJZ+(E3IXKYJZ-E3IXJYJZ)*DFLOAT(IY-JY)/DFLOAT(KY-JY) &
                +E(J,IX,IY,JZ)-E(J,IX,JY,JZ)-(E(J,IX,KY,JZ)-E(J,IX,JY,JZ))*DFLOAT(IY-JY)/DFLOAT(KY-JY)
        E3IXIYKZ=E3IXJYKZ+(E3IXKYKZ-E3IXJYKZ)*DFLOAT(IY-JY)/DFLOAT(KY-JY) &
                +E(J,IX,IY,KZ)-E(J,IX,JY,KZ)-(E(J,IX,KY,KZ)-E(J,IX,JY,KZ))*DFLOAT(IY-JY)/DFLOAT(KY-JY)
       ELSE
        E3IXIYJZ=E3IXJYJZ
        E3IXIYKZ=E3IXJYKZ
       ENDIF
       IF (KZ /= JZ) THEN
        E3(J,IX,IY,IZ)=E3(J,IX,JY,JZ)+(E3IXIYJZ-E3IXIYKZ)*DFLOAT(IZ-JZ)/DFLOAT(KZ-JZ) &
                +E(J,IX,IY,IZ)-E(J,IX,IY,JZ)-(E(J,IX,IY,KZ)-E(J,IX,IY,JZ))*DFLOAT(IZ-JZ)/DFLOAT(KZ-JZ)
       ELSE
        E3(J,IX,IY,JZ)=E3IXJYJZ
       ENDIF
! ---- INTERPOLATE POLE STRENGTH
       IF (KX /= JX) THEN
        PIXJYJZ=P(J,JX,JY,JZ)+(P(J,KX,JY,JZ)-P(J,JX,JY,JZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX) 
        PIXKYJZ=P(J,JX,KY,JZ)+(P(J,KX,KY,JZ)-P(J,JX,KY,JZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX)
        PIXJYKZ=P(J,JX,JY,KZ)+(P(J,KX,JY,KZ)-P(J,JX,JY,KZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX)
        PIXKYKZ=P(J,JX,KY,KZ)+(P(J,KX,KY,KZ)-P(J,JX,KY,KZ))*DFLOAT(IX-JX)/DFLOAT(KX-JX)
       ELSE
        PIXJYJZ=P(J,JX,JY,JZ)
        PIXKYJZ=P(J,JX,KY,JZ)
        PIXJYKZ=P(J,JX,JY,KZ)
        PIXKYKZ=P(J,JX,KY,KZ)
       ENDIF
       IF (KY /= JY) THEN
        PIXIYJZ=PIXJYJZ+(PIXKYJZ-PIXJYJZ)*DFLOAT(IY-JY)/DFLOAT(KY-JY)
        PIXIYKZ=PIXJYKZ+(PIXKYKZ-PIXJYKZ)*DFLOAT(IY-JY)/DFLOAT(KY-JY)
       ELSE
        PIXIYJZ=PIXJYJZ
        PIXIYKZ=PIXJYKZ
       ENDIF
       IF (KZ /= JZ) THEN
        P(J,IX,IY,IZ)=P(J,IX,JY,JZ)+(PIXIYJZ-PIXIYKZ)*DFLOAT(IZ-JZ)/DFLOAT(KZ-JZ)
       ELSE
        P(J,IX,IY,JZ)=PIXJYJZ
       ENDIF
      ENDIF
     ENDDO
     ENDDO
     ENDDO
    ENDDO
   ENDIF

! ---- WRITE ENERGY BANDS TO FILE fort.100+
   IF (MYID == 0) THEN
    DO J=1,REDUCED_NCGS
     REWIND(100+J) ! OVERWRITE fort.100+
     DO IX=-KVCX*MESH,KVCX*MESH
     DO IY=-KVCY*MESH,KVCY*MESH
     DO IZ=-KVCZ*MESH,KVCZ*MESH
      IF (IOPTN(9) == 3) WRITE(*,*) 'LOOP 2:',IX,IY,IZ
      IF ((IALG == 2).AND.(J >= IOCC+1-IOPTN(44)).AND.(J <= IOCC+IOPTN(44))) THEN
       WRITE(100+J,'(3I5,2F25.15)') IX,IY,IZ,E(J,IX,IY,IZ),E2(J,IX,IY,IZ)
      ELSE IF ((IALG == 3).AND.(J >= IOCC+1-IOPTN(44)).AND.(J <= IOCC+IOPTN(44))) THEN
       WRITE(100+J,'(3I5,3F25.15,F10.5)') IX,IY,IZ,E(J,IX,IY,IZ),E2(J,IX,IY,IZ),E3(J,IX,IY,IZ),P(J,IX,IY,IZ)
      ELSE
       WRITE(100+J,'(3I5,F25.15)') IX,IY,IZ,E(J,IX,IY,IZ)
      ENDIF
     ENDDO
     ENDDO
     ENDDO
    ENDDO
   ENDIF

! ---- ACCUMULATE DOS
   DOSHF=0.0D0
   DOSMP=0.0D0
   DOSGF=0.0D0
   DOSGFW=0.0D0
   DO J=1,REDUCED_NCGS
    IF (IOPTN(9) == 3) WRITE(*,*) 'LOOP 3:',J
    DO IX=-KVCX*MESH,KVCX*MESH
    DO IY=-KVCY*MESH,KVCY*MESH
    DO IZ=-KVCZ*MESH,KVCZ*MESH
! ---- HF DOS
     IF (IX < KVCX*MESH) THEN
      ES=INT(E(J,IX,IY,IZ)*MESH2)
      EE=INT(E(J,IX+1,IY,IZ)*MESH2)
      IF (ES < EMIN+10) CYCLE
      IF (EE < EMAX-10) CYCLE
      IF (IOPTN(9) == 3) WRITE(*,*) 'ES,EE,EMIN,EMAX:',ES,EE,EMIN,EMAX
      IF (EE > ES) THEN
       DO K=ES,EE-1
        DOSHF(K,J)=DOSHF(K,J)+1.0D0/DABS(E(J,IX+1,IY,IZ)-E(J,IX,IY,IZ))
       ENDDO
      ELSE IF (EE < ES) THEN
       DO K=ES,EE+1,-1
        DOSHF(K,J)=DOSHF(K,J)+1.0D0/DABS(E(J,IX+1,IY,IZ)-E(J,IX,IY,IZ))
       ENDDO
      ELSE
       DOSHF(ES,J)=DOSHF(ES,J)+1.0D0/DABS(E(J,IX+1,IY,IZ)-E(J,IX,IY,IZ))
      ENDIF
     ENDIF
     IF (IY < KVCY*MESH) THEN
      ES=INT(E(J,IX,IY,IZ)*MESH2)
      EE=INT(E(J,IX,IY+1,IZ)*MESH2)
      IF (ES < EMIN+10) CYCLE
      IF (EE < EMAX-10) CYCLE
      IF (EE > ES) THEN
       DO K=ES,EE-1
        DOSHF(K,J)=DOSHF(K,J)+1.0D0/DABS(E(J,IX,IY,IZ)-E(J,IX,IY+1,IZ))
       ENDDO
      ELSE IF (EE < ES) THEN
       DO K=ES,EE+1,-1
        DOSHF(K,J)=DOSHF(K,J)+1.0D0/DABS(E(J,IX,IY,IZ)-E(J,IX,IY+1,IZ))
       ENDDO
      ELSE
       DOSHF(ES,J)=DOSHF(ES,J)+1.0D0/DABS(E(J,IX,IY,IZ)-E(J,IX,IY+1,IZ))
      ENDIF
     ENDIF
     IF (IZ < KVCZ*MESH) THEN
      ES=INT(E(J,IX,IY,IZ)*MESH2)
      EE=INT(E(J,IX,IY,IZ+1)*MESH2)
      IF (ES < EMIN+10) CYCLE
      IF (EE < EMAX-10) CYCLE
      IF (EE > ES) THEN
       DO K=ES,EE-1
        DOSHF(K,J)=DOSHF(K,J)+1.0D0/DABS(E(J,IX,IY,IZ)-E(J,IX,IY,IZ+1))
       ENDDO
      ELSE IF (EE < ES) THEN
       DO K=ES,EE+1,-1
        DOSHF(K,J)=DOSHF(K,J)+1.0D0/DABS(E(J,IX,IY,IZ)-E(J,IX,IY,IZ+1))
       ENDDO
      ELSE
       DOSHF(ES,J)=DOSHF(ES,J)+1.0D0/DABS(E(J,IX,IY,IZ)-E(J,IX,IY,IZ+1))
      ENDIF
     ENDIF
! ---- MP2 DOS
     IF ((IALG >= 2).AND.(J >= IOCC+1-IOPTN(44)).AND.(J <= IOCC+IOPTN(44))) THEN
     IF (IX < KVCX*MESH) THEN
      ES=INT(E2(J,IX,IY,IZ)*MESH2)
      EE=INT(E2(J,IX+1,IY,IZ)*MESH2)
      IF (ES < EMIN+10) CYCLE
      IF (EE < EMAX-10) CYCLE
      IF (EE > ES) THEN
       DO K=ES,EE-1
        DOSMP(K,J)=DOSMP(K,J)+1.0D0/DABS(E2(J,IX+1,IY,IZ)-E2(J,IX,IY,IZ))
       ENDDO
      ELSE IF (EE < ES) THEN
       DO K=ES,EE+1,-1
        DOSMP(K,J)=DOSMP(K,J)+1.0D0/DABS(E2(J,IX+1,IY,IZ)-E2(J,IX,IY,IZ))
       ENDDO
      ELSE
       DOSMP(ES,J)=DOSMP(ES,J)+1.0D0/DABS(E2(J,IX+1,IY,IZ)-E2(J,IX,IY,IZ))
      ENDIF
     ENDIF
     IF (IY < KVCY*MESH) THEN
      ES=INT(E2(J,IX,IY,IZ)*MESH2)
      EE=INT(E2(J,IX,IY+1,IZ)*MESH2)
      IF (ES < EMIN+10) CYCLE
      IF (EE < EMAX-10) CYCLE
      IF (EE > ES) THEN
       DO K=ES,EE-1
        DOSMP(K,J)=DOSMP(K,J)+1.0D0/DABS(E2(J,IX,IY,IZ)-E2(J,IX,IY+1,IZ))
       ENDDO
      ELSE IF (EE < ES) THEN
       DO K=ES,EE+1,-1
        DOSMP(K,J)=DOSMP(K,J)+1.0D0/DABS(E2(J,IX,IY,IZ)-E2(J,IX,IY+1,IZ))
       ENDDO
      ELSE
       DOSMP(ES,J)=DOSMP(ES,J)+1.0D0/DABS(E2(J,IX,IY,IZ)-E2(J,IX,IY+1,IZ))
      ENDIF
     ENDIF
     IF (IZ < KVCZ*MESH) THEN
      ES=INT(E2(J,IX,IY,IZ)*MESH2)
      EE=INT(E2(J,IX,IY,IZ+1)*MESH2)
      IF (ES < EMIN+10) CYCLE
      IF (EE < EMAX-10) CYCLE
      IF (EE > ES) THEN
       DO K=ES,EE-1
        DOSMP(K,J)=DOSMP(K,J)+1.0D0/DABS(E2(J,IX,IY,IZ)-E2(J,IX,IY,IZ+1))
       ENDDO
      ELSE IF (EE < ES) THEN
       DO K=ES,EE+1,-1
        DOSMP(K,J)=DOSMP(K,J)+1.0D0/DABS(E2(J,IX,IY,IZ)-E2(J,IX,IY,IZ+1))
       ENDDO
      ELSE
       DOSMP(ES,J)=DOSMP(ES,J)+1.0D0/DABS(E2(J,IX,IY,IZ)-E2(J,IX,IY,IZ+1))
      ENDIF
     ENDIF
     ENDIF
! ---- GF2 DOS
     IF ((IALG == 3).AND.(J >= IOCC+1-IOPTN(44)).AND.(J <= IOCC+IOPTN(44))) THEN
     IF (IX < KVCX*MESH) THEN
      ES=INT(E3(J,IX,IY,IZ)*MESH2)
      EE=INT(E3(J,IX+1,IY,IZ)*MESH2)
      IF (ES < EMIN+10) CYCLE
      IF (EE < EMAX-10) CYCLE
      IF (EE > ES) THEN
       DO K=ES,EE-1
        DOSGF(K,J)=DOSGF(K,J)+1.0D0/DABS(E3(J,IX+1,IY,IZ)-E3(J,IX,IY,IZ))
        DOSGFW(K,J)=DOSGFW(K,J)+1.0D0/DABS(E3(J,IX+1,IY,IZ)-E3(J,IX,IY,IZ))*P(J,IX,IY,IZ)
       ENDDO
      ELSE IF (EE < ES) THEN
       DO K=ES,EE+1,-1
        DOSGF(K,J)=DOSGF(K,J)+1.0D0/DABS(E3(J,IX+1,IY,IZ)-E3(J,IX,IY,IZ))
        DOSGFW(K,J)=DOSGFW(K,J)+1.0D0/DABS(E3(J,IX+1,IY,IZ)-E3(J,IX,IY,IZ))*P(J,IX,IY,IZ)
       ENDDO
      ELSE
       DOSGF(ES,J)=DOSGF(ES,J)+1.0D0/DABS(E3(J,IX+1,IY,IZ)-E3(J,IX,IY,IZ))
       DOSGFW(ES,J)=DOSGFW(ES,J)+1.0D0/DABS(E3(J,IX+1,IY,IZ)-E3(J,IX,IY,IZ))*P(J,IX,IY,IZ)
      ENDIF
     ENDIF
     IF (IY < KVCY*MESH) THEN
      ES=INT(E3(J,IX,IY,IZ)*MESH2)
      EE=INT(E3(J,IX,IY+1,IZ)*MESH2)
      IF (ES < EMIN+10) CYCLE
      IF (EE < EMAX-10) CYCLE
      IF (EE > ES) THEN
       DO K=ES,EE-1
        DOSGF(K,J)=DOSGF(K,J)+1.0D0/DABS(E3(J,IX,IY,IZ)-E3(J,IX,IY+1,IZ))
        DOSGFW(K,J)=DOSGFW(K,J)+1.0D0/DABS(E3(J,IX,IY,IZ)-E3(J,IX,IY+1,IZ))*P(J,IX,IY,IZ)
       ENDDO
      ELSE IF (EE < ES) THEN
       DO K=ES,EE+1,-1
        DOSGF(K,J)=DOSGF(K,J)+1.0D0/DABS(E3(J,IX,IY,IZ)-E3(J,IX,IY+1,IZ))
        DOSGFW(K,J)=DOSGFW(K,J)+1.0D0/DABS(E3(J,IX,IY,IZ)-E3(J,IX,IY+1,IZ))*P(J,IX,IY,IZ)
       ENDDO
      ELSE
       DOSGF(ES,J)=DOSGF(ES,J)+1.0D0/DABS(E3(J,IX,IY,IZ)-E3(J,IX,IY+1,IZ))
       DOSGFW(ES,J)=DOSGFW(ES,J)+1.0D0/DABS(E3(J,IX,IY,IZ)-E3(J,IX,IY+1,IZ))*P(J,IX,IY,IZ)
      ENDIF
     ENDIF
     IF (IZ < KVCZ*MESH) THEN
      ES=INT(E3(J,IX,IY,IZ)*MESH2)
      EE=INT(E3(J,IX,IY,IZ+1)*MESH2)
      IF (ES < EMIN+10) CYCLE
      IF (EE < EMAX-10) CYCLE
      IF (EE > ES) THEN
       DO K=ES,EE-1
        DOSGF(K,J)=DOSGF(K,J)+1.0D0/DABS(E3(J,IX,IY,IZ)-E3(J,IX,IY,IZ+1))
        DOSGFW(K,J)=DOSGFW(K,J)+1.0D0/DABS(E3(J,IX,IY,IZ)-E3(J,IX,IY,IZ+1))*P(J,IX,IY,IZ)
       ENDDO
      ELSE IF (EE < ES) THEN
       DO K=ES,EE+1,-1
        DOSGF(K,J)=DOSGF(K,J)+1.0D0/DABS(E3(J,IX,IY,IZ)-E3(J,IX,IY,IZ+1))
        DOSGFW(K,J)=DOSGFW(K,J)+1.0D0/DABS(E3(J,IX,IY,IZ)-E3(J,IX,IY,IZ+1))*P(J,IX,IY,IZ)
       ENDDO
      ELSE
       DOSGF(ES,J)=DOSGF(ES,J)+1.0D0/DABS(E3(J,IX,IY,IZ)-E3(J,IX,IY,IZ+1))
       DOSGFW(ES,J)=DOSGFW(ES,J)+1.0D0/DABS(E3(J,IX,IY,IZ)-E3(J,IX,IY,IZ+1))*P(J,IX,IY,IZ)
      ENDIF
     ENDIF
     ENDIF
    ENDDO
    ENDDO
    ENDDO
   ENDDO
   
! ---- WRITE DOS TO FILE fort.200
   IF (MYID == 0) THEN
    IF (IALG == 1) THEN
     REWIND(200)
     DO I=EMIN,EMAX
      WRITE(200,'(F9.4,100E11.3)') DFLOAT(I)/MESH2,(DOSHF(I,J),J=ICORE+1,IALLMAX-IVIRTCORE)
     ENDDO
    ELSE IF (IALG == 2) THEN
     REWIND(201)
     DO I=EMIN,EMAX
      WRITE(201,'(F9.4,100E11.3)') DFLOAT(I)/MESH2,(DOSMP(I,J),J=MAX(ICORE+1,IOCC+1-IOPTN(44)),MIN(IALLMAX,IOCC+IOPTN(44)))
     ENDDO
    ELSE IF (IALG == 3) THEN
     REWIND(201)
     REWIND(202)
     REWIND(203)
     DO I=EMIN,EMAX
      WRITE(201,'(F9.4,100E11.3)') DFLOAT(I)/MESH2,(DOSMP(I,J),J=MAX(ICORE+1,IOCC+1-IOPTN(44)),MIN(IALLMAX,IOCC+IOPTN(44)))
      WRITE(202,'(F9.4,100E11.3)') DFLOAT(I)/MESH2,(DOSGF(I,J),J=MAX(ICORE+1,IOCC+1-IOPTN(44)),MIN(IALLMAX,IOCC+IOPTN(44)))
      WRITE(203,'(F9.4,100E11.3)') DFLOAT(I)/MESH2,(DOSGFW(I,J),J=MAX(ICORE+1,IOCC+1-IOPTN(44)),MIN(IALLMAX,IOCC+IOPTN(44)))
     ENDDO
    ENDIF
   ENDIF

! -- debug satellite
!  do kix=-kvcx,kvcx
!  do kiy=-kvcy,kvcy
!  do kiz=-kvcz,kvcz
!  do moi=iocc-1,iocc
!  do kjx=-kvcx,kvcx
!  do kjy=-kvcy,kvcy
!  do kjz=-kvcz,kvcz
!  do moj=iocc-1,iocc
!  do kax=-kvcx,kvcx
!  do kay=-kvcy,kvcy
!  do kaz=-kvcz,kvcz
!  do moa=iocc+1,iocc+2
!  kbx=kix+kjx-kax
!  kby=kiy+kjy-kay
!  kbz=kiz+kjz-kaz
!  if (kbx < -kvcx) kbx=kbx+2*kvcx
!  if (kbx > kvcx) kbx=kbx-2*kvcx
!  if (kby < -kvcy) kby=kby+2*kvcy
!  if (kby > kvcy) kby=kby-2*kvcy
!  if (kbz < -kvcz) kbz=kbz+2*kvcz
!  if (kbz > kvcz) kbz=kbz-2*kvcz
!  if (myid == 0) write(400,*) kbx*mesh,kby*mesh,kbz*mesh, &
!  e(moi,kix*mesh,kiy*mesh,kiz*mesh)+e(moj,kjx*mesh,kjy*mesh,kjz*mesh)-e(moa,kax*mesh,kay*mesh,kaz*mesh)
!  enddo
!  enddo
!  enddo
!  enddo
!  enddo
!  enddo
!  enddo
!  enddo
!  enddo
!  enddo
!  enddo
!  enddo
! -- debug end

   DEALLOCATE(DOSHF,DOSMP,DOSGF,DOSGFW)
   IF (IALG >= 2) DEALLOCATE(E2)
   IF (IALG == 3) DEALLOCATE(E3,P)
   DEALLOCATE(LSWAPPED,SIMILARITY,WS,WF,WC,WC_PREVIOUS,E,WK)

   RETURN
END SUBROUTINE



SUBROUTINE DUMP_ENERGYBANDS_FCC
! DUMP THE HF/DFT ENERGY BANDS ON FINE MESH FOR FCC LATTICE

   USE MPI_F08
   USE CONSTANTS
   USE CONTROL
   USE STRUCTURE
   USE BASISSET
   USE INTEGRAL
   USE GRADIENT
   USE MULTIPOLE
   USE OEP

   IMPLICIT NONE
!  INCLUDE "mpif.h"
   INTEGER :: QX,QY,QZ,IX,IY,IZ,I,J,K,L,N
   INTEGER :: REDUCED_NCGS
   DOUBLE PRECISION :: ASYMMETRY,W,A
   DOUBLE PRECISION,ALLOCATABLE :: E(:)
   DOUBLE COMPLEX :: PX,PY,PZ,PQX,PQY,PQZ,AVE
   DOUBLE COMPLEX,ALLOCATABLE :: WS(:,:),WF(:,:),WC(:,:)
   DOUBLE COMPLEX,ALLOCATABLE :: WK(:,:)
   INTEGER,ALLOCATABLE :: PATHX(:),PATHY(:),PATHZ(:)
   DOUBLE PRECISION,ALLOCATABLE :: PATH(:)
   DOUBLE PRECISION :: PATH_TMP
   INTEGER :: PATHLENGTH

   IF ((KVCX /= KVCY).OR.(KVCX /= KVCZ)) THEN
    CALL WARNING('SET KVCX = KVCY = KVCZ')
    RETURN
   ENDIF

   ALLOCATE(WS(NCGS,NCGS),WF(NCGS,NCGS),WC(NCGS,NCGS),WK(NCGS,NCGS))
   ALLOCATE(E(NCGS))
   ASYMMETRY=0.0D0

   PATHLENGTH=0
   ! L(1/2,1/2,1/2) to Gamma(0,0,0) (Lambda)
   IF (MYID == 0) WRITE(6,'(A,I5)') 'L    ',PATHLENGTH+1
   DO IX=KVCX*16/2,0,-1
   IY=IX
   IZ=IX
    PATHLENGTH=PATHLENGTH+1
   ENDDO
   ! Gamma(0,0,0) to X(0,1/2,1/2) (Delta)
   IF (MYID == 0) WRITE(6,'(A,I5)') 'Gamma',PATHLENGTH+1
   DO IX=0,0
   DO IY=0,KVCY*16/2
   IZ=IY
    PATHLENGTH=PATHLENGTH+1
   ENDDO
   ENDDO
   ! X(0,1/2,1/2) to U(1/4,5/8,5/8)
   IF (MYID == 0) WRITE(6,'(A,I5)') 'X    ',PATHLENGTH+1
   DO IX=0,KVCX*16/4
   IY=KVCY*16/2+IX/2
   IZ=KVCZ*16/2+IX/2
    PATHLENGTH=PATHLENGTH+1
   ENDDO
   ! U(1/4,5/8,5/8) to Gamma(0,0,0) (Sigma)
   IF (MYID == 0) WRITE(6,'(A,I5)') 'U    ',PATHLENGTH+1
   DO IX=KVCX*16/4,0,-1
   IY=IX*5/2
   IZ=IX*5/2
    PATHLENGTH=PATHLENGTH+1
   ENDDO
   IF (MYID == 0) WRITE(6,'(A,I5)') 'Gamma',PATHLENGTH
   ALLOCATE(PATHX(PATHLENGTH),PATHY(PATHLENGTH),PATHZ(PATHLENGTH),PATH(PATHLENGTH))
   PATHLENGTH=0
   PATH_TMP=0.0D0
   ! L(1/2,1/2,1/2) to Gamma(0,0,0) (Lambda)
   DO IX=KVCX*16/2,0,-1
   IY=IX
   IZ=IX
    PATHLENGTH=PATHLENGTH+1
    PATHX(PATHLENGTH)=IX
    PATHY(PATHLENGTH)=IY
    PATHZ(PATHLENGTH)=IZ
    PATH_TMP=PATH_TMP+DSQRT(0.5D0**2+0.5D0**2+0.5D0**2)
    PATH(PATHLENGTH)=PATH_TMP
   ENDDO
   ! Gamma(0,0,0) to X(0,1/2,1/2) (Delta)
   DO IX=0,0
   DO IY=0,KVCY*16/2
   IZ=IY
    PATHLENGTH=PATHLENGTH+1
    PATHX(PATHLENGTH)=IX
    PATHY(PATHLENGTH)=IY
    PATHZ(PATHLENGTH)=IZ
    PATH_TMP=PATH_TMP+DSQRT(0.5D0**2+0.5D0**2)
    PATH(PATHLENGTH)=PATH_TMP
   ENDDO
   ENDDO
   ! X(0,1/2,1/2) to U(1/4,5/8,5/8)
   DO IX=0,KVCX*16/4
   IY=KVCY*16/2+IX/2
   IZ=KVCZ*16/2+IX/2
    PATHLENGTH=PATHLENGTH+1
    PATHX(PATHLENGTH)=IX
    PATHY(PATHLENGTH)=IY
    PATHZ(PATHLENGTH)=IZ
    PATH_TMP=PATH_TMP+DSQRT(0.25D0**2+0.125D0**2+0.125D0**2)
    PATH(PATHLENGTH)=PATH_TMP
   ENDDO
   ! U(1/4,5/8,5/8) to Gamma(0,0,0) (Sigma)
   DO IX=KVCX*16/4,0,-1
   IY=KVCY*16/2+IX/2
   IZ=KVCZ*16/2+IX/2
    PATHLENGTH=PATHLENGTH+1
    PATHX(PATHLENGTH)=IX
    PATHY(PATHLENGTH)=IY
    PATHZ(PATHLENGTH)=IZ
    PATH_TMP=PATH_TMP+DSQRT(0.25D0**2+0.625D0**2+0.625D0**2)
    PATH(PATHLENGTH)=PATH_TMP
   ENDDO

   IF (MYID == 0) WRITE(6,'(A)') 'DUMPING HF/DFT FCC ENERGY BANDS (L-G-X-U-G) IN fort.200+'

!  DUMP FOCK MATRIX
   IF ((IOPTN(9) >= 2).AND.(MYID == 0)) THEN
    WRITE(6,'(A)') 'FOCK MATRIX FOR CONTRACTED GAUSSIANS'
    CALL DUMP1(F_C,NCGS,CEL1X,CEL1Y,CEL1Z)
   ENDIF

!  ZERO SCRATCH DENSITY MATRIX AND ENERGY BAND
   P_C_OUT=0.0D0
   W_C=0.0D0
   E=0.0D0

!  LOOP OVER K VECTORS
   DO I=1,PATHLENGTH
 
    IX=PATHX(I)
    IY=PATHY(I)
    IZ=PATHZ(I)

!   CONSTRUCT DYNAMICAL OVERLAP & FOCK MATRICES
    IF (KVCX == 0) THEN
     PX=DCMPLX(0.0D0,0.0D0)
    ELSE
     PX=DCMPLX(0.0D0,PI*DFLOAT(IX)/DFLOAT(KVCX*16))
    ENDIF
    IF (KVCY == 0) THEN
     PY=DCMPLX(0.0D0,0.0D0)
    ELSE
     PY=DCMPLX(0.0D0,PI*DFLOAT(IY)/DFLOAT(KVCY*16))
    ENDIF
    IF (KVCZ == 0) THEN
     PZ=DCMPLX(0.0D0,0.0D0)
    ELSE
     PZ=DCMPLX(0.0D0,PI*DFLOAT(IZ)/DFLOAT(KVCZ*16))
    ENDIF
    WS=DCMPLX(0.0D0,0.0D0)
    WF=DCMPLX(0.0D0,0.0D0)
    DO J=1,NCGS
     DO K=1,NCGS
      DO QX=-CEL1X,CEL1X
      DO QY=-CEL1Y,CEL1Y
      DO QZ=-CEL1Z,CEL1Z
       PQX=PX*DFLOAT(QX)
       PQY=PY*DFLOAT(QY)
       PQZ=PZ*DFLOAT(QZ)
       WS(K,J)=WS(K,J)+S_C(K,J,QX,QY,QZ)*CDEXP(PQX+PQY+PQZ)
       WF(K,J)=WF(K,J)+F_C(K,J,QX,QY,QZ)*CDEXP(PQX+PQY+PQZ)
      ENDDO
      ENDDO
      ENDDO
     ENDDO
    ENDDO
    DO J=1,NCGS
     DO K=1,NCGS
      IF (CDABS(WF(K,J)-DCONJG(WF(J,K))) > ASYMMETRY) THEN
       ASYMMETRY=CDABS(WF(K,J)-DCONJG(WF(J,K)))
      ENDIF
      AVE=(WF(K,J)+DCONJG(WF(J,K)))/2.0D0
      WF(K,J)=AVE
      WF(J,K)=DCONJG(AVE)
     ENDDO
    ENDDO

!   CONVERT DYNAMICAL FOCK & OVERLAP FROM CARTESIAN TO SPHERICAL GAUSSIANS
    DO J=1,NCGS
     DO K=1,NCGS
      WK(K,J)=DCMPLX(0.0D0,0.0D0)
      DO L=1,NCGS
       WK(K,J)=WK(K,J)+C2S(K,L)*WS(L,J)
      ENDDO
     ENDDO
    ENDDO
    DO J=1,NCGS
     DO K=1,NCGS
      WS(K,J)=DCMPLX(0.0D0,0.0D0)
      DO L=1,NCGS
       WS(K,J)=WS(K,J)+C2S(J,L)*WK(K,L)
      ENDDO
     ENDDO
    ENDDO
    DO J=1,NCGS
     DO K=1,NCGS
      WK(K,J)=DCMPLX(0.0D0,0.0D0)
      DO L=1,NCGS
       WK(K,J)=WK(K,J)+C2S(K,L)*WF(L,J)
      ENDDO
     ENDDO
    ENDDO
    DO J=1,NCGS
     DO K=1,NCGS
      WF(K,J)=DCMPLX(0.0D0,0.0D0)
      DO L=1,NCGS
       WF(K,J)=WF(K,J)+C2S(J,L)*WK(K,L)
      ENDDO
     ENDDO
    ENDDO

!   CANONICAL ORTHOGONALIZATION OF DYNAMICAL OVERLAP MATRIX

    CALL HHBS(NCGS,NCGS,WS,E,WC)

    IF (IOPTN(9) == 3) WRITE(6,'(A,100F10.5:)') 'METRIC  = ',(E,J=1,NCGS)
    REDUCED_NCGS=0
    DO J=1,NCGS
     IF (E(J) > DOPTN(85)) THEN
      REDUCED_NCGS=REDUCED_NCGS+1
      DO K=1,NCGS
       WS(K,REDUCED_NCGS)=WC(K,J)/DSQRT(E(J))
      ENDDO
     ENDIF
    ENDDO

!   DIAGONALIZATION OF DYNAMICAL FOCK MATRIX
    DO J=1,REDUCED_NCGS
     DO K=1,NCGS
      WC(K,J)=DCMPLX(0.0D0,0.0D0)
      DO L=1,NCGS
       WC(K,J)=WC(K,J)+WF(K,L)*WS(L,J)
      ENDDO
     ENDDO
    ENDDO
    DO J=1,REDUCED_NCGS
     DO K=1,REDUCED_NCGS
      WF(K,J)=DCMPLX(0.0D0,0.0D0)
      DO L=1,NCGS
       WF(K,J)=WF(K,J)+DCONJG(WS(L,K))*WC(L,J)
      ENDDO
     ENDDO
    ENDDO

    CALL HHBS(NCGS,REDUCED_NCGS,WF,E,WC)

    DO J=1,REDUCED_NCGS
     IF (MYID == 0) WRITE(200+J,'(I5,2F25.15)') I,PATH(I),E(J)
    ENDDO

   ENDDO

   DEALLOCATE(WS,WF,WC,E,PATH,PATHX,PATHY,PATHZ,WK)

   RETURN
END SUBROUTINE
