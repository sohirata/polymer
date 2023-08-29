SUBROUTINE SHELL1
! CLASSIFY ORBITAL BASIS FUNCTIONS INTO SHELLS.
! A SHELL CONSISTS OF BASIS FUNCTIONS WHICH SHARE THE ONE AND THE SAME CENTER AND THE EXPONENT.

   USE CONTROL
   USE BASISSET
   
   IMPLICIT NONE
   INTEGER :: I,J,K,L,FLAG
   
   IF (IOPTN(7) >= 1) THEN
    L=1
   ELSE
    L=0
   ENDIF

   J=0
   DO I=1,NPGS
    K=PAX(I)+PAY(I)+PAZ(I)
    IF (K == 0) THEN
     J=J+1
     P_SHELL(J,0,0,0)=I
     P_SHELL_ANG(J)=0
    ENDIF
   ENDDO
   IF (J /= NPSHELL) CALL PABORT('AN INTERNAL PROGRAM ERROR IS DETECTED')
   
   DO I=1,NPGS
    K=PAX(I)+PAY(I)+PAZ(I)
    FLAG=0
    IF (K+L > BASIS_LMAX) THEN
     CALL PABORT('BASIS ANGULAR MOMENTUM TOO HIGH')
    ELSE IF (K > 0) THEN
     DO J=1,NPSHELL
      IF ((PGSX(P_SHELL(J,0,0,0)) == PGSX(I)).AND.(PGSY(P_SHELL(J,0,0,0)) == PGSY(I)).AND. &
      (PGSZ(P_SHELL(J,0,0,0)) == PGSZ(I)).AND.(ZT(P_SHELL(J,0,0,0)) == ZT(I))) THEN
       FLAG=1
       P_SHELL(J,PAX(I),PAY(I),PAZ(I))=I
       IF (K > P_SHELL_ANG(J)) P_SHELL_ANG(J)=K
      ENDIF
     ENDDO
     IF (FLAG == 0) CALL PABORT('P AND HIGHER FUNCTIONS MUST HAVE LOWER ANGULAR MOMENTUM COUNTERPARTS IN THE BASIS SET')
    ENDIF
   ENDDO
   
   IF (IOPTN(9) >= 2) THEN
    WRITE(6,'(A)') 'BASIS SET SHELL STRUCTURE'
    DO I=1,NPSHELL
     SELECT CASE(P_SHELL_ANG(I))
      CASE(0)
       WRITE(6,'(I3,A3,I3)') I,'S',P_SHELL(I,0,0,0)
      CASE(1)
       WRITE(6,'(I3,A3,I3,A3,3I3)') I,'S',P_SHELL(I,0,0,0), &
       'P',P_SHELL(I,1,0,0),P_SHELL(I,0,1,0),P_SHELL(I,0,0,1)
      CASE(2)
       WRITE(6,'(I3,A3,I3,A3,3I3,A3,6I3)') I,'S',P_SHELL(I,0,0,0), &
       'P',P_SHELL(I,1,0,0),P_SHELL(I,0,1,0),P_SHELL(I,0,0,1), &
       'D',P_SHELL(I,2,0,0),P_SHELL(I,0,2,0),P_SHELL(I,0,0,2), &
       P_SHELL(I,1,1,0),P_SHELL(I,0,1,1),P_SHELL(I,1,0,1)
      CASE(3)
       WRITE(6,'(I3,A3,I3,A3,3I3,A3,6I3,A3,10I3)') I,'S',P_SHELL(I,0,0,0), &
       'P',P_SHELL(I,1,0,0),P_SHELL(I,0,1,0),P_SHELL(I,0,0,1), &
       'D',P_SHELL(I,2,0,0),P_SHELL(I,0,2,0),P_SHELL(I,0,0,2), &
       P_SHELL(I,1,1,0),P_SHELL(I,0,1,1),P_SHELL(I,1,0,1), &
       'F',P_SHELL(I,3,0,0),P_SHELL(I,0,3,0),P_SHELL(I,0,0,3), &
       P_SHELL(I,2,1,0),P_SHELL(I,0,2,1),P_SHELL(I,1,0,2), &
       P_SHELL(I,1,2,0),P_SHELL(I,0,1,2),P_SHELL(I,2,0,1), &
       P_SHELL(I,1,1,1)
     END SELECT
    ENDDO
   ENDIF
   
   RETURN
END SUBROUTINE



SUBROUTINE SHELL2
! CLASSIFY AUXILIARY BASIS FUNCTIONS INTO SHELLS.
! A SHELL CONSISTS OF BASIS FUNCTIONS WHICH SHARE THE ONE AND THE SAME CENTER AND THE EXPONENT.

   USE CONTROL
   USE AUXILIARY
   
   IMPLICIT NONE
   INTEGER :: I,J,K,L,FLAG
   
   IF (IOPTN(7) >= 1) THEN
    L=1
   ELSE
    L=0
   ENDIF
   
   J=0
   DO I=1,NAGS
    K=AAX(I)+AAY(I)+AAZ(I)
    IF (K == 0) THEN
     J=J+1
     A_SHELL(J,0,0,0)=I
     A_SHELL_ANG(J)=0
    ENDIF
   ENDDO
   IF (J /= NASHELL) CALL PABORT('AN INTERNAL PROGRAM ERROR IS DETECTED')
   
   DO I=1,NAGS
    K=AAX(I)+AAY(I)+AAZ(I)
    FLAG=0
    IF (K+L > 4) THEN
     CALL PABORT('A H- OR HIGHER ANGULAR MOMENTUM FUNCTIONS CANNOT BE USED AS ORBITAL BASIS FUNCTIONS')
    ELSE IF (K > 0) THEN
     DO J=1,NASHELL
      IF ((AGSX(A_SHELL(J,0,0,0)) == AGSX(I)).AND.(AGSY(A_SHELL(J,0,0,0)) == AGSY(I)).AND. &
      (AGSZ(A_SHELL(J,0,0,0)) == AGSZ(I)).AND.(AZT(A_SHELL(J,0,0,0)) == AZT(I))) THEN
       FLAG=1
       A_SHELL(J,AAX(I),AAY(I),AAZ(I))=I
       IF (K > A_SHELL_ANG(J)) A_SHELL_ANG(J)=K
      ENDIF
     ENDDO
     IF (FLAG == 0) CALL PABORT('P-, D-, F-, AND G-TYPE FUNCTIONS MUST HAVE LOWER ANGULAR MOMENTUM COUNTERPARTS IN THE BASIS SET')
    ENDIF
   ENDDO
   
   IF (IOPTN(9) >= 2) THEN
    WRITE(6,'(A)') 'BASIS SET SHELL STRUCTURE'
    DO I=1,NASHELL
     SELECT CASE(A_SHELL_ANG(I))
      CASE(0)
       WRITE(6,'(I3,A3,I3)') I,'S',A_SHELL(I,0,0,0)
      CASE(1)
       WRITE(6,'(I3,A3,I3,A3,3I3)') I,'S',A_SHELL(I,0,0,0), &
       'P',A_SHELL(I,1,0,0),A_SHELL(I,0,1,0),A_SHELL(I,0,0,1)
      CASE(2)
       WRITE(6,'(I3,A3,I3,A3,3I3,A3,6I3)') I,'S',A_SHELL(I,0,0,0), &
       'P',A_SHELL(I,1,0,0),A_SHELL(I,0,1,0),A_SHELL(I,0,0,1), &
       'D',A_SHELL(I,2,0,0),A_SHELL(I,0,2,0),A_SHELL(I,0,0,2), &
       A_SHELL(I,1,1,0),A_SHELL(I,0,1,1),A_SHELL(I,1,0,1)
      CASE(3)
       WRITE(6,'(I3,A3,I3,A3,3I3,A3,6I3,A3,10I3)') I,'S',A_SHELL(I,0,0,0), &
       'P',A_SHELL(I,1,0,0),A_SHELL(I,0,1,0),A_SHELL(I,0,0,1), &
       'D',A_SHELL(I,2,0,0),A_SHELL(I,0,2,0),A_SHELL(I,0,0,2), &
       A_SHELL(I,1,1,0),A_SHELL(I,0,1,1),A_SHELL(I,1,0,1), &
       'F',A_SHELL(I,3,0,0),A_SHELL(I,0,3,0),A_SHELL(I,0,0,3), &
       A_SHELL(I,2,1,0),A_SHELL(I,0,2,1),A_SHELL(I,1,0,2), &
       A_SHELL(I,1,2,0),A_SHELL(I,0,1,2),A_SHELL(I,2,0,1), &
       A_SHELL(I,1,1,1)
      CASE(4)
       WRITE(6,'(I3,A3,I3,A3,3I3,A3,6I3,A3,10I3,A7)') I,'S',A_SHELL(I,0,0,0), &
       'P',A_SHELL(I,1,0,0),A_SHELL(I,0,1,0),A_SHELL(I,0,0,1), &
       'D',A_SHELL(I,2,0,0),A_SHELL(I,0,2,0),A_SHELL(I,0,0,2), &
       A_SHELL(I,1,1,0),A_SHELL(I,0,1,1),A_SHELL(I,1,0,1), &
       'F',A_SHELL(I,3,0,0),A_SHELL(I,0,3,0),A_SHELL(I,0,0,3), &
       A_SHELL(I,2,1,0),A_SHELL(I,0,2,1),A_SHELL(I,1,0,2), &
       A_SHELL(I,1,2,0),A_SHELL(I,0,1,2),A_SHELL(I,2,0,1), &
       A_SHELL(I,1,1,1),'G ...'
     END SELECT
    ENDDO
   ENDIF
   
   RETURN
END SUBROUTINE
