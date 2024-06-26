SUBROUTINE LEGEND

   USE MPI_F08
   USE CONSTANTS

   IMPLICIT NONE
!  INCLUDE "mpif.h"

   IF (MYID == 0) THEN
     WRITE(6,'(A)') ' '
     WRITE(6,'(A)') '**************************'
     WRITE(6,'(A)') '* POLYMER 3D Version 1.0 *'
     WRITE(6,'(A)') '**************************'
     WRITE(6,'(A)') ' '
     WRITE(6,'(A)') '(c) So Hirata, University of Illinois at Urbana-Champaign'
     WRITE(6,'(A)') ' '
   ENDIF

   RETURN
END SUBROUTINE
