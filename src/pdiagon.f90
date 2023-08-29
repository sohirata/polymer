SUBROUTINE TRED2(A,N,NP,D,E)
   INTEGER N,NP
   DOUBLE PRECISION A(NP,NP),D(NP),E(NP)
   INTEGER I,J,K,L
   DOUBLE PRECISION F,G,H,HH,SCALE
   DO 18 I=N,2,-1
     L=I-1
     H=0.0D0
     SCALE=0.0D0
     IF(L.GT.1)THEN
       DO 11 K=1,L
         SCALE=SCALE+DABS(A(I,K))
11     CONTINUE
       IF(SCALE.EQ.0.0D0)THEN
         E(I)=A(I,L)
       ELSE
         DO 12 K=1,L
           A(I,K)=A(I,K)/SCALE
           H=H+A(I,K)**2
12       CONTINUE
         F=A(I,L)
         G=-DSIGN(DSQRT(H),F)
         E(I)=SCALE*G
         H=H-F*G
         A(I,L)=F-G
         F=0.0D0
         DO 15 J=1,L
!  OMIT FOLLOWING LINE IF FINDING ONLY EIGENVALUES
           A(J,I)=A(I,J)/H
           G=0.0D0
           DO 13 K=1,J
             G=G+A(J,K)*A(I,K)
13         CONTINUE
           DO 14 K=J+1,L
             G=G+A(K,J)*A(I,K)
14         CONTINUE
           E(J)=G/H
           F=F+E(J)*A(I,J)
15       CONTINUE
         HH=F/(H+H)
         DO 17 J=1,L
           F=A(I,J)
           G=E(J)-HH*F
           E(J)=G
           DO 16 K=1,J
             A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16         CONTINUE
17       CONTINUE
       ENDIF
     ELSE
       E(I)=A(I,L)
     ENDIF
     D(I)=H
18 CONTINUE
!  OMIT FOLLOWING LINE IF FINDING ONLY EIGENVALUES.
   D(1)=0.0D0
   E(1)=0.0D0
   DO 24 I=1,N
!  DELETE LINES FROM HERE ...
     L=I-1
     IF(D(I).NE.0.0D0)THEN
       DO 22 J=1,L
         G=0.0D0
         DO 19 K=1,L
           G=G+A(I,K)*A(K,J)
19       CONTINUE
         DO 21 K=1,L
           A(K,J)=A(K,J)-G*A(K,I)
21       CONTINUE
22     CONTINUE
     ENDIF
!  ... TO HERE WHEN FINDING ONLY EIGENVALUES.
     D(I)=A(I,I)
!  ALSO DELETE LINES FROM HERE ...
     A(I,I)=1.0D0
     DO 23 J=1,L
       A(I,J)=0.0D0
       A(J,I)=0.0D0
23   CONTINUE
!  ... TO HERE WHEN FINDING ONLY EIGENVALUES.
24 CONTINUE
   RETURN
!  (C) COPR. 1986-92 NUMERICAL RECIPES SOFTWARE T.)-5I.
END SUBROUTINE



SUBROUTINE TQLI(D,E,N,NP,Z)
   INTEGER N,NP
   DOUBLE PRECISION D(NP),E(NP),Z(NP,NP)
!U USES PYTHAG
   INTEGER I,ITER,K,L,M
   DOUBLE PRECISION B,C,DD,F,G,P,R,S,PYTHAG
   DO 11 I=2,N
     E(I-1)=E(I)
11 CONTINUE
   E(N)=0.0D0
   DO 15 L=1,N
     ITER=0
1    DO 12 M=L,N-1
       DD=DABS(D(M))+DABS(D(M+1))
       IF (DABS(E(M))+DD.EQ.DD) GOTO 2
12   CONTINUE
     M=N
2    IF(M.NE.L)THEN
       IF(ITER.EQ.30) CALL PABORT('TOO MANY ITERATIONS IN TQLI')
       ITER=ITER+1
       G=(D(L+1)-D(L))/(2.0D0*E(L))
       R=PYTHAG(G,1.0D0)
       G=D(M)-D(L)+E(L)/(G+DSIGN(R,G))
       S=1.0D0
       C=1.0D0
       P=0.0D0
       DO 14 I=M-1,L,-1
         F=S*E(I)
         B=C*E(I)
         R=PYTHAG(F,G)
         E(I+1)=R
         IF(R.EQ.0.0D0)THEN
           D(I+1)=D(I+1)-P
           E(M)=0.0D0
           GOTO 1
         ENDIF
         S=F/R
         C=G/R
         G=D(I+1)-P
         R=(D(I)-G)*S+2.0D0*C*B
         P=S*R
         D(I+1)=G+P
         G=C*R-B
!  OMIT LINES FROM HERE ...
         DO 13 K=1,N
           F=Z(K,I+1)
           Z(K,I+1)=S*Z(K,I)+C*F
           Z(K,I)=C*Z(K,I)-S*F
13       CONTINUE
!  ... TO HERE WHEN FINDING ONLY EIGENVALUES.
14     CONTINUE
       D(L)=D(L)-P
       E(L)=G
       E(M)=0.0D0
       GOTO 1
     ENDIF
15 CONTINUE
   RETURN
!  (C) COPR. 1986-92 NUMERICAL RECIPES SOFTWARE T.)-5I.
END SUBROUTINE



SUBROUTINE TRED2_EVALONLY(A,N,NP,D,E)
   INTEGER N,NP
   DOUBLE PRECISION A(NP,NP),D(NP),E(NP)
   INTEGER I,J,K,L
   DOUBLE PRECISION F,G,H,HH,SCALE
   DO 18 I=N,2,-1
     L=I-1
     H=0.0D0
     SCALE=0.0D0
     IF(L.GT.1)THEN
       DO 11 K=1,L
         SCALE=SCALE+DABS(A(I,K))
11     CONTINUE
       IF(SCALE.EQ.0.0D0)THEN
         E(I)=A(I,L)
       ELSE
         DO 12 K=1,L
           A(I,K)=A(I,K)/SCALE
           H=H+A(I,K)**2
12       CONTINUE
         F=A(I,L)
         G=-DSIGN(DSQRT(H),F)
         E(I)=SCALE*G
         H=H-F*G
         A(I,L)=F-G
         F=0.0D0
         DO 15 J=1,L
           G=0.0D0
           DO 13 K=1,J
             G=G+A(J,K)*A(I,K)
13         CONTINUE
           DO 14 K=J+1,L
             G=G+A(K,J)*A(I,K)
14         CONTINUE
           E(J)=G/H
           F=F+E(J)*A(I,J)
15       CONTINUE
         HH=F/(H+H)
         DO 17 J=1,L
           F=A(I,J)
           G=E(J)-HH*F
           E(J)=G
           DO 16 K=1,J
             A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16         CONTINUE
17       CONTINUE
       ENDIF
     ELSE
       E(I)=A(I,L)
     ENDIF
     D(I)=H
18 CONTINUE
   E(1)=0.0D0
   DO 24 I=1,N
     D(I)=A(I,I)
24 CONTINUE
   RETURN
!  (C) COPR. 1986-92 NUMERICAL RECIPES SOFTWARE T.)-5I.
END SUBROUTINE



SUBROUTINE TQLI_EVALONLY(D,E,N,NP)
   INTEGER N,NP
   DOUBLE PRECISION D(NP),E(NP)
!U USES PYTHAG
   INTEGER I,ITER,K,L,M
   DOUBLE PRECISION B,C,DD,F,G,P,R,S,PYTHAG
   DO 11 I=2,N
     E(I-1)=E(I)
11 CONTINUE
   E(N)=0.0D0
   DO 15 L=1,N
     ITER=0
1    DO 12 M=L,N-1
       DD=DABS(D(M))+DABS(D(M+1))
       IF (DABS(E(M))+DD.EQ.DD) GOTO 2
12   CONTINUE
     M=N
2    IF(M.NE.L)THEN
       IF(ITER.EQ.30) CALL PABORT('TOO MANY ITERATIONS IN TQLI')
       ITER=ITER+1
       G=(D(L+1)-D(L))/(2.0D0*E(L))
       R=PYTHAG(G,1.0D0)
       G=D(M)-D(L)+E(L)/(G+DSIGN(R,G))
       S=1.0D0
       C=1.0D0
       P=0.0D0
       DO 14 I=M-1,L,-1
         F=S*E(I)
         B=C*E(I)
         R=PYTHAG(F,G)
         E(I+1)=R
         IF(R.EQ.0.0D0)THEN
           D(I+1)=D(I+1)-P
           E(M)=0.0D0
           GOTO 1
         ENDIF
         S=F/R
         C=G/R
         G=D(I+1)-P
         R=(D(I)-G)*S+2.0D0*C*B
         P=S*R
         D(I+1)=G+P
         G=C*R-B
14     CONTINUE
       D(L)=D(L)-P
       E(L)=G
       E(M)=0.0D0
       GOTO 1
     ENDIF
15 CONTINUE
   RETURN
!  (C) COPR. 1986-92 NUMERICAL RECIPES SOFTWARE T.)-5I.
END SUBROUTINE



FUNCTION PYTHAG(A,B)
   DOUBLE PRECISION A,B,PYTHAG
   DOUBLE PRECISION ABSA,ABSB
   ABSA=DABS(A)
   ABSB=DABS(B)
   IF(ABSA.GT.ABSB)THEN
     PYTHAG=ABSA*DSQRT(1.0D0+(ABSB/ABSA)**2)
   ELSE
     IF(ABSB.EQ.0.0D0)THEN
       PYTHAG=0.0D0
     ELSE
       PYTHAG=ABSB*DSQRT(1.0D0+(ABSA/ABSB)**2)
     ENDIF
   ENDIF
   RETURN
!  (C) COPR. 1986-92 NUMERICAL RECIPES SOFTWARE T.)-5I.
END FUNCTION



SUBROUTINE SVBKSB(A,N,NP,B)
   INTEGER N,NP
   DOUBLE PRECISION A(NP,NP),B(NP),C(N)
   INTEGER I,J

   DO I=1,N
    C(I)=0.0D0
    DO J=1,N
     C(I)=C(I)+A(I,J)*B(J)
    ENDDO
   ENDDO
   DO I=1,N
    B(I)=C(I)
   ENDDO
   RETURN
!  (C) COPR. 1986-92 NUMERICAL RECIPES SOFTWARE T.)-5I.
END SUBROUTINE



SUBROUTINE PIKSRT_EVALONLY(N,ARR)
   INTEGER N
   DOUBLE PRECISION ARR(N)
   INTEGER I,J
   DOUBLE PRECISION A

   DO 12 J=2,N
    A=ARR(J)
    DO 11 I=J-1,1,-1
     IF(ARR(I).LE.A)GOTO 10
      ARR(I+1)=ARR(I)
11   CONTINUE
     I=0
10   ARR(I+1)=A
12 CONTINUE
   RETURN
!  (C) COPR. 1986-92 NUMERICAL RECIPES SOFTWARE T.)-5I.
END SUBROUTINE



SUBROUTINE PIKSRT(N,NP,ARR,V,W)
   INTEGER N,NP
   DOUBLE PRECISION ARR(NP),V(NP,NP),W(NP)
   INTEGER I,J,K
   DOUBLE PRECISION A

   DO 12 J=2,N
    A=ARR(J)
    DO 14 K=1,N
     W(K)=V(K,J)
14  CONTINUE
    DO 11 I=J-1,1,-1
     IF(ARR(I).LE.A)GOTO 10
      ARR(I+1)=ARR(I)
      DO 15 K=1,N
       V(K,I+1)=V(K,I)
15    CONTINUE
11   CONTINUE
     I=0
10   ARR(I+1)=A
     DO 13 K=1,N
      V(K,I+1)=W(K)
13   CONTINUE
12 CONTINUE
   RETURN
END SUBROUTINE



SUBROUTINE PIKSRT2(N,NP,E1,E2,V,W)
   INTEGER N,NP
   DOUBLE PRECISION E1(NP),E2(NP),V(NP,NP),W(NP)
   INTEGER I,J,K
   DOUBLE PRECISION A,B

   DO 12 J=2,N
    A=E1(J)
    B=E2(J)
    DO 14 K=1,N
     W(K)=V(K,J)
14  CONTINUE
    DO 11 I=J-1,1,-1
     IF(E1(I).LE.A)GOTO 10
      E1(I+1)=E1(I)
      E2(I+1)=E2(I)
      DO 15 K=1,N
       V(K,I+1)=V(K,I)
15    CONTINUE
11   CONTINUE
     I=0
10   E1(I+1)=A
     E2(I+1)=B
     DO 13 K=1,N
      V(K,I+1)=W(K)
13   CONTINUE
12 CONTINUE
   RETURN
END SUBROUTINE




SUBROUTINE PIKSRT3(N,NP,E1,E2,E3)
   INTEGER N,NP
   DOUBLE PRECISION E1(NP),E2(NP),E3(NP)
   INTEGER I,J
   DOUBLE PRECISION A,B,C

   DO 12 J=2,N
    A=E1(J)
    B=E2(J)
    C=E3(J)
    DO 11 I=J-1,1,-1
     IF(E1(I).LE.A)GOTO 10
      E1(I+1)=E1(I)
      E2(I+1)=E2(I)
      E3(I+1)=E3(I)
11   CONTINUE
     I=0
10   E1(I+1)=A
     E2(I+1)=B
     E3(I+1)=C
12 CONTINUE
   RETURN
END SUBROUTINE



SUBROUTINE PIKSRT4(N,NP,E1,E2,E3,E4)
   INTEGER N,NP
   DOUBLE PRECISION E1(NP),E2(NP),E3(NP),E4(NP)
   INTEGER I,J
   DOUBLE PRECISION A,B,C,D

   DO 12 J=2,N
    A=E1(J)
    B=E2(J)
    C=E3(J)
    D=E4(J)
    DO 11 I=J-1,1,-1
     IF(E1(I).LE.A)GOTO 10
      E1(I+1)=E1(I)
      E2(I+1)=E2(I)
      E3(I+1)=E3(I)
      E4(I+1)=E4(I)
11   CONTINUE
     I=0
10   E1(I+1)=A
     E2(I+1)=B
     E3(I+1)=C
     E4(I+1)=D
12 CONTINUE
   RETURN
END SUBROUTINE




SUBROUTINE PIKSRT5(N,NP,E1,E2,E3,E4,E5)
   INTEGER N,NP
   DOUBLE PRECISION E1(NP),E2(NP),E3(NP),E4(NP),E5(NP)
   INTEGER I,J
   DOUBLE PRECISION A,B,C,D,E

   DO 12 J=2,N
    A=E1(J)
    B=E2(J)
    C=E3(J)
    D=E4(J)
    E=E5(J)
    DO 11 I=J-1,1,-1
     IF(E1(I).LE.A)GOTO 10
      E1(I+1)=E1(I)
      E2(I+1)=E2(I)
      E3(I+1)=E3(I)
      E4(I+1)=E4(I)
      E5(I+1)=E5(I)
11   CONTINUE
     I=0
10   E1(I+1)=A
     E2(I+1)=B
     E3(I+1)=C
     E4(I+1)=D
     E5(I+1)=E
12 CONTINUE
   RETURN
END SUBROUTINE

