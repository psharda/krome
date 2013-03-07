
!!$*DECK DGEFA

subroutine dgefa (a,lda,n,ipvt,info)
!!$ C***BEGIN PROLOGUE DGEFA
!!$ C***PURPOSE FACTOR A MATRIX USING GAUSSIAN ELIMINATION.
!!$ C***CATEGORY D2A1
!!$ C***TYPE   real*8::(SGEFA-S,DGEFA-D,CGEFA-C)
!!$ C***KEYWORDS GENERAL MATRIX,LINEAR ALGEBRA,LINPACK,
!!$ C       MATRIX FACTORIZATION
!!$ C***AUTHOR MOLER,C. B.,(U. OF NEW MEXICO)
!!$ C***DESCRIPTION
!!$ C
!!$ C   DGEFA FACTORS A real*8::MATRIX BY GAUSSIAN ELIMINATION.
!!$ C
!!$ C   DGEFA IS USUALLY CALLED BY DGECO,BUT IT CAN BE CALLED
!!$ C   DIRECTLY WITH A SAVING IN TIME IF RCOND IS NOT NEEDED.
!!$ C   (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .
!!$ C
!!$ C   ON ENTRY
!!$ C
!!$ C    A    DOUBLE PRECISION(LDA,N)
!!$ C        THE MATRIX TO BE FACTORED.
!!$ C
!!$ C    LDA   INTEGER
!!$ C        THE LEADING DIMENSION OF THE ARRAY A .
!!$ C
!!$ C    N    INTEGER
!!$ C        THE ORDER OF THE MATRIX A .
!!$ C
!!$ C   ON RETURN
!!$ C
!!$ C    A    AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
!!$ C        WHICH WERE USED TO OBTAIN IT.
!!$ C        THE FACTORIZATION CAN BE WRITTEN A = L*U WHERE
!!$ C        L IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!!$ C        TRIANGULAR MATRICES AND U IS UPPER TRIANGULAR.
!!$ C
!!$ C    IPVT  INTEGER(N)
!!$ C        AN integer::VECTOR OF PIVOT INDICES.
!!$ C
!!$ C    INFO  INTEGER
!!$ C        = 0 NORMAL VALUE.
!!$ C        = K IF U(K,K) .EQ. 0.0 . THIS IS NOT AN ERROR
!!$ C           CONDITION FOR THIS SUBROUTINE,BUT IT DOES
!!$ C           INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO
!!$ C           IF CALLED. USE RCOND IN DGECO FOR A RELIABLE
!!$ C           INDICATION OF SINGULARITY.
!!$ C
!!$ C***REFERENCES J. J. DONGARRA,J. R. BUNCH,C. B. MOLER,AND G. W.
!!$ C         STEWART,LINPACK USERS' GUIDE,SIAM,1979.
!!$ C***ROUTINES CALLED DAXPY,DSCAL,IDAMAX
!!$ C***REVISION HISTORY (YYMMDD)
!!$ C  780814 DATE WRITTEN
!!$ C  890831 MODIFIED ARRAY DECLARATIONS. (WRB)
!!$ C  890831 REVISION DATE FROM VERSION 3.2
!!$ C  891214 PROLOGUE CONVERTED TO VERSION 4.0 FORMAT. (BAB)
!!$ C  900326 REMOVED DUPLICATE INFORMATION FROM DESCRIPTION SECTION.
!!$ C      (WRB)
!!$ C  920501 REFORMATTED THE REFERENCES SECTION. (WRB)
!!$ C***END PROLOGUE DGEFA
  integer::lda,n,ipvt(*),info
  real*8::a(lda,*)
!!$ C
  real*8::t
  integer::idamax,j,k,kp1,l,nm1
!!$ C
!!$ C   GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!!$ C
!!$ C***FIRST EXECUTABLE STATEMENT DGEFA
  info = 0
  nm1 = n - 1
  if (nm1 .lt. 1) go to 70
  do k = 1,nm1
     kp1 = k + 1
!!$ C
!!$ C    FIND L = PIVOT INDEX
!!$ C
     l = idamax(n-k+1,a(k,k),1) + k - 1
     ipvt(k) = l
!!$ C
!!$ C    ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!!$ C
     if (a(l,k) .eq. 0.0d0) go to 40
!!$ C
!!$ C      INTERCHANGE IF NECESSARY
!!$ C
     if (l .eq. k) go to 10
     t = a(l,k)
     a(l,k) = a(k,k)
     a(k,k) = t
10   continue
!!$ C
!!$ C      COMPUTE MULTIPLIERS
!!$ C
     t = -1.0d0/a(k,k)
     call dscal(n-k,t,a(k+1,k),1)
!!$ C
!!$ C      ROW ELIMINATION WITH COLUMN INDEXING
!!$ C
     do j = kp1,n
        t = a(l,j)
        if (l .eq. k) go to 20
        a(l,j) = a(k,j)
        a(k,j) = t
20      continue
        call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
30      continue
     end do
     go to 50
40   continue
     info = k
50   continue
60   continue
  end do
70 continue
  ipvt(n) = n
  if (a(n,n) .eq. 0.0d0) info = n
  return
end subroutine dgefa

!!$*DECK DGESL

subroutine dgesl (a,lda,n,ipvt,b,job)
!!$ C***BEGIN PROLOGUE DGESL
!!$ C***PURPOSE SOLVE THE REAL SYSTEM A*X=B OR TRANS(A)*X=B USING THE
!!$ C      FACTORS COMPUTED BY DGECO OR DGEFA.
!!$ C***CATEGORY D2A1
!!$ C***TYPE   real*8::(SGESL-S,DGESL-D,CGESL-C)
!!$ C***KEYWORDS LINEAR ALGEBRA,LINPACK,MATRIX,SOLVE
!!$ C***AUTHOR MOLER,C. B.,(U. OF NEW MEXICO)
!!$ C***DESCRIPTION
!!$ C
!!$ C   DGESL SOLVES THE real*8::SYSTEM
!!$ C   A * X = B OR TRANS(A) * X = B
!!$ C   USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
!!$ C
!!$ C   ON ENTRY
!!$ C
!!$ C    A    DOUBLE PRECISION(LDA,N)
!!$ C        THE OUTPUT FROM DGECO OR DGEFA.
!!$ C
!!$ C    LDA   INTEGER
!!$ C        THE LEADING DIMENSION OF THE ARRAY A .
!!$ C
!!$ C    N    INTEGER
!!$ C        THE ORDER OF THE MATRIX A .
!!$ C
!!$ C    IPVT  INTEGER(N)
!!$ C        THE PIVOT VECTOR FROM DGECO OR DGEFA.
!!$ C
!!$ C    B    DOUBLE PRECISION(N)
!!$ C        THE RIGHT HAND SIDE VECTOR.
!!$ C
!!$ C    JOB   INTEGER
!!$ C        = 0     TO SOLVE A*X = B ,
!!$ C        = NONZERO  TO SOLVE TRANS(A)*X = B WHERE
!!$ C              TRANS(A) IS THE TRANSPOSE.
!!$ C
!!$ C   ON RETURN
!!$ C
!!$ C    B    THE SOLUTION VECTOR X .
!!$ C
!!$ C   ERROR CONDITION
!!$ C
!!$ C    A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
!!$ C    ZERO ON THE DIAGONAL. TECHNICALLY THIS INDICATES SINGULARITY
!!$ C    BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
!!$ C    SETTING OF LDA . IT WILL NOT OCCUR IF THE SUBROUTINES ARE
!!$ C    CALLED CORRECTLY AND IF DGECO HAS SET RCOND .GT. 0.0
!!$ C    OR DGEFA HAS SET INFO .EQ. 0 .
!!$ C
!!$ C   TO COMPUTE INVERSE(A) * C WHERE C IS A MATRIX
!!$ C   WITH P COLUMNS
!!$ C      CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
!!$ C      IF (RCOND IS TOO SMALL) GO TO ...
!!$ C      DO 10 J = 1,P
!!$ C       CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
!!$ C    10 CONTINUE
!!$ C
!!$ C***REFERENCES J. J. DONGARRA,J. R. BUNCH,C. B. MOLER,AND G. W.
!!$ C         STEWART,LINPACK USERS' GUIDE,SIAM,1979.
!!$ C***ROUTINES CALLED DAXPY,DDOT
!!$ C***REVISION HISTORY (YYMMDD)
!!$ C  780814 DATE WRITTEN
!!$ C  890831 MODIFIED ARRAY DECLARATIONS. (WRB)
!!$ C  890831 REVISION DATE FROM VERSION 3.2
!!$ C  891214 PROLOGUE CONVERTED TO VERSION 4.0 FORMAT. (BAB)
!!$ C  900326 REMOVED DUPLICATE INFORMATION FROM DESCRIPTION SECTION.
!!$ C      (WRB)
!!$ C  920501 REFORMATTED THE REFERENCES SECTION. (WRB)
!!$ C***END PROLOGUE DGESL
  integer::lda,n,ipvt(*),job
  real*8::a(lda,*),b(*)
!!$ C
  real*8::ddot,t
  integer::k,kb,l,nm1
!!$ C***FIRST EXECUTABLE STATEMENT DGESL
  nm1 = n - 1
  if (job .ne. 0) go to 50
!!$ C
!!$ C    JOB = 0 ,SOLVE A * X = B
!!$ C    FIRST SOLVE L*Y = B
!!$ C
  if (nm1 .lt. 1) go to 30
  do k = 1,nm1
     l = ipvt(k)
     t = b(l)
     if (l .eq. k) go to 10
     b(l) = b(k)
     b(k) = t
10   continue
     call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
20   continue
  end do
30 continue
!!$ C
!!$ C    NOW SOLVE U*X = Y
!!$ C
  do kb = 1,n
     k = n + 1 - kb
     b(k) = b(k)/a(k,k)
     t = -b(k)
     call daxpy(k-1,t,a(1,k),1,b(1),1)
40   continue
  end do
  go to 100
50 continue
!!$ C
!!$ C    JOB = NONZERO,SOLVE TRANS(A) * X = B
!!$ C    FIRST SOLVE TRANS(U)*Y = B
!!$ C
  do k = 1,n
     t = ddot(k-1,a(1,k),1,b(1),1)
     b(k) = (b(k) - t)/a(k,k)
60   continue
  end do
!!$ C
!!$ C    NOW SOLVE TRANS(L)*X = Y
!!$ C
  if (nm1 .lt. 1) go to 90
  do kb = 1,nm1
     k = n - kb
     b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
     l = ipvt(k)
     if (l .eq. k) go to 70
     t = b(l)
     b(l) = b(k)
     b(k) = t
70   continue
80   continue
  end do
90 continue
100 continue
  return
end subroutine dgesl

!!$*DECK DGBFA

subroutine dgbfa (abd,lda,n,ml,mu,ipvt,info)
!!$ C***BEGIN PROLOGUE DGBFA
!!$ C***PURPOSE FACTOR A BAND MATRIX USING GAUSSIAN ELIMINATION.
!!$ C***CATEGORY D2A2
!!$ C***TYPE   real*8::(SGBFA-S,DGBFA-D,CGBFA-C)
!!$ C***KEYWORDS BANDED,LINEAR ALGEBRA,LINPACK,MATRIX FACTORIZATION
!!$ C***AUTHOR MOLER,C. B.,(U. OF NEW MEXICO)
!!$ C***DESCRIPTION
!!$ C
!!$ C   DGBFA FACTORS A real*8::BAND MATRIX BY ELIMINATION.
!!$ C
!!$ C   DGBFA IS USUALLY CALLED BY DGBCO,BUT IT CAN BE CALLED
!!$ C   DIRECTLY WITH A SAVING IN TIME IF RCOND IS NOT NEEDED.
!!$ C
!!$ C   ON ENTRY
!!$ C
!!$ C    ABD   DOUBLE PRECISION(LDA,N)
!!$ C        CONTAINS THE MATRIX IN BAND STORAGE. THE COLUMNS
!!$ C        OF THE MATRIX ARE STORED IN THE COLUMNS OF ABD AND
!!$ C        THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
!!$ C        ML+1 THROUGH 2*ML+MU+1 OF ABD .
!!$ C        SEE THE COMMENTS BELOW FOR DETAILS.
!!$ C
!!$ C    LDA   INTEGER
!!$ C        THE LEADING DIMENSION OF THE ARRAY ABD .
!!$ C        LDA MUST BE .GE. 2*ML + MU + 1 .
!!$ C
!!$ C    N    INTEGER
!!$ C        THE ORDER OF THE ORIGINAL MATRIX.
!!$ C
!!$ C    ML   INTEGER
!!$ C        NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
!!$ C        0 .LE. ML .LT. N .
!!$ C
!!$ C    MU   INTEGER
!!$ C        NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!!$ C        0 .LE. MU .LT. N .
!!$ C        MORE EFFICIENT IF ML .LE. MU .
!!$ C   ON RETURN
!!$ C
!!$ C    ABD   AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
!!$ C        THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
!!$ C        THE FACTORIZATION CAN BE WRITTEN A = L*U WHERE
!!$ C        L IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!!$ C        TRIANGULAR MATRICES AND U IS UPPER TRIANGULAR.
!!$ C
!!$ C    IPVT  INTEGER(N)
!!$ C        AN integer::VECTOR OF PIVOT INDICES.
!!$ C
!!$ C    INFO  INTEGER
!!$ C        = 0 NORMAL VALUE.
!!$ C        = K IF U(K,K) .EQ. 0.0 . THIS IS NOT AN ERROR
!!$ C           CONDITION FOR THIS SUBROUTINE,BUT IT DOES
!!$ C           INDICATE THAT DGBSL WILL DIVIDE BY ZERO IF
!!$ C           CALLED. USE RCOND IN DGBCO FOR A RELIABLE
!!$ C           INDICATION OF SINGULARITY.
!!$ C
!!$ C   BAND STORAGE
!!$ C
!!$ C      IF A IS A BAND MATRIX,THE FOLLOWING PROGRAM SEGMENT
!!$ C      WILL SET UP THE INPUT.
!!$ C
!!$ C          ML = (BAND WIDTH BELOW THE DIAGONAL)
!!$ C          MU = (BAND WIDTH ABOVE THE DIAGONAL)
!!$ C          M = ML + MU + 1
!!$ C          DO 20 J = 1,N
!!$ C           I1 = MAX(1,J-MU)
!!$ C           I2 = MIN(N,J+ML)
!!$ C           DO 10 I = I1,I2
!!$ C             K = I - J + M
!!$ C             ABD(K,J) = A(I,J)
!!$ C        10  CONTINUE
!!$ C        20 CONTINUE
!!$ C
!!$ C      THIS USES ROWS ML+1 THROUGH 2*ML+MU+1 OF ABD .
!!$ C      IN ADDITION,THE FIRST ML ROWS IN ABD ARE USED FOR
!!$ C      ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
!!$ C      THE TOTAL NUMBER OF ROWS NEEDED IN ABD IS 2*ML+MU+1 .
!!$ C      THE ML+MU BY ML+MU UPPER LEFT TRIANGLE AND THE
!!$ C      ML BY ML LOWER RIGHT TRIANGLE ARE NOT REFERENCED.
!!$ C
!!$ C***REFERENCES J. J. DONGARRA,J. R. BUNCH,C. B. MOLER,AND G. W.
!!$ C         STEWART,LINPACK USERS' GUIDE,SIAM,1979.
!!$ C***ROUTINES CALLED DAXPY,DSCAL,IDAMAX
!!$ C***REVISION HISTORY (YYMMDD)
!!$ C  780814 DATE WRITTEN
!!$ C  890531 CHANGED ALL SPECIFIC INTRINSICS TO GENERIC. (WRB)
!!$ C  890831 MODIFIED ARRAY DECLARATIONS. (WRB)
!!$ C  890831 REVISION DATE FROM VERSION 3.2
!!$ C  891214 PROLOGUE CONVERTED TO VERSION 4.0 FORMAT. (BAB)
!!$ C  900326 REMOVED DUPLICATE INFORMATION FROM DESCRIPTION SECTION.
!!$ C      (WRB)
!!$ C  920501 REFORMATTED THE REFERENCES SECTION. (WRB)
!!$ C***END PROLOGUE DGBFA
  integer::lda,n,ml,mu,ipvt(*),info
  real*8::abd(lda,*)
!!$ C
  real*8::t
  integer::i,idamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
!!$ C
!!$ C***FIRST EXECUTABLE STATEMENT DGBFA
  m = ml + mu + 1
  info = 0
!!$ C
!!$ C   ZERO INITIAL FILL-IN COLUMNS
!!$ C
  j0 = mu + 2
  j1 = min(n,m) - 1
  if (j1 .lt. j0) go to 30
  do jz = j0,j1
     i0 = m + 1 - jz
     do i = i0,ml
        abd(i,jz) = 0.0d0
10      continue
     end do
20   continue
  end do
30 continue
  jz = j1
  ju = 0
!!$ C
!!$ C   GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!!$ C
  nm1 = n - 1
  if (nm1 .lt. 1) go to 130
  do k = 1,nm1
     kp1 = k + 1
!!$ C
!!$ C    ZERO NEXT FILL-IN COLUMN
!!$ C
     jz = jz + 1
     if (jz .gt. n) go to 50
     if (ml .lt. 1) go to 50
     do i = 1,ml
        abd(i,jz) = 0.0d0
40      continue
     end do
50   continue
!!$ C
!!$ C    FIND L = PIVOT INDEX
!!$ C
     lm = min(ml,n-k)
     l = idamax(lm+1,abd(m,k),1) + m - 1
     ipvt(k) = l + k - m
!!$ C
!!$ C    ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!!$ C
     if (abd(l,k) .eq. 0.0d0) go to 100
!!$ C
!!$ C      INTERCHANGE IF NECESSARY
!!$ C
     if (l .eq. m) go to 60
     t = abd(l,k)
     abd(l,k) = abd(m,k)
     abd(m,k) = t
60   continue
!!$ C
!!$ C      COMPUTE MULTIPLIERS
!!$ C
     t = -1.0d0/abd(m,k)
     call dscal(lm,t,abd(m+1,k),1)
!!$ C
!!$ C      ROW ELIMINATION WITH COLUMN INDEXING
!!$ C
     ju = min(max(ju,mu+ipvt(k)),n)
     mm = m
     if (ju .lt. kp1) go to 90
     do j = kp1,ju
        l = l - 1
        mm = mm - 1
        t = abd(l,j)
        if (l .eq. mm) go to 70
        abd(l,j) = abd(mm,j)
        abd(mm,j) = t
70      continue
        call daxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
80      continue
     end do
90   continue
     go to 110
100  continue
     info = k
110  continue
120  continue
  end do
130 continue
  ipvt(n) = n
  if (abd(m,n) .eq. 0.0d0) info = n
  return
end subroutine dgbfa

!!$*DECK DGBSL

subroutine dgbsl (abd,lda,n,ml,mu,ipvt,b,job)
!!$ C***BEGIN PROLOGUE DGBSL
!!$ C***PURPOSE SOLVE THE REAL BAND SYSTEM A*X=B OR TRANS(A)*X=B USING
!!$ C      THE FACTORS COMPUTED BY DGBCO OR DGBFA.
!!$ C***CATEGORY D2A2
!!$ C***TYPE   real*8::(SGBSL-S,DGBSL-D,CGBSL-C)
!!$ C***KEYWORDS BANDED,LINEAR ALGEBRA,LINPACK,MATRIX,SOLVE
!!$ C***AUTHOR MOLER,C. B.,(U. OF NEW MEXICO)
!!$ C***DESCRIPTION
!!$ C
!!$ C   DGBSL SOLVES THE real*8::BAND SYSTEM
!!$ C   A * X = B OR TRANS(A) * X = B
!!$ C   USING THE FACTORS COMPUTED BY DGBCO OR DGBFA.
!!$ C
!!$ C   ON ENTRY
!!$ C
!!$ C    ABD   DOUBLE PRECISION(LDA,N)
!!$ C        THE OUTPUT FROM DGBCO OR DGBFA.
!!$ C
!!$ C    LDA   INTEGER
!!$ C        THE LEADING DIMENSION OF THE ARRAY ABD .
!!$ C
!!$ C    N    INTEGER
!!$ C        THE ORDER OF THE ORIGINAL MATRIX.
!!$ C
!!$ C    ML   INTEGER
!!$ C        NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
!!$ C
!!$ C    MU   INTEGER
!!$ C        NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!!$ C
!!$ C    IPVT  INTEGER(N)
!!$ C        THE PIVOT VECTOR FROM DGBCO OR DGBFA.
!!$ C
!!$ C    B    DOUBLE PRECISION(N)
!!$ C        THE RIGHT HAND SIDE VECTOR.
!!$ C
!!$ C    JOB   INTEGER
!!$ C        = 0     TO SOLVE A*X = B ,
!!$ C        = NONZERO  TO SOLVE TRANS(A)*X = B ,WHERE
!!$ C              TRANS(A) IS THE TRANSPOSE.
!!$ C
!!$ C   ON RETURN
!!$ C
!!$ C    B    THE SOLUTION VECTOR X .
!!$ C
!!$ C   ERROR CONDITION
!!$ C
!!$ C    A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
!!$ C    ZERO ON THE DIAGONAL. TECHNICALLY THIS INDICATES SINGULARITY
!!$ C    BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
!!$ C    SETTING OF LDA . IT WILL NOT OCCUR IF THE SUBROUTINES ARE
!!$ C    CALLED CORRECTLY AND IF DGBCO HAS SET RCOND .GT. 0.0
!!$ C    OR DGBFA HAS SET INFO .EQ. 0 .
!!$ C
!!$ C   TO COMPUTE INVERSE(A) * C WHERE C IS A MATRIX
!!$ C   WITH P COLUMNS
!!$ C      CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
!!$ C      IF (RCOND IS TOO SMALL) GO TO ...
!!$ C      DO 10 J = 1,P
!!$ C       CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
!!$ C    10 CONTINUE
!!$ C
!!$ C***REFERENCES J. J. DONGARRA,J. R. BUNCH,C. B. MOLER,AND G. W.
!!$ C         STEWART,LINPACK USERS' GUIDE,SIAM,1979.
!!$ C***ROUTINES CALLED DAXPY,DDOT
!!$ C***REVISION HISTORY (YYMMDD)
!!$ C  780814 DATE WRITTEN
!!$ C  890531 CHANGED ALL SPECIFIC INTRINSICS TO GENERIC. (WRB)
!!$ C  890831 MODIFIED ARRAY DECLARATIONS. (WRB)
!!$ C  890831 REVISION DATE FROM VERSION 3.2
!!$ C  891214 PROLOGUE CONVERTED TO VERSION 4.0 FORMAT. (BAB)
!!$ C  900326 REMOVED DUPLICATE INFORMATION FROM DESCRIPTION SECTION.
!!$ C      (WRB)
!!$ C  920501 REFORMATTED THE REFERENCES SECTION. (WRB)
!!$ C***END PROLOGUE DGBSL
  integer::lda,n,ml,mu,ipvt(*),job
  real*8::abd(lda,*),b(*)
!!$ C
  real*8::ddot,t
  integer::k,kb,l,la,lb,lm,m,nm1
!!$ C***FIRST EXECUTABLE STATEMENT DGBSL
  m = mu + ml + 1
  nm1 = n - 1
  if (job .ne. 0) go to 50
!!$ C
!!$ C    JOB = 0 ,SOLVE A * X = B
!!$ C    FIRST SOLVE L*Y = B
!!$ C
  if (ml .eq. 0) go to 30
  if (nm1 .lt. 1) go to 30
  do k = 1,nm1
     lm = min(ml,n-k)
     l = ipvt(k)
     t = b(l)
     if (l .eq. k) go to 10
     b(l) = b(k)
     b(k) = t
10   continue
     call daxpy(lm,t,abd(m+1,k),1,b(k+1),1)
20   continue
  end do
30 continue
!!$ C
!!$ C    NOW SOLVE U*X = Y
!!$ C
  do kb = 1,n
     k = n + 1 - kb
     b(k) = b(k)/abd(m,k)
     lm = min(k,m) - 1
     la = m - lm
     lb = k - lm
     t = -b(k)
     call daxpy(lm,t,abd(la,k),1,b(lb),1)
40   continue
  end do
  go to 100
50 continue
!!$ C
!!$ C    JOB = NONZERO,SOLVE TRANS(A) * X = B
!!$ C    FIRST SOLVE TRANS(U)*Y = B
!!$ C
  do k = 1,n
     lm = min(k,m) - 1
     la = m - lm
     lb = k - lm
     t = ddot(lm,abd(la,k),1,b(lb),1)
     b(k) = (b(k) - t)/abd(m,k)
60   continue
  end do
!!$ C
!!$ C    NOW SOLVE TRANS(L)*X = Y
!!$ C
  if (ml .eq. 0) go to 90
  if (nm1 .lt. 1) go to 90
  do kb = 1,nm1
     k = n - kb
     lm = min(ml,n-k)
     b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1)
     l = ipvt(k)
     if (l .eq. k) go to 70
     t = b(l)
     b(l) = b(k)
     b(k) = t
70   continue
80   continue
  end do
90 continue
100 continue
  return
end subroutine dgbsl

!!$*DECK DAXPY

subroutine daxpy (n,da,dx,incx,dy,incy)
!!$ C***BEGIN PROLOGUE DAXPY
!!$ C***PURPOSE COMPUTE A CONSTANT TIMES A VECTOR PLUS A VECTOR.
!!$ C***CATEGORY D1A7
!!$ C***TYPE   real*8::(SAXPY-S,DAXPY-D,CAXPY-C)
!!$ C***KEYWORDS BLAS,LINEAR ALGEBRA,TRIAD,VECTOR
!!$ C***AUTHOR LAWSON,C. L.,(JPL)
!!$ C      HANSON,R. J.,(SNLA)
!!$ C      KINCAID,D. R.,(U. OF TEXAS)
!!$ C      KROGH,F. T.,(JPL)
!!$ C***DESCRIPTION
!!$ C
!!$ C        B L A S SUBPROGRAM
!!$ C  DESCRIPTION OF PARAMETERS
!!$ C
!!$ C   --INPUT--
!!$ C    N NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!!$ C    DA real*8::SCALAR MULTIPLIER
!!$ C    DX real*8::VECTOR WITH N ELEMENTS
!!$ C   INCX STORAGE SPACING BETWEEN ELEMENTS OF DX
!!$ C    DY real*8::VECTOR WITH N ELEMENTS
!!$ C   INCY STORAGE SPACING BETWEEN ELEMENTS OF DY
!!$ C
!!$ C   --OUTPUT--
!!$ C    DY real*8::RESULT (UNCHANGED IF N .LE. 0)
!!$ C
!!$ C   OVERWRITE real*8::DY WITH real*8::DA*DX + DY.
!!$ C   FOR I = 0 TO N-1,REPLACE DY(LY+I*INCY) WITH DA*DX(LX+I*INCX) +
!!$ C    DY(LY+I*INCY),
!!$ C   WHERE LX = 1 IF INCX .GE. 0,ELSE LX = 1+(1-N)*INCX,AND LY IS
!!$ C   DEFINED IN A SIMILAR WAY USING INCY.
!!$ C
!!$ C***REFERENCES C. L. LAWSON,R. J. HANSON,D. R. KINCAID AND F. T.
!!$ C         KROGH,BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN
!!$ C         USAGE,ALGORITHM NO. 539,TRANSACTIONS ON MATHEMATICAL
!!$ C         SOFTWARE 5,3 (SEPTEMBER 1979),PP. 308-323.
!!$ C***ROUTINES CALLED (NONE)
!!$ C***REVISION HISTORY (YYMMDD)
!!$ C  791001 DATE WRITTEN
!!$ C  890831 MODIFIED ARRAY DECLARATIONS. (WRB)
!!$ C  890831 REVISION DATE FROM VERSION 3.2
!!$ C  891214 PROLOGUE CONVERTED TO VERSION 4.0 FORMAT. (BAB)
!!$ C  920310 CORRECTED DEFINITION OF LX IN DESCRIPTION. (WRB)
!!$ C  920501 REFORMATTED THE REFERENCES SECTION. (WRB)
!!$ C***END PROLOGUE DAXPY
  real*8::dx(*),dy(*),da
  integer::n,incx,incy,ix
  integer::iy,i,m,mp1,ns
!!$ C***FIRST EXECUTABLE STATEMENT DAXPY
  if (n.le.0 .or. da.eq.0.0d0) return
  if (incx .eq. incy) if (incx-1) 5,20,60
!!$ C
!!$ C   CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!!$ C
5 ix = 1
  iy = 1
  if (incx .lt. 0) ix = (-n+1)*incx + 1
  if (incy .lt. 0) iy = (-n+1)*incy + 1
  do i = 1,n
     dy(iy) = dy(iy) + da*dx(ix)
     ix = ix + incx
     iy = iy + incy
10   continue
  end do
  return
!!$ C
!!$ C   CODE FOR BOTH INCREMENTS EQUAL TO 1.
!!$ C
!!$ C   CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
!!$ C
20 m = mod(n,4)
  if (m .eq. 0) go to 40
  do i = 1,m
     dy(i) = dy(i) + da*dx(i)
30   continue
  end do
  if (n .lt. 4) return
40 mp1 = m + 1
  do i = mp1,n,4
     dy(i) = dy(i) + da*dx(i)
     dy(i+1) = dy(i+1) + da*dx(i+1)
     dy(i+2) = dy(i+2) + da*dx(i+2)
     dy(i+3) = dy(i+3) + da*dx(i+3)
50   continue
  end do
  return
!!$ C
!!$ C   CODE FOR EQUAL,POSITIVE,NON-UNIT INCREMENTS.
!!$ C
60 ns = n*incx
  do i = 1,ns,incx
     dy(i) = da*dx(i) + dy(i)
70   continue
  end do
  return
end subroutine daxpy

!!$*DECK DCOPY

subroutine dcopy (n,dx,incx,dy,incy)
!!$ C***BEGIN PROLOGUE DCOPY
!!$ C***PURPOSE COPY A VECTOR.
!!$ C***CATEGORY D1A5
!!$ C***TYPE   real*8::(SCOPY-S,DCOPY-D,CCOPY-C,ICOPY-I)
!!$ C***KEYWORDS BLAS,COPY,LINEAR ALGEBRA,VECTOR
!!$ C***AUTHOR LAWSON,C. L.,(JPL)
!!$ C      HANSON,R. J.,(SNLA)
!!$ C      KINCAID,D. R.,(U. OF TEXAS)
!!$ C      KROGH,F. T.,(JPL)
!!$ C***DESCRIPTION
!!$ C
!!$ C        B L A S SUBPROGRAM
!!$ C  DESCRIPTION OF PARAMETERS
!!$ C
!!$ C   --INPUT--
!!$ C    N NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!!$ C    DX real*8::VECTOR WITH N ELEMENTS
!!$ C   INCX STORAGE SPACING BETWEEN ELEMENTS OF DX
!!$ C    DY real*8::VECTOR WITH N ELEMENTS
!!$ C   INCY STORAGE SPACING BETWEEN ELEMENTS OF DY
!!$ C
!!$ C   --OUTPUT--
!!$ C    DY COPY OF VECTOR DX (UNCHANGED IF N .LE. 0)
!!$ C
!!$ C   COPY real*8::DX TO real*8::DY.
!!$ C   FOR I = 0 TO N-1,COPY DX(LX+I*INCX) TO DY(LY+I*INCY),
!!$ C   WHERE LX = 1 IF INCX .GE. 0,ELSE LX = 1+(1-N)*INCX,AND LY IS
!!$ C   DEFINED IN A SIMILAR WAY USING INCY.
!!$ C
!!$ C***REFERENCES C. L. LAWSON,R. J. HANSON,D. R. KINCAID AND F. T.
!!$ C         KROGH,BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN
!!$ C         USAGE,ALGORITHM NO. 539,TRANSACTIONS ON MATHEMATICAL
!!$ C         SOFTWARE 5,3 (SEPTEMBER 1979),PP. 308-323.
!!$ C***ROUTINES CALLED (NONE)
!!$ C***REVISION HISTORY (YYMMDD)
!!$ C  791001 DATE WRITTEN
!!$ C  890831 MODIFIED ARRAY DECLARATIONS. (WRB)
!!$ C  890831 REVISION DATE FROM VERSION 3.2
!!$ C  891214 PROLOGUE CONVERTED TO VERSION 4.0 FORMAT. (BAB)
!!$ C  920310 CORRECTED DEFINITION OF LX IN DESCRIPTION. (WRB)
!!$ C  920501 REFORMATTED THE REFERENCES SECTION. (WRB)
!!$ C***END PROLOGUE DCOPY
  real*8::dx(*),dy(*)
  integer::n,incx,incy,ix,iy
  integer::i,m,mp1,nd,ns
!!$ C***FIRST EXECUTABLE STATEMENT DCOPY
  if (n .le. 0) return
  if (incx .eq. incy) if (incx-1) 5,20,60
!!$ C
!!$ C   CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!!$ C
5 ix = 1
  iy = 1
  if (incx .lt. 0) ix = (-n+1)*incx + 1
  if (incy .lt. 0) iy = (-n+1)*incy + 1
  do i = 1,n
     dy(iy) = dx(ix)
     ix = ix + incx
     iy = iy + incy
10   continue
  end do
  return
!!$ C
!!$ C   CODE FOR BOTH INCREMENTS EQUAL TO 1.
!!$ C
!!$ C   CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.
!!$ C
20 m = mod(n,7)
  if (m .eq. 0) go to 40
  do i = 1,m
     dy(i) = dx(i)
30   continue
  end do
  if (n .lt. 7) return
40 mp1 = m + 1
  do i = mp1,n,7
     dy(i) = dx(i)
     dy(i+1) = dx(i+1)
     dy(i+2) = dx(i+2)
     dy(i+3) = dx(i+3)
     dy(i+4) = dx(i+4)
     dy(i+5) = dx(i+5)
     dy(i+6) = dx(i+6)
50   continue
  end do
  return
!!$ C
!!$ C   CODE FOR EQUAL,POSITIVE,NON-UNIT INCREMENTS.
!!$ C
60 ns = n*incx
  do i = 1,ns,incx
     dy(i) = dx(i)
70   continue
  end do
  return
end subroutine dcopy

!!$*DECK DDOT

function ddot (n,dx,incx,dy,incy)
  real*8::ddot
!!$ C***BEGIN PROLOGUE DDOT
!!$ C***PURPOSE COMPUTE THE INNER PRODUCT OF TWO VECTORS.
!!$ C***CATEGORY D1A4
!!$ C***TYPE   real*8::(SDOT-S,DDOT-D,CDOTU-C)
!!$ C***KEYWORDS BLAS,INNER PRODUCT,LINEAR ALGEBRA,VECTOR
!!$ C***AUTHOR LAWSON,C. L.,(JPL)
!!$ C      HANSON,R. J.,(SNLA)
!!$ C      KINCAID,D. R.,(U. OF TEXAS)
!!$ C      KROGH,F. T.,(JPL)
!!$ C***DESCRIPTION
!!$ C
!!$ C        B L A S SUBPROGRAM
!!$ C  DESCRIPTION OF PARAMETERS
!!$ C
!!$ C   --INPUT--
!!$ C    N NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!!$ C    DX real*8::VECTOR WITH N ELEMENTS
!!$ C   INCX STORAGE SPACING BETWEEN ELEMENTS OF DX
!!$ C    DY real*8::VECTOR WITH N ELEMENTS
!!$ C   INCY STORAGE SPACING BETWEEN ELEMENTS OF DY
!!$ C
!!$ C   --OUTPUT--
!!$ C   DDOT real*8::DOT PRODUCT (ZERO IF N .LE. 0)
!!$ C
!!$ C   RETURNS THE DOT PRODUCT OF real*8::DX AND DY.
!!$ C   DDOT = SUM FOR I = 0 TO N-1 OF DX(LX+I*INCX) * DY(LY+I*INCY),
!!$ C   WHERE LX = 1 IF INCX .GE. 0,ELSE LX = 1+(1-N)*INCX,AND LY IS
!!$ C   DEFINED IN A SIMILAR WAY USING INCY.
!!$ C
!!$ C***REFERENCES C. L. LAWSON,R. J. HANSON,D. R. KINCAID AND F. T.
!!$ C         KROGH,BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN
!!$ C         USAGE,ALGORITHM NO. 539,TRANSACTIONS ON MATHEMATICAL
!!$ C         SOFTWARE 5,3 (SEPTEMBER 1979),PP. 308-323.
!!$ C***ROUTINES CALLED (NONE)
!!$ C***REVISION HISTORY (YYMMDD)
!!$ C  791001 DATE WRITTEN
!!$ C  890831 MODIFIED ARRAY DECLARATIONS. (WRB)
!!$ C  890831 REVISION DATE FROM VERSION 3.2
!!$ C  891214 PROLOGUE CONVERTED TO VERSION 4.0 FORMAT. (BAB)
!!$ C  920310 CORRECTED DEFINITION OF LX IN DESCRIPTION. (WRB)
!!$ C  920501 REFORMATTED THE REFERENCES SECTION. (WRB)
!!$ C***END PROLOGUE DDOT
  real*8::dx(*),dy(*)
  integer::n,incx,incy,ix,iy
  integer::i,m,mp1,ns
!!$ C***FIRST EXECUTABLE STATEMENT DDOT
  ddot = 0.0d0
  if (n .le. 0) return
  if (incx .eq. incy) if (incx-1) 5,20,60
!!$ C
!!$ C   CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!!$ C
5 ix = 1
  iy = 1
  if (incx .lt. 0) ix = (-n+1)*incx + 1
  if (incy .lt. 0) iy = (-n+1)*incy + 1
  do i = 1,n
     ddot = ddot + dx(ix)*dy(iy)
     ix = ix + incx
     iy = iy + incy
10   continue
  end do
  return
!!$ C
!!$ C   CODE FOR BOTH INCREMENTS EQUAL TO 1.
!!$ C
!!$ C   CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
!!$ C
20 m = mod(n,5)
  if (m .eq. 0) go to 40
  do i = 1,m
     ddot = ddot + dx(i)*dy(i)
30   continue
  end do
  if (n .lt. 5) return
40 mp1 = m + 1
  do i = mp1,n,5
     ddot = ddot + dx(i)*dy(i) + dx(i+1)*dy(i+1) + dx(i+2)*dy(i+2) +&
          dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
50   continue
  end do
  return
!!$ C
!!$ C   CODE FOR EQUAL,POSITIVE,NON-UNIT INCREMENTS.
!!$ C
60 ns = n*incx
  do i = 1,ns,incx
     ddot = ddot + dx(i)*dy(i)
70   continue
  end do
  return
end function ddot

!!****************************************
!!$*DECK DNRM2

function dnrm2 (n,dx,incx)
  real*8::dnrm2
!!$ C***BEGIN PROLOGUE DNRM2
!!$ C***PURPOSE COMPUTE THE EUCLIDEAN LENGTH (L2 NORM) OF A VECTOR.
!!$ C***CATEGORY D1A3B
!!$ C***TYPE   real*8::(SNRM2-S,DNRM2-D,SCNRM2-C)
!!$ C***KEYWORDS BLAS,EUCLIDEAN LENGTH,EUCLIDEAN NORM,L2,
!!$ C       LINEAR ALGEBRA,UNITARY,VECTOR
!!$ C***AUTHOR LAWSON,C. L.,(JPL)
!!$ C      HANSON,R. J.,(SNLA)
!!$ C      KINCAID,D. R.,(U. OF TEXAS)
!!$ C      KROGH,F. T.,(JPL)
!!$ C***DESCRIPTION
!!$ C
!!$ C        B L A S SUBPROGRAM
!!$ C  DESCRIPTION OF PARAMETERS
!!$ C
!!$ C   --INPUT--
!!$ C    N NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!!$ C    DX real*8::VECTOR WITH N ELEMENTS
!!$ C   INCX STORAGE SPACING BETWEEN ELEMENTS OF DX
!!$ C
!!$ C   --OUTPUT--
!!$ C  DNRM2 real*8::RESULT (ZERO IF N .LE. 0)
!!$ C
!!$ C   EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX WITH STORAGE
!!$ C   INCREMENT INCX.
!!$ C   IF N .LE. 0,RETURN WITH RESULT = 0.
!!$ C   IF N .GE. 1,THEN INCX MUST BE .GE. 1
!!$ C
!!$ C   FOUR PHASE METHOD USING TWO BUILT-IN CONSTANTS THAT ARE
!!$ C   HOPEFULLY APPLICABLE TO ALL MACHINES.
!!$ C     CUTLO = MAXIMUM OF SQRT(U/EPS) OVER ALL KNOWN MACHINES.
!!$ C     CUTHI = MINIMUM OF SQRT(V)   OVER ALL KNOWN MACHINES.
!!$ C   WHERE
!!$ C     EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
!!$ C     U  = SMALLEST POSITIVE NO.  (UNDERFLOW LIMIT)
!!$ C     V  = LARGEST NO.      (OVERFLOW LIMIT)
!!$ C
!!$ C   BRIEF OUTLINE OF ALGORITHM.
!!$ C
!!$ C   PHASE 1 SCANS ZERO COMPONENTS.
!!$ C   MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
!!$ C   MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
!!$ C   MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
!!$ C   WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
!!$ C
!!$ C   VALUES FOR CUTLO AND CUTHI.
!!$ C   FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
!!$ C   DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS:
!!$ C   CUTLO,S.P.  U/EPS = 2**(-102) FOR HONEYWELL. CLOSE SECONDS ARE
!!$ C          UNIVAC AND DEC AT 2**(-103)
!!$ C          THUS CUTLO = 2**(-51) = 4.44089E-16
!!$ C   CUTHI,S.P.  V = 2**127 FOR UNIVAC,HONEYWELL,AND DEC.
!!$ C          THUS CUTHI = 2**(63.5) = 1.30438E19
!!$ C   CUTLO,D.P.  U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
!!$ C          THUS CUTLO = 2**(-33.5) = 8.23181D-11
!!$ C   CUTHI,D.P.  SAME AS S.P. CUTHI = 1.30438D19
!!$ C   DATA CUTLO,CUTHI /8.232D-11,1.304D19/
!!$ C   DATA CUTLO,CUTHI /4.441E-16,1.304E19/
!!$ C
!!$ C***REFERENCES C. L. LAWSON,R. J. HANSON,D. R. KINCAID AND F. T.
!!$ C         KROGH,BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN
!!$ C         USAGE,ALGORITHM NO. 539,TRANSACTIONS ON MATHEMATICAL
!!$ C         SOFTWARE 5,3 (SEPTEMBER 1979),PP. 308-323.
!!$ C***ROUTINES CALLED (NONE)
!!$ C***REVISION HISTORY (YYMMDD)
!!$ C  791001 DATE WRITTEN
!!$ C  890531 CHANGED ALL SPECIFIC INTRINSICS TO GENERIC. (WRB)
!!$ C  890831 MODIFIED ARRAY DECLARATIONS. (WRB)
!!$ C  890831 REVISION DATE FROM VERSION 3.2
!!$ C  891214 PROLOGUE CONVERTED TO VERSION 4.0 FORMAT. (BAB)
!!$ C  920501 REFORMATTED THE REFERENCES SECTION. (WRB)
!!$ C***END PROLOGUE DNRM2
  integer::next
  real*8::dx(*),cutlo,cuthi,hitest,sum,xmax,zero,        one
  integer::n,nn,incx,i,j
  save cutlo,cuthi,zero,one
  data zero,one /0.0d0,1.0d0/
!!$ C
  data cutlo,cuthi /8.232d-11,1.304d19/
!!$ C***FIRST EXECUTABLE STATEMENT DNRM2
  if (n .le. 0) then
     dnrm2 = zero
     go to 300
  end if
!!$ C
  next = 30 !###@@
  assign 30 to next
  sum = zero
  nn = n * incx
!!$ C
!!$ C                         BEGIN MAIN LOOP
!!$ C
  i = 1
20 go to next,(30,50,70,110)
30 if (abs(dx(i)) .gt. cutlo) go to 85
  assign 50 to next
  xmax = zero
!!$ C
!!$ C            PHASE 1. SUM IS ZERO
!!$ C
50 if (dx(i) .eq. zero) go to 200
  if (abs(dx(i)) .gt. cutlo) go to 85
!!$ C
!!$ C                PREPARE FOR PHASE 2.
!!$ C
  assign 70 to next
  go to 105
!!$ C
!!$ C                PREPARE FOR PHASE 4.
!!$ C
100 i = j
  assign 110 to next
  sum = (sum / dx(i)) / dx(i)
105 xmax = abs(dx(i))
  go to 115
!!$ C
!!$ C          PHASE 2. SUM IS SMALL.
!!$ C               SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
!!$ C
70 if (abs(dx(i)) .gt. cutlo) go to 75
!!$ C
!!$ C           COMMON CODE FOR PHASES 2 AND 4.
!!$ C           IN PHASE 4 SUM IS LARGE. SCALE TO AVOID OVERFLOW.
!!$ C
110 if (abs(dx(i)) .le. xmax) go to 115
  sum = one + sum * (xmax / dx(i))**2
  xmax = abs(dx(i))
  go to 200
!!$ C
115 sum = sum + (dx(i)/xmax)**2
  go to 200
!!$ C
!!$ C         PREPARE FOR PHASE 3.
!!$ C
75 sum = (sum * xmax) * xmax
!!$ C
!!$ C   FOR REAL OR D.P. SET HITEST = CUTHI/N
!!$ C   FOR COMPLEX   SET HITEST = CUTHI/(2*N)
!!$ C
85 hitest = cuthi / n
!!$ C
!!$ C          PHASE 3. SUM IS MID-RANGE. NO SCALING.
!!$ C
  do j = i,nn,incx
     if (abs(dx(j)) .ge. hitest) go to 100
     sum = sum + dx(j)**2
  end do
  dnrm2 = sqrt(sum)
  go to 300
!!$ C
200 continue
  i = i + incx
  if (i .le. nn) go to 20
!!$ C
!!$ C       END OF MAIN LOOP.
!!$ C
!!$ C       COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
!!$ C
  dnrm2 = xmax * sqrt(sum)
300 continue
  return
end function dnrm2

!!$*DECK DSCAL

subroutine dscal (n,da,dx,incx)
!!$ C***BEGIN PROLOGUE DSCAL
!!$ C***PURPOSE MULTIPLY A VECTOR BY A CONSTANT.
!!$ C***CATEGORY D1A6
!!$ C***TYPE   real*8::(SSCAL-S,DSCAL-D,CSCAL-C)
!!$ C***KEYWORDS BLAS,LINEAR ALGEBRA,SCALE,VECTOR
!!$ C***AUTHOR LAWSON,C. L.,(JPL)
!!$ C      HANSON,R. J.,(SNLA)
!!$ C      KINCAID,D. R.,(U. OF TEXAS)
!!$ C      KROGH,F. T.,(JPL)
!!$ C***DESCRIPTION
!!$ C
!!$ C        B L A S SUBPROGRAM
!!$ C  DESCRIPTION OF PARAMETERS
!!$ C
!!$ C   --INPUT--
!!$ C    N NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!!$ C    DA real*8::SCALE FACTOR
!!$ C    DX real*8::VECTOR WITH N ELEMENTS
!!$ C   INCX STORAGE SPACING BETWEEN ELEMENTS OF DX
!!$ C
!!$ C   --OUTPUT--
!!$ C    DX real*8::RESULT (UNCHANGED IF N.LE.0)
!!$ C
!!$ C   REPLACE real*8::DX BY real*8::DA*DX.
!!$ C   FOR I = 0 TO N-1,REPLACE DX(IX+I*INCX) WITH DA * DX(IX+I*INCX),
!!$ C   WHERE IX = 1 IF INCX .GE. 0,ELSE IX = 1+(1-N)*INCX.
!!$ C
!!$ C***REFERENCES C. L. LAWSON,R. J. HANSON,D. R. KINCAID AND F. T.
!!$ C         KROGH,BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN
!!$ C         USAGE,ALGORITHM NO. 539,TRANSACTIONS ON MATHEMATICAL
!!$ C         SOFTWARE 5,3 (SEPTEMBER 1979),PP. 308-323.
!!$ C***ROUTINES CALLED (NONE)
!!$ C***REVISION HISTORY (YYMMDD)
!!$ C  791001 DATE WRITTEN
!!$ C  890831 MODIFIED ARRAY DECLARATIONS. (WRB)
!!$ C  890831 REVISION DATE FROM VERSION 3.2
!!$ C  891214 PROLOGUE CONVERTED TO VERSION 4.0 FORMAT. (BAB)
!!$ C  900821 MODIFIED TO CORRECT PROBLEM WITH A NEGATIVE INCREMENT.
!!$ C      (WRB)
!!$ C  920501 REFORMATTED THE REFERENCES SECTION. (WRB)
!!$ C***END PROLOGUE DSCAL
  real*8::da,dx(*)
  integer::i,incx,ix,m,mp1,n
!!$ C***FIRST EXECUTABLE STATEMENT DSCAL
  if (n .le. 0) return
  if (incx .eq. 1) goto 20
!!$ C
!!$ C   CODE FOR INCREMENT NOT EQUAL TO 1.
!!$ C
  ix = 1
  if (incx .lt. 0) ix = (-n+1)*incx + 1
  do i = 1,n
     dx(ix) = da*dx(ix)
     ix = ix + incx
10   continue
  end do
  return
!!$ C
!!$ C   CODE FOR INCREMENT EQUAL TO 1.
!!$ C
!!$ C   CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
!!$ C
20 m = mod(n,5)
  if (m .eq. 0) goto 40
  do i = 1,m
     dx(i) = da*dx(i)
30   continue
  end do
  if (n .lt. 5) return
40 mp1 = m + 1
  do i = mp1,n,5
     dx(i) = da*dx(i)
     dx(i+1) = da*dx(i+1)
     dx(i+2) = da*dx(i+2)
     dx(i+3) = da*dx(i+3)
     dx(i+4) = da*dx(i+4)
50   continue
  end do
  return
end subroutine dscal

!!$*DECK IDAMAX

function idamax (n,dx,incx)
  integer::idamax
!!$ C***BEGIN PROLOGUE IDAMAX
!!$ C***PURPOSE FIND THE SMALLEST INDEX OF THAT COMPONENT OF A VECTOR
!!$ C      HAVING THE MAXIMUM MAGNITUDE.
!!$ C***CATEGORY D1A2
!!$ C***TYPE   real*8::(ISAMAX-S,IDAMAX-D,ICAMAX-C)
!!$ C***KEYWORDS BLAS,LINEAR ALGEBRA,MAXIMUM COMPONENT,VECTOR
!!$ C***AUTHOR LAWSON,C. L.,(JPL)
!!$ C      HANSON,R. J.,(SNLA)
!!$ C      KINCAID,D. R.,(U. OF TEXAS)
!!$ C      KROGH,F. T.,(JPL)
!!$ C***DESCRIPTION
!!$ C
!!$ C        B L A S SUBPROGRAM
!!$ C  DESCRIPTION OF PARAMETERS
!!$ C
!!$ C   --INPUT--
!!$ C    N NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!!$ C    DX real*8::VECTOR WITH N ELEMENTS
!!$ C   INCX STORAGE SPACING BETWEEN ELEMENTS OF DX
!!$ C
!!$ C   --OUTPUT--
!!$ C  IDAMAX SMALLEST INDEX (ZERO IF N .LE. 0)
!!$ C
!!$ C   FIND SMALLEST INDEX OF MAXIMUM MAGNITUDE OF real*8::DX.
!!$ C   IDAMAX = FIRST I,I = 1 TO N,TO MAXIMIZE ABS(DX(IX+(I-1)*INCX)),
!!$ C   WHERE IX = 1 IF INCX .GE. 0,ELSE IX = 1+(1-N)*INCX.
!!$ C
!!$ C***REFERENCES C. L. LAWSON,R. J. HANSON,D. R. KINCAID AND F. T.
!!$ C         KROGH,BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN
!!$ C         USAGE,ALGORITHM NO. 539,TRANSACTIONS ON MATHEMATICAL
!!$ C         SOFTWARE 5,3 (SEPTEMBER 1979),PP. 308-323.
!!$ C***ROUTINES CALLED (NONE)
!!$ C***REVISION HISTORY (YYMMDD)
!!$ C  791001 DATE WRITTEN
!!$ C  890531 CHANGED ALL SPECIFIC INTRINSICS TO GENERIC. (WRB)
!!$ C  890531 REVISION DATE FROM VERSION 3.2
!!$ C  891214 PROLOGUE CONVERTED TO VERSION 4.0 FORMAT. (BAB)
!!$ C  900821 MODIFIED TO CORRECT PROBLEM WITH A NEGATIVE INCREMENT.
!!$ C      (WRB)
!!$ C  920501 REFORMATTED THE REFERENCES SECTION. (WRB)
!!$ C***END PROLOGUE IDAMAX
  real*8::dx(*),dmax,xmag
  integer::i,incx,ix,n
!!$ C***FIRST EXECUTABLE STATEMENT IDAMAX
  idamax = 0
  if (n .le. 0) return
  idamax = 1
  if (n .eq. 1) return
!!$ C
  if (incx .eq. 1) goto 20
!!$ C
!!$ C   CODE FOR INCREMENTS NOT EQUAL TO 1.
!!$ C
  ix = 1
  if (incx .lt. 0) ix = (-n+1)*incx + 1
  dmax = abs(dx(ix))
  ix = ix + incx
  do i = 2,n
     xmag = abs(dx(ix))
     if (xmag .gt. dmax) then
        idamax = i
        dmax = xmag
     endif
     ix = ix + incx
10   continue
  end do
  return
!!$ C
!!$ C   CODE FOR INCREMENTS EQUAL TO 1.
!!$ C
20 dmax = abs(dx(1))
  do i = 2,n
     xmag = abs(dx(i))
     if (xmag .gt. dmax) then
        idamax = i
        dmax = xmag
     endif
30   continue
  end do
  return
end function idamax

!!$*DECK XERRWD

subroutine xerrwd (msg,nmes,nerr,level,ni,i1,i2,nr,r1,r2)
!!$ C***BEGIN PROLOGUE XERRWD
!!$ C***SUBSIDIARY
!!$ C***PURPOSE WRITE ERROR MESSAGE WITH VALUES.
!!$ C***CATEGORY R3C
!!$ C***TYPE   real*8::(XERRWV-S,XERRWD-D)
!!$ C***AUTHOR HINDMARSH,ALAN C.,(LLNL)
!!$ C***DESCRIPTION
!!$ C
!!$ C SUBROUTINES XERRWD,XSETF,XSETUN,AND THE FUNCTION ROUTINE IXSAV,
!!$ C AS GIVEN HERE,CONSTITUTE A SIMPLIFIED VERSION OF THE SLATEC ERROR
!!$ C HANDLING PACKAGE.
!!$ C
!!$ C ALL ARGUMENTS ARE INPUT ARGUMENTS.
!!$ C
!!$ C MSG  = THE MESSAGE (CHARACTER ARRAY).
!!$ C NMES  = THE LENGTH OF MSG (NUMBER OF CHARACTERS).
!!$ C NERR  = THE ERROR NUMBER (NOT USED).
!!$ C LEVEL = THE ERROR LEVEL..
!!$ C      0 OR 1 MEANS RECOVERABLE (CONTROL RETURNS TO CALLER).
!!$ C      2 MEANS FATAL (RUN IS ABORTED--SEE NOTE BELOW).
!!$ C NI   = NUMBER OF INTEGERS (0,1,OR 2) TO BE PRINTED WITH MESSAGE.
!!$ C I1,I2 = INTEGERS TO BE PRINTED,DEPENDING ON NI.
!!$ C NR   = NUMBER OF REALS (0,1,OR 2) TO BE PRINTED WITH MESSAGE.
!!$ C R1,R2 = REALS TO BE PRINTED,DEPENDING ON NR.
!!$ C
!!$ C NOTE.. THIS ROUTINE IS MACHINE-DEPENDENT AND SPECIALIZED FOR USE
!!$ C IN LIMITED CONTEXT,IN THE FOLLOWING WAYS..
!!$ C 1. THE ARGUMENT MSG IS ASSUMED TO BE OF TYPE CHARACTER,AND
!!$ C   THE MESSAGE IS PRINTED WITH A FORMAT OF (1X,A).
!!$ C 2. THE MESSAGE IS ASSUMED TO TAKE ONLY ONE LINE.
!!$ C   MULTI-LINE MESSAGES ARE GENERATED BY REPEATED CALLS.
!!$ C 3. IF LEVEL = 2,CONTROL PASSES TO THE STATEMENT  STOP
!!$ C   TO ABORT THE RUN. THIS STATEMENT MAY BE MACHINE-DEPENDENT.
!!$ C 4. R1 AND R2 ARE ASSUMED TO BE IN real*8::AND ARE PRINTED
!!$ C   IN D21.13 FORMAT.
!!$ C
!!$ C***ROUTINES CALLED IXSAV
!!$ C***REVISION HISTORY (YYMMDD)
!!$ C  920831 DATE WRITTEN
!!$ C  921118 REPLACED MFLGSV/LUNSAV BY IXSAV. (ACH)
!!$ C  930329 MODIFIED PROLOGUE TO SLATEC FORMAT. (FNF)
!!$ C  930407 CHANGED MSG FROM CHARACTER*1 ARRAY TO VARIABLE. (FNF)
!!$ C  930922 MINOR COSMETIC CHANGE. (FNF)
!!$ C***END PROLOGUE XERRWD
!!$ C
!!$ C*INTERNAL NOTES:
!!$ C
!!$ C FOR A DIFFERENT DEFAULT LOGICAL UNIT NUMBER,IXSAV (OR A SUBSIDIARY
!!$ C ROUTINE THAT IT CALLS) WILL NEED TO BE MODIFIED.
!!$ C FOR A DIFFERENT RUN-ABORT COMMAND,CHANGE THE STATEMENT FOLLOWING
!!$ C STATEMENT 100 AT THE END.
!!$ C-----------------------------------------------------------------------
!!$ C SUBROUTINES CALLED BY XERRWD.. NONE
!!$ C FUNCTION ROUTINE CALLED BY XERRWD.. IXSAV
!!$ C-----------------------------------------------------------------------
!!$ C**END
!!$ C
!!$ C DECLARE ARGUMENTS.
!!$ C
  real*8::r1,r2
  integer::nmes,nerr,level,ni,i1,i2,nr
  character*(*) msg
!!$ C
!!$ C DECLARE LOCAL VARIABLES.
!!$ C
  integer::lunit,ixsav,mesflg
!!$ C
!!$ C GET LOGICAL UNIT NUMBER AND MESSAGE PRINT FLAG.
!!$ C
!!$ C***FIRST EXECUTABLE STATEMENT XERRWD
  lunit = ixsav (1,0,.false.)
  mesflg = ixsav (2,0,.false.)
  if (mesflg .eq. 0) go to 100
!!$ C
!!$ C WRITE THE MESSAGE.
!!$ C
  write (lunit,10) msg
10 format(1x,a)
  if (ni .eq. 1) write (lunit,20) i1
20 format(6x,'in above message,i1 =',i10)
  if (ni .eq. 2) write (lunit,30) i1,i2
30 format(6x,'in above message,i1 =',i10,3x,'i2 =',i10)
  if (nr .eq. 1) write (lunit,40) r1
40 format(6x,'in above message,r1 =',d21.13)
  if (nr .eq. 2) write (lunit,50) r1,r2
50 format(6x,'in above,r1 =',d21.13,3x,'r2 =',d21.13)
!!$ C
!!$ C ABORT THE RUN IF LEVEL = 2.
!!$ C
100 if (level .ne. 2) return
  stop
!!$ C----------------------- END OF SUBROUTINE XERRWD ----------------------
end subroutine xerrwd

!!$*DECK XSETF

subroutine xsetf (mflag)
!!$ C***BEGIN PROLOGUE XSETF
!!$ C***PURPOSE RESET THE ERROR PRINT CONTROL FLAG.
!!$ C***CATEGORY R3A
!!$ C***TYPE   ALL (XSETF-A)
!!$ C***KEYWORDS ERROR CONTROL
!!$ C***AUTHOR HINDMARSH,ALAN C.,(LLNL)
!!$ C***DESCRIPTION
!!$ C
!!$ C  XSETF SETS THE ERROR PRINT CONTROL FLAG TO MFLAG:
!!$ C   MFLAG=1 MEANS PRINT ALL MESSAGES (THE DEFAULT).
!!$ C   MFLAG=0 MEANS NO PRINTING.
!!$ C
!!$ C***SEE ALSO XERRWD,XERRWV
!!$ C***REFERENCES (NONE)
!!$ C***ROUTINES CALLED IXSAV
!!$ C***REVISION HISTORY (YYMMDD)
!!$ C  921118 DATE WRITTEN
!!$ C  930329 ADDED SLATEC FORMAT PROLOGUE. (FNF)
!!$ C  930407 CORRECTED SEE ALSO SECTION. (FNF)
!!$ C  930922 MADE USER-CALLABLE,AND OTHER COSMETIC CHANGES. (FNF)
!!$ C***END PROLOGUE XSETF
!!$ C
!!$ C SUBROUTINES CALLED BY XSETF.. NONE
!!$ C FUNCTION ROUTINE CALLED BY XSETF.. IXSAV
!!$ C-----------------------------------------------------------------------
!!$ C**END
  integer::mflag,junk,ixsav
!!$ C
!!$ C***FIRST EXECUTABLE STATEMENT XSETF
  if (mflag .eq. 0 .or. mflag .eq. 1) junk = ixsav (2,mflag,.true.)
  return
!!$ C----------------------- END OF SUBROUTINE XSETF -----------------------
end subroutine xsetf

!!$*DECK XSETUN

subroutine xsetun (lun)
!!$ C***BEGIN PROLOGUE XSETUN
!!$ C***PURPOSE RESET THE LOGICAL UNIT NUMBER FOR ERROR MESSAGES.
!!$ C***CATEGORY R3B
!!$ C***TYPE   ALL (XSETUN-A)
!!$ C***KEYWORDS ERROR CONTROL
!!$ C***DESCRIPTION
!!$ C
!!$ C  XSETUN SETS THE LOGICAL UNIT NUMBER FOR ERROR MESSAGES TO LUN.
!!$ C
!!$ C***AUTHOR HINDMARSH,ALAN C.,(LLNL)
!!$ C***SEE ALSO XERRWD,XERRWV
!!$ C***REFERENCES (NONE)
!!$ C***ROUTINES CALLED IXSAV
!!$ C***REVISION HISTORY (YYMMDD)
!!$ C  921118 DATE WRITTEN
!!$ C  930329 ADDED SLATEC FORMAT PROLOGUE. (FNF)
!!$ C  930407 CORRECTED SEE ALSO SECTION. (FNF)
!!$ C  930922 MADE USER-CALLABLE,AND OTHER COSMETIC CHANGES. (FNF)
!!$ C***END PROLOGUE XSETUN
!!$ C
!!$ C SUBROUTINES CALLED BY XSETUN.. NONE
!!$ C FUNCTION ROUTINE CALLED BY XSETUN.. IXSAV
!!$ C-----------------------------------------------------------------------
!!$ C**END
  integer::lun,junk,ixsav
!!$ C
!!$ C***FIRST EXECUTABLE STATEMENT XSETUN
  if (lun .gt. 0) junk = ixsav (1,lun,.true.)
  return
!!$ C----------------------- END OF SUBROUTINE XSETUN ----------------------
end subroutine xsetun

!!$*DECK IXSAV

function ixsav (ipar,ivalue,iset)
  integer::ixsav
!!$ C***BEGIN PROLOGUE IXSAV
!!$ C***SUBSIDIARY
!!$ C***PURPOSE SAVE AND RECALL ERROR MESSAGE CONTROL PARAMETERS.
!!$ C***CATEGORY R3C
!!$ C***TYPE   ALL (IXSAV-A)
!!$ C***AUTHOR HINDMARSH,ALAN C.,(LLNL)
!!$ C***DESCRIPTION
!!$ C
!!$ C IXSAV SAVES AND RECALLS ONE OF TWO ERROR MESSAGE PARAMETERS:
!!$ C  LUNIT,THE LOGICAL UNIT NUMBER TO WHICH MESSAGES ARE PRINTED,AND
!!$ C  MESFLG,THE MESSAGE PRINT FLAG.
!!$ C THIS IS A MODIFICATION OF THE SLATEC LIBRARY ROUTINE J4SAVE.
!!$ C
!!$ C SAVED LOCAL VARIABLES..
!!$ C  LUNIT = LOGICAL UNIT NUMBER FOR MESSAGES. THE DEFAULT IS OBTAINED
!!$ C      BY A CALL TO IUMACH (MAY BE MACHINE-DEPENDENT).
!!$ C  MESFLG = PRINT CONTROL FLAG..
!!$ C      1 MEANS PRINT ALL MESSAGES (THE DEFAULT).
!!$ C      0 MEANS NO PRINTING.
!!$ C
!!$ C ON INPUT..
!!$ C  IPAR  = PARAMETER INDICATOR (1 FOR LUNIT,2 FOR MESFLG).
!!$ C  IVALUE = THE VALUE TO BE SET FOR THE PARAMETER,IF ISET = .TRUE.
!!$ C  ISET  = LOGICAL FLAG TO INDICATE WHETHER TO READ OR WRITE.
!!$ C       IF ISET = .TRUE.,THE PARAMETER WILL BE GIVEN
!!$ C       THE VALUE IVALUE. IF ISET = .FALSE.,THE PARAMETER
!!$ C       WILL BE UNCHANGED,AND IVALUE IS A DUMMY ARGUMENT.
!!$ C
!!$ C ON RETURN..
!!$ C  IXSAV = THE (OLD) VALUE OF THE PARAMETER.
!!$ C
!!$ C***SEE ALSO XERRWD,XERRWV
!!$ C***ROUTINES CALLED IUMACH
!!$ C***REVISION HISTORY (YYMMDD)
!!$ C  921118 DATE WRITTEN
!!$ C  930329 MODIFIED PROLOGUE TO SLATEC FORMAT. (FNF)
!!$ C  930915 ADDED IUMACH CALL TO GET DEFAULT OUTPUT UNIT. (ACH)
!!$ C  930922 MINOR COSMETIC CHANGES. (FNF)
!!$ C  010425 TYPE DECLARATION FOR IUMACH ADDED. (ACH)
!!$ C***END PROLOGUE IXSAV
!!$ C
!!$ C SUBROUTINES CALLED BY IXSAV.. NONE
!!$ C FUNCTION ROUTINE CALLED BY IXSAV.. IUMACH
!!$ C-----------------------------------------------------------------------
!!$ C**END
  logical iset
  integer::ipar,ivalue
!!$ C-----------------------------------------------------------------------
  integer::iumach,lunit,mesflg
!!$ C-----------------------------------------------------------------------
!!$ C THE FOLLOWING FORTRAN-77 DECLARATION IS TO CAUSE THE VALUES OF THE
!!$ C LISTED (LOCAL) VARIABLES TO BE SAVED BETWEEN CALLS TO THIS ROUTINE.
!!$ C-----------------------------------------------------------------------
  save lunit,mesflg
  data lunit/-1/,mesflg/1/
!!$ C
!!$ C***FIRST EXECUTABLE STATEMENT IXSAV
  if (ipar .eq. 1) then
     if (lunit .eq. -1) lunit = iumach()
     ixsav = lunit
     if (iset) lunit = ivalue
  endif
!!$ C
  if (ipar .eq. 2) then
     ixsav = mesflg
     if (iset) mesflg = ivalue
  endif
!!$ C
  return
!!$ C----------------------- END OF FUNCTION IXSAV -------------------------
end function ixsav

!!$*DECK IUMACH

function iumach()
  integer::iumach
!!$ C***BEGIN PROLOGUE IUMACH
!!$ C***PURPOSE PROVIDE STANDARD OUTPUT UNIT NUMBER.
!!$ C***CATEGORY R1
!!$ C***TYPE   integer::(IUMACH-I)
!!$ C***KEYWORDS MACHINE CONSTANTS
!!$ C***AUTHOR HINDMARSH,ALAN C.,(LLNL)
!!$ C***DESCRIPTION
!!$ C *USAGE:
!!$ C    integer:: LOUT,IUMACH
!!$ C    LOUT = IUMACH()
!!$ C
!!$ C *FUNCTION RETURN VALUES:
!!$ C   LOUT : THE STANDARD LOGICAL UNIT FOR FORTRAN OUTPUT.
!!$ C
!!$ C***REFERENCES (NONE)
!!$ C***ROUTINES CALLED (NONE)
!!$ C***REVISION HISTORY (YYMMDD)
!!$ C  930915 DATE WRITTEN
!!$ C  930922 MADE USER-CALLABLE,AND OTHER COSMETIC CHANGES. (FNF)
!!$ C***END PROLOGUE IUMACH
!!$ C
!!$ C*INTERNAL NOTES:
!!$ C THE BUILT-IN VALUE OF 6 IS STANDARD ON A WIDE RANGE OF FORTRAN
!!$ C SYSTEMS. THIS MAY BE MACHINE-DEPENDENT.
!!$ C**END
!!$ C***FIRST EXECUTABLE STATEMENT IUMACH
  iumach = 6
!!$ C
  return
!!$ C----------------------- END OF FUNCTION IUMACH ------------------------
end function iumach
