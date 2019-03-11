* COPYRIGHT (c) 1988 AEA Technology and
* Council for the Central Laboratory of the Research Councils
C Original date 14 June 2001
C  June 2001: threadsafe version of MC41
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC71A(N,KASE,X,EST,W,IW,KEEP)
C
C      MC71A/AD ESTIMATES THE 1-NORM OF A SQUARE MATRIX A.
C      REVERSE COMMUNICATION IS USED FOR EVALUATING
C      MATRIX-VECTOR PRODUCTS.
C
C
C         N       INTEGER
C                 THE ORDER OF THE MATRIX.  N .GE. 1.
C
C         KASE    INTEGER
C                 SET INITIALLY TO ZERO . IF N .LE. 0 SET TO -1
C                 ON INTERMEDIATE RETURN
C                 = 1 OR 2.
C                ON FINAL RETURN
C                 =  0  ,IF SUCCESS
C                 = -1  ,IF N .LE.0
C
C         X       REAL ARRAY OF DIMENSION (N)
C                 IF 1-NORM IS REQUIRED
C                 MUST BE OVERWRITTEN BY
C
C                      A*X,             IF KASE=1,
C                      TRANSPOSE(A)*X,  IF KASE=2,
C
C                 AND MC71 MUST BE RE-CALLED, WITH ALL THE OTHER
C                 PARAMETERS UNCHANGED.
C                 IF INFINITY-NORM IS REQUIRED
C                 MUST BE OVERWRITTEN BY
C
C                      TRANSPOSE(A)*X,  IF KASE=1,
C                      A*X,             IF KASE=2,
C
C                 AND MC71 MUST BE RE-CALLED, WITH ALL THE OTHER
C                 PARAMETERS UNCHANGED.
C
C         EST     REAL
C                 CONTAINS AN ESTIMATE (A LOWER BOUND) FOR NORM(A).
C
C         W       REAL ARRAY OF DIMENSION (N)
C                 = A*V,   WHERE  EST = NORM(W)/NORM(V)
C                          (V  IS NOT RETURNED).
C         IW      INTEGER(N) USED AS WORKSPACE.
C
C         KEEP    INTEGER ARRAY LENGTH 5 USED TO PRESERVE PRIVATE
C                 DATA, JUMP, ITER, J AND JLAST BETWEEN CALLS,
C                 KEEP(5) IS SPARE.
C
C      REFERENCE
C      N.J. HIGHAM (1987) FORTRAN CODES FOR ESTIMATING
C      THE ONE-NORM OF A
C      REAL OR COMPLEX MATRIX, WITH APPLICATIONS
C      TO CONDITION  ESTIMATION, NUMERICAL ANALYSIS REPORT NO. 135,
C      UNIVERSITY OF MANCHESTER, MANCHESTER M13 9PL, ENGLAND.
C
C      SUBROUTINES AND FUNCTIONS
C
C
C      INTERNAL VARIABLES
C
C     .. Parameters ..
      INTEGER ITMAX
      PARAMETER (ITMAX=5)
      REAL ZERO,ONE
      PARAMETER (ZERO=0.0,ONE=1.0)
C     ..
C     .. Scalar Arguments ..
      REAL EST
      INTEGER KASE,N
C     ..
C     .. Array Arguments ..
      REAL W(*),X(*)
      INTEGER IW(*),KEEP(5)
C     ..
C     .. Local Scalars ..
      REAL ALTSGN,TEMP
      INTEGER I,ITER,J,JLAST,JUMP
C     ..
C     .. External Functions ..
      INTEGER ISAMAX
      EXTERNAL ISAMAX
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN,NINT,REAL
C     ..
C     .. Executable Statements ..
C
      IF (N.LE.0) THEN
        KASE = -1
        RETURN

      END IF

      IF (KASE.EQ.0) THEN
        DO 10 I = 1,N
          X(I) = ONE/REAL(N)
   10   CONTINUE
        KASE = 1
        JUMP = 1
        KEEP(1) = JUMP
        KEEP(2) = 0
        KEEP(3) = 0
        KEEP(4) = 0
        RETURN

      END IF
C
      JUMP  = KEEP(1)
      ITER  = KEEP(2)
      J     = KEEP(3)
      JLAST = KEEP(4)
C
      GO TO (100,200,300,400,500) JUMP
C
C      ................ ENTRY   (JUMP = 1)
C
  100 CONTINUE
      IF (N.EQ.1) THEN
        W(1) = X(1)
        EST = ABS(W(1))
C         ... QUIT
        GO TO 510

      END IF
C
      DO 110 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  110 CONTINUE
      KASE = 2
      JUMP = 2
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 2)
C
  200 CONTINUE
      J = ISAMAX(N,X,1)
      ITER = 2
C
C      MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
C
  220 CONTINUE
      DO 230 I = 1,N
        X(I) = ZERO
  230 CONTINUE
      X(J) = ONE
      KASE = 1
      JUMP = 3
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 3)
C
  300 CONTINUE
C
C      COPY X INTO W
C
      DO 310 I = 1,N
        W(I) = X(I)
  310 CONTINUE
      DO 320 I = 1,N
        IF (NINT(SIGN(ONE,X(I))).NE.IW(I)) GO TO 330
  320 CONTINUE
C
C      REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 410
C
  330 CONTINUE
      DO 340 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  340 CONTINUE
      KASE = 2
      JUMP = 4
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 4)
C
  400 CONTINUE
      JLAST = J
      J = ISAMAX(N,X,1)
      IF ((ABS(X(JLAST)).NE.ABS(X(J))) .AND. (ITER.LT.ITMAX)) THEN
        ITER = ITER + 1
        GO TO 220

      END IF
C
C      ITERATION COMPLETE.  FINAL STAGE.
C
  410 CONTINUE
      EST = ZERO
      DO 420 I = 1,N
        EST = EST + ABS(W(I))
  420 CONTINUE
C
      ALTSGN = ONE
      DO 430 I = 1,N
        X(I) = ALTSGN* (ONE+REAL(I-1)/REAL(N-1))
        ALTSGN = -ALTSGN
  430 CONTINUE
      KASE = 1
      JUMP = 5
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 5)
C
  500 CONTINUE
      TEMP = ZERO
      DO 520 I = 1,N
        TEMP = TEMP + ABS(X(I))
  520 CONTINUE
      TEMP = 2.0*TEMP/REAL(3*N)
      IF (TEMP.GT.EST) THEN
C
C      COPY X INTO W
C
        DO 530 I = 1,N
          W(I) = X(I)
  530   CONTINUE
        EST = TEMP
      END IF
C
  510 KASE = 0
C
 1010 CONTINUE
      KEEP(1) = JUMP
      KEEP(2) = ITER
      KEEP(3) = J
      KEEP(4) = JLAST
      RETURN
C
      END
* COPYRIGHT (c) 1988 AEA Technology
* Original date 17 Feb 2005

C 17th February 2005 Version 1.0.0. Replacement for FD05.

      REAL FUNCTION FD15A(T)
C----------------------------------------------------------------
C  Fortran 77 implementation of the Fortran 90 intrinsic
C    functions: EPSILON, TINY, HUGE and RADIX.  Note that
C    the RADIX result is returned as REAL.
C
C  The CHARACTER argument specifies the type of result:
C       
C   'E'  smallest positive real number: 1.0 + FD15A > 1.0, i.e.
C          EPSILON(REAL)
C   'T'  smallest full precision positive real number, i.e.
C          TINY(REAL)
C   'H'  largest finite positive real number, i.e.
C          HUGE(REAL)
C   'R'  the base of the floating point arithematic, i.e.
C          RADIX(REAL)
C
C    any other value gives a result of zero.
C----------------------------------------------------------------
      CHARACTER T

      IF ( T.EQ.'E' ) THEN
         FD15A = EPSILON(1.0)
      ELSE IF ( T.EQ.'T' ) THEN
         FD15A = TINY(1.0)
      ELSE IF ( T.EQ.'H' ) THEN
         FD15A = HUGE(1.0)
      ELSE IF ( T.EQ.'R' ) THEN
         FD15A = REAL(RADIX(1.0))
      ELSE
         FD15A = 0.0
      ENDIF
      RETURN
      END
! COPYRIGHT (c) 1995 Council for the Central Laboratory
!                    of the Research Councils
! Original date 23 March 2001
! Threadsafe version of IM01
!
! Version 1.3.0
! See ChangeLog for version history
!
      SUBROUTINE MI21I( ICNTL, CNTL, ISAVE, RSAVE )
      REAL CNTL( 5 )
      INTEGER ICNTL( 8 )
      INTEGER ISAVE(10)
      REAL RSAVE(6)
C
C  If A is symmetric, positive definite, MI21 solves the linear system
C
C      A x = b
C
C  using the
C
C  ===================
C  Conjugate Gradients
C  ===================
C
C  iterative method  optionally using preconditioning. If P
C  is a symmetric positive definite preconditioner and
C  P = P_L * P_L(transpose), the algorithm actually solves
C
C                             ^       ^                      ^
C        P_L A P_L(transpose) x = P_L b,  x = P_L(transpose) x
C
C  Each time a matrix-vector product with A or P is required,
C  control is passed back to the user.
C
C  MI21I/ID is the initialisation routine for MI21A/AD and should
C  be called once prior to calls to MI21A/AD.
C
C  Argument list.
C
C  ICNTL   (output) INTEGER control array, dimension 8.
C          ICNTL(1) is the stream number for error messages.
C             On exit, ICNTL(2) = 6.
C          ICNTL(2) is the stream number for warning messages.
C             On exit, ICNTL(2) = 6.
C          ICNTL(3) indicates whether the user wishes to use a
C             preconditioner. If ICNTL(3) is nonzero, preconditioning.
C             On exit, ICNTL(3) = 0
C          ICNTL(4) indicates whether
C             the convergence test offered by MI21A/AD is to be used.
C             If ICNTL(4) = 0, the computed solution x is accepted
C             if ||Ax - b|| < max( ||A x_0 - b || * CNTL(1), CNTL(2) ).
C             Otherwise, the user most perform his/her
C             own test for convergence when IACT = 1 is returned.
C             On exit, ICNTL(4) = 0
C          ICNTL(5) indicates whether the user wishes to supply an
C             initial guess for the solution vector x.
C             If ICNTL(5) = 0, the user does not wish to supply
C             an initial guess and x = (0,0,...,0) will be used
C             as the initial guess. Otherwise, the user
C             must supply an intial guess on the first call to
C             MI21A/AD. On exit, ICNTL(5) = 0.
C          ICNTL(6) determines the maximum number of iterations
C             allowed. It has default value -1 and, in this case,
C             the maximum number will be N. If the user does
C             not want the maximum number to be N, ICNTL(6) should
C             be set to the maximum  number the user wishes
C             to allow.
C          ICNTL(7) determines whether the normalized curvature is used
C             when testing for breakdown.  The
C             normalized curvature is only used if ICNTL(7)=1.
C             The default value is 0.
C          ICNTL(8) is spare and set to zero.
C
C  CNTL    (output)  REAL (DOUBLE PRECISION) control array, dimension 5.
C          CNTL(1) is a convergence tolerance.
C             On exit, set to the square root of machine precision.
C             If ICNTL(4) is nonzero, CNTL(1) not accessed by MI21A/AD.
C          CNTL(2) is a convergence tolerance. On exit, set to zero.
C             If ICNTL(4) is nonzero, CNTL(2) not accessed by MI21A/AD.
C          CNTL(3) is tolerance used to check whether the algorithm
C             has broken down. On exit, set to the machine precision
C          CNTL(4) and CNTL(5) are spare and set to zero.
C
C  ISAVE   (output) INTEGER ARRAY, length 10, used to hold the
C          routine's persistent integer data.  The array contents
C          must not be altered by the user.
C
C  RSAVE   (output) REAL ARRAY, length 6, used to hold the
C          routine's persistent real data.  The array contents
C          must not be altered by the user.
C
C  Local variables
C
      INTEGER I
C
C  Intrinsic Functions
C
      INTRINSIC SQRT
      ICNTL( 1 ) = 6
      ICNTL( 2 ) = 6
      ICNTL( 3 ) = 0
      ICNTL( 4 ) = 0
      ICNTL( 5 ) = 0
      ICNTL( 6 ) = - 1

      ICNTL(7) = 0
      ICNTL(8) = 0

      CNTL( 1 ) = SQRT( EPSILON(CNTL) )
      CNTL( 2 ) = 0.0
      CNTL( 3 ) = - 1.0

      CNTL(4) = 0.0
      CNTL(5) = 0.0

C  Initialize persistent data to avoid undefined assignment in MI21A
      DO 10 I = 1, 10
      ISAVE(I) = 0
   10 CONTINUE
      DO 20 I = 1, 6
      RSAVE(I) = 0.0
   20 CONTINUE
      RETURN
C
C  End of MI21I/MI21ID
C
      END

      SUBROUTINE MI21A( IACT, N, W, LDW, LOCY, LOCZ, RESID,
     *                   ICNTL, CNTL, INFO, ISAVE, RSAVE )
      REAL RESID
      INTEGER          IACT, LDW, LOCY, LOCZ, N
      REAL CNTL( 5 ), W( LDW, 4 )
      INTEGER          ICNTL( 8 ), INFO( 4 )
      INTEGER ISAVE(10)
      REAL RSAVE(6)
C
C  Argument list.
C
C  IACT    (input) INTEGER.
C          IACT must be set to 0 prior to first call.
C          On each exit, IACT indicates the action required by
C          the user. Possible values of IACT and the action
C          required are as follows:
C  -1      Fatal error (see INFO(1)). Terminate computation.
C   1      If ICNTL(4) = 0 (the default), convergence has been
C          achieved and the user should terminate the computation.
C          If ICNTL(4) is nonzero, the user should test the norm of
C          the residual in RESID for convergence and recall MI21A/AD
C          if convergence has not been achieved.
C   2      The user must perform the matrix-vector product
C          y := Az,
C          and recall MI21A/AD. The vectors  y and z are held in the
C          first N entries of columns LOCZ and LOCZ of array W,
C          respectively.
C   3      The user must perform the preconditioning operation
C          y := Pz, where P is the preconditioner,
C          and recall MI21A/AD. The vectors y and z
C          are  held in the first N entries of
C          columns LOCZ and LOCZ of array W, respectively.
C
C  N       (input) INTEGER.
C          On entry, the dimension of the matrix.
C          Unchanged on exit.
C
C  W       (right-hand side, solution and workspace, input/output)
C          REAL (DOUBLE PRECISION) array, dimension (LDW,4).
C          Prior to the first call, the first N entries of column
C          1 must be set to hold the right-hand side vector b and,
C          if ICNTL(5) is nonzero, the first N entries of column 2
C          must be set to the initial estimate of the solution vector
C          x.  On exit with IACT = 1, the first N entries of column 1
C          hold the current residual vector r = b - Ax, and the
C          current estimate of the solution x is held in the
C          first N entries of column 2.  On exit
C          with IACT > 1, the user is required to perform a
C          computation with columns LOCZ and LOCZ of W
C          (see argument IACT).  The remaining columns of
C          W must not be altered by the user between calls to
C          MI21A/AD.
C
C  LDW     (input) INTEGER
C          The leading dimension of the array W. LDW >= max(1,N).
C
C  LOCY, LOCZ (output) INTEGER
C          On exit with IACT > 1, LOCY and LOCZ define the columns
C          of the array W which hold y and z.
C          (see IACT).
C
C  RESID   (output) REAL (DOUBLE PRECISION)
C          On exit with IACT = 1,
C          RESID holds ||b - Ax||, where x is the
C          iterated solution.
C          If ICNTL(4) is nonzero, on exit with IACT = 1,
C          the user should carry out his/her
C          own test for convergnce at this point.
C
C  ICNTL   (input) INTEGER control array of dimension 8.
C          ICNTL may be initalised by calling MI21I/ID.
C          See MI21I/ID for details.
C
C  CNTL    (input) REAL (DOUBLE PRECISION) control array of dimension 5.
C          CNTL may be initalised by calling MI21I/ID.
C          See MI21I/ID for details.
C
C  INFO    (output) INTEGER ARRAY , length 4.
C          If INFO(1) = 0 on exit, no errors or warnings issued.
C          If INFO(1) > 0 on exit, warning to user.
C          INFO(1) = 1, value of CNTL(1) is out-of-range (u,1.0)
C          The default sqrt(u) is used, u=machine precision.
C          If INFO(1) < 0 on exit, illegal input parameter,
C               or breakdown occured during iteration.
C
C                Error input data:
C
C                   -1: matrix dimension N < 1
C                   -2: LDW < N
C
C                BREAKDOWN: If curvature smaller than CNTL(3) is
C                   encountered, the program will terminate.
C
C                   -3: CURV < CNTL(3):
C          On each exit, INFO(2) holds the number of
C          iterations performed so far.
C
C  ISAVE   (input/output) INTEGER ARRAY, length 10, used to hold the
C          routine's persistent integer data.  The array contents
C          must not be altered by the user.
C
C  RSAVE   (input/output) REAL ARRAY, length 6, used to hold the
C          routine's persistent real data.  The array contents
C          must not be altered by the user.
C
C
C  BLAS calls:    SAXPY, SCOPY, SDOT, SNRM2, SSCAL
C
      REAL ONE, ZERO
      PARAMETER ( ONE = 1.0E+0, ZERO = 0.0E+0 )
C
C  Local variables
C
      INTEGER          I, IPOS, ITMAX, B, P, Q, R, X, Z
      REAL ALPHA, BETA, BNRM2, CURV, RHO, RHO1
C
C  External Functions
C
      REAL SDOT, SNRM2
      EXTERNAL         SDOT, SNRM2
C
C  Intrinsic Functions
C
      INTRINSIC        SQRT
C
C  External Subroutines
C
      EXTERNAL         SAXPY, SCOPY, SSCAL
C
C  Restore persistent data
C
      IPOS  = ISAVE(1)
      ITMAX = ISAVE(2)
      B     = ISAVE(3)
      P     = ISAVE(4)
      Q     = ISAVE(5)
      R     = ISAVE(6)
      X     = ISAVE(7)
      Z     = ISAVE(8)

      BNRM2 = RSAVE(1)
      RHO   = RSAVE(2)
      RHO1  = RSAVE(3)
      CURV  = RSAVE(4)
C
C  Jump to appropriate place in code
C
      IF ( IACT .NE. 0 ) THEN
C
C  Immediate return if error on a previous call
C
         IF ( IACT .LT. 0 ) GO TO 1000
C
C  Immediate return if convergence already achieved
C
         IF ( IACT .EQ. 1 .AND. ICNTL( 4 ) .EQ. 0 ) GO TO 1000
         IF ( IACT .EQ. 1 .AND. BNRM2 .EQ. ZERO ) GO TO 1000
C
C  Branch
C
         GO TO ( 30, 50, 60, 70 ), IPOS
      END IF
C
C  Initial call
C
      INFO( 1 ) = 0
C
C  Test the input parameters
C
      IF ( N .LT. 1 ) THEN
         INFO( 1 ) = - 1
      ELSE IF ( LDW .LT. N ) THEN
         INFO( 1 ) = - 2
      END IF
      IF ( INFO( 1 ) .LT. 0 ) THEN
         IACT = - 1
         IF ( ICNTL( 1 ) .GT. 0 ) WRITE( ICNTL( 1 ), 2000 ) INFO( 1 )
         GO TO 1000
      END IF
C
C  Alias workspace columns
C
      B = 1
      R = B
      X = 2
      P = 3
      Q = 4
      IF ( ICNTL( 3 ) .NE. 0 ) THEN
         Z = Q
      ELSE
         Z = R
      END IF
C
C  Set INFO(2) and ITMAX
C
      INFO( 2 ) = 0
      IF ( ICNTL( 6 ) .GT. 0 ) THEN
         ITMAX = ICNTL( 6 )
      ELSE
         ITMAX = N
      END IF
C
C  Set CNTL(3)
C
      IF ( CNTL( 3 ) .LE. ZERO ) CNTL( 3 ) = N * EPSILON(CNTL)
C
C  Compute ||b||
C
      BNRM2 = SNRM2( N, W( 1, B ), 1 )
C
C  Immediate return if ||b|| = 0
C
      IF ( BNRM2. EQ. ZERO ) THEN
         IACT = 1
         DO 10 I = 1, N
            W( I, X ) = ZERO
            W( I, B ) = ZERO
   10    CONTINUE
         RESID = ZERO
         GO TO 1000
      END IF
C
C  Check value of CNTL(1)
C
      IF ( ICNTL( 4 ) .EQ. 0 ) THEN
         IF ( CNTL( 1 ).LT.EPSILON(CNTL) .OR. CNTL( 1 ).GT.ONE ) THEN
            INFO( 1 ) = 1
            IF (ICNTL( 2 ) .GT. 0 ) THEN
               WRITE( ICNTL( 2 ), 2010 ) INFO( 1 )
               WRITE( ICNTL( 2 ), 2020 )
            END IF
            CNTL( 1 ) = SQRT( EPSILON(CNTL) )
         END IF
      END IF
C
C  Compute initial residual
C
C  If the user has not supplied an initial guess, set X = 0
C  as the initial guess
C
      IF ( ICNTL( 5 ) .EQ. 0 ) THEN
         DO 20 I = 1, N
            W( I, X ) = ZERO
   20    CONTINUE
         GO TO 40
      ELSE
C
C  Initial guess supplied by user
C  If initial guess for solution is x = 0 no action is required (r = b)
C
         IF ( SNRM2( N, W( 1, X ), 1 ) .EQ. ZERO ) GO TO 40
C
C  Otherwise, return to user to compute Ax.
C  Column P can be used temporarily to hold Ax
C
         IPOS = 1
         IACT = 2
         LOCY = P
         LOCZ = X
         GO TO 1000
      END IF
C
C  Compute r = b - Ax
C
   30 CONTINUE
      CALL SAXPY( N, - ONE, W( 1, P ), 1, W( 1, R ), 1 )
C
C  Compute ||r||
C
      BNRM2 = SNRM2( N, W( 1, R ), 1 )
C
C  Main iteration loop
C
   40 CONTINUE
C
C  Update iteration count
C
      INFO( 2 ) = INFO( 2 ) + 1
C
C  Check maximum number of iterations has not been exceeded
C
      IF ( INFO( 2 ) .GT. ITMAX ) THEN
         INFO( 1 ) = - 4
         IACT = - 1
         IF (ICNTL( 1 ) .GT. 0 ) THEN
            WRITE( ICNTL( 1 ), 2000 ) INFO( 1 )
            WRITE( ICNTL( 1 ), 2030 ) ITMAX
         END IF
         GO TO 1000
      END IF
C
C  Return to user to obtain preconditioner Z = P^-1 R
C
      IF ( ICNTL( 3 ) .NE. 0 ) THEN
         IPOS = 2
         IACT = 3
         LOCY = Z
         LOCZ = R
         GO TO 1000
      END IF
   50 CONTINUE
C
C  Compute the inner product R^T Z
C
      RHO = SDOT( N, W( 1, R ), 1, W( 1, Z ), 1 )
C
C  Compute search direction vector P
C
      IF ( INFO( 2 ) .EQ. 1 ) THEN
C
C  Special case: first iteration
C
         CALL SCOPY( N, W( 1, Z ), 1, W( 1, P ), 1 )
      ELSE
         BETA = RHO / RHO1
C
C  later iterations
C
         CALL SSCAL( N, BETA, W( 1, P ), 1 )
         CALL SAXPY( N, ONE, W( 1, Z ), 1, W( 1, P ), 1 )
      END IF
C
C  Return to user for matrix-vector product Q = A P
C
      IPOS = 3
      IACT = 2
      LOCY = Q
      LOCZ = P
      GO TO 1000
   60 CONTINUE
C
C  Obtain the curvature along P
C
      CURV = SDOT( N, W( 1, P ), 1, W( 1, Q ), 1 )
C
C  If the curvature is negative, A is indefinite
C
      IF (ICNTL(7).EQ.1) THEN
        IF ( CURV .LT. CNTL( 3 )*SDOT( N, W( 1, P ), 1, W( 1, P ), 1 ))
     +     INFO( 1 ) = - 3
      ELSE IF ( CURV .LT. CNTL( 3 ) ) THEN
         INFO( 1 ) = - 3
      END IF
      IF (INFO(1).EQ.-3) THEN
         IACT = - 1
         IF ( ICNTL( 1 ) .GT. 0 ) WRITE( ICNTL( 1 ), 2000 ) INFO( 1 )
         GO TO 1000
      END IF
C
C  Compute the stepsize
C
      ALPHA = RHO / CURV
C
C  Update the estimate of the solution
C
      CALL SAXPY( N, ALPHA, W( 1, P ), 1, W( 1, X ), 1 )
C
C  Update the residual
C
      CALL SAXPY( N, - ALPHA, W( 1, Q ), 1, W( 1, R ), 1 )
C
C  Iteration complete. Check convergence
C
      RESID = SNRM2( N, W( 1, R ), 1 )
      IPOS = 4
      IF ( ICNTL( 4 ) .NE. 0 ) THEN
C
C  Return the residual to the user for convergence testing
C
         IACT = 1
         GO TO 1000
      ELSE
C
C  Test the scaled residual for convergence
C
         IF ( RESID .LE. MAX( BNRM2 * CNTL( 1 ), CNTL( 2 ) ) ) THEN
C
C  Convergence achieved
C
            IACT = 1
            GO TO 1000
         END IF
      END IF
   70 CONTINUE
      RHO1 = RHO
C
C  Next iteration
C
      GO TO 40
C
C  Save persistent data and return to caller
C
 1000 CONTINUE
      ISAVE(1) = IPOS
      ISAVE(2) = ITMAX
      ISAVE(3) = B
      ISAVE(4) = P
      ISAVE(5) = Q
      ISAVE(6) = R
      ISAVE(7) = X
      ISAVE(8) = Z
      RSAVE(1) = BNRM2
      RSAVE(2) = RHO
      RSAVE(3) = RHO1
      RSAVE(4) = CURV
      RETURN
C
C  Non-executable statements
C
 2000 FORMAT( / ' Error message from MI21A/AD. INFO(1) = ', I4 )
 2010 FORMAT( / ' Warning message from MI21A/AD. INFO(1) = ', I4 )
 2020 FORMAT( ' Convergence tolerance out of range.' )
 2030 FORMAT( ' Number of iterations required exceeds the maximum of ',
     *        I8, / ' allowed by ICNTL(6)' )
C
C  End of MI21A/AD
C
      END
! COPYRIGHT (c) 1995 Council for the Central Laboratory
!                    of the Research Councils
!
! Version 1.2.0
! See ChangeLog for version history.
!
      SUBROUTINE MI24I(ICNTL,CNTL,ISAVE,RSAVE,LSAVE)
      INTEGER          ICNTL( 8 )
      REAL CNTL( 4 )
      INTEGER ISAVE(17)
      REAL RSAVE(9)
      LOGICAL LSAVE(4)
      INTEGER I
C
C  MI24 solves the linear system Ax = b using the
C  Generalized Minimal Residual with restarts (GMRES) iterative method
C  optionally using preconditioning
C                    ^                   ^
C           P(1)AP(2)x = P(1)b,  x = P(2)x

C  P(1), P(2) are the preconditioners, which are not passed to the code,
C  but each time a matrix-vector product with P(1) or P(2) is required,
C  control is passed back to the user.
C
C  Similarly, the matrix A is not passed to the code, but when
C  a matrix-vector with A is required, control is
C  passed back to the user.
C
C  MI24I/ID is the initialisation routine for MI24A/AD and should
C  be called once prior to calls to MI24A/AD.
C
C  Argument list.
C
C  ICNTL   (output) INTEGER control array, dimension 8.
C          ICNTL(1) is the stream number for error messages.
C             On exit, ICNTL(2) = 6.
C          ICNTL(2) is the stream number for warning messages.
C             On exit, ICNTL(2) = 6.
C          ICNTL(3) indicates whether the user wishes to use a
C             preconditioner. If ICNTL(3) is 1 preconditioning will
C             be performed on the left, if it is 2 preconditioning will
C             be performed on the right, and if it is zero, no
C             preconditioning is required.
C             On exit, ICNTL(3) = 0
C          ICNTL(4) indicates whether
C             the convergence test offered by MI24A/AD is to be used.
C             If ICNTL(4) = 0, the computed solution x is accepted
C             if ||Ax - b||/||b|| < CNTL(1).
C             Otherwise, the user most perform his/her
C             own test for convergence when IACT = 1 is returned.
C             On exit, ICNTL(4) = 0
C          ICNTL(5) indicates whether the user wishes to supply an
C             initial guess for the solution vector x.
C             If ICNTL(5) = 0, the user does not wish to supply
C             an initial guess and x = (0,0,...,0) will be used
C             as the initial guess. Otherwise, the user
C             must supply an intial guess on the first call to
C             MI24A/AD. On exit, ICNTL(5) = 0.
C          ICNTL(6) determines the maximum number of iterations
C             allowed. It has default value -1 and, in this case,
C             the maximum number will be N. If the user does
C             not want the maximum number to be N, ICNTL(6) should
C             be set to the maximum  number the user wishes
C             to allow.
C          ICNTRL(7) and ICNTRL(8) are spare and set to zero.
C  CNTL    (output)  REAL (DOUBLE PRECISION) control array, dimension 4.
C          CNTL(1) is a convergence tolerance.
C             On exit, set to square root of machine precision.
C             If ICNTL(4) is nonzero, CNTL(1) not accessed by MI24A/AD.
C          CNTL(2) is a convergence tolerance.
C             On exit, set to zero.
C             If ICNTL(4) is nonzero, CNTL(2) not accessed by MI24A/AD.
C          CNTRL(3) and CNTRL(4) are spare and set to zero.
C
C  ISAVE   (output) INTEGER ARRAY, length 17, used to hold the
C          routine's persistent integer data.  The array contents
C          must not be altered by the user.
C
C  RSAVE   (output) REAL ARRAY, length 9, used to hold the
C          routine's persistent real data.  The array contents
C          must not be altered by the user.
C
C  LSAVE   (output) LOGICAL ARRAY, length 4, used to hold the
C          routine's persistent logical data.  The array contents
C          must not be altered by the user.
C
C     Intrinsic Functions
C
      INTRINSIC SQRT
      ICNTL( 1 ) = 6
      ICNTL( 2 ) = 6
      ICNTL( 3 ) = 0
      ICNTL( 4 ) = 0
      ICNTL( 5 ) = 0
      ICNTL( 6 ) = - 1

      ICNTL(7) = 0
      ICNTL(8) = 0

      CNTL( 1 ) = SQRT( EPSILON(CNTL) )
      CNTL( 2 ) = 0.0

      CNTL(3) = 0.0
      CNTL(4) = 0.0

C  Initialize persistent data to avoid undefined assignment in MI24A
      DO 10 I = 1, 17
      ISAVE(I) = 0
   10 CONTINUE
      DO 20 I = 1, 9
      RSAVE(I) = 0.0
   20 CONTINUE
      DO 30 I = 1, 4
      LSAVE(I) = .FALSE.
   30 CONTINUE

      RETURN
C
C  End of MI24I/MI24ID
C
      END

      SUBROUTINE MI24A( IACT, N, M, W, LDW, LOCY, LOCZ, H, LDH, RESID,
     *                   ICNTL, CNTL, INFO, ISAVE, RSAVE, LSAVE )
      REAL RESID
      INTEGER          IACT, N, M, LDW, LOCY, LOCZ, LDH
      REAL CNTL( 4 ), W( LDW, M + 7 ), H( LDH, M + 2 )
      INTEGER          ICNTL( 8 ), INFO( 4 )
      INTEGER ISAVE(17)
      REAL RSAVE(9)
      LOGICAL LSAVE(4)
C
C  Argument list.
C
C  IACT    (input) INTEGER.
C          IACT must be set to 0 prior to first call.
C          On each exit, IACT indicates the action required by
C          the user. Possible values of IACT and the action
C          required are as follows:
C  -1      Fatal error (see INFO(1)). Terminate computation.
C   1      If ICNTL(4) = 0 (the default), convergence has been
C          achieved and the user should terminate the computation.
C          If ICNTL(4) is nonzero, the user should test the norm of
C          the residual in RESID for convergence and recall MI24A/AD
C          if convergence has not been achieved.
C   2      The user must perform the matrix-vector product
C          y := Az,
C          and recall MI24A/AD. The vectors  y and z are held in the
C          first N entries of columns LOCZ and LOCZ of array W,
C          respectively.
C   3      The user must perform the left preconditioning operation
C          y := P_L z,  where P_L is the left preconditioner
C          and recall MI24A/AD. The vectors y and z
C          are  held in the first N entries of
C          columns LOCZ and LOCZ of array W, respectively.
C   4      The user must perform the right preconditioning operation
C          y := P_R z,  where P_R is the right preconditioner
C          and recall MI24A/AD. The vectors y and z
C          are  held in the first N entries of
C          columns LOCZ and LOCZ of array W, respectively.
C
C  N       (input) INTEGER.
C          On entry, the dimension of the matrix.
C          Unchanged on exit.
C
C  W       (right-hand side, solution and workspace, input/output)
C          REAL (DOUBLE PRECISION) array, dimension (LDW,M+7).
C          Prior to the first call, the first N entries of column
C          1 must be set to hold the right-hand side vector b and,
C          if ICNTL(5) is nonzero, the first N entries of column 2
C          must be set to the initial estimate of the solution vector
C          x.  On exit with IACT = 1, the first N entries of column 1
C          hold the current residual vector r = b - Ax, and the
C          current estimate of the solution x is held in the
C          first N entries of column 2.  On exit
C          with IACT > 1, the user is required to perform a
C          computation with columns LOCZ and LOCZ of W
C          (see argument IACT).  The remaining columns of
C          W must not be altered by the user between calls to
C          MI24A/AD.
C
C  LDW     (input) INTEGER
C          The leading dimension of the array W. LDW >= max(1,N).
C
C  LOCY, LOCZ (output) INTEGER
C          On exit with IACT > 1, LOCY and LOCZ define the columns
C          of the array W which hold y and z.
C          (see IACT).
C
C  H       (workspace, output)
C          DOUBLE PRECISION array, dimension (LDH,M+2).
C          This workspace is used for constructing and storing the
C          upper Hessenberg matrix. The two extra columns are used to
C          store the Givens rotation matrices.
C
C  LDH    (input) INTEGER
C          The leading dimension of the array H. LDH >= M+1.
C
C  RESID   (output) REAL (DOUBLE PRECISION)
C          On exit with IACT = 1,
C          RESID holds ||b - Ax||, where x is the
C          iterated solution.
C          If ICNTL(4) is nonzero, on exit with IACT = 1,
C          the user should carry out his/her
C          own test for convergnce at this point.
C
C  ICNTL   (input) INTEGER control array of dimension 8.
C          ICNTL may be initalised by calling MI24I/ID.
C          See MI24I/ID for details.
C
C  CNTL    (input) REAL (DOUBLE PRECISION) control array of dimension 4.
C          CNTL may be initalised by calling MI24I/ID.
C          See MI24I/ID for details.
C
C  INFO    (output) INTEGER ARRAY , length 4.
C          If INFO(1) = 0 on exit, no errors or warnings issued.
C          If INFO(1) > 0 on exit, warning to user.
C
C             Warning message:
C
C             1: value of CNTL(1) is out-of-range (u,1.0)
C              The default sqrt(u) is used, u=machine precision.
C
C          If INFO(1) < 0 on exit, an error has occurred.
C
C                Error input data:
C
C                   -1: matrix dimension N < 1
C                   -2: restart interval M < 1
C                   -3: LDW < N
C                   -4: LDH < M + 1
C                   -5: Too many iterations have been performed
C
C          On each exit, INFO(2) holds the number of
C          iterations performed so far.
C
C  ISAVE   (output) INTEGER ARRAY, length 17, used to hold the
C          routine's persistent integer data.  The array contents
C          must not be altered by the user.
C
C  RSAVE   (output) REAL ARRAY, length 9, used to hold the
C          routine's persistent real data.  The array contents
C          must not be altered by the user.
C
C  LSAVE   (output) LOGICAL ARRAY, length 4, used to hold the
C          routine's persistent logical data.  The array contents
C          must not be altered by the user.
C
C  BLAS CALLS:   SAXPY, SCOPY, SDOT, SNRM2, SROT, SROTG, SSCAL,
C                STRSV, SGEMV
C
      REAL    ZERO, ONE, POINT1
      PARAMETER         ( ZERO = 0.0E+0, ONE = 1.0E+0, POINT1 = 1.0E-1 )
      INTEGER             I, K, ITMAX, CS, SN, R, S, V, UU, Y, RES
      INTEGER             B, X, IPOS, U
      LOGICAL             LEFT, RIGHT
      REAL    AA, BB, BNRM2 , RNORM , SDOT, SNRM2, RSTOP
      REAL    PRESID, PRSTOP
      EXTERNAL            SAXPY, SCOPY, SDOT, SNRM2, SROT, SROTG, SSCAL
      EXTERNAL            STRSV, SGEMV
C
C  Restore persistent data
C
      IPOS   = ISAVE(1)
      ITMAX  = ISAVE(2)
      B      = ISAVE(3)
      I      = ISAVE(4)
      K      = ISAVE(5)
      R      = ISAVE(6)
      X      = ISAVE(7)
      U      = ISAVE(8)
      V      = ISAVE(9)
      S      = ISAVE(10)
      Y      = ISAVE(11)
      CS     = ISAVE(12)
      SN     = ISAVE(13)
      UU     = ISAVE(14)
      RES    = ISAVE(15)

      BNRM2  = RSAVE(1)
      AA     = RSAVE(2)
      BB     = RSAVE(3)
      RNORM  = RSAVE(4)
      PRESID = RSAVE(5)
      RSTOP  = RSAVE(6)
      PRSTOP = RSAVE(7)

      LEFT   = LSAVE(1)
      RIGHT  = LSAVE(2)
C
C  Jump to appropriate place in code
C
      IF ( IACT .NE. 0 ) THEN
C
C  Immediate return if error on a previous call
C
         IF ( IACT .LT. 0 ) GO TO 1000
C
C  Immediate return if convergence already achieved
C
         IF ( IACT .EQ. 1 .AND. ICNTL( 4 ) .EQ. 0 ) GO TO 1000
         IF ( IACT .EQ. 1 .AND. BNRM2 .EQ. ZERO ) GO TO 1000
C
C  Branch
C
         GO TO ( 40, 60, 70, 100, 110, 120, 160 ), IPOS
      END IF
C
C  Initial call
C
      INFO( 1 ) = 0
C
C     Test the input parameters.
C
      IF ( N .LT. 1 ) THEN
         INFO( 1 ) = - 1
      ELSE IF ( M .LT. 1 ) THEN
         INFO( 1 ) = - 2
      ELSE IF ( LDW .LT. MAX( 1, N ) ) THEN
         INFO( 1 ) = - 3
      ELSE IF ( LDH .LT. M + 1 ) THEN
         INFO( 1 ) = - 4
      ENDIF
      IF ( INFO( 1 ) .LT. 0 ) THEN
         IACT = - 1
         IF ( ICNTL( 1 ) .GT. 0 ) WRITE( ICNTL( 1 ), 2000 ) INFO( 1 )
         GO TO 1000
      END IF
C
C  Set INFO(2) and ITMAX
C
      INFO( 2 ) = 0
      IF ( ICNTL( 6 ) .GT. 0 ) THEN
         ITMAX = ICNTL( 6 )
      ELSE
         ITMAX = 2 * N
      END IF
C
C  Alias workspace columns
C
      RES = 1
      X   = 2
      S   = 3
      B   = 4
      UU  = 5
      Y   = 6
      V   = 7
C
C  Store the Givens parameters in matrix H
C
      CS = M + 1
      SN = CS + 1
C
C  Compute ||b||
C
      BNRM2 = SNRM2( N, W( 1, RES ), 1 )
C
C  Immediate return if ||b|| = 0
C
      IF ( BNRM2. EQ. ZERO ) THEN
         IACT = 1
         DO 10 I = 1, N
            W( I, X ) = ZERO
            W( I, RES ) = ZERO
   10    CONTINUE
         RESID = ZERO
         GO TO 1000
      END IF
C
C  Check value of CNTL(1)
C
      IF ( ICNTL( 4 ) .EQ. 0 ) THEN
         IF ( CNTL( 1 ).LT.EPSILON(CNTL) .OR. CNTL( 1 ).GT.ONE ) THEN
            INFO( 1 ) = 1
            IF (ICNTL( 2 ) .GT. 0 ) THEN
               WRITE( ICNTL( 2 ), 2010 ) INFO( 1 )
               WRITE( ICNTL( 2 ), 2020 )
            END IF
            CNTL( 1 ) = SQRT( EPSILON(CNTL) )
         END IF
      END IF
C
C  LEFT indicates that a left preconditioner will be used
C  RIGHT indicates that a right preconditioner will be used
C
      LEFT  = ICNTL( 3 ) .EQ. 1 .OR. ICNTL( 3 ) .EQ. 3
      RIGHT = ICNTL( 3 ) .EQ. 2 .OR. ICNTL( 3 ) .EQ. 3
C
C  If no initial guess is provided, set x to zero
C
      IF ( ICNTL( 5 ) .EQ. 0 ) THEN
         DO 20 I = 1, N
            W( I, X ) = ZERO
   20    CONTINUE
      END IF
C
C  Start computing the residual
C
      CALL SCOPY( N, W( 1, RES ), 1, W( 1, B ), 1 )
      IF ( SNRM2( N, W( 1, X ), 1 ) .EQ. ZERO ) GO TO 50
C
C  Start of the main iteration
C
   30 CONTINUE
C
C  Return to user for matrix-vector product Y = A * X
C
         IPOS = 1
         IACT = 2
         LOCY = Y
         LOCZ = X
         GO TO 1000
   40    CONTINUE
C
C  Finalize the residual
C
         CALL SAXPY( N, - ONE, W( 1, Y ), 1, W( 1, RES ), 1 )
   50    CONTINUE
C
C  Compute the norm of the residual
C
         RESID = SNRM2( N, W( 1, RES ), 1 )
C
C  Assign the stopping tolerence on the first iteration
C
         IF ( INFO( 2 ) .EQ. 0 ) THEN
            RSTOP  = MAX( RESID * CNTL( 1 ), CNTL( 2 ) )
            PRSTOP = RSTOP
         END IF
C
C  If too many iterations have occured, exit
C
         IF ( INFO( 1 ) .LT. 0 ) THEN
            IACT = - 1
            IF ( ICNTL( 1 ) .GT. 0 ) WRITE( ICNTL( 1 ), 2000 ) INFO( 1 )
            GO TO 1000
         END IF
C
C  Return to user to obtain preconditioner R = P_L^-1 RES
C
         IF ( LEFT ) THEN
            R = UU
            IPOS = 2
            IACT = 3
            LOCY = R
            LOCZ = RES
            GO TO 1000
         ELSE
            R = RES
         END IF
   60    CONTINUE
C
C  Check for convergence
C
         IF ( ICNTL( 4 ) .NE. 0 .OR. ( ICNTL( 4 ) .EQ. 0 .AND.
     *        RESID .LE. RSTOP ) ) THEN
            IACT = 1
            IPOS = 3
            GO TO 1000
         END IF
   70    CONTINUE
C
C  Construct the first column of V
C
         CALL SCOPY( N, W( 1, R ), 1, W( 1, V ), 1 )
         RNORM = SNRM2( N, W( 1, V ), 1 )
         CALL SSCAL( N, ONE / RNORM, W( 1, V ), 1 )
C
C  Initialize S to the elementary vector E1 scaled by RNORM
C
         W( 1, S ) = RNORM
         DO 80 K = 2, N
            W( K, S ) = ZERO
   80    CONTINUE
C
C  Start of the inner iteration
C
         I = 0
   90    CONTINUE
            I = I + 1
C
C  Update iteration count
C
            INFO( 2 ) = INFO( 2 ) + 1
C
C  Check maximum number of iterations has not been exceeded
C
            IF ( INFO( 2 ) .GT. ITMAX ) THEN
               I = I - 1
               INFO( 1 ) = - 5
               IF ( ICNTL( 1 ) .GT. 0 ) WRITE( ICNTL( 1 ), 2030 ) ITMAX
               IF ( I .NE. 0 ) GO TO 150
               IACT = - 1
               IF ( ICNTL( 1 ) .GT. 0 )
     *            WRITE( ICNTL( 1 ), 2000 ) INFO( 1 ) 
               GO TO 1000
            END IF
C
C  Return to user to obtain preconditioner Y = P_R^-1 V
C
            IF ( RIGHT ) THEN
               IPOS = 4
               IACT = 4
               LOCY = Y
               LOCZ = V + I - 1
               GO TO 1000
            END IF
  100       CONTINUE
C
C  Return to user for matrix-vector product
C
            IPOS = 5
            IACT = 2
            IF ( RIGHT ) THEN
               LOCY = RES
               LOCZ = Y
            ELSE
               LOCY = RES
               LOCZ = V + I - 1
            END IF
            GO TO 1000
  110       CONTINUE
C
C  Return to user to obtain preconditioner W = P_L^-1 RES
C
            IF ( LEFT ) THEN
               U = UU
               IPOS = 6
               IACT = 3
               LOCY = UU
               LOCZ = RES
               GO TO 1000
            ELSE
               U = RES
            END IF
  120       CONTINUE
C
C  Construct I-th column of H orthonormal to the previous I-1 columns
C  using the Gram-Schmidt process on V and U
C
C           CALL MI24B( I, N, H( 1, I ), W( 1, V ), LDW, W( 1, U ) )
            DO 130 K = 1, I
               H( K, I ) = SDOT( N, W( 1, U ), 1,
     *                              W( 1, V + K - 1 ), 1 )
               CALL SAXPY( N, - H( K, I ), W( 1, V + K - 1 ), 1,
     *                                     W( 1, U ), 1 )
  130       CONTINUE
            H( I + 1, I ) = SNRM2( N, W( 1, U ), 1 )
            CALL SCOPY( N, W( 1, U ), 1, W( 1, V + I ), 1 )
            CALL SSCAL( N, ONE / H( I + 1, I ), W( 1, V + I ), 1 )
C
C  Apply Givens rotations to the I-th column of H. This "updating" of
C  the QR factorization effectively reduces the Hessenberg matrix to
C  upper triangular form during the M iterations.
C
            DO 140 K = 1, I - 1
               CALL SROT( 1, H( K, I ), LDH, H( K + 1, I ), LDH,
     *                    H( K, CS ), H( K, SN ) )
  140       CONTINUE
C
C  Construct the I-th rotation matrix, and apply it to H so that
C  H(I+1,I) = 0
C
            AA = H( I, I )
            BB = H( I + 1, I )
            CALL SROTG( AA, BB, H( I, CS ), H( I, SN ) )
            CALL SROT( 1, H( I, I ), LDH, H( I + 1, I ), LDH,
     *            H( I, CS ), H( I, SN ) )
C
C  Apply the I-th rotation matrix to [ S(I), S(I+1) ]'. This
C  gives an approximation of the residual norm. If less than
C  tolerance, update the approximation vector X and quit
C
            IF ( I .LT. N ) THEN
               CALL SROT( 1, W( I, S ), LDW, W( I + 1, S ),
     *                    LDW, H( I, CS ), H( I, SN ) )
               PRESID = ABS( W( I + 1, S ) )
               IF ( PRESID .LE. PRSTOP .AND. ICNTL( 4 ) .EQ. 0 ) THEN
                  PRSTOP = PRSTOP * POINT1
                  GO TO 150
               END IF
               IF ( I .LT. M )  GO TO 90
            END IF
C
C  Compute current solution vector X
C
  150    CONTINUE
         CALL SCOPY( I, W( 1, S ), 1, W( 1, Y ), 1 )
         CALL STRSV( 'UPPER', 'NOTRANS', 'NONUNIT', I, H, LDH,
     *               W( 1, Y ), 1 )
C
C  Compute current update vector UU = V*Y
C
         CALL SGEMV( 'NOTRANS', N, I, ONE, W( 1, V ), LDW, W( 1, Y ), 1,
     *               ZERO, W( 1, UU ), 1 )
C
C  Return to user to obtain preconditioner Y = P_R^-1 UU
C
         IF ( RIGHT ) THEN
            IPOS = 7
            IACT = 4
            LOCY = Y
            LOCZ = UU
            GO TO 1000
         END IF
  160    CONTINUE
C
C  Update X
C
         IF ( RIGHT ) THEN
            CALL SAXPY( N, ONE, W( 1, Y  ), 1, W( 1, X ), 1 )
         ELSE
            CALL SAXPY( N, ONE, W( 1, UU ), 1, W( 1, X ), 1 )
         END IF
C
C  Start computing the residual
C
         CALL SCOPY( N, W( 1, B ), 1, W( 1, RES ), 1 )
C
C  Restart
C
      GO TO 30
C
C  Save persistent data and return to caller
C
 1000 CONTINUE
      ISAVE(1) = IPOS
      ISAVE(2) = ITMAX
      ISAVE(3) = B
      ISAVE(4) = I
      ISAVE(5) = K
      ISAVE(6) = R
      ISAVE(7) = X
      ISAVE(8) = U
      ISAVE(9) = V
      ISAVE(10) = S
      ISAVE(11) = Y
      ISAVE(12) = CS
      ISAVE(13) = SN
      ISAVE(14) = UU
      ISAVE(15) = RES

      RSAVE(1) = BNRM2
      RSAVE(2) = AA
      RSAVE(3) = BB
      RSAVE(4) = RNORM
      RSAVE(5) = PRESID
      RSAVE(6) = RSTOP
      RSAVE(7) = PRSTOP

      LSAVE(1) = LEFT
      LSAVE(2) = RIGHT
      RETURN
C
C  Non-executable statements
C
 2000 FORMAT( / ' Error message from MI24A/AD. INFO(1) = ', I4 )
 2010 FORMAT( / ' Warning message from MI24A/AD. INFO(1) = ', I4 )
 2020 FORMAT( ' Convergence tolerance out of range.' )
 2030 FORMAT( /, ' # iterations required exceeds the maximum of ',
     *        I8, ' allowed by ICNTL(6)' )
C
C  End of MI24A/AD
C
      END
C COPYRIGHT (c) 1995 Council for the Central Laboratory
*                    of the Research Councils
C Original date 29 March 2001
C  March 2001: threadsafe version of MI06 and extend control arguments

C 12th July 2004 Version 1.0.0. Version numbering added.
C 22nd February 2005 Version 1.1.0. FD05 dependence changed to FD15.

      SUBROUTINE MI26I(ICNTL,CNTL,ISAVE,RSAVE)
C
C  MI26 solves the linear system Ax = b using the
C  BiConjugate Gradient Stabilized iterative method with
C  optionally using preconditioning
C                    ^                   ^
C           P(L)AP(R)x = P(L)b,  x = P(R)x

C  P(L), P(R) are the preconditioners, which are not passed to the code,
C  but each time a matrix-vector product with P = P(L)P(R)
C  is required, control is passed back to the user.
C
C  Similarly, the matrix A is not passed to the code, but when
C  a matrix-vector with A is required, control is
C  passed back to the user.
C
C  MI26I/ID is the initialisation routine for MI26A/AD and should
C  be called once prior to calls to MI26A/AD.
C
C  Argument list.
C
C  ICNTL   (output) INTEGER control array, dimension 8.
C          ICNTL(1) is the stream number for error messages.
C          On exit, ICNTL(2) = 6.
C          ICNTL(2) is the stream number for warning messages.
C          On exit, ICNTL(2) = 6.
C          ICNTL(3) indicates whether the user wishes to use a
C          preconditioner. If ICNTL(3) is nonzero, preconditioning.
C          On exit, ICNTL(3) = 0
C          ICNTL(4) indicates whether
C          the convergence test offered by MI26A/AD is to be used.
C          If ICNTL(4) = 0, the computed solution x is accepted
C          if ||Ax - b||/||r sub 0|| < CNTL(1) (r sub 0 is
C          initial residual).
C          Otherwise, the user may perform his/her
C          own test for convergence when IACT = 1 is returned.
C          On exit, ICNTL(4) = 0
C          ICNTL(5) indicates whether the user wishes to supply an
C          initial guess for the solution vector x.
C          If ICNTL(5) = 0, the user does not wish to supply
C          an initial guess and x = (0,0,...,0) will be used
C          as the initial guess. Otherwise, the user
C          must supply an intial guess on the first call to
C          MI26A/AD. On exit, ICNTL(5) = 0.
C          ICNTL(6) determines the maximum number of iterations
C          allowed. It has default value -1 and, in this case,
C          the maximum number will be N. If the user does
C          not want the maximum number to be N, ICNTL(6) should
C          be set to the maximum  number the user wishes
C          to allow.
C          ICNTRL(7) and ICNTRL(8) are spare and set to zero.
C
C  CNTL    (output)  REAL (DOUBLE PRECISION) control array, dimension 5.
C          CNTL(1) and CNTL(2) are convergence tolerances.
C          On exit, set to square root of machine precision and 0,
C          respectively.
C          If ICNTL(4) is nonzero, CNTL(1) and CNTL(2) are
C          not accessed by MI26A/AD.
C          CNTL(3) is tolerance used to check whether the algorithm
C          has broken down. On exit, set to machine precision
C          CNTRL(4) and CNTRL(5) are spare and set to zero.
C
C  ISAVE   (output) INTEGER ARRAY, length 14, used to hold the
C          routine's persistent integer data.  The array contents
C          must not be altered by the user.
C
C  RSAVE   (output) REAL ARRAY, length 9, used to hold the
C          routine's persistent real data.  The array contents
C          must not be altered by the user.
C
C     .. Array Arguments ..
      REAL CNTL(5)
      INTEGER ICNTL(8)
      INTEGER ISAVE(14)
      REAL RSAVE(9)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E+0)
C     ..
C     .. External Functions ..
      REAL FD15A
      EXTERNAL FD15A
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = 0
      ICNTL(4) = 0
      ICNTL(5) = 0
      ICNTL(6) = -1

      ICNTL(7) = 0
      ICNTL(8) = 0

      CNTL(1) = SQRT(FD15A('E'))
      CNTL(2) = ZERO
      CNTL(3) = FD15A('E')

      CNTL(4) = ZERO
      CNTL(5) = ZERO

C  Initialize persistent data to avoid undefined assignment
      DO 10 I = 1, 14
      ISAVE(I) = 0
   10 CONTINUE
      DO 20 I = 1, 9
      RSAVE(I) = 0.0
   20 CONTINUE

      RETURN
      END
C
      SUBROUTINE MI26A(IACT,N,W,LDW,LOCY,LOCZ,RESID,ICNTL,CNTL,INFO,
     +                 ISAVE,RSAVE)
C
C     .. Scalar Arguments ..
      REAL RESID
      INTEGER IACT,LDW,LOCY,LOCZ,N
C     ..
C     .. Array Arguments ..
      REAL CNTL(5),W(LDW,8)
      INTEGER ICNTL(8),INFO(4)
      INTEGER ISAVE(14)
      REAL RSAVE(9)
C     ..
C
C  Argument list.
C
C  IACT    (input) INTEGER.
C          IACT must be set to 0 prior to first call.
C          On each exit, IACT indicates the action required by
C          the user. Possible values of IACT and the action
C          required are as follows:
C  -1      Fatal error (see INFO(1)). Terminate computation.
C   1      If ICNTL(4) = 0 (the default), convergence has been
C          achieved and the user should terminate the computation.
C          If ICNTL(4) is nonzero, the user should test the norm of
C          the residual in RESID for convergence and recall MI26A/AD
C          if convergence has not been achieved.
C   2      The user must perform the matrix-vector product
C          y := Az ,
C          and recall MI26A/AD. The vectors y and z
C          are  held in the first N entries of
C          columns LOCY and LOCZ of array W, respectively.
C   3      The user must perform the preconditioning operations
C          y := Pz,
C          where P = P sub L P sub R  is the preconditioner
C          and recall MI26A/AD. The vectors y and z
C          are  held in the first N entries of
C          columns LOCY and LOCZ of array W, respectively.
C
C  N       (input) INTEGER.
C          On entry, the dimension of the matrix.
C          Unchanged on exit.
C
C  W       (right-hand side, solution and workspace, input/output)
C          REAL (DOUBLE PRECISION) array, dimension (LDW,8).
C          Prior to the first call, the first N entries of column
C          1 must be set to hold the right-hand side vector b and,
C          if ICNTL(5) is nonzero, the first N entries of column 2
C          must be set to the initial estimate of the solution vector
C          x.  On exit with IACT = 1, the first N entries of column 1
C          hold the current residual vector r = b - Ax, and the
C          current estimate of the solution x is held in the
C          first N entries of column 2.  On exit
C          with IACT > 1, the user is required to perform a
C          computation with columns LOCZ and LOCZ of W
C          (see argument IACT).  The remaining columns of
C          W must not be altered by the user between calls to
C          MI26A/AD.
C
C  LDW     (input) INTEGER
C          The leading dimension of the array W. LDW >= max(1,N).
C
C LOCY, LOCZ (output) INTEGER  variables
C          On exit with IACT > 1, LOCY and LOCZ define
C          the columns of the array W which hold y and z.
C          (see IACT).
C
C  RESID   (output) REAL (DOUBLE PRECISION)
C          On exit with IACT = 1,
C          RESID holds ||b - Ax||, where x is the
C          iterated solution.
C          If ICNTL(4) is nonzero, on exit with IACT = 1,
C          the user may carry out his/her
C          own test for convergnce at this point.
C
C  ICNTL   (input) INTEGER control array of dimension 8.
C          ICNTL may be initialised by calling MI26I/ID.
C          See MI26I/ID for details.
C
C  CNTL    (input) REAL (DOUBLE PRECISION) control array of dimension 5.
C          CNTL may be initialised by calling MI26I/ID.
C          See MI26I/ID for details.
C
C  INFO    (output) INTEGER ARRAY , length 4.
C          If INFO(1) = 0 on exit, no errors or warnings issued.
C          If INFO(1) > 0 on exit, warning to user.
C          INFO(1) = 1, value of CNTL(1) is out-of-range (u,1.0)
C          The default sqrt(u) is used, u=machine precision.
C          If INFO(1) < 0 on exit, illegal input parameter,
C               or breakdown occured during iteration.
C
C                Error input data:
C
C                   -1: matrix dimension N < 0
C                   -2: LDW < N
C
C                BREAKDOWN: If RHO becomes too small,
C                   the program will terminate.
C
C                   -3: RHO < CNTL(3)*||R|*||RTLD||
C          R and RTLD have become  orthogonal.
C          Error -3 also returned if OMEGA becomes too small
C          (S and T have become orthogonal relative to T'*T).
C
C          On each exit, INFO(2) holds the number of
C          iterations performed so far.
C
C  ISAVE   (output) INTEGER ARRAY, length 14, used to hold the
C          routine's persistent integer data.  The array contents
C          must not be altered by the user.
C
C  RSAVE   (output) REAL ARRAY, length 9, used to hold the
C          routine's persistent real data.  The array contents
C          must not be altered by the user.
C
C  BLAS CALLS: SAXPY, SCOPY, SDOT, SNRM2, SSCAL
C***************************************************
C
C     .. Parameters ..
      REAL ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
C     ..
C     .. Local Scalars ..
      REAL ALPHA,BETA,BNRM2,OMEGA,RHO,RHO1,RNRM2,RTNRM2,RSTOP,SNORM2,
     +     TNORM2
      INTEGER B,I,IPOS,ITMAX,P,PHAT,R,RTLD,S,SHAT,T,V,X
C     ..
C     .. External Functions ..
      REAL SDOT,SNRM2,FD15A
      EXTERNAL SDOT,SNRM2,FD15A
C     ..
C     .. Intrinsic Functions ..

      INTRINSIC ABS,MAX,SQRT
C     ..
C     .. External Subroutines ..
      EXTERNAL SAXPY,SCOPY,SSCAL
C
C  Restore persistent data
C
      IPOS   = ISAVE(1)
      ITMAX  = ISAVE(2)
      B      = ISAVE(3)
      R      = ISAVE(4)
      X      = ISAVE(5)
      P      = ISAVE(6)
      S      = ISAVE(7)
      T      = ISAVE(8)
      V      = ISAVE(9)
      PHAT   = ISAVE(10)
      RTLD   = ISAVE(11)
      SHAT   = ISAVE(12)

      BNRM2  = RSAVE(1)
      ALPHA  = RSAVE(2)
      BETA   = RSAVE(3)
      RHO    = RSAVE(4)
      RHO1   = RSAVE(5)
      RSTOP  = RSAVE(6)
      OMEGA  = RSAVE(7)

C Jump to appropriate place in code
      IF (IACT.EQ.0) GO TO 10
C Immediate return if error on a previous call
      IF (IACT.LT.0) GO TO 1000
C Immediate return if convergence already achieved
      IF (IACT.EQ.1 .AND. ICNTL(4).EQ.0) GO TO 1000
      IF (IACT.EQ.1 .AND. BNRM2.EQ.ZERO) GO TO 1000
C
      IF (IPOS.EQ.1) GO TO 40
      IF (IPOS.EQ.2) GO TO 70
      IF (IPOS.EQ.3) GO TO 80
      IF (IPOS.EQ.4) GO TO 90
      IF (IPOS.EQ.5) GO TO 100
      IF (IPOS.EQ.6) GO TO 110
      IF (IPOS.EQ.7) GO TO 120


   10 CONTINUE
C
C  Initial call.
      INFO(1) = 0
C  Test the input parameters.
      IF (N.LE.0) THEN
         INFO(1) = -1
      ELSE IF (LDW.LT.MAX(1,N)) THEN
         INFO(1) = -2
      END IF
      IF (INFO(1).LT.0) THEN
         IACT = -1
         IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9000) INFO(1)
         GO TO 1000
      END IF

C
C Alias workspace columns.
C
      B = 1
      X = 2
      R = 1
      RTLD = 3
      P = 4
      V = 5
      T = 6
      PHAT = 7
      SHAT = 8
      S = 1
C Set INFO(2) and ITMAX.

      INFO(2) = 0
      ITMAX = N
      IF (ICNTL(6).GT.0) ITMAX = ICNTL(6)

C Compute ||b||
C
      BNRM2 = SNRM2(N,W(1,B),1)
C Immediate return if ||b|| = 0.
      IF (BNRM2.EQ.ZERO) THEN
         IACT = 1
         DO 20 I = 1,N
            W(I,X) = ZERO
            W(I,B) = ZERO
   20    CONTINUE
         RESID = ZERO
         GO TO 1000
      END IF
C
      IF (ICNTL(4).EQ.0) THEN
C Check value of CNTL(1)
         IF (CNTL(1).LT.FD15A('E') .OR. CNTL(1).GT.ONE) THEN
            INFO(1) = 1
            IF (ICNTL(2).GT.0) THEN
               WRITE (ICNTL(2),FMT=9010) INFO(1)
               WRITE (ICNTL(2),FMT=9020)
            END IF
            CNTL(1) = SQRT(FD15A('E'))
         END IF
      END IF
C
C Compute initial residual.
C
C If the user has not supplied an initial guess, set x = 0
C as the initial guess.
      IF (ICNTL(5).EQ.0) THEN
         DO 30 I = 1,N
            W(I,X) = ZERO
   30    CONTINUE
         GO TO 50
      ELSE
C Initial guess supplied by user
C If initial guess for solution is x = 0 no action is required.
C (r = b)
         IF (SNRM2(N,W(1,X),1).EQ.ZERO) GO TO 50
C
C Otherwise, return to user to compute Ax.
C Column P can be used temporarily to hold Ax.
         IPOS = 1
         IACT = 2
         LOCY = P
         LOCZ = X
         GO TO 1000
      END IF
C
C Compute r = b - Ax
   40 CALL SAXPY(N,-ONE,W(1,P),1,W(1,R),1)

   50 CONTINUE
C
C  Compute the norm of the initial residual
C
      RSTOP = SNRM2(N,W(1,R),1)
C
C Choose RTLD such that initially, (R,RTLD) = RHO is not equal to 0.
C Here we choose RTLD = R.
C
      CALL SCOPY(N,W(1,R),1,W(1,RTLD),1)

C Perform BiConjugate Gradient Stabilized iteration.

   60 CONTINUE

C Update iteration count

      INFO(2) = INFO(2) + 1
C
Check maximum number of iterations has not been exceeded.
      IF (INFO(2).GT.ITMAX) THEN
         INFO(1) = -4
         IACT = -1
         IF (ICNTL(1).GT.0) THEN
            WRITE (ICNTL(1),FMT=9000) INFO(1)
            WRITE (ICNTL(1),FMT=9030) ITMAX
         END IF
         GO TO 1000
      END IF


      RHO = SDOT(N,W(1,RTLD),1,W(1,R),1)
C Check for breakdown.
      IF (ABS(RHO).LT.CNTL(3)*N) THEN
C RHO getting small. Carry out more rigorous test
         RNRM2 = SNRM2(N,W(1,R),1)
         RTNRM2 = SNRM2(N,W(1,RTLD),1)
         IF (ABS(RHO).LT.CNTL(3)*RNRM2*RTNRM2) THEN
            INFO(1) = -3
            IACT = -1
            IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9000) INFO(1)
            GO TO 1000
         END IF
      END IF
C
C Compute vector P.
C
      IF (INFO(2).GT.1) THEN
         BETA = (RHO/RHO1)* (ALPHA/OMEGA)
         CALL SAXPY(N,-OMEGA,W(1,V),1,W(1,P),1)
         CALL SSCAL(N,BETA,W(1,P),1)
         CALL SAXPY(N,ONE,W(1,R),1,W(1,P),1)
      ELSE
         CALL SCOPY(N,W(1,R),1,W(1,P),1)
      END IF
C
C Compute direction adjusting vector PHAT and scalar ALPHA.
C
      IF (ICNTL(3).NE.0) THEN
C Return to user for preconditioning.
         IPOS = 2
         IACT = 3
         LOCY = PHAT
         LOCZ = P
         GO TO 1000
      ELSE
C No preconditioning (i.e identity preconditioning)
         CALL SCOPY(N,W(1,P),1,W(1,PHAT),1)
      END IF

   70 CONTINUE
C
C Return to user for matrix-vector product
      IPOS = 3
      IACT = 2
      LOCY = V
      LOCZ = PHAT
      GO TO 1000

   80 CONTINUE

      ALPHA = RHO/SDOT(N,W(1,RTLD),1,W(1,V),1)
C
C Early check for tolerance.
C
      CALL SAXPY(N,-ALPHA,W(1,V),1,W(1,R),1)

C Note: R=1 and S=1, so no need to copy R into S
C        CALL SCOPY( N, W(1,R), 1, W(1,S), 1 )

      CALL SAXPY(N,ALPHA,W(1,PHAT),1,W(1,X),1)

      RESID = SNRM2(N,W(1,S),1)
      IPOS = 4
      IF (ICNTL(4).NE.0) THEN
C Return the residual to the user for convergence testing.
         IACT = 1
         GO TO 1000
      ELSE
C Test the scaled residual for convergence.
         IF (RESID.LE.MAX(CNTL(2),RSTOP*CNTL(1))) THEN
C Convergence achieved
            IACT = 1
            GO TO 1000
         END IF
      END IF
C
   90 CONTINUE

C
C  Compute stabilizer vector SHAT and scalar OMEGA.
C
      IF (ICNTL(3).NE.0) THEN
C Return to user for preconditioning.
         IPOS = 5
         IACT = 3
         LOCY = SHAT
         LOCZ = S
         GO TO 1000
      ELSE
C No preconditioning (i.e identity preconditioning)
         CALL SCOPY(N,W(1,S),1,W(1,SHAT),1)
      END IF

  100 CONTINUE

C Return to user for matrix-vector product
      IPOS = 6
      IACT = 2
      LOCY = T
      LOCZ = SHAT
      GO TO 1000

  110 CONTINUE

      OMEGA = SDOT(N,W(1,T),1,W(1,S),1)/SDOT(N,W(1,T),1,W(1,T),1)

C Check OMEGA is not too small.
      IF (ABS(OMEGA).LT.CNTL(3)*N) THEN
C OMEGA getting small. Carry out more rigorous test
         SNORM2 = SNRM2(N,W(1,S),1)
         TNORM2 = SNRM2(N,W(1,T),1)
         IF (ABS(RHO).LT.CNTL(3)*SNORM2/TNORM2) THEN
            INFO(1) = -3
            IACT = -1
            IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9000) INFO(1)
            GO TO 1000
         END IF
      END IF
C
C Compute new solution approximation vector X.
C
      CALL SAXPY(N,OMEGA,W(1,SHAT),1,W(1,X),1)

C  Compute residual R, check for tolerance.
C
      CALL SAXPY(N,-OMEGA,W(1,T),1,W(1,R),1)

      RESID = SNRM2(N,W(1,R),1)
      IPOS = 7
      IF (ICNTL(4).NE.0) THEN
C Return the residual to the user for convergence testing.
         IACT = 1
         GO TO 1000
      ELSE
C Test the scaled residual for convergence.
         IF (RESID.LE.MAX(CNTL(2),RSTOP*CNTL(1))) THEN
C Convergence achieved
            IACT = 1
            GO TO 1000
         END IF
      END IF
C
  120 CONTINUE

      RHO1 = RHO

C Next iteration
      GO TO 60
C
C  Save persistent data and return to caller
C
 1000 CONTINUE
      ISAVE(1)  = IPOS
      ISAVE(2)  = ITMAX
      ISAVE(3)  = B
      ISAVE(4)  = R
      ISAVE(5)  = X
      ISAVE(6)  = P
      ISAVE(7)  = S
      ISAVE(8)  = T
      ISAVE(9)  = V
      ISAVE(10) = PHAT
      ISAVE(11) = RTLD
      ISAVE(12) = SHAT

      RSAVE(1)  = BNRM2
      RSAVE(2)  = ALPHA
      RSAVE(3)  = BETA
      RSAVE(4)  = RHO
      RSAVE(5)  = RHO1
      RSAVE(6)  = RSTOP
      RSAVE(7)  = OMEGA
      RETURN
C
C End of MI26A/AD
C
 9000 FORMAT (/' Error message from MI26A/AD. INFO(1) = ',I4)
 9010 FORMAT (/' Warning message from MI26A/AD. INFO(1) = ',I4)
 9020 FORMAT (' Convergence tolerance out of range.')
 9030 FORMAT (' Number of iterations required exceeds the maximum of ',
     +       I8,/' allowed by ICNTL(6)')
      END
