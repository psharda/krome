!THESE ARE THE SUBROUTINES FOR NLEQ1 STORED IN A SINGLE FILE
!NLEQ1 SOLVES A NON-LINEAR SYSTEM OF EQUATIONS
!COMMENTED SUBROUTINES CONFLICTING WITH DLSODES
!FURTHER DETAILS: !http://www.zib.de/en/numerik/software/ant/nleq1.html

      SUBROUTINE NLEQ1(N,FCN,JAC,X,XSCAL,RTOL,IOPT,IERR,
     $LIWK,IWK,LRWK,RWK)
C*    Begin Prologue NLEQ1
      INTEGER N
      EXTERNAL FCN,JAC
      DOUBLE PRECISION X(N),XSCAL(N)
      DOUBLE PRECISION RTOL
      INTEGER IOPT(50)
      INTEGER IERR
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LRWK
      DOUBLE PRECISION RWK(LRWK)
C     ------------------------------------------------------------
C
C*  Title
C
C     Numerical solution of nonlinear (NL) equations (EQ)
C     especially designed for numerically sensitive problems.
C
C*  Written by        U. Nowak, L. Weimann 
C*  Purpose           Solution of systems of highly nonlinear equations
C*  Method            Damped affine invariant Newton method
C                     (see references below)
C*  Category          F2a. - Systems of nonlinear equations
C*  Keywords          Nonlinear equations, Newton methods
C*  Version           2.4
C*  Revision          May 2009
C*  Latest Change     May 2009
C*  Library           CodeLib
C*  Code              Fortran 77, Double Precision
C*  Environment       Standard Fortran 77 environment on PC's,
C                     workstations and hosts.
C*  Copyright     (c) Konrad-Zuse-Zentrum fuer
C                     Informationstechnik Berlin (ZIB)
C                     Takustrasse 7, D-14195 Berlin-Dahlem
C                     phone : + 49/30/84185-0
C                     fax   : + 49/30/84185-125
C*  Contact           Lutz Weimann
C                     ZIB, Division Scientific Computing, 
C                          Department Numerical Analysis and Modelling
C                     phone : + 49/30/84185-185
C                     fax   : + 49/30/84185-107
C                     e-mail: weimann@zib.de
C
C*    References:
C
C     /1/ P. Deuflhard:
C         Newton Methods for Nonlinear Problems. -
C         Affine Invariance and Adaptive Algorithms.
C         Series Computational Mathematics 35, Springer (2004)
C
C     /2/ U. Nowak, L. Weimann:
C         A Family of Newton Codes for Systems of Highly Nonlinear
C         Equations - Algorithm, Implementation, Application.
C         ZIB, Technical Report TR 90-10 (December 1990)
C
C  ---------------------------------------------------------------
C
C* Licence
C    You may use or modify this code for your own non commercial
C    purposes for an unlimited time. 
C    In any case you should not deliver this code without a special 
C    permission of ZIB.
C    In case you intend to use the code commercially, we oblige you
C    to sign an according licence agreement with ZIB.
C
C* Warranty 
C    This code has been tested up to a certain level. Defects and
C    weaknesses, which may be included in the code, do not establish
C    any warranties by ZIB. ZIB does not take over any liabilities
C    which may follow from acquisition or application of this code.
C
C* Software status 
C    This code is under care of ZIB and belongs to ZIB software class 1.
C
C     ------------------------------------------------------------
C
C*    Summary:
C     ========
C     Damped Newton-algorithm for systems of highly nonlinear
C     equations - damping strategy due to Ref. (1).
C
C     (The iteration is done by subroutine N1INT currently. NLEQ1
C      itself does some house keeping and builds up workspace.)
C
C     Jacobian approximation by numerical differences or user
C     supplied subroutine JAC.
C
C     The numerical solution of the arising linear equations is
C     done by means of the subroutines *GETRF and *GETRS ( Gauss-
C     algorithm with column-pivoting and row-interchange ) in the
C     dense matrix case, or by the subroutines *GBTRF and *GBTRS in
C     the band matrix case from LAPACK (replace '*' by 'S' or 'D'
C     for single or double precision version respectively).
C     For special purposes these routines may be substituted.
C
C     This is a driver routine for the core solver N1INT.
C
C     ------------------------------------------------------------
C
C*    Parameters list description (* marks inout parameters)
C     ======================================================
C
C*    External subroutines (to be supplied by the user)
C     =================================================
C 
C     (Caution: Arguments declared as (input) must not
C               be altered by the user subroutines ! )
C
C     FCN(N,X,F,IFAIL) Ext    Function subroutine
C       N              Int    Number of vector components (input)
C       X(N)           Dble   Vector of unknowns (input)
C       F(N)           Dble   Vector of function values (output)
C       IFAIL          Int    FCN evaluation-failure indicator. (output)
C                             On input:  Has always value 0 (zero).
C                             On output: Indicates failure of FCN eval-
C                                uation, if having a value <= 2.
C                             If <0: NLEQ1 will be terminated with 
C                                    error code = 82, and IFAIL stored
C                                    to IWK(23).
C                             If =1: A new trial Newton iterate will
C                                    computed, with the damping factor
C                                    reduced to it's half.
C                             If =2: A new trial Newton iterate will
C                                    computed, with the damping factor
C                                    reduced by a reduct. factor, which
C                                    must be output through F(1) by FCN,
C                                    and it's value must be >0 and < 1.
C                             Note, that if IFAIL = 1 or 2, additional
C                             conditions concerning the damping factor,
C                             e.g. the minimum damping factor or the
C                             bounded damping strategy may also influ-
C                             ence the value of the reduced damping 
C                             factor.
C
C     JAC(N,LDJAC,X,DFDX,IFAIL) 
C                       Ext    Jacobian matrix subroutine
C       N                 Int    Number of vector components (input)
C       LDJAC             Int    Leading dimension of Jacobian array 
C                                (input)
C       X(N)              Dble   Vector of unknowns (input)
C       DFDX(LDJAC,N)     Dble   DFDX(i,k): partial derivative of
C                                I-th component of FCN with respect
C                                to X(k) (output)
C       IFAIL             Int    JAC evaluation-failure indicator.
C                                (output)
C                                Has always value 0 (zero) on input.
C                                Indicates failure of JAC evaluation
C                                and causes termination of NLEQ1,
C                                if set to a negative value on output
C
C
C*    Input parameters of NLEQ1
C     =========================
C
C     N              Int    Number of unknowns
C   * X(N)           Dble   Initial estimate of the solution
C   * XSCAL(N)       Dble   User scaling (lower threshold) of the 
C                           iteration vector X(N)
C   * RTOL           Dble   Required relative precision of
C                           solution components -
C                           RTOL.GE.EPMACH*TEN*N
C   * IOPT(50)       Int    Array of run-time options. Set to zero
C                           to get default values (details see below)
C
C*    Output parameters of NLEQ1
C     ==========================
C
C   * X(N)           Dble   Solution values ( or final values,
C                           respectively )
C   * XSCAL(N)       Dble   After return with IERR.GE.0, it contains
C                           the latest internal scaling vector used
C                           After return with IERR.EQ.-1 in onestep-
C                           mode it contains a possibly adapted 
C                           (as described below) user scaling vector:
C                           If (XSCAL(I).LT. SMALL) XSCAL(I) = SMALL ,
C                           If (XSCAL(I).GT. GREAT) XSCAL(I) = GREAT .
C                           For SMALL and GREAT, see section machine
C                           constants below  and regard note 1.
C   * RTOL           Dble   Finally achieved (relative) accuracy.
C                           The estimated absolute error of component i
C                           of x_out is approximately given by
C                             abs_err(i) = RTOL * XSCAL_out(i) ,
C                           where (approximately)
C                             XSCAL_out(i) = 
C                                max(abs(X_out(i)),XSCAL_in(i)).
C     IERR           Int    Return value parameter
C                           =-1 sucessfull completion of one iteration
C                               step, subsequent iterations are needed 
C                               to get a solution. (stepwise mode only) 
C                           = 0 successfull completion of the iteration,
C                               solution has been computed
C                           > 0 see list of error messages below
C
C     Note 1.
C        The machine dependent values SMALL and EPMACH are
C        gained from calls of the ZIBCONST.
C
C*    Workspace parameters of NLEQ1
C     =============================
C
C     LIWK           Int    Declared dimension of integer workspace.
C                           Required minimum (for standard linear system
C                           solver) : N+50
C  *  IWK(LIWK)      Int    Integer Workspace
C     LRWK           Int    Declared dimension of real workspace.
C                           Required minimum (for standard linear system
C                           solver and Jacobian computed by numerical
C                           approximation - if the Jacobian is computed
C                           by a user subroutine JAC, decrease the 
C                           expressions noted below by N):
C                           for full case Jacobian: (N+NBROY+13)*N+61
C                           for a band-matrix Jacobian:
C                              (2*ML+MU+NBROY+14)*N+61 with
C                           ML = lower bandwidth , MU = upper bandwidth
C                           NBROY = Maximum number of Broyden steps
C                           (Default: if Broyden steps are enabled, e.g.
C                                                IOPT(32)=1            -
C                                       NBROY=N (full Jacobian), 
C                                            =ML+MU+1 (band Jacobian),
C                                       but at least NBROY=10 
C                                     else (if IOPT(32)=0) - 
C                                       NBROY=0 ;
C                            see equally named IOPT and IWK-fields below)
C   * RWK(LRWK)      Dble   Real Workspace
C
C     Note 2a.  A test on sufficient workspace is made. If this
C               test fails, IERR is set to 10 and an error-message
C               is issued from which the minimum of required
C               workspace size can be obtained.
C
C     Note 2b.  The first 50 elements of IWK and RWK are partially 
C               used as input for internal algorithm parameters (for
C               details, see below). In order to set the default values
C               of these parameters, the fields must be set to zero.
C               Therefore, it's recommanded always to initialize the
C               first 50 elements of both workspaces to zero.
C
C*   Options IOPT:
C    =============
C
C     Pos. Name   Default  Meaning
C
C       1  QSUCC  0        =0 (.FALSE.) initial call:
C                             NLEQ1 is not yet initialized, i.e. this is
C                             the first call for this nonlinear system.
C                             At successfull return with MODE=1,
C                             QSUCC is set to 1.
C                          =1 (.TRUE.) successive call:
C                             NLEQ1 is initialized already and is now
C                             called to perform one or more following
C                             Newton-iteration steps.
C                             ATTENTION:
C                                Don't destroy the contents of
C                                IOPT(i) for 1 <= i <= 50 ,
C                                IWK(j)  for 1 <= j < NIWKFR and
C                                RWK(k)  for 1 <= k < NRWKFR.
C                                (Nevertheless, some of the options, e.g.
C                                 FCMIN, SIGMA, MPR..., can be modified
C                                 before successive calls.)
C       2  MODE   0        =0 Standard mode initial call:
C                             Return when the required accuracy for the
C                             iteration vector is reached. User defined
C                             parameters are evaluated and checked.
C                             Standard mode successive call:
C                             If NLEQ1 was called previously with MODE=1,
C                             it performs all remaining iteration steps.
C                          =1 Stepwise mode:
C                             Return after one Newton iteration step.
C       3  JACGEN 0        Method of Jacobian generation
C                          =0 Standard method is JACGEN=2
C                          =1 User supplied subroutine JAC will be 
C                             called to generate Jacobian matrix
C                          =2 Jacobian approximation by numerical
C                             differentation (see subroutines N1JAC
C                             and N1JACB)
C                          =3 Jacobian approximation by numerical
C                             differentation with feedback control
C                             (see subroutines N1JCF and N1JCFB)
C       4  MSTOR  0        =0 The Jacobian A is a dense matrix
C                          =1 A is a band matrix
C       5                  Reserved
C       6  ML     0        Lower bandwidth of A (excluding the
C                          diagonal);
C                          IOPT(6) ignored, if IOPT(4).NE.1
C       7  MU     0        Upper bandwidth of A (excluding the
C                          diagonal);
C                          IOPT(7) ignored, if IOPT(4).NE.1
C       8                  Reserved
C       9  ISCAL  0        Determines how to scale the iterate-vector:
C                          =0 The user supplied scaling vector XSCAL is
C                             used as a (componentwise) lower threshold
C                             of the current scaling vector
C                          =1 The vector XSCAL is always used as the
C                             current scaling vector
C      10                  Reserved
C      11  MPRERR 0        Print error messages
C                          =0 No output
C                          =1 Error messages
C                          =2 Warnings additionally
C                          =3 Informal messages additionally
C      12  LUERR  6        Logical unit number for error messages
C      13  MPRMON 0        Print iteration Monitor
C                          =0 No output
C                          =1 Standard output
C                          =2 Summary iteration monitor additionally
C                          =3 Detailed iteration monitor additionally
C                          =4,5,6 Outputs with increasing level addi-
C                             tional increasing information for code
C                             testing purposes. Level 6 produces
C                             in general extremely large output!
C      14  LUMON  6        Logical unit number for iteration monitor
C      15  MPRSOL 0        Print solutions
C                          =0 No output
C                          =1 Initial values and solution values
C                          =2 Intermediate iterates additionally
C      16  LUSOL  6        Logical unit number for solutions
C      17..18              Reserved
C      19  MPRTIM 0        Output level for the time monitor
C                          = 0 : no time measurement and no output
C                          = 1 : time measurement will be done and
C                                summary output will be written -
C                                regard note 5a.
C      20  LUTIM  6        Logical output unit for time monitor
C      21..30              Reserved
C      31  NONLIN 3        Problem type specification
C                          =1 Linear problem
C                             Warning: If specified, no check will be
C                             done, if the problem is really linear, and
C                             NLEQ1 terminates unconditionally after one
C                             Newton-iteration step.
C                          =2 Mildly nonlinear problem
C                          =3 Highly nonlinear problem
C                          =4 Extremely nonlinear problem
C      32  QRANK1 0        =0 (.FALSE.) Rank-1 updates by Broyden-
C                             approximation are inhibited.
C                          =1 (.TRUE.) Rank-1 updates by Broyden-
C                             approximation are allowed.
C      33  QORDI  0        =0 (.FALSE.) Standard program mode 
C                          =1 (.TRUE.)  Special program mode:
C                             Ordinary Newton iteration is done, e.g.:
C                             No damping strategy and no monotonicity
C                             test is applied
C      34  QSIMPL 0        =0 (.FALSE.) Standard program mode
C                          =1 (.TRUE.)  Special program mode:
C                             Simplified Newton iteration is done, e.g.:
C                             The Jacobian computed at the starting
C                             point is fixed for all subsequent
C                             iteration steps, and
C                             no damping strategy and no monotonicity
C                             test is applied.
C      35  QNSCAL 0        Inhibit automatic row scaling: 
C                          =0 (.FALSE.) Automatic row scaling of
C                             the linear system is activ: 
C                             Rows i=1,...,N will be divided by
C                             max j=1,...,N (abs(a(i,j))) 
C                          =1 (.TRUE.) No row scaling of the linear
C                             system. Recommended only for well row-
C                             scaled nonlinear systems.
C      36..37              Reserved
C      38  IBDAMP          Bounded damping strategy switch:
C                          =0 The default switch takes place, dependent
C                             on the setting of NONLIN (=IOPT(31)):
C                             NONLIN = 0,1,2,3 -> IBDAMP = off ,
C                             NONLIN = 4 -> IBDAMP = on
C                          =1 means always IBDAMP = on 
C                          =2 means always IBDAMP = off
C      39  IORMON          Convergence order monitor 
C                          =0 Standard option is IORMON=2 
C                          =1 Convergence order is not checked,
C                             the iteration will be always proceeded
C                             until the solution has the required 
C                             precision RTOL (or some error condition
C                             occured)
C                          =2 Use additional 'weak stop' criterion:
C                             Convergence order is monitored
C                             and termination due to slowdown of the
C                             convergence may occur.
C                          =3 Use additional 'hard stop' criterion:
C                             Convergence order is monitored
C                             and termination due to superlinear 
C                             convergence slowdown may occur. 
C                          In case of termination due to convergence
C                          slowdown, the warning code IERR=4 will be
C                          set.
C                          In cases, where the Newton iteration con-
C                          verges but superlinear convergence order has
C                          never been detected, the warning code IERR=5 
C                          is returned.
C      40..45              Reserved
C      46..50              User options (see note 5b)
C
C     Note 3:
C         If NLEQ1 terminates with IERR=2 (maximum iterations)
C         or  IERR=3 (small damping factor), you may try to continue
C         the iteration by increasing NITMAX or decreasing FCMIN
C         (see RWK) and setting QSUCC to 1.
C
C     Note 4 : Storage of user supplied banded Jacobian
C        In the band matrix case, the following lines may build
C        up the analytic Jacobian A;
C        Here AFL denotes the quadratic matrix A in dense form,
C        and ABD the rectangular matrix A in banded form :
C
C                   ML = IOPT(6)
C                   MU = IOPT(7)
C                   MH = MU+1
C                   DO 20 J = 1,N
C                     I1 = MAX0(1,J-MU)
C                     I2 = MIN0(N,J+ML)
C                     DO 10 I = I1,I2
C                       K = I-J+MH
C                       ABD(K,J) = AFL(I,J)
C           10        CONTINUE
C           20      CONTINUE
C
C         The total number of rows needed in  ABD  is  ML+MU+1 .
C         The  MU by MU  upper left triangle and the
C         ML by ML  lower right triangle are not referenced.
C
C     Note 5a:
C        The integrated time monitor calls the machine dependent
C        subroutine SECOND to get the current time stamp in form
C        of a real number (Single precision). As delivered, this
C        subroutine always return 0.0 as time stamp value. Refer
C        to the compiler- or library manual of the FORTRAN compiler
C        which you currently use to find out how to get the current
C        time stamp on your machine.
C
C     Note 5b:
C         The user options may be interpreted by the user replacable
C         routines N1SOUT, N1FACT, N1SOLV - the distributed version
C         of N1SOUT currently uses IOPT(46) as follows:
C         0 = standard plotdata output (may be postprocessed by a user-
C             written graphical program)
C         1 = plotdata output is suitable as input to the graphical
C             package GRAZIL (based on GKS), which has been developed
C             at ZIB. 
C
C
C*   Optional INTEGER input/output in IWK:
C    =======================================
C
C     Pos. Name          Meaning
C
C      1   NITER  IN/OUT Number of Newton-iterations
C      2                 reserved
C      3   NCORR  IN/OUT Number of corrector steps
C      4   NFCN   IN/OUT Number of FCN-evaluations
C      5   NJAC   IN/OUT Number of Jacobian generations or
C                        JAC-calls
C      6                 reserved
C      7                 reserved
C      8   NFCNJ  IN/OUT Number of FCN-evaluations for Jacobian
C                        approximation
C      9   NREJR1 IN/OUT Number of rejected Newton iteration steps
C                        done with a rank-1 approximated Jacobian
C     10..11             Reserved
C     12   IDCODE IN/OUT Output: The 8 decimal digits program identi-
C                        fication number ppppvvvv, consisting of the
C                        program code pppp and the version code vvvv.
C                        Input: If containing a negative number,
C                        it will only be overwritten by the identi-
C                        fication number, immediately followed by
C                        a return to the calling program.      
C     13..15             Reserved
C     16   NIWKFR OUT    First element of IWK which is free to be used
C                        as workspace between Newton iteration steps
C                        for standard linear solvers: 51
C     17   NRWKFR OUT    First element of RWK which is free to be used
C                        as workspace between Newton iteration steps.
C                        For standard linear solvers and numerically 
C                        approximated Jacobian computed by one of the 
C                        expressions:
C                        (N+7+NBROY)*N+61        for a full Jacobian
C                        (2*ML+MU+8+NBROY)*N+61  for a banded Jacobian
C                        If the Jacobian is computed by a user routine
C                        JAC, subtract N in both expressions.
C     18   LIWKA  OUT    Length of IWK currently required
C     19   LRWKA  OUT    Length of RWK currently required
C     20..22             Reserved
C     23   IFAIL  OUT    Set in case of failure of N1FACT (IERR=80),
C                        N1SOLV (IERR=81), FCN (IERR=82) or JAC(IERR=83)
C                        to the nonzero IFAIL value returned by the 
C                        routine indicating the failure .
C     24   ICONV  OUT    Current status of of the convergence monitor
C                        (only if convergence order monitor is on - 
C                         see IORMON(=IOPT(39)))
C                        =0: No convergence indicated yet
C                        =1: Damping factor is 1.0d0
C                        =2: Superlinear convergence in progress
C                        =3: Quadratic convergence in progress
C     25..30             Reserved
C     31   NITMAX IN     Maximum number of permitted iteration
C                        steps (default: 50)
C     32                 Reserved
C     33   NEW    IN/OUT Count of consecutive rank-1 updates
C     34..35             Reserved
C     36   NBROY  IN     Maximum number of possible consecutive 
C                        iterative Broyden steps. The total real 
C                        workspace needed (RWK) depends on this value
C                        (see LRWK above).
C                        Default is N (see parameter N),
C                        if MSTOR=0 (=IOPT(4)), 
C                        and ML+MU+1 (=IOPT(6)+IOPT(7)+1), if MSTOR=1
C                        (but minimum is always 10) - 
C                        provided that Broyden is allowed. 
C                        If Broyden is inhibited, NBROY is always set to
C                        zero.
C     37..50             Reserved
C
C*   Optional REAL input/output in RWK:
C    ====================================
C
C     Pos. Name          Meaning
C
C      1..16             Reserved
C     17   CONV   OUT    The achieved relative accuracy after the  
C                        current step
C     18   SUMX   OUT    Natural level (not Normx of printouts)
C                        of the current iterate, i.e. Sum(DX(i)**2),
C                        where DX = scaled Newton correction.
C     19   DLEVF  OUT    Standard level (not Normf of printouts)
C                        of the current iterate, i.e. Norm2(F(X)),
C                        where F =  nonlinear problem function.
C     20   FCBND  IN     Bounded damping strategy restriction factor
C                        (Default is 10)
C     21   FCSTRT IN     Damping factor for first Newton iteration -
C                        overrides option NONLIN, if set (see note 6)
C     22   FCMIN  IN     Minimal allowed damping factor (see note 6)
C     23   SIGMA  IN     Broyden-approximation decision parameter
C                        Required choice: SIGMA.GE.1. Increasing this
C                        parameter make it less probable that the algo-
C                        rith performs rank-1 updates.
C                        Rank-1 updates are inhibited, if 
C                        SIGMA.GT.1/FCMIN is set. (see note 6)
C     24   SIGMA2 IN     Decision parameter about increasing damping
C                        factor to corrector if predictor is small.
C                        Required choice: SIGMA2.GE.1. Increasing this
C                        parameter make it less probable that the algo-
C                        rith performs rank-1 updates.
C     25                 Reserved
C     26   AJDEL  IN     Jacobian approximation without feedback:
C                        Relative pertubation for components
C                        (Default: sqrt(epmach*10), epmach: relative
C                         machine precision) 
C     27   AJMIN  IN     Jacobian approximation without feedback:
C                        Threshold value (Default: 0.0d0)
C                          The absolute pertubation for component k is
C                          computed by 
C                          DELX := AJDEL*max(abs(Xk),AJMIN)
C     28  ETADIF  IN     Jacobian approximation with feedback:
C                        Target value for relative pertubation ETA of X
C                        (Default: 1.0d-6)
C     29  ETAINI  IN     Jacobian approximation with feedback:
C                        Initial value for denominator differences
C                        (Default: 1.0d-6)
C     30..50             Reserved
C
C     Note 6:
C       The default values of the internal parameters may be obtained
C       from the monitor output with at least IOPT field MPRMON set to 2
C       and by initializing the corresponding RWK-fields to zero. 
C
C*   Error and warning messages:
C    ===========================
C
C      1    Termination, since jacobian matrix became singular
C      2    Termination after NITMAX iterations ( as indicated by
C           input parameter NITMAX=IWK(31) )
C      3    Termination, since damping factor became to small
C      4    Warning: Superlinear or quadratic convergence slowed down
C           near the solution.
C           Iteration has been stopped therefore with an approximation
C           of the solution not such accurate as requested by RTOL,
C           because possibly the RTOL requirement may be too stringent
C           (i.e. the nonlinear problem is ill-conditioned)
C      5    Warning: Iteration stopped with termination criterion 
C           (using RTOL as requested precision) satisfied, but no 
C           superlinear or quadratic convergence has been indicated yet.
C           Therefore, possibly the error estimate for the solution may
C           not match good enough the really achieved accuracy.
C     10    Integer or real workspace too small
C     20    Bad input to dimensional parameter N
C     21    Nonpositive value for RTOL supplied
C     22    Negative scaling value via vector XSCAL supplied
C     30    One or more fields specified in IOPT are invalid
C           (for more information, see error-printout)
C     80    Error signalled by linear solver routine N1FACT,
C           for more detailed information see IFAIL-value
C           stored to IWK(23)
C     81    Error signalled by linear solver routine N1SOLV,
C           for more detailed information see IFAIL-value
C           stored to IWK(23)
C           (not used by standard routine N1SOLV)
C     82    Error signalled by user routine FCN (Nonzero value
C           returned via IFAIL-flag; stored to IWK(23) )
C     83    Error signalled by user routine JAC (Nonzero value
C           returned via IFAIL-flag; stored to IWK(23) )
C
C     Note 7 : in case of failure:
C        -    use non-standard options
C        -    or turn to Newton-algorithm with rank strategy
C        -    use another initial guess
C        -    or reformulate model
C        -    or apply continuation techniques
C
C*    Machine dependent constants used:
C     =================================
C
C     DOUBLE PRECISION EPMACH  in  N1PCHK, N1INT
C     DOUBLE PRECISION GREAT   in  N1PCHK
C     DOUBLE PRECISION SMALL   in  N1PCHK, N1INT, N1SCAL
C
C*    Subroutines called: N1PCHK, N1INT
C
C     ------------------------------------------------------------
C*    End Prologue
C
C*    Summary of changes:
C     ===================
C      
C     2.2.1  91, May 30     Time monitor included
C     2.2.2  91, May 30     Bounded damping strategy implemented
C     2.2.3  91, June 19    AJDEL, AJMIN as RWK-options for JACGEN.EQ.2,
C                           ETADIF, ETAINI as RWK-opts. for JACGEN.EQ.3
C                           FCN-count changed for anal. Jacobian
C     2.2.4  91, August  9  Convergence order monitor included
C     2.2.5  91, August 13  Standard Broyden updates replaced by
C                           iterative Broyden
C     2.2.6  91, Sept.  16  Damping factor reduction by FCN-fail imple-
C                           mented
C     2.3    91, Dec.   20  New Release for CodeLib
C            92, March  11  Level of accepted simplified correction
C                           stored to RWK(IRWKI+4)
C            00, July   12  RTOL output-value bug fixed
C            06, Jan.   24  IERR=5 no longer returned if residuum of
C                           final iterate is exactly zero
C     2.4    09, May    29  Changed routines N1FACT and N1SOLV to use
C                           LAPACK Routines DGETRF, DGETRS, DGBTRF and
C                           DGBTRS instead of LINPACK routines to solve 
C                           linear systems.
C            10, July   26  Subroutine N1INT: Initialization of unitialized
C                           Variable FCMON fixed.
C   
C     ------------------------------------------------------------
C
C     PARAMETER (IRWKI=xx, LRWKI=yy)  
C     IRWKI: Start position of internally used RWK part
C     LRWKI: Length of internally used RWK part
C     (current values see parameter statement below)
C
C     INTEGER L4,L5,L51,L6,L61,L62,L63,L7,L71,L8,L9,L10,L11,L12,L121,
C             L13,L14,L20
C     Starting positions in RWK of formal array parameters of internal
C     routine N1INT (dynamically determined in driver routine NLEQ1,
C     dependent on N and options setting)
C
C     Further RWK positions (only internally used)
C
C     Position  Name     Meaning
C
C     IRWKI     FCKEEP   Damping factor of previous successfull iter.
C     IRWKI+1   FCA      Previous damping factor
C     IRWKI+2   FCPRI    A priori estimate of damping factor
C     IRWKI+3   DMYCOR   Number My of latest corrector damping factor
C                        (kept for use in rank-1 decision criterium)
C     IRWKI+4   SUMXS    natural level of accepted simplified correction
C     IRWKI+(5..LRWKI-1) Free
C
C     Internal arrays stored in RWK (see routine N1INT for descriptions)
C
C     Position  Array         Type   Remarks
C
C     L4        A(M1,N)       Perm   M1=N (full mode) or 
C                                    M1=2*IOPT(6)+IOPT(7)+1 (band mode)
C     L41       DXSAVE(N,NBROY)
C                             Perm   NBROY=IWK(36) (Default: N or 0)
C     L5        DX(N)         Perm  
C     L51       DXQ(N)        Perm 
C     L6        XA(N)         Perm
C     L61       F(N)          Perm
C     L62       FW(N)         Perm
C     L63       XWA(N)        Perm
C     L7        FA(N)         Perm
C     L71       ETA(N)        Perm   Only used for JACGEN=IOPT(3)=3
C     L9        XW(N)         Temp
C     L11       DXQA(N)       Temp
C     L12       T1(N)         Temp
C     L121      T2(N)         Temp
C     L13       T3(N)         Temp
C     L14                     Temp   Start position of array workspace 
C                                    needed for linear solver  
C     
C
      EXTERNAL N1INT
      INTRINSIC DBLE
      INTEGER IRWKI, LRWKI
      PARAMETER (IRWKI=51, LRWKI=10)  
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=1.0D1)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER NITMAX,LUERR,LUMON,LUSOL,MPRERR,MPRMON,MPRSOL,
     $MSTOR,M1,M2,NRWKFR,NRFRIN,NRW,NIWKFR,NIFRIN,NIW,NONLIN,JACGEN
      INTEGER L4,L41,L5,L51,L6,L61,L62,L63,L7,L71,L8,L9,L11,L12,L121,
     $L13,L14,L20
      DOUBLE PRECISION FC,FCMIN,PERCI,PERCR
      LOGICAL QINIMO,QRANK1,QFCSTR,QSUCC,QBDAMP,QSIMPL
      CHARACTER CHGDAT*20, PRODCT*8
C     Which version ?
      LOGICAL QVCHK
      INTEGER IVER
      PARAMETER( IVER=21112401 )
C
C     Version: 2.4.0.1           Latest change:
C     -----------------------------------------
C
      DATA      CHGDAT      /'July 26, 2010       '/
      DATA      PRODCT      /'NLEQ1   '/
C*    Begin
      IERR = 0
      QVCHK = IWK(12).LT.0
      IWK(12) = IVER
      IF (QVCHK) RETURN
C        Print error messages?
      MPRERR = IOPT(11)
      LUERR = IOPT(12)
      IF (LUERR .EQ. 0) THEN
        LUERR = 6
        IOPT(12)=LUERR
      ENDIF
C        Print iteration monitor?
      MPRMON = IOPT(13)
      LUMON = IOPT(14)
      IF (LUMON .LE. 0 .OR. LUMON .GT. 99) THEN
        LUMON = 6
        IOPT(14)=LUMON
      ENDIF
C        Print intermediate solutions?
      MPRSOL = IOPT(15)
      LUSOL = IOPT(16)
      IF (LUSOL .EQ. 0) THEN
        LUSOL = 6
        IOPT(16)=LUSOL
      ENDIF
C        Print time summary statistics?
      MPRTIM = IOPT(19)
      LUTIM = IOPT(20)
      IF (LUTIM .EQ. 0) THEN
        LUTIM = 6
        IOPT(20)=LUTIM
      ENDIF
      QSUCC = IOPT(1).EQ.1
      QINIMO = MPRMON.GE.1.AND..NOT.QSUCC
C     Print NLEQ1 heading lines
      IF(QINIMO)THEN
10000   FORMAT('   N L E Q 1  *****  V e r s i o n  ',
     $         '2 . 3 ***',//,1X,'Newton-Method ',
     $         'for the solution of nonlinear systems',//)
        WRITE(LUMON,10000)
      ENDIF
C     Check input parameters and options
      CALL N1PCHK(N,X,XSCAL,RTOL,IOPT,IERR,LIWK,IWK,LRWK,RWK)
C     Exit, if any parameter error was detected till here
      IF (IERR.NE.0) RETURN 
C
      MSTOR=IOPT(4)
      IF (MSTOR.EQ.0) THEN
        M1=N
        M2=N
      ELSE IF (MSTOR.EQ.1) THEN
        ML=IOPT(6)
        MU=IOPT(7)
        M1=2*ML+MU+1
        M2=ML+MU+1
      ENDIF
      JACGEN=IOPT(3)
      IF (JACGEN.EQ.0) JACGEN=2
      IOPT(3)=JACGEN
      QRANK1=IOPT(32).EQ.1
      QSIMPL=IOPT(34).EQ.1
      IF (QRANK1) THEN
        NBROY=IWK(36)
        IF (NBROY.EQ.0) NBROY=MAX(M2,10)
        IWK(36)=NBROY
      ELSE
        NBROY=0
      ENDIF
C     WorkSpace: RWK
      L4=IRWKI+LRWKI
      L41=L4+M1*N
      L5=L41+NBROY*N
      L51=L5+N
      L6=L51+N
      L61=L6+N
      L62=L61+N
      L63=L62+N
      L7=L63+N
      L71=L7+N
      IF (JACGEN.NE.3) THEN
        L8=L71
      ELSE
        L8=L71+N
      ENDIF
      NRWKFR = L8
      L9=L8
      L11=L9+N
      L12=L11+N
      L121=L12+N
      L13=L121+N
      L14=L13+N
      NRW=L14-1
C     End WorkSpace at NRW
C     WorkSpace: IWK
      L20=51
      NIWKFR = L20
      IF (QRANK1.OR.QSIMPL) NIWKFR = NIWKFR+N
      NIW=L20-1
C     End WorkSpace at NIW
      IWK(16) = NIW+1
      IWK(17) = NRW+1
      NIFRIN = NIW+1
      NRFRIN = NRW+1
C
      IF(NRW.GT.LRWK.OR.NIW.GT.LIWK)THEN
        IERR=10
      ELSE
        IF(QINIMO)THEN
          PERCR = DBLE(NRW)/DBLE(LRWK)*100.0D0
          PERCI = DBLE(NIW)/DBLE(LIWK)*100.0D0
C         Print statistics concerning workspace usage
10050     FORMAT(' Real    Workspace declared as ',I9,
     $    ' is used up to ',I9,' (',F5.1,' percent)',//,
     $    ' Integer Workspace declared as ',I9,
     $    ' is used up to ',I9,' (',F5.1,' percent)',//)
          WRITE(LUMON,10050)LRWK,NRW,PERCR,LIWK,NIW,PERCI
        ENDIF
        IF(QINIMO)THEN
10051     FORMAT(/,' N =',I4,//,' Prescribed relative ',
     $    'precision',D10.2,/)
          WRITE(LUMON,10051)N,RTOL
10052     FORMAT(' The Jacobian is supplied by ',A)
          IF (JACGEN.EQ.1) THEN
            WRITE(LUMON,10052) 'a user subroutine'
          ELSE IF (JACGEN.EQ.2) THEN
             WRITE(LUMON,10052) 
     $        'numerical differentiation (without feedback strategy)'
          ELSE IF (JACGEN.EQ.3) THEN
             WRITE(LUMON,10052) 
     $        'numerical differentiation (feedback strategy included)'
          ENDIF
10055     FORMAT(' The Jacobian will be stored in ',A,' mode')
          IF (MSTOR.EQ.0) THEN
            WRITE(LUMON,10055) 'full'
          ELSE IF (MSTOR.EQ.1) THEN
            WRITE(LUMON,10055) 'banded'
10056       FORMAT(' Lower bandwidth : ',I3,'   Upper bandwidth : ',I3)
            WRITE(LUMON,10056) ML,MU
          ENDIF
10057     FORMAT(' Automatic row scaling of the Jacobian is ',A,/)
          IF (IOPT(35).EQ.1) THEN
            WRITE(LUMON,10057) 'inhibited'
          ELSE
            WRITE(LUMON,10057) 'allowed'
          ENDIF
        ENDIF
        NONLIN=IOPT(31)
        IF (IOPT(38).EQ.0) QBDAMP = NONLIN.EQ.4
        IF (IOPT(38).EQ.1) QBDAMP = .TRUE.
        IF (IOPT(38).EQ.2) QBDAMP = .FALSE.
        IF (QBDAMP) THEN
          IF (RWK(20).LT.ONE) RWK(20) = TEN
        ENDIF
10064   FORMAT(' Rank-1 updates are ',A)
        IF (QINIMO) THEN
          IF (QRANK1) THEN
            WRITE(LUMON,10064) 'allowed'
          ELSE
            WRITE(LUMON,10064) 'inhibited'
          ENDIF
10065     FORMAT(' Problem is specified as being ',A)
          IF (NONLIN.EQ.1) THEN
            WRITE(LUMON,10065) 'linear'
          ELSE IF (NONLIN.EQ.2) THEN
            WRITE(LUMON,10065) 'mildly nonlinear'
          ELSE IF (NONLIN.EQ.3) THEN
            WRITE(LUMON,10065) 'highly nonlinear'
          ELSE IF (NONLIN.EQ.4) THEN
            WRITE(LUMON,10065) 'extremely nonlinear'
          ENDIF
10066     FORMAT(' Bounded damping strategy is ',A,:,/, 
     $           ' Bounding factor is ',D10.3)
          IF (QBDAMP) THEN
            WRITE(LUMON,10066) 'active', RWK(20)
          ELSE
            WRITE(LUMON,10066) 'off'
          ENDIF
10067     FORMAT(' Special mode: ',A,' Newton iteration will be done')
          IF (IOPT(33).EQ.1) WRITE(LUMON,10067) 'Ordinary'
          IF (IOPT(34).EQ.1) WRITE(LUMON,10067) 'Simplified'
        ENDIF
C       Maximum permitted number of iteration steps
        NITMAX=IWK(31)
        IF (NITMAX.LE.0) NITMAX=50
        IWK(31)=NITMAX
10068   FORMAT(' Maximum permitted number of iteration steps : ',
     $         I6)
        IF (QINIMO) WRITE(LUMON,10068) NITMAX
C       Initial damping factor for highly nonlinear problems
        QFCSTR=RWK(21).GT.ZERO
        IF (.NOT.QFCSTR) THEN
          RWK(21)=1.0D-2
          IF (NONLIN.EQ.4) RWK(21)=1.0D-4
        ENDIF
C       Minimal permitted damping factor
        IF (RWK(22).LE.ZERO) THEN
          RWK(22)=1.0D-4
          IF (NONLIN.EQ.4) RWK(22)=1.0D-8
        ENDIF
        FCMIN=RWK(22)
C       Rank1 decision parameter SIGMA
        IF (RWK(23).LT.ONE) RWK(23)=3.0D0
        IF (.NOT.QRANK1) RWK(23)=10.0D0/FCMIN
C       Decision parameter about increasing too small predictor
C       to greater corrector value
        IF (RWK(24).LT.ONE) RWK(24)=10.0D0/FCMIN       
C       Starting value of damping factor (FCMIN.LE.FC.LE.1.0)
        IF(NONLIN.LE.2.AND..NOT.QFCSTR)THEN
C         for linear or mildly nonlinear problems
          FC = ONE
        ELSE
C         for highly or extremely nonlinear problems
          FC = RWK(21)
        ENDIF
C       Simplied Newton iteration implies ordinary Newton it. mode
        IF (IOPT(34).EQ.1) IOPT(33)=1
C       If ordinary Newton iteration, factor is always one
        IF (IOPT(33).EQ.1) FC=1.0D0
        RWK(21)=FC
        IF (MPRMON.GE.2.AND..NOT.QSUCC) THEN
10069     FORMAT(//,' Internal parameters:',//,
     $      ' Starting value for damping factor FCSTART = ',D9.2,/,
     $      ' Minimum allowed damping factor FCMIN = ',D9.2,/,
     $      ' Rank-1 updates decision parameter SIGMA = ',D9.2)
          WRITE(LUMON,10069) RWK(21),FCMIN,RWK(23)
        ENDIF
C       Store lengths of currently required workspaces
        IWK(18) = NIFRIN-1
        IWK(19) = NRFRIN-1
C
C       Initialize and start time measurements monitor
C
        IF ( IOPT(1).EQ.0 .AND. MPRTIM.NE.0 ) THEN
          CALL MONINI (' NLEQ1',LUTIM)
          CALL MONDEF (0,'NLEQ1')
          CALL MONDEF (1,'FCN')
          CALL MONDEF (2,'Jacobi')
          CALL MONDEF (3,'Lin-Fact')
          CALL MONDEF (4,'Lin-Sol')
          CALL MONDEF (5,'Output')
          CALL MONSTR (IERR)
        ENDIF
C
C
C
        IERR=-1
C       If IERR is unmodified on exit, successive steps are required
C       to complete the Newton iteration
        IF (NBROY.EQ.0) NBROY=1
        CALL N1INT(N,FCN,JAC,X,XSCAL,RTOL,NITMAX,NONLIN,IOPT,IERR,
     $  LRWK,RWK,NRFRIN,LIWK,IWK,NIFRIN,
     $  M1,M2,NBROY,
     $  RWK(L4),RWK(L41),RWK(L5),RWK(L51),RWK(L6),RWK(L63),RWK(L61),
     $  RWK(L7),
     $  RWK(L71),RWK(L9),RWK(L62),RWK(L11),RWK(L12),RWK(L121),RWK(L13),
     $  RWK(21),RWK(22),RWK(23),RWK(24),RWK(IRWKI+1),RWK(IRWKI),
     $  RWK(IRWKI+2),RWK(IRWKI+3),RWK(17),RWK(18),RWK(IRWKI+4),RWK(19),
     $  MSTOR,MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,IWK(1),IWK(3),
     $  IWK(4),IWK(5),IWK(8),IWK(9),IWK(33),IWK(24),QBDAMP)
C
        IF (MPRTIM.NE.0.AND.IERR.NE.-1.AND.IERR.NE.10) THEN
          CALL MONHLT
          CALL MONPRT
        ENDIF
C
C       Free workspaces, so far not used between steps
        IWK(16) = NIWKFR
        IWK(17) = NRWKFR
      ENDIF
C     Print statistics
      IF (MPRMON.GE.1.AND.IERR.NE.-1.AND.IERR.NE.10) THEN
10080   FORMAT(/, '   ******  Statistics * ', A8, ' *******', /,
     $            '   ***  Newton iterations : ', I7,'  ***', /,
     $            '   ***  Corrector steps   : ', I7,'  ***', /,
     $            '   ***  Rejected rk-1 st. : ', I7,'  ***', /,
     $            '   ***  Jacobian eval.    : ', I7,'  ***', /,
     $            '   ***  Function eval.    : ', I7,'  ***', /,
     $            '   ***  ...  for Jacobian : ', I7,'  ***', /,
     $            '   *************************************', /)
        WRITE (LUMON,10080) PRODCT,IWK(1),IWK(3),IWK(9),IWK(5),
     $  IWK(4),IWK(8)
      ENDIF
C     Print workspace requirements, if insufficient
      IF (IERR.EQ.10) THEN
10090   FORMAT(//,' ',20('*'),'Workspace Error',20('*'))
        IF (MPRERR.GE.1) WRITE(LUERR,10090)
        IF(NRW.GT.LRWK)THEN
10091     FORMAT(/,' Real Workspace dimensioned as',1X,I9,
     $    1X,'must be enlarged at least up to ',
     $    I9,//)
          IF (MPRERR.GE.1) WRITE(LUERR,10091)LRWK,NRFRIN-1
        ENDIF
        IF(NIW.GT.LIWK)THEN
10092     FORMAT(/,' Integer Workspace dimensioned as ',
     $    I9,' must be enlarged at least up ',
     $    'to ',I9,//)
          IF (MPRERR.GE.1) WRITE(LUERR,10092)LIWK,NIFRIN-1
        ENDIF
      ENDIF
C     End of subroutine NLEQ1
      RETURN
      END
C
      SUBROUTINE N1PCHK(N,X,XSCAL,RTOL,IOPT,IERR,LIWK,IWK,LRWK,RWK)
C*    Begin Prologue N1PCHK
      INTEGER N
      DOUBLE PRECISION X(N),XSCAL(N)
      DOUBLE PRECISION RTOL
      INTEGER IOPT(50)
      INTEGER IERR
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LRWK
      DOUBLE PRECISION RWK(LRWK)
C     ------------------------------------------------------------
C
C*    Summary :
C
C     N 1 P C H K : Checking of input parameters and options
C                   for NLEQ1.
C
C*    Parameters:
C     ===========
C
C     See parameter description in driver routine.
C
C*    Subroutines called: ZIBCONST
C
C*    Machine dependent constants used:
C     =================================
C
C     EPMACH = relative machine precision
C     GREAT = squareroot of maxreal divided by 10
C     SMALL = squareroot of "smallest positive machine number
C             divided by relative machine precision"
      DOUBLE PRECISION EPMACH,GREAT,SMALL
C
C     ------------------------------------------------------------
C*    End Prologue
C
      INTRINSIC DBLE
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=1.0D1)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
      PARAMETER (NUMOPT=50)
      INTEGER MSTOR,IOPTL(NUMOPT),IOPTU(NUMOPT)
      DOUBLE PRECISION TOLMIN,TOLMAX,DEFSCL
C
      DATA IOPTL /0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,1,
     $            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     $            0,0,0,0,0,0,0,0,0,0,
     $            -9999999,-9999999,-9999999,-9999999,-9999999/
      DATA IOPTU /1,1,3,1,0,9999999,9999999,0,1,0,3,99,6,99,3,99,0,0,1,
     $            99,0,0,0,0,0,0,0,0,0,0,4,1,1,1,1,
     $            0,0,2,3,0,0,0,0,0,0,
     $            9999999,9999999,9999999,9999999,9999999/
C
      CALL ZIBCONST(EPMACH,SMALL)
      GREAT  = 1.0D0/SMALL
      IERR = 0
C        Print error messages?
      MPRERR = IOPT(11)
      LUERR = IOPT(12)
      IF (LUERR .LE. 0 .OR. LUERR .GT. 99) THEN
        LUERR = 6
        IOPT(12)=LUERR
      ENDIF
C
C     Checking dimensional parameter N
      IF ( N.LE.0 ) THEN
        IF (MPRERR.GE.1)  WRITE(LUERR,10011) N
10011   FORMAT(/,' Error: Bad input to dimensional parameter N supplied'
     $         ,/,8X,'choose N positive, your input is: N = ',I5)
        IERR = 20
      ENDIF
C
C     Problem type specification by user
      NONLIN=IOPT(31)
      IF (NONLIN.EQ.0) NONLIN=3
      IOPT(31)=NONLIN
C
C     Checking and conditional adaption of the user-prescribed RTOL
      IF (RTOL.LE.ZERO) THEN
        IF (MPRERR.GE.1) 
     $      WRITE(LUERR,'(/,A)') ' Error: Nonpositive RTOL supplied'
        IERR = 21
      ELSE
        TOLMIN = EPMACH*TEN*DBLE(N)
        IF(RTOL.LT.TOLMIN) THEN
          RTOL = TOLMIN
          IF (MPRERR.GE.2) 
     $      WRITE(LUERR,10012) 'increased ','smallest',RTOL
        ENDIF
        TOLMAX = 1.0D-1
        IF(RTOL.GT.TOLMAX) THEN
          RTOL = TOLMAX
          IF (MPRERR.GE.2) 
     $      WRITE(LUERR,10012) 'decreased ','largest',RTOL
        ENDIF
10012   FORMAT(/,' Warning: User prescribed RTOL ',A,'to ',
     $         'reasonable ',A,' value RTOL = ',D11.2)
      ENDIF
C     
C     Test user prescribed accuracy and scaling on proper values
      IF (N.LE.0) RETURN 
      IF (NONLIN.GE.3) THEN
        DEFSCL = RTOL
      ELSE
        DEFSCL = ONE
      ENDIF
      DO 10 I=1,N
        IF (XSCAL(I).LT.ZERO) THEN
          IF (MPRERR.GE.1) THEN 
            WRITE(LUERR,10013) I
10013       FORMAT(/,' Error: Negative value in XSCAL(',I5,') supplied')
          ENDIF
          IERR = 22
        ENDIF
        IF (XSCAL(I).EQ.ZERO) XSCAL(I) = DEFSCL
        IF ( XSCAL(I).GT.ZERO .AND. XSCAL(I).LT.SMALL ) THEN
          IF (MPRERR.GE.2) THEN
            WRITE(LUERR,10014) I,XSCAL(I),SMALL
10014       FORMAT(/,' Warning: XSCAL(',I5,') = ',D9.2,' too small, ',
     $             'increased to',D9.2)
          ENDIF
          XSCAL(I) = SMALL
        ENDIF
        IF (XSCAL(I).GT.GREAT) THEN
          IF (MPRERR.GE.2) THEN
            WRITE(LUERR,10015) I,XSCAL(I),GREAT
10015       FORMAT(/,' Warning: XSCAL(',I5,') = ',D9.2,' too big, ',
     $             'decreased to',D9.2)
          ENDIF
          XSCAL(I) = GREAT
        ENDIF
10    CONTINUE
C     Special dependence on Jacobian's storage mode for ML and MU
      MSTOR = IOPT(4)
      IF (MSTOR.EQ.0) THEN
        IOPTU(6)=0
        IOPTU(7)=0
      ELSE IF (MSTOR.EQ.1) THEN
        IOPTU(6)=N-1
        IOPTU(7)=N-1
      ENDIF
C     Checks options
      DO 20 I=1,30
        IF (IOPT(I).LT.IOPTL(I) .OR. IOPT(I).GT.IOPTU(I)) THEN
          IERR=30
          IF (MPRERR.GE.1) THEN
            WRITE(LUERR,20001) I,IOPT(I),IOPTL(I),IOPTU(I)
20001       FORMAT(' Invalid option specified: IOPT(',I2,')=',I12,';',
     $             /,3X,'range of permitted values is ',I8,' to ',I8)
          ENDIF
        ENDIF
20    CONTINUE
C     End of subroutine N1PCHK
      RETURN
      END
C
      SUBROUTINE N1INT(N,FCN,JAC,X,XSCAL,RTOL,NITMAX,NONLIN,IOPT,IERR,
     $LRWK,RWK,NRWKFR,LIWK,IWK,NIWKFR,M1,M2,NBROY,
     $A,DXSAVE,DX,DXQ,XA,XWA,F,FA,ETA,XW,FW,DXQA,T1,T2,T3,FC,FCMIN,
     $SIGMA,SIGMA2,FCA,FCKEEP,FCPRI,DMYCOR,CONV,SUMX,SUMXS,DLEVF,MSTOR,
     $MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,NITER,NCORR,NFCN,NJAC,
     $NFCNJ,NREJR1,NEW,ICONV,QBDAMP)
C*    Begin Prologue N1INT
      INTEGER N
      EXTERNAL FCN,JAC
      DOUBLE PRECISION X(N),XSCAL(N)
      DOUBLE PRECISION RTOL
      INTEGER NITMAX,NONLIN
      INTEGER IOPT(50)
      INTEGER IERR
      INTEGER LRWK
      DOUBLE PRECISION RWK(LRWK)
      INTEGER NRWKFR,LIWK
      INTEGER IWK(LIWK)
      INTEGER NIWKFR,M1,M2,NBROY
      DOUBLE PRECISION A(M1,N),DXSAVE(N,NBROY)
      DOUBLE PRECISION DX(N),DXQ(N),XA(N),XWA(N),F(N),FA(N),ETA(N)
      DOUBLE PRECISION XW(N),FW(N),DXQA(N),T1(N),T2(N),T3(N)
      DOUBLE PRECISION FC,FCMIN,SIGMA,SIGMA2,FCA,FCKEEP,CONV,SUMX,SUMXS,
     $                 DLEVF,FCPRI,DMYCOR
      INTEGER MSTOR,MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,NITER,
     $NCORR,NFCN,NJAC,NFCNJ,NREJR1,NEW,ICONV
      LOGICAL QBDAMP
C     ------------------------------------------------------------
C
C*    Summary :
C
C     N 1 I N T : Core routine for NLEQ1 .
C     Damped Newton-algorithm for systems of highly nonlinear
C     equations especially designed for numerically sensitive
C     problems.
C
C*    Parameters:
C     ===========
C
C       N,FCN,JAC,X,XSCAL,RTOL   
C                         See parameter description in driver routine
C
C       NITMAX      Int    Maximum number of allowed iterations
C       NONLIN      Int    Problem type specification
C                          (see IOPT-field NONLIN)
C       IOPT        Int    See parameter description in driver routine
C       IERR        Int    See parameter description in driver routine
C       LRWK        Int    Length of real workspace
C       RWK(LRWK)   Dble   Real workspace array
C       NRWKFR      Int    First free position of RWK on exit 
C       LIWK        Int    Length of integer workspace
C       IWK(LIWK)   Int    Integer workspace array
C       NIWKFR      Int    First free position of IWK on exit 
C       M1          Int    Leading dimension of Jacobian array A
C                          for full case Jacobian: N
C                          for banded Jacobian: 2*ML+MU+1;
C                          ML, MU see IOPT-description in driver routine
C       M2          Int    for full case Jacobian: N
C                          for banded Jacobian: ML+MU+1
C       NBROY       Int    Maximum number of possible consecutive
C                          iterative Broyden steps. (See IWK(36))
C       A(M1,N)     Dble   Holds the Jacobian matrix (decomposed form
C                          after call of linear decomposition
C                          routine)
C       DXSAVE(X,NBROY)
C                   Dble   Used to save the quasi Newton corrections of
C                          all previously done consecutive Broyden
C                          steps.
C       DX(N)       Dble   Current Newton correction
C       DXQ(N)      Dble   Simplified Newton correction J(k-1)*X(k)
C       XA(N)       Dble   Previous Newton iterate
C       XWA(N)      Dble   Scaling factors used for latest decomposed
C                          Jacobian for column scaling - may differ
C                          from XW, if Broyden updates are performed
C       F(N)        Dble   Function (FCN) value of current iterate
C       FA(N)       Dble   Function (FCN) value of previous iterate
C       ETA(N)      Dble   Jacobian approximation: updated scaled
C                          denominators
C       XW(N)       Dble   Scaling factors for iteration vector
C       FW(N)       Dble   Scaling factors for rows of the system
C       DXQA(N)     Dble   Previous Newton correction
C       T1(N)       Dble   Workspace for linear solvers and internal
C                          subroutines
C       T2(N)       Dble   Workspace array for internal subroutines
C       T3(N)       Dble   Workspace array for internal subroutines
C       FC          Dble   Current Newton iteration damping factor.
C       FCMIN       Dble   Minimum permitted damping factor. If
C                          FC becomes smaller than this value, one
C                          of the following may occur:
C                          a.    Recomputation of the Jacobian
C                                matrix by means of difference
C                                approximation (instead of Rank1
C                                update), if Rank1 - update
C                                previously was used
C                          b.    Fail exit otherwise
C       SIGMA       Dble   Decision parameter for rank1-updates
C       SIGMA2      Dble   Decision parameter for damping factor
C                          increasing to corrector value
C       FCA         Dble   Previous Newton iteration damping factor.
C       FCKEEP      Dble   Keeps the damping factor as it is at start
C                          of iteration step.
C       CONV        Dble   Scaled maximum norm of the Newton-
C                          correction. Passed to RWK-field on output.
C       SUMX        Dble   Square of the natural level (see equal-
C                          named RWK-output field)
C       SUMXS       Dble   Square of the "simplified" natural level
C                          (see equal-named RWK-internal field)
C       DLEVF       Dble   The standard level (see equal-
C                          named RWK-output field)
C       MSTOR       Int    see description of IOPT-field MSTOR
C       MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,
C       NITER,NCORR,NFCN,NJAC,NFCNJ,NREJR1,NEW :
C                          See description of equal named IWK-fields
C                          in the driver subroutine
C       QBDAMP      Logic  Flag, that indicates, whether bounded damping
C                          strategy is active:
C                          .true.  = bounded damping strategy is active
C                          .false. = normal damping strategy is active
C
C*    Internal double variables
C     =========================
C
C       AJDEL    See RWK(26) (num. diff. without feedback)
C       AJMIN    See RWK(27) (num. diff. without feedback)
C       CONVA    Holds the previous value of CONV .
C       DMUE     Temporary value used during computation of damping 
C                factors predictor.
C       EPDIFF   sqrt(10*epmach) (num. diff. with feedback)
C       ETADIF   See description of RWK(28) (num. diff. with feedback)
C       ETAINI   Initial value for all ETA-components (num. diff. fb.)
C       ETAMAX   Maximum allowed pertubation (num. diff. with feedback)
C       ETAMIN   Minimum allowed pertubation (num. diff. with feedback)
C       FCDNM    Used to compute the denominator of the damping 
C                factor FC during computation of it's predictor,
C                corrector and aposteriori estimate (in the case of
C                performing a Rank1 update) .
C       FCK2     Aposteriori estimate of FC.
C       FCH      Temporary storage of new corrector FC.
C       FCMIN2   FCMIN**2 . Used for FC-predictor computation.
C       FCNUMP   Gets the numerator of the predictor formula for FC.
C       FCNMP2   Temporary used for predictor numerator computation.
C       FCNUMK   Gets the numerator of the corrector computation 
C                of FC .
C       SUMXA    Natural level of the previous iterate.
C       TH       Temporary variable used during corrector- and 
C                aposteriori computations of FC.
C
C*    Internal integer variables
C     ==========================
C
C     IFAIL      Gets the return value from subroutines called from
C                N1INT (N1FACT, N1SOLV, FCN, JAC) 
C     ISCAL      Holds the scaling option from the IOPT-field ISCAL      
C     MODE       Matrix storage mode (see IOPT-field MODE) 
C     NRED       Count of successive corrector steps
C     NILUSE     Gets the amount of IWK used by the linear solver
C     NRLUSE     Gets the amount of RWK used by the linear solver
C     NIWLA      Index of first element of IWK provided to the
C                linear solver
C     NRWLA      Index of first element of RWK provided to the
C                linear solver
C     LIWL       Holds the maximum amount of integer workspace
C                available to the linear solver
C     LRWL       Holds the maximum amount of real workspace
C                available to the linear solver
C
C
C*    Internal logical variables
C     ==========================
C
C     QGENJ      Jacobian updating technique flag:
C                =.TRUE.  : Call of analytical subroutine JAC or
C                           numerical differentiation
C                =.FALSE. : rank1- (Broyden-) update
C     QINISC     Iterate initial-scaling flag:
C                =.TRUE.  : at first call of N1SCAL
C                =.FALSE. : at successive calls of N1SCAL
C     QSUCC      See description of IOPT-field QSUCC.
C     QJCRFR     Jacobian refresh flag:
C                set to .TRUE. if damping factor gets too small
C                and Jacobian was computed by rank1-update. 
C                Indicates, that the Jacobian needs to be recomputed
C                by subroutine JAC or numerical differentation.
C     QLINIT     Initialization state of linear solver workspace:
C                =.FALSE. : Not yet initialized
C                =.TRUE.  : Initialized - N1FACT has been called at
C                           least one time.
C     QSCALE     Holds the value of .NOT.QNSCAL. See description
C                of IOPT-field QNSCAL.
C
C*    Subroutines called:
C     ===================
C
C       N1FACT, N1SOLV, N1JAC,  N1JACB, N1JCF,  N1JCFB, N1LVLS, 
C       N1SCRF, N1SCRB, N1SOUT, N1PRV1, N1PRV2, N1SCAL,
C       MONON,  MONOFF
C
C*    Functions called:
C     =================
C
C       ZIBCONST, WNORM
C
C*    Machine constants used
C     ======================
C
      DOUBLE PRECISION EPMACH,SMALL
C 
C     ------------------------------------------------------------
C*    End Prologue
      EXTERNAL N1FACT, N1SOLV, N1JAC,  N1JACB, N1JCF,  N1JCFB, N1LVLS,
     $         N1SCRF, N1SCRB, N1SOUT, N1PRV1, N1PRV2, N1SCAL,
     $         MONON,  MONOFF, WNORM
      INTRINSIC DSQRT,DMIN1,MAX0,MIN0
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=10.0D0)
      INTEGER IFAIL,ILOOP,ISCAL,K,MODE,NRED,NILUSE,NRLUSE,NIWLA,
     $NRWLA,LIWL,LRWL,L1,L2,JACGEN
      DOUBLE PRECISION AJDEL,AJMIN,ALFA1,ALFA2,ALFA,BETA,CONVA,
     $DLEVXA,DMYPRI,DXANRM,DXNRM,WNORM,EPDIFF,ETAMIN,ETAMAX,
     $ETAINI,ETADIF,FCDNM,FCMIN2,FCNUMP,FCCOR,FCNMP2,FCH,FCBND,
     $FCBH,FCK2,FCNUMK,FCREDU,DLEVFN,SUMXA,SUM1,SUM2,S1,TH,RSMALL,
     $APREC
      LOGICAL QGENJ,QINISC,QSUCC,QJCRFR,QLINIT,QSCALE,QNEXT,QREP,
     $        QRANK1,QMIXIO,QORDI,QSIMPL,QLU
CWEI
      INTRINSIC DLOG
      DOUBLE PRECISION CLIN0,CLIN1,CALPHA,CALPHK,ALPHAE,ALPHAK,ALPHAA,
     $                 SUMXA0,SUMXA1,SUMXA2,SUMXTE,FCMON,DLOG
      INTEGER IORMON
      LOGICAL QMSTOP
      SAVE CLIN0,CLIN1,CALPHA,ALPHAE,ALPHAK,ALPHAA,SUMXA0,SUMXA1,SUMXA2,
     $     QMSTOP,FCMON
C
      CALL ZIBCONST(EPMACH,SMALL)
C*    Begin
C       ----------------------------------------------------------
C       1 Initialization
C       ----------------------------------------------------------
C       1.1 Control-flags and -integers
        QSUCC = IOPT(1).EQ.1
        QSCALE = .NOT. IOPT(35).EQ.1
        QORDI  = IOPT(33).EQ.1
        QSIMPL = IOPT(34).EQ.1
        QRANK1 = IOPT(32).EQ.1
        IORMON = IOPT(39)
        IF (IORMON.EQ.0) IORMON=2
        ISCAL = IOPT(9)
        MODE = IOPT(2)
        JACGEN = IOPT(3)
        QMIXIO = LUMON.EQ.LUSOL .AND. MPRMON.NE.0 .AND. MPRSOL.NE.0
        QLU    = .NOT. QSIMPL
        MPRTIM = IOPT(19)
C       ----------------------------------------------------------
C       1.2 Derivated dimensional parameters
        IF (MSTOR.EQ.0) THEN
          ML=0
        ELSE IF (MSTOR.EQ.1) THEN
          ML=M1-M2
          MU=M2-1-ML
        ENDIF
C       ----------------------------------------------------------
C       1.3 Derivated internal parameters
        FCMIN2 = FCMIN*FCMIN
        RSMALL = DSQRT(TEN*RTOL)
C       ----------------------------------------------------------
C       1.4 Adaption of input parameters, if necessary
        IF(FC.LT.FCMIN) FC = FCMIN
        IF(FC.GT.ONE) FC = ONE
C       ----------------------------------------------------------
C       1.5 Initial preparations
        QJCRFR = .FALSE.
        QLINIT = .FALSE.
        IFAIL = 0
        FCBND = ZERO
        IF (QBDAMP) FCBND = RWK(20)
C       ----------------------------------------------------------
C       1.5.1 Numerical differentiation related initializations
        IF (JACGEN.EQ.2) THEN
          AJDEL = RWK(26)
          IF (AJDEL.LE.SMALL) AJDEL = DSQRT(EPMACH*TEN)
          AJMIN = RWK(27)
        ELSE IF (JACGEN.EQ.3) THEN
          ETADIF = RWK(28)
          IF (ETADIF .LE. SMALL) ETADIF = 1.0D-6
          ETAINI = RWK(29)
          IF (ETAINI .LE. SMALL) ETAINI = 1.0D-6
          EPDIFF = DSQRT(EPMACH*TEN)
          ETAMAX = DSQRT(EPDIFF)
          ETAMIN = EPDIFF*ETAMAX
        ENDIF
C       ----------------------------------------------------------
C       1.5.2 Miscellaneous preparations of first iteration step
        IF (.NOT.QSUCC) THEN
          NITER = 0
          NCORR = 0
          NREJR1 = 0
          NFCN = 0
          NJAC = 0
          NFCNJ = 0
          QGENJ = .TRUE.
          QINISC = .TRUE.
          FCKEEP = FC
          FCA = FC
          FCPRI = FC
          FCK2 = FC
          FCMON = FC
          CONV = ZERO
          IF (JACGEN.EQ.3) THEN
            DO 1520 L1=1,N
              ETA(L1)=ETAINI
1520        CONTINUE
          ENDIF
          DO 1521 L1=1,N
            XA(L1)=X(L1)
1521      CONTINUE
CWEI      
          ICONV = 0
          ALPHAE = ZERO
          SUMXA1 = ZERO
          SUMXA0 = ZERO
          CLIN0  = ZERO
          QMSTOP = .FALSE.
C         ------------------------------------------------------
C         1.6 Print monitor header
          IF(MPRMON.GE.2 .AND. .NOT.QMIXIO)THEN
16003       FORMAT(///,2X,66('*'))
            WRITE(LUMON,16003)
16004       FORMAT(/,8X,'It',7X,'Normf ',10X,'Normx ',8X,
     $             'Damp.Fct.',3X,'New')
            WRITE(LUMON,16004)
          ENDIF
C         --------------------------------------------------------
C         1.7 Startup step
C         --------------------------------------------------------
C         1.7.1 Computation of the residual vector
          IF (MPRTIM.NE.0) CALL MONON(1)
          CALL FCN(N,X,F,IFAIL)
          IF (MPRTIM.NE.0) CALL MONOFF(1)
          NFCN = NFCN+1
C     Exit, if ...
          IF (IFAIL.NE.0) THEN
            IERR = 82
            GOTO 4299
          ENDIF
        ELSE
          QINISC = .FALSE.
        ENDIF
C
C       Main iteration loop
C       ===================
C
C       Repeat
2       CONTINUE
C         --------------------------------------------------------
C         2 Startup of iteration step
          IF (.NOT.QJCRFR) THEN
C           ------------------------------------------------------
C           2.1 Scaling of variables X(N)
            CALL N1SCAL(N,X,XA,XSCAL,XW,ISCAL,QINISC,IOPT,LRWK,RWK)
            QINISC = .FALSE.
            IF(NITER.NE.0)THEN
C             ----------------------------------------------------
C             2.2 Aposteriori estimate of damping factor
              DO 2200 L1=1,N
                DXQA(L1)=DXQ(L1)
2200          CONTINUE
              IF (.NOT.QORDI) THEN
                FCNUMP = ZERO
                DO 2201 L1=1,N
                  FCNUMP=FCNUMP+(DX(L1)/XW(L1))**2
2201            CONTINUE
                TH = FC-ONE
                FCDNM = ZERO
                DO 2202 L1=1,N
                  FCDNM=FCDNM+((DXQA(L1)+TH*DX(L1))/XW(L1))**2
2202            CONTINUE
C               --------------------------------------------------
C               2.2.2 Decision criterion for Jacobian updating
C                     technique:
C                     QGENJ.EQ..TRUE. numerical differentation,
C                     QGENJ.EQ..FALSE. rank1 updating
                QGENJ = .TRUE.
                IF (FC.EQ.FCPRI) THEN
                  QGENJ = FC.LT.ONE.OR.FCA.LT.ONE.OR.DMYCOR.LE.FC*SIGMA
     $                    .OR. .NOT.QRANK1 .OR. NEW+2.GT.NBROY 
                  FCA = FC
                ELSE
                  DMYCOR = FCA*FCA*HALF*DSQRT(FCNUMP/FCDNM)
                  IF (NONLIN.LE.3) THEN
                    FCCOR = DMIN1(ONE,DMYCOR)
                  ELSE
                    FCCOR = DMIN1(ONE,HALF*DMYCOR)
                  ENDIF
                  FCA = DMAX1(DMIN1(FC,FCCOR),FCMIN)
C$Test-begin
                  IF (MPRMON.GE.5) THEN
                    WRITE(LUMON,22201) FCCOR, FC, DMYCOR, FCNUMP,
     $                                 FCDNM
22201               FORMAT (/, ' +++ aposteriori estimate +++', /,
     $                    ' FCCOR  = ', D18.10, '  FC     = ', D18.10, /,
     $                    ' DMYCOR = ', D18.10, '  FCNUMP = ', D18.10, /,
     $                    ' FCDNM  = ', D18.10, /,
     $                       ' ++++++++++++++++++++++++++++', /)
                  ENDIF
C$Test-end 
                ENDIF
                FCK2 = FCA
C               ------------------------------------------------------
C               2.2.1 Computation of the numerator of damping
C                     factor predictor
                FCNMP2 = ZERO
                DO 221 L1=1,N
                  FCNMP2=FCNMP2+(DXQA(L1)/XW(L1))**2
221             CONTINUE
                FCNUMP = FCNUMP*FCNMP2
              ENDIF
            ENDIF
          ENDIF
          QJCRFR =.FALSE.
C         --------------------------------------------------------
C         2.3 Jacobian matrix (stored to array A(M1,N))
C         --------------------------------------------------------
C         2.3.1 Jacobian generation by routine JAC or
C               difference approximation (If QGENJ.EQ..TRUE.)
C               - or -
C               Rank-1 update of Jacobian (If QGENJ.EQ..FALSE.)
          IF (QGENJ .AND. (.NOT.QSIMPL .OR. NITER.EQ.0)) THEN
            NEW = 0
            IF (JACGEN.EQ.1) THEN
               IF (MPRTIM.NE.0) CALL MONON(2)
               CALL JAC(N,M1,X,A,IFAIL)
               IF (MPRTIM.NE.0) CALL MONOFF(2)
            ELSE
              IF (MSTOR.EQ.0) THEN
                IF (MPRTIM.NE.0) CALL MONON(2)
                IF (JACGEN.EQ.3) 
     $            CALL N1JCF(FCN,N,N,X,F,A,XW,ETA,ETAMIN,ETAMAX,
     $                       ETADIF,CONV,NFCNJ,T1,IFAIL)
                IF (JACGEN.EQ.2) 
     $            CALL N1JAC(FCN, N, N, X, F, A, XW, AJDEL, AJMIN,
     $                       NFCNJ, T1, IFAIL)
                IF (MPRTIM.NE.0) CALL MONOFF(2)
              ELSE IF (MSTOR.EQ.1) THEN
                IF (MPRTIM.NE.0) CALL MONON(2)
                IF (JACGEN.EQ.3) 
     $            CALL N1JCFB(FCN,N,M1,ML,X,F,A,XW,ETA,ETAMIN,
     $                        ETAMAX,ETADIF,CONV,NFCNJ,T1,T2,T3,IFAIL)
                IF (JACGEN.EQ.2) 
     $            CALL N1JACB (FCN, N, M1, ML, X, F, A, XW,  AJDEL,
     $                         AJMIN, NFCNJ, T1, T2, T3, IFAIL)
                IF (MPRTIM.NE.0) CALL MONOFF(2)
              ENDIF
            ENDIF
            NJAC = NJAC + 1
C     Exit, If ...
            IF (JACGEN.EQ.1 .AND. IFAIL.LT.0) THEN
              IERR = 83
              GOTO 4299
            ENDIF
            IF (JACGEN.NE.1 .AND. IFAIL.NE.0) THEN
              IERR = 82
              GOTO 4299
            ENDIF
          ELSE IF (.NOT.QSIMPL) THEN
            NEW = NEW+1
          ENDIF
          IF ( NEW.EQ.0 .AND. (QLU.OR.NITER.EQ.0) ) THEN
C             ------------------------------------------------------
C             2.3.2.1 Save scaling values
              DO 2321 L1=1,N
                XWA(L1) = XW(L1)
2321          CONTINUE
C             ------------------------------------------------------
C             2.3.2.2 Prepare Jac. for use by band-solver DGBFA/DGBSL
              IF (MSTOR.EQ.1) THEN
                DO 2322 L1=1,N
                  DO 2323 L2=M2,1,-1
                    A(L2+ML,L1)=A(L2,L1)
2323              CONTINUE
2322            CONTINUE
              ENDIF
C             ------------------------------------------------------
C             2.4 Prepare solution of the linear system
C             ------------------------------------------------------
C             2.4.1 internal column scaling of matrix A
              IF (MSTOR.EQ.0) THEN
                DO 2410 K=1,N
                  S1 =-XW(K)
                  DO 2412 L1=1,N
                    A(L1,K)=A(L1,K)*S1
2412              CONTINUE
2410            CONTINUE
              ELSE IF (MSTOR.EQ.1) THEN
                DO 2413 K=1,N
                  L2=MAX0(1+M2-K,ML+1)
                  L3=MIN0(N+M2-K,M1)
                  S1 =-XW(K)
                  DO 2415 L1=L2,L3
                    A(L1,K)=A(L1,K)*S1
2415              CONTINUE
2413            CONTINUE
              ENDIF
C             ----------------------------------------------------
C             2.4.2 Row scaling of matrix A
              IF (QSCALE) THEN
                IF (MSTOR.EQ.0) THEN
                  CALL N1SCRF(N,N,A,FW)
                ELSE IF (MSTOR.EQ.1) THEN
                  CALL N1SCRB(N,M1,ML,MU,A,FW)
                ENDIF
              ELSE
                DO 242 K=1,N
                  FW(K)=ONE
242             CONTINUE
              ENDIF
          ENDIF
C         --------------------------------------------------------
C         2.4.3 Save and scale values of F(N)
          DO 243 L1=1,N
            FA(L1)=F(L1)
            T1(L1)=F(L1)*FW(L1)
243       CONTINUE
C         --------------------------------------------------------
C         3 Central part of iteration step
C         --------------------------------------------------------
C         3.1 Solution of the linear system
C         --------------------------------------------------------
C         3.1.1 Decomposition of (N,N)-matrix A
          IF (.NOT.QLINIT) THEN
            NIWLA = IWK(18)+1
            NRWLA = IWK(19)+1
            LIWL = LIWK-NIWLA+1
            LRWL = LRWK-NRWLA+1
          ENDIF
          IF ( NEW.EQ.0 .AND. (QLU.OR.NITER.EQ.0)) THEN
            IF (MPRTIM.NE.0) CALL MONON(3)
            CALL N1FACT(N,M1,ML,MU,A,IOPT,IFAIL,LIWL,IWK(NIWLA),
     $                  NILUSE,LRWL,RWK(NRWLA),NRLUSE)
            IF (MPRTIM.NE.0) CALL MONOFF(3)
            IF (.NOT.QLINIT) THEN
              NIWKFR = NIWKFR+NILUSE
              NRWKFR = NRWKFR+NRLUSE
C             Store lengths of currently required workspaces
              IWK(18) = NIWKFR-1
              IWK(19) = NRWKFR-1
            ENDIF
C       Exit Repeat If ...
            IF(IFAIL.NE.0) THEN
              IF (IFAIL.EQ.1) THEN
                IERR = 1
              ELSE
                IERR = 80
              ENDIF
              GOTO 4299
            ENDIF
          ENDIF
          QLINIT = .TRUE.
C         --------------------------------------------------------
C         3.1.2 Solution of linear (N,N)-system
          IF(NEW.EQ.0) THEN 
            IF (MPRTIM.NE.0) CALL MONON(4)
            CALL N1SOLV(N,M1,ML,MU,A,T1,IOPT,IFAIL,LIWL,IWK(NIWLA),
     $                 IDUMMY,LRWL,RWK(NRWLA),IDUMMY)
            IF (MPRTIM.NE.0) CALL MONOFF(4)
C     Exit Repeat If ...
            IF(IFAIL.NE.0)  THEN
              IERR = 81
              GOTO 4299
            ENDIF
          ELSE  
            ALFA1=ZERO
            ALFA2=ZERO
            DO 3121 I=1,N
              ALFA1=ALFA1+DX(I)*DXQ(I)/XW(I)**2
              ALFA2=ALFA2+DX(I)**2/XW(I)**2
3121        CONTINUE
            ALFA=ALFA1/ALFA2
            BETA=ONE-ALFA
            DO 3122 I=1,N
              T1(I)=(DXQ(I)+(FCA-ONE)*ALFA*DX(I))/BETA
3122        CONTINUE
            IF(NEW.EQ.1) THEN
              DO 3123 I=1,N
                DXSAVE(I,1)=DX(I)
3123          CONTINUE
            ENDIF
            DO 3124 I=1,N
              DXSAVE(I,NEW+1)=T1(I)
              DX(I)=T1(I)
              T1(I)=T1(I)/XW(I)
3124        CONTINUE
          ENDIF
C         --------------------------------------------------------
C         3.2 Evaluation of scaled natural level function SUMX
C             scaled maximum error norm CONV
C             evaluation of (scaled) standard level function
C             DLEVF ( DLEVF only, if MPRMON.GE.2 )
C             and computation of ordinary Newton corrections 
C             DX(N)
          IF (.NOT. QSIMPL) THEN
            CALL N1LVLS(N,T1,XW,F,DX,CONV,SUMX,DLEVF,MPRMON,NEW.EQ.0)
          ELSE
            CALL N1LVLS(N,T1,XWA,F,DX,CONV,SUMX,DLEVF,MPRMON,NEW.EQ.0)
          ENDIF
          DO 32 L1=1,N
            XA(L1)=X(L1)
32        CONTINUE
          SUMXA = SUMX
          DLEVXA = DSQRT(SUMXA/DBLE(FLOAT(N)))
          CONVA = CONV
          DXANRM = WNORM(N,DX,XW)
C         --------------------------------------------------------
C         3.3 A - priori estimate of damping factor FC
          IF(NITER.NE.0.AND.NONLIN.NE.1.AND.NEW.EQ.0.AND.
     $       .NOT. QORDI)THEN
C           ------------------------------------------------------
C           3.3.1 Computation of the denominator of a-priori
C                 estimate
            FCDNM = ZERO
            DO 331 L1=1,N
              FCDNM=FCDNM+((DX(L1)-DXQA(L1))/XW(L1))**2
331         CONTINUE
            FCDNM = FCDNM*SUMX
C           ------------------------------------------------------
C           3.3.2 New damping factor
            IF(FCDNM.GT.FCNUMP*FCMIN2 .OR.
     $        (NONLIN.EQ.4 .AND. FCA**2*FCNUMP .LT. 4.0D0*FCDNM)) THEN
              DMYPRI = FCA*DSQRT(FCNUMP/FCDNM)
              FCPRI = DMIN1(DMYPRI,ONE)
              IF (NONLIN.EQ.4) FCPRI = DMIN1(HALF*DMYPRI,ONE)
            ELSE
              FCPRI = ONE
C$Test-begin
              DMYPRI = -1.0D0
C$Test-end
            ENDIF
C$Test-begin
            IF (MPRMON.GE.5) THEN
              WRITE(LUMON,33201) FCPRI, FC, FCA, DMYPRI, FCNUMP,
     $                           FCDNM
33201         FORMAT (/, ' +++ apriori estimate +++', /,
     $                ' FCPRI  = ', D18.10, '  FC     = ', D18.10, /,
     $                ' FCA    = ', D18.10, '  DMYPRI = ', D18.10, /,
     $                ' FCNUMP = ', D18.10, '  FCDNM  = ', D18.10, /,
     $                   ' ++++++++++++++++++++++++', /)
            ENDIF
C$Test-end 
            FC = DMAX1(FCPRI,FCMIN)
            IF (QBDAMP) THEN
              FCBH = FCA*FCBND
              IF (FC.GT.FCBH) THEN
                FC = FCBH
                IF (MPRMON.GE.4)
     $            WRITE(LUMON,*) ' *** incr. rest. act. (a prio) ***'
              ENDIF
              FCBH = FCA/FCBND
              IF (FC.LT.FCBH) THEN
                FC = FCBH
                IF (MPRMON.GE.4)
     $            WRITE(LUMON,*) ' *** decr. rest. act. (a prio) ***'
              ENDIF
            ENDIF
          ENDIF
CWEI
          IF (IORMON.GE.2) THEN
            SUMXA2=SUMXA1
            SUMXA1=SUMXA0
            SUMXA0=DLEVXA
            IF (SUMXA0.EQ.ZERO) SUMXA0=SMALL
C           Check convergence rates (linear and superlinear)
C           ICONV : Convergence indicator
C                   =0: No convergence indicated yet
C                   =1: Damping factor is 1.0d0
C                   =2: Superlinear convergence detected (alpha >=1.2)
C                   =3: Quadratic convergence detected (alpha > 1.8)
            FCMON = DMIN1(FC,FCMON)
            IF (FCMON.LT.ONE) THEN
              ICONV = 0
              ALPHAE = ZERO
            ENDIF
            IF (FCMON.EQ.ONE .AND. ICONV.EQ.0) ICONV=1
            IF (NITER.GE.1) THEN
              CLIN1 = CLIN0
              CLIN0 = SUMXA0/SUMXA1
            ENDIF
            IF (ICONV.GE.1.AND.NITER.GE.2) THEN
              ALPHAK = ALPHAE
              ALPHAE = ZERO
              IF (CLIN1.LE.0.95D0) ALPHAE = DLOG(CLIN0)/DLOG(CLIN1)
              IF (ALPHAK.NE.ZERO) ALPHAK =0.5D0*(ALPHAE+ALPHAK)
              ALPHAA = DMIN1(ALPHAK,ALPHAE)
              CALPHK = CALPHA
              CALPHA = ZERO
              IF (ALPHAE.NE.ZERO) CALPHA = SUMXA1/SUMXA2**ALPHAE
              SUMXTE = DSQRT(CALPHA*CALPHK)*SUMXA1**ALPHAK-SUMXA0
              IF (ALPHAA.GE.1.2D0 .AND. ICONV.EQ.1) ICONV = 2
              IF (ALPHAA.GT.1.8D0) ICONV = 3
              IF (MPRMON.GE.4)  WRITE(LUMON,32001) ICONV, ALPHAE, 
     $                            CALPHA, CLIN0, ALPHAK, SUMXTE
32001         FORMAT(' ** ICONV: ',I1,'  ALPHA: ',D9.2,
     $               '  CONST-ALPHA: ',D9.2,'  CONST-LIN: ',D9.2,' **',
     $               /,' **',11X,'ALPHA-POST: ',D9.2,' CHECK: ',D9.2,
     $               25X,'**')
              IF ( ICONV.GE.2 .AND. ALPHAA.LT.0.9D0 ) THEN
                 IF (IORMON.EQ.3) THEN
                   IERR = 4
                   GOTO 4299
                 ELSE
                   QMSTOP = .TRUE.
                 ENDIF 
              ENDIF
            ENDIF
          ENDIF
          FCMON = FC
C
C         --------------------------------------------------------
C         3.4 Save natural level for later computations of
C             corrector and print iterate
          FCNUMK = SUMX
          IF (MPRMON.GE.2) THEN
            IF (MPRTIM.NE.0) CALL MONON(5)
            CALL N1PRV1(DLEVF,DLEVXA,FCKEEP,NITER,NEW,MPRMON,LUMON,
     $                  QMIXIO)
            IF (MPRTIM.NE.0) CALL MONOFF(5)
          ENDIF
          NRED = 0
          QNEXT = .FALSE.
          QREP  = .FALSE.   
C         QREP = ITER .GT. ITMAX   or  QREP = ITER .GT. 0
C
C         Damping-factor reduction loop
C         ================================
C         DO (Until)
34        CONTINUE
C           ------------------------------------------------------
C           3.5 Preliminary new iterate
            DO 35 L1=1,N
              X(L1)=XA(L1)+DX(L1)*FC
35          CONTINUE
C           -----------------------------------------------------
C           3.5.2 Exit, if problem is specified as being linear
C     Exit Repeat If ...
            IF( NONLIN.EQ.1 )THEN
              IERR = 0
              GOTO 4299
            ENDIF
C           ------------------------------------------------------
C           3.6.1 Computation of the residual vector
            IF (MPRTIM.NE.0) CALL MONON(1)
            CALL FCN(N,X,F,IFAIL)
            IF (MPRTIM.NE.0) CALL MONOFF(1)
            NFCN = NFCN+1
C     Exit, If ...
            IF(IFAIL.LT.0)THEN
              IERR = 82
              GOTO 4299
            ENDIF
            IF(IFAIL.EQ.1 .OR. IFAIL.EQ.2) THEN
              IF (QORDI) THEN
                IERR = 82
                GOTO 4299
              ENDIF
              IF (IFAIL.EQ.1) THEN
                FCREDU = HALF
              ELSE
                FCREDU = F(1)
C     Exit, If ...
                IF (FCREDU.LE.0 .OR. FCREDU.GE.1) THEN
                  IERR = 82
                  GOTO 4299
                ENDIF
              ENDIF
              IF (MPRMON.GE.2) THEN
36101           FORMAT(8X,I2,' FCN could not be evaluated  ',
     $                 8X,F7.5,4X,I2)
                WRITE(LUMON,36101)NITER,FC,NEW
              ENDIF
              FCH = FC
              FC = FCREDU*FC
              IF (FCH.GT.FCMIN) FC = DMAX1(FC,FCMIN)
              IF (QBDAMP) THEN
                FCBH = FCH/FCBND
                IF (FC.LT.FCBH) THEN
                  FC = FCBH
                  IF (MPRMON.GE.4) WRITE(LUMON,*)
     $               ' *** decr. rest. act. (FCN redu.) ***'
                ENDIF
              ENDIF
              IF (FC.LT.FCMIN) THEN
                IERR = 3
                GOTO 4299
              ENDIF  
C     Break DO (Until) ...
              GOTO 3109
            ENDIF
            IF (QORDI) THEN
C             -------------------------------------------------------
C             3.6.2 Convergence test for ordinary Newton iteration
C     Exit Repeat If ...
              IF( DXANRM.LE.RTOL )THEN
                IERR = 0
                GOTO 4299
              ENDIF
            ELSE
              DO 361 L1=1,N
               T1(L1)=F(L1)*FW(L1)
361           CONTINUE
C             --------------------------------------------------
C             3.6.3 Solution of linear (N,N)-system
              IF (MPRTIM.NE.0) CALL MONON(4)
              CALL N1SOLV(N,M1,ML,MU,A,T1,IOPT,IFAIL,LIWL,IWK(NIWLA),
     $                   IDUMMY,LRWL,RWK(NRWLA),IDUMMY)
              IF (MPRTIM.NE.0) CALL MONOFF(4)
C     Exit Repeat If ...
              IF(IFAIL.NE.0)  THEN
                IERR = 81
                GOTO 4299
              ENDIF
              IF(NEW.GT.0) THEN 
                DO 3630 I=1,N
                  DXQ(I) = T1(I)*XWA(I)
3630            CONTINUE                   
                DO 363 ILOOP=1,NEW 
                  SUM1=ZERO
                  SUM2=ZERO
                  DO 3631 I=1,N
                    SUM1=SUM1+(DXQ(I)*DXSAVE(I,ILOOP))/ XW(I)**2
                    SUM2=SUM2+(DXSAVE(I,ILOOP)/XW(I))**2
3631              CONTINUE
                  BETA=SUM1/SUM2
                  DO 3632 I=1,N
                    DXQ(I)=DXQ(I)+BETA*DXSAVE(I,ILOOP+1)
                    T1(I) = DXQ(I)/XW(I)
3632              CONTINUE
363             CONTINUE
              ENDIF
C             ----------------------------------------------------
C             3.6.4 Evaluation of scaled natural level function
C                   SUMX
C                   scaled maximum error norm CONV and evaluation
C                   of (scaled) standard level function DLEVFN
              IF (.NOT. QSIMPL) THEN
                CALL N1LVLS(N,T1,XW,F,DXQ,CONV,SUMX,DLEVFN,MPRMON,
     $                      NEW.EQ.0)
              ELSE
                CALL N1LVLS(N,T1,XWA,F,DXQ,CONV,SUMX,DLEVFN,MPRMON,
     $                      NEW.EQ.0)
              ENDIF  
              DXNRM = WNORM(N,DXQ,XW)
C             -----------------------------------------------------
C             3.6.5 Convergence test
C     Exit Repeat If ...
              IF ( DXNRM.LE.RTOL .AND. DXANRM.LE.RSMALL .AND. 
     $            FC.EQ.ONE ) THEN
                IERR = 0
                GOTO 4299
              ENDIF
C           
              FCA = FC
C             ----------------------------------------------------
C             3.6.6 Evaluation of reduced damping factor
              TH = FCA-ONE
              FCDNM = ZERO
              DO 39 L1=1,N
                FCDNM=FCDNM+((DXQ(L1)+TH*DX(L1))/XW(L1))**2
39            CONTINUE
              IF (FCDNM.NE.ZERO) THEN
                DMYCOR = FCA*FCA*HALF*DSQRT(FCNUMK/FCDNM)
              ELSE
                DMYCOR = 1.0D+35
              ENDIF
              IF (NONLIN.LE.3) THEN
                FCCOR = DMIN1(ONE,DMYCOR)
              ELSE
                FCCOR = DMIN1(ONE,HALF*DMYCOR)
              ENDIF
C$Test-begin
              IF (MPRMON.GE.5) THEN
                WRITE(LUMON,39001) FCCOR, FC, DMYCOR, FCNUMK,
     $                             FCDNM, FCA
39001           FORMAT (/, ' +++ corrector computation +++', /,
     $                ' FCCOR  = ', D18.10, '  FC     = ', D18.10, /,
     $                ' DMYCOR = ', D18.10, '  FCNUMK = ', D18.10, /,
     $                ' FCDNM  = ', D18.10, '  FCA    = ', D18.10, /,
     $                   ' +++++++++++++++++++++++++++++', /)
              ENDIF
C$Test-end 
            ENDIF
C           ------------------------------------------------------
C           3.7 Natural monotonicity test
            IF(SUMX.GT.SUMXA .AND. .NOT.QORDI)THEN
C             ----------------------------------------------------
C             3.8 Output of iterate
              IF(MPRMON.GE.3) THEN
                IF (MPRTIM.NE.0) CALL MONON(5)
                CALL N1PRV2(DLEVFN,DSQRT(SUMX/DBLE(FLOAT(N))),FC,
     $                      NITER,MPRMON,LUMON,QMIXIO,'*')
                IF (MPRTIM.NE.0) CALL MONOFF(5)
              ENDIF
              IF (QMSTOP) THEN
                IERR = 4
                GOTO 4299
              ENDIF
              FCH = DMIN1(FCCOR,HALF*FC)
              IF (FC.GT.FCMIN) THEN
                FC=DMAX1(FCH,FCMIN)
              ELSE
                FC=FCH
              ENDIF
              IF (QBDAMP) THEN
                FCBH = FCA/FCBND
                IF (FC.LT.FCBH) THEN
                  FC = FCBH
                  IF (MPRMON.GE.4)
     $              WRITE(LUMON,*) ' *** decr. rest. act. (a post) ***'
                ENDIF
              ENDIF
CWEI
              FCMON = FC
C
C$Test-begin
                IF (MPRMON.GE.5) THEN
                  WRITE(LUMON,39002) FC
39002             FORMAT (/, ' +++ corrector setting 1 +++', /,
     $                    ' FC     = ', D18.10, /,
     $                       ' +++++++++++++++++++++++++++', /)
                ENDIF
C$Test-end 
              QREP = .TRUE.
              NCORR = NCORR+1
              NRED = NRED+1
C             ----------------------------------------------------
C             3.10 If damping factor is too small:
C                  Refresh Jacobian,if current Jacobian was computed
C                  by a Rank1-update, else fail exit
              QJCRFR  = FC.LT.FCMIN.OR.NEW.GT.0.AND.NRED.GT.1
C     Exit Repeat If ...
              IF(QJCRFR.AND.NEW.EQ.0)THEN
                IERR = 3
                GOTO 4299
              ENDIF
            ELSE
              IF (.NOT.QORDI.AND..NOT.QREP .AND. FCCOR.GT.SIGMA2*FC)
     $        THEN
                IF(MPRMON.GE.3) THEN
                  IF (MPRTIM.NE.0) CALL MONON(5)
                  CALL N1PRV2(DLEVFN,DSQRT(SUMX/DBLE(FLOAT(N))),FC,
     $                        NITER,MPRMON,LUMON,QMIXIO,'+')
                  IF (MPRTIM.NE.0) CALL MONOFF(5)
                ENDIF
                FC = FCCOR
C$Test-begin
                IF (MPRMON.GE.5) THEN
                  WRITE(LUMON,39003) FC
39003             FORMAT (/, ' +++ corrector setting 2 +++', /,
     $                    ' FC     = ', D18.10, /,
     $                       ' +++++++++++++++++++++++++++', /)
                ENDIF
C$Test-end 
                QREP = .TRUE.
              ELSE
                QNEXT = .TRUE.
              ENDIF
            ENDIF
3109      CONTINUE
          IF(.NOT.(QNEXT.OR.QJCRFR)) GOTO  34
C         UNTIL ( expression - negated above)
C         End of damping-factor reduction loop
C         =======================================
          IF(QJCRFR)THEN
C           ------------------------------------------------------
C           3.11 Restore former values for repeting iteration
C                step
            NREJR1 = NREJR1+1
            DO 3111 L1=1,N
              X(L1)=XA(L1)
3111        CONTINUE
            DO 3112 L1=1,N
              F(L1)=FA(L1)
3112        CONTINUE
            IF(MPRMON.GE.2)THEN
31130           FORMAT(8X,I2,' Not accepted damping factor ',
     $                 8X,F7.5,4X,I2)
                WRITE(LUMON,31130)NITER,FC,NEW
            ENDIF
            FC = FCKEEP
            FCA = FCK2
            IF(NITER.EQ.0)THEN
              FC = FCMIN
            ENDIF
            QGENJ = .TRUE.
          ELSE
C           ------------------------------------------------------
C           4 Preparations to start the following iteration step
C           ------------------------------------------------------
C           4.1 Print values
            IF(MPRMON.GE.3 .AND. .NOT.QORDI) THEN
              IF (MPRTIM.NE.0) CALL MONON(5)
              CALL N1PRV2(DLEVFN,DSQRT(SUMX/DBLE(FLOAT(N))),FC,NITER+1,
     $                    MPRMON,LUMON,QMIXIO,'*')
              IF (MPRTIM.NE.0) CALL MONOFF(5)
            ENDIF
C           Print the natural level of the current iterate and return
C           it in one-step mode
            SUMXS = SUMX
            SUMX = SUMXA
            IF(MPRSOL.GE.2.AND.NITER.NE.0) THEN
              IF (MPRTIM.NE.0) CALL MONON(5)
              CALL N1SOUT(N,XA,2,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
              IF (MPRTIM.NE.0) CALL MONOFF(5)
            ELSE IF(MPRSOL.GE.1.AND.NITER.EQ.0)THEN
              IF (MPRTIM.NE.0) CALL MONON(5)
              CALL N1SOUT(N,XA,1,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
              IF (MPRTIM.NE.0) CALL MONOFF(5)
            ENDIF
            NITER = NITER+1
            DLEVF = DLEVFN
C     Exit Repeat If ...
            IF(NITER.GE.NITMAX)THEN
              IERR = 2
              GOTO 4299
            ENDIF
            FCKEEP = FC
C           ------------------------------------------------------
C           4.2 Return, if in one-step mode
C
C Exit Subroutine If ...
            IF (MODE.EQ.1) THEN
              IWK(18)=NIWLA-1
              IWK(19)=NRWLA-1
              IOPT(1)=1
              RETURN
            ENDIF
          ENDIF
        GOTO 2
C       End Repeat
4299    CONTINUE
C       End of main iteration loop
C       ==========================
C       ----------------------------------------------------------
C       9 Exits
C       ----------------------------------------------------------
C       9.1 Solution exit
        APREC = -1.0D0
C 
        IF(IERR.EQ.0 .OR. IERR.EQ.4)THEN
          IF (NONLIN.NE.1) THEN
            IF (.NOT.QORDI) THEN
              IF ( IERR.EQ.0 ) THEN
                APREC = DSQRT(SUMX/DBLE(FLOAT(N)))
                DO 91 L1=1,N
                  X(L1)=X(L1)+DXQ(L1)
91              CONTINUE
              ELSE 
                APREC = DSQRT(SUMXA/DBLE(FLOAT(N)))
                IF (ALPHAA.GT.ZERO .AND. IORMON.EQ.3) THEN
                  DO 92 L1=1,N
                    X(L1)=X(L1)+DX(L1)
92                CONTINUE
                ENDIF
              ENDIF
C             Print final monitor output
              IF(MPRMON.GE.2) THEN
                IF (IERR.EQ.0) THEN
                  IF (MPRTIM.NE.0) CALL MONON(5)
                  CALL N1PRV2(DLEVFN,DSQRT(SUMX/DBLE(FLOAT(N))),FC,
     $                        NITER+1,MPRMON,LUMON,QMIXIO,'*')
                  IF (MPRTIM.NE.0) CALL MONOFF(5)
                ELSE IF (IORMON.EQ.3) THEN
                  IF (MPRTIM.NE.0) CALL MONON(5)
                  CALL N1PRV1(DLEVFN,DSQRT(SUMXA/DBLE(FLOAT(N))),FC,
     $                        NITER,NEW,MPRMON,LUMON,QMIXIO)
                  IF (MPRTIM.NE.0) CALL MONOFF(5)
                ENDIF
              ENDIF
              IF (  IORMON.GE.2 ) THEN
                IF ( ICONV.LE.1 .AND. ALPHAE .NE. ZERO 
     $                          .AND. ALPHAK .NE. ZERO ) IERR = 5
              ENDIF
            ELSE
C           IF (QORDI) THEN
              APREC = DSQRT(SUMXA/DBLE(FLOAT(N)))
            ENDIF
            IF(MPRMON.GE.1) THEN
91001         FORMAT(///' Solution of nonlinear system ',
     $        'of equations obtained within ',I3,
     $        ' iteration steps',//,' Achieved relative accuracy',D10.3)
              IF (QORDI .OR. IERR.EQ.4) THEN
                WRITE(LUMON,91001) NITER,APREC
              ELSE
                WRITE(LUMON,91001) NITER+1,APREC
              ENDIF 
            ENDIF
          ELSE
            IF(MPRMON.GE.1) THEN
91002         FORMAT(///' Solution of linear system ',
     $        'of equations obtained by NLEQ1',//,' No estimate ',
     $        'available for the achieved relative accuracy')
                WRITE(LUMON,91002)
            ENDIF
          ENDIF
        ENDIF
C       ----------------------------------------------------------
C       9.2 Fail exit messages
C       ----------------------------------------------------------
C       9.2.1 Termination, since jacobian matrix became singular
        IF(IERR.EQ.1.AND.MPRERR.GE.1)THEN
92101     FORMAT(/,' Iteration terminates due to ',
     $    'singular jacobian matrix',/)
          WRITE(LUERR,92101)
        ENDIF
C       ----------------------------------------------------------
C       9.2.2 Termination after more than NITMAX iterations
        IF(IERR.EQ.2.AND.MPRERR.GE.1)THEN
92201     FORMAT(/,' Iteration terminates after NITMAX ',
     $    '=',I3,'  Iteration steps')
          WRITE(LUERR,92201)NITMAX
        ENDIF
C       ----------------------------------------------------------
C       9.2.3 Damping factor FC became too small
        IF(IERR.EQ.3.AND.MPRERR.GE.1)THEN
92301     FORMAT(/,' Damping factor has become too ',
     $    'small: lambda =',D10.3,2X,/)
          WRITE(LUERR,92301)FC
        ENDIF
CWEI
C       ----------------------------------------------------------
C       9.2.4.1 Superlinear convergence slowed down
        IF(IERR.EQ.4.AND.MPRERR.GE.1)THEN
92401     FORMAT(/,' Warning: Monotonicity test failed after ',A,
     $           ' convergence was already checked;',/,
     $    ' RTOL requirement may be too stringent',/)
92402     FORMAT(/,' Warning: ',A,' convergence slowed down;',/,
     $    ' RTOL requirement may be too stringent',/)
          IF (QMSTOP) THEN
            IF (ICONV.EQ.2) WRITE(LUERR,92401) 'superlinear'
            IF (ICONV.EQ.3) WRITE(LUERR,92401) 'quadratic'
          ELSE
            IF (ICONV.EQ.2) WRITE(LUERR,92402) 'superlinear'
            IF (ICONV.EQ.3) WRITE(LUERR,92402) 'quadratic'
          ENDIF
        ENDIF
C       ----------------------------------------------------------
C       9.2.4.2 Convergence criterion satisfied before superlinear
C               convergence has been established
        IF (IERR.EQ.5.AND.DLEVFN.EQ.ZERO) IERR=0
        IF(IERR.EQ.5.AND.MPRERR.GE.1)THEN
92410     FORMAT(/,' Warning: No quadratic or superlinear convergence ',
     $           'established yet',/,
     $           10X,'your solution may perhaps may be less accurate ',
     $           /,10X,'as indicated by the standard error estimate')
          WRITE(LUERR,92410)
        ENDIF
C       ----------------------------------------------------------
C       9.2.5 Error exit due to linear solver routine N1FACT
        IF(IERR.EQ.80.AND.MPRERR.GE.1)THEN
92501     FORMAT(/,' Error ',I5,' signalled by linear solver N1FACT')
          WRITE(LUERR,92501) IFAIL
        ENDIF
C       ----------------------------------------------------------
C       9.2.6 Error exit due to linear solver routine N1SOLV
        IF(IERR.EQ.81.AND.MPRERR.GE.1)THEN
92601     FORMAT(/,' Error ',I5,' signalled by linear solver N1SOLV')
          WRITE(LUERR,92601) IFAIL
        ENDIF
C       ----------------------------------------------------------
C       9.2.7 Error exit due to fail of user function FCN
        IF(IERR.EQ.82.AND.MPRERR.GE.1)THEN
92701     FORMAT(/,' Error ',I5,' signalled by user function FCN')
          WRITE(LUERR,92701) IFAIL
        ENDIF
C       ----------------------------------------------------------
C       9.2.8 Error exit due to fail of user function JAC
        IF(IERR.EQ.83.AND.MPRERR.GE.1)THEN
92801     FORMAT(/,' Error ',I5,' signalled by user function JAC')
          WRITE(LUERR,92801) IFAIL
        ENDIF
        IF(IERR.GE.80.AND.IERR.LE.83) IWK(23) = IFAIL
        IF ((IERR.EQ.82.OR.IERR.EQ.83).AND.NITER.LE.1.AND.MPRERR.GE.1)
     $  THEN
          WRITE (LUERR,92810)
92810     FORMAT(' Try to find a better initial guess for the solution')
        ENDIF
C       ----------------------------------------------------------
C       9.3 Common exit
        IF (MPRERR.GE.3.AND.IERR.NE.0.AND.IERR.NE.4.AND.NONLIN.NE.1)
     $    THEN
93100     FORMAT(/,'    Achieved relative accuracy',D10.3,2X)
          WRITE(LUERR,93100)CONVA
          APREC = CONVA
        ENDIF
        RTOL = APREC
        SUMX = SUMXA
        IF(MPRSOL.GE.2.AND.NITER.NE.0) THEN
           MODE=2
           IF (QORDI) MODE=3
           IF (MPRTIM.NE.0) CALL MONON(5)
           CALL N1SOUT(N,XA,MODE,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
           IF (MPRTIM.NE.0) CALL MONOFF(5)
        ELSE IF(MPRSOL.GE.1.AND.NITER.EQ.0)THEN
           IF (MPRTIM.NE.0) CALL MONON(5)
           CALL N1SOUT(N,XA,1,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
           IF (MPRTIM.NE.0) CALL MONOFF(5)
        ENDIF
        IF (.NOT.QORDI) THEN
          IF (IERR.NE.4) NITER = NITER+1
          DLEVF = DLEVFN
          IF(MPRSOL.GE.1)THEN
C           Print Solution or final iteration vector
            IF(IERR.EQ.0)THEN
               MODEFI = 3
            ELSE
               MODEFI = 4
            ENDIF
            IF (MPRTIM.NE.0) CALL MONON(5)
            CALL N1SOUT(N,X,MODEFI,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
            IF (MPRTIM.NE.0) CALL MONOFF(5)
          ENDIF
        ENDIF
C       Return the latest internal scaling to XSCAL
        DO 93 I=1,N
          XSCAL(I)=XW(I)
93      CONTINUE
C       End of exits
C       End of subroutine N1INT
      RETURN
      END
C
      SUBROUTINE N1SCAL(N,X,XA,XSCAL,XW,ISCAL,QINISC,IOPT,LRWK,RWK)
C*    Begin Prologue SCALE
      INTEGER N
      DOUBLE PRECISION X(N),XSCAL(N),XA(N),XW(N)
      INTEGER ISCAL
      LOGICAL QINISC
      INTEGER IOPT(50),LRWK
      DOUBLE PRECISION RWK(LRWK)
C     ------------------------------------------------------------
C
C*    Summary :
C    
C     S C A L E : To be used in connection with NLEQ1 .
C       Computation of the internal scaling vector XW used for the
C       Jacobian matrix, the iterate vector and it's related
C       vectors - especially for the solution of the linear system
C       and the computations of norms to avoid numerical overflow.
C
C*    Input parameters
C     ================
C
C     N         Int     Number of unknowns
C     X(N)      Dble    Current iterate
C     XA(N)     Dble    Previous iterate
C     XSCAL(N)  Dble    User scaling passed from parameter XSCAL
C                       of interface routine NLEQ1
C     ISCAL     Int     Option ISCAL passed from IOPT-field
C                       (for details see description of IOPT-fields)
C     QINISC    Logical = .TRUE.  : Initial scaling
C                       = .FALSE. : Subsequent scaling
C     IOPT(50)  Int     Options array passed from NLEQ1 parameter list
C     LRWK      Int     Length of real workspace
C     RWK(LRWK) Dble    Real workspace (see description above)
C
C*    Output parameters
C     =================
C
C     XW(N)     Dble   Scaling vector computed by this routine
C                      All components must be positive. The follow-
C                      ing relationship between the original vector
C                      X and the scaled vector XSCAL holds:
C                      XSCAL(I) = X(I)/XW(I) for I=1,...N
C
C*    Subroutines called: ZIBCONST
C
C*    Machine constants used
C     ======================
C
      DOUBLE PRECISION EPMACH, SMALL
C
C     ------------------------------------------------------------
C*    End Prologue
      INTRINSIC DABS,DMAX1
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
      INTEGER MPRMON,LUMON
      CALL ZIBCONST(EPMACH, SMALL)
C*    Begin
      DO 1 L1=1,N
        IF (ISCAL.EQ.1) THEN
          XW(L1) = XSCAL(L1)
        ELSE
          XW(L1)=DMAX1(XSCAL(L1),(DABS(X(L1))+DABS(XA(L1)))*HALF,SMALL)
        ENDIF
1     CONTINUE
C$Test-Begin
      MPRMON = IOPT(13)
      IF (MPRMON.GE.6) THEN
        LUMON = IOPT(14)
        WRITE(LUMON,*) ' '
        WRITE(LUMON,*) ' ++++++++++++++++++++++++++++++++++++++++++'
        WRITE(LUMON,*) '      X-components   Scaling-components    '
        WRITE(LUMON,10) (X(L1), XW(L1), L1=1,N)
10      FORMAT('  ',D18.10,'  ',D18.10)
        WRITE(LUMON,*) ' ++++++++++++++++++++++++++++++++++++++++++'
        WRITE(LUMON,*) ' '
      ENDIF
C$Test-End
C     End of subroutine N1SCAL
      RETURN
      END
C
      SUBROUTINE N1SCRF(M,N,A,FW)
C*    Begin Prologue SCROWF
      INTEGER M,N
      DOUBLE PRECISION A(M,N),FW(M)
C     ------------------------------------------------------------
C
C*    Summary :
C
C     S C R O W F : Row Scaling of a (M,N)-matrix in full storage
C                   mode
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C       M           Int    Number of rows of the matrix
C       N           Int    Number of columns of the matrix
C     * A(M,N)      Dble   Matrix to be scaled
C
C*    Output parameters
C     =================
C
C       FW(M)       Dble   Row scaling factors - FW(i) contains
C                          the factor by which the i-th row of A
C                          has been multiplied
C
C     ------------------------------------------------------------
C*    End Prologue
      INTRINSIC DABS
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER J,K
      DOUBLE PRECISION S1,S2
C*    Begin
      DO 1 K=1,M
        S1=ZERO
        DO 2 J=1,N
          S2=DABS(A(K,J))
          IF (S2.GT.S1) S1=S2
2       CONTINUE
        IF (S1.GT.ZERO) THEN
          S1=ONE/S1
          FW(K)=S1
          DO 3 J=1,N
            A(K,J)=A(K,J)*S1
3         CONTINUE
        ELSE
          FW(K)=ONE
        ENDIF
1     CONTINUE
C     End of subroutine N1SCRF
      RETURN
      END
C
      SUBROUTINE N1SCRB(N,LDA,ML,MU,A,FW)
C*    Begin Prologue SCROWB
      INTEGER N,LDA,ML,MU
      DOUBLE PRECISION A(LDA,N),FW(N)
C     ------------------------------------------------------------
C
C*    Summary :
C
C     S C R O W B : Row Scaling of a (N,N)-matrix in band storage
C                   mode
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C       N           Int    Number of rows and columns of the matrix
C       LDA         Int    Leading dimension of the matrix array
C       ML          Int    Lower bandwidth of the matrix (see IOPT)
C       MU          Int    Upper bandwidth of the matrix (see IOPT)
C     * A(LDA,N)    Dble   Matrix to be scaled (see Note in routine
C                          DGBFA for details on storage technique)
C
C*    Output parameters
C     =================
C
C       FW(N)       Dble   Row scaling factors - FW(i) contains
C                          the factor by which the i-th row of
C                          the matrix has been multiplied
C
C     ------------------------------------------------------------
C*    End Prologue
      INTRINSIC DABS,MAX0,MIN0
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER K,K1,L1,L2,L3,M2
      DOUBLE PRECISION S1,S2
C*    Begin
      M2=ML+MU+1
      DO 1 K=1,N
        S1=ZERO
        L2=MAX0(1,K-ML)
        L3=MIN0(N,K+MU)
        K1=M2+K
        DO 2 L1=L2,L3
          S2=DABS(A(K1-L1,L1))
          IF (S2.GT.S1) S1=S2
2       CONTINUE
        IF (S1.GT.ZERO) THEN
          S1=ONE/S1
          FW(K)=S1
          DO 3 L1=L2,L3
            A(K1-L1,L1)=A(K1-L1,L1)*S1
3         CONTINUE
        ELSE
          FW(K)=ONE
        ENDIF
1     CONTINUE
C     End of subroutine N1SCRB
      RETURN
      END
C
      SUBROUTINE N1FACT(N,LDA,ML,MU,A,IOPT,IFAIL,LIWK,IWK,LAIWK,LRWK,
     $RWK,LARWK)
C*    Begin Prologue FACT
      INTEGER N,LDA,ML,MU
      DOUBLE PRECISION A(LDA,N)
      INTEGER IOPT(50)
      INTEGER IFAIL
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LAIWK,LRWK
      DOUBLE PRECISION RWK(LRWK)
      INTEGER LARWK
C     ------------------------------------------------------------
C
C*    Summary :
C
C     F A C T : Call linear algebra subprogram for factorization of
C               a (N,N)-matrix
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C     N             Int    Order of the linear system
C     LDA           Int    Leading dimension of the matrix array A
C     ML            Int    Lower bandwidth of the matrix (only for
C                          banded systems)
C     MU            Int    Upper bandwidth of the matrix (only for
C                          banded systems)
C   * A(LDA,N)      Dble   Matrix storage. See main routine NLEQ1,
C                          note 4 for how to store a banded Jacobian
C     IOPT(50)      Int    Option vector passed from NLEQ1
C
C*    Output parameters
C     =================
C
C     IFAIL         Int    Error indicator returned by this routine:
C                          = 0 matrix decomposition successfull
C                          = 1 decomposition failed - 
C                              matrix is numerically singular
C                          =10 supplied (integer) workspace too small
C
C*    Workspace parameters
C     ====================
C
C     LIWK          Int    Length of integer workspace passed to this
C                          routine (In)
C     IWK(LIWK)     Int    Integer Workspace supplied for this routine
C     LAIWK         Int    Length of integer Workspace used by this 
C                          routine (out)       
C     LRWK          Int    Length of real workspace passed to this
C                          routine (In)                  
C     RWK(LRWK)     Dble   Real Workspace supplied for this routine
C     LARWK         Int    Length of real Workspace used by this 
C                          routine (out)
C
C*    Subroutines called:  DGETRF, DGBTRF
C
C     ------------------------------------------------------------
C*    End Prologue
      EXTERNAL DGETRF, DGBTRF
C*    Begin
      LAIWK = N
      LARWK = 0
C     ( Real workspace starts at RWK(NRWKFR) )
C     IF (LIWK.GE.LAIWK.AND.LRWK.GE.LARWK) THEN
      IF (LIWK.GE.LAIWK) THEN
        MSTOR = IOPT(4)
        IF (MSTOR.EQ.0) THEN
          CALL DGETRF(N,N,A,LDA,IWK,IFAIL)
        ELSE IF (MSTOR.EQ.1) THEN
          CALL DGBTRF(N,N,ML,MU,A,LDA,IWK,IFAIL)
        ENDIF
        IF(IFAIL.NE.0)THEN
          IFAIL = 1
        ENDIF
      ELSE
        IFAIL = 10
        MPRERR=IOPT(11)
        LUERR=IOPT(12)
10001   FORMAT(/,' Insuffient workspace for linear solver,',
     $         ' at least needed more needed : ',/,
     $         ' ',A,' workspace : ',I4)
        IF (LIWK.LT.LAIWK.AND.MPRERR.GT.0) 
     $    WRITE(LUERR,10001) 'Integer',LAIWK-LIWK
      ENDIF
      RETURN
      END
C
      SUBROUTINE N1SOLV(N,LDA,ML,MU,A,B,IOPT,IFAIL,LIWK,IWK,LAIWK,
     $LRWK,RWK,LARWK)
C*    Begin Prologue SOLVE
      INTEGER N,LDA,ML,MU
      DOUBLE PRECISION A(LDA,N)
      DOUBLE PRECISION B(N)
      INTEGER IOPT(50)
      INTEGER IFAIL
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LRWK,LAIWK
      DOUBLE PRECISION RWK(LRWK)
      INTEGER LARWK
C     ------------------------------------------------------------
C
C*    Summary :
C
C     S O L V E : Call linear algebra subprogram for solution of
C                 the linear system A*Z = B
C
C*    Parameters
C     ==========
C
C     N,LDA,ML,MU,A,IOPT,IFAIL,LIWK,IWK,LAIWK,LRWK,RWK,LARWK :
C                        See description for subroutine N1FACT.          
C     B          Dble    In:  Right hand side of the linear system
C                        Out: Solution of the linear system
C
C     Subroutines called: DGETRS, DGBTRS
C
C     ------------------------------------------------------------
C*    End Prologue
      EXTERNAL DGETRS, DGBTRS
C*    Begin
      MSTOR = IOPT(4)
      IF (MSTOR.EQ.0) THEN
        CALL DGETRS('N',N,1,A,LDA,IWK,B,N,IFAIL)
      ELSE IF (MSTOR.EQ.1) THEN
        CALL DGBTRS('N',N,ML,MU,1,A,LDA,IWK,B,N,IFAIL)
      ENDIF
      RETURN
      END
C
      SUBROUTINE N1LVLS(N,DX1,XW,F,DXQ,CONV,SUMX,DLEVF,MPRMON,QDSCAL)
C*    Begin Prologue LEVELS
      INTEGER N,MPRMON
      DOUBLE PRECISION DX1(N),XW(N),F(N),DXQ(N)
      DOUBLE PRECISION CONV,SUMX,DLEVF
      LOGICAL QDSCAL
C
C     ------------------------------------------------------------
C
C*    Summary :
C
C     L E V E L S : To be used in connection with NLEQ1 .
C     provides descaled solution, error norm and level functions
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C       N              Int    Number of parameters to be estimated
C       DX1(N)         Dble   array containing the scaled Newton
C                             correction
C       XW(N)          Dble   Array containing the scaling values
C       F(N)           Dble   Array containing the residuum
C
C*    Output parameters
C     =================
C
C       DXQ(N)         Dble   Array containing the descaled Newton
C                             correction
C       CONV           Dble   Scaled maximum norm of the Newton
C                             correction
C       SUMX           Dble   Scaled natural level function value
C       DLEVF          Dble   Standard level function value (only
C                             if needed for print)
C       MPRMON         Int    Print information parameter (see
C                             driver routine NLEQ1 )
C       QDSCAL         Logic  .TRUE., if descaling of DX1 required,
C                             else .FALSE.
C
C     ------------------------------------------------------------
C*    End Prologue
      INTRINSIC DABS
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER L1
      DOUBLE PRECISION S1
C*    Begin
      IF (QDSCAL) THEN
C       ------------------------------------------------------------
C       1.2 Descaling of solution DX1 ( stored to DXQ )
        DO 12 L1=1,N
          DXQ(L1)=DX1(L1)*XW(L1)
12      CONTINUE
      ENDIF
C     ------------------------------------------------------------
C     2 Evaluation of scaled natural level function SUMX and
C       scaled maximum error norm CONV
      CONV = ZERO
      DO 20 L1=1,N
        S1 = DABS(DX1(L1))
        IF(S1.GT.CONV) CONV=S1
20    CONTINUE
      SUMX = ZERO
      DO 21 L1=1,N
        SUMX = SUMX+DX1(L1)**2
21    CONTINUE
C     ------------------------------------------------------------
C     3 Evaluation of (scaled) standard level function DLEVF
      DLEVF = ZERO
      DO 3 L1=1,N
        DLEVF = DLEVF+F(L1)**2
3     CONTINUE
      DLEVF = DSQRT(DLEVF/DBLE(FLOAT(N)))
C     End of subroutine N1LVLS
      RETURN
      END
C
      SUBROUTINE N1JAC (FCN, N, LDA, X, FX, A, YSCAL, AJDEL, AJMIN,
     $                  NFCN, FU, IFAIL)
C* Begin Prologue N1JAC
      EXTERNAL FCN
      INTEGER N, LDA
      DOUBLE PRECISION X(N), FX(N), A(LDA,N), YSCAL(N), AJDEL, AJMIN
      INTEGER NFCN
      DOUBLE PRECISION FU(N)
      INTEGER IFAIL
C
C  ---------------------------------------------------------------------
C
C* Title
C
C  Evaluation of a dense Jacobian matrix using finite difference
C  approximation adapted for use in nonlinear systems solver NLEQ1
C
C* Environment       Fortran 77
C                    Double Precision
C                    Sun 3/60, Sun OS
C* Latest Revision   May 1990
C
C
C* Parameter list description
C  --------------------------
C
C* External subroutines (to be supplied by the user)
C  -------------------------------------------------
C
C  FCN        Ext     FCN (N, X, FX, IFAIL)
C                     Subroutine in order to provide right-hand
C                     side of first-order differential equations
C    N        Int     Number of rows and columns of the Jacobian
C    X(N)     Dble    The current scaled iterates
C    FX(N)    Dble    Array containing FCN(X)
C    IFAIL    Int     Return code
C                     Whenever a negative value is returned by FCN
C                     routine N1JAC is terminated immediately.
C
C
C* Input parameters (* marks inout parameters)
C  ----------------
C
C  N          Int     Number of rows and columns of the Jacobian
C  LDA        Int     Leading Dimension of array A
C  X(N)       Dble    Array containing the current scaled
C                     iterate
C  FX(N)      Dble    Array containing FCN(X)
C  YSCAL(N)   Dble    Array containing the scaling factors
C  AJDEL      Dble    Perturbation of component k: abs(Y(k))*AJDEL
C  AJMIN      Dble    Minimum perturbation is AJMIN*AJDEL
C  NFCN       Int  *  FCN - evaluation count
C
C* Output parameters (* marks inout parameters)
C  -----------------
C
C  A(LDA,N)   Dble    Array to contain the approximated
C                     Jacobian matrix ( dF(i)/dx(j)in A(i,j))
C  NFCN       Int  *  FCN - evaluation count adjusted
C  IFAIL      Int     Return code non-zero if Jacobian could not
C                     be computed.
C
C* Workspace parameters
C  --------------------
C
C  FU(N)      Dble    Array to contain FCN(x+dx) for evaluation of
C                     the numerator differences
C
C* Called
C  ------
C
      INTRINSIC DABS, DMAX1, DSIGN
C  ---------------------------------------------------------------------
C
C* End Prologue
C
C* Local variables
C  ---------------
C
      INTEGER I, K
      DOUBLE PRECISION U, W
C
C* Begin
C
      IFAIL = 0
      DO 1 K = 1,N
         W = X(K)
         U = DSIGN(DMAX1(DABS(X(K)),AJMIN,YSCAL(K))*AJDEL, X(K))
         X(K) = W + U
C
         CALL FCN (N, X, FU, IFAIL)
         NFCN = NFCN + 1
         IF (IFAIL .NE. 0) GOTO 99
C
         X(K) = W
         DO 11 I = 1,N
            A(I,K) = (FU(I) - FX(I)) / U  
 11      CONTINUE
 1    CONTINUE
C
99    CONTINUE
      RETURN
C
C
C* End of N1JAC
C
      END
      SUBROUTINE N1JACB (FCN, N, LDA, ML, X, FX, A, YSCAL, AJDEL, AJMIN,
     $     NFCN, FU, U, W, IFAIL)
C* Begin Prologue N1JACB
      EXTERNAL FCN
      INTEGER N, LDA, ML
      DOUBLE PRECISION X(N), FX(N), A(LDA,N), YSCAL(N), AJDEL, AJMIN
      INTEGER NFCN
      DOUBLE PRECISION FU(N), U(N), W(N)
      INTEGER IFAIL
C
C  ---------------------------------------------------------------------
C
C* Title
C
C  Evaluation of a banded Jacobian matrix using finite difference
C  approximation adapted for use in nonlinear systems solver NLEQ1
C
C* Environment       Fortran 77
C                    Double Precision
C                    Sun 3/60, Sun OS
C* Latest Revision   May 1990
C
C
C* Parameter list description
C  --------------------------
C
C* External subroutines (to be supplied by the user)
C  -------------------------------------------------
C
C  FCN        Ext     FCN (N, X, FX, IFAIL)
C                     Subroutine in order to provide right-hand
C                     side of first-order differential equations
C    N        Int     Number of rows and columns of the Jacobian
C    X(N)     Dble    The current scaled iterates
C    FX(N)    Dble    Array containing FCN(X)
C    IFAIL    Int     Return code
C                     Whenever a negative value is returned by FCN
C                     routine N1JACB is terminated immediately.
C
C
C* Input parameters (* marks inout parameters)
C  ----------------
C
C  N          Int     Number of rows and columns of the Jacobian
C  LDA        Int     Leading Dimension of array A
C  ML         Int     Lower bandwidth of the Jacobian matrix
C  X(N)       Dble    Array containing the current scaled
C                     iterate
C  FX(N)      Dble    Array containing FCN(X)
C  YSCAL(N)   Dble    Array containing the scaling factors
C  AJDEL      Dble    Perturbation of component k: abs(Y(k))*AJDEL
C  AJMIN      Dble    Minimum perturbation is AJMIN*AJDEL
C  NFCN       Int  *  FCN - evaluation count
C
C* Output parameters (* marks inout parameters)
C  -----------------
C
C  A(LDA,N)   Dble    Array to contain the approximated
C                     Jacobian matrix ( dF(i)/dx(j)in A(i,j))
C  NFCN       Int  *  FCN - evaluation count adjusted
C  IFAIL      Int     Return code non-zero if Jacobian could not
C                     be computed.
C
C* Workspace parameters
C  --------------------
C
C  FU(N)      Dble    Array to contain FCN(x+dx) for evaluation of
C                     the numerator differences
C  U(N)       Dble    Array to contain dx(i)
C  W(N)       Dble    Array to save original components of X
C
C* Called
C  ------
C
      INTRINSIC DABS, DMAX1, MAX0, MIN0
C
C* Constants
C  ---------
C
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0 D0)
C
C  ---------------------------------------------------------------------
C
C* End Prologue
C
C* Local variables
C  ---------------
C
      INTEGER I, I1, I2, JJ, K, MH, MU
C
C* Begin
C
      IFAIL = 0
      MU = LDA - 2*ML - 1
      LDAB = ML + MU + 1
      DO 1 I = 1,LDAB
         DO 11 K = 1,N
            A(I,K) = ZERO
 11      CONTINUE
 1    CONTINUE
C
      DO 2 JJ = 1,LDAB
         DO 21 K = JJ,N,LDAB
            W(K) = X(K)
            U(K) = DSIGN(DMAX1(DABS(X(K)),AJMIN,YSCAL(K))*AJDEL, X(K))
            X(K) = W(K) + U(K)
 21      CONTINUE
C
         CALL FCN (N, X, FU, IFAIL)
         NFCN = NFCN + 1
         IF (IFAIL .NE. 0) GOTO 99
C
         DO 22 K = JJ,N,LDAB 
            X(K) = W(K)
            I1 = MAX0 (1, K - MU)
            I2 = MIN0 (N, K + ML)
            MH = MU + 1 - K
            DO 221 I = I1,I2
               A(MH+I,K) = (FU(I) - FX(I)) / U(K)
 221        CONTINUE
 22      CONTINUE
C
 2    CONTINUE
C
99    CONTINUE
      RETURN
C
C
C* End of N1JACB
C
      END
      SUBROUTINE N1JCF (FCN, N, LDA, X, FX, A, YSCAL, ETA, ETAMIN,
     $     ETAMAX, ETADIF, CONV, NFCN, FU, IFAIL)
C* Begin Prologue N1JCF
      EXTERNAL FCN
      INTEGER N, LDA
      DOUBLE PRECISION X(N), FX(N), A(LDA,N), YSCAL(N), ETA(N),
     $     ETAMIN, ETAMAX, ETADIF, CONV
      INTEGER NFCN
      DOUBLE PRECISION FU(N)
      INTEGER IFAIL
C
C  ---------------------------------------------------------------------
C
C* Title
C
C  Approximation of dense Jacobian matrix for nonlinear systems solver
C  NLEQ1 with feed-back control of discretization and rounding errors
C
C* Environment       Fortran 77
C                    Double Precision
C                    Sun 3/60, Sun OS
C* Latest Revision   May 1990
C
C
C* Parameter list description
C  --------------------------
C
C* External subroutines (to be supplied by the user)
C  -------------------------------------------------
C
C  FCN        Ext     FCN (N, X, FX, IFAIL)
C                     Subroutine in order to provide right-hand
C                     side of first-order differential equations
C    N        Int     Number of rows and columns of the Jacobian
C    X(N)     Dble    The current scaled iterates
C    FX(N)    Dble    Array containing FCN(X)
C    IFAIL    Int     Return code
C                     Whenever a negative value is returned by FCN
C                     routine N1JCF is terminated immediately.
C
C
C* Input parameters (* marks inout parameters)
C  ----------------
C
C  N          Int     Number of rows and columns of the Jacobian
C  LDA        Int     Leading dimension of A (LDA .GE. N)
C  X(N)       Dble    Array containing the current scaled
C                     iterate
C  FX(N)      Dble    Array containing FCN(X)
C  YSCAL(N)   Dble    Array containing the scaling factors
C  ETA(N)     Dble *  Array containing the scaled denominator
C                     differences
C  ETAMIN     Dble    Minimum allowed scaled denominator
C  ETAMAX     Dble    Maximum allowed scaled denominator
C  ETADIF     Dble    DSQRT (1.1*EPMACH)
C                     EPMACH = machine precision
C  CONV       Dble    Maximum norm of last (unrelaxed) Newton correction
C  NFCN       Int  *  FCN - evaluation count
C
C* Output parameters (* marks inout parameters)
C  -----------------
C
C  A(LDA,N)   Dble    Array to contain the approximated
C                     Jacobian matrix ( dF(i)/dx(j)in A(i,j))
C  ETA(N)     Dble *  Scaled denominator differences adjusted
C  NFCN       Int  *  FCN - evaluation count adjusted
C  IFAIL      Int     Return code non-zero if Jacobian could not
C                     be computed.
C
C* Workspace parameters
C  --------------------
C
C  FU(N)      Dble    Array to contain FCN(x+dx) for evaluation of
C                     the numerator differences
C
C* Called
C  ------
C
      INTRINSIC DABS, DMAX1, DMIN1, DSIGN, DSQRT
C
C* Constants
C  ---------
C
      DOUBLE PRECISION SMALL2, ZERO
      PARAMETER (SMALL2 = 0.1D0,
     $           ZERO   = 0.0D0)
C
C  ---------------------------------------------------------------------
C
C* End Prologue
C
C* Local variables
C  ---------------
C
      INTEGER I, K, IS
      DOUBLE PRECISION FHI, HG, U, SUMD, W
      LOGICAL QFINE
C
C* Begin
C
      DO 1 K = 1,N
         IS = 0
C        DO (Until)
 11         CONTINUE
            W = X(K)
            U = DSIGN (ETA(K)*YSCAL(K), X(K))
            X(K) = W + U
            CALL FCN (N, X, FU, IFAIL)
            NFCN = NFCN + 1
C           Exit, If ...
            IF (IFAIL .NE. 0) GOTO 99
            X(K) = W
            SUMD = ZERO
            DO 111 I = 1,N
               HG = DMAX1 (DABS (FX(I)), DABS (FU(I)))
               FHI = FU(I) - FX(I)
               IF (HG .NE. ZERO) SUMD = SUMD + (FHI/HG)**2
               A(I,K) = FHI / U
 111        CONTINUE
            SUMD = DSQRT (SUMD / DBLE(N))
            QFINE = .TRUE.
            IF (SUMD .NE. ZERO .AND. IS .EQ. 0)THEN
               ETA(K) = DMIN1 (ETAMAX,
     $              DMAX1 (ETAMIN, DSQRT (ETADIF / SUMD)*ETA(K)))
               IS = 1
               QFINE = CONV .LT. SMALL2 .OR. SUMD .GE. ETAMIN
            ENDIF
            IF (.NOT.(QFINE)) GOTO  11
C        UNTIL ( expression - negated above)
 1    CONTINUE
C
C     Exit from DO-loop
 99   CONTINUE
C
      RETURN
C
C* End of subroutine N1JCF
C
      END
      SUBROUTINE N1JCFB (FCN, N, LDA, ML, X, FX, A, YSCAL, ETA,
     $     ETAMIN, ETAMAX, ETADIF, CONV, NFCN, FU, U, W, IFAIL)
C* Begin Prologue N1JCFB
      EXTERNAL FCN
      INTEGER N, LDA, ML
      DOUBLE PRECISION X(N), FX(N), A(LDA,N), YSCAL(N), ETA(N),
     $     ETAMIN, ETAMAX, ETADIF, CONV
      INTEGER NFCN
      DOUBLE PRECISION FU(N), U(N), W(N)
      INTEGER IFAIL
C     
C     ---------------------------------------------------------------------
C     
C* Title
C  
C  Approximation of banded Jacobian matrix for nonlinear systems solver
C  NLEQ1 with feed-back control of discretization and rounding errors
C  
C* Environment       Fortran 77
C                    Double Precision
C                    Sun 3/60, Sun OS
C* Latest Revision   May 1990
C  
C  
C* Parameter list description
C  --------------------------
C  
C* External subroutines (to be supplied by the user)
C  -------------------------------------------------
C  
C  FCN        Ext     FCN (N, X, FX, IFAIL)
C                     Subroutine in order to provide right-hand
C                     side of first-order differential equations
C    N        Int     Number of rows and columns of the Jacobian
C    X(N)     Dble    The current scaled iterates
C    FX(N)    Dble    Array containing FCN(X)
C    IFAIL    Int     Return code
C                     Whenever a negative value is returned by FCN
C                     routine N1JCFB is terminated immediately.
C  
C  
C* Input parameters (* marks inout parameters)
C  ----------------
C  
C  N          Int     Number of rows and columns of the Jacobian
C  LDA        Int     Leading dimension of A (LDA .GE. ML+MU+1)
C  ML         Int     Lower bandwidth of Jacobian matrix
C  X(N)       Dble    Array containing the current scaled
C                     iterate
C  FX(N)      Dble    Array containing FCN(X)
C  YSCAL(N)   Dble    Array containing the scaling factors
C  ETA(N)     Dble *  Array containing the scaled denominator
C                     differences
C  ETAMIN     Dble    Minimum allowed scaled denominator
C  ETAMAX     Dble    Maximum allowed scaled denominator
C  ETADIF     Dble    DSQRT (1.1*EPMACH)
C                     EPMACH = machine precision
C  CONV       Dble    Maximum norm of last (unrelaxed) Newton correction
C  NFCN       Int  *  FCN - evaluation count
C  
C* Output parameters (* marks inout parameters)
C  -----------------
C  
C  A(LDA,N)   Dble    Array to contain the approximated
C                     Jacobian matrix ( dF(i)/dx(j)in A(i,j))
C  ETA(N)     Dble *  Scaled denominator differences adjusted
C  NFCN       Int  *  FCN - evaluation count adjusted
C  IFAIL      Int     Return code non-zero if Jacobian could not
C                     be computed.
C  
C* Workspace parameters
C  --------------------
C  
C  FU(N)      Dble    Array to contain FCN(x+dx) for evaluation of
C                     the numerator differences
C  U(N)       Dble    Array to contain dx(i)
C  W(N)       Dble    Array to save original components of X
C  
C* Called
C  ------
C  
      INTRINSIC DABS, DMAX1, DMIN1, DSIGN, DSQRT, MAX0, MIN0
C  
C* Constants
C  ---------
C  
      DOUBLE PRECISION SMALL2, ZERO
      PARAMETER (SMALL2 = 0.1D0,
     $           ZERO   = 0.0D0)
C  
C  ---------------------------------------------------------------------
C  
C* End Prologue
C  
C* Local variables
C  ---------------
C  
      INTEGER I, IS, I1, I2, JJ, K, MH, MU, LDAB
      DOUBLE PRECISION FHI, HG, SUMD
      LOGICAL QFINE
C  
C* Begin
C  
      MU = LDA - 2*ML - 1
      LDAB = ML + MU + 1
      DO 1 I = 1,LDAB
         DO 11 K = 1,N
            A(I,K) = ZERO
 11      CONTINUE
 1    CONTINUE
C  
      DO 2 JJ = 1,LDAB
         IS = 0
C        DO (Until)
 21         CONTINUE
            DO 211  K = JJ,N,LDAB
               W(K) = X(K)
               U(K) = DSIGN (ETA(K)*YSCAL(K), X(K))
               X(K) = W(K) + U(K)
 211        CONTINUE
C     
            CALL FCN (N, X, FU, IFAIL)
            NFCN = NFCN + 1
C           Exit, If ...
            IF (IFAIL .NE. 0) GOTO 99
C     
            DO 212 K = JJ,N,LDAB
               X(K) = W(K)
               SUMD = ZERO
               I1 = MAX0 (1, K - MU)
               I2 = MIN0 (N, K + ML)
               MH = MU + 1 - K
               DO 213 I = I1,I2
                  HG = DMAX1 (DABS(FX(I)), DABS(FU(I)))
                  FHI = FU(I) - FX(I)
                  IF (HG .NE. ZERO) SUMD = SUMD + (FHI/HG)**2
                  A(MH+I, K) = FHI / U(K)
 213           CONTINUE
               SUMD = DSQRT (SUMD / DBLE(N))
               QFINE = .TRUE.
               IF (SUMD .NE. ZERO .AND. IS .EQ. 0) THEN
                  ETA(K) = DMIN1 (ETAMAX,
     $                 DMAX1 (ETAMIN, DSQRT (ETADIF / SUMD)*ETA(K)))
                  IS = 1
                  QFINE = CONV .LT. SMALL2 .OR. SUMD .GE. ETAMIN
               ENDIF
 212        CONTINUE
            IF (.NOT.(QFINE)) GOTO 21
C        UNTIL ( expression - negated above)
 2    CONTINUE
C
C     Exit from DO-loop
 99   CONTINUE
C
      RETURN
C
C* End of subroutine N1JCFB
C
      END
C
      SUBROUTINE N1PRV1(DLEVF,DLEVX,FC,NITER,NEW,MPRMON,LUMON,QMIXIO)
C*    Begin Prologue N1PRV1
      DOUBLE PRECISION DLEVF,DLEVX,FC
      INTEGER NITER,NEW,MPRMON,LUMON
      LOGICAL QMIXIO
C     ------------------------------------------------------------
C
C*    Summary :
C
C     N 1 P R V 1 : Printing of intermediate values (Type 1 routine)
C
C*    Parameters
C     ==========
C
C     DLEVF, DLEVX   See descr. of internal double variables of N1INT
C     FC,NITER,NEW,MPRMON,LUMON
C                  See parameter descr. of subroutine N1INT
C     QMIXIO Logical  = .TRUE.  , if LUMON.EQ.LUSOL
C                     = .FALSE. , if LUMON.NE.LUSOL
C
C     ------------------------------------------------------------
C*    End Prologue
C     Print Standard - and natural level
      IF(QMIXIO)THEN
1       FORMAT(2X,66('*'))
        WRITE(LUMON,1)
2       FORMAT(8X,'It',7X,'Normf ',10X,'Normx ',20X,'New')
        IF (MPRMON.GE.3) WRITE(LUMON,2)
3       FORMAT(8X,'It',7X,'Normf ',10X,'Normx ',8X,'Damp.Fct.',3X,'New')
        IF (MPRMON.EQ.2) WRITE(LUMON,3)
      ENDIF
4     FORMAT(6X,I4,5X,D10.3,2X,4X,D10.3,17X,I2)
      IF (MPRMON.GE.3.OR.NITER.EQ.0) 
     $  WRITE(LUMON,4) NITER,DLEVF,DLEVX,NEW
5     FORMAT(6X,I4,5X,D10.3,6X,D10.3,6X,F7.5,4X,I2)
      IF (MPRMON.EQ.2.AND.NITER.NE.0) 
     $  WRITE(LUMON,5) NITER,DLEVF,DLEVX,FC,NEW
      IF(QMIXIO)THEN
6       FORMAT(2X,66('*'))
        WRITE(LUMON,6)
      ENDIF
C     End of subroutine N1PRV1
      RETURN
      END
C
      SUBROUTINE N1PRV2(DLEVF,DLEVX,FC,NITER,MPRMON,LUMON,QMIXIO,
     $                  CMARK)
C*    Begin Prologue N1PRV2
      DOUBLE PRECISION DLEVF,DLEVX,FC
      INTEGER NITER,MPRMON,LUMON
      LOGICAL QMIXIO
      CHARACTER*1 CMARK
C     ------------------------------------------------------------
C
C*    Summary :
C
C     N 1 P R V 2 : Printing of intermediate values (Type 2 routine)
C
C*    Parameters
C     ==========
C
C     DLEVF, DLEVX   See descr. of internal double variables of N2INT
C     FC,NITER,MPRMON,LUMON
C                  See parameter descr. of subroutine N2INT
C     QMIXIO Logical  = .TRUE.  , if LUMON.EQ.LUSOL
C                     = .FALSE. , if LUMON.NE.LUSOL
C     CMARK Char*1    Marker character to be printed before DLEVX
C
C     ------------------------------------------------------------
C*    End Prologue
C     Print Standard - and natural level, and damping
C     factor
      IF(QMIXIO)THEN
1       FORMAT(2X,66('*'))
        WRITE(LUMON,1)
2       FORMAT(8X,'It',7X,'Normf ',10X,'Normx ',8X,'Damp.Fct.')
        WRITE(LUMON,2)
      ENDIF
3     FORMAT(6X,I4,5X,D10.3,4X,A1,1X,D10.3,2X,4X,F7.5)
      WRITE(LUMON,3)NITER,DLEVF,CMARK,DLEVX,FC
      IF(QMIXIO)THEN
4       FORMAT(2X,66('*'))
        WRITE(LUMON,4)
      ENDIF
C     End of subroutine N1PRV2
      RETURN
      END
C
      SUBROUTINE N1SOUT(N,X,MODE,IOPT,RWK,NRW,IWK,NIW,MPRINT,LUOUT)
C*    Begin Prologue SOLOUT
      INTEGER N
      DOUBLE PRECISION X(N)
      INTEGER NRW
      INTEGER MODE
      INTEGER IOPT(50)
      DOUBLE PRECISION RWK(NRW)
      INTEGER NIW
      INTEGER IWK(NIW)
      INTEGER MPRINT,LUOUT
C     ------------------------------------------------------------
C
C*    Summary :
C
C     S O L O U T : Printing of iterate (user customizable routine)
C
C*    Input parameters
C     ================
C
C     N         Int Number of equations/unknowns
C     X(N)   Double iterate vector
C     MODE          =1 This routine is called before the first
C                      Newton iteration step
C                   =2 This routine is called with an intermedi-
C                      ate iterate X(N)
C                   =3 This is the last call with the solution
C                      vector X(N)
C                   =4 This is the last call with the final, but
C                      not solution vector X(N)
C     IOPT(50)  Int The option array as passed to the driver
C                   routine(elements 46 to 50 may be used
C                   for user options)
C     MPRINT    Int Solution print level 
C                   (see description of IOPT-field MPRINT)
C     LUOUT     Int the solution print unit 
C                   (see description of see IOPT-field LUSOL)
C
C
C*    Workspace parameters
C     ====================
C
C     NRW, RWK, NIW, IWK    see description in driver routine
C
C*    Use of IOPT by this routine
C     ===========================
C
C     Field 46:       =0 Standard output
C                     =1 GRAZIL suitable output
C
C     ------------------------------------------------------------
C*    End Prologue
      LOGICAL QGRAZ,QNORM
C*    Begin
      QNORM = IOPT(46).EQ.0
      QGRAZ = IOPT(46).EQ.1
      IF(QNORM) THEN
1        FORMAT('  ',A,' data:',/)
         IF (MODE.EQ.1) THEN
101        FORMAT('  Start data:',/,'  N =',I5,//,
     $            '  Format: iteration-number, (x(i),i=1,...N), ',
     $            'Normf , Normx ',/)
           WRITE(LUOUT,101) N
           WRITE(LUOUT,1) 'Initial'
         ELSE IF (MODE.EQ.3) THEN
           WRITE(LUOUT,1) 'Solution'
         ELSE IF (MODE.EQ.4) THEN
           WRITE(LUOUT,1) 'Final'
         ENDIF
2        FORMAT(' ',I5)
C        WRITE          NITER
         WRITE(LUOUT,2) IWK(1)
3        FORMAT((12X,3(D18.10,1X)))
         WRITE(LUOUT,3)(X(L1),L1=1,N)
C        WRITE          DLEVF,   DLEVX
         WRITE(LUOUT,3) RWK(19),DSQRT(RWK(18)/DBLE(FLOAT(N)))
         IF(MODE.EQ.1.AND.MPRINT.GE.2) THEN
           WRITE(LUOUT,1) 'Intermediate'
         ELSE IF(MODE.GE.3) THEN
           WRITE(LUOUT,1) 'End'
         ENDIF
      ENDIF
      IF(QGRAZ) THEN
        IF(MODE.EQ.1) THEN
10        FORMAT('&name com',I3.3,:,255(7(', com',I3.3,:),/))
          WRITE(LUOUT,10)(I,I=1,N+2)
15        FORMAT('&def  com',I3.3,:,255(7(', com',I3.3,:),/))
          WRITE(LUOUT,15)(I,I=1,N+2)
16        FORMAT(6X,': X=1, Y=',I3)
          WRITE(LUOUT,16) N+2
        ENDIF
20      FORMAT('&data ',I5)
C        WRITE          NITER
        WRITE(LUOUT,20) IWK(1) 
21      FORMAT((6X,4(D18.10)))
        WRITE(LUOUT,21)(X(L1),L1=1,N)
C        WRITE          DLEVF,   DLEVX
        WRITE(LUOUT,21) RWK(19),DSQRT(RWK(18)/DBLE(FLOAT(N)))
        IF(MODE.GE.3) THEN
30        FORMAT('&wktype 3111',/,'&atext x ''iter''')
          WRITE(LUOUT,30)
35        FORMAT('&vars = com',I3.3,/,'&atext y ''x',I3,'''',
     $           /,'&run')
          WRITE(LUOUT,35) (I,I,I=1,N)
36        FORMAT('&vars = com',I3.3,/,'&atext y ''',A,'''',
     $           /,'&run')
          WRITE(LUOUT,36) N+1,'Normf ',N+2,'Normx '
C39       FORMAT('&stop')
C         WRITE(LUOUT,39)
        ENDIF
      ENDIF
C     End of subroutine N1SOUT
      RETURN
      END
C*    End package
      DOUBLE PRECISION FUNCTION WNORM(N,Z,XW)
      INTEGER N
      DOUBLE PRECISION Z(N), XW(N)
C     ------------------------------------------------------------
C
C*    Summary :
C
C     E N O R M : Return the norm to be used in exit (termination)
C                 criteria
C
C*    Input parameters
C     ================
C
C     N         Int Number of equations/unknowns
C     Z(N)     Dble  The vector, of which the norm is to be computed
C     XW(N)    Dble  The scaling values of Z(N)
C
C*    Output
C     ======
C
C     WNORM(N,Z,XW)  Dble  The mean square root norm of Z(N) subject
C                          to the scaling values in XW(N):
C                          = Sqrt( Sum(1,...N)((Z(I)/XW(I))**2) / N )
C
C     ------------------------------------------------------------
C*    End Prologue
      INTEGER I
      DOUBLE PRECISION S
C*    Begin
      S = 0.0D0
      DO 10 I=1,N
        S = S + ( Z(I)/XW(I) ) ** 2
10    CONTINUE
      WNORM = DSQRT( S / DBLE(FLOAT(N)) )
C     End of function WNORM
      RETURN
      END
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of DGETRF
*
      END
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETF2 computes an LU factorization of a general m-by-n matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN 
      INTEGER            I, J, JP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH      
      INTEGER            IDAMAX
      EXTERNAL           DLAMCH, IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Compute machine safe minimum 
* 
      SFMIN = DLAMCH('S')  
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M ) THEN 
               IF( ABS(A( J, J )) .GE. SFMIN ) THEN 
                  CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 ) 
               ELSE 
                 DO 20 I = 1, M-J 
                    A( J+I, J ) = A( J+I, J ) / A( J, J ) 
   20            CONTINUE 
               END IF 
            END IF 
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGETF2
*
      END
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASWP performs a series of row interchanges on the matrix A.
*  One row interchange is initiated for each of rows K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (K2*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*
*  Further Details
*  ===============
*
*  Modified by
*   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
*
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
*
      RETURN
*
*     End of DLASWP
*
      END
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by DGETRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASWP, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( NOTRAN ) THEN
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand sides.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
*
*        Solve L*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        Solve A' * X = B.
*
*        Solve U'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve L'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,
     $               A, LDA, B, LDB )
*
*        Apply row interchanges to the solution vectors.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
*
      RETURN
*
*     End of DGETRS
*
      END
      SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGBTF2 computes an LU factorization of a real m-by-n band matrix A
*  using partial pivoting with row interchanges.
*
*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
*          On entry, the matrix A in band storage, in rows KL+1 to
*          2*KL+KU+1; rows 1 to KL of the array need not be set.
*          The j-th column of A is stored in the j-th column of the
*          array AB as follows:
*          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
*
*          On exit, details of the factorization: U is stored as an
*          upper triangular band matrix with KL+KU superdiagonals in
*          rows 1 to KL+KU+1, and the multipliers used during the
*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
*          See below for further details.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  Further Details
*  ===============
*
*  The band storage scheme is illustrated by the following example, when
*  M = N = 6, KL = 2, KU = 1:
*
*  On entry:                       On exit:
*
*      *    *    *    +    +    +       *    *    *   u14  u25  u36
*      *    *    +    +    +    +       *    *   u13  u24  u35  u46
*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
*     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
*     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
*
*  Array elements marked * are not used by the routine; elements marked
*  + need not be set on entry, but are required by the routine to store
*  elements of U, because of fill-in resulting from the row
*  interchanges.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, JP, JU, KM, KV
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     KV is the number of superdiagonals in the factor U, allowing for
*     fill-in.
*
      KV = KU + KL
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Gaussian elimination with partial pivoting
*
*     Set fill-in elements in columns KU+2 to KV to zero.
*
      DO 20 J = KU + 2, MIN( KV, N )
         DO 10 I = KV - J + 2, KL
            AB( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
*
*     JU is the index of the last column affected by the current stage
*     of the factorization.
*
      JU = 1
*
      DO 40 J = 1, MIN( M, N )
*
*        Set fill-in elements in column J+KV to zero.
*
         IF( J+KV.LE.N ) THEN
            DO 30 I = 1, KL
               AB( I, J+KV ) = ZERO
   30       CONTINUE
         END IF
*
*        Find pivot and test for singularity. KM is the number of
*        subdiagonal elements in the current column.
*
         KM = MIN( KL, M-J )
         JP = IDAMAX( KM+1, AB( KV+1, J ), 1 )
         IPIV( J ) = JP + J - 1
         IF( AB( KV+JP, J ).NE.ZERO ) THEN
            JU = MAX( JU, MIN( J+KU+JP-1, N ) )
*
*           Apply interchange to columns J to JU.
*
            IF( JP.NE.1 )
     $         CALL DSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1,
     $                     AB( KV+1, J ), LDAB-1 )
*
            IF( KM.GT.0 ) THEN
*
*              Compute multipliers.
*
               CALL DSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )
*
*              Update trailing submatrix within the band.
*
               IF( JU.GT.J )
     $            CALL DGER( KM, JU-J, -ONE, AB( KV+2, J ), 1,
     $                       AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ),
     $                       LDAB-1 )
            END IF
         ELSE
*
*           If pivot is zero, set INFO to the index of the pivot
*           unless a zero pivot has already been found.
*
            IF( INFO.EQ.0 )
     $         INFO = J
         END IF
   40 CONTINUE
      RETURN
*
*     End of DGBTF2
*
      END
      SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGBTRF computes an LU factorization of a real m-by-n band matrix A
*  using partial pivoting with row interchanges.
*
*  This is the blocked version of the algorithm, calling Level 3 BLAS.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
*          On entry, the matrix A in band storage, in rows KL+1 to
*          2*KL+KU+1; rows 1 to KL of the array need not be set.
*          The j-th column of A is stored in the j-th column of the
*          array AB as follows:
*          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
*
*          On exit, details of the factorization: U is stored as an
*          upper triangular band matrix with KL+KU superdiagonals in
*          rows 1 to KL+KU+1, and the multipliers used during the
*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
*          See below for further details.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  Further Details
*  ===============
*
*  The band storage scheme is illustrated by the following example, when
*  M = N = 6, KL = 2, KU = 1:
*
*  On entry:                       On exit:
*
*      *    *    *    +    +    +       *    *    *   u14  u25  u36
*      *    *    +    +    +    +       *    *   u13  u24  u35  u46
*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
*     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
*     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
*
*  Array elements marked * are not used by the routine; elements marked
*  + need not be set on entry, but are required by the routine to store
*  elements of U because of fill-in resulting from the row interchanges.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP,
     $                   JU, K2, KM, KV, NB, NW
      DOUBLE PRECISION   TEMP
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   WORK13( LDWORK, NBMAX ),
     $                   WORK31( LDWORK, NBMAX )
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX, ILAENV
      EXTERNAL           IDAMAX, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGBTF2, DGEMM, DGER, DLASWP, DSCAL,
     $                   DSWAP, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     KV is the number of superdiagonals in the factor U, allowing for
*     fill-in
*
      KV = KU + KL
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment
*
      NB = ILAENV( 1, 'DGBTRF', ' ', M, N, KL, KU )
*
*     The block size must not exceed the limit set by the size of the
*     local arrays WORK13 and WORK31.
*
      NB = MIN( NB, NBMAX )
*
      IF( NB.LE.1 .OR. NB.GT.KL ) THEN
*
*        Use unblocked code
*
         CALL DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
      ELSE
*
*        Use blocked code
*
*        Zero the superdiagonal elements of the work array WORK13
*
         DO 20 J = 1, NB
            DO 10 I = 1, J - 1
               WORK13( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
*
*        Zero the subdiagonal elements of the work array WORK31
*
         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE
*
*        Gaussian elimination with partial pivoting
*
*        Set fill-in elements in columns KU+2 to KV to zero
*
         DO 60 J = KU + 2, MIN( KV, N )
            DO 50 I = KV - J + 2, KL
               AB( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE
*
*        JU is the index of the last column affected by the current
*        stage of the factorization
*
         JU = 1
*
         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )
*
*           The active part of the matrix is partitioned
*
*              A11   A12   A13
*              A21   A22   A23
*              A31   A32   A33
*
*           Here A11, A21 and A31 denote the current block of JB columns
*           which is about to be factorized. The number of rows in the
*           partitioning are JB, I2, I3 respectively, and the numbers
*           of columns are JB, J2, J3. The superdiagonal elements of A13
*           and the subdiagonal elements of A31 lie outside the band.
*
            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )
*
*           J2 and J3 are computed after JU has been updated.
*
*           Factorize the current block of JB columns
*
            DO 80 JJ = J, J + JB - 1
*
*              Set fill-in elements in column JJ+KV to zero
*
               IF( JJ+KV.LE.N ) THEN
                  DO 70 I = 1, KL
                     AB( I, JJ+KV ) = ZERO
   70             CONTINUE
               END IF
*
*              Find pivot and test for singularity. KM is the number of
*              subdiagonal elements in the current column.
*
               KM = MIN( KL, M-JJ )
               JP = IDAMAX( KM+1, AB( KV+1, JJ ), 1 )
               IPIV( JJ ) = JP + JJ - J
               IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                  IF( JP.NE.1 ) THEN
*
*                    Apply interchange to columns J to J+JB-1
*
                     IF( JP+JJ-1.LT.J+KL ) THEN
*
                        CALL DSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1,
     $                              AB( KV+JP+JJ-J, J ), LDAB-1 )
                     ELSE
*
*                       The interchange affects columns J to JJ-1 of A31
*                       which are stored in the work array WORK31
*
                        CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                              WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                        CALL DSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1,
     $                              AB( KV+JP, JJ ), LDAB-1 )
                     END IF
                  END IF
*
*                 Compute multipliers
*
                  CALL DSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ),
     $                        1 )
*
*                 Update trailing submatrix within the band and within
*                 the current block. JM is the index of the last column
*                 which needs to be updated.
*
                  JM = MIN( JU, J+JB-1 )
                  IF( JM.GT.JJ )
     $               CALL DGER( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1,
     $                          AB( KV, JJ+1 ), LDAB-1,
     $                          AB( KV+1, JJ+1 ), LDAB-1 )
               ELSE
*
*                 If pivot is zero, set INFO to the index of the pivot
*                 unless a zero pivot has already been found.
*
                  IF( INFO.EQ.0 )
     $               INFO = JJ
               END IF
*
*              Copy current column of A31 into the work array WORK31
*
               NW = MIN( JJ-J+1, I3 )
               IF( NW.GT.0 )
     $            CALL DCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1,
     $                        WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            IF( J+JB.LE.N ) THEN
*
*              Apply the row interchanges to the other blocks.
*
               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )
*
*              Use DLASWP to apply the row interchanges to A12, A22, and
*              A32.
*
               CALL DLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB,
     $                      IPIV( J ), 1 )
*
*              Adjust the pivot indices.
*
               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE
*
*              Apply the row interchanges to A13, A23, and A33
*              columnwise.
*
               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO 100 II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     IF( IP.NE.II ) THEN
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     END IF
  100             CONTINUE
  110          CONTINUE
*
*              Update the relevant part of the trailing submatrix
*
               IF( J2.GT.0 ) THEN
*
*                 Update A12
*
                  CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     $                        JB, J2, ONE, AB( KV+1, J ), LDAB-1,
     $                        AB( KV+1-JB, J+JB ), LDAB-1 )
*
                  IF( I2.GT.0 ) THEN
*
*                    Update A22
*
                     CALL DGEMM( 'No transpose', 'No transpose', I2, J2,
     $                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     $                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,
     $                           AB( KV+1, J+JB ), LDAB-1 )
                  END IF
*
                  IF( I3.GT.0 ) THEN
*
*                    Update A32
*
                     CALL DGEMM( 'No transpose', 'No transpose', I3, J2,
     $                           JB, -ONE, WORK31, LDWORK,
     $                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,
     $                           AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  END IF
               END IF
*
               IF( J3.GT.0 ) THEN
*
*                 Copy the lower triangle of A13 into the work array
*                 WORK13
*
                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE
*
*                 Update A13 in the work array
*
                  CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     $                        JB, J3, ONE, AB( KV+1, J ), LDAB-1,
     $                        WORK13, LDWORK )
*
                  IF( I2.GT.0 ) THEN
*
*                    Update A23
*
                     CALL DGEMM( 'No transpose', 'No transpose', I2, J3,
     $                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     $                           WORK13, LDWORK, ONE, AB( 1+JB, J+KV ),
     $                           LDAB-1 )
                  END IF
*
                  IF( I3.GT.0 ) THEN
*
*                    Update A33
*
                     CALL DGEMM( 'No transpose', 'No transpose', I3, J3,
     $                           JB, -ONE, WORK31, LDWORK, WORK13,
     $                           LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                  END IF
*
*                 Copy the lower triangle of A13 back into place
*
                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE
*
*              Adjust the pivot indices.
*
               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            END IF
*
*           Partially undo the interchanges in the current block to
*           restore the upper triangular form of A31 and copy the upper
*           triangle of A31 back into place
*
            DO 170 JJ = J + JB - 1, J, -1
               JP = IPIV( JJ ) - JJ + 1
               IF( JP.NE.1 ) THEN
*
*                 Apply interchange to columns J to JJ-1
*
                  IF( JP+JJ-1.LT.J+KL ) THEN
*
*                    The interchange does not affect A31
*
                     CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                           AB( KV+JP+JJ-J, J ), LDAB-1 )
                  ELSE
*
*                    The interchange does affect A31
*
                     CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                           WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                  END IF
               END IF
*
*              Copy the current column of A31 back into place
*
               NW = MIN( I3, JJ-J+1 )
               IF( NW.GT.0 )
     $            CALL DCOPY( NW, WORK31( 1, JJ-J+1 ), 1,
     $                        AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
  180    CONTINUE
      END IF
*
      RETURN
*
*     End of DGBTRF
*
      END
      SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
     $                   INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGBTRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general band matrix A using the LU factorization computed
*  by DGBTRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations.
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
*          Details of the LU factorization of the band matrix A, as
*          computed by DGBTRF.  U is stored as an upper triangular band
*          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
*          the multipliers used during the factorization are stored in
*          rows KL+KU+2 to 2*KL+KU+1.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices; for 1 <= i <= N, row i of the matrix was
*          interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LNOTI, NOTRAN
      INTEGER            I, J, KD, L, LM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER, DSWAP, DTBSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.( 2*KL+KU+1 ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      KD = KU + KL + 1
      LNOTI = KL.GT.0
*
      IF( NOTRAN ) THEN
*
*        Solve  A*X = B.
*
*        Solve L*X = B, overwriting B with X.
*
*        L is represented as a product of permutations and unit lower
*        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
*        where each transformation L(i) is a rank-one modification of
*        the identity matrix.
*
         IF( LNOTI ) THEN
            DO 10 J = 1, N - 1
               LM = MIN( KL, N-J )
               L = IPIV( J )
               IF( L.NE.J )
     $            CALL DSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
               CALL DGER( LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ),
     $                    LDB, B( J+1, 1 ), LDB )
   10       CONTINUE
         END IF
*
         DO 20 I = 1, NRHS
*
*           Solve U*X = B, overwriting B with X.
*
            CALL DTBSV( 'Upper', 'No transpose', 'Non-unit', N, KL+KU,
     $                  AB, LDAB, B( 1, I ), 1 )
   20    CONTINUE
*
      ELSE
*
*        Solve A'*X = B.
*
         DO 30 I = 1, NRHS
*
*           Solve U'*X = B, overwriting B with X.
*
            CALL DTBSV( 'Upper', 'Transpose', 'Non-unit', N, KL+KU, AB,
     $                  LDAB, B( 1, I ), 1 )
   30    CONTINUE
*
*        Solve L'*X = B, overwriting B with X.
*
         IF( LNOTI ) THEN
            DO 40 J = N - 1, 1, -1
               LM = MIN( KL, N-J )
               CALL DGEMV( 'Transpose', LM, NRHS, -ONE, B( J+1, 1 ),
     $                     LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB )
               L = IPIV( J )
               IF( L.NE.J )
     $            CALL DSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
   40       CONTINUE
         END IF
      END IF
      RETURN
*
*     End of DGBTRS
*
      END
      SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
*    .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER INCX,INCY,LDA,M,N
*    ..
*    .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*    ..
*
* Purpose
* =======
*
* DGER performs the rank 1 operation
*
*    A := alpha*x*y' + A,
*
* where alpha is a scalar, x is an m element vector, y is an n element
* vector and A is an m by n matrix.
*
* Arguments
* ==========
*
* M       - INTEGER.
*             On entry, M specifies the number of rows of the matrix A.
*             M must be at least zero.
*             Unchanged on exit.
*
* N       - INTEGER.
*             On entry, N specifies the number of columns of the matrix A.
*             N must be at least zero.
*             Unchanged on exit.
*
* ALPHA - DOUBLE PRECISION.
*             On entry, ALPHA specifies the scalar alpha.
*             Unchanged on exit.
*
* X       - DOUBLE PRECISION array of dimension at least
*             ( 1 + ( m - 1 )*abs( INCX ) ).
*             Before entry, the incremented array X must contain the m
*             element vector x.
*             Unchanged on exit.
*
* INCX - INTEGER.
*             On entry, INCX specifies the increment for the elements of
*             X. INCX must not be zero.
*             Unchanged on exit.
*
* Y       - DOUBLE PRECISION array of dimension at least
*             ( 1 + ( n - 1 )*abs( INCY ) ).
*             Before entry, the incremented array Y must contain the n
*             element vector y.
*             Unchanged on exit.
*
* INCY - INTEGER.
*             On entry, INCY specifies the increment for the elements of
*             Y. INCY must not be zero.
*             Unchanged on exit.
*
* A       - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*             Before entry, the leading m by n part of the array A must
*             contain the matrix of coefficients. On exit, A is
*             overwritten by the updated matrix.
*
* LDA - INTEGER.
*             On entry, LDA specifies the first dimension of A as declared
*             in the calling (sub) program. LDA must be at least
*             max( 1, m ).
*             Unchanged on exit.
*
*
* Level 2 Blas routine.
*
* -- Written on 22-October-1986.
*    Jack Dongarra, Argonne National Lab.
*    Jeremy Du Croz, Nag Central Office.
*    Sven Hammarling, Nag Central Office.
*    Richard Hanson, Sandia National Labs.
*
*
*    .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
*    ..
*    .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JY,KX
*    ..
*    .. External Subroutines ..
      EXTERNAL XERBLA
*    ..
*    .. Intrinsic Functions ..
      INTRINSIC MAX
*    ..
*
*    Test the input parameters.
*
      INFO = 0
      IF (M.LT.0) THEN
         INFO = 1
      ELSE IF (N.LT.0) THEN
         INFO = 2
      ELSE IF (INCX.EQ.0) THEN
         INFO = 5
      ELSE IF (INCY.EQ.0) THEN
         INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = 9
      END IF
      IF (INFO.NE.0) THEN
         CALL XERBLA('DGER ',INFO)
         RETURN
      END IF
*
*    Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
*
*    Start the operations. In this version the elements of A are
*    accessed sequentially with one pass through A.
*
      IF (INCY.GT.0) THEN
         JY = 1
      ELSE
         JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN
         DO 20 J = 1,N
            IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  DO 10 I = 1,M
                        A(I,J) = A(I,J) + X(I)*TEMP
10             CONTINUE
            END IF
            JY = JY + INCY
20    CONTINUE
      ELSE
         IF (INCX.GT.0) THEN
            KX = 1
         ELSE
            KX = 1 - (M-1)*INCX
         END IF
         DO 40 J = 1,N
            IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  IX = KX
                  DO 30 I = 1,M
                        A(I,J) = A(I,J) + X(IX)*TEMP
                        IX = IX + INCX
30             CONTINUE
            END IF
            JY = JY + INCY
40    CONTINUE
      END IF
*
      RETURN
*
*    End of DGER .
*
      END
      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
*    .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*    ..
*    .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*    ..
*
* Purpose
* =======
*
*    interchanges two vectors.
*    uses unrolled loops for increments equal one.
*    jack dongarra, linpack, 3/11/78.
*    modified 12/3/93, array(1) declarations changed to array(*)
*
*
*    .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
*    ..
*    .. Intrinsic Functions ..
      INTRINSIC MOD
*    ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*       code for unequal increments or equal increments not equal
*       to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         DTEMP = DX(IX)
         DX(IX) = DY(IY)
         DY(IY) = DTEMP
         IX = IX + INCX
         IY = IY + INCY
10    CONTINUE
      RETURN
*
*       code for both increments equal to 1
*
*
*       clean-up loop
*
20    M = MOD(N,3)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
         DTEMP = DX(I)
         DX(I) = DY(I)
         DY(I) = DTEMP
30    CONTINUE
      IF (N.LT.3) RETURN
40    MP1 = M + 1
      DO 50 I = MP1,N,3
         DTEMP = DX(I)
         DX(I) = DY(I)
         DY(I) = DTEMP
         DTEMP = DX(I+1)
         DX(I+1) = DY(I+1)
         DY(I+1) = DTEMP
         DTEMP = DX(I+2)
         DX(I+2) = DY(I+2)
         DY(I+2) = DTEMP
50    CONTINUE
      RETURN
      END
      SUBROUTINE XERBLA(SRNAME,INFO)
*
* -- LAPACK auxiliary routine (preliminary version) --
*    Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*    November 2006
*
*    .. Scalar Arguments ..
      INTEGER INFO
      CHARACTER*6 SRNAME
*    ..
*
* Purpose
* =======
*
* XERBLA is an error handler for the LAPACK routines.
* It is called by an LAPACK routine if an input parameter has an
* invalid value. A message is printed and execution stops.
*
* Installers may consider modifying the STOP statement in order to
* call system-specific exception-handling facilities.
*
* Arguments
* =========
*
* SRNAME (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
* INFO (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
*
      WRITE (*,FMT=9999) SRNAME,INFO
*
      STOP
*
9999  FORMAT (' ** On entry to ',A6,' parameter number ',I2,' had ',
     +       'an illegal value')
*
*    End of XERBLA
*
      END
      SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*    .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*    ..
*    .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
*    ..
*
* Purpose
* =======
*
* DTRSM solves one of the matrix equations
*
*    op( A )*X = alpha*B, or X*op( A ) = alpha*B,
*
* where alpha is a scalar, X and B are m by n matrices, A is a unit, or
* non-unit, upper or lower triangular matrix and op( A ) is one of
*
*    op( A ) = A or op( A ) = A'.
*
* The matrix X is overwritten on B.
*
* Arguments
* ==========
*
* SIDE - CHARACTER*1.
*             On entry, SIDE specifies whether op( A ) appears on the left
*             or right of X as follows:
*
*             SIDE = 'L' or 'l' op( A )*X = alpha*B.
*
*             SIDE = 'R' or 'r' X*op( A ) = alpha*B.
*
*             Unchanged on exit.
*
* UPLO - CHARACTER*1.
*             On entry, UPLO specifies whether the matrix A is an upper or
*             lower triangular matrix as follows:
*
*             UPLO = 'U' or 'u' A is an upper triangular matrix.
*
*             UPLO = 'L' or 'l' A is a lower triangular matrix.
*
*             Unchanged on exit.
*
* TRANSA - CHARACTER*1.
*             On entry, TRANSA specifies the form of op( A ) to be used in
*             the matrix multiplication as follows:
*
*             TRANSA = 'N' or 'n' op( A ) = A.
*
*             TRANSA = 'T' or 't' op( A ) = A'.
*
*             TRANSA = 'C' or 'c' op( A ) = A'.
*
*             Unchanged on exit.
*
* DIAG - CHARACTER*1.
*             On entry, DIAG specifies whether or not A is unit triangular
*             as follows:
*
*             DIAG = 'U' or 'u' A is assumed to be unit triangular.
*
*             DIAG = 'N' or 'n' A is not assumed to be unit
*                                     triangular.
*
*             Unchanged on exit.
*
* M       - INTEGER.
*             On entry, M specifies the number of rows of B. M must be at
*             least zero.
*             Unchanged on exit.
*
* N       - INTEGER.
*             On entry, N specifies the number of columns of B. N must be
*             at least zero.
*             Unchanged on exit.
*
* ALPHA - DOUBLE PRECISION.
*             On entry, ALPHA specifies the scalar alpha. When alpha is
*             zero then A is not referenced and B need not be set before
*             entry.
*             Unchanged on exit.
*
* A       - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*             when SIDE = 'L' or 'l' and is n when SIDE = 'R' or 'r'.
*             Before entry with UPLO = 'U' or 'u', the leading k by k
*             upper triangular part of the array A must contain the upper
*             triangular matrix and the strictly lower triangular part of
*             A is not referenced.
*             Before entry with UPLO = 'L' or 'l', the leading k by k
*             lower triangular part of the array A must contain the lower
*             triangular matrix and the strictly upper triangular part of
*             A is not referenced.
*             Note that when DIAG = 'U' or 'u', the diagonal elements of
*             A are not referenced either, but are assumed to be unity.
*             Unchanged on exit.
*
* LDA - INTEGER.
*             On entry, LDA specifies the first dimension of A as declared
*             in the calling (sub) program. When SIDE = 'L' or 'l' then
*             LDA must be at least max( 1, m ), when SIDE = 'R' or 'r'
*             then LDA must be at least max( 1, n ).
*             Unchanged on exit.
*
* B       - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*             Before entry, the leading m by n part of the array B must
*             contain the right-hand side matrix B, and on exit is
*             overwritten by the solution matrix X.
*
* LDB - INTEGER.
*             On entry, LDB specifies the first dimension of B as declared
*             in the calling (sub) program. LDB must be at least
*             max( 1, m ).
*             Unchanged on exit.
*
*
* Level 3 Blas routine.
*
*
* -- Written on 8-February-1989.
*    Jack Dongarra, Argonne National Laboratory.
*    Iain Duff, AERE Harwell.
*    Jeremy Du Croz, Numerical Algorithms Group Ltd.
*    Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*    .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*    ..
*    .. External Subroutines ..
      EXTERNAL XERBLA
*    ..
*    .. Intrinsic Functions ..
      INTRINSIC MAX
*    ..
*    .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
*    ..
*    .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*    ..
*
*    Test the input parameters.
*
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
         INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
         INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +       (.NOT.LSAME(TRANSA,'T')) .AND.
     +       (.NOT.LSAME(TRANSA,'C'))) THEN
         INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
         INFO = 4
      ELSE IF (M.LT.0) THEN
         INFO = 5
      ELSE IF (N.LT.0) THEN
         INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
         INFO = 11
      END IF
      IF (INFO.NE.0) THEN
         CALL XERBLA('DTRSM ',INFO)
         RETURN
      END IF
*
*    Quick return if possible.
*
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
*
*    And when alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
         DO 20 J = 1,N
            DO 10 I = 1,M
                  B(I,J) = ZERO
10       CONTINUE
20    CONTINUE
         RETURN
      END IF
*
*    Start the operations.
*
      IF (LSIDE) THEN
         IF (LSAME(TRANSA,'N')) THEN
*
*             Form B := alpha*inv( A )*B.
*
            IF (UPPER) THEN
                  DO 60 J = 1,N
                        IF (ALPHA.NE.ONE) THEN
                              DO 30 I = 1,M
                                 B(I,J) = ALPHA*B(I,J)
30                         CONTINUE
                        END IF
                        DO 50 K = M,1,-1
                              IF (B(K,J).NE.ZERO) THEN
                                 IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                                 DO 40 I = 1,K - 1
                                    B(I,J) = B(I,J) - B(K,J)*A(I,K)
40                            CONTINUE
                              END IF
50                   CONTINUE
60             CONTINUE
            ELSE
                  DO 100 J = 1,N
                        IF (ALPHA.NE.ONE) THEN
                              DO 70 I = 1,M
                                 B(I,J) = ALPHA*B(I,J)
70                         CONTINUE
                        END IF
                        DO 90 K = 1,M
                              IF (B(K,J).NE.ZERO) THEN
                                 IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                                 DO 80 I = K + 1,M
                                    B(I,J) = B(I,J) - B(K,J)*A(I,K)
80                            CONTINUE
                              END IF
90                   CONTINUE
100             CONTINUE
            END IF
         ELSE
*
*             Form B := alpha*inv( A' )*B.
*
            IF (UPPER) THEN
                  DO 130 J = 1,N
                        DO 120 I = 1,M
                              TEMP = ALPHA*B(I,J)
                              DO 110 K = 1,I - 1
                                 TEMP = TEMP - A(K,I)*B(K,J)
110                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/A(I,I)
                              B(I,J) = TEMP
120                   CONTINUE
130             CONTINUE
            ELSE
                  DO 160 J = 1,N
                        DO 150 I = M,1,-1
                              TEMP = ALPHA*B(I,J)
                              DO 140 K = I + 1,M
                                 TEMP = TEMP - A(K,I)*B(K,J)
140                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/A(I,I)
                              B(I,J) = TEMP
150                   CONTINUE
160             CONTINUE
            END IF
         END IF
      ELSE
         IF (LSAME(TRANSA,'N')) THEN
*
*             Form B := alpha*B*inv( A ).
*
            IF (UPPER) THEN
                  DO 210 J = 1,N
                        IF (ALPHA.NE.ONE) THEN
                              DO 170 I = 1,M
                                 B(I,J) = ALPHA*B(I,J)
170                         CONTINUE
                        END IF
                        DO 190 K = 1,J - 1
                              IF (A(K,J).NE.ZERO) THEN
                                 DO 180 I = 1,M
                                    B(I,J) = B(I,J) - A(K,J)*B(I,K)
180                            CONTINUE
                              END IF
190                   CONTINUE
                        IF (NOUNIT) THEN
                              TEMP = ONE/A(J,J)
                              DO 200 I = 1,M
                                 B(I,J) = TEMP*B(I,J)
200                         CONTINUE
                        END IF
210             CONTINUE
            ELSE
                  DO 260 J = N,1,-1
                        IF (ALPHA.NE.ONE) THEN
                              DO 220 I = 1,M
                                 B(I,J) = ALPHA*B(I,J)
220                         CONTINUE
                        END IF
                        DO 240 K = J + 1,N
                              IF (A(K,J).NE.ZERO) THEN
                                 DO 230 I = 1,M
                                    B(I,J) = B(I,J) - A(K,J)*B(I,K)
230                            CONTINUE
                              END IF
240                   CONTINUE
                        IF (NOUNIT) THEN
                              TEMP = ONE/A(J,J)
                              DO 250 I = 1,M
                                 B(I,J) = TEMP*B(I,J)
250                         CONTINUE
                        END IF
260             CONTINUE
            END IF
         ELSE
*
*             Form B := alpha*B*inv( A' ).
*
            IF (UPPER) THEN
                  DO 310 K = N,1,-1
                        IF (NOUNIT) THEN
                              TEMP = ONE/A(K,K)
                              DO 270 I = 1,M
                                 B(I,K) = TEMP*B(I,K)
270                         CONTINUE
                        END IF
                        DO 290 J = 1,K - 1
                              IF (A(J,K).NE.ZERO) THEN
                                 TEMP = A(J,K)
                                 DO 280 I = 1,M
                                    B(I,J) = B(I,J) - TEMP*B(I,K)
280                            CONTINUE
                              END IF
290                   CONTINUE
                        IF (ALPHA.NE.ONE) THEN
                              DO 300 I = 1,M
                                 B(I,K) = ALPHA*B(I,K)
300                         CONTINUE
                        END IF
310             CONTINUE
            ELSE
                  DO 360 K = 1,N
                        IF (NOUNIT) THEN
                              TEMP = ONE/A(K,K)
                              DO 320 I = 1,M
                                 B(I,K) = TEMP*B(I,K)
320                         CONTINUE
                        END IF
                        DO 340 J = K + 1,N
                              IF (A(J,K).NE.ZERO) THEN
                                 TEMP = A(J,K)
                                 DO 330 I = 1,M
                                    B(I,J) = B(I,J) - TEMP*B(I,K)
330                            CONTINUE
                              END IF
340                   CONTINUE
                        IF (ALPHA.NE.ONE) THEN
                              DO 350 I = 1,M
                                 B(I,K) = ALPHA*B(I,K)
350                         CONTINUE
                        END IF
360             CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*    End of DTRSM .
*
      END
C      SUBROUTINE DSCAL(N,DA,DX,INCX)
*    .. Scalar Arguments ..
C      DOUBLE PRECISION DA
C      INTEGER INCX,N
*    ..
*    .. Array Arguments ..
C      DOUBLE PRECISION DX(*)
*    ..
*
* Purpose
* =======
**
*    scales a vector by a constant.
*    uses unrolled loops for increment equal to one.
*    jack dongarra, linpack, 3/11/78.
*    modified 3/93 to return if incx .le. 0.
*    modified 12/3/93, array(1) declarations changed to array(*)
*
*
*    .. Local Scalars ..
C      INTEGER I,M,MP1,NINCX
*    ..
*    .. Intrinsic Functions ..
C      INTRINSIC MOD
*    ..
C      IF (N.LE.0 .OR. INCX.LE.0) RETURN
C      IF (INCX.EQ.1) GO TO 20
*
*       code for increment not equal to 1
*
C      NINCX = N*INCX
C      DO 10 I = 1,NINCX,INCX
C         DX(I) = DA*DX(I)
C10    CONTINUE
C      RETURN
*
*       code for increment equal to 1
*
*
*       clean-up loop
*
C20    M = MOD(N,5)
C      IF (M.EQ.0) GO TO 40
C      DO 30 I = 1,M
C         DX(I) = DA*DX(I)
C30    CONTINUE
C      IF (N.LT.5) RETURN
C40    MP1 = M + 1
C      DO 50 I = MP1,N,5
C         DX(I) = DA*DX(I)
C         DX(I+1) = DA*DX(I+1)
C         DX(I+2) = DA*DX(I+2)
C         DX(I+3) = DA*DX(I+3)
C         DX(I+4) = DA*DX(I+4)
C50    CONTINUE
C      RETURN
C      END
C      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
*    .. Scalar Arguments ..
C      INTEGER INCX,INCY,N
*    ..
*    .. Array Arguments ..
C      DOUBLE PRECISION DX(*),DY(*)
*    ..
*
* Purpose
* =======
*
*    copies a vector, x, to a vector, y.
*    uses unrolled loops for increments equal to one.
*    jack dongarra, linpack, 3/11/78.
*    modified 12/3/93, array(1) declarations changed to array(*)
*
*
*    .. Local Scalars ..
C      INTEGER I,IX,IY,M,MP1
*    ..
*    .. Intrinsic Functions ..
C      INTRINSIC MOD
*    ..
C      IF (N.LE.0) RETURN
C      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*       code for unequal increments or equal increments
*          not equal to 1
*
C      IX = 1
C      IY = 1
C      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
C      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
C      DO 10 I = 1,N
C         DY(IY) = DX(IX)
C         IX = IX + INCX
C         IY = IY + INCY
C10    CONTINUE
C      RETURN
*
*       code for both increments equal to 1
*
*
*       clean-up loop
*
C20    M = MOD(N,7)
C      IF (M.EQ.0) GO TO 40
C      DO 30 I = 1,M
C         DY(I) = DX(I)
C30    CONTINUE
C      IF (N.LT.7) RETURN
C40    MP1 = M + 1
C      DO 50 I = MP1,N,7
C         DY(I) = DX(I)
C         DY(I+1) = DX(I+1)
C         DY(I+2) = DX(I+2)
C         DY(I+3) = DX(I+3)
C         DY(I+4) = DX(I+4)
C         DY(I+5) = DX(I+5)
C         DY(I+6) = DX(I+6)
C50    CONTINUE
C      RETURN
C      END
      SUBROUTINE DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
*    .. Scalar Arguments ..
      INTEGER INCX,K,LDA,N
      CHARACTER DIAG,TRANS,UPLO
*    ..
*    .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*)
*    ..
*
* Purpose
* =======
*
* DTBSV solves one of the systems of equations
*
*    A*x = b, or A'*x = b,
*
* where b and x are n element vectors and A is an n by n unit, or
* non-unit, upper or lower triangular band matrix, with ( k + 1 )
* diagonals.
*
* No test for singularity or near-singularity is included in this
* routine. Such tests must be performed before calling this routine.
*
* Arguments
* ==========
*
* UPLO - CHARACTER*1.
*             On entry, UPLO specifies whether the matrix is an upper or
*             lower triangular matrix as follows:
*
*             UPLO = 'U' or 'u' A is an upper triangular matrix.
*
*             UPLO = 'L' or 'l' A is a lower triangular matrix.
*
*             Unchanged on exit.
*
* TRANS - CHARACTER*1.
*             On entry, TRANS specifies the equations to be solved as
*             follows:
*
*             TRANS = 'N' or 'n' A*x = b.
*
*             TRANS = 'T' or 't' A'*x = b.
*
*             TRANS = 'C' or 'c' A'*x = b.
*
*             Unchanged on exit.
*
* DIAG - CHARACTER*1.
*             On entry, DIAG specifies whether or not A is unit
*             triangular as follows:
*
*             DIAG = 'U' or 'u' A is assumed to be unit triangular.
*
*             DIAG = 'N' or 'n' A is not assumed to be unit
*                                     triangular.
*
*             Unchanged on exit.
*
* N       - INTEGER.
*             On entry, N specifies the order of the matrix A.
*             N must be at least zero.
*             Unchanged on exit.
*
* K       - INTEGER.
*             On entry with UPLO = 'U' or 'u', K specifies the number of
*             super-diagonals of the matrix A.
*             On entry with UPLO = 'L' or 'l', K specifies the number of
*             sub-diagonals of the matrix A.
*             K must satisfy 0 .le. K.
*             Unchanged on exit.
*
* A       - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*             Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
*             by n part of the array A must contain the upper triangular
*             band part of the matrix of coefficients, supplied column by
*             column, with the leading diagonal of the matrix in row
*             ( k + 1 ) of the array, the first super-diagonal starting at
*             position 2 in row k, and so on. The top left k by k triangle
*             of the array A is not referenced.
*             The following program segment will transfer an upper
*             triangular band matrix from conventional full matrix storage
*             to band storage:
*
*                   DO 20, J = 1, N
*                      M = K + 1 - J
*                      DO 10, I = MAX( 1, J - K ), J
*                         A( M + I, J ) = matrix( I, J )
*             10 CONTINUE
*             20 CONTINUE
*
*             Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
*             by n part of the array A must contain the lower triangular
*             band part of the matrix of coefficients, supplied column by
*             column, with the leading diagonal of the matrix in row 1 of
*             the array, the first sub-diagonal starting at position 1 in
*             row 2, and so on. The bottom right k by k triangle of the
*             array A is not referenced.
*             The following program segment will transfer a lower
*             triangular band matrix from conventional full matrix storage
*             to band storage:
*
*                   DO 20, J = 1, N
*                      M = 1 - J
*                      DO 10, I = J, MIN( N, J + K )
*                         A( M + I, J ) = matrix( I, J )
*             10 CONTINUE
*             20 CONTINUE
*
*             Note that when DIAG = 'U' or 'u' the elements of the array A
*             corresponding to the diagonal elements of the matrix are not
*             referenced, but are assumed to be unity.
*             Unchanged on exit.
*
* LDA - INTEGER.
*             On entry, LDA specifies the first dimension of A as declared
*             in the calling (sub) program. LDA must be at least
*             ( k + 1 ).
*             Unchanged on exit.
*
* X       - DOUBLE PRECISION array of dimension at least
*             ( 1 + ( n - 1 )*abs( INCX ) ).
*             Before entry, the incremented array X must contain the n
*             element right-hand side vector b. On exit, X is overwritten
*             with the solution vector x.
*
* INCX - INTEGER.
*             On entry, INCX specifies the increment for the elements of
*             X. INCX must not be zero.
*             Unchanged on exit.
*
*
* Level 2 Blas routine.
*
* -- Written on 22-October-1986.
*    Jack Dongarra, Argonne National Lab.
*    Jeremy Du Croz, Nag Central Office.
*    Sven Hammarling, Nag Central Office.
*    Richard Hanson, Sandia National Labs.
*
*
*    .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
*    ..
*    .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KPLUS1,KX,L
      LOGICAL NOUNIT
*    ..
*    .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*    ..
*    .. External Subroutines ..
      EXTERNAL XERBLA
*    ..
*    .. Intrinsic Functions ..
      INTRINSIC MAX,MIN
*    ..
*
*    Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
         INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +       .NOT.LSAME(TRANS,'C')) THEN
         INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
         INFO = 3
      ELSE IF (N.LT.0) THEN
         INFO = 4
      ELSE IF (K.LT.0) THEN
         INFO = 5
      ELSE IF (LDA.LT. (K+1)) THEN
         INFO = 7
      ELSE IF (INCX.EQ.0) THEN
         INFO = 9
      END IF
      IF (INFO.NE.0) THEN
         CALL XERBLA('DTBSV ',INFO)
         RETURN
      END IF
*
*    Quick return if possible.
*
      IF (N.EQ.0) RETURN
*
      NOUNIT = LSAME(DIAG,'N')
*
*    Set up the start point in X if the increment is not unity. This
*    will be ( N - 1 )*INCX too small for descending loops.
*
      IF (INCX.LE.0) THEN
         KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
         KX = 1
      END IF
*
*    Start the operations. In this version the elements of A are
*    accessed by sequentially with one pass through A.
*
      IF (LSAME(TRANS,'N')) THEN
*
*       Form x := inv( A )*x.
*
         IF (LSAME(UPLO,'U')) THEN
            KPLUS1 = K + 1
            IF (INCX.EQ.1) THEN
                  DO 20 J = N,1,-1
                        IF (X(J).NE.ZERO) THEN
                              L = KPLUS1 - J
                              IF (NOUNIT) X(J) = X(J)/A(KPLUS1,J)
                              TEMP = X(J)
                              DO 10 I = J - 1,MAX(1,J-K),-1
                                 X(I) = X(I) - TEMP*A(L+I,J)
10                         CONTINUE
                        END IF
20             CONTINUE
            ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 40 J = N,1,-1
                        KX = KX - INCX
                        IF (X(JX).NE.ZERO) THEN
                              IX = KX
                              L = KPLUS1 - J
                              IF (NOUNIT) X(JX) = X(JX)/A(KPLUS1,J)
                              TEMP = X(JX)
                              DO 30 I = J - 1,MAX(1,J-K),-1
                                 X(IX) = X(IX) - TEMP*A(L+I,J)
                                 IX = IX - INCX
30                         CONTINUE
                        END IF
                        JX = JX - INCX
40             CONTINUE
            END IF
         ELSE
            IF (INCX.EQ.1) THEN
                  DO 60 J = 1,N
                        IF (X(J).NE.ZERO) THEN
                              L = 1 - J
                              IF (NOUNIT) X(J) = X(J)/A(1,J)
                              TEMP = X(J)
                              DO 50 I = J + 1,MIN(N,J+K)
                                 X(I) = X(I) - TEMP*A(L+I,J)
50                         CONTINUE
                        END IF
60             CONTINUE
            ELSE
                  JX = KX
                  DO 80 J = 1,N
                        KX = KX + INCX
                        IF (X(JX).NE.ZERO) THEN
                              IX = KX
                              L = 1 - J
                              IF (NOUNIT) X(JX) = X(JX)/A(1,J)
                              TEMP = X(JX)
                              DO 70 I = J + 1,MIN(N,J+K)
                                 X(IX) = X(IX) - TEMP*A(L+I,J)
                                 IX = IX + INCX
70                         CONTINUE
                        END IF
                        JX = JX + INCX
80             CONTINUE
            END IF
         END IF
      ELSE
*
*       Form x := inv( A')*x.
*
         IF (LSAME(UPLO,'U')) THEN
            KPLUS1 = K + 1
            IF (INCX.EQ.1) THEN
                  DO 100 J = 1,N
                        TEMP = X(J)
                        L = KPLUS1 - J
                        DO 90 I = MAX(1,J-K),J - 1
                              TEMP = TEMP - A(L+I,J)*X(I)
90                   CONTINUE
                        IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                        X(J) = TEMP
100             CONTINUE
            ELSE
                  JX = KX
                  DO 120 J = 1,N
                        TEMP = X(JX)
                        IX = KX
                        L = KPLUS1 - J
                        DO 110 I = MAX(1,J-K),J - 1
                              TEMP = TEMP - A(L+I,J)*X(IX)
                              IX = IX + INCX
110                   CONTINUE
                        IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                        X(JX) = TEMP
                        JX = JX + INCX
                        IF (J.GT.K) KX = KX + INCX
120             CONTINUE
            END IF
         ELSE
            IF (INCX.EQ.1) THEN
                  DO 140 J = N,1,-1
                        TEMP = X(J)
                        L = 1 - J
                        DO 130 I = MIN(N,J+K),J + 1,-1
                              TEMP = TEMP - A(L+I,J)*X(I)
130                   CONTINUE
                        IF (NOUNIT) TEMP = TEMP/A(1,J)
                        X(J) = TEMP
140             CONTINUE
            ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 160 J = N,1,-1
                        TEMP = X(JX)
                        IX = KX
                        L = 1 - J
                        DO 150 I = MIN(N,J+K),J + 1,-1
                              TEMP = TEMP - A(L+I,J)*X(IX)
                              IX = IX - INCX
150                   CONTINUE
                        IF (NOUNIT) TEMP = TEMP/A(1,J)
                        X(JX) = TEMP
                        JX = JX - INCX
                        IF ((N-J).GE.K) KX = KX - INCX
160             CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*    End of DTBSV .
*
      END
      SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*    .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
*    ..
*    .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*    ..
*
* Purpose
* =======
*
* DGEMV performs one of the matrix-vector operations
*
*    y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y,
*
* where alpha and beta are scalars, x and y are vectors and A is an
* m by n matrix.
*
* Arguments
* ==========
*
* TRANS - CHARACTER*1.
*             On entry, TRANS specifies the operation to be performed as
*             follows:
*
*             TRANS = 'N' or 'n' y := alpha*A*x + beta*y.
*
*             TRANS = 'T' or 't' y := alpha*A'*x + beta*y.
*
*             TRANS = 'C' or 'c' y := alpha*A'*x + beta*y.
*
*             Unchanged on exit.
*
* M       - INTEGER.
*             On entry, M specifies the number of rows of the matrix A.
*             M must be at least zero.
*             Unchanged on exit.
*
* N       - INTEGER.
*             On entry, N specifies the number of columns of the matrix A.
*             N must be at least zero.
*             Unchanged on exit.
*
* ALPHA - DOUBLE PRECISION.
*             On entry, ALPHA specifies the scalar alpha.
*             Unchanged on exit.
*
* A       - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*             Before entry, the leading m by n part of the array A must
*             contain the matrix of coefficients.
*             Unchanged on exit.
*
* LDA - INTEGER.
*             On entry, LDA specifies the first dimension of A as declared
*             in the calling (sub) program. LDA must be at least
*             max( 1, m ).
*             Unchanged on exit.
*
* X       - DOUBLE PRECISION array of DIMENSION at least
*             ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*             and at least
*             ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*             Before entry, the incremented array X must contain the
*             vector x.
*             Unchanged on exit.
*
* INCX - INTEGER.
*             On entry, INCX specifies the increment for the elements of
*             X. INCX must not be zero.
*             Unchanged on exit.
*
* BETA - DOUBLE PRECISION.
*             On entry, BETA specifies the scalar beta. When BETA is
*             supplied as zero then Y need not be set on input.
*             Unchanged on exit.
*
* Y       - DOUBLE PRECISION array of DIMENSION at least
*             ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*             and at least
*             ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*             Before entry with BETA non-zero, the incremented array Y
*             must contain the vector y. On exit, Y is overwritten by the
*             updated vector y.
*
* INCY - INTEGER.
*             On entry, INCY specifies the increment for the elements of
*             Y. INCY must not be zero.
*             Unchanged on exit.
*
*
* Level 2 Blas routine.
*
* -- Written on 22-October-1986.
*    Jack Dongarra, Argonne National Lab.
*    Jeremy Du Croz, Nag Central Office.
*    Sven Hammarling, Nag Central Office.
*    Richard Hanson, Sandia National Labs.
*
*
*    .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*    ..
*    .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
*    ..
*    .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*    ..
*    .. External Subroutines ..
      EXTERNAL XERBLA
*    ..
*    .. Intrinsic Functions ..
      INTRINSIC MAX
*    ..
*
*    Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     + .NOT.LSAME(TRANS,'C')) THEN
         INFO = 1
      ELSE IF (M.LT.0) THEN
         INFO = 2
      ELSE IF (N.LT.0) THEN
         INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = 6
      ELSE IF (INCX.EQ.0) THEN
         INFO = 8
      ELSE IF (INCY.EQ.0) THEN
         INFO = 11
      END IF
      IF (INFO.NE.0) THEN
         CALL XERBLA('DGEMV ',INFO)
         RETURN
      END IF
*
*    Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     + ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
*
*    Set LENX and LENY, the lengths of the vectors x and y, and set
*    up the start points in X and Y.
*
      IF (LSAME(TRANS,'N')) THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF (INCX.GT.0) THEN
         KX = 1
      ELSE
         KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
         KY = 1
      ELSE
         KY = 1 - (LENY-1)*INCY
      END IF
*
*    Start the operations. In this version the elements of A are
*    accessed sequentially with one pass through A.
*
*    First form y := beta*y.
*
      IF (BETA.NE.ONE) THEN
         IF (INCY.EQ.1) THEN
            IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,LENY
                        Y(I) = ZERO
10             CONTINUE
            ELSE
                  DO 20 I = 1,LENY
                        Y(I) = BETA*Y(I)
20             CONTINUE
            END IF
         ELSE
            IY = KY
            IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,LENY
                        Y(IY) = ZERO
                        IY = IY + INCY
30             CONTINUE
            ELSE
                  DO 40 I = 1,LENY
                        Y(IY) = BETA*Y(IY)
                        IY = IY + INCY
40             CONTINUE
            END IF
         END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
*
*       Form y := alpha*A*x + y.
*
         JX = KX
         IF (INCY.EQ.1) THEN
            DO 60 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                        TEMP = ALPHA*X(JX)
                        DO 50 I = 1,M
                              Y(I) = Y(I) + TEMP*A(I,J)
50                   CONTINUE
                  END IF
                  JX = JX + INCX
60       CONTINUE
         ELSE
            DO 80 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                        TEMP = ALPHA*X(JX)
                        IY = KY
                        DO 70 I = 1,M
                              Y(IY) = Y(IY) + TEMP*A(I,J)
                              IY = IY + INCY
70                   CONTINUE
                  END IF
                  JX = JX + INCX
80       CONTINUE
         END IF
      ELSE
*
*       Form y := alpha*A'*x + y.
*
         JY = KY
         IF (INCX.EQ.1) THEN
            DO 100 J = 1,N
                  TEMP = ZERO
                  DO 90 I = 1,M
                        TEMP = TEMP + A(I,J)*X(I)
90             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
100       CONTINUE
         ELSE
            DO 120 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  DO 110 I = 1,M
                        TEMP = TEMP + A(I,J)*X(IX)
                        IX = IX + INCX
110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*    End of DGEMV .
*
      END
      SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*    .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*    ..
*    .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*    ..
*
* Purpose
* =======
*
* DGEMM performs one of the matrix-matrix operations
*
*    C := alpha*op( A )*op( B ) + beta*C,
*
* where op( X ) is one of
*
*    op( X ) = X or op( X ) = X',
*
* alpha and beta are scalars, and A, B and C are matrices, with op( A )
* an m by k matrix, op( B ) a k by n matrix and C an m by n matrix.
*
* Arguments
* ==========
*
* TRANSA - CHARACTER*1.
*             On entry, TRANSA specifies the form of op( A ) to be used in
*             the matrix multiplication as follows:
*
*             TRANSA = 'N' or 'n', op( A ) = A.
*
*             TRANSA = 'T' or 't', op( A ) = A'.
*
*             TRANSA = 'C' or 'c', op( A ) = A'.
*
*             Unchanged on exit.
*
* TRANSB - CHARACTER*1.
*             On entry, TRANSB specifies the form of op( B ) to be used in
*             the matrix multiplication as follows:
*
*             TRANSB = 'N' or 'n', op( B ) = B.
*
*             TRANSB = 'T' or 't', op( B ) = B'.
*
*             TRANSB = 'C' or 'c', op( B ) = B'.
*
*             Unchanged on exit.
*
* M       - INTEGER.
*             On entry, M specifies the number of rows of the matrix
*             op( A ) and of the matrix C. M must be at least zero.
*             Unchanged on exit.
*
* N       - INTEGER.
*             On entry, N specifies the number of columns of the matrix
*             op( B ) and the number of columns of the matrix C. N must be
*             at least zero.
*             Unchanged on exit.
*
* K       - INTEGER.
*             On entry, K specifies the number of columns of the matrix
*             op( A ) and the number of rows of the matrix op( B ). K must
*             be at least zero.
*             Unchanged on exit.
*
* ALPHA - DOUBLE PRECISION.
*             On entry, ALPHA specifies the scalar alpha.
*             Unchanged on exit.
*
* A       - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*             k when TRANSA = 'N' or 'n', and is m otherwise.
*             Before entry with TRANSA = 'N' or 'n', the leading m by k
*             part of the array A must contain the matrix A, otherwise
*             the leading k by m part of the array A must contain the
*             matrix A.
*             Unchanged on exit.
*
* LDA - INTEGER.
*             On entry, LDA specifies the first dimension of A as declared
*             in the calling (sub) program. When TRANSA = 'N' or 'n' then
*             LDA must be at least max( 1, m ), otherwise LDA must be at
*             least max( 1, k ).
*             Unchanged on exit.
*
* B       - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*             n when TRANSB = 'N' or 'n', and is k otherwise.
*             Before entry with TRANSB = 'N' or 'n', the leading k by n
*             part of the array B must contain the matrix B, otherwise
*             the leading n by k part of the array B must contain the
*             matrix B.
*             Unchanged on exit.
*
* LDB - INTEGER.
*             On entry, LDB specifies the first dimension of B as declared
*             in the calling (sub) program. When TRANSB = 'N' or 'n' then
*             LDB must be at least max( 1, k ), otherwise LDB must be at
*             least max( 1, n ).
*             Unchanged on exit.
*
* BETA - DOUBLE PRECISION.
*             On entry, BETA specifies the scalar beta. When BETA is
*             supplied as zero then C need not be set on input.
*             Unchanged on exit.
*
* C       - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*             Before entry, the leading m by n part of the array C must
*             contain the matrix C, except when beta is zero, in which
*             case C need not be set on entry.
*             On exit, the array C is overwritten by the m by n matrix
*             ( alpha*op( A )*op( B ) + beta*C ).
*
* LDC - INTEGER.
*             On entry, LDC specifies the first dimension of C as declared
*             in the calling (sub) program. LDC must be at least
*             max( 1, m ).
*             Unchanged on exit.
*
*
* Level 3 Blas routine.
*
* -- Written on 8-February-1989.
*    Jack Dongarra, Argonne National Laboratory.
*    Iain Duff, AERE Harwell.
*    Jeremy Du Croz, Numerical Algorithms Group Ltd.
*    Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*    .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*    ..
*    .. External Subroutines ..
      EXTERNAL XERBLA
*    ..
*    .. Intrinsic Functions ..
      INTRINSIC MAX
*    ..
*    .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL NOTA,NOTB
*    ..
*    .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*    ..
*
*    Set NOTA and NOTB as true if A and B respectively are not
*    transposed and set NROWA, NCOLA and NROWB as the number of rows
*    and columns of A and the number of rows of B respectively.
*
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      IF (NOTA) THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF (NOTB) THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*    Test the input parameters.
*
      INFO = 0
      IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND.
     + (.NOT.LSAME(TRANSA,'T'))) THEN
         INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND.
     +       (.NOT.LSAME(TRANSB,'T'))) THEN
         INFO = 2
      ELSE IF (M.LT.0) THEN
         INFO = 3
      ELSE IF (N.LT.0) THEN
         INFO = 4
      ELSE IF (K.LT.0) THEN
         INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
         INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
         INFO = 13
      END IF
      IF (INFO.NE.0) THEN
         CALL XERBLA('DGEMM ',INFO)
         RETURN
      END IF
*
*    Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     + (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
*
*    And if alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
         IF (BETA.EQ.ZERO) THEN
            DO 20 J = 1,N
                  DO 10 I = 1,M
                        C(I,J) = ZERO
10             CONTINUE
20       CONTINUE
         ELSE
            DO 40 J = 1,N
                  DO 30 I = 1,M
                        C(I,J) = BETA*C(I,J)
30             CONTINUE
40       CONTINUE
         END IF
         RETURN
      END IF
*
*    Start the operations.
*
      IF (NOTB) THEN
         IF (NOTA) THEN
*
*             Form C := alpha*A*B + beta*C.
*
            DO 90 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                        DO 50 I = 1,M
                              C(I,J) = ZERO
50                   CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                        DO 60 I = 1,M
                              C(I,J) = BETA*C(I,J)
60                   CONTINUE
                  END IF
                  DO 80 L = 1,K
                        IF (B(L,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(L,J)
                              DO 70 I = 1,M
                                 C(I,J) = C(I,J) + TEMP*A(I,L)
70                         CONTINUE
                        END IF
80             CONTINUE
90       CONTINUE
         ELSE
*
*             Form C := alpha*A'*B + beta*C
*
            DO 120 J = 1,N
                  DO 110 I = 1,M
                        TEMP = ZERO
                        DO 100 L = 1,K
                              TEMP = TEMP + A(L,I)*B(L,J)
100                   CONTINUE
                        IF (BETA.EQ.ZERO) THEN
                              C(I,J) = ALPHA*TEMP
                        ELSE
                              C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                        END IF
110             CONTINUE
120       CONTINUE
         END IF
      ELSE
         IF (NOTA) THEN
*
*             Form C := alpha*A*B' + beta*C
*
            DO 170 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                        DO 130 I = 1,M
                              C(I,J) = ZERO
130                   CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                        DO 140 I = 1,M
                              C(I,J) = BETA*C(I,J)
140                   CONTINUE
                  END IF
                  DO 160 L = 1,K
                        IF (B(J,L).NE.ZERO) THEN
                              TEMP = ALPHA*B(J,L)
                              DO 150 I = 1,M
                                 C(I,J) = C(I,J) + TEMP*A(I,L)
150                         CONTINUE
                        END IF
160             CONTINUE
170       CONTINUE
         ELSE
*
*             Form C := alpha*A'*B' + beta*C
*
            DO 200 J = 1,N
                  DO 190 I = 1,M
                        TEMP = ZERO
                        DO 180 L = 1,K
                              TEMP = TEMP + A(L,I)*B(J,L)
180                   CONTINUE
                        IF (BETA.EQ.ZERO) THEN
                              C(I,J) = ALPHA*TEMP
                        ELSE
                              C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                        END IF
190             CONTINUE
200       CONTINUE
         END IF
      END IF
*
      RETURN
*
*    End of DGEMM .
*
      END
C      INTEGER FUNCTION IDAMAX(N,DX,INCX)
*     .. Scalar Arguments ..
C      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
C      DOUBLE PRECISION DX(*)
*     ..
*
*  Purpose
*  =======
*
*     finds the index of element having max. absolute value.
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
C      DOUBLE PRECISION DMAX
C      INTEGER I,IX
*     ..
*     .. Intrinsic Functions ..
C      INTRINSIC DABS
*     ..
C      IDAMAX = 0
C      IF (N.LT.1 .OR. INCX.LE.0) RETURN
C      IDAMAX = 1
C      IF (N.EQ.1) RETURN
C      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
C      IX = 1
C      DMAX = DABS(DX(1))
C      IX = IX + INCX
C      DO 10 I = 2,N
C          IF (DABS(DX(IX)).LE.DMAX) GO TO 5
C          IDAMAX = I
C          DMAX = DABS(DX(IX))
C    5     IX = IX + INCX
C   10 CONTINUE
C      RETURN
*
*        code for increment equal to 1
*
C   20 DMAX = DABS(DX(1))
C      DO 30 I = 2,N
C          IF (DABS(DX(I)).LE.DMAX) GO TO 30
C          IDAMAX = I
C          DMAX = DABS(DX(I))
C   30 CONTINUE
C      RETURN
C      END
      LOGICAL FUNCTION LSAME(CA,CB)
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER CA,CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER INTA,INTB,ZCODE
*     ..
*
*     Test if the characters are equal
*
      LSAME = CA .EQ. CB
      IF (LSAME) RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR('Z')
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR(CA)
      INTB = ICHAR(CB)
*
      IF (ZCODE.EQ.90 .OR. ZCODE.EQ.122) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
          IF (INTA.GE.97 .AND. INTA.LE.122) INTA = INTA - 32
          IF (INTB.GE.97 .AND. INTB.LE.122) INTB = INTB - 32
*
      ELSE IF (ZCODE.EQ.233 .OR. ZCODE.EQ.169) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
          IF (INTA.GE.129 .AND. INTA.LE.137 .OR.
     +        INTA.GE.145 .AND. INTA.LE.153 .OR.
     +        INTA.GE.162 .AND. INTA.LE.169) INTA = INTA + 64
          IF (INTB.GE.129 .AND. INTB.LE.137 .OR.
     +        INTB.GE.145 .AND. INTB.LE.153 .OR.
     +        INTB.GE.162 .AND. INTB.LE.169) INTB = INTB + 64
*
      ELSE IF (ZCODE.EQ.218 .OR. ZCODE.EQ.250) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
          IF (INTA.GE.225 .AND. INTA.LE.250) INTA = INTA - 32
          IF (INTB.GE.225 .AND. INTB.LE.250) INTB = INTB - 32
      END IF
      LSAME = INTA .EQ. INTB
*
*     RETURN
*
*     End of LSAME
*
      END
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*  Purpose
*  =======
*
*  DLAMCH determines double precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by DLAMCH:
*          = 'E' or 'e',   DLAMCH := eps
*          = 'S' or 's ,   DLAMCH := sfmin
*          = 'B' or 'b',   DLAMCH := base
*          = 'P' or 'p',   DLAMCH := eps*base
*          = 'N' or 'n',   DLAMCH := t
*          = 'R' or 'r',   DLAMCH := rnd
*          = 'M' or 'm',   DLAMCH := emin
*          = 'U' or 'u',   DLAMCH := rmin
*          = 'L' or 'l',   DLAMCH := emax
*          = 'O' or 'o',   DLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC2
*     ..
*     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,
     $                   EMAX, RMAX, PREC
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
*
      DLAMCH = RMACH
      RETURN
*
*     End of DLAMCH
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  Purpose
*  =======
*
*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         B = 1
         C = DLAMC3( A, B )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         LT = 0
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
*
*     End of DLAMC1
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
*     ..
*
*  Purpose
*  =======
*
*  DLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) DOUBLE PRECISION
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) DOUBLE PRECISION
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
*
*        Start to find EPS.
*
         B = LBETA
         A = B**( -LT )
         LEPS = A
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS
*
         LEPS = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( A.LT.LEPS )
     $      LEPS = A
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
*
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
***
* Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         IEEE = IEEE .OR. LIEEE1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
*
*        Finally, call DLAMC5 to compute EMAX and RMAX.
*
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', I8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
*
*     End of DLAMC2
*
      END
*
************************************************************************
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
*  Purpose
*  =======
*
*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A, B    (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC4( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
*     ..
*
*  Purpose
*  =======
*
*  DLAMC4 is a service routine for DLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) EMIN
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) DOUBLE PRECISION
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Executable Statements ..
*
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of DLAMC4
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
*     ..
*
*  Purpose
*  =======
*
*  DLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest machine floating-point number.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
         EMAX = EMAX - 1
      END IF
*
      IF( IEEE ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
         EMAX = EMAX - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
*
      RMAX = Y
      RETURN
*
*     End of DLAMC5
*
      END
      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines (DEPRECATED)
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR method
*               for nonsymmetric eigenvalue problems (DEPRECATED)
*          = 9: maximum size of the subproblems at the bottom of the
*               computation tree in the divide-and-conquer algorithm
*               (used by xGELSD and xGESDD)
*          =10: ieee NaN arithmetic can be trusted not to trap
*          =11: infinity arithmetic can be trusted not to trap
*          12 <= ISPEC <= 16:
*               xHSEQR or one of its subroutines,
*               see IPARMQ for detailed explanation
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
* (ILAENV) (output) INTEGER
*          >= 0: the value of the parameter specified by ISPEC
*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IC, IZ, NB, NBMIN, NX
      LOGICAL            CNAME, SNAME
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. External Functions ..
      INTEGER            IEEECK, IPARMQ
      EXTERNAL           IEEECK, IPARMQ
*     ..
*     .. Executable Statements ..
*
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120,
     $        130, 140, 150, 160, 160, 160, 160, 160 )ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
   10 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC+64 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:
     $             I ) = CHAR( IC+64 )
   30       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   40       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1: 1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
*
      GO TO ( 50, 60, 70 )ISPEC
*
   50 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
   60 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
   70 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
   80 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
   90 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  100 CONTINUE
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  110 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  120 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
  130 CONTINUE
*
*     ISPEC = 9:  maximum size of the subproblems at the bottom of the
*                 computation tree in the divide-and-conquer algorithm
*                 (used by xGELSD and xGESDD)
*
      ILAENV = 25
      RETURN
*
  140 CONTINUE
*
*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
*
  150 CONTINUE
*
*     ISPEC = 11: infinity arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
*
  160 CONTINUE
*
*     12 <= ISPEC <= 16: xHSEQR or one of its subroutines. 
*
      ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      RETURN
*
*     End of ILAENV
*
      END
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
*     ..
*
*  Purpose
*  =======
*
*  IEEECK is called from the ILAENV to verify that Infinity and
*  possibly NaN arithmetic is safe (i.e. will not trap).
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies whether to test just for inifinity arithmetic
*          or whether to test for infinity and NaN arithmetic.
*          = 0: Verify infinity arithmetic only.
*          = 1: Verify infinity and NaN arithmetic.
*
*  ZERO    (input) REAL
*          Must contain the value 0.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  ONE     (input) REAL
*          Must contain the value 1.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  RETURN VALUE:  INTEGER
*          = 0:  Arithmetic failed to produce the correct answers
*          = 1:  Arithmetic produced the correct answers
*
*     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
*     ..
*     .. Executable Statements ..
      IEEECK = 1
*
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
*
*
*
*     Return if we were only asked to check infinity arithmetic
*
      IF( ISPEC.EQ.0 )
     $   RETURN
*
      NAN1 = POSINF + NEGINF
*
      NAN2 = POSINF / NEGINF
*
      NAN3 = POSINF / POSINF
*
      NAN4 = POSINF*ZERO
*
      NAN5 = NEGINF*NEGZRO
*
      NAN6 = NAN5*0.0
*
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      RETURN
      END
      INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*     
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, ISPEC, LWORK, N
      CHARACTER          NAME*( * ), OPTS*( * )
*
*  Purpose
*  =======
*
*       This program sets problem and machine dependent parameters
*       useful for xHSEQR and its subroutines. It is called whenever 
*       ILAENV is called with 12 <= ISPEC <= 16
*
*  Arguments
*  =========
*
*       ISPEC  (input) integer scalar
*              ISPEC specifies which tunable parameter IPARMQ should
*              return.
*
*              ISPEC=12: (INMIN)  Matrices of order nmin or less
*                        are sent directly to xLAHQR, the implicit
*                        double shift QR algorithm.  NMIN must be
*                        at least 11.
*
*              ISPEC=13: (INWIN)  Size of the deflation window.
*                        This is best set greater than or equal to
*                        the number of simultaneous shifts NS.
*                        Larger matrices benefit from larger deflation
*                        windows.
*
*              ISPEC=14: (INIBL) Determines when to stop nibbling and
*                        invest in an (expensive) multi-shift QR sweep.
*                        If the aggressive early deflation subroutine
*                        finds LD converged eigenvalues from an order
*                        NW deflation window and LD.GT.(NW*NIBBLE)/100,
*                        then the next QR sweep is skipped and early
*                        deflation is applied immediately to the
*                        remaining active diagonal block.  Setting
*                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
*                        multi-shift QR sweep whenever early deflation
*                        finds a converged eigenvalue.  Setting
*                        IPARMQ(ISPEC=14) greater than or equal to 100
*                        prevents TTQRE from skipping a multi-shift
*                        QR sweep.
*
*              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
*                        a multi-shift QR iteration.
*
*              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
*                        following meanings.
*                        0:  During the multi-shift QR sweep,
*                            xLAQR5 does not accumulate reflections and
*                            does not use matrix-matrix multiply to
*                            update the far-from-diagonal matrix
*                            entries.
*                        1:  During the multi-shift QR sweep,
*                            xLAQR5 and/or xLAQRaccumulates reflections and uses
*                            matrix-matrix multiply to update the
*                            far-from-diagonal matrix entries.
*                        2:  During the multi-shift QR sweep.
*                            xLAQR5 accumulates reflections and takes
*                            advantage of 2-by-2 block structure during
*                            matrix-matrix multiplies.
*                        (If xTRMM is slower than xGEMM, then
*                        IPARMQ(ISPEC=16)=1 may be more efficient than
*                        IPARMQ(ISPEC=16)=2 despite the greater level of
*                        arithmetic work implied by the latter choice.)
*
*       NAME    (input) character string
*               Name of the calling subroutine
*
*       OPTS    (input) character string
*               This is a concatenation of the string arguments to
*               TTQRE.
*
*       N       (input) integer scalar
*               N is the order of the Hessenberg matrix H.
*
*       ILO     (input) INTEGER
*       IHI     (input) INTEGER
*               It is assumed that H is already upper triangular
*               in rows and columns 1:ILO-1 and IHI+1:N.
*
*       LWORK   (input) integer scalar
*               The amount of workspace available.
*
*  Further Details
*  ===============
*
*       Little is known about how best to choose these parameters.
*       It is possible to use different values of the parameters
*       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
*
*       It is probably best to choose different parameters for
*       different matrices and different parameters at different
*       times during the iteration, but this has not been
*       implemented --- yet.
*
*
*       The best choices of most of the parameters depend
*       in an ill-understood way on the relative execution
*       rate of xLAQR3 and xLAQR5 and on the nature of each
*       particular eigenvalue problem.  Experiment may be the
*       only practical way to determine which choices are most
*       effective.
*
*       Following is a list of default values supplied by IPARMQ.
*       These defaults may be adjusted in order to attain better
*       performance in any particular computational environment.
*
*       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
*                        Default: 75. (Must be at least 11.)
*
*       IPARMQ(ISPEC=13) Recommended deflation window size.
*                        This depends on ILO, IHI and NS, the
*                        number of simultaneous shifts returned
*                        by IPARMQ(ISPEC=15).  The default for
*                        (IHI-ILO+1).LE.500 is NS.  The default
*                        for (IHI-ILO+1).GT.500 is 3*NS/2.
*
*       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
*
*       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
*                        a multi-shift QR iteration.
*
*                        If IHI-ILO+1 is ...
*
*                        greater than      ...but less    ... the
*                        or equal to ...      than        default is
*
*                                0               30       NS =   2+
*                               30               60       NS =   4+
*                               60              150       NS =  10
*                              150              590       NS =  **
*                              590             3000       NS =  64
*                             3000             6000       NS = 128
*                             6000             infinity   NS = 256
*
*                    (+)  By default matrices of this order are
*                         passed to the implicit double shift routine
*                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
*                         values of NS are used only in case of a rare
*                         xLAHQR failure.
*
*                    (**) The asterisks (**) indicate an ad-hoc
*                         function increasing from 10 to 64.
*
*       IPARMQ(ISPEC=16) Select structured matrix multiply.
*                        (See ISPEC=16 above for details.)
*                        Default: 3.
*
*     ================================================================
*     .. Parameters ..
      INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14,
     $                   ISHFTS = 15, IACC22 = 16 )
      INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14,
     $                   NIBBLE = 14, KNWSWP = 500 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
*     ..
*     .. Local Scalars ..
      INTEGER            NH, NS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LOG, MAX, MOD, NINT, REAL
*     ..
*     .. Executable Statements ..
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR.
     $    ( ISPEC.EQ.IACC22 ) ) THEN
*
*        ==== Set the number simultaneous shifts ====
*
         NH = IHI - ILO + 1
         NS = 2
         IF( NH.GE.30 )
     $      NS = 4
         IF( NH.GE.60 )
     $      NS = 10
         IF( NH.GE.150 )
     $      NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
         IF( NH.GE.590 )
     $      NS = 64
         IF( NH.GE.3000 )
     $      NS = 128
         IF( NH.GE.6000 )
     $      NS = 256
         NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
*
      IF( ISPEC.EQ.INMIN ) THEN
*
*
*        ===== Matrices of order smaller than NMIN get sent
*        .     to xLAHQR, the classic double shift algorithm.
*        .     This must be at least 11. ====
*
         IPARMQ = NMIN
*
      ELSE IF( ISPEC.EQ.INIBL ) THEN
*
*        ==== INIBL: skip a multi-shift qr iteration and
*        .    whenever aggressive early deflation finds
*        .    at least (NIBBLE*(window size)/100) deflations. ====
*
         IPARMQ = NIBBLE
*
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
*
*        ==== NSHFTS: The number of simultaneous shifts =====
*
         IPARMQ = NS
*
      ELSE IF( ISPEC.EQ.INWIN ) THEN
*
*        ==== NW: deflation window size.  ====
*
         IF( NH.LE.KNWSWP ) THEN
            IPARMQ = NS
         ELSE
            IPARMQ = 3*NS / 2
         END IF
*
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
*
*        ==== IACC22: Whether to accumulate reflections
*        .     before updating the far-from-diagonal elements
*        .     and whether to use 2-by-2 block structure while
*        .     doing it.  A small amount of work could be saved
*        .     by making this choice dependent also upon the
*        .     NH=IHI-ILO+1.
*
         IPARMQ = 0
         IF( NS.GE.KACMIN )
     $      IPARMQ = 1
         IF( NS.GE.K22MIN )
     $      IPARMQ = 2
*
      ELSE
*        ===== invalid value of ispec =====
         IPARMQ = -1
*
      END IF
*
*     ==== End of IPARMQ ====
*
      END
      SUBROUTINE ZIBCONST(EPMACH,SMALL)
      DOUBLE PRECISION EPMACH,SMALL
C
C*********************************************************************
C  Set Approximations to machine constants.
C  This routine is machine dependent.
C*********************************************************************
C
C  Output parameters
C    EPMACH    DOUBLE     Relative machine precision
C    SMALL     DOUBLE     SQRT of smallest positive machine number
C 
C  The proposed values should work on Intel compatible cpus,
C  PowerPCs and Sun sparcs
C
C*********************************************************************
C
      EPMACH = 1.0D-17
      SMALL  = 1.0D-150
      RETURN
      END
      SUBROUTINE ZIBSEC (CPTIM, IFAIL)
C
C*********************************************************************
C  Set CPTIM to cpu time in seconds.
C  This routine is machine dependent.
C*********************************************************************
C
C  Output parameters
C    CPTIM    REAL        Cpu time in seconds
C    IFAIL    INTEGER     Errorcode
C 
C  CSUN code works on sun and MacOSX (g77).
C  CXLF code works on MacOSX (xlf)
C  CVIS code works with Visual Fortran (Windows)
C
C*********************************************************************
C
      REAL CPTIM
      INTEGER IFAIL
C
CSUN      REAL RTIME(2)
CXLF      REAL(4) etime_
CXLF      TYPE TB_TYPE
CXLF      SEQUENCE
CXLF      REAL(4) USRTIME
CXLF      REAL(4) SYSTIME
CXLF      END TYPE
CXLF      TYPE (TB_TYPE) ETIME_STRUCT
C
      CPTIM = 0
      IFAIL = 0
C
CSUN      CPTIM = ETIME(RTIME)
CSUN      CPTIM = RTIME(1)
C
CXLF      CPTIM = etime_(ETIME_STRUCT)
C
CVIS      CALL CPU_TIME(CPTIM)
C
      RETURN
      END
      SUBROUTINE MON
C
C* Begin Prologue MON
C
C  ---------------------------------------------------------------------
C
C* Title
C
C  Time monitor to discover cpu time consumption of interesting parts
C  of a code.
C
C* Written by        U. Nowak, U. Poehle, L. Weimann
C* Purpose           Cpu time measuring
C* Category          ...
C* File              mon.f
C* Version           2.1
C* Latest Change     96/03/04 (1.3)
C* Library           CodeLib
C* Code              Fortran 77
C                    Single precision
C
C  ---------------------------------------------------------------------
C
C* Summary
C  -------
C 
C  MON provides different entries 
C   - to initialize the time monitor, 
C   - to start and to stop the measuring, 
C   - to start and to stop a number of stop-watches which may be nested,
C   - to print a summary table,
C   - to define text labels for the summary table, and
C   - to store the vector of average cpu times
C
C  Errors detected during the measuring will disable the time monitor
C  but will not affect the calling program.
C
C  ---------------------------------------------------------------------
C
C* Parameter list description
C  --------------------------
C
C* Input parameters
C  ----------------
C
C  TEXTH      Char    Text to identify the measuring (up to 31
C                     characters will be printed in the headline of the
C                     summary table)
C                     Used by MONINI
C
C  IOUNIT     Int     Output unit for error messages and summary table
C                     Used by MONINI
C
C  INDX       Int     Number of a range of source code lines (a distinct
C                     part of the algorithm normally) for which the
C                     cpu time will be measured by a kind of stop-watch
C                     Used by MONDEF, MONON, and MONOFF
C
C  NAMEH      Char    Text to label the summary lines (up to 17
C                     characters will be printed)
C                     Used by MONDEF
C
C* Output parameters
C  -----------------
C
C  IRET       Int     Return code of ZIBSEC
C                     Used by MONSTR
C
C  AVER       Dble    Vector of average cpu time for all measurings
C                     Used by MONGET
C
C
C* End Prologue
C  ------------
C
C
C* Constants
C  ---------
C
C  MAXTAB     Int     Maximum number of stop-watches
C  MNEST      Int     Maximum number of nested measurings
C
C
C* Local Variables
C  ---------------
C
C  ASEC       Real    Array of averages per call
C  INDACT     Int     Array for indices of nested measurings
C  NCALL      Int     Array of call counts
C  PC1        Real    Array of per cent values with respect to all
C                     measurings
C  PC2        Real    Array of per cent values with respect to the sum
C                     of all stop-watches
C  QDISAB     Log     Time monitor disabled
C  QON        Log     Array reflecting if an index is active
C  QSTART     Log     Time monitor started
C  SEC        Real    Array of cpu time measurings
C
C
C
      INTEGER MAXTAB, MNEST
C
      PARAMETER (MAXTAB = 25, MNEST=20)
C
      INTEGER INDACT(MNEST), NCALL(0:MAXTAB)
      REAL ASEC(0:MAXTAB), PC1(0:MAXTAB), PC2(0:MAXTAB), SEC(0:MAXTAB)
      DOUBLE PRECISION AVER(0:MAXTAB)
      LOGICAL QDISAB, QON(0:MAXTAB), QSTART
      CHARACTER NAME(0:MAXTAB)*17, TEXT*30 
      CHARACTER*(*) NAMEH, TEXTH
C
      SAVE ASEC, CPTIM, INDACT, IONCNT, MAXIND, MONI, NAME, NCALL, PC1,
     $     PC2, QDISAB, QON, QSTART, SEC, TEXT
C
C
      DATA CPTIM /0.0E0/, MONI /6/, INFO /1/, IONCNT /-1/,
     $     QDISAB /.FALSE./, QSTART /.FALSE./
C
      RETURN
C
C
      ENTRY MONINI (TEXTH, IOUNIT)
C     Initialize monitor.
C     Has to be called first. May be called again after MONHLT.
C
      TEXT = TEXTH
      MONI = IOUNIT
C
      IF (IONCNT .GT. 0 .AND. .NOT. QDISAB) GOTO 1070
C
      MAXIND = 0
      IONCNT = 0
      QDISAB = .FALSE.
C
      DO 1000 I = 0,MAXTAB
         SEC(I) = 0.
         ASEC(I) = 0.
         NCALL(I) = 0
         QON(I) = .FALSE.
         WRITE (NAME(I), '(A, I2)') 'Part ', I
 1000 CONTINUE
      NAME(0) = 'General'
C
      DO 1010 I = 1,MNEST
         INDACT(I) = 0
 1010 CONTINUE
C
      RETURN
C
C
      ENTRY MONDEF(INDX, NAMEH)
C     Define one monitor entry.
C     May be called at any time before MONPRT.
C
      IF (QDISAB) RETURN
      IF (INDX.LT.0 .OR. INDX.GT.MAXTAB) GOTO 1080
C
      NAME(INDX) = NAMEH
C
      RETURN
C
C
      ENTRY MONSTR (IRET)
C     Start monitor measurings.
C     Has to be called once after initialization.
C
      IF (QDISAB) RETURN
      IF (IONCNT .LT. 0) GOTO 1090
      IF (IONCNT .GT. 0) GOTO 1100
      IF (QON(0)) GOTO 1110
C
      IFAIL = 0
      CALL ZIBSEC (CPTIM, IFAIL)
C
      IF (IFAIL .EQ. 0) THEN
C        Switch on general stop-watch
         SEC(0) = -CPTIM
         QON(0) = .TRUE.
         IONCNT = 1
         QSTART = .TRUE.
      ENDIF
      IRET = IFAIL
C
      RETURN
C
C
      ENTRY MONON (INDX)
C     Start one measuring.
C     A running stop-watch will be deactivated until the new stop-watch
C     stops.
C
      IF (.NOT. QSTART) RETURN
      IF (QDISAB) RETURN
      IF (IONCNT .LT. 1) GOTO 1120
      IF (INDX .GT. MAXTAB .OR. INDX .LE. 0) GOTO 1130
      IF (QON(INDX)) GOTO 1140
C
      MAXIND = MAX(MAXIND, INDX)
      CALL ZIBSEC (CPTIM, IFAIL)
C
C     Hold actual stop-watch
      SEC(INDACT(IONCNT)) = SEC(INDACT(IONCNT)) + CPTIM
C
C     Switch on new stop-watch
      IF (INFO .GT. 1) WRITE (MONI,*) ' Enter ', NAME(INDX), SEC(INDX)
C
      IONCNT = IONCNT + 1
      IF (IONCNT .GT. MNEST) GOTO 1150
C
      INDACT(IONCNT) = INDX
      SEC(INDX) = SEC(INDX) - CPTIM
      QON(INDX) = .TRUE.
C
      RETURN
C
C
      ENTRY MONOFF (INDX)
C     Stop one measuring.
C     May be called for the stop-watch started most recently.
C     The stop-watch deactivated most recently will be activated.
C
      IF (.NOT. QSTART) RETURN
      IF (QDISAB) RETURN
      IF (INDX .GT. MAXTAB .OR. INDX .LE. 0) GOTO 1160
      IF (INDACT(IONCNT) .NE. INDX) GOTO 1170
C
      CALL ZIBSEC (CPTIM, IFAIL)
C
C     Switch off actual stop-watch
      QON(INDX) = .FALSE.
      SEC(INDX) = SEC(INDX) + CPTIM
      NCALL(INDX) = NCALL(INDX) + 1
      IONCNT = IONCNT - 1
      IF (INFO .GT. 1) WRITE (MONI,*) ' Exit ', NAME(INDX), SEC(INDX)
C
C     Continue previous stop-watch
      SEC(INDACT(IONCNT)) = SEC(INDACT(IONCNT)) - CPTIM
C
      RETURN
C
C
      ENTRY MONHLT
C     Terminate monitor.
C     Stops all active stop-watches.
C
      IF (.NOT. QSTART) RETURN
      IF (QDISAB) RETURN
C
      CALL ZIBSEC (CPTIM, IFAIL)
C
      DO 1020 I=IONCNT,1,-1
         QON(INDACT(I)) = .FALSE.
         SEC(INDACT(I)) = SEC(INDACT(I)) + CPTIM
         NCALL(INDACT(I)) = NCALL(INDACT(I)) + 1
 1020 CONTINUE
C
      IONCNT = 0
C     This means that the time monitor has to be started by MONSTR
C     before another measuring.
C
      RETURN
C
C
C
      ENTRY MONPRT
C     Print statistics.
C     May be called after MONHLT only.
C
      IF (IONCNT .GT. 0) GOTO 1180
      IF (.NOT. QSTART) GOTO 1190
C
      SUM = 1.E - 10
      DO 1030 I = 1,MAXIND
         SUM = SUM + SEC(I)
         IF (NCALL(I) .GT. 0) ASEC(I) = SEC(I)/FLOAT(NCALL(I))
 1030 CONTINUE
      SUM0 = SUM + SEC(0)
      IF (NCALL(0) .GT. 0) ASEC(0) = SEC(0)/FLOAT(NCALL(0))
C
      DO 1040 I = 1,MAXIND
         PC1(I) = 100.*SEC(I)/SUM0
         PC2(I) = 100.*SEC(I)/SUM
 1040 CONTINUE
      PC1(0) = 100.*SEC(0)/SUM0
C
      WRITE (MONI,9000)
      WRITE (MONI,9010)
      WRITE (MONI,9020) ' '
 9000 FORMAT (///)
 9010 FORMAT (1X, 77('#'))
 9020 FORMAT (' #   ', A, T75,  '   #')
C
      IF (QDISAB) THEN
         WRITE (MONI,9020) ' '
         WRITE (MONI,9020)
     $        'Warning  The following results may be misleading'
         WRITE (MONI,9020)
     $        'because an error occured and disabled the time monitor'
      ENDIF
C
      WRITE (MONI,9020) ' '
      WRITE (MONI,9030) 'Results from time monitor program for:', TEXT
 9030 FORMAT (' #   ', A40, 1X, A31, '#')
C
      WRITE (MONI,9020) ' '
      WRITE (MONI,9040) 'Total time:', SUM0, 'Sum of parts:', SUM
 9040 FORMAT (' #   ', A11, F11.3, 5X, A13, F11.3, 21X, '#')
C
      WRITE (MONI,9020) ' '
      WRITE (MONI,9050)
 9050 FORMAT (' #   ', 2X, 'Name', 14X, 'Calls', 7X, 'Time', 4X,
     $     'Av-time', 4X, '% Total', 6X, '% Sum   #')
C
      WRITE (MONI,9060) NAME(0), NCALL(0), SEC(0), ASEC(0), PC1(0)
 9060 FORMAT (' #   ', A17, I8, F11.3, F11.4, F11.2, 14X, '#')
C
      DO 1050 I = 1,MAXIND
         WRITE (MONI,9070) NAME(I), NCALL(I), SEC(I), ASEC(I), PC1(I),
     $        PC2(I)
 1050 CONTINUE
C
 9070 FORMAT (' #   ', A17, I8, F11.3, F11.4, F11.2, F11.2, 3X, '#')
C
C
      WRITE (MONI,9020) ' '
      WRITE (MONI,9010)
      WRITE (MONI,9000)
C
      RETURN
C
C
      ENTRY MONGET(AVER)
C     Store average cpu times per call.
C     May be called at any time after MONPRT.
C
      IF (.NOT. QSTART) RETURN
      IF (QDISAB) RETURN
C
      DO 1060 I=0,MAXIND
         AVER(I) = DBLE(ASEC(I))
 1060 CONTINUE
C
      RETURN
C
C  Error exits
C
 1070 CONTINUE
      WRITE (MONI,9080) 'MONINI', 'Time monitor is running already.'
      GOTO 1200
C
 1080 CONTINUE
      WRITE (MONI,9090) 'MONDEF', 'Index out of range', INDX
      GOTO 1200
C
 1090 CONTINUE
      WRITE (MONI,9080) 'MONSTR',
     $     'Time monitor has to be initialized by MONINI first.'
      GOTO 1200
C
 1100 CONTINUE
      WRITE (MONI,9080) 'MONSTR', 'Time monitor is running already.'
      GOTO 1200
C
 1110 CONTINUE
      WRITE (MONI,9080) 'MONSTR',
     $     'Time monitor has been started already.'
      GOTO 1200
C
 1120 CONTINUE
      WRITE (MONI,9080) 'MONON', 'Time monitor is not yet started.'
      GOTO 1200
C
 1130 CONTINUE
      WRITE (MONI,9090) 'MONON', 'Index out of range', INDX
      GOTO 1200
C
 1140 CONTINUE
      WRITE (MONI,9090) 'MONON',
     $     'Measuring is running already for this INDX', INDX
      GOTO 1200
C
 1150 CONTINUE
      WRITE (MONI,9100) 'MONON', 'Nesting is too deep.',
     $     'The following indices are active', (INDACT(I),I=0,IONCNT)
      GOTO 1200
C
 1160 CONTINUE
      WRITE (MONI,9090) 'MONOFF', 'Index out of range', INDX
      GOTO 1200
C
 1170 CONTINUE
      WRITE (MONI,9110) 'MONOFF', 'Measuring ', INDX,
     $     'cannot be stopped.',
     $     'The following indices are active', (INDACT(I),I=0,IONCNT)
      GOTO 1200
C
 1180 CONTINUE
      WRITE (MONI,9080) 'MONPRT', 'Time monitor is still running.'
      GOTO 1200
C
 1190 CONTINUE
      WRITE (MONI,9080) 'MONPRT', 'Time monitor was not started.'
      GOTO 1200
C
 1200 QDISAB = .TRUE.
      RETURN
C
 9080 FORMAT (/, ' ++ Error in subroutine ', A, ' ++',/, 4X, A)
C
 9090 FORMAT (/, ' ++ Error in subroutine ', A, ' ++',/, 4X, A,
     $     ' (', I6, ').')
C
 9100 FORMAT (/, ' ++ Error in subroutine ', A, ' ++', 4X, A,/, 4X, A,
     $     (10I4))
C
 9110 FORMAT (/, ' ++ Error in subroutine ', A, ' ++', 4X, A, I3, 1X, A,
     $     /, 4X, A, (10I4))
C
C  End subroutine monitor
C
      END
