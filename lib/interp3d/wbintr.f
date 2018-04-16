C+------------------------------------------------------------------------------
C
      SUBROUTINE WBINTR (IDIMW, JDIMW, I1W, I2W, J1W, J2W, XW, YW, ZW,
     >                   IDIMB, JDIMB, I1B, I2B, J1B, J2B, XB, YB, ZB,
     >                   IDEGW, IDEGB, UB, VB,
     >                   XINTR, YINTR, ZINTR, UINTR, VINTR, TINTR,
     >                   IINTR, JINTR, KINTR, FIRST, LUNWR, IFAIL, FAIL)
C
C     WBINTR calculates a "wing"/"body" intersection given a regular mesh for
C     each surface.  The surfaces are quite arbitrary.  Here, "wing" refers to
C     the surface which is treated one "I" line at a time from I1W to I2W in
C     order along the intersection (increasing or decreasing "I"), while the
C     "body" is treated as a parametric bilinear or bicubic surface.
C
C     Provision is made for saving and reusing surface/line indices/(u,v,t)
C     values for repeated application to related cases. The option to retry
C     each "I" multiple times may be suppressed if necessary for efficient
C     handling of cases where the intersection curve hits a surface boundary.
C
C     More than one call per intersection may be appropriate, as in the case
C     of proceeding from leading to trailing edge for both the upper and lower
C     surfaces of a wing/body intersection.  This forces UB and VB to be
C     arguments, not automatic arrays.  A given application of WBINTR should
C     determine the active surface/line indices, manage the starting guess
C     inputs, and handle any off-the-surface parts of the intersection curve
C     once IFAIL is returned.
C
C     Environment: Fortran 90
C
C     History:
C
C     11/21/97  DAS  Adaptation of the more specialized WFINTR wing/body
C                    intersection routine from the SYN87 design code, in
C                    which James Reuther had wrapped a systematic retry
C                    scheme around the underlying INTSEC4 utility.
C                    WBINTR should be applicable to all types of intersection.
C
C     Author: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER IDIMW, JDIMW         ! I   Dimensions of "wing" section arrays
      INTEGER I1W, I2W             ! I   Range of "I" along desired intersection
      INTEGER J1W, J2W             ! I   Range of "J" used for each "I" line
                                   !     N.B.: J1W < J2W, but I1W > I2W is OK.
      REAL    XW(IDIMW,JDIMW),     ! I   Regularized "wing" sections
     >        YW(IDIMW,JDIMW),
     >        ZW(IDIMW,JDIMW)
      INTEGER IDIMB, JDIMB         ! I   Dimensions of "body" surface arrays
      INTEGER I1B, I2B,            ! I   Indices defining the active surface
     >        J1B, J2B
      REAL    XB(IDIMB,JDIMB),     ! I   Regularized "body" sections
     >        YB(IDIMB,JDIMB),
     >        ZB(IDIMB,JDIMB)
      INTEGER IDEGW                ! I   1 = piecewise linear "I" lines;
                                   !     2 = monotonic ("tight") cubic splines
                                   !     3 = "loose" local cubic splines
      INTEGER IDEGB                ! I   1 = bilinear surface interpolations;
                                   !     3 = bicubic     "       "       "
      REAL    UB(IDIMB,JDIMB),     ! IOS (u,v)s for the "body" surface;
     >        VB(IDIMB,JDIMB)      !     UB(I1B,J1B) = -1. on input forces
                                   !     INTSEC4/5 to calculate them; they may
                                   !     be reqd. on another call with new I1/2W
      REAL    XINTR(IDIMW),        !  O  Intersection coordinates
     >        YINTR(IDIMW),
     >        ZINTR(IDIMW)
      REAL    UINTR(IDIMW),        !  O  Corresp. surface "u" & "v" values for
     >        VINTR(IDIMW),        !     the "body" & "t" values for the "wing"
     >        TINTR(IDIMW)
      INTEGER IINTR(IDIMW),        ! IO  Corresp. cell indices for the "body"
     >        JINTR(IDIMW),        !     (I, J) and the "wing" lines (K);
     >        KINTR(IDIMW)         !     see FIRST
      LOGICAL FIRST                !     T means I/J/K/U/V/TINTR(I1W) ONLY are
                                   !       used as inputs (then values from the
                                   !       previous I are used up to I2W);
                                   !     F means I/J/K/U/V/TINTR(I) are used as
                                   !       starting guesses for ALL I lines
      INTEGER LUNWR                ! I   Logical unit for diagnostics if > 0
      INTEGER IFAIL                !  O  Trouble I index if FAIL = T on return
      LOGICAL FAIL                 ! IO  T on input means DON'T use the retry
                                   !       scheme upon a failure for an I after
                                   !       I = I1W; the I is returned as IFAIL
                                   !     T on output means failure after either
                                   !       1 try or many tries according to the
                                   !       input value of FAIL
C-------------------------------------------------------------------------------

C     Local constants:

      INTEGER, PARAMETER ::
     >   TEN = 10
      REAL, PARAMETER ::
     >   ONE = 1., TENTH = 0.1, TOLER = 1.E-7, UNDEFINED = -999.,
     >   ZERO = 0.
      LOGICAL, PARAMETER ::
     >   CLOSED = .FALSE.

C     Local variables:

      INTEGER
     >   I, IBODY, IER, INC, ISMAT, ITOTER, ITNT, IUNT, IVNT,
     >   J, JBODY, JSMAT, KCURV, LUNERR
      INTEGER, SAVE ::
     >   NPRINT
      REAL
     >   EPS, TINT, TSAVE, TTOTAL, UINT, USAVE, VINT, VSAVE,
     >   XINT, YINT, ZINT, SMATRIX(4,4,3),
     >   TT(J2W), XT(J2W), YT(J2W), ZT(J2W)
      LOGICAL
     >   RETRY
      CHARACTER * 1
     >   METHOD
      CHARACTER * 1, SAVE ::
     >   METHODS(3)

C     N.B.: TT(J1W:J2W), etc., would hamper wing surface gridding - don't do it.

C     Storage:

      DATA
     >   METHODS /'L', 'M', 'B'/, NPRINT /0/

C     Execution:

      METHOD = METHODS(IDEGW)
      LUNERR = SIGN (LUNWR, -1)  ! Force -ve to suppress INTSEC* iterations
      EPS    = MAX (5.* EPSILON (TOLER), TOLER)  ! For p, q in [0, 1], etc.
      RETRY  = .NOT. FAIL
      FAIL   = .FALSE.

      IF (FIRST) THEN  ! Guesses for first I must be input
         I     = I1W
         UINT  = UINTR(I)
         VINT  = VINTR(I)
         TINT  = TINTR(I)
         IBODY = IINTR(I)
         JBODY = JINTR(I)
         KCURV = KINTR(I)

         IF (NPRINT == 0) THEN  ! First may be true more than once
            NPRINT = 1
            IF (LUNWR .GT. 0) WRITE (LUNWR, '(/,A,1P,E10.3)')
     >         ' WBINTR: Tolerance passed to INTSEC4/5 =', EPS
         END IF

      END IF

      XINT  = UNDEFINED  ! In case of printing for I = I1W
      YINT  = UNDEFINED
      ZINT  = UNDEFINED
      ISMAT = I1B - 1    ! Surface cell matrix is initially unusable
      JSMAT = J1B - 1

      IF (I2W .GE. I1W) THEN
         INC = 1
      ELSE
         INC = -1
      END IF

      DO 100 I = I1W, I2W, INC

         DO J = J1W, J2W
            XT(J) = XW(I,J)
            YT(J) = YW(I,J)
            ZT(J) = ZW(I,J)
         END DO

         ITOTER = 0
         IUNT   = 0
         IVNT   = 0
         ITNT   = 0

         IF (.NOT. FIRST) THEN
            UINT  = UINTR(I)
            VINT  = VINTR(I)
            TINT  = TINTR(I)
            IBODY = IINTR(I)
            JBODY = JINTR(I)
            KCURV = KINTR(I)
         END IF

         USAVE  = UINT
         VSAVE  = VINT
         TSAVE  = TINT
         TTOTAL = ZERO  ! INTSEC* parameterizes the line

   50    CONTINUE

         IF (IDEGB .EQ. 3) THEN

            CALL INTSEC4
     >        (IDIMB, JDIMB, I1B, I2B, J1B, J2B, XB, YB, ZB, UB, VB,
     >         IBODY, JBODY, ISMAT, JSMAT, SMATRIX, EPS,
     >         J1W, J2W, XT, YT, ZT, TT, KCURV, TTOTAL, METHOD, CLOSED,
     >         TINT, UINT, VINT, XINT, YINT, ZINT, LUNERR, IER)
         ELSE

            CALL INTSEC5
     >        (IDIMB, JDIMB, I1B, I2B, J1B, J2B, XB, YB, ZB, UB, VB,
     >         IBODY, JBODY, EPS,
     >         J1W, J2W, XT, YT, ZT, TT, KCURV, TTOTAL, METHOD, CLOSED,
     >         TINT, UINT, VINT, XINT, YINT, ZINT, LUNERR, IER)
         END IF

         IF (IER .GT. 2) THEN

            IF (I .NE. I1W .AND. .NOT. RETRY) THEN

C              Normally this should mean the intersection dropped off the edge,
C              which needs to be allowed for at the higher level.

               IFAIL = I
               FAIL  = .TRUE.
               GO TO 999
            END IF

            IF (ITOTER .LT. 3) THEN

               IF (LUNWR .GT. 0) THEN
                  WRITE (LUNWR, '(/,1X,A,I4,3X,A,I4,3X,A)')
     >               'WBINTR: Intersection failure. IER:', IER,
     >               'Curve I:', I, 'Curve X,Y,Z:'
                  WRITE (LUNWR, '(3F14.6)')
     >               (XT(J), YT(J), ZT(J), J = J1W, J2W)
                  WRITE (LUNWR, '(A,3I5)')
     >               ' Surface I, J & curve K: ', IBODY, JBODY, KCURV
                  WRITE (LUNWR, '(A,3F10.6,/,A,3E13.4,/,A)')
     >               ' U, V, T: ', UINT, VINT, TINT,
     >               ' X, Y, Z: ', XINT, YINT, ZINT,
     >               ' Change initial guess & try again.'
               END IF

               IF (ITOTER .EQ. 0) THEN
                  IUNT = IUNT + 1
                  UINT = REAL (IUNT) * TENTH
                  VINT = VSAVE
                  TINT = TSAVE
                  IF (IUNT .EQ. TEN) ITOTER = ITOTER + 1
                  GO TO 50
               ELSE IF (ITOTER .EQ. 1) THEN
                  IVNT = IVNT + 1
                  VINT = REAL (IVNT) * TENTH
                  UINT = USAVE
                  TINT = TSAVE
                  IF (IVNT .EQ. TEN) ITOTER = ITOTER + 1
                  GO TO 50
               ELSE IF (ITOTER .EQ. 2) THEN
                  ITNT = ITNT + 1
                  TINT = REAL (ITNT) * TENTH
                  UINT = USAVE
                  VINT = VSAVE
                  IF (ITNT .EQ. TEN) ITOTER = ITOTER + 1
                  GO TO 50
               END IF
            ELSE
               IFAIL = I
               FAIL  = .TRUE.
               IF (LUNWR .GT. 0) THEN
                  WRITE (LUNWR,'(A)') ' WBINTR: Aborting.'
               END IF
               GO TO 999

            END IF

         END IF

         XINTR(I) = XINT
         YINTR(I) = YINT
         ZINTR(I) = ZINT
         UINTR(I) = UINT
         VINTR(I) = VINT
         TINTR(I) = TINT
         IINTR(I) = IBODY
         JINTR(I) = JBODY
         KINTR(I) = KCURV

  100 CONTINUE

  999 RETURN
      END
