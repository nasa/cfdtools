C+----------------------------------------------------------------------
C
      SUBROUTINE OBJEPS (N, X, H, FUN, FUNTYPE, KTH, LUNOUT, FEPS, IER)
C
C  PURPOSE:
C
C        OBJEPS estimates the precision of a function of N variables
C     (one estimate per variable), according to the method of Hamming
C     involving the table of differences.  Such an estimate is desirable
C     for effective application of numerical optimization methods.
C
C        The function is supplied as an external module, FUN, for which
C     some calling sequence must be assumed here.  Two likely choices
C     are supported (with others readily added).  They are:
C
C        FUNTYPE = 'Q' (for QNMDIF):       FUN (N, X, F)
C        FUNTYPE = 'N' (for NPSOL/NZSOL):  FUN (MODE, N, X, F, G, NSTATE).
C
C        This implementation is suited to functions which are not too
C     expensive.  It gives an indication of how precise the function
C     appears as each variable X (I) is varied while remaining elements
C     of X (*) are held constant.  The largest element of FEPS (*) (i.e.,
C     the measure of least precision) should be passed to the optimizer.
C
C        Very expensive functions require a more flexible implementation
C     (forthcoming).
C
C  METHOD:
C
C        A table is constructed for each of the variables in turn, and
C     all of the differences up to K=KTH are calculated, with an estimate
C     of function precision derived from each column of differences.  Only
C     final estimates are returned in FEPS (*), since successive columns
C     of the table use the same local work-space.  The columns need to be
C     printed as rows for this reason.  (Printing may be suppressed.)
C
C        The formula for the estimate from column K for variable I is:
C
C           FEPS (I) = MAX (DEL (K) F) / SQRT (CK)   where J refers to rows
C                       J                            in the Kth column, and
C           CK = (2K!) / (K!) ** 2
C
C        Seven elements are used in the K=KTH difference column where the
C     expected alternating-sign pattern should appear (although no search
C     for such alternation is made - only absolute values are used in the
C     estimate).  KTH is expected to be 4 or 5, but may be as high as the
C     KMAX parameter below (presently 7).
C
C  ARGUMENTS:
C   NAME    DIM   TYPE I/O/S DESCRIPTION
C   N        -     I   I     Number of variables; N >= 1
C   X        N     R   I     Point in N-space at which function precision
C                            is to be estimated
C   H        N     R   I     Finite difference intervals (should be small
C                            so that high-order differences of the function
C                            (involving H**K) are small
C   FUN                I     Name of function routine in a form given above.
C                            In the case of FUNTYPE = 'N', we use
C                            MODE = 0 (no gradients required); and
C                            NSTATE = 1 on first call here; 0 elsewhere.
C   FUNTYPE  -    C*1  I     Uppercase character controlling function type
C                            as indicated above
C   KTH      -     I   I     Difference columns K = 1 : KTH are calculated
C   LUNOUT   -     I   I     LUNOUT > 0 prints the difference tables and
C                                       other results
C                            LUNOUT < 0 suppresses printing (but any error
C                                       message goes to |LUNOUT|)
C   FEPS     N     R     O   FEPS (I) contains the last (KTH) estimate of
C                            precision with respect to variable I at X (*).
C                            The LARGEST element of FEPS (*) is intended to
C                            serve as the function precision parameter input
C                            to the optimization package.
C   IER      -     I     O   IER = 0 means no problems;
C                                  1 means KTH > KMAX (not enough local
C                                    storage)
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77 with IMPLICIT NONE extension
C
C  HISTORY:
C
C  05/09/91  D.A.Saunders  Initial implementation of Hamming's algorithm
C                          for N-dimensional functions.  The algorithm
C                          is applied N times, once for each variable.
C                          It is not suited to very expensive functions.
C  05/10/91    "    "      Adapted to allow FUN in either of the QNMDIF
C                          or NPSOL/NZSOL forms.
C                        
C  AUTHOR:  David Saunders, with Robert Kennelly
C           NASA/Ames Research Center, Moffett Field, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER
     >   IER, KTH, LUNOUT, N
      REAL
     >   FEPS (N), H (N), X (N)
      CHARACTER
     >   FUNTYPE * 1
      EXTERNAL
     >   FUN

C ... Local constants:

      INTEGER
     >   KMAX, NJKTH, NJ0MAX
      REAL
     >   TWO, ZERO
      PARAMETER
     >  (KMAX  = 7,     ! Max. difference column allowed for by local storage
     >   NJKTH = 3,     ! No. of difference rows either side of target for K=KTH
     >   NJ0MAX= NJKTH + (KMAX + 1) / 2,     ! Derived half-length of F column
     >   TWO   = 2.E+0,
     >   ZERO  = 0.E+0)

C ... Local variables:

      INTEGER
     >   I, ISHIFT, J, JP1, JP2, J1, J2, K, LUNERR, MODE, NJ0, NSTATE
      REAL
     >   CK, DELKMAX, DXI, F (-NJ0MAX : NJ0MAX), F0, GRAD, XI

C     Execution:

      IER = 0
      LUNERR = ABS (LUNOUT)
      IF (KTH .GT. KMAX) THEN
         IER = 1
         WRITE (LUNERR, 1000) KMAX
         GO TO 999
      END IF

      IF (LUNOUT .GT. 0) THEN
         WRITE (LUNOUT, 1001) X
         WRITE (LUNOUT, 1002) H
      END IF

C     One function evaluation at X provides the middle value for each table.

      IF (FUNTYPE .EQ. 'Q') THEN   ! QNMDIF application

         CALL FUN (N, X, F0)

      ELSE
         MODE = 0        ! No gradients needed by NPSOL/NZSOL
         NSTATE = 1      ! First call to FUN   "    "    "

         CALL FUN (MODE, N, X, F0, GRAD, NSTATE)
      END IF

      NSTATE = 0
      NJ0 = NJKTH + (KTH + 1) / 2   ! Half-length of function column

C     For each variable...

      DO 500, I = 1, N
         XI  = X (I)
         DXI = H (I)

C        For each row of the function column in the difference table...

         DO 200, J = -NJ0, NJ0
            IF (J .EQ. 0) THEN
               F (J) = F0
            ELSE
               X (I) = XI + REAL (J) * DXI

               IF (FUNTYPE .EQ. 'Q') THEN
                  CALL FUN (N, X, F (J))
               ELSE
                  CALL FUN (MODE, N, X, F (J), GRAD, NSTATE)
               END IF
            END IF
  200    CONTINUE

         X (I) = XI

         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1003)
     >      I, (XI + REAL (J) * DXI, J = -NJKTH, NJKTH),
     >      (F (J), J = -NJKTH, NJKTH)

C        For each column K of the difference table...

         CK = TWO
         J1 = -NJ0
         J2 =  NJ0 - 1   ! Difference column shortens by 1 with each K
         JP1 = -NJKTH    ! Window on printed portion of table is centered
         JP2 =  NJKTH    ! around the input X (I), but it zig-zags...

         DO 400, K = 1, KTH
            DELKMAX = ZERO

C           For each row of the table...

            DO 300, J = J1, J2
               F (J) = F (J + 1) - F (J)             ! Has the effect of making
               DELKMAX = MAX (ABS (F (J)), DELKMAX)  ! the columns drift "up"
  300       CONTINUE           

            FEPS (I) = DELKMAX / SQRT (CK)
            CK = TWO * CK * REAL (K + K + 1) / REAL (K + 1)

            IF (LUNOUT .GT. 0) THEN
               ISHIFT = MOD (K, 2)  ! The print window must shift "up" with
               JP1 = JP1 - ISHIFT   ! the difference columns in order to
               JP2 = JP2 - ISHIFT   ! stay centered ... (Draw a picture!)
               WRITE (LUNOUT, 1004) K, (F (J), J = JP1, JP2), FEPS (I)
            END IF

            J2 = J2 - 1
  400    CONTINUE
  500 CONTINUE

      IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1005) FEPS

  999 RETURN

C     Formats:

 1000 FORMAT ('0FEPSILON: Too high a difference column requested.', /,
     >        ' Present limit: ', I2)
 1001 FORMAT (//, ' Current (central) values of each variable:',
     >        //, (1P, 1X, 7E14.6))
 1002 FORMAT (//, ' Finite differencing intervals for each variable:',
     >        //, (1P, 1X, 7E10.2))
 1003 FORMAT (//, ' Difference table for variable I =', I2, ':',
     >        //, ' X (I):', 1P, 7E14.6,
     >        /,  ' F (J):', 7E14.6)
 1004 FORMAT (' D (', I1, '):', 1P, 7E14.6, '   Eps:', E10.2)
 1005 FORMAT (//, ' Function precision estimates for each variable:',
     >        //, (1P, 1X, 7E10.2))
      END
