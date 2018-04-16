C+----------------------------------------------------------------------
C
      SUBROUTINE NMLIZX ( DOWN, NX, XIN, XOUT, XMIN, XMAX )
C
C  PURPOSE: NMLIZX normalizes the given elements to the interval [0, 1],
C           with options to reverse the process or to determine the data
C           range first.   Put another way,  it  scales the given points
C           XIN (*) "up" or "down" according to the following formulas:
C
C           DOWN = TRUE:    XOUT = ( XIN - XMIN ) / ( XMAX - XMIN )
C
C                           where XMIN and XMAX are used as inputs
C                           IF THEY ARE NOT EQUAL, else they are
C                           determined by NMLIZX and returned as
C                           outputs. (See NOTES.)
C
C           DOWN = FALSE:   XOUT = XMIN + ( XMAX - XMIN ) * XIN
C
C                           where XMIN and XMAX are inputs (possibly
C                           determined from a prior call to NMLIZX).
C
C  NOTES:   (1) The option to use supplied XMIN and XMAX values even
C               if DOWN=TRUE permits normalizing different datasets in
C               the same way (apart from permitting some efficiency
C               if the data range is already known).
C           (2) The input elements are not assumed to be ordered -
C               NMLIZX was prompted by the least squares problem of
C               polynomial fitting, where the data may be scattered.
C           (3) Scaling in-place is permissible.
C           (4) Normalizing by a zero range makes no sense, and is
C               signalled to the calling program as output XMIN=XMAX.
C
C  PARAMETERS:
C   ARG    TYPE  I/O/S   DIM     DESCRIPTION
C   DOWN    L      I      -    DOWN=T means "normalize";
C                              DOWN=F means "denormalize"
C   NX      I      I      -    Number of points to be scaled; NX >= 1
C   XIN     R      I     NX    Data to be normalized/de-normalized;
C                              need not be monotonic
C   XOUT    R      O     NX    Scaled results; may be same array as XIN;
C                              may be undefined - see XMIN, XMAX
C   XMIN,   R     I/O     -    See PURPOSE and NOTES. A summary:
C   XMAX                  -    DOWN=T means if input XMIN=XMAX then
C                                     determine data range as output,
C                                     else use input XMIN and XMAX
C                                     to "normalize" the data;
C                              DOWN=F means use input XMIN and XMAX to
C                                     "denormalize" the data;
C                              XMIN=XMAX on return means XOUT(*) is not
C                              defined because of implied divide by zero
C
C  ERROR HANDLING:  See NOTE (4) and XMIN, XMAX description.
C
C  ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C  HISTORY:
C   June '86  Original implementation, to deal with fitting polynomials
C             to data involving abcissas like 1985., 1986., ..., which
C             can cause single-precision failure even in robust methods.
C   Dec. '86  Name changed from SCALEX to NMLIZX after need for similar
C             routine XFORMX arose.
C
C  AUTHOR: David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

C ... Arguments:

      INTEGER    NX
      REAL       XIN (NX), XOUT (NX), XMIN, XMAX
      LOGICAL    DOWN

C ... Local Variables:

      INTEGER    I
      REAL       RANGE

C ... Execution:

      RANGE = XMAX - XMIN
     
      IF ( DOWN ) THEN

C ...    "Normalize".  Can we use the given data range?

         IF ( RANGE .EQ. 0.E+0 ) THEN

C ...       No - calculate data range, assuming no particular ordering:

            XMIN = XIN (1)
            XMAX = XMIN

            DO 20, I = 2, NX
               XMIN = MIN ( XIN (I), XMIN )
               XMAX = MAX ( XIN (I), XMAX )
   20       CONTINUE

            RANGE = XMAX - XMIN
         END IF

         IF ( XMIN .NE. XMAX ) THEN

C ...       The range is valid - it may or may not be derived from
C           the current dataset - hence the option to have it input:

            RANGE = 1.E+0 / RANGE
            DO 30, I = 1, NX
               XOUT (I) = ( XIN (I) - XMIN ) * RANGE
   30       CONTINUE

C        ELSE
C ...       Drop out with XMAX = XMIN as a signal that something is
C           abnormal.
         END IF

      ELSE

C ...    DOWN is FALSE, meaning scale "up" or "denormalize":

         DO 40, I = 1, NX
            XOUT (I) = XIN (I) * RANGE + XMIN
   40    CONTINUE

      END IF
C
      RETURN
      END
