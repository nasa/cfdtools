C+----------------------------------------------------------------------
      SUBROUTINE XSORT(NPTS,X,IXPTR)
C
C  PURPOSE:
C	Sorts an array of pointers such that the X array is in ascending
C	order when accessed indirectly through the pointers.  I.e.,
C	X(IXPTR(1)) <= X(IXPTR(2)) ... <= X(IXPTR(NPTS)).  Useful for
C	referring to several arrays to be sorted according to one array.
C
C  METHOD:
C	Bubble sort without work space.  References to X are all of form
C	X(IXPTR(I)), which is up to 15% slower than using work space to
C	hold a copy of the X array.  However, this method obviates the
C	need for another argument or for an artificial limit on the
C	number of points (from a local array).
C
C  PARAMETERS:
C  ARGUMENT  TYPE  I/S/O  DIMENSION  DESCRIPTION
C  --------  ----  -----  ---------  -----------
C    NPTS     I      I        -      Number of points in the X and IXPTR arrays.
C    X        R      I      (NPTS)   X array for which IXPTR is to be sorted.
C    IXPTR    I      O      (NPTS)   Pointers, sorted such that Xs so referred
C                                    are sorted in ascending order.
C
C  COMMONS USED:  (None)
C
C  FILES USED:  (None)
C
C  NOTES:  (None)
C
C  NON-STANDARD CODE:  (None)
C
C  EXTERNAL REFERENCES:  (None)
C
C  ENVIRONMENT:  VAX 11/780, VMS V4.4, FORTRAN V4.6
C
C  AUTHOR:  Nancy Karp, Informatics PMI
C
C  REVISIONS:	DATE  PERSON	DESCRIPTION
C		----  ------	-----------
C		 8/87	DLH	Removed limit of number of points
C		 6/83	AVB	Up XTEMP to 300
C		 1/83	DLH	UP XTEMP TO 200
C		10/78	NLK	Original release
C
C-----------------------------------------------------------------------
C
      DIMENSION X(NPTS), IXPTR(NPTS)
C
C ... Initialize pointers as if X were monotically ascending
      DO 50 I=1,NPTS
         IXPTR(I)=I
   50 CONTINUE
C
C ... Check for enough points to sort
      IF(NPTS.LT.2)  GO TO 9999
C
C ... Sort into ascending order
      DO 200 I=1,NPTS-1
         SMALL=X(IXPTR(I))
         DO 100 J=I+1,NPTS
            IF(X(IXPTR(J)).GE.SMALL)  GO TO 100
C ...          SMALL is smallest value; ISMALL is pointer to smallest value
               SMALL=X(IXPTR(J))
               ISMALL=J
  100    CONTINUE
         IF(SMALL.EQ.X(IXPTR(I)))  GO TO 200
C ...       Else need to interchange points in IXPTR
            ITEMP=IXPTR(I)
            IXPTR(I)=IXPTR(ISMALL)
            IXPTR(ISMALL)=ITEMP
  200 CONTINUE
C
 9999 CONTINUE
      RETURN
      END
