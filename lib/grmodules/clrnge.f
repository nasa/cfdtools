C+**********************************************************************
C
      SUBROUTINE CLRNGE ( MXCOLR, COLARY, NVALS, VALS )
C
C ACRONYM: Set CoLor RaNGE
C              - -   - ---
C PURPOSE:
C    CLRNGE assigns real valued color indices over a range to the
C    elements of an array. It handles the case where the size of the
C    array is greater than the maximum number of colors available.
C
C ARGUMENTS:
C     ARG   TYPE I/O/S  DIM     DESCRIPTION
C    MXCOLR   I    I     -      Maximum number of colors available.
C    COLARY   R    S  MXCOLR    Scratch array holds base indices.
C    NVALS    I    I     -      Number of array elements to be assigned.
C    VALS     R    O  NVALS     Array of color indices from BLUE=7.
C                               down to MAGENTA=2.
C EXTERNAL REFERENCES: None
C
C STANDARDS VIOLATIONS: None.
C
C ENVIRONMENT:  VAX/VMS -- FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE     INITIALS   DESCRIPTION 
C  26 Mar.1987  RCL       Original design and implementation.
C
C AUTHOR: Rosalie Lefkowitz, Sterling Software, Palo Alto, CA.
C
C-**********************************************************************

      IMPLICIT   NONE

      INTEGER    INDEX, LEV, MXCOLR, NVALS

      REAL       CINC, COLARY(1), SCALE, VALS(1)


      IF ( NVALS.LE.MXCOLR ) THEN
            CINC = 5./REAL(NVALS-1)
            DO 100 LEV=1,NVALS
               VALS(NVALS-LEV+1) = 2. + (LEV-1)*CINC
 100        CONTINUE
         ELSE
            CINC = 5./REAL(MXCOLR-1)
            DO 200 LEV=1,MXCOLR
               COLARY(MXCOLR-LEV+1) = 2. + (LEV-1)*CINC
 200        CONTINUE
            SCALE = REAL(MXCOLR)/REAL(NVALS)
            DO 300 LEV = 1,NVALS
               INDEX = SCALE*LEV
               INDEX = MAX ( INDEX, 1 )
               VALS(LEV) = COLARY(INDEX)
 300        CONTINUE
      END IF

      RETURN
      END
