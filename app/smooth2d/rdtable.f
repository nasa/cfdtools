C+----------------------------------------------------------------------

      SUBROUTINE RDTABLE (MAXN, NX, NY, X, Y, F, LUNIN, LUNOUT, LUNCRT)
C
C  Acronym:  ReaD rectangular TABLE of data points (X,Y,F)
C  --------  -  -             -----
C
C  Description:
C  ------------
C
C        RDTABLE reads a formatted bivariate dataset defined on a
C     rectangular mesh.  This version determines which variable changes
C     most rapidly by examining the first two data points.  All reads
C     are list-directed.
C
C        This version does not deal with the data title, and assumes that
C     NX and NY have been read already, in order to be usable on the multi-
C     dataset files now handled by SMOOTH2D.
C
C  INPUT DATASET FORMAT:
C
C       <TITLE                ! title line - no longer read by RDTABLE>
C       <NX   NY              ! no. of columns and rows resp. (>=2)>
C        X    Y    F          ! first data point
C        X    Y    F          ! second data point  (either same X or same Y)
C        :    :    :          :
C        :    :    :          :
C        X    Y    F          ! last data point
C       <NX   NY              ! repeat for further tables>
C        :    :    :          :
C
C  Error Handling:
C  ---------------
C
C        RDTABLE writes diagnostics to LUNOUT and LUNCRT (which may be the
C     same unit) and terminates if X, Y values are not compatible with
C     rectangularity.  This version no longer checks for bad NX, NY.
C     EOF encountered anywhere is abnormal.
C
C  Environment:
C  ------------
C
C     DEC VAX/VMS, FORTRAN
C     (Any mainframe should be OK.  IMPLICIT NONE and ! comments are
C     only extensions used.)
C     
C  History:
C  --------
C
C     83-11-01   Fritsch, et al., Livermore:
C                Initial READER implementation for BIMOND package.
C     83-11-02   Renumbered statement labels in numerical order.
C     83-11-03   1. Added tests for rectangular mesh.
C                2. Added coding to determine which variable changes
C                   most rapidly.
C     83-11-08   1. Added 1000 to all statement labels, to avoid possible
C                   conflict with numbers used in BIMOND.
C                2. Minor cosmetic changes.
C     87-03-20   Minor changes to output formats.
C
C     87-11-27   David Saunders, Sterling Software:
C                Adapted READER as RDTABLE for use by SMOOTH2D:
C                1. Went to list-directed reads; used EOF rather than
C                   blank line to indicate "no-more tables."
C                2. Introduced DO ... LUN = LUN1, LUN2, LUN2 - LUN1
C                   to allow LUN1 = LUN2 and condense error handling, which
C                   is not in-line now (more trouble than it was worth).
C     89-06-12   TITLE and NX, NY, should now be read (and range-checked)
C                externally, as needed for SMOOTH2D's handling of more than
C                one dataset per file and distinguishing between rectangular
C                and random data.
C
C ----------------------------------------------------------------------

      IMPLICIT   NONE

C  Arguments:
C  ----------

      INTEGER    MAXN           ! Max. no. of rows in calling prog.  (I)
      INTEGER    NX             ! No. of columns found by caller     (I)
      INTEGER    NY             !  "   "   rows    "   "   "         (I)
      REAL       X (*)          ! X-coords. for any row  in X (1:NX) (O)
      REAL       Y (*)          ! Y-coords.  "  all rows in Y (1:NY) (O)
      REAL       F (MAXN, *)    ! Function values in  F (1:NX, 1:NY) (O)
      INTEGER    LUNIN          ! Logical unit for input dataset     (I)
      INTEGER    LUNOUT         !  "   "   "   "   print file if any (I)
      INTEGER    LUNCRT         !  "   "   "   "   screen if any     (I)
                                ! Note: LUNCRT=LUNOUT is OK.

C  Local Variables:
C  ----------------

      INTEGER    I, INC, J, LUN, LUN1, LUN2

      REAL       FIN, XIN, YIN

C-----------------------------------------------------------------------


C  Execution:
C  ----------

      LUN1 = LUNOUT
      LUN2 = LUNCRT
      INC = LUN2 - LUN1            ! INC may be negative,
      IF (INC .EQ. 0) INC = 1      ! but it must not be zero.
C*****      NX = 0
C*****      NY = 0
C*****      READ (LUNIN, *, END=999, ERR=810) NX, NY

C*****      IF ((NX .LT. 2) .OR. (NX .GT. MAXN)) GO TO 820
C*****      IF ((NY .LT. 2) .OR. (NY .GT. MAXN)) GO TO 830

C     Read first two points.

      READ (LUNIN, *, ERR=840)  X (1), Y (1), F (1, 1)
      READ (LUNIN, *, ERR=840)  XIN, YIN, FIN

      IF (XIN .EQ. X (1)) THEN   ! Y runs most rapidly.

         Y (2) = YIN
         F (1, 2) = FIN
         DO 610, I = 1, NX
            DO 600, J = 1, NY

C              Skip two points already read.

               IF ((I .EQ. 1) .AND. (J .LE. 2))  GO TO 600

               READ (LUNIN, *, ERR=840)  XIN, YIN, F (I, J)

               IF (J .EQ. 1)  THEN
                  X (I) = XIN
               ELSE IF (XIN .NE. X (I)) THEN
                  GO TO 850
               END IF

               IF ( I.EQ.1 )  THEN
                  Y (J) = YIN
               ELSE IF (YIN .NE. Y (J))  THEN
                  GO TO 850
               ENDIF

  600       CONTINUE
  610    CONTINUE

      ELSE IF (YIN .EQ. Y (1)) THEN   ! X  runs most rapidly.

         X (2) = XIN
         F (2, 1) = FIN
         DO 710, J = 1, NY
            DO 700, I = 1, NX

C              Skip two points already read.

               IF ((J .EQ. 1) .AND. (I .LE. 2)) GO TO 700

               READ (LUNIN, *, ERR=840)  XIN, YIN, F (I, J)

               IF (I .EQ. 1)  THEN
                  Y (J) = YIN
               ELSE IF (YIN .NE. Y(J))  THEN
                  GO TO 850
               ENDIF

               IF (J .EQ. 1)  THEN
                  X (I) = XIN
               ELSE IF (XIN .NE. X (I))  THEN
                  GO TO 850
               ENDIF

  700       CONTINUE
  710    CONTINUE

      ELSE   ! The data are nonrectangular.

         I = 0
         J = 0
         GO TO 850

      ENDIF

      GO TO 999


C  Error Handling:
C  ---------------

C*****  810 DO 812, LUN = LUN1, LUN2, INC
C*****         WRITE (LUN, 1000) 'Error reading NX, NY.'
C*****  812 CONTINUE
C*****      GO TO 900

C*****  820 DO 822, LUN = LUN1, LUN2, INC
C*****         WRITE (LUN, 1000) 'NX is too big or less than 2.'
C*****         WRITE (LUN, 1010) 'NX: ', NX, 'MAXN: ', MAXN
C*****  822 CONTINUE
C*****      GO TO 900

C*****  830 DO 832, LUN = LUN1, LUN2, INC
C*****         WRITE (LUN, 1000) 'NY is too big or less than 2.'
C*****         WRITE (LUN, 1010) 'NY: ', NY, 'MAXN: ', MAXN
C*****  832 CONTINUE
C*****      GO TO 900

  840 DO 842, LUN = LUN1, LUN2, INC
         WRITE (LUN, 1000) 'Error reading X, Y, F.'
  842 CONTINUE
      GO TO 900

  850 DO 852, LUN = LUN1, LUN2, INC
         WRITE (LUNOUT, 1020)  I, J, XIN, YIN
  852 CONTINUE
      GO TO 900

  900 STOP 'Abort.'

  999 RETURN


C  Formats:
C  --------

 1000 FORMAT (//, '   *** RDTABLE: ', A)
C***** 1010 FORMAT (1X, A, I5)
 1020 FORMAT (//, '   *** RDTABLE: Non-rectangular data found where:',
     >        /, '  I =', I3, ',  J =', I3, ',  X =', 1P, E15.8,
     >        ',  Y =', E15.8)

      END                ! End of RDTABLE
