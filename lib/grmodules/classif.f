C+------------------------------------------------------------------------------
C
      SUBROUTINE CLASSIF (ORIENT, CLASS, EXPLAN, LENTXT)
C
C     One liner:  Writes classification headings on plots using DISSPLA.
C
C     Purpose:    CLASSIF is intended to write special headings for plots
C                 in either landscape or portrait mode.  The format consists
C                 of both a classification status in large print and an
C                 explanation, both to appear at the top and bottom of
C                 the page.  The headings are suitable for 8.5 x 11" paper.
C
C     Method:     The titles are written using DISSPLA routines.  They
C                 are centered, and positioned at hard-coded distances
C                 from the top and bottom of the page.  Text heights are
C                 also hard-coded.
C                 
C                 CAUTION:  Care must be taken since a new subplot is
C                 defined by the routine, and the page size, plot area and
C                 origin are changed.  CLASSIF should be used after all
C                 other plotting is done.
C
C     Arguments: 
C        ARG         DIM    TYPE   I/O/S    DESCRIPTION
C      ORIENT       LENDIC   C*1     I     'Landscape' or 'Portrait' mode
C      CLASS          *       C      I      Classification title
C      EXPLAN         *       C      I      Explanation title
C      LENTXT         2       I      I      Lengths of titles. Pass 100 to
C                                           invoke self-counting option.
C
C     Procedures (all DISSPLA utilities):
C         AREA2D    Sets the plot area in inches
C         BSHIFT    Used with (0.,0.) to ensure origin is at edge of page
C         ENDGR     Terminates previous subplot
C         HEIGHT    Sets character height
C         MESSAG    Prints text string at specified position
C         PAGE      Sets the page size
C         PHYSOR    Sets the origin
C         XMESS     Returns the length of a character string
C
C     Environment:  VAX/VMS; FORTRAN 77; DISSPLA V11-9003
C
C     History:  03/02/89  M.D.Wong     Initial design and coding.
C               08/27/89   "   "       Updated for DISSPLA V11.0.
C               02/22/90   "   "       Added LENTXT to argument list.
C               08/20/90  D.A.Saunders Recovered DISSPLA version from
C                                      SMDLIB version accidentally installed.
C               03/26/91   "   "       PostScript/Landscape combination was
C                                      giving trouble.  Page size cannot
C                                      exceed 8.5 x 11 it seems.  Use of
C                                      BSHIFT instead (in QPLOT and here)
C                                      got around the problem.
C
C     Author:  Michael Wong, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   LENTXT (2)
      CHARACTER
     >   ORIENT * 1, CLASS * (*), EXPLAN * (*)

C     Local constants.

      REAL
     >   HALF, HTIT, HEXP, SAFETY, SPACE
      PARAMETER
     >  (HALF = 0.5, HTIT = .15, HEXP =.12, SAFETY = .20, SPACE = .10)
      CHARACTER
     >   BLANK * 1
      PARAMETER
     >  (BLANK = ' ')

C     Local variables.

      REAL
     >   BOTTOM, TOP, XCENTER, XPOSN, YPOSN

C     Procedures.

      EXTERNAL
     >   AREA2D, BSHIFT, ENDGR, HEIGHT, MESSAG, PAGE, PHYSOR, XMESS
      REAL
     >   XMESS

C     Execution.

      IF (CLASS .EQ. BLANK .AND. EXPLAN .EQ. BLANK) GO TO 99
   
C     Terminate current subplot.

      CALL ENDGR (0)

C     Set page values according to landscape or portrait mode.

      IF (ORIENT .EQ. 'L') THEN
         CALL BSHIFT (0., 0.)  ! Reset any earlier fudge used to lower the plot
         CALL PAGE (11.0, 8.5)
         XCENTER = 5.75        ! 1/4" fudge for nicer centering of typical plots
         BOTTOM = 0.125        ! 1/8" fudge needed for landscape/Postscript
         TOP = 8.5 - BOTTOM    ! (SAFETY below isn't enough)
      ELSE
         CALL PAGE (8.5, 11.0)      
         XCENTER = 4.50
         BOTTOM = 0.0
         TOP = 11.0
      END IF

C     Define physical origin and new subplot plot area.  Is 6.5 arbitrary?

      CALL PHYSOR (0., 0.)
      CALL AREA2D (6.5, 6.5)

C     Write classification titles in large letters.

      IF (CLASS .NE. BLANK) THEN

         CALL HEIGHT (HTIT)

C        Center first title at top of page, then at bottom of page.
  
         XPOSN = XCENTER - HALF * XMESS (CLASS, LENTXT (1))
         YPOSN = TOP - SAFETY - HTIT
         CALL MESSAG (CLASS, LENTXT (1), XPOSN, YPOSN)

         YPOSN = BOTTOM + SAFETY
         CALL MESSAG (CLASS, LENTXT (1), XPOSN, YPOSN)
      END IF

C     Write explanation titles in smaller letters. 

      IF (EXPLAN .NE. BLANK) THEN

         CALL HEIGHT (HEXP)

         XPOSN = XCENTER - HALF * XMESS (EXPLAN, LENTXT (2))
         YPOSN = TOP - SAFETY - HTIT - SPACE - HEXP
         CALL MESSAG (EXPLAN, LENTXT (2), XPOSN, YPOSN)

         YPOSN = BOTTOM + SAFETY + HTIT + SPACE
         CALL MESSAG (EXPLAN, LENTXT (2), XPOSN, YPOSN)
      END IF

   99 RETURN
      END
