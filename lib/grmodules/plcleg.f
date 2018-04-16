C+**********************************************************************
C
      SUBROUTINE PLCLEG ( NLEVS, RLEVS, NPRLEV, NDGLBL,
     +                    XLEGND, YLEGND, TXLGND, HTLGND )
C
C  PURPOSE:  PLCLEG produces the sort of legend needed by contour plots.
C	     See PLLEVS for its usage.
C
C  ARGUMENTS:
C    NAME   TYPE  I/O/S   DIM    DESCRIPTION
C    NLEVS   I      I      -     Number of contour levels.
C    RLEVS   R      I    NLEVS   Actual contour levels.
C    NPRLEV  I      I    NLEVS   Number of contour lines per level.
C    NDGLBL  I      I      -     Number of digits for contour labels.
C    XLEGND  R      I      -     Horizontal position of legend.
C    YLEGND  R      I      -     Vertical position of legend.
C    TXLGND  C*10   I      -     Legend text ( 10 characters ).
C    HTLGND  R      I      -     Height of legend and labels.
C                                0.0 means take default of 0.07.
C  EXTERNAL REFERENCES:
C     NAME    DESCRIPTION  AND  SOURCE
C    ANGLE  - For ensuring horizontal legend lines.           DISSPLA
C    HEIGHT - Sets height of text.                            DISSPLA
C    INTNO  - Writes integer number (from physical origin).   DISSPLA
C    MESSAG - Writes message (from physical origin).          DISSPLA
C    REALNO - Writes real number (from physical origin).      DISSPLA
C    XMESS  - Outputs length of message.                      DISSPLA
C
C  STANDARDS VIOLATIONS: None.
C
C  ENVIRONMENT:  VAX/VMS, CRAY-1S -- FORTRAN 77
C
C  DEVELOPMENT HISTORY:
C    DATE     PERSON    STATAMENT OF CHANGES
C  09/01/89    MDW      Updated to run on DISSPLA version 11.0.  
C  11/17/82    RGL      Introduced FORTRAN 77, esp. CHARACTER variables.
C  01/09/81    CJG      Original coding.
C
C  AUTHOR: Carol Green, Ron Langhi, Informatics Inc.
C
C-**********************************************************************
C
      CHARACTER  TXLGND*10
      DIMENSION  RLEVS(NLEVS), NPRLEV(NLEVS)
      LOGICAL    NOLINE
C
C  *  Print legend headings:
C
      HLEG = HTLGND
      IF ( HLEG.LE.0. ) HLEG = .07
      CALL HEIGHT ( HLEG )
      CALL ANGLE ( 0.0 )
      CALL MESSAG ( 7HCONTOUR, 7, XLEGND, YLEGND )
      ALENG = XMESS( 7HCONTOUR, 7 )
      XPOS  = XLEGND + 1.2 * ALENG
      CALL MESSAG ( TXLGND, 10, XPOS, YLEGND )
      YINC = HLEG + HLEG
      YPOS = YLEGND
C
C  *  Fill in legend values:
C
      NOLINE = .FALSE.
C
      DO 10 LEVEL = 1, NLEVS
         YPOS = YPOS - YINC
         CALL INTNO ( LEVEL, XLEGND, YPOS )
         CALL REALNO ( RLEVS( LEVEL ), NDGLBL, XPOS, YPOS )
         IF ( NPRLEV( LEVEL ).EQ.0 ) THEN
            X = XPOS + ALENG / 7. * ( NDGLBL + 3 )
            CALL MESSAG( ' *', 2, X, YPOS )
            NOLINE = .TRUE.
         END IF
   10 CONTINUE
C
      IF ( NOLINE ) THEN
         YPOS = YPOS - 1.5 * YINC
         CALL MESSAG ( '*  NO CONTOUR LINES FOUND', 25, XLEGND, YPOS )
      END IF
C
      RETURN
      END
