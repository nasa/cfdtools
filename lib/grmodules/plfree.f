C+----------------------------------------------------------------------
C
      SUBROUTINE PLFREE ( NFREE, FRETXT, FREVAL, NDGFRE, HITFRE, YTXT,
     +                    PAGEX )
C
C  PURPOSE:  PLFREE plots (horizontally centered) tabulated "free"
C            values with headings at the specified vertical position.
C            It assumes prior <sub>plot initialization (page size
C            and physical origin). Text is printed in the current
C            character font.
C
C  ARGUMENTS:
C    VAR    TYPE I/O/S DIM   DESCRIPTION
C    NFREE   I     I    -    Number of free values to be tabulated.
C                            NFREE must be greater than 0, and 10 is
C                            the usual upper limit.
C    FRETXT  C*10  I  NFREE  Descriptive text for free values.
C    FREVAL  R     I  NFREE  Free values.
C    NDGFRE  I     I  NFREE  >0 means floating point with NDGFRE decimal
C                             places after decimal point.
C                            >100 means free point, NDGFRE-100 digits total.
C                            =0 means INTEGER.
C                            <0 means exponent form with NDGFRE decimal
C                             places after decimal point.
C    HITFRE  R     I    -    Approximate height (in inches) of free
C                            values and free text. 0.0 means take default.
C    YTXT    R     I    -    Vertical plot position for FRETXT headings.
C                            Free values appear below FRETXT.
C    PAGEX   R     I    -    Specifies length (inches) measured from
C                            current physical origin of free text/values
C                            display.
C
C  ENVIRONMENT:  VAX/VMS, CRAY-1 -- FORTRAN-77
C
C  HISTORY:
C  09/01/89    MDW    Updated to run on DISSPLA version 11.0.
C  12/15/82    RCL    Expanded use of NDGFRE.
C  11/17/82    RGL    Introduced CHARACTER variables.
C  01/04/82    PJT    Original design and coding.
C
C  AUTHOR:  Jeff Trosin, NASA Ames/Informatics, Inc.
C
C-----------------------------------------------------------------------

      CHARACTER  FRETXT( NFREE ) * 10

      DIMENSION  FREVAL( NFREE ), NDGFRE( NFREE )


      HFREE = HITFRE
      IF ( HFREE.LE.0. ) HFREE = .07
      CALL HEIGHT ( HFREE )
      XLEN   = XMESS ( '     ', 5 )
      PAGELN = PAGEX - 0.3
      XINC   = PAGELN / NFREE
      XTXT   = PAGELN*0.5 - XINC*( NFREE - 1 )/2 - XLEN
      XVAL   = XTXT + 0.3*XLEN
      YVAL   = YTXT - 1.5*HFREE

      DO 400 J = 1, NFREE
         CALL MESSAG ( FRETXT(J), 10, XTXT, YTXT )
         IF ( NDGFRE(J).NE.0 ) THEN
            CALL REALNO ( FREVAL(J), NDGFRE(J), XVAL, YVAL )
         ELSE
            IVAL = INT ( FREVAL(J) + SIGN ( 0.5, FREVAL(J) ) )
            CALL INTNO ( IVAL, XVAL + 0.25, YVAL )
         END IF
         XTXT = XTXT + XINC
         XVAL = XVAL + XINC
 400  CONTINUE

      CALL RESET ( 'HEIGHT' )

      RETURN
      END
