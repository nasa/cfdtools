C+----------------------------------------------------------------------
      SUBROUTINE MAPFIX (WXLEFT, WXRGHT, WYBOT, WYTOP, VXLFTD, VXRGTD,
     1                   VYBOTD, VYTOPD, VXLFT, VXRGT, VYBOT, VYTOP,IER)
C  PURPOSE -
C    Map a window (user, data space) onto all or part of a viewport
C    (program, plot space) such that the shape of any data, plotted after
C    the window-viewport mapping is established, is preserved.
C    The actual viewport returned by this routine will be centered in
C    the viewport passed to this routine.  If necessary, the viewport
C    will be shortened on one axis.  This routine does not make any
C    graphics calls.  It only sets viewport values for the caller.
C
C  PARAMETERS -
C    ARG  TYPE I/O/S  DIM         DESCRIPTION
C   WXLEFT R     I     -       The X-values on the sides of the window
C   WXRGHT R     I     -	(WXLEFT may be greater or less than WXRGHT).
C   WYBOT  R     I     -       The Y-values at bottom/top of the window
C   WYTOP  R     I     -	(WYBOT  may be greater or less than WYTOP).
C   VXLFTD R     I     -       Default X-values for left/right edges of
C   VXRGTD R     I     -         the viewport (VXRGTD>VXLFTD).
C   VYBOTD R     I     -       Default Y-values for bottom/top edges of
C   VYTOPD R     I     -         the viewport (VYTOPD>VYBOTD).
C   VXLFT  R     O     -       Actual  X-values for left/right edges of
C   VXRGT  R     O     -         the viewport (may be adjusted closer).
C   VYBOT  R     O     -       Actual  Y-values for bottom/top edges of
C   VYTOP  R     O     -         the viewport (may be adjusted closer).
C   IER    I     O     -       Error return code
C				0 = OK
C				1 = window is degenerate (e.g., WYBOT=WYTOP).
C				2 = default viewport is degenerate or inverted
C				    (e.g., VYBOTD>=VYTOPD).
C
C  COMMONS USED -
C    VAR  TYPE I/O/S  DIM   OFFSET  DESCRIPTION
C   (None)
C  EXTERNAL REFERENCES -
C    (None)
C
C  AUTHOR - Dexter L. Hermstad		Informatics General Corporation
C		01/11/84		Palo Alto, California
C
C  REVISIONS -		LAST MODIFIED:  DLH 01/11/84
C    DATE      PERSON  STATEMENT OF CHANGES
C    01/16/84	DLH	Allows reverse direction in either window axis.
C    01/11/84	DLH	Original coding
C    01/11/84	DLH	Design
C
C  NOTES -
C-----------------------------------------------------------------------
C
C     ... Check for a zero length window boundary
      IF (WXLEFT.EQ.WXRGHT .OR. WYBOT.EQ.WYTOP) THEN
         IER= 1
         GO TO 800
      ENDIF
C
C     ... Check for a zero or negative length viewport boundary
      IF (VXLFTD.GE.VXRGTD .OR. VYBOTD.GE.VYTOPD) THEN
         IER= 2
         GO TO 800
      ENDIF
C
      IER= 0
C
C     ... Compute aspect ratio of window (subject space or data)
      RATWND= (WYTOP-  WYBOT)/ (WXRGHT-WXLEFT)
C
C     ... Compute aspect ratio of viewport (object space or plot area)
      RATVWP= (VYTOPD-VYBOTD)/ (VXRGTD-VXLFTD)
C
C     ... See which viewport axis (if any) needs to be shortened to
C     ... make the aspect ratios equal.  ABS allows the software to
C     ... return correct values even if one or both sides of the window
C     ... go from higher to lower.
      IF (ABS(RATWND/RATVWP)-1.) 100, 200, 300
C
C     ... Y axis needs to be shortened
  100 CONTINUE
      VXLFT= VXLFTD
      VXRGT= VXRGTD
      YDIF=  ((VYTOPD-VYBOTD)-((VXRGTD-VXLFTD)*RATWND))*0.5
      VYBOT= VYBOTD+ YDIF
      VYTOP= VYTOPD- YDIF
C     ... The viewport values must be ascending - this correction is
C     ... necessary if the window values were descending
      IF (VYBOT.GT.VYTOP) THEN
         TEMP= VYBOT
         VYBOT= VYTOP
         VYTOP= TEMP
      ENDIF
      GO TO 800
C
C     ... No adjustment necessary
  200 CONTINUE
      VXLFT= VXLFTD
      VXRGT= VXRGTD
      VYBOT= VYBOTD
      VYTOP= VYTOPD
      GO TO 800
C
C     ... X axis needs to be shortened
  300 CONTINUE
      XDIF=  ((VXRGTD-VXLFTD)-((VYTOPD-VYBOTD)/RATWND))*0.5
      VXLFT= VXLFTD+ XDIF
      VXRGT= VXRGTD- XDIF
      VYBOT= VYBOTD
      VYTOP= VYTOPD
C     ... The viewport values must be ascending - this correction is
C     ... necessary if the window values were descending
      IF (VXLFT.GT.VXRGT) THEN
         TEMP= VXLFT
         VXLFT= VXRGT
         VXRGT= TEMP
      ENDIF
C
  800 CONTINUE
      RETURN
      END
