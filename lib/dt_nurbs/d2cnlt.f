      SUBROUTINE D2CNLT (CMEM, IMEM, DMEM, IXLP, ILP, ITS, IER)

C     PURPOSE:
C        ConNect Loop to Trimmed Surface.  Record Loop ILP as the IXLPth
C        Loop in Trimmed Surface ITS.
C
C        If ILP is zero, the ITS(IXLP) is set to zero.
C
C        If ILP has a non-zero Trimmed Surface pointer, or ITS already
C        has a Loop defined in its IXLPth index, an attempt will be made
C        to remove those connections.
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IXLP     Index value of the Loop.
C        ILP      Pointer to the Loop.
C        ITS      Pointer to the Trimmed Surface.
C
C     OUTPUT:
C        IER      Error flag.
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = ITS does not point to a Trimmed Surface entity
C                 -3  = IXLP not between 1 and MAXLP
C                 -4  = ILP does not point to a Loop entity
C                 -999= Unexpected error in a lower level routine
C
C     CALLS:
C        D1MBAD
C        D0LPFP
C        D0TSFP
C        D2FEEI
C        D2STEI
C        DTERR
C
C     HISTORY:
C        09Jul92  D. Parsons  Created.
C        21Jul92  D. Parsons  Added zeroing of old connections.
C
C     ------

C     Long name alias:
C        ENTRY D2_CONNECT_LOOP_TRIMMED (CMEM, IMEM, DMEM, IXLP, ILP,
C    +         ITS, IER)

      EXTERNAL D1MBAD
      LOGICAL  D1MBAD

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IXLP, ILP, ITS, IER
      CHARACTER         LABELX
      INTEGER           IDMY, MAXLP, ILPX, IXLPX, ITSX, IERX

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2CNLT'/

C     ------

      IER = 0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

C     Verify the validity of the Trimmed Surface and get MAXLP.  Try
C     to disconnect any old connections

      IF (ITS .NE. 0) THEN
         CALL D0TSFP (CMEM, IMEM, DMEM, ITS, LABELX, IDMY, IDMY, IDMY,
     +         IDMY, MAXLP, IDMY, IERX)
         IF (IERX .NE. 0) THEN
            IER = -2
            GOTO 9000
         ENDIF
         IF ((IXLP .LT. 1) .OR. (IXLP .GT. MAXLP)) THEN
            IER = -3
            GOTO 9000
         ENDIF

C        Look for old connections

         CALL D2FEEI (CMEM, IMEM, DMEM, ITS, 5+IXLP, ILPX, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
         IF (ILPX .NE. 0) THEN
            CALL D0LPFP (CMEM, IMEM, DMEM, ILPX, LABELX, ITSX, IDMY,
     +            IDMY, IDMY, IERX)
            IF (ITSX .EQ. ITS) THEN
C              Zero the Trimmed Surface pointer
               CALL D2STEI (CMEM, IMEM, DMEM, 0, 1, ILPX, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -999
                  GOTO 9000
               ENDIF
C              Zero the index
               CALL D2STEI (CMEM, IMEM, DMEM, 0, 2, ILPX, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -999
                  GOTO 9000
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      IF (ILP .NE. 0) THEN

C        Verify the validity of the Loop and get its ITS pointer

         CALL D0LPFP (CMEM, IMEM, DMEM, ILP, LABELX, ITSX, IXLPX,
     +      IDMY, IDMY, IERX)
         IF (IERX .NE. 0) THEN
            IER = -4
            GOTO 9000
         ENDIF

C        If the Loop is already pointing to a Trimmed Surface,
C        try to disconnect the link

         IF (ITSX .NE. 0) THEN
            CALL D0TSFP (CMEM, IMEM, DMEM, ITSX, LABELX, IDMY, IDMY,
     +            IDMY, IDMY, IDMY, IDMY, IERX)
            IF (IERX .EQ. 0) THEN
               CALL D2FEEI (CMEM, IMEM, DMEM, ITSX, 5+IXLPX, ILPX, IERX)
               IF (ILPX .EQ. ILP)
     +               CALL D2STEI (CMEM, IMEM, DMEM, 0, 5+IXLPX,
     +               ITSX, IERX)
            ENDIF
         ENDIF
      ENDIF


C     Replace the IXLPth Loop in ITS with ILP

      IF (ITS .NE. 0) THEN
         CALL D2STEI (CMEM, IMEM, DMEM, ILP, 5+IXLP, ITS, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
      ENDIF

C     Set the Loop to point back to this ITS

      IF (ILP .NE. 0) THEN
C        Store the Trimmed Surface pointer into ILP
         CALL D2STEI (CMEM, IMEM, DMEM, ITS, 1, ILP, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
C        Store the index IXLP pointer into ILP
         IF (ITS .NE. 0) THEN
            CALL D2STEI (CMEM, IMEM, DMEM, IXLP, 2, ILP, IERX)
         ELSE
C           Zero the index pointer when zeroing the ITS
            CALL D2STEI (CMEM, IMEM, DMEM, 0, 2, ILP, IERX)
         ENDIF
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF

      ENDIF

      GOTO 9999

 9000 CONTINUE

C     Error Handling

      IF (IER .EQ. -999) THEN
         CALL DTERR (5, SUBNAM, IER, 0)
      ELSE
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

 9999 CONTINUE

      RETURN
      END
