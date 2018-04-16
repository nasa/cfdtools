      SUBROUTINE D2CNEL (CMEM, IMEM, DMEM, IXEG, IEG, ILP, IER)

C     PURPOSE:
C        ConNect Edge to Loop.  Record Edge IEG as the IXEGth
C        Edge in Loop ILP.
C
C        If IEG is zero, the ILP(IXEG) is set to zero.
C        If ILP is zero, the IEG pointer to ILP is set to zero.
C
C        If IEG has a non-zero Loop pointer, or ILP already has an
C        Edge defined in its IXEGth index, an attempt will be made
C        to remove those connections.
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IXEG     Index value of the Edge.
C        IEG      Pointer to the Edge.
C        ILP      Pointer to the Loop.
C
C     OUTPUT:
C        IER      Error flag.
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = ILP does not point to a Loop entity
C                 -3  = IXEG not between 1 and MAXEG
C                 -4  = IEG does not point to an Edge entity
C                 -999= Unexpected error in a lower level routine
C
C     CALLS:
C        D1MBAD
C        D0EGFP
C        D0LPFP
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
C        ENTRY D2_CONNECT_EDGE_LOOP (CMEM, IMEM, DMEM, IXEG, IEG,
C    +         ILP, IER)

      EXTERNAL D1MBAD
      LOGICAL  D1MBAD

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IXEG, IEG, ILP, IER
      CHARACTER         LABELX
      INTEGER           IDMY, MAXEG, IEGX, ILPX, IXEGX, IERX

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2CNEL'/

C     ------

      IER = 0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

C     Verify the validity of the Loop and get MAXEG.  Try to
C     disconnect any old connections

      IF (ILP .NE. 0) THEN
         CALL D0LPFP (CMEM, IMEM, DMEM, ILP, LABELX, IDMY, IDMY,
     +         MAXEG, IDMY, IERX)
         IF (IERX .NE. 0) THEN
            IER = -2
            GOTO 9000
         ENDIF
         IF ((IXEG .LT. 1) .OR. (IXEG .GT. MAXEG)) THEN
            IER = -3
            GOTO 9000
         ENDIF

C        Look for old connections

         CALL D2FEEI (CMEM, IMEM, DMEM, ILP, 3+IXEG, IEGX, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
         IF (IEGX .NE. 0) THEN
            CALL D0EGFP (CMEM, IMEM, DMEM, IEGX, LABELX, IDMY, ILPX,
     +            IDMY, IDMY, IDMY, IERX)
            IF (ILPX .EQ. ILP) THEN
C              Zero the Loop pointer
               CALL D2STEI (CMEM, IMEM, DMEM, 0, 2, IEGX, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -999
                  GOTO 9000
               ENDIF
C              Zero the index
               CALL D2STEI (CMEM, IMEM, DMEM, 0, 3, IEGX, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -999
                  GOTO 9000
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      IF (IEG .NE. 0) THEN

C        Verify the validity of the Edge and get its ILP pointer

         CALL D0EGFP (CMEM, IMEM, DMEM, IEG, LABELX, IDMY, ILPX,
     +         IXEGX, IDMY, IDMY, IERX)
         IF (IERX .NE. 0) THEN
            IER = -4
            GOTO 9000
         ENDIF

C        If the Edge is already pointing to a loop, try to disconnect
C        the link

         IF (ILPX .NE. 0) THEN
            CALL D0LPFP (CMEM, IMEM, DMEM, ILPX, LABELX, IDMY, IDMY,
     +            IDMY, IDMY, IERX)
            IF (IERX .EQ. 0) THEN
               CALL D2FEEI (CMEM, IMEM, DMEM, ILPX, 3+IXEGX, IEGX, IERX)
               IF (IEGX .EQ. IEG)
     +               CALL D2STEI (CMEM, IMEM, DMEM, 0, 3+IXEGX,
     +               ILPX, IERX)
            ENDIF
         ENDIF
      ENDIF

C     Replace the IXEGth Edge in ILP with IEG

      IF (ILP .NE. 0) THEN
         CALL D2STEI (CMEM, IMEM, DMEM, IEG, 3+IXEG, ILP, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
      ENDIF

C     Set the Edge IEG to point back to this ILP

      IF (IEG .NE. 0) THEN
C        Store the Loop pointer into IEG
         CALL D2STEI (CMEM, IMEM, DMEM, ILP, 2, IEG, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
C        Store the index pointer IXEG into IEG
         IF (ILP .NE. 0) THEN
            CALL D2STEI (CMEM, IMEM, DMEM, IXEG, 3, IEG, IERX)
         ELSE
C           Zero the index pointer when zeroing the ILP
            CALL D2STEI (CMEM, IMEM, DMEM, 0, 3, IEG, IERX)
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


