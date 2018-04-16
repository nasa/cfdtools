      SUBROUTINE D2CNEE (CMEM, IMEM, DMEM, IEG1, IEG2, IER)

C     PURPOSE:
C        ConNect Edge to Edge.  Join the two Edges IEG1 and IEG2
C        of a Joined Surface.
C
C        If either edge is zero, the JEG pointer for the other is
C        set to zero.
C
C        If either Edge already has a JEG pointer, an attempt will
C        be made to remove that connection.
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IEG1     Pointer to the first Edge.
C        IEG2     Pointer to the second Edge.
C
C     OUTPUT:
C        IER      Error flag.
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = IEG1 does not point to an Edge entity
C                 -3  = IEG2 does not point to an Edge entity
C                 -999= Unexpected error in a lower level routine
C
C     CALLS:
C        D1MBAD
C        D0EGFP
C        D2STEI
C        DTERR
C
C     HISTORY:
C        09Jul92  D. Parsons  Created.
C        21Jul92  D. Parsons  Added zeroing of old connections.
C
C     ------

C     Long name alias:
C        ENTRY D2_CONNECT_EDGE_EDGE (CMEM, IMEM, DMEM, IEG1, IEG2, IER)

      EXTERNAL D1MBAD
      LOGICAL  D1MBAD

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IEG1, IEG2, IER
      CHARACTER         LABELX
      INTEGER           JEG, IDMY, IERX

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2CNEE'/

C     ------

      IER = 0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

C     Verify the validity of the 1st Edge and disconnect any old
C     connections, if possible

      IF (IEG1 .NE. 0) THEN
         CALL D0EGFP (CMEM, IMEM, DMEM, IEG1, LABELX, IDMY, IDMY,
     +         IDMY, JEG, IDMY, IERX)
         IF (IERX .NE. 0) THEN
            IER = -2
            GOTO 9000
         ENDIF

         IF (JEG .NE. 0) THEN
C           Try to set this third Edge's JEG pointer to zero.
C           1st verify that it is an Edge entity
            CALL D0EGFP (CMEM, IMEM, DMEM, JEG, LABELX, IDMY, IDMY,
     +            IDMY, IDMY, IDMY, IERX)
            IF (IERX .EQ. 0) THEN
C              2nd store 0 into its JEG pointer
               CALL D2STEI (CMEM, IMEM, DMEM, 0, 4, JEG, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -999
                  GOTO 9000
               ENDIF
            ENDIF
         ENDIF
      ENDIF

C     Verify the validity of the 2nd Edge and disconnect any old
C     connections, if possible

      IF (IEG2 .NE. 0) THEN
         CALL D0EGFP (CMEM, IMEM, DMEM, IEG2, LABELX, IDMY, IDMY,
     +         IDMY, JEG, IDMY, IERX)
         IF (IERX .NE. 0) THEN
            IER = -3
            GOTO 9000
         ENDIF

         IF (JEG .NE. 0) THEN
C           Try to set this third Edge's JEG pointer to zero.
C           1st verify that it is an Edge entity
            CALL D0EGFP (CMEM, IMEM, DMEM, JEG, LABELX, IDMY, IDMY,
     +            IDMY, IDMY, IDMY, IERX)
            IF (IERX .EQ. 0) THEN
C              2nd store 0 into its JEG pointer
               CALL D2STEI (CMEM, IMEM, DMEM, 0, 4, JEG, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -999
                  GOTO 9000
               ENDIF
            ENDIF
         ENDIF
      ENDIF

C     Replace the JEG value in IEG1 with IEG2

      IF (IEG1 .NE. 0) THEN
         CALL D2STEI (CMEM, IMEM, DMEM, IEG2, 4, IEG1, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
      ENDIF

C     Replace the JEG value in IEG2 with IEG1

      IF (IEG2 .NE. 0) THEN
         CALL D2STEI (CMEM, IMEM, DMEM, IEG1, 4, IEG2, IERX)
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
