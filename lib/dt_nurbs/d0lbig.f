      SUBROUTINE D0LBIG (CMEM, IMEM, DMEM, ISRC, IGE)

C     Create an entity label and/or entity subscript for an IGES Entity
C     from a B-spline Function, Edge, Loop, Trimmed Surface or Joined
C     Surface entity's label.  Default behavior in the event of error is
C     to blank the entity label (directory entry 18) and zero the entity
C     subscript (directory entry 19)
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     ISRC    MEM pointer to entity whose character data is a label
C   INPUT-OUTPUT:
C     IGE     MEM pointer to an IGES Entity whose directory entry 18 and 19
C             are to be given a value extracted from the label in ISRC.
C
C   CALLS:
C     D2FEAC
C     D2FELC
C     D2IGSC
C     D2IGSI
C     D2LENC
C     DTERPT
C
C   HISTORY:
C     08Sep93  P Kraushar   Created
C*****
      CHARACTER CMEM*(*)
      INTEGER IMEM(*), ISRC, IGE
      DOUBLE PRECISION DMEM(*)

      CHARACTER STRG*16
      INTEGER LENC, IERX, N, I, J, K
C****
      STRG = ' '
      N = 0
      I = 8
      CALL DTERPT (0)
C     Get length of label
      CALL D2FELC (CMEM, IMEM, DMEM, ISRC, LENC, IERX)
      IF (IERX .NE. 0) GOTO 9000
      IF (LENC .LE. 16) THEN
C        Fetch whole label from source
         CALL D2FEAC (CMEM, IMEM, DMEM, ISRC, 1, LENC, STRG, IERX)
         IF (IERX .NE. 0) GOTO 9000
      ELSE
C        Fetch first eight and last eight characters from source label
         CALL D2FEAC (CMEM, IMEM, DMEM, ISRC, 1, 8, STRG(1:8), IERX)
         IF (IERX .NE. 0) GOTO 9000
         CALL D2FEAC (CMEM, IMEM, DMEM, ISRC, LENC-7, LENC, STRG(9:16),
     +                IERX)
         IF (IERX .NE. 0) GOTO 9000
         LENC = 16
      ENDIF
C     If trailing digits are numeric, translate into entity subscript
      N = 0
      K = 1
      DO 100 I=LENC,LENC-7,-1
         J = INDEX ('0123456789', STRG(I:I)) - 1
         IF (J .LT. 0) GOTO 110
         N = N + K*J
         K = K*10
  100 CONTINUE
      I = LENC - 8
  110 CONTINUE
C     First up-to-8 characters becomes label entity
      IF (I .GT. 8) THEN
         I = 8
      ELSE IF (STRG(I:I) .EQ. '_') THEN
C        Delete any underscore between label and subscript
         I = I - 1
      ENDIF
C   Error handling rejoin point (couldn't get source info)
 9000 CONTINUE
C     Now store the values in IGE, ignoring any errors.
      CALL D2IGSC (CMEM, IMEM, DMEM, 'C', STRG(1:I), IGE, 'D', 18, IERX)
      CALL D2IGSI (CMEM, IMEM, DMEM, 'I', N, IGE, 'D', 19, IERX)
      CALL DTERPT (1)
      RETURN            
      END
