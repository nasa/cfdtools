      INTEGER FUNCTION D0TRMC (STRING)

C     Very simple routine to find the last non-blank character
C     in a string.
C
C     INPUT
C        STRING
C
C     OUTPUT
C        D0TRMC   Position of the last non--blank character in the
C                 string.
C
C     HISTORY
C        20May92     D. Parsons  Created
C

      CHARACTER*(*) STRING
      INTEGER I

      DO 10 I = LEN(STRING), 1, -1
         IF (STRING(I:I) .NE. ' ') THEN
            D0TRMC = I
            GOTO 20
         ENDIF
  10  CONTINUE

      D0TRMC = 0

  20  CONTINUE

      RETURN
      END
