      SUBROUTINE D0FITS (STRGIN, STRGOU, LENOU, IER)

C   Copy a string to another, truncating or padding on the right with
C   spaces to make it fit, as usual, but indicate the "used" length of
C   the output string in LENOU.  When truncating is necessary, report the
C   true length of the input string as a positive value in IER.
C
C   INPUT:
C     STRGIN  Input string.
C   OUTPUT:
C     STRGOU  Output string variable.
C     LENOU   Length used in STRGOU to hold input from STRGIN.  Remainder,
C             if any, of STRGOU represents padding.
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             a positive value indicates truncation occurred, and gives
C             the true length of the input string.  Negative values are
C             not generated.
C
C   CALLS:    none
C   HISTORY:
C     02Jun92 P Kraushar    Created.
C*****
C
      CHARACTER STRGIN*(*), STRGOU*(*)
      INTEGER LENOU, IER
C****
      STRGOU = STRGIN
      LENOU = LEN (STRGIN)
      IF (LENOU .LE. LEN (STRGOU)) THEN
        IER = 0
      ELSE
        IER = LENOU
        LENOU = LEN (STRGOU)
      END IF
      RETURN
      END
