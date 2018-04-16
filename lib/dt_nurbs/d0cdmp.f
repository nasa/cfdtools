      SUBROUTINE D0CDMP (CMEMX, OUTSTR, OUTCTL, SPECAL)

C   Format characters for D1EDMP

C****
      CHARACTER CMEMX*(*), OUTSTR*56, OUTCTL*56
      LOGICAL SPECAL

      INTEGER I, J, K, L, KCH
      CHARACTER CTLCH*33
      CHARACTER*1 C0NML, C0CTL, C1NML, C1CTL, CNOT

      DATA CTLCH / '@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_ ' /
      DATA C0NML, C0CTL, C1NML, C1CTL, CNOT
     +    / ' ', '^', '''', '"', ' ' /
C****
      L = LEN (CMEMX)
      IF (L .GT. 50) STOP 'First argument too long - D0CDMP'
      SPECAL = .FALSE.
      J = 2
      K = 10
      DO 100 I=1,L
        KCH = ICHAR (CMEMX(I:I))
        IF (KCH .LT. 32) THEN
          OUTSTR(J:J) = CTLCH(KCH+1:KCH+1)
          OUTCTL(J:J) = C0CTL
          SPECAL = .TRUE.
        ELSE IF (KCH .LT. 127) THEN
          OUTSTR(J:J) = CMEMX(I:I)
          OUTCTL(J:J) = C0NML
        ELSE IF (KCH .LT. 128) THEN
          OUTSTR(J:J) = CTLCH(33:33)
          OUTCTL(J:J) = C0CTL
          SPECAL = .TRUE.
        ELSE IF (KCH .LT. 160) THEN
          OUTSTR(J:J) = CTLCH(KCH-127:KCH-127)
          OUTCTL(J:J) = C1CTL
          SPECAL = .TRUE.
        ELSE IF (KCH .LT. 255) THEN
          OUTSTR(J:J) = CHAR(KCH-128)
          OUTCTL(J:J) = C1NML
          SPECAL = .TRUE.
        ELSE
          OUTSTR(J:J) = CTLCH(33:33)
          OUTCTL(J:J) = C1CTL
          SPECAL = .TRUE.
        END IF
        J = J + 1
        K = K - 1
        IF (K .EQ. 0) THEN
          K = 10
          J = J + 1
        END IF
  100 CONTINUE
      DO 110 I=L+1,50
        OUTSTR(J:J) = CNOT
        OUTCTL(J:J) = ' '
        J = J + 1
        K = K - 1
        IF (K .EQ. 0) THEN
          K = 10
          J = J + 1
        END IF
  110 CONTINUE
      RETURN
      END
