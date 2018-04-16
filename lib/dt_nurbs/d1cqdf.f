      SUBROUTINE D1CQDF (CMEM, IMEM, DMEM, NSTR, NCHR, INSTR, ICQ, IER)

C   Shadow routine for D2CQDF.
C
C   Define a character string sequence entity.  The space to allocate
C   may be given explicitly using NSTR and NCHR, or implicitly via the
C   initialization string INSTR.  INSTR contains a sequence of strings
C   marked off by delimiter characters.  Any character not present in
C   any of the strings may be used as a delimiter.  The subroutine
C   uses the first character in INSTR as the string delimiter.
C
C        Character String Sequence Entity:
C
C           CMEM  [1..len()]  Strings (with delimiters removed)
C
C           DMEM  (not used)
C
C           IMEM  [1]         Number of strings
C                 [2..1+No.]  Ending position of each string
C
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     NSTR    Number of items in sequence.  The integer space allocated
C             is the maximum of this number and the number of strings
C             found in INSTR.
C     NCHR    Length of character space to allocate.  The character
C             space allocated is the maximum of this number and the
C             number of non-delimiter characters in INSTR.
C     INSTR   Initialization string.  First character is delimiter
C             used to separate successive strings.  Last character
C             should also be delimiter, but may be omitted.
C   OUTPUT:
C     ICQ     Pointer to character string sequence (CQ) entity
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             -2  = Insufficient space available for character data
C             -3  = Insufficient space available for integer data
C
C   CALLS:    D1DEFE
C   HISTORY:
C     25May93 D Parsons     Created shadow routine from D2CQDF.  Error
C                           checking should be done by D2CQDF.
C
C*****

      EXTERNAL D0TRMC
      INTEGER  D0TRMC

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)

      INTEGER NSTR, NCHR, ICQ, IER
      CHARACTER INSTR*(*)

      INTEGER ENTCQ
      PARAMETER (ENTCQ=249)

      INTEGER NST, NCH, NIN, NSTRA, NCHRA, JB, JE, I
      INTEGER IHDR, IXC0, IXCB, IXCE, IXI1
      CHARACTER DELIM*1
C****
      IER = 0

      NIN = D0TRMC (INSTR)

C     Find number of strings and non-delimiter characters in INSTR
      NST = 0
      IF (NIN .GT. 1) THEN
        DELIM = INSTR(1:1)
        DO 100 I=2,NIN
          IF (INSTR(I:I) .EQ. DELIM) NST = NST + 1
  100   CONTINUE
        NCH = NIN - NST - 1
        IF (INSTR(NIN:NIN) .NE. DELIM) NST = NST + 1
      ELSE
        NCH = 0
      END IF

C     Allocate the new entity
      NSTRA = 1 + MAX (NSTR, NST)
      NCHRA = MAX (NCHR, NCH)
      CALL D1DEFE (CMEM, IMEM, DMEM, ENTCQ, NCHRA, NSTRA, 0, .FALSE.,
     +    ICQ, IER)
      IF (IER .NE. 0) GOTO 9900

C     Initialize it from INSTR
      IHDR = ICQ/256
      IXC0 = IMEM(IHDR+1) - 1
      IXI1 = IMEM(IHDR+2)
C     Store current number of strings as first integer datum
      IMEM(IXI1) = NST
C     Extract strings from INSTR and store in CMEM, saving ending character
C     (relative) positions in successive locations of IMEM
      JE = 1
      IXCE = IXC0
      DO 200 I=1,NST
        JB = JE + 1
C       Start search loop for end of next string
  150     CONTINUE
          JE = JE + 1
          IF (JE .GT. NIN) GO TO 160
          IF (INSTR(JE:JE) .NE. DELIM) GO TO 150
C       End loop
  160   CONTINUE
        IF (JE .GT. JB) THEN
          IXCB = IXCE + 1
          IXCE = IXCE + JE - JB
          CMEM(IXCB:IXCE) = INSTR(JB:JE-1)
        END IF
        IMEM(IXI1+I) = IXCE - IXC0
  200 CONTINUE

C     Exit normally

      IER = 0

C   (Any error reporting should be handled by D2CQDF)

 9900 CONTINUE

      RETURN

      END
