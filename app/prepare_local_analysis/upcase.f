C+----------------------------------------------------------------------
C
      SUBROUTINE UPCASE( STRING )
C
C PURPOSE:  UPCASE changes all lower case letters in the given
C           character string to upper case.
C
C METHOD:   Each character in STRING is treated in turn.  The intrinsic
C           function INDEX effectively allows a table lookup, with
C           the local strings LOW and UPP acting as two tables.
C           This method avoids the use of CHAR and ICHAR, which appear
C           to behave differently on ASCII and EBCDIC machines.
C
C ARGUMENTS
C    ARG       DIM     TYPE I/O/S DESCRIPTION
C  STRING       *       C   I/O   Character string possibly containing
C                                 some lower-case letters on input;
C                                 strictly upper-case letters on output
C                                 with no change to any non-alphabetic
C                                 characters.
C
C EXTERNAL REFERENCES:
C  LEN    - Returns the declared length of a CHARACTER variable.
C  INDEX  - Returns the position of second string within first.
C
C ENVIRONMENT:  ANSI FORTRAN 77
C
C AUTHOR: Michael Saunders, Systems Optimization Lab., Stanford, 9/10/85
C         (rewrite of version by Hooper/Kennelly, Informatics, 1983)
C
C-----------------------------------------------------------------------

      CHARACTER      STRING * (*)
      INTEGER        I, J
      CHARACTER      C*1, LOW*26, UPP*26
      DATA           LOW /'abcdefghijklmnopqrstuvwxyz'/,
     $               UPP /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      DO 10 J = 1, LEN(STRING)
         C = STRING(J:J)
         IF (C .GE. 'a'  .AND.  C .LE. 'z') THEN
            I = INDEX( LOW, C )
            IF (I .GT. 0) STRING(J:J) = UPP(I:I)
         END IF
   10 CONTINUE

      RETURN
      END
