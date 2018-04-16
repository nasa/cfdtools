C+------------------------------------------------------------------------------
C
      SUBROUTINE LOCASE (STRING)
C
C PURPOSE:  LOCASE changes all upper case letters in the given
C           character string to lower case.
C
C METHOD:   Each character in STRING is treated in turn.  The intrinsic
C           function INDEX effectively allows a table lookup, with
C           the local strings LOW and UPP acting as two tables.
C           This method avoids the use of CHAR and ICHAR, which appear
C           to behave differently on ASCII and EBCDIC machines.
C
C ARGUMENTS:
C    ARG    TYPE I/O/S DESCRIPTION
C  STRING    C   I/O   Character string possibly containing some
C                      uppercase letters on input;
C                      strictly lowercase letters on output with no
C                      change to any non-alphabetic characters.
C
C HISTORY:
C
C  09/10/1985  M. Saunders   Rewrite of version by Hooper/Kennelly, 1983.
C  July 1999   M. Rimlinger  Adaptation of existing UPCASE.
C
C AUTHOR (UPCASE): Michael Saunders, Systems Optimization Lab., Stanford U.
C
C-------------------------------------------------------------------------------

C   Arguments:

      CHARACTER, INTENT (INOUT) :: STRING * (*)

C   Local constants:

      CHARACTER, PARAMETER :: LOW*26 = 'abcdefghijklmnopqrstuvwxyz',
     .                        UPP*26 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
C   Local variables:

      INTEGER    I, J
      CHARACTER  C*1

C   Execution:

      DO J = 1, LEN (STRING)
         C = STRING(J:J)
         IF (C >= 'A' .AND. C <= 'Z') THEN
            I = INDEX (UPP, C)
            IF (I > 0) STRING(J:J) = LOW(I:I)
         END IF
      END DO

      END SUBROUTINE LOCASE
