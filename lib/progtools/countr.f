C+----------------------------------------------------------------------
C
      SUBROUTINE COUNTR (LREAD, BUFFER, NUMTOK, OK, EOF)
C
C
C     Description and usage:
C
C           COUNTR determines the number of "non-blank fields" (tokens)
C        present in the first significant line read from the specified
C        logical unit, and then backspaces the file one record so that
C        it is ready for further use.   The delimiters between tokens
C        are taken to be those of subroutine SCANNR, namely blank, tab,
C        comma, colon, or equal sign.  This version skips blank lines
C        and ignores '!' comments through its use of subroutine GETLINE.
C
C           COUNTR was originally written for use with subroutine RDCOLS,
C        where specified columns are to be read in a free-format manner.
C
C     Arguments:
C
C        Name   Dimension  Type  I/O/S  Description
C        LREAD              I    I      Logical unit number for input.
C        BUFFER             C    I      Workspace to hold a line of data.
C        NUMTOK             I      O    Number of non-blank fields found.
C                                       Normally greater than 0 (unless
C                                       EOF encountered).
C        OK                 L      O    Error flag.  Set .TRUE. only if
C                                       the READ and BACKSPACE operations
C                                       were error-free or if EOF = .TRUE.
C        EOF                L      O    End-of-file flag.
C
C
C     Procedures:
C
C        GETLINE  Gets next significant line; suppresses '!' comments.
C        SCANNR   Looks for the position of a non-blank field.
C
C
C     External files:
C
C        LREAD    Input file.
C
C
C     Environment:  Digital VAX/VMS, FORTRAN 77
C                   IMPLICIT NONE is non-standard.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     History:
C
C        31 Dec. 1983   RAK   Initial design and coding.
C        14 Jan. 1984   RAK   Rewritten as a SUBROUTINE.
C        31 Mar. 1984   RAK   Skips blank lines (cf. NLCHEK).
C        19 Jul. 1984   RGL   Added workspace argument BUFFER.
C        26 Dec. 1991   DAS   Installed GETLINE when non-use of
C                             LEN (BUFFER) seemed to be giving trouble.
C                             IOCHEK proved redundant as a result.
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      LOGICAL
     >   OK, EOF
      INTEGER
     >   LREAD, NUMTOK
      CHARACTER
     >   BUFFER * (*)

C     Local variables.

      INTEGER
     >   FIRST, IOS, LAST, MARK

C     Procedures.

      EXTERNAL
     >   GETLINE, SCANNR


C     Execution.
C     ----------

      EOF = .FALSE.
      OK = .TRUE.
   10 CONTINUE

C        Look for a (significant) line of input.

         CALL GETLINE (LREAD, '!', BUFFER, LAST, IOS)
         IF (IOS .LT. 0) THEN
            EOF = .TRUE.
         ELSE IF (IOS .GT. 0) THEN  ! Read error
            OK = .FALSE.
         ELSE IF (LAST .EQ. 0) THEN  ! Insignificant line
            GO TO 10
         ELSE

            FIRST = 1
            NUMTOK = 0

C           Count the number of "non-blank fields" read.

   20       CONTINUE
               CALL SCANNR (BUFFER, FIRST, LAST, MARK)
               IF (FIRST .GT. 0) NUMTOK = NUMTOK + 1
               FIRST = MARK + 2
               IF (FIRST .LE. LAST)
     >      GO TO 20

         END IF

      BACKSPACE (LREAD, IOSTAT = IOS)
      OK = OK .AND. (IOS .EQ. 0)


C     Termination.
C     ------------

      RETURN
      END
