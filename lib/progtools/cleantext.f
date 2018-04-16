C+----------------------------------------------------------------------
C
      SUBROUTINE CLEANTEXT (STRING1, STRING2, LAST)
C
C  One-liner:  Remove line feeds, commas, etc., from list-directed output
C
C  Description and usage:
C
C        CLEANTEXT was prompted by the fact that list-directed output of
C     long I/O lists appears to split them at 132-column boundaries.  For
C     multi-column plottable data, this may be intolerable.  Formatted
C     output could overcome the problem, but the reduced bulk of list-
C     directed I/O (on some systems) may be preferable.
C
C        This utility simply strips out all likely separators and repacks
C     the tokens with just one blank between each pair.  The repacking
C     may be done in place if the same string is used for both arguments.
C
C  Environment:  Fortran 90
C
C  History:
C
C     10/03/01  DAS  Initial implementation prompted by the need to
C                    regularize the steps in multi-column trajectory
C                    time histories.
C
C  Author:  David Saunders, ELORET/NASA Ames Research Center
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      CHARACTER, INTENT (INOUT) ::
     >   STRING1 * (*), ! Input string; STRING1 (1 : LAST) is processed
     >   STRING2 * (*)  ! Output string with single blanks between tokens,
                        ! ending at character LAST (appropriately updated);
                        ! the same variable may be passed as STRING1 & -2

      INTEGER,   INTENT (INOUT) ::
     >   LAST           ! Input as the last significant character of STRING1;
                        ! output ................................... STRING2

C     Constants:

      CHARACTER, PARAMETER ::
     >   BLANK * 1 = ' '

C     Local variables:

      INTEGER
     >   FINAL, FIRST, MARK, START

      CHARACTER
     >   DELIMITERS * 7

C     Procedures:

      EXTERNAL
     >   SCAN2 ! Token-finding utility

C     Execution:

C     Set up the delimiters for SCAN2:

      DELIMITERS (1 : 1) = BLANK
      DELIMITERS (2 : 2) = ','
      DELIMITERS (3 : 3) = CHAR (8)   ! Back space
      DELIMITERS (4 : 4) = CHAR (9)   ! Horizontal tab
      DELIMITERS (5 : 5) = CHAR (10)  ! Line feed
      DELIMITERS (6 : 6) = CHAR (13)  ! Carriage return
      DELIMITERS (7 : 7) = CHAR (127) ! Delete

      FIRST = 1
      START = 1 ! Start of next token to be placed in STRING2
      FINAL = 0 ! End of it (initialized in case STRING1 is empty)

      DO ! Until end of STRING1, look for tokens to transcribe

         CALL SCAN2 (STRING1, DELIMITERS, FIRST, LAST, MARK)

         IF (MARK == 0) EXIT

         FINAL = START + MARK - FIRST ! End of this token within STRING2

         STRING2 (START : FINAL) = STRING1 (FIRST : MARK)

         START = FINAL + 1 ! For the single space, but avoid a trailing blank

         IF (START < LAST) STRING2 (START : START) = BLANK

         START = START + 1
         FIRST = MARK  + 2

      END DO

      LAST = FINAL ! In STRING2

      END SUBROUTINE CLEANTEXT
