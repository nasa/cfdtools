C+----------------------------------------------------------------------
C
      SUBROUTINE GETLINE (LUN, COMMENT, LINE, LAST, IOS)
C
C  One-liner:  Low-level data reader; suppresses trailing comments/blanks
C
C  Description and usage:
C
C        GETLINE is intended to be a standard low level text input utility,
C     providing a uniform input format which permits free use of blank
C     lines, comment lines, trailing comments, and "commented-out" lines.
C     It reads a record from logical unit LUN and returns the "significant"
C     portion in LINE, with only trailing COMMENTs, blanks, or tabs after
C     LAST.  It returns LAST = 0 if the line is effectively empty.
C
C        Double-COMMENTs are replaced with single ones so that COMMENT may
C     still be used in a string if required.  SPECIAL CASE: if COMMENT is
C     the FIRST significant character, the line is considered empty.  This
C     covers the common case of "commenting-out" lines with more than one
C     COMMENT character, as in  !!!! 0.800000   1.
C
C  Arguments:
C
C     Name    Type/Dimension  I/O/S  Description
C
C     LUN     I               I      Logical unit for reading text.
C
C     COMMENT C*1             I      Character signalling the end of the
C                                    "significant" part of a text string
C                                    and the beginning of comments. If
C                                    blank, this feature is disabled and
C                                    much less scanning is needed.
C
C     LINE    C*(*)             O/S  Buffer for reading one record;
C                                    returned with "significant" text in
C                                    LINE (1 : LAST) unless LAST is zero.
C
C     LAST    I                 O    Index of last significant character
C                                    in LINE.  If LINE is null or all
C                                    blanks, tabs, and/or comment, LAST
C                                    is returned as zero.
C
C     IOS     I                 O    Error status of read;  = 0 means no
C                                    read error, negative is EOF, and the
C                                    meaning of positive IOS is system-
C                                    dependent.
C
C  Environment:  Fortran 90
C
C  Notes:
C
C     (1)  Note that GETLINE does NOT keep reading if it encounters an
C          empty line.  Part of the purpose of GETLINE is to avoid this
C          very drawback of FORTRAN's list-directed I/O.
C
C  Authors:  Ronald Langhi/Robert Kennelly,  Sterling Software, Palo Alto, CA
C
C  History:
C
C     10 Mar. 1987  RGL/RAK  Initial design and code as GETSTRING.
C      4 Apr. 1987    RAK    Cosmetics, documentation revised.
C      8 May  1987    DAS    Name changed to GETLINE.
C     22 Apr. 1988  DAS/RAK  Trim trailing tabs as well as the blanks.
C      7 July 1990  DAS/RAK  Handle special case of commenting-out lines
C                            with more than one leading COMMENT.
C     19 May  1999    DAS    Fortran 90's SIZE=N replaces nonstandard
C                            '(Q, A)' read for counting # characters.
C     20 May  1999    DAS    NO!  SIZE=... raises EOR/EOF issues that
C                            appear insoluble in a portable way.
C                            Use LEN_TRIM (LINE) and avoid SIZE=...
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Constants.

      CHARACTER
     >   BLANK * 1, HTAB * 1
      PARAMETER
     >   (BLANK = ' ',
     >    HTAB  = CHAR (9))   ! HTAB = 9 for Absoft FORTRAN on 68000 machines

C     Arguments.

      INTEGER
     >   LUN, LAST, IOS
      CHARACTER
     >   COMMENT * 1, LINE * (*)

C     Local variables.

      INTEGER
     >   I, J, N              ! N can be LAST throughout, but may be inefficient
      LOGICAL
     >   SUPPRESD

C     Execution.
C     ----------

C     Read the next line.  The ADVANCE = 'NO' must be present with SIZE = N.
C     NO:  Even with PAD = 'YES' on the OPEN, we get an end-of-record value
C          for IOS (-2 for DEC, -4006 for SGI) which cannot be portably
C          distinguished from EOF (IOS < 0).  It seems SIZE = ... is unusable.

C***  READ (LUN, '(A)', SIZE = N, IOSTAT = IOS, ADVANCE = 'NO') LINE
C***  IF (IOS /= 0) GO TO 99

C***  IF (N == 0) THEN
C***     GO TO 99
C***  ELSE
C***     N = MIN (N, LEN (LINE))
C***  END IF

      N = 0
      READ (LUN, '(A)', IOSTAT = IOS) LINE
      IF (IOS /= 0) GO TO 99

      N = LEN_TRIM (LINE)

C     Perform search for comments?
C     ----------------------------

      IF (COMMENT /= BLANK) THEN

C        Examine one character at a time in LINE (1 : N) for a
C        transition between significant characters and comments.

         I = 0
   10    CONTINUE
            I = I + 1
            IF (LINE (I : I) == COMMENT) THEN

C              Handle a special case here rather than impact all lines.
C              If the FIRST significant character is a COMMENT, consider
C              the line suppressed.  This covers the common case where
C              more than one COMMENT character is used to "comment-out"
C              an input line.
C
C              Search for preceding significant text:

               SUPPRESD = .TRUE.
               DO 20, J = 1, I-1
                  IF (LINE (J : J) /= BLANK .AND.
     >                LINE (J : J) /= HTAB) SUPPRESD = .FALSE.
   20          CONTINUE

               IF (SUPPRESD) THEN
                  N = 0
                  GO TO 40       ! For consistency; GO TO 99 is more direct.
               ELSE IF (I == N) THEN

C                 Single comment at the very end of the text - done.

                  N = I - 1
                  GO TO 40
               ELSE IF (LINE (I+1 : I+1) == COMMENT) THEN

C                 Double comment - strip the first one from the text,
C                 decrement N, and keep searching.  The DO-loop is
C                 clunky but required by the standard.

                  DO 30, J = I, N - 1
                     LINE(J : J) = LINE(J+1 : J+1)
   30             CONTINUE
                  N = N - 1
               ELSE

C                 Single comment embedded in the text - done.

                  N = I - 1
                  GO TO 40
               END IF
            END IF
            IF (I < N) GO TO 10

   40    CONTINUE
      END IF

C     Remove trailing blanks and tabs.
C     --------------------------------

C     We're taking advantage of the fact that if the DO-loop runs to
C     completion, then N = 0 on exit.

      DO 50, N = N, 1, -1
         IF (LINE (N : N) /= BLANK .AND.
     >       LINE (N : N) /= HTAB) GO TO 99
   50 CONTINUE

C     Termination.
C     ------------

   99 LAST = N
      RETURN
      END
