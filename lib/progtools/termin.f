C+----------------------------------------------------------------------
C
      SUBROUTINE TERMIN (STRING, MARKER)
C
C
C     Description and usage:
C
C           Adds the specified marker to the end of a string.  It is
C        positioned either following the last significant character or
C        packed at the end of the string.  A blank string is a special
C        case:  to distinguish it from a string with the marker already
C        at the beginning, the new marker begins in position 2 with a
C        leading blank.  Thus a "null" string begins with the marker,
C        while a "blank" string begins with a blank.  This feature is
C        intended for handling defaults more flexibly.
C
C           TERMIN was written for QPLOT, to be used for adding the '$'
C        terminator required by DISSPLA for its self-counting option.
C        As generalized, it may find other uses as well.
C
C
C     Parameters:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING              H    I/O    Character string to be terminated.
C        MARKER              H    I      Terminal character string.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C        (2)  The marker may actually be a string - the positioning goes
C             through as for a single character, but STRING must be longer
C             than MARKER if the whole marker is to be fitted in !  No
C             check is made for this; instead, the marker is simply
C             truncated if it won't fit.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C        16 Jan. 1984    RAK    Initial design and coding.
C        24 Jan. 1984    RAK    Generalized for arbitrary markers.
C
C-----------------------------------------------------------------------


C     Variable declarations.
C     ----------------------

      IMPLICIT NONE

C     Parameters.

      CHARACTER
     >   BLANK, HT
      PARAMETER
     >   (BLANK = ' ', HT = CHAR (9))

C     Variables.

      INTEGER
     >   CURSOR, LENGTH, TAIL
      CHARACTER
     >   MARKER * (*), STRING * (*)


C     Executable statements.
C     ----------------------

      LENGTH = LEN (STRING)

C     Position the cursor at the last non-blank, or at the beginning of
C     the string whether blank or not.

      CURSOR = LENGTH + 1
   10 CONTINUE
         CURSOR = CURSOR - 1
         IF ((STRING (CURSOR : CURSOR) .EQ. BLANK .OR.
     >        STRING (CURSOR : CURSOR) .EQ. HT) .AND.
     >        CURSOR .GT. 1) GO TO 10

C     If the string isn't marked, do so.  Note that an entirely blank string
C     is left with a leading blank, e.g. ' $' ("blank") is distinguished
C     from '$' ("null").

      IF (STRING (CURSOR : LENGTH) .NE. MARKER) THEN
         TAIL = MIN (CURSOR + 1,
     >               MAX (2, LENGTH - LEN (MARKER) + 1))
         STRING (TAIL : LENGTH) = MARKER
      END IF


C     Termination.
C     ------------

      RETURN
      END
