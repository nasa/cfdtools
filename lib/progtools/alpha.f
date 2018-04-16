C+----------------------------------------------------------------------
C
      FUNCTION ALPHA (STRING)
C
C
C     Description and usage:
C
C           A simple(-minded) test for numeric data is implemented by
C        searching an input string for disqualifying characters.  Most,
C        but not all, disqualifiers are included in the definition of
C        statement function LETTER - included are some punctuation marks
C        and the alphabet, except for the letters D and E which may be
C        part of a single or double precision quantity.  Some insurance
C        is provided by requiring that any numeric string have at least
C        one digit.  Note that a few ambiguities remain:
C
C           (a)  A string might have the form of numerical data but be
C                intended as text - no general test can hope to detect
C                such cases.
C        
C           (b)  This routine does not check for correctness of the
C                data format.  A meaningless string such as 'E1E2E3'
C                will not be tagged as ALPHA, but is not numeric either.
C
C        Despite these weaknesses, this method should work in the vast
C        majority of ordinary cases.
C
C           This routine was written for use by QPLOT but may find other
C        applications where loosely-formatted input data is required.
C
C
C     Parameters:
C
C        Name    Dimension  Type  I/O/S  Description
C        ALPHA               L      O    Set .TRUE. if STRING is not
C                                        numerical data.
C        STRING              C    I      Input data to be tested.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C        (2)  The test and conversion of lowercase letters is meaningless
C             for systems which do not recognize lowercase - it may be
C             deleted.
C
C        (3)  The ASCII character set and collating sequence is assumed in
C             the choice of some symbols used and in the use of comparison
C             functions .LE. and .GE. (rather than the more general lexical
C             comparison library functions).
C
C        (4)  For QPLOT, any title line may be guaranteed to be alphabetic
C             by enclosing it in matching quotes (single or double).  The
C             quotes are stripped out by STRIPPER.
C
C        (5)  COMPLEX data with parentheses will look alphabetic.
C
C
C     Development history:
C
C        23 Feb. 1984      RAK     Initial design and coding.
C        27 Feb. 1984      RAK     Added more punctuation, including '$',
C                                  to the list of disqualifiers.
C        14 June 1984      RAK     Amended header.
C        16 Nov. 1991  D.Saunders  Lone dates like 11/16/91 are now deemed
C                                  alphabetic.  (Added '/' to LETTER;
C                                  also added ATOM .GE. '{'.)
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      LOGICAL
     >   ALPHA
      CHARACTER
     >   STRING * (*)

C     Local variables.

      LOGICAL
     >   DIGIT, LETTER, MAYBE
      INTEGER
     >   I, LENGTH
      CHARACTER
     >   ATOM * 1

C     Statement functions.

      LETTER (ATOM) = ATOM .GE. ':' .AND. ATOM .LE. '_' .AND.
     >                ATOM .NE. 'D' .AND. ATOM .NE. 'E' .OR.
     >                ATOM .GE. '!' .AND. ATOM .LE. '*' .OR.
     >                ATOM .GE. '{'  .OR. ATOM .EQ. '/'

      DIGIT (ATOM) =  ATOM .GE. '0' .AND. ATOM .LE. '9'


C     Execution.
C     ----------

      LENGTH = LEN (STRING)
      MAYBE = .FALSE.

      I = 0
   10 CONTINUE
         I = I + 1
         ATOM = STRING (I : I)

C        Convert any lowercase letters to uppercase.

         IF (ATOM .GE. 'a' .AND. ATOM .LE. 'z')
     >      ATOM = CHAR (ICHAR (ATOM) - 32)

C        The presence of any "letter" means the string is NOT numeric.  For
C        insurance, keep track of whether any honest digits have been found.

         ALPHA = LETTER (ATOM)
         MAYBE = MAYBE .OR. DIGIT (ATOM)
         IF (.NOT.ALPHA .AND. I .LT. LENGTH) GO TO 10

      ALPHA = ALPHA .OR. .NOT.MAYBE


C     Termination.
C     ------------

      RETURN
      END
