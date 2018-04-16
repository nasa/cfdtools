C+----------------------------------------------------------------------
C
      SUBROUTINE RDTUPLES (SCREEN, PROMPT, KEYBRD, DISK, NDIM, NTUPLES,
     >   X1, X2, X3)
C
C
C     Description and usage:
C
C        RDTUPLES reads an indefinite number of real pairs or triples (or
C     singles), one per line or record, with an optional prompt.  Standard
C     delimiters (comma, blank, or tab) are assumed between values.  An
C     indirect option is provided (entering @filename rather than entering
C     values at the keyboard).
C
C        RDTUPLES is appropriate for entering (X,Y) or (X,Y,Z) values into
C     distinct arrays, as opposed to entering multiple values into a single
C     array (for which RDREALS is appropriate).  The degenerate case (one
C     real X per record) is retained because it may be handy for picking off
C     the first column of a file.
C
C        For the indirect option, embedded "!" comments are ignored.  Blank
C     lines are also skipped since they're a poor indicator of end-of-data.
C     Reading continues to end-of-file.  For interactive input, interpretation
C     of carriage-return and control-Z is left to the application, with either
C     normally indicating end-of-data.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        SCREEN              I    I      Logical unit number to which the
C                                        prompt and diagnostics are written.
C        PROMPT              C    I      Prompt string with carriage control
C                                        in PROMPT (1:1).  (Use ' ', '0', '1',
C                                        or '$' in the usual way.)
C                                        If PROMPT is blank it is suppressed.
C        KEYBRD              I    I      Logical unit number from which data
C                                        is to be read interactively.
C        DISK                I    I      Logical unit number to use if the
C                                        indirect option (@filename) is invoked.
C                                        DISK = KEYBRD on input implies file
C                                        is already open (rare possibility).
C                                        Tuples are read until EOF; embedded "!"
C                                        comments or blank lines are permitted.
C        NDIM                I    I      NDIM = 1, 2, or 3 controls entry of
C                                        values in the arrays X1, X2, X3.
C                                        E.g. NDIM = 2 means X3 is ignored;
C                                        further, if only one value is read,
C                                        X2 is left undisturbed for that read.
C                                        This facilitates defaulting.
C        NTUPLES             I    I/O    Input with the maximum number of
C                                        tuples provided for.  On output,
C                                        NTUPLES > 0 is the no. of tuples read;
C                                        NTUPLES = 0 corresponds to null (CR);
C                                        NTUPLES =-1 corresponds to end-of-data;
C                                        NTUPLES =-2 means an invalid value was
C                                                    found (indirect mode only);
C                                        NTUPLES =-3 means a system-dependent
C                                                    I/O error was encountered.
C        X1      NTUPLES     R      O    X1, X2, X3 contain the tuples entered,
C        X2         "        "      "    if any.  X2 is ignored if NDIM < 2;
C        X3         "        "      "    X3 is ignored if NDIM < 3.
C
C
C     Implementation Notes:
C
C        Limiting NDIM to 3 is surely a reasonable compromise.  Full generality
C     here would require a matrix of outputs, defeating the utility of this
C     routine's multi-vector (or -scalar) handling over the single vector of
C     RDREALS.  The latter can still be used directly for more than 3 values
C     per input line.
C
C        Error handling is mostly performed by the lower level RDREALS.
C     Fewer than NDIM values on a record is not considered an error because
C     some flexibility would be lost:  defaulting would be impeded.  More than
C     NDIM values or NTUPLES tuples is also not considered an error:  RDREALS
C     and RDTUPLES terminate normally if they fill their arrays without finding
C     end of data.  Also, any indirect file is left open (for possible further
C     input).
C
C
C     Disclaimer:
C
C        We have here a case of diminishing returns:  the functionality
C     provided is nice to have, but it may not warrant yet another gold-plated
C     module.  On the other hand, RDREALS answers only part of the problem.
C     Where do we draw the line?  The hope remains that keeping this type of
C     code reusable and out of new applications will enhance productivity in
C     the long run.
C
C
C     Procedures:
C
C        OPENER   Opens a disk file (for the indirect option).
C        RDREALS  Decodes one record at a time (with optional prompt).
C
C
C     Environment:  DEC VAX/VMS FORTRAN, including:
C
C        >  IMPLICIT NONE
C        >  Trailing ! comments
C        >  Some 8-character variable names
C        >  '$' carriage control
C
C
C     Author:  David Saunders, Sterling Software, Palo Alto, CA.
C
C     02/28/89  DAS  (With Robert Kennelly:)  Initial implementation, making
C                    use of the newly-available RDREALS for the distinct case
C                    of tuples from one or more records, with the indirect
C                    option thrown in to add functionality.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      CHARACTER
     >   PROMPT * (*)
      INTEGER
     >   DISK, KEYBRD, NDIM, NTUPLES, SCREEN
      REAL
     >   X1 (*), X2 (*), X3 (*)

C     Internal COMMON (required by the RDTUPLES application of RDREALS). (Ugh!)

      INTEGER    MAXBUF
      PARAMETER (MAXBUF = 132)       ! Arbitrary limit on input record length
      CHARACTER  BUFFER * (MAXBUF)
      COMMON /TUPLES/ BUFFER

C     Local constants.

      CHARACTER
     >   BLANK * 1, INDIRECT * 1

      PARAMETER
     >  (BLANK    = ' ',
     >   INDIRECT = '@')

C     Local variables.

      INTEGER
     >   COUNT, LUN, ND
      REAL
     >   TUPLE (3)
      LOGICAL
     >   INTERACT

C     Procedures.

      EXTERNAL
     >   OPENER, RDREALS

C     Execution.


      LUN = KEYBRD
      INTERACT = LUN .NE. DISK   ! KEYBRD .EQ. DISK allows an already-open file.
      COUNT = 0

  200 CONTINUE

C        Look for a tuple.  Shorten the prompt after the first tuple if
C        mode is interactive; suppress it if a disk file is being read.

         ND = NDIM

         IF (COUNT .EQ. 0 .AND. INTERACT) THEN
            CALL RDREALS (SCREEN, PROMPT,    LUN, ND, TUPLE)
         ELSE IF (INTERACT) THEN
            CALL RDREALS (SCREEN, '$More? ', LUN, ND, TUPLE)
         ELSE
            CALL RDREALS (SCREEN, BLANK,     LUN, ND, TUPLE)
         END IF

         IF (BUFFER (1 : 1) .EQ. INDIRECT) THEN  ! Input has been redirected
            INTERACT = .FALSE.
            LUN = DISK

            CALL OPENER (SCREEN, BLANK, KEYBRD, BUFFER (2:), LUN, 'OLD')

         ELSE IF (ND .GT. 0) THEN                ! Valid tuple found.
            COUNT = COUNT + 1
            X1 (COUNT) = TUPLE (1)
            IF (ND .GT. 1 .AND. NDIM .GT. 1) X2 (COUNT) = TUPLE (2)
            IF (ND .GT. 2 .AND. NDIM .GT. 2) X3 (COUNT) = TUPLE (3)
         END IF

C        If room for more, and not end-of-data or line was blank...

         IF (COUNT .LT. NTUPLES .AND.
     >      (ND .GT. 0 .OR.
     >      (ND .EQ. 0 .AND. .NOT.INTERACT)))
     >GO TO 200                                  ! ... keep reading


C     Termination.

      IF (ND .GE. -1 .AND. COUNT .GT. 0) THEN
         NTUPLES = COUNT     ! ND = -1 or 0 means EOF or CR (OK if COUNT > 0)
      ELSE
         NTUPLES = ND        ! Distinguish EOF or CR if COUNT = 0; also, flag
      END IF                 ! a read error in indirect mode (PROMPT=BLANK).

      RETURN
      END
