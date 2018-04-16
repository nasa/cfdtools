C+----------------------------------------------------------------------
C
      SUBROUTINE LOOKUP (NDICT, DICTRY, ALPHA, KEY, ENTRY)
C
C
C     Description and usage:
C
C           Performs dictionary lookups.  A pointer is returned if a
C        match is found between the input key and the corresponding
C        initial characters of one of the elements of the dictionary.
C        If a "synonym" has been provided for an entry, the search is
C        continued until a match to a primary dictionary entry is found.
C        Cases of no match, or multiple matches, are also provided for.
C
C           Dictionary entries must be left-justified, and may be alphabetized
C        for faster searches.  Secondary entries, if any, are composed of
C        two words separated by one or more characters such as blank, tab,
C        comma, colon, or equal sign which are treated as non-significant
C        by SCANNR.  The first entry of each such pair serves as a synonym
C        for the second, more fundamental keyword.
C
C           The ordered search stops after the section of the dictionary
C        having the same first letters as the key has been checked, or
C        after a specified number of entries have been examined.  A special
C        dictionary entry, the vertical bar '|', will also terminate the
C        search.  This will speed things up if an appropriate dictionary
C        length parameter cannot be determined.  Both types of search are
C        sequential.  See "Notes" below for some suggestions if efficiency
C        is an issue.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        NDICT               I    I      Number of dictionary entries to be
C                                        examined.
C        DICTRY  NDICT       C    I      Array of dictionary entries,
C                                        left-justified in their fields.
C                                        May be alphabetized for efficiency,
C                                        in which case ALPHA should be .TRUE.
C                                        Entries with synonyms are of the form
C                                        'ENTRY:SYNONYM', where 'SYNONYM'
C                                        is a more fundamental entry in the
C                                        same dictionary.  NOTE: Don't build
C                                        "circular" dictionaries!
C        ALPHA               L    I      Indicates whether the dictionary
C                                        is in alphabetical order, in which
C                                        case the search can be terminated
C                                        sooner.
C        KEY                 C    I/O    String to be compared against the
C                                        dictionary.  Abbreviations are legal
C                                        provided they correspond to a unique
C                                        entry in the dictionary.  KEY is
C                                        replaced on termination by its most
C                                        fundamental equivalent dictionary
C                                        entry (uppercase, left-justified) if
C                                        a match was found.
C        ENTRY               I      O    Dictionary pointer.  If > 0, it
C                                        indicates which entry matched KEY.
C                                        In case of trouble, a negative value
C                                        means that a UNIQUE match was not
C                                        found - the absolute value of ENTRY
C                                        points to the second dictionary entry
C                                        which matched KEY.  Zero means that
C                                        NO match could be found.  ENTRY
C                                        always refers to the last search
C                                        performed - in searching a chain of
C                                        synonyms a non-positive value will be
C                                        returned if there is any break, even
C                                        if the original input key was found.
C
C
C     External references:
C
C        Name    Description
C        SCANNR  Finds first and last significant characters.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C        (2)  We have assumed that the dictionary is not too big.  If
C             many searches are to be done or if the dictionary has more
C             than a dozen or so entries, it may be advantageous to build
C             an index array of pointers to the beginning of the section
C             of the dictionary containing each letter, then pass in the
C             portion of the dictionary beginning with DICTRY (INDEX).
C             (This won't generally work for dictionaries with synonyms.)
C             For very large problems, a completely different approach may
C             be advisable, e.g. a binary search for ordered dictionaries.
C
C        (3)  LOOKUP is case sensitive.  In most applications it will be
C             necessary to use an uppercase dictionary, and to convert the
C             input key to uppercase before calling LOOKUP.  Companion
C             routines TOKENS and PAIRS, available from the author, already
C             take care of this.
C
C        (4)  The key need not be left-justified.  Any leading (or
C             trailing) characters which are "non-significant" to SCANNR
C             will be ignored.  These include blanks, horizontal tabs,
C             commas, colons, and equal signs.  See SCANNR for details.
C
C        (5)  The ASCII collating sequence for character data is assumed.
C             (Note that this means the numerals precede the alphabet, unlike
C             common practice!)  On some machines, it may be necessary to
C             use the FORTRAN lexical library routines to force use of the
C             ASCII sequence.
C
C        (6)  Parameter NUMSIG sets a limit on the length of significant
C             dictionary entries.  Special applications may require that this
C             be increased.  (It is 16 in the present version.)
C
C        (7)  No protection against "circular" dictionaries is provided: don't
C             claim that A is B, and that B is A.  All synonym chains must
C             terminate!  Other potential errors not checked for include
C             duplicate or mis-ordered entries.
C
C        (8)  The handling of ambiguities introduces some ambiguity:
C
C                ALPHA = .TRUE.  A potential problem, when one entry
C                                looks like an abbreviation for another
C                                (eg. does 'A' match 'A' or 'AB'?) was
C                                resolved by dropping out of the search
C                                immediately when an "exact" match is found.
C
C                ALPHA = .FALSE. The programmer must ensure that the above
C                                situation does not arise: each dictionary
C                                entry must be recognizable, at least when
C                                specified to full length.  Otherwise, the
C                                result of a search will depend on the
C                                order of entries.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C        24 Feb. 1984  RAK/DAS  Initial design and coding.
C        25 Feb. 1984    RAK    Combined the two searches by suitable
C                               choice of terminator FLAG.
C        28 Feb. 1984    RAK    Optional synonyms in dictionary, no
C                               longer update KEY.
C        29 Mar. 1984    RAK    Put back replacement of KEY by its
C                               corresponding entry.
C        21 June 1984    RAK    Corrected bug in error handling for cases
C                               where no match was found.
C        23 Apr. 1985    RAK    Introduced test for exact matches, which
C                               permits use of dictionary entries which
C                               would appear to be ambiguous (for ordered
C                               case).  Return -I to point to the entry
C                               which appeared ambiguous (had been -1).
C                               Repaired loop termination - had to use
C                               equal length strings or risk quitting too
C                               soon when one entry is an abbreviation
C                               for another.  Eliminated HIT, reduced
C                               NUMSIG to 16.
C        16 May  1988    DAS    Had to use LEN in definition of LAST to
C                               suit revised SCANNR; early termination if
C                               KEY length exceeds LEN (DICTRY (1)).
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   ENTRY, NDICT
      LOGICAL
     >   ALPHA
      CHARACTER
     >   DICTRY (NDICT) * (*), KEY * (*)

C     Local constants.

      INTEGER
     >   NUMSIG
      CHARACTER
     >   BLANK, CURLY
      PARAMETER
     >   (BLANK = ' ', CURLY = '{', NUMSIG = 16)

C     Local variables.

      INTEGER
     >   FIRST, I, LAST, LENDIC, LENGTH, MARK
      CHARACTER
     >   FLAG * (NUMSIG), TARGET * (NUMSIG)

C     Procedures.

      EXTERNAL
     >   SCANNR


C     Execution.
C     ----------

      ENTRY = 0

C     Isolate the significant portion of the input key (if any).

      FIRST = 1
      LAST = LEN (KEY)
      CALL SCANNR (KEY, FIRST, LAST, MARK)
      IF (MARK .EQ. 0) GO TO 99

C     Can't hope to find a match if the key is longer than dictionary entries.

      LENGTH = MARK - FIRST + 1
      LENDIC = LEN (DICTRY (1))
      IF (LENGTH .GT. LENDIC) GO TO 99


C     The search starts with the input key, but may be repeated if that
C     target is just a synonym for a more fundamental dictionary entry.
C     NUMSIG = LEN (TARGET) is assumed to be plenty big enough.

      TARGET = KEY (FIRST:MARK)

   10 CONTINUE

C        Select search strategy by cunning choice of termination test
C        flag.  The left curly bracket follows all the alphabetic
C        characters in the ASCII collating sequence, but precedes the
C        vertical bar.

         IF (ALPHA) THEN
            FLAG = TARGET
         ELSE
            FLAG = CURLY
         END IF


C        Perform search.
C        ---------------

         I = 0
   20    CONTINUE
            I = I + 1
            IF (TARGET (1:LENGTH) .EQ. DICTRY (I) (1:LENGTH)) THEN
               IF (ENTRY .EQ. 0) THEN

C                 First "hit" - must still guard against ambiguities
C                 by searching until we've gone beyond the key (ordered
C                 dictionary), or until the end-of-dictionary mark is
C                 reached (exhaustive search).

                  ENTRY = I

C                 Special handling if match is exact - terminate search.
C                 We thus avoid confusion if one dictionary entry looks
C                 like an abbreviation of another.  This fix won't
C                 generally work for un-ordered dictionaries!

                  FIRST = 1
                  LAST = LENDIC
                  CALL SCANNR (DICTRY (ENTRY), FIRST, LAST, MARK)
                  IF (MARK .EQ. LENGTH) I = NDICT
               ELSE


C                 Oops - two hits!  Abnormal termination.
C                 ---------------------------------------

                  ENTRY = -I
                  GO TO 99
               END IF
            END IF

C           Check whether we've gone past the appropriate section of the 
C           dictionary.  The test on the index provides insurance and an
C           optional means for limiting the extent of the search.

            IF (DICTRY (I) (1:LENGTH) .LE. FLAG .AND. I .LT. NDICT)
     >   GO TO 20


C        Check for a synonym.
C        --------------------

         IF (ENTRY .GT. 0) THEN

C           Look for a second entry "behind" the first entry.  (FIRST
C           and MARK were determined above when the hit was detected.)

            FIRST = MARK + 2
            CALL SCANNR (DICTRY (ENTRY), FIRST, LAST, MARK)
            IF (MARK .GT. 0) THEN

C              Reset the target and dictionary pointer and repeat the
C              search for the synonym instead of the original key.

               TARGET = DICTRY (ENTRY) (FIRST:MARK)
               LENGTH = MARK - FIRST + 1
               ENTRY = 0
               GO TO 10

            ELSE

C              Expand the key to the full dictionary entry as a possible aid
C              to the calling program (which may prefer to avoid dealing with
C              entry numbers).

               KEY = DICTRY (ENTRY)

            END IF
         END IF


C     Normal termination.
C     -------------------

   99 RETURN
      END
