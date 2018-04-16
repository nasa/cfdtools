C+----------------------------------------------------------------------
C
      SUBROUTINE LOOKUP2 (NDICT, DICTRY, KEY, ENTRY)
C
C     Description:
C
C           This is a simplified version of LOOKUP for the case when exact
C        matches are essential and the abbreviation and synonym ideas are
C        inappropriate.  The alphabetic option is also abandoned.
C
C           LOOKUP2 performs dictionary lookups.  A pointer is returned if a
C        match is found between the input key and the corresponding initial
C        characters of one of the elements of the dictionary.
C
C           Dictionary entries must be left-justified.
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        NDICT               I    I      Number of dictionary entries to be
C                                        examined.
C        DICTRY  NDICT       C    I      Array of dictionary entries,
C                                        left-justified in their fields.
C        KEY                 C    I/O    String to be compared against the
C                                        dictionary.
C        ENTRY               I      O    Dictionary pointer.  If ENTRY > 0, it
C                                        which matched KEY exactly (except for
C                                        possible trailing blanks in the
C                                        dictionary elements.  Zero means that
C                                        no match could be found.
C
C     12/16/05  DAS  Initial hack of LOOKUP for a real gas species application.
C
C     Author:  David Saunders, ELORET/NASA Ames Research Center.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   ENTRY, NDICT
      CHARACTER
     >   DICTRY (NDICT) * (*), KEY * (*)

      INTEGER
     >   I, LENDIC, LENKEY

C     Execution.

      ENTRY  = 0
      LENDIC = LEN (DICTRY (1))
      LENKEY = LEN_TRIM (KEY)

      IF (LENDIC >= LENKEY) THEN  ! Else an exact match is impossible

         DO I = 1, NDICT

            LENDIC = LEN_TRIM (DICTRY (I))

            IF (LENDIC == LENKEY) THEN
               IF (KEY (1:LENKEY) == DICTRY (I) (1:LENKEY)) THEN
                  ENTRY = I
                  EXIT
               END IF
            END IF

         END DO

      END IF

      END SUBROUTINE LOOKUP2
