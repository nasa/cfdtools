C+----------------------------------------------------------------------
C
      SUBROUTINE PARSENAM (STRING, NAMFIRST, NAMLAST, EXTFIRST, EXTLAST)
C
C     One-liner:
C
C        Isolate file name and extension if any (VMS or Unix).
C
C     Description and usage:
C
C           PARSNAME identifies the basic file name and extension (if any)
C        within a string containing a file specification.  Since it uses a
C        backward search, it may apply to VAX/VMS file names or to Unix
C        file names, and is thus portable.  STRING (NAMFIRST : NAMLAST) is
C        determined as the basic file name; STRING (EXTFIRST : EXTLAST) is
C        determined as the file extension if one is present.  For example:
C
C           dev:[dir.subdir]some_file.type;2    (VMS) and
C           /d1/dx/whatever/some.file.type      (Unix)
C                           ^       ^ ^  ^
C        would produce the same pointers as indicated.  Note that the LAST
C        "extension" is isolated under Unix, while ;n version numbers are
C        not valid for Unix.
C
C           In order for STRING (NAMFIRST : EXTLAST) always to delimit the
C        file name apart from directory and version number even if the
C        basic name or an extension is missing, it was necessary to use
C        NAMLAST = 0 and EXTFIRST = 0 as the signals for missing elements
C        - user beware.  (If BOTH are zero, at most a directory was present.)
C
C           Any device/directory or version number can thus be deduced if
C        desired, but the initial intent here is to modularize extraction
C        of the basic name for the common case of deriving a related file
C        name.  Delimination of a possible file extension, and of the full
C        name other than directory and version number, were included in
C        case they prove handy.
C
C     Arguments:
C
C        Name      Type  I/O/S  Description
C        STRING     C    I      STRING (1 : LEN (STRING)) will be scanned
C                               for a file name and extension.
C
C        NAMFIRST,  I      O    First & last characters of basic file name.
C          NAMLAST              NAMLAST = 0 signals a missing basic name as
C                               for /dir/.alias (say).  NAMFIRST = 6 here.
C
C        EXTFIRST,  I      O    First & last characters of any extension.
C          EXTLAST              EXTFIRST = 0 signals no extension, as for
C                               /dir/mydata or mydata. (say).  EXTLAST is
C                               11 or 12 in these cases respectively.
C
C     Procedures:
C
C        SCAN3      Backward-scan utility
C
C
C     Environment:  VAX/VMS or Unix, FORTRAN 77, with
C                   IMPLICIT NONE, 8-character names, and ! comments
C
C     History:
C
C        06 Jan. 1991   DAS   Initial implementation.
C
C     Author:  David Saunders, NASA Ames/Sterling Software, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      CHARACTER
     >   STRING * (*)
      INTEGER
     >   NAMFIRST, NAMLAST, EXTFIRST, EXTLAST

C     Local variables.

      INTEGER
     >   FIRST, LAST, MARK

C     Procedures.

      EXTERNAL
     >   SCAN3

C     Execution.

C     First, dismiss any device/directory spec.
C     Allow for trailing blanks, but not undefined characters.

      LAST = LEN (STRING)

      CALL SCAN3 (STRING, ' ]/', 1, LAST, MARK)  ! Backward search

      NAMFIRST = 1
      NAMLAST = 0
      EXTFIRST = 0
      EXTLAST = LAST

      IF (MARK .EQ. 0) GO TO 99          ! Remove the most degenerate case

C     Now dismiss any version number.

      FIRST = MARK                       ! Start of significant text
      NAMFIRST = FIRST

      CALL SCAN3 (STRING, ';', FIRST, LAST, MARK)

      IF (MARK .EQ. 0) GO TO 99          ! ; is a grubby possibility
      IF (MARK .EQ. FIRST + 1) GO TO 99  ! ;n is too

      IF (MARK .GT. FIRST + 1) THEN      ! More than just a version
         LAST = MARK - 2
         EXTLAST = LAST
      END IF                             ! Else MARK = FIRST (no version #)

C     Directory and version have been cleared out of the way.
C     Now distinguish basic name from extension.  Either may be missing.

      IF (STRING (LAST : LAST) .EQ. '.') THEN  ! No extension proper
         IF (LAST .GT. FIRST) THEN
            NAMLAST = LAST - 1
         END IF                             ! And we're done either way

      ELSE                               ! Search backwards for '.'

         CALL SCAN3 (STRING, '.', FIRST, LAST, MARK)

         IF (MARK .GT. FIRST) THEN          ! Valid extension
            EXTFIRST = MARK                 ! Normal case
            IF (STRING (FIRST : FIRST) .NE. '.') THEN
               NAMLAST = MARK - 2              ! Normal case
            END IF                             ! Else name is missing
         ELSE                               ! Name only
            NAMLAST = LAST                  ! Since extension is missing
         END IF

      END IF
        
   99 RETURN
      END
