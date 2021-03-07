C+----------------------------------------------------------------------
C
      SUBROUTINE GETDIS (LUNCRT, LUNKBD, MODE, NP, P, NPTS)
C
C  PURPOSE:  GETDIS does the prompting needed for most uses of DSTRIB,
C            and possibly other 1-D grid utilities.
C
C  ARGUMENTS:
C    ARG    TYPE  I/O/S   DIM     DESCRIPTION
C    LUNCRT   I     I      -      Logical unit for screen prompts.
C    LUNKBD   I     I      -        "   "   "   "  keyboard responses.
C    MODE     I     O      -      MODE with which DSTRIB is to be used
C                                 (-1, 0, 1, 2, or 3 - see DSTRIB - or
C                                 4 for the Vinokur distribution, or
C                                 5 for the sinusoid + quadratic). If
C                                 MODE = -99, calling program can quit.
C    NP       I    I/O     -      Number of distribution params. needed.
C                                 Input with maximum room provided (at
C                                 least 1);  output with actual length
C                                 of P (*) if MODE > 0, else not set.
C    P        R     O      NP     Parameters prompted for.  See NP.
C    NPTS     I    I/O      -     Input with max. # pts. allowed by the
C                                 application; output with # points in
C                                 desired distribution.
C
C  PROCEDURES:
C    READER      Prompting utility
C
C  ENVIRONMENT:  FORTRAN 77
C
C  HISTORY:
C  05/26/85    DAS    Original design and code (DSTRIB options only).
C  06/17/85    DAS    Switched to MODE=-99, not NPTS<0, to mean "quit".
C  05/29/90    DAS    Clarified internal point prompt (X, not I).
C  05/04/92    DAS    Provided for the Vinokur distribution.
C  04/01/95    DAS    Provided for FOILGRID.
C  10/18/96    DAS    Replaced FOILGRID with FOILGRD.
C  07/15/97    DAS    Slight adaptation of the PROFILE utility for
C                     use by BSPROFILE (64-bit arithmetic).
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, MODE, NP, NPTS
      REAL
     >   P (NP)

C     Local variables:

      INTEGER
     >   MAXPTS
      CHARACTER
     >   YESNO * 1
      LOGICAL
     >   DEFAULT, QUIT

C     Procedures:

      EXTERNAL
     >   READI, READD

C     Execution:

      WRITE (LUNCRT, 1001) ' ',
     >  'Discretization options:',
     >  '   0:  Uniform distribution',
     >  '   1:  Sinusoidal bunching towards the lower end',
     >  '   2:  Symmetric sinusoidal bunching towards both ends',
     >  '   3:  Sinusoidal bunching around an internal point',
     >  '   4:  Vinokur distribution (1st & last increment specified)',
     >  '   5:  Linear + Quadratic + Sine + Cosine combination',
     >  ' '

  100 MODE = 5
      CALL READI (LUNCRT,
     >   'Enter grid choice. <CR> = 5; ^Z (or ^D) = quit: ',
     >   LUNKBD, MODE, DEFAULT, QUIT)

      IF (MODE .LT. 0 .OR. MODE .GT. 5) GO TO 100

      IF (QUIT) THEN
         MODE = -99
         GO TO 999
      END IF

      MAXPTS = NPTS
  200 NPTS = 97
      CALL READI (LUNCRT, 'How many points per surface?  <CR> = 97: ',
     >            LUNKBD, NPTS, DEFAULT, QUIT)

      IF (NPTS .GT. MAXPTS .OR. MAXPTS .LT. 2) THEN
         WRITE (LUNCRT, '(A, 2I5)') ' Present limits: ', 2, MAXPTS
         GO TO 200
      END IF

      IF (MODE .EQ. 0) THEN

C  *     Uniform distribution - no other parameters needed:

      ELSE IF (MODE .GE. 1 .AND. MODE .LE. 3) THEN

C  *     Sinusoidal-type distributions require a "WIDTH"-type exponent:

         WRITE (LUNCRT, 1001) ' ',
     >      'Exponents > 1.0 give bunching in further regions;',
     >      'fractional exponents in the range (0.,1.] do not.',
     >      ' '

         P (1) = 1.0
         CALL READD (LUNCRT, 'Enter exponent. (Default is 1.) ',
     >               LUNKBD, P (1), DEFAULT, QUIT)

         IF (MODE .EQ. 3) THEN

            NP = 2

C           Default here is misleading - data not necessarily in [0.,1.]

            P (2) = 0.5
            CALL READD (LUNCRT,
     >     'Internal pt. (X or T) about which pts. are to be bunched: ',
     >         LUNKBD, P (2), DEFAULT, QUIT)

         ELSE

            NP = 1

         END IF

      ELSE IF (MODE .EQ. 4) THEN   ! Vinokur distribution

         NP = 2
  400    P (1) = -0.2
         WRITE (LUNCRT, 1001)
     >      'First increment?  +ve is absolute; -ve is relative;'
         CALL READD (LUNCRT,
     >      '-r means r% of range; default = 0.2% of range: ',
     >      LUNKBD, P (1), DEFAULT, QUIT)
         IF (P (1) .EQ. 0.) GO TO 400

         P (2) = P (1) + P (1)
         CALL READD (LUNCRT,
     >      'Last increment?  Default = twice the first: ',
     >      LUNKBD, P (2), DEFAULT, QUIT)
         IF (P (2) .EQ. 0.) GO TO 400

      ELSE IF (MODE .EQ. 5) THEN   ! Linear + quadratic + sine + cosine

         NP = 4
         P (1) = 0.04
         P (2) = 0.0
         P (3) = 0.3
         P (4) = 0.66

         WRITE (LUNCRT, 1001)
     >      'Defaults for L, Q, S, C terms are 0.04, 0.0, 0.3, 0.66.'
         YESNO = 'Y'
         CALL READC (LUNCRT, 'Take the defaults? (Y/N; <CR>=Y): ',
     >      LUNKBD, YESNO, DEFAULT, QUIT)

         IF (YESNO .NE. 'Y') THEN
            CALL READD (LUNCRT, 'Weight on    LINEAR term?  [0.04]: ',
     >         LUNKBD, P (1), DEFAULT, QUIT)
            CALL READD (LUNCRT, 'Weight on QUADRATIC term?  [0.00]: ',
     >         LUNKBD, P (2), DEFAULT, QUIT)
            CALL READD (LUNCRT, 'Weight on      SINE term?  [0.30]: ',
     >         LUNKBD, P (3), DEFAULT, QUIT)
            CALL READD (LUNCRT, 'Weight on    COSINE term?  [0.66]: ',
     >         LUNKBD, P (4), DEFAULT, QUIT)
         END IF

      END IF

  999 RETURN

C     Formats:

 1001 FORMAT (1X, A)

      END
