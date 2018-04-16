C+----------------------------------------------------------------------
C
      SUBROUTINE TOGGLE (NITEMS, MENU, ONOFF, LUNCRT, LUNKBD)
C
C
C     One-liner:  Utility for displaying and switching on/off options
C
C     Purpose:
C
C           TOGGLE modularizes the situation of turning on or off any of
C        a number of options.  It displays the current settings, prompts
C        for switching a setting, redisplays, and so on until a simple
C        carriage return is detected, meaning "proceed."  Encountering
C        EOF (CTRL Z) is interpreted as meaning "stop here."
C
C           Dealing strictly with yes/no options makes for an easy
C        implementation.  Whether TOGGLE provides functionality of broad
C        enough applicability to be warranted remains to be seen.
C
C     Sample Usage:
C
C        DATA
C           MENU
C             /'1: Plot Cl vs. span station?',
C              '2: Plot Cm vs. span station?',
C              '3: Plot Cd vs. span station?',
C                :   :   :   :   :   :   :  /,
C           ONOFF
C             /MXMENU * .TRUE./
C        
C                :   :   :   :
C
C        CALL TOGGLE (MXMENU, MENU, ONOFF, LUNCRT, LUNKBD)
C
C                :   :   :   :
C
C        IF (ONOFF (1)) THEN
C                :   :   :   :
C        END IF
C
C                :   :   :   :
C
C
C           Note that inclusion of the item number in the menu text (as
C        opposed to its insertion by TOGGLE) makes for better clarity
C        in the application program when ONOFF (*) is referenced.  Also,
C        inclusion of the ?s saves worrying where the last nonblank is here.
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C
C        NITEMS      -        I     I    Number of items in the menu
C        MENU     NITEMS      C     I    Menu, with each line in the form
C                                        shown by the above example, including
C                                        the item numbers 1 through NITEMS.
C        ONOFF    NITEMS      L    I/O   Option switches, input with default
C                                        values; output with settings possibly
C                                        switched
C        LUNCRT      -        I     I    Logical unit number for the screen.
C        LUNKBD      -        I     I    Logical unit number for the keyboard.
C
C     External references:
C
C        Name    Description
C        READER  Prompting utility
C
C     Environment:
C
C        VAX/VMS FORTRAN 77
C
C     Known System Dependencies:
C
C        1.  IMPLICIT NONE is non-standard.
C        3.  Use of escape sequences for inverse video is a frill that could
C            be eliminated easily if necessary.  (What happens on hard copy
C            terminals?)
C
C     Author:
C
C        David Saunders, Sterling Software, Palo Alto, CA.
C
C     Development History:
C
C        06/25/87  DAS  Initial implementation (adapted from SELECT)
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   NITEMS, LUNCRT, LUNKBD
      LOGICAL
     >   ONOFF (NITEMS)
      CHARACTER
     >   MENU (NITEMS) * (*)

C     Local Constants:

      CHARACTER
     >   ESC, INVERSE * 4, RESET * 4

      PARAMETER
     >  (ESC     = CHAR (27),
     >   INVERSE = ESC // '[7m',
     >   RESET   = ESC // '[0m')

C     Local Variables:

      INTEGER
     >   ITEM
      LOGICAL
     >   PROCEED, QUIT
      CHARACTER
     >   YESNO * 3

C     Procedures:

      EXTERNAL
     >   READI

C     Execution.

 200  CONTINUE

C        Choose the simple approach: display the menu every time.

         WRITE (LUNCRT, 1001)
         DO 300, ITEM = 1, NITEMS
            IF (ONOFF (ITEM)) THEN
               YESNO = 'Yes'
            ELSE
               YESNO = 'No '
            END IF
            WRITE (LUNCRT, 1001) MENU (ITEM), INVERSE, YESNO, RESET
  300    CONTINUE

         WRITE (LUNCRT, 1001)
         CALL READI (LUNCRT,
     >      'Enter item number to switch, or <CR> to proceed: ',
     >      LUNKBD, ITEM, PROCEED, QUIT)

         IF (QUIT) STOP 'Stopping as requested.'
         IF (PROCEED) GO TO 999
         IF (ITEM .GE. 1 .AND. ITEM .LE. NITEMS) THEN
            ONOFF (ITEM) = .NOT. ONOFF (ITEM)
         ELSE
C           Bad entry - just redisplay and prompt again.
         END IF

      GO TO 200

C     Normal termination.

  999 RETURN

C     Format statements.

 1001 FORMAT (4X, A, 1X, A4, A3, A4)

      END
