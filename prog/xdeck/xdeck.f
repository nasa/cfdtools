C+----------------------------------------------------------------------
C
      PROGRAM xdeck
C
C     XDECK is a program originally intended to remove the deck-names or
C     sequence numbers in columns 73-80 of files derived from card decks
C     or from the UPDATE facility on the CDC 7600,  and to eliminate any
C     trailing blanks.  (This halves the size of typical files.)   XDECK
C     has since been generalized to operate beyond any specified column.
C
C     The user is prompted for the file name and the column beyond which
C     truncation is desired.  (The default is column 72.)  A next higher
C     version of the file is produced,  in which all trailing blanks  as
C     well as the requested columns are suppressed.   The original  file
C     is left intact.   [Under Linux: no version number means the output
C     file name is the input name with '.xdeck' appended.]
C
C     Also,  XDECK  will  now  handle  the  removal  of  leading columns
C     (designed for use with files disposed from the CRAY via the  CYBER
C     front-end).
C
C     Logical units:      1  for input  file
C                         2  for output file
C                         5  for terminal input
C                         6  for terminal output
C
C     History: Original code on PDP-11:  David Saunders/Dan McCoy, 1/78.
C              Added leading column removal option:    Greg Howe, 11/82.
C              Patched "leading > length of short line" case: DAS 01/85.
C              Patched to handle long records gracefully:     DAS 08/85.
C              Added warning before first prompt:             RGL  1/86.
C              Under Linux, version number treatment was not  DAS 10/13.
C              appropriate.  Add .xdeck to the file name now.
C              The style encountered here is no longer recommended!
C              Carriage control in particular no longer uses '$'.
C
C     Author:  Dan McKernan, Informatics, 1/82: Generalized; FORTRAN 77.
C
C-----------------------------------------------------------------------

      IMPLICIT    NONE

      INTEGER     ioerr_open     , ioerr_close     ,
     >            read_file      , write_file      ,
     >            sys_input      , sys_output

      PARAMETER ( ioerr_open = 1 , ioerr_close = 2 ,
     >            read_file  = 1 , write_file  = 2 ,
     >            sys_input  = 5 , sys_output  = 6 )

      CHARACTER   errmsg(2)*14 , line*512 , file(2)*4 , message*60
      LOGICAL     eof
      INTEGER     ifile , ioerr , ipos , itrunc , leading , n , status

      DATA        itrunc    / 72  /              ,
     >            leading   / 0 /                ,
     >            eof       / .FALSE. /          ,
     >            errmsg(1) / ' open error on' / ,
     >            errmsg(2) / 'close error on' / ,
     >            file(1)   / ' old' /           ,
     >            file(2)   / ' new' /

C     -- Warn the user about stripping the leading blank and
C     -- truncating after 72 in the same step:

      WRITE ( UNIT = sys_output,
     >        FMT  = '(A)'     )
     >        'Beware:  Stripping the leading blank and truncating',
     >        'beyond column 72 leaves only 71 columns intact.',
     >        'In this case, truncate after 73 to yield a full 72.',
     >        ' '

C     -- Get input file name:

      n = 0

      DO WHILE ( n .LE. 0 .OR. n .GT. 80 )

         WRITE ( UNIT = sys_output ,
     >           FMT  = '(A)'      ,
     >           ADVANCE = 'no'    ) 'Enter input file name: '

         READ  ( UNIT = sys_input  ,
     >           FMT  = '(A)'      ,
     >           END  = 99         ) line
      n = len_trim(line)

      END DO

C

   5  CONTINUE

C     -- Get column # to truncate beyond:

      CALL accept_integer( 'Enter column to truncate beyond (default=72): ', itrunc )

      IF ( itrunc.EQ.-1 ) itrunc = 72

C     -- Get # of leading columns to remove:

      CALL accept_integer('Enter number of leading columns to remove (default=0): ', leading )

      IF ( leading.EQ.-1 ) leading = 0
      leading = leading+1

C

      IF ( leading.GE.itrunc ) THEN
         WRITE ( UNIT = sys_output ,
     >           FMT  = '(A)'      )
     >           ' ',
     >           '%%_ERROR, leading columns removed greater than or',
     >           '          equal to truncated columns -  re-enter.',
     >           ' '
         goto 5
      ENDIF

C     -- Input file name may or may not have version #. Either is OK:

      ifile = read_file

      ioerr = ioerr_open
      OPEN ( UNIT   = read_file ,
     >       FILE   = line(1:n) ,
     >       STATUS = 'old'     ,
     >       ERR    = 20        )

C     -- Output file is to be the next version number.
C     -- Ensure that any version number in the command is ignored:

      ipos = index( line , ';' )
      IF ( ipos. EQ. 0 ) then
         line(n+1:n+6) = '.xdeck'
         n = n + 6
      ELSE
         n = ipos
      ENDIF

C     -- Open output file:

      ifile = write_file
      OPEN ( UNIT   = write_file ,
     >       FILE   = line(1:n)  ,
     >       STATUS = 'new'      ,
     >       ERR    = 20         )

C     -- Process files line-by-line:

      DO WHILE ( .NOT. eof )

         READ  ( UNIT = read_file ,
     >           FMT  = '(A)'   ,
     >           IOSTAT = status  ) line

         IF ( status .EQ. 0 )

     >      THEN

C              -- Suppress trailing blanks:

               n = min ( len_trim(line) , itrunc )
               DO WHILE ( n .GT. 0  .AND. line(n:n) .EQ. ' ' )
                  n = n - 1
               END DO

C              -- Check for empty line. (Patch for n < leading could
C                 probably be done more efficiently. But it's safe.)

               IF ( n .LE. 0  .OR.  n .LT. leading )
     >            THEN
                     WRITE ( UNIT = write_file , FMT = '(A)' )
                  ELSE
                     WRITE ( UNIT = write_file , FMT = '(A)' )
     >                  line(leading:n)
               END IF

            ELSE IF ( status .LT. 0 )

     >         THEN

                  eof = .TRUE.

               ELSE

C                 -- This code may have been clobbered if n exceeded the
C                    buffer length - but try to die gracefully anyway:

                  WRITE ( UNIT = sys_output,
     >                    FMT  = '( '' Error reading input file -'',
     >                              '' record length may exceed 512.''/
     >                              '' Value found MAY be '', I8 )' ) n
                  eof = .TRUE.

            END IF

      END DO


      ioerr = ioerr_close

      ifile = read_file
      CLOSE ( UNIT = read_file  , ERR = 20 )

      ifile = write_file
      CLOSE ( UNIT = write_file , ERR = 20 )

      goto 99

   20 WRITE( UNIT = sys_output             ,
     >       FMT  = '(/,A,A,''file'')' ) errmsg(ioerr) , file(ifile)

   99 CALL exit

      END
C+----------------------------------------------------------------------
C
      SUBROUTINE accept_integer( message, integer_variable )
C
C     This routine reads in a character string and converts
C     it to an integer. The integer must be in the range of
C     0 <= integer < 10^10 .  If it is not then the user is
C     prompted again for a value.  If an end of file ( ^Z )
C     is read then the program is stopped.
C
C     Author: Daniel S. McKernan, Informatics, March 4 1982
C
C-----------------------------------------------------------------------

      INTEGER   successful_conversion , sys_input ,
     >          sys_output            , status

      CHARACTER buffer*80 , integer_format*4 , null*1
      CHARACTER(len=*) message

      LOGICAL   conversion_error

C     -- Define I/O return codes:

      PARAMETER ( successful_conversion = 0 )

C     -- Define I/O units:

      PARAMETER ( sys_input = 5 , sys_output = 6 )

      null = char( 0 )

C     -- Read character string until it is converted
C     -- to an integer or an end of file is read:

      conversion_error = .TRUE.
      DO WHILE ( conversion_error )

         WRITE ( UNIT = sys_output ,
     >           FMT  = '(A)'      )
     >           message

         READ ( UNIT = sys_input ,
     >          FMT  = '(A)'     ,
     >          END  = 20        )  buffer

C        -- Find the first nonblank and nonnull character:

         i = 1
         DO WHILE (( buffer(i:i) .EQ. ' ' ) .OR.
     >               buffer(i:i) .EQ. null )
            i = i + 1
         END DO

C        -- Find the last nonblank and nonnull character:

         j = 80
         DO WHILE (( buffer(j:j) .EQ. ' ' ) .OR.
     >               buffer(j:j) .EQ. null )
            j = j - 1
         END DO

         IF ( i .LE. j )
     >      THEN

C              -- Convert the characters to an integer:

               integer_format = '(I' // char( j - i + 1 + ichar( '0' ) )
     >                               // ')'
               READ ( UNIT   = buffer(i:j)    ,
     >                FMT    = integer_format ,
     >                IOSTAT = status         ) integer_variable

C              -- Check for an error in converting:

               IF ( status .EQ. successful_conversion )
     >              conversion_error = .FALSE.

            ELSE

C              -- User entered a carriage return use default:

               conversion_error = .FALSE.
               integer_variable = -1

            END IF

      END DO

      RETURN

C     -- A control Z was entered so stop the program:

   20 CALL exit
      END
