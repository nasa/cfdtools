!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program a2b  ! ASCII to Binary

!  Convert an ASCII file (probably large) to its unformatted equivalent.
!  The all-text case is much simpler, so it is prompted for.  If reals,
!  integers, and logicals have to be converted, each line has to be tokenized,
!  and each token has to be tested for type.  Writing each token as it is
!  decoded avoids the problem of assembling an equivalent list of correct
!  data types to be written as one record.  Default types are assumed, so
!  we don't attempt to distinguish 2-byte or 8-byte integers, etc.
!  The output file name has the string '.unformatted' (or possibly '.u')
!  appended to the input file name.  Reading this file is bound to be
!  application-specific.
!
!  History:
!
!     12/20/2017  D.A.Saunders  Initial implementation, with a FIAT material
!                               database (7,700+ lines) in mind.
!     12/22/2017    "     "     Pondered the hard case.  Jeff Hill suggested
!                               writing tokens as we go, using advance='no',
!                               but that keyword is illegal for both unformatted
!                               and stream file writes.
!     12/28/2017    "     "     Tried stream file output with the pos=n keyword.
!     12/29/2017    "     "     Stream files work, but alphnumeric tokens need
!                               to be preceded by their length in the output.
!     01/03/2018    "     "     Hadn't meant to overlook logical variables.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      lunin   = 1,       &  ! Input ASCII file
      lunout  = 2,       &  ! Output file, unformatted
      lunkbd  = 5,       &  ! Keyboard inputs
      luncrt  = 6,       &  ! Screen prompts and diagnostics
      maxbuf  = 256,     &  ! Limit on the length of any input line
      maxname = 132         ! Limit on file name lengths

!  Variables:

   integer :: &
      ios, l, length, n, nibytes, nrbytes

   logical :: &
      cr, eof, pure_text, r4

   real :: &
      eps

   character (3) :: &
      seps
   character (maxbuf) :: &
      buffer
   character (maxname) :: &
      filein, fileout

!  Execution:

   seps = ' ' // ',' // char (9)  ! Include tab; do it once at the higher level
   r4   = epsilon (eps) > 1.e-10  ! Else 8-byte reals
   if (r4) then
      nrbytes = 4
   else
      nrbytes = 8
   end if
   nibytes = 4  ! Default integers only

   filein = 'ASCII.txt'
   call reads (luncrt, 'Input ASCII file: ', &
               lunkbd, filein, cr, eof)
   if (eof) go to 99
   l = len_trim (filein)
   open (lunin, file=filein(1:l), status='old', iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(2a)') 'Trouble opening file ', filein(1:l)
      go to 99
   end if

   if (l < maxname-1) then
      fileout = filein(1:l) // '.unformatted'
   else
      write (luncrt, '(a, i4)') 'File name limit encountered:', l
      go to 99
   end if

   pure_text = .false.
   call ready (luncrt, 'Treat the file as plain text? [y|n; <cr>=no]: ', &
               lunkbd, pure_text, cr, eof)
   if (eof) go to 99

   if (pure_text) then
      open (lunout, file=fileout, status='unknown', form='unformatted')
   else
      open (lunout, file=fileout, status='unknown', form='unformatted', &
            access='stream')
   end if

   l = 0
   n = 1  ! Stream file position

   do  ! Until EOF
      l = l + 1
      read (lunin, '(a)', iostat=ios) buffer
      if (ios < 0) exit  ! EOF

      length = len_trim (buffer)

      if (length == maxbuf) then
         write (luncrt, '(a, i12, i6)') &
            'Trouble: Buffer length encountered. Line #, limit:', l, maxbuf, &
            'Aborting.'
         go to 99
      end if

      if (pure_text) then
         write (lunout) buffer(1:length)  ! Plain character data
      else
         call treat_data_types ()  ! Treats integers, reals, logicals too
      end if

   end do

   write (luncrt, '(a, i12)') 'Number of lines processed:', l

99 close (lunin,  iostat=ios)
   close (lunout, iostat=ios)

!  Internal procedure for program a2b:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine treat_data_types ()
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: first, inumber, ios, last, mark, nchar
      logical :: TorF
      real    :: rnumber
      character (7) :: TFbuffer
      logical, external :: alpha  ! Finds if a token is unlikely to be a number

      first = 1
      last  = length
      do  ! Until no more tokens

         call scan2 (buffer, seps, first, last, mark)  ! Isolate next token
         if (mark == 0) exit

         if (alpha (buffer(first:mark))) then  ! Not a number
            nchar = mark - first + 1
            if (nchar == 1 .or. nchar == 6 .or. nchar == 7) then  ! Logical?
               TFbuffer(1:nchar) = buffer(first:mark)  ! Avoid changing text
               call upcase (TFbuffer(1:nchar))
               select case (nchar)
                  case (1)
                     TorF = TFbuffer(1:1) == 'T' .or. &
                            TFbuffer(1:1) == 'F'
                  case (6)
                     TorF = TFbuffer(1:nchar) == '.TRUE.'
                  case (7)
                     TorF = TFbuffer(1:nchar) == '.FALSE.'
               end select
               if (TorF) then
                   TorF = index (TFbuffer(1:nchar), 'T') /= 0
                   write (lunout, pos=n) TorF
                   n = n + nibytes
               else  ! As for the not-a-logical case
                  write (lunout, pos=n) nchar, buffer(first:mark)
                  n = n + nibytes + nchar
               end if
            else  ! Not a logical value; treat as text
               write (lunout, pos=n) nchar, buffer(first:mark)
               n = n + nibytes + nchar
            end if
         else
            call decode_number (buffer(first:mark), inumber, rnumber, ios)

            select case (ios)
               case (0)  ! Neither integer nor real; logical handled above
                  nchar = mark - first + 1
                  write (lunout, pos=n) nchar, buffer(first:mark)  ! Text
                  n = n + nibytes + nchar
               case (1)  ! Real, not integer
                  write (lunout, pos=n) rnumber
                  n = n + nrbytes
               case (2)  ! Integer or real, but treat as integer
                  write (lunout, pos=n) inumber
                  n = n + nibytes
            end select

         end if
         first = mark + 2
      end do

      end subroutine treat_data_types

   end program a2b
