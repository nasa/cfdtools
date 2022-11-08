!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine transcribe (lunin, lunout, ios)

!  Transcribe one text file to another, as commonly needed for echoing a control
!  file to a log file. Both files are assumed to be open on the indicated units,
!  ready for the transcription. The input file is rewound upon return, ready for
!  reading as though there was no transcription.
!
!  History:
!     02/20/2020  D.A.Saunders  Belated implementation of an obvious utility.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: lunin, lunout  ! Input and output logical units
   integer, intent (out) :: ios            ! 0 means no problem;
                                           ! 1 means a read error
!  Local constants:

   integer, parameter :: maxchar = 132     ! Sensible limit on input lines

!  Local variables:

   character (maxchar) :: buffer

!  Execution:

   do ! Until EOF
      read (lunin, '(a)', iostat=ios) buffer
      if (ios /= 0) exit

      write (lunout, '(1x, a)') trim (buffer)
   end do

   if (ios < 0) then  ! Normal EOF
       ios = 0
   else
       ios = 1
   end if

   rewind (lunin)

   end subroutine transcribe
