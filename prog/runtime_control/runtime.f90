
   module runtime

!     Just trying out some ideas ...

      implicit none

      integer, public, parameter :: mxcontrols = 50

      integer,   public :: ncontrols
      integer,   public :: action_value
      character, public :: action_keyword * 16  ! Upper case; e.g., @ITER[ATION]
      character, public :: control_names(mxcontrols) * 16
      character, public :: control_values(mxcontrols) * 16

      public :: look_for_controls

      private

      integer,   parameter :: mx_token_length = 16
      integer,   parameter :: mx_tokens_per_line = 16
      character, parameter :: comment * 1 = '#'
      character, parameter :: delimiters * 5 = ' :,=' // char (9)  ! tab

      character, dimension (mx_tokens_per_line) :: tokens * (mx_token_length)

      contains

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine look_for_controls (lun, filename)

!        Open the indicated file, scan for an action keyword line containing
!        an associated value not less than action_value, and return all
!        keyword/value string pairs up to the next action keyword if any.
!        The action keyword initially envisaged is '@iteration', and its value
!        is assumed to be an integer.

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Arguments:

         integer,   intent (in)  :: lun
         character, intent (in)  :: filename * (*)

!        Local constants:

         character, parameter :: routine * 26 = 'Run-time control warning: '

!        Local variables:

         integer   :: ios, itoken, ivalue, last, n, n1, ntokens
         logical   :: active
         character :: buffer * 132, keyword * 16, keyvalue * 16

!        Procedures:

         external :: getline ! Handles leading and trailing comments
         external :: token2  ! Tokenizes a string

!        Execution:

         ncontrols = 0
         active = .false.

         open (lun, file=filename, status='old', iostat=ios)

         do while (ios == 0 .and. ncontrols < mxcontrols)

            call getline (lun, comment, buffer, last, ios)
!! write (*,*) 'last: ', last
            if (ios < 0) exit  ! EOF

            if (ios > 0 ) then  ! Shouldn't be possible
            end if

            if (last == 0) cycle  ! Blank line or comment

            ntokens = mx_tokens_per_line  ! Most allowed for

            call token2 (buffer(1:last), delimiters, ntokens, tokens)

!! write (*,*) 'ntokens: ', ntokens
            n1 = 1  ! First token number

!           Special handling of the action keyword:

            if (tokens(1)(1:1) == action_keyword(1:1)) then

!              May want to handle, say, '@ iter 1000' as well as '@iter 1000'
!              Don't bother yet.

               if (tokens(1)(1:4) == action_keyword(1:4)) then

                  if (active) exit  ! Ignore a set of control intended for later

                  active = .true.

                  if (ntokens < 2) exit  ! Should issue a warning

                  read (tokens(2), *, iostat=ios) ivalue
!! write (*,*) 'ivalue: ', ivalue
                  if (ios /= 0) then
                     write (*, '(/, 2a)') &
                        routine, 'Trouble parsing this line:', buffer(1:last)
                     exit
                  end if
!! write (*,*) 'action-value: ', action_value
                  if (ivalue < action_value) then  ! Keep looking
                     active = .false.
                     cycle
                  end if

                  n1 = 3

               end if

            end if

            if (.not. active) cycle  ! Skip earlier action items

!           Transfer any [remaining] controls on this line as (keyword, value)s:

            n = ntokens  ! In case they don't come in pairs

            do itoken = n1, n, 2
               ncontrols = ncontrols + 1
!! write (*,*) 'ncontrols: ', ncontrols
               if (ncontrols <= mxcontrols) then
                  control_names(ncontrols) = tokens(itoken)
                  if (itoken < ntokens) then
                     control_values(ncontrols) = tokens(itoken + 1)
                  else
                     ! Warn of missing keyword value?
                     ntokens = ntokens - 1
                  end if
               end if
            end do

         end do  ! Next input line

         close (lun, iostat=ios)

         end subroutine look_for_controls

   end module runtime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program runtime_controls

!  Try some ideas for how to modularize periodic run-time control changes.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use runtime

   implicit none

   integer, parameter :: lun = 1
   character, parameter :: filename * 15 = 'runtime.control'

   integer :: ios, ivalue, l1, l2, n
   real    :: rvalue

   action_keyword = '@ITERATION'
   action_value   = 1000

   call look_for_controls (lun, filename)

   write (*, '(a, i3)') 'ncontrols:', ncontrols

   do n = 1, ncontrols


      l1 = len_trim (control_names(n))
      l2 = len_trim (control_values(n))

      write (*, '(4a)') &
         'Keyword and value strings: ', control_names(n)(1:l1), &
         '  ', control_values(n)(1:l2)

      select case (control_names(n)(1:l1))

         case ('KEY1')

            read (control_values(n)(1:l2), *, iostat=ios) ivalue

            if (ios /= 0) then
               write (*, '(2a)') ' *** Bad value: ', control_values(n)(1:l2)
            else
               write (*, '(a, i6)', iostat=ios) 'Integer value: ', ivalue
            end if

         case ('KEY2')

            read (control_values(n)(1:l2), *, iostat=ios) rvalue

            if (ios /= 0) then
               write (*, '(2a)') ' *** Bad value: ', control_values(n)(1:l2)
            else
               write (*, '(a, f10.6)', iostat=ios) 'Real value: ', rvalue
            end if

         case ('KEY3')

            read (control_values(n)(1:l2), *, iostat=ios) ivalue

            if (ios /= 0) then
               write (*, '(2a)') ' *** Bad value: ', control_values(n)(1:l2)
            else
               write (*, '(a, i6)', iostat=ios) 'Integer value: ', ivalue
            end if

         case default

            write (*, '(2a)') 'Unknown keyword: ', control_names(n)(1:l1)

      end select

   end do

   end program runtime_controls
