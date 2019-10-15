!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program merge_zones
!
!  Description:
!
!     This utility was prompted by the need to mix parts of results from two
!     turbulent flow calculations, namely a Baldwin-Lomax solution on a capsule
!     forebody and an SST solution on the aft body.  In general, it merges the
!     indicated zones from a first file with zones indicated from a second file.
!     If the second file is suppressed, it could also be used to reorder the
!     zones of the first file.  The files are expected to be in Tecplot format,
!     and ASCII, since we cannot read Tecplot binary (plt) files.
!
!  Control:
!     A few prompts suffice (no control file).  Following each file name, any
!     reasonable list of integers on a single line define the zones to be trans-
!     ferred in the appropriate order to the output file, which is also in ASCII
!     Tecplot form.  The list(s) can be of any form accepted by the RDLIST
!     utility (see the progtools library).  E.g. 1:8 for the first 8 zones.
!
!  Method:
!     Since the zone order(s) may be shuffled, it is expedient to read all of
!     each file rather than try to process one zone at a time.  The output file
!     is written header first then one zone at a time, though.  For this to make
!     sense, the variable lists must be identical in the two files.
!
!  History:
!     10/13/2018  D.A.Saunders  Initial implementation.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_header_structure  ! See tecplot_io.f90
   use grid_block_structure   !   "     "     "
   use tecplot_io_module      !   "     "     "

   implicit none

!  Local constants:

   integer, parameter :: lunin1 = 1   ! First input file
   integer, parameter :: lunin2 = 2   ! Second  "    "
   integer, parameter :: lunout = 3   ! Output result
   integer, parameter :: lunkbd = 5   ! Keyboard inputs
   integer, parameter :: luncrt = 6   ! Screen prompts and diagnostics
   integer, parameter :: mxlist = 512 ! List length limit

!  Local variables:

   integer :: i, ios, iv, iz, n1, n2, nzone1
   integer, dimension (mxlist) :: list1, list2
   logical :: cr, eof
   type (grid_header) :: header1, header2
   type (grid_type), pointer, dimension (:) :: file1, file2

!  Execution:

   write (luncrt, '(/, a)', advance='no') ' 1st file to transcribe zones from: '
   read  (lunkbd, '(a)') header1%filename
   header1%ndim = -1   ! Automates detection of 2D or 3D
   header1%formatted = .true.
   ios = 0  ! 1 means verbose read
   call Tecplot_read (lunin1, header1, file1, ios)
   nzone1 = header1%nblocks  ! We'll reuse header1 for the output file
   if (ios /= 0) go to 99

   n1 = mxlist
   call rdlist (luncrt, 'List of zones to transcribe (any order): ', &
                lunkbd, n1, list1)
   if (n1 <= 0 .or. n1 > nzone1) then
      write (luncrt, '(a, 2i5)') ' Bad list length.  Limits: ', 1, nzone1
      go to 99
   end if

   write (luncrt, '(a)', advance='no') &
      ' 2nd file to transcribe zones from: [<cr> = ctrl D = none] '
   call reads (luncrt, '2nd file to transcribe zones from: [ctrl D = none] ', &
               lunkbd, header2%filename, cr, eof)
   if (cr .or. eof) then
      n2 = 0
   else
      header2%formatted = .true.
      header2%ndim = -1  ! Automate detection
      ios = 0  ! 1 means verbose read
      call Tecplot_read (lunin2, header2, file2, ios)
      if (ios /= 0) go to 99

!     Make sure the variables are the same:

      if (header2%ndim /= header1%ndim .or. &
          header2%numq /= header1%numq) then
          write (luncrt, '(a)') &
             ' File dimensions or variable counts don''t match.'
          go to 99
      end if

      do iv = 1, header1%numq
         if (header2%varname(iv) /= header1%varname(iv)) then
             write (luncrt, '(a, i4)') ' Name mismatch: variable #:', iv
             go to 99
         end if
      end do

      n2 = mxlist
      call rdlist (luncrt, 'List of zones to transcribe (any order): ', &
                   lunkbd, n2, list2)
      if (n2 <= 0 .or. n2 > header2%nblocks) then
         write (luncrt, '(a, 2i5)') ' Bad list length.  Limits: ', 1, &
            header2%nblocks
         go to 99
      end if
   end if

   write (luncrt, '(a)', advance='no') ' Output file name: '
   read  (lunkbd, '(a)', iostat=ios) header1%filename
   if (ios < 0) go to 99

   call readi (luncrt, 'Output datapacking? [0=point order, else 1:] ', &
               lunkbd, header1%datapacking, cr, eof)

   header1%nblocks = n1 + n2
   call Tec_header_write (lunout, header1, ios)  ! Open and write header
   if (ios /= 0) go to 99

!  Transcribe the indicated zones:

   do i = 1, n1
      iz = list1(i)
      call Tec_block_write (lunout, header1, file1(iz), ios)
      if (ios /= 0) go to 99
   end do

   do i = 1, n2
      iz = list2(i)
      call Tec_block_write (lunout, header2, file2(iz), ios)
      if (ios /= 0) go to 99
   end do

   close (lunout)

99 continue

   end program merge_zones
