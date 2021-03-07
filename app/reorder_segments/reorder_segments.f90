!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program reorder_segments
!
!  This is a rather specialized utility prompted by the need to manipulate a
!  centerline dataset from a full-body flow solution on a CAPSULE_GRID grid.
!  The dataset is first extracted with POSTFLOW as 8 zones then turned into a
!  single zone via the XLINE utility, which is used to remove all the " zone"
!  lines).  Actually, this utility now has the option to remove the "zone" lines
!  itself.
!
!  The problem with merging as a single zone is that the line segments don't fit
!  end-to-end as needed for plotting as a single curve.  Simple lists of zone
!  dimensions, reverse-order switches, and output order of the adjusted segments
!  serve to drive the utility.
!
!  The control file name is hard-coded as 'reorder_segments.inp' as opposed to
!  being read from standard input so that the file to be processed can be
!  prompted for.  The input file name is prepended with 'reordered.' for the
!  output file name.
!
!  Sample reorder_segments.inp:
!
!     1   ! nheader    = # file header lines to transmit unchanged
!     1   ! nzonelines = # zone header lines per zone to suppress in the output
!     5   ! ncolumns   = # columns in all remaining lines
!     8   ! nzones     = # line segments (zones) to be shuffled
!     17 17 90 90 111 111 17 17  ! Input segment lengths
!     T  F  F  T  F   T   T  F   ! T = reverse segment order
!     2  3  6  7  8   5   4  1   ! Output order of segments
!
!  05/07/13  D.A.Saunders   Initial implementation for Mars InSight work.
!  05/09/13    "      "     Added the option to suppress zone header lines here.
!                           Enter nzonelines = 0 if there are none to delete.
!
!  Author:  David Saunders  ERC, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      lunin  = 1,        &
      lunout = 2

!  Variables:

   integer :: &
      i, i1, i2, inc, lmax, m, n, ncolumns, nheader, nzonelines, nzones
   character (80) :: &
      filename
   integer, allocatable, dimension (:) :: &
      length, output_order
   logical, allocatable, dimension (:) :: &
      reverse
   character (256), allocatable, dimension (:,:) :: &
      dataset

!  Execution:

   open  (lunin, file='reorder_segments.inp', status='old')
   read  (lunin, *) nheader
   read  (lunin, *) nzonelines
   read  (lunin, *) ncolumns
   read  (lunin, *) nzones
   allocate (length(nzones), reverse(nzones), output_order(nzones))
   read  (lunin, *) length
   read  (lunin, *) reverse
   read  (lunin, *) output_order
   close (lunin)

   lmax = maxval (length)
   allocate (dataset(lmax,nzones))  ! Storage for nzones of size biggest zone

   write (*, '(a)', advance='no') 'Dataset file name: '
   read  (*, *) filename
   open  (lunin, file=filename, status='old')
   open  (lunout, file='reordered.'//trim (filename), status='unknown')

   do i = 1, nheader
      read (lunin, '(a)') dataset(1,1)  ! Use as a buffer
      write (lunout, '(a)') trim (dataset(1,1))
   end do

   do n = 1, nzones
      do i = 1, nzonelines
         read (lunin, '(a)') ! Skip zone header lines
      end do
      do i = 1, length(n)
         read (lunin, '(a)') dataset(i,n)
      end do
   end do

   close (lunin)

   do n = 1, nzones
      m = output_order(n)
      i1  = 1
      i2  = length(m)
      inc = 1
      if (reverse(m)) then
         i1 = i2
         i2 = 1
         inc = -1
      end if
      do i = i1, i2, inc
         write (lunout, '(a)') trim (dataset(i,m))
      end do
   end do

   close (lunout)
   deallocate (length, reverse, output_order, dataset)

   end program reorder_segments
