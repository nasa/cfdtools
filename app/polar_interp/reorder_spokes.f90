!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program reorder_spokes

!  Reorder a dataset that follows the POLAR_INTERP convention except that it
!  starts with the 6 o'clock spoke rather than the 12 o'clock spoke that the
!  earlier utility expects.  This is easier and less error-prone than modifying
!  POLAR_INTERP to allow for starting at 6 o'clock.
!
!  Assumptions:
!
!  (0) The symmetry plane is assumed to be y = 0.
!  (1) The first spoke contains the nose point that is expected with just one
!      spoke.  It should have y = 0.
!  (2) Input (and output) spokes are ordered from nose to rim.
!  (3) Input spokes are ordered counterclockwise from 6 o'clock;
!      output spokes remain counterclockwise, but from 12 o'clock.
!      The same statements would apply to clockwise ordering.
!  (4) After the first (nose) point, the spokes can be reordered in two copies
!      for a whole body.  A half-body is not so straightforward and is NOT
!      treated initially.  (It would require looking at all the azimuthal
!      angles, along with sorting, to treat one spoke at a time.)
!      A half-body dataset can be detected by whether spokes 2 and nspoke-1
!      have the same azimuthal angle sign, and will cause early termination
!      in this initial implementation.
!  (5) To avoid assuming equal numbers of points per spoke (excluding point 1),
!      azimuthal angles are used to detect the last index in the same half as
!      spoke 2.
!  (5) Any number of functions may accompany the x/y/z coordinates that are
!      expected as columns 1-3.  (No header records are assumed.)
!
!  History:
!
!     08/14/2019  Dinesh had produced a bunch of datasets starting with the
!                 6 o'clock spoke before becoming aware of the 12 o'clock
!                 stipulation.  This is the workaround for the whole-body case.
!
!  Author:  David Saunders, AMA, Inc., NASA Ames Research Center.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use table_type_module  ! Both are in table_io.f90
   use table_io

   implicit none

!  Constants:

   integer, parameter :: lunin  = 1  ! Input table
   integer, parameter :: lunout = 2  ! Output  "
   real,    parameter :: eps = 1.e-7 ! For detecting half or whole body

!  Variables:

   integer :: i, iguess, ios, isplit, ncols, nrows
   logical :: halfbody
   real    :: sign2
   type (table_type) :: tablein, tableout

!  Execution:

   write (*, '(a)', advance='no') 'input dataset: '
   read  (*, '(a)') tablein%filename

   call table_io_read_real (lunin, tablein, ios)  ! Read as reals, not alphas
   if (ios /= 0) then
      write (*, '(2a)') '*** Trouble reading dataset ', trim (tablein%filename)
   end if
   nrows = tablein%nrows
   ncols = tablein%ncols
   write (*, '(3x, a, i4)') &
      'nheader:', tablein%nheader, &
      '# rows: ', nrows,   &
      '# cols: ', ncols
   if (ios /= 0) go to 99

!  Check that the point common to all spokes is point 1 on the symmetry plane.
!  Allow for its not being at z = 0 necessarily, but then the test is not
!  conclusive ...

   if (abs (tablein%values(2,1)) > eps) then  ! y of point 1
      write (*, '(a)') '*** Spoke point 1 does not seem to be at the nose.'
      go to 99
   end if

!  Do the spokes span a half body or a whole body?  We don't know how many
!  spokes we have, but don't want to have to figure that out.

   iguess = nrows/4  ! Beyond spoke 1
   halfbody = sign (1., tablein%values(2,iguess)) == &
              sign (1., tablein%values(2,nrows-iguess))
   if (halfbody) then
      write (*, '(a)') '*** This appears to be a half-body dataset.', &
                       'Handle that differently if it''s ever needed.'
      go to 99
   end if

!  Clone the input table, as the output table has the same dimensions.
!  The logical arguments indicate that the copy's arrays have not been allocated
!  yet and that the table values are to be treated as real, not character:

   tableout%filename = 'reordered.' // trim (tablein%filename)
   call table_io_copy (tablein, tableout, .true., .true., ios)
   if (ios /= 0) then
      write (*, '(a)') '*** Trouble cloning the input table.'
      go to 99
   end if

!  Locate the last index of a spoke point that is in the same half as spoke 2:

   sign2 = sign (1., tablein%values(2,iguess))
   do i = iguess, nrows
      isplit = i
      if (sign (1., tablein%values(2,i)) == sign2 .and. &
               abs (tablein%values(2,i)) > eps) then  ! Avoid y ~ 0.
         cycle
      else
         isplit = isplit - 1
         exit  ! This could surely be done better
      end if
   end do

   write (*, '(a, i5)')  &
      'isplit:', isplit

!  Transcribe the spoke points in two groups:

!! tableout%values(:,1) = tablein%values(:,1) ! Nose; lready done by the cloning
   tableout%values(:,2:isplit) = tablein%values(:,isplit+1:nrows)
   tableout%values(:,isplit+1:nrows) = tablein%values(:,2:isplit)

!  Save the reordering:

   table_io_format = '(3f12.7, 20es15.7)'

   call table_io_write_real (lunout, tableout, table_io_format_sep, ios)
   if (ios /= 0) write (*, '(a)') '*** Trouble writing the output dataset.'

99 continue

   end program reorder_spokes
