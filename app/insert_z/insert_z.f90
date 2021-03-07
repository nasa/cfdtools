!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program insert_z
!
!  Description:
!
!  This is a specialized utility prompted by a parametric study of axisymmetric
!  flow solutions where one parameter is varied at a time and can serve as the
!  third dimension for turning line plots into contour plots. It reads an x-y
!  line dataset (Tecplot ASCII, e.g., x, y, si, qw along a surface) and inserts
!  variable 3 as "z" with a constant parameter value that is prompted for.
!  Use of the name z facilitates subsequent Tecplot contouring of concatenated
!  datasets.
!
!  Part of the motivation is to turn line data into surface data so that the
!  BUMP_FACTORS utility can be applied to judiciously concatenated and edited
!  datasets (in the absence of a line form of BUMP_FACTORS).  The associated
!  baseline "surface" would be constructed by replicating the baseline line
!  dataset and inserting "z" values in the copies as for the datasets with
!  the parameter truly varying.  Then concatenating and a little editing
!  produces two surface datasets suited to BUMP_FACTORS.
!
!  04/23/12  David Saunders  Initial implementation, for a sphere-cone study.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_header_structure  ! All of these are part of Tecplot_io.f90
   use grid_block_structure
   use tecplot_io_module

   implicit none

!  Constants:

   integer, parameter :: &
      lunin  = 1,  &     ! Input xy-space dataset (Tecplot ASCII, 1+ zones)
      lunout = 2,  &     ! Output dataset (Tecplot ASCII) with "z" substituted
      lunkbd = 5,  &     ! Keyboard inputs
      luncrt = 6         ! Prompts and diagnostics

   real, parameter ::    &
      one  = 1.0,        &
      zero = 0.0

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

!  Variables:

   integer :: &
      ib, ios, nblocks, numf, nvars

   real :: &
      z                  ! Desired constant value to substitute for variable 3

   logical :: &
      formatted

   character :: &
      filein*64, fileout*64

!  Composite data types:

   type (grid_header) :: &
      Tecplot_header_in, Tecplot_header_out

   type (grid_type), pointer, dimension (:) :: &
      dataset_in,     &  ! (x,y,f) for one or more input zones
      dataset_out        ! (x,y,z,f) for the output zone(s)

!  Execution:

   write (luncrt, '(a)', advance='no') &
      'Input x-y line dataset (Tecplot ASCII, 1 or more zones): '
   read (luncrt, *) filein

   write (luncrt, '(a)', advance='no') 'Output dataset name (Tecplot ASCII): '
   read (luncrt, *) fileout

   write (luncrt, '(a)', advance='no') &
      'Constant "z" value to insert as variable 3: '
   read (luncrt, *) z

   Tecplot_header_in%filename  = filein
   Tecplot_header_in%formatted = true
   Tecplot_header_in%ndim      = 2

   ios = 1  ! Turn on verbose mode

   call Tecplot_read (lunin, Tecplot_header_in, dataset_in, ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Trouble reading input Tecplot dataset ', &
         trim (filein)
      go to 99
   end if

   nblocks = Tecplot_header_in%nblocks
   numf    = Tecplot_header_in%numq

!  Set up the output dataset:

   Tecplot_header_out%filename  = fileout
   Tecplot_header_out%formatted = true
   Tecplot_header_out%ndim      = 3
   Tecplot_header_out%numq      = numf
   Tecplot_header_out%nblocks   = nblocks
   Tecplot_header_out%title     = Tecplot_header_in%title

   nvars = 3 + numf
   allocate (Tecplot_header_out%varname(nvars))
   Tecplot_header_out%varname(1:2) = Tecplot_header_in%varname(1:2)
   Tecplot_header_out%varname(3)   = 'z'
   Tecplot_header_out%varname(4:)  = Tecplot_header_in%varname(3:)

   allocate (dataset_out(nblocks))

   do ib = 1, nblocks

      call clone_zone (dataset_in(ib), 3, numf, dataset_out(ib), ios)

      dataset_out(ib)%x = dataset_in(ib)%x
      dataset_out(ib)%y = dataset_in(ib)%y
      dataset_out(ib)%z = z
      dataset_out(ib)%q = dataset_in(ib)%q

   end do

!  Save the adjusted dataset:

   call Tecplot_write (lunout, Tecplot_header_out, dataset_out, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble writing the adjusted dataset.'
   end if

   call deallocate_header (Tecplot_header_in)
   call deallocate_header (Tecplot_header_out)
   call deallocate_blocks (1, nblocks, 2, numf, dataset_in,  ios)
   call deallocate_blocks (1, nblocks, 3, numf, dataset_out, ios)

99 continue

   end program insert_z
