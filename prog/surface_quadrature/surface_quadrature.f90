!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program surface_quadrature
!
!  Description:
!
!     Read a structured surface dataset (Tecplot ASCII) of one or more zones
!     with one or more function values per grid point.  Integrate each function
!     over the surface via trapezoidal quadrature, one zone at a time.  Tabulate
!     results for each zone and for the total surface.
!
!     Prompts suffice for the input file and the output file.
!
!  Sample Input Dataset (BLOCK or POINT order):
!
!     TITLE = "sample.dat"
!     VARIABLES = "x" "y" "z"
!     "function  1"
!     "function  2"
!     "function  3"
!     ZONE T = "Gridded data"
!     I = 49, J = 49, K=1, ZONETYPE=ORDERED, DATAPACKING=BLOCK
!     0.000000E+00   3.724260E-04   1.418648E-03   3.032083E-03   5.106146E-03
!     7.534256E-03   1.020982E-02   1.302626E-02   1.627440E-02   2.059799E-02
!      :              :              :              :              :
!
!  Sample Output (Screen and Output File):
!
!     Zone   Function 1   Function 2   Function 3  ...
!        1      1.23456      2.98765      3.45678  ...
!       [2      2.34567       :            :
!     Total     3.58023       :            :]
!
!  05/02/2017  D.A.Saunders  Initial implementation, to deal with mass-loss
!                            data from NEQAIR via POLAR_INTERP.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_header_structure  ! See Tecplot_io.f90 for these modules
   use grid_block_structure
   use tecplot_io_module

   implicit none

!  Constants:

   integer, parameter :: &
      lunsurf = 1,       &  ! Input structured surface solution  (Tecplot ASCII)
      luntab  = 2,       &  ! Output tabulation
      lunkbd  = 5,       &  ! Keyboard inputs
      luncrt  = 6,       &  ! Screen outputs
      ndim    = 3

   real, parameter :: &
      zero = 0.

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

!  Variables:

   integer :: &
      ifun, ios, iz, lun, nf, ni, nj, nzone

   real :: &
      total_area

   logical :: &
      cr, eof

   character (80) :: &
      filename

   real, allocatable, dimension (:) :: &
      area, total_quad

   real, allocatable, dimension (:,:) :: &
      quad

!  Derived data types:

   type (grid_header) :: &
      header_surface

   type (grid_type), pointer, dimension (:) :: &
      surface               ! (x,y,z,f) for the input surface (1 zone or more)

!  Execution:

   call reads (luncrt, 'Input surface file (Tecplot structured dataset): ', &
               lunkbd, header_surface%filename, cr, eof)
   if (eof) go to 99

   call reads (luncrt, 'Output file of tabulated quadratures: ', &
               lunkbd, filename, cr, eof)
   if (eof) go to 99

   ios = 1  ! Verbose mode
   header_surface%ndim = ndim
   header_surface%formatted = true

   call Tec_header_read (lunsurf, header_surface, surface, ios)
   if (ios /= 0) go to 99

   nzone = header_surface%nblocks
   nf    = header_surface%numq

   open (luntab, file=filename, status='unknown')

!  Process one zone at a time:

   do lun = luntab, luncrt, luncrt - luntab
      write (lun, '(a)') &
         '   Zone            Area      Function 1      Function 2  ...'
   end do

   allocate (area(nzone), quad(nf,nzone), total_quad(nf))
   area(:) = zero;  quad(:,:) = zero;  total_quad(:) = zero
   total_area = zero
   
   do iz = 1, nzone

      call Tec_block_allocate (surface(iz), ndim, nf, ios)
      if (ios /= 0) go to 99

      call Tec_block_read (lunsurf, header_surface, surface(iz), ios)
      if (ios /= 0) go to 99

      ni = surface(iz)%ni
      nj = surface(iz)%nj

      call trapezoidal_quadrature_3d (ni, nj, nf, &
                                      surface(iz)%x, surface(iz)%y, &
                                      surface(iz)%z, surface(iz)%q, &
                                      area(iz), quad(:,iz))
      total_area = area(iz) + total_area
      total_quad(:) = quad(:,iz) + total_quad(:)

      do lun = luntab, luncrt, luncrt - luntab
         write (lun, '(i7, 40es16.8)') iz, area(iz), quad(:,iz)
      end do

   end do

   if (nzone > 1) then
      do lun = luntab, luncrt, luncrt - luntab
         write (lun, '(a, 40es16.8)') 'Totals:', total_area, total_quad(:)
      end do
   end if

99 continue

   end program surface_quadrature
