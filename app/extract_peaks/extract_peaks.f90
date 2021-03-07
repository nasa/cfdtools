!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program extract_peaks

!  For a list of files of column data, extract peak values and tabulate them.
!  More precisely:  if the first ndim columns are coordinates (time or spatial)
!  followed by nfun columns of function data, write nfun peak values to one
!  line per input file.  This could easily be extended to (say) writing the
!  coordinate(s) along with the peak value to nfun output files if the peak
!  locations are also of interest.  Initially, they're not.
!
!  Some prompts handle the number of header lines to skip, and the names of the
!  functions to include in the output files.
!
!  Sample Input File of Files to Process:
!
!     DPLR-AoA=0-pqtau-laminar-MHR-t655.dat
!     DPLR-AoA=0-pqtau-laminar-MHR-t665.dat
!     DPLR-AoA=0-pqtau-laminar-MHR-t674.dat
!     DPLR-AoA=0-pqtau-laminar-MHR-t680.dat
!     DPLR-AoA=0-pqtau-laminar-MHR-t687.dat
!     DPLR-AoA=0-pqtau-laminar-MHR-t695.dat
!
!  Representative Format of One Such File (1 Header Line, 2 Dims., 3 Functions):
!
!     variables=x,y,p,qw,tau
!      0.000000000  0.000000000  1.55635363E+00  4.22335080E+01  9.34081158E-01
!      0.000122464  0.008980126  1.55635363E+00  4.22335322E+01  9.34081693E-01
!      0.000610847  0.026908363  1.55164269E+00  4.19596878E+01  2.79743491E+00
!      0.001584453  0.044774007  1.54491638E+00  4.15442151E+01  4.43967952E+00
!       :            :            :               :               :
!
!  Corresponding Output (All Function Peaks/No Coordinates), Tab-delimited:
!
!     x   y   p   qdot   tau                     kPa      W/cm^2          Pa
!     DPLR-AoA=0-pqtau-laminar-MHR-t655   2.3456E+01  1.1234E+00  1.4142E+02
!     DPLR-AoA=0-pqtau-laminar-MHR-t665   2.5678E+01  1.5432E+00  2.1234E+02
!                                     :        :           :           :
!
!  History:
!
!  04/12/13  D.A.Saunders  Initial implementation, for analyzing centerline data
!                          from 48 flow solutions for the Mars InSight project.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      luncases = 1,      &  ! File of file names to process
      lundat   = 2,      &  ! One such file
      lunout   = 3,      &  ! Output tabulation
      lunkbd   = 5,      &  ! Keyboard inputs
      luncrt   = 6,      &  ! Screen prompts, error messages
      maxcol   = 12,     &  ! Limit on # columns in 3 write formats below
      maxlen   = 64         ! Limit on file and variable name lengths

!  Variables:

   integer :: &
      ios, l, m, n, ndata, ndim, nfun, nheader

   integer :: &
      location(1)

   integer, allocatable, dimension (:) :: &
      location_peak

   character (1)  :: &
      tab

   character (maxlen) :: &
      casename, filename

   real, allocatable, dimension (:) :: &
      peak_function

   real, allocatable, dimension (:,:) :: &
      coordinate, function, peak_coordinate

   character (maxlen), allocatable, dimension (:) :: &
      coordinate_name, function_name, unit_name

!  Execution:

   write (luncrt, '(/, a)', advance='no') 'File of file names to process: '
   read  (lunkbd, '(a)') filename

   ios = -1
   do while (ios /= 0)
      open (luncases, file=filename, status='old', iostat=ios)
      if (ios /= 0) write (luncrt, '(a)') 'File not found; try again.'
   end do

   write (luncrt, '(a)', advance='no') 'nheader, ndim and nfun: '
   read  (lunkbd, *) nheader, ndim, nfun

   n = ndim + nfun
   if (n > maxcol) then
      write (luncrt, '(a, i4)') 'Please recompile with maxcol >=', n
      go to 99
   end if

   allocate (coordinate_name(ndim), function_name(nfun), unit_name(nfun))
   write (luncrt, '(a)', advance='no') 'Coordinate & function names: '
   read  (lunkbd, *) coordinate_name, function_name
   write (luncrt, '(a)', advance='no') 'Units of functions (only): '
   read  (lunkbd, *) unit_name

   write (luncrt, '(a)', advance='no') 'Output file of peak values: '
   read  (lunkbd, '(a)') filename
   open  (lunout, file=filename, status='unknown')

   tab = char (9)

   allocate (peak_coordinate(ndim,nfun), peak_function(nfun))

   write (lunout, '(a, 12(3x,a))', advance='no') &
      (trim (coordinate_name(l)), l = 1, ndim), &
      (trim (function_name(m)), m = 1, nfun)
   write (lunout, '(12(2a))') (tab, trim (unit_name(m)), m = 1, nfun)

!  Loop over cases:

   do  ! Until case list EOF

      read (luncases, *, iostat=ios) casename
      if (ios < 0) exit

      open (lundat, file=casename, status='old', iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(2a)') 'Cannot open file ', trim (casename)
         exit
      end if

!     Count the number of lines and allocate storage accordingly:

      do l = 1, nheader
         read (lundat, '(a)')
      end do

      ndata = 0
      do  ! Until EOF
         read (lundat, '(a)', iostat=ios)
         if (ios < 0) exit
         ndata = ndata + 1
      end do

      rewind (lundat)
      do l = 1, nheader
         read (lundat, '(a)')
      end do

      allocate (coordinate(ndim,ndata), function(nfun,ndata))

      read  (lundat, *) (coordinate(:,n), function(:,n), n = 1, ndata)
      close (lundat)

      do m = 1, nfun
         location = maxloc (function(m,:))      ! Index of peak for mth function
         l = location(1)                        ! Location is a 1-D array
         peak_coordinate(:,m) = coordinate(:,l) ! Coordinates of this peak
         peak_function(m)     = function(m,l)   ! Peak value of this function
      end do

      l = index (casename, '.') - 1 ! Suppress any file extension in case name
      if (l < 0) l = len_trim (casename)

      write (lunout, '(a, 12(a1, es20.12))') &
         casename(1:l), (tab, peak_function(m), m = 1, nfun)

      deallocate (coordinate, function)

   end do  ! Next case

99 deallocate (coordinate_name, peak_coordinate, peak_function)
   close (luncases)
   close (lunout)

   end program extract_peaks
