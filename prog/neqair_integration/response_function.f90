!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program response_function
!
!
!  Apply a radiometer response function to spectral radiance output from NEQAIR.
!  Save the result in the same form as the NEQAIR output (that of intensity.out
!  or intensity_scanned.out but with just two columns).
!
!  Any header lines in either of the input datasets are ignored.  The response
!  function is expected to span a relatively small wavelength range, and assumed
!  to be zero elsewhere.  Wavelengths should be angstrom units for both inputs.
!  The column for spectral radiance from NEQAIR is prompted for, with wavelength
!  assumed to be in column 1.
!
!  06/12/2018  D.A.Saunders  Initial interpretation of a request by Todd White,
!                            with Schiaparelli data in mind.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use table_type_module  ! Both are in table_io.f90
   use table_io

   implicit none

!  Constants:

   integer, parameter :: lunin1 = 1  ! Spectral radiance input from (say) NEQAIR
   integer, parameter :: lunin2 = 2  ! Response function to be applied
   integer, parameter :: lunout = 3  ! Output function vs. wavelength
   integer, parameter :: lunkbd = 5  ! Keyboard inputs
   integer, parameter :: luncrt = 6  ! Screen outputs
   real,    parameter :: arrow  = 1. ! Ascending wavelengths assumed
   character (1), parameter :: blank  = ' '
   character (1), parameter :: method = 'M' ! Monotonic spline interpolation of
                                            ! the response function

!  Variables:

   integer :: icol, ier, j, nr
   real    :: flambda, lambda, unused
   real, allocatable, dimension (:) :: fx, fy
   logical :: cr, eof, new
   type (table_type) :: response, &  ! Response fn. vs. wavelength
                        spectrum, &  ! Predicted spectral radiance data
                        result       ! Output (lambda, response) pairs
!  Execution:

   spectrum%filename = 'intensity_scanned.out'
   call reads (luncrt, 'Spectral radiance file? <cr>=intensity_scanned.out: ', &
               lunkbd, spectrum%filename, cr, eof)
   if (eof) go to 99

   icol = 2
   call readi (luncrt, &
            'Col. 1 = wavelength assumed; spectral radiance col.? (<cr>=2): ', &
               lunkbd, icol, cr, eof)
   if (eof) go to 99

   response%filename = 'response_function.dat'
   call reads (luncrt, 'Response function file? <cr>=response_function.dat: ', &
               lunkbd, response%filename, cr, eof)
   if (eof) go to 99

   call table_io_read_real (lunin2, response, ier)
   if (ier /= 0) go to 99

   write (luncrt, '(/, 2a)') ' Reading ', trim (spectrum%filename)

   call table_io_read_real (lunin1, spectrum, ier)
   if (ier /= 0) go to 99

   write (luncrt, '(a, i9)') &
      ' # rows found:', spectrum%nrows, &
      ' # cols found:', spectrum%ncols

!  Set up the results as a table in order to use the table write utility:

   result%filename = 'response.dat'
   result%nrows    = spectrum%nrows
   result%ncols    = 2
   table_io_format = '(f14.4, es13.6)'

!  Apply the response function to each spectral radiance.
!  Searching the response function requires transposed data.

   nr = response%nrows
   allocate (fx(nr), fy(nr), result%values(2,result%nrows))
   fx(:) = response%values(1,:)
   fy(:) = response%values(2,:)
   deallocate (response%values)

   new = .true.
   do j = 1, spectrum%nrows
      lambda = spectrum%values(1,j)
      call lcsfit (nr, fx, fy, new, method, 1, lambda, flambda, unused)
      new = .false.
      result%values(1,j) = lambda
      result%values(2,j) = flambda*spectrum%values(icol,j)
   end do

   write (luncrt, '(2a)') ' Writing results to ', trim (result%filename)

   call table_io_write_real (lunout, result, blank, ier)

99 continue

   end program response_function
