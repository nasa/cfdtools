   program test_interface

!  Test the Fortran-callable interface to Gary Allen's equilibrium gas
!  composition C code.
!  Also, by generating tables of outputs that can be contoured, ascertain
!  the precision available.
!
!  12/12/05  D. A. Saunders  Initial implementation.
!  09/09/11   "  "  "  "     Revived to test CO2 cases.

   implicit none

!  Constants:

   integer, parameter :: &
      luncrt = 6, &
      lunkbd = 5, &
      lunin  = 1, &  ! Input file of the type read by Gary's code
      lunout = 2     ! Structured 2-space results in ASCII Tecplot "POINT" form.

   character, parameter :: &
      file_out * 7 = 'gas.dat'

!  Variables:

   integer :: &
      i, j, ni, nj, nspecies, numf, status

   real :: &
      dv1, dv2, v1, v1min, v2, v2min, &
      pressure, density, temperature, enthalpy, entropy, speed_of_sound

   real, allocatable, dimension (:) :: &
      mole_fractions

   character :: &
      file_in * 32, &
      specified_state * 5                ! One more input to the C code ...

   character, allocatable, dimension (:) :: &
      names * 4                          ! ... and another

!  Execution:

   status = -1
   do while (status /= 0)
      write (luncrt, '(/, a)', advance='no') ' Input file: '
      read  (lunkbd, *) file_in
      open  (lunin, file=file_in, status='old', iostat=status)
   end do

   read (lunin, *) specified_state, v1min, v2min
   nspecies = 0
   do while (status == 0)
      nspecies = nspecies + 1
      read (lunin, '(a)', iostat=status)
   end do
   nspecies = nspecies - 1

   write (luncrt, '(a, i3)') ' # species found: ', nspecies

   allocate (names(nspecies), mole_fractions(nspecies))
   rewind (lunin)
   read   (lunin, *)
   do i = 1, nspecies  ! Avoid trailing comments
      read (lunin, *) names(i), mole_fractions(i)
   end do
   close  (lunin)

   write (luncrt, '(2a, 1p, 2e13.6)') &
      ' Specified state: ', specified_state, v1min, v2min
   write (luncrt, '(a)', advance='no') &
      ' Enter n1, dv1 and n2, dv2 for a grid starting as (v1,v2) [2 pairs]: '
   read  (lunkbd, *) ni, dv1, nj, dv2

!  Set up the output file in Tecplot ASCII format, POINT order:

   open (lunout, file=file_out, status='unknown')

   write (lunout, '(a)') &
      'Gas mixture calculations via Gary Allen''s lewis.c', &
      'variables=v1,v2,pressure,density,temperature,enthalpy,entropy,a'
   write (lunout, '(a, i4, a, i4)') &
      'zone t="table" F=point, i=', ni, ' j=', nj

!  Evaluate all the implied cases of gas properties:

   do j = 1, nj

      v2 = v2min + real (j-1) * dv2

      do i = 1, ni

         v1 = v1min + real (i-1) * dv1

         call equilibrium_composition (specified_state, v1, v2,                &
                                       nspecies, names, mole_fractions,        &
                                       pressure, density, temperature,         &
                                       enthalpy, entropy, speed_of_sound,      &
                                       status)
         if (status /= 0) then
            write (luncrt, '(a, i4, a, 2i4, a, 1p, 2e13.6)') &
               ' lewis status:', status, ' i, j:', i, j,    &
               ' v1, v2:', v1, v2
            go to 99
         end if

         write (lunout, '(1p, 8e13.6)') v1, v2, &
            pressure, density, temperature,     &
            enthalpy, entropy, speed_of_sound

      end do

   end do

   close (lunout)

99 continue

   end program test_interface
