!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program point_source

!  This extra step after Stardust_Integration processes many intensity.out files
!  from NEQAIR was requested by the Hayabusa 2 observation team to remove the
!  square of the observation distance that has been incorporated in the output
!  integrations.dat file.  This file contains integrations (one per wavelength)
!  w.r.t. solid angle subtended at the sensor location over the elements in the
!  pixel array used to define the lines of sight through the flow field that
!  NEQAIR has processed.  Its units are spectral irradiance, W/cm^2-micron.
!
!  The slant range used for the integration is read from the same 'lines.inp'
!  control expected by Stardust_Integrations (on the viewing angle line).  Those
!  units should be meters, so in view of the W/cm^2-micron units above, the
!  distance is converted to cm before being applied, and the output units are
!  W/micron.  [Later:] No: Brett Cruden says the output after multiplying by
!  distance squared (Spectral Radiant Intensity) has units W/sr-micron.
!
!  Note that the viewing angle remains inherent in these predicted values
!  (which also ignore atmospheric absorption over the viewing distance)
!  so a range of datasets from the entry path and a representative observation
!  flight path is needed to provide upper bounds for some other flight path.
!
!  See logical unit descriptions below for input & out file descriptions (two
!  columns with no header records).
!
!     07/06/2020  DAS  Initial implementation.
!     07/13/2020   "   Nomenclature corrections.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      line_to_read = 8,  & ! ... within lines.inp control file for view/distance
      lunin        = 1,  & ! Input integrations.dat (angstrom, W/cm^2-micron)
      lunout       = 2     ! Output point_source.dat (nm, W/sr-micron) pairs

!  Variables:

   integer :: i, ios
   real    :: angle, angstrom, distance, dsq, nm, spectral_irradiance

!  Execution:

!  The control file for Stardust_Lines is named lines.inp and
!  line 8 with the viewing angle of elevation above the horizontal
!  now should have a viewing distance in meters following the angle.

   open (lunin, file='lines.inp', iostat=ios)
   if (ios /= 0) then
      write (*, '(a)') '*** Trouble opening lines.inp control file.'
      go to 99
   end if

   do i = 1, line_to_read - 1
      read (lunin, *)  ! Skip the first so many lines
   end do

   read (lunin, *) angle, distance  ! Degrees and meters;  close (lunin)
   distance = distance*100.  ! Meters -> cm
   dsq = distance**2

!  Process integrations.dat one line at a time because we don't know how many:

   open (lunin, file='integrations.dat', iostat=ios)
   if (ios /= 0) then
      write (*, '(a)') '*** Trouble opening integrations.dat control file.'
      go to 99
   end if

   open (lunout, file='point_source.dat', status='unknown')

   do  ! Until EOF
      read (lunin, *, iostat=ios) angstrom, spectral_irradiance
      if (ios < 0) exit

      nm = 0.1*angstrom
      write (lunout, '(f12.4, es12.5)') nm, spectral_irradiance*dsq
   end do

   close (lunin)
   close (lunout)

99 continue

   end program point_source
