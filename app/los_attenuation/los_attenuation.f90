!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program los_attenuation
!
!  Description:
!
!     For the input line of sight data (los.g and los.f in Plot3D format) and
!     mission-specific radio frequency or frequencies, calculate a measure of
!     the expected radio attenuation.  The formulation used initially here is
!     that employed by program Radio_Attenuation, c. 2016, which exploited the
!     existence of a substantial aerothermal database for the Orion spacecraft.
!     It generated a secondary database for all signal directions from multiple
!     antennas.
!
!     This version can also read los.dat (point-order Tecplot file), and it
!     allows for variable function numbers for p, T, eND.
!
!  Reference:
!
!     Mitigation of Communication Blackout During Re-entry Using Static Magnetic
!     Field, Neha Mehra, Rajesh K. Singh, and Subhash C. Bera,
!     Progress in Electromagnetic Research B, Vol. 63, pp. 161-172, 2015
!
!  Formulation Correction From Radio_Attenuation:
!
!     03/17/2016    Rick Thompson determined that plasma frequency
!                   omegap is an angular frequency in radians/sec.
!                   Therefore, the radio frequencies should be in
!                   the same units, so multiply Hz by 2 pi.
!
!  Possible Issue:
!
!     Is this formulation Earth-specific?
!     This version includes an electron number density-only method that should
!     be applicable to any atmosphere, and is invoked via atmosphere = 'any'.
!
!  Namelist Controls (attenuation.inp; omit trailing comments)
!
!     $attenuation
!     los_form                                   (1|2 = los.g + los.f | los.dat)
!     time = 50.5                                            (BET time, seconds)
!     bp = 1.23 2.34 3.45 30. 0. (x,y,z,cone,clock coords. of radio transmitter)
!     receiver = 'MRO'                (receiving asset; 'DTO' = direct to Earth)
!     atmosphere = 'Earth'    (any planet or 'Titan' or 'any', case-insensitive)
!     frequency  = 400.                                                    (MHz)
!     ip = 1                                    (function numbers for p, T, eND)
!     iT = 2
!     ie = 3
!     $end
!
!  LOS data:
!
!     The initial formulation expects p, T, electron number density as functions
!     1, 2, 3 in los.f, in SI units.  Any further functions will be ignored.
!     [This has now been generalized via ip, iT, ie which default to 1, 2, 3.]
!     Other forms of the attenuation calculations may require a different set of
!     functions, including eND alone.  Since a given atmosphere may be treated
!     with more than one formulation, use of a planet name to define the formu-
!     lation is problematic.
!        Suggestion:
!     Introduce names 'Earth1', 'Earth2', ... if the need arises.
!
!     Since the plottable outputs contain p and T, they should be included in
!     the LOS dataset whether the attenuation calculation uses them or now.
!
!  Output:
!
!     o  Teclotable file of LOS data and the integrand, integrands.dat:
!
!        # time (s): 50.5
!        # bp x/y/z (m): 1.23 2.34 3.45 30. 0.
!        # receiver: 'MRO'
!        # atmosphere: 'Mars'
!        # frequency (MHz): 400.
!        # attenuation 1(dB): 50.
!        # eNDmax  (per m^3): 1.234e20 
!        # eNDcrit (per m^3): 1.012e20 
!        variables=arc[m],p[Pa],T[K],eND[m^-3],integrand[?]
!        0.000000E+00 4.123456E+01 3.456789E+03 1.234567E+18
!        1.123456E-05  :            :            :
!         :            :            :            :
!
!  History:
!
!     Jan-Jul 2016  D.A.Saunders  Radio_Attenuation implementation for Orion.
!     08-07-2022      "     "     Design of single-LOS adaptation.
!     08-10-2022      "     "     Allow for atmosphere specifics; initial code.
!     08-12-2022      "     "     Realized that a planet name may not be
!                                 sufficient to define what is expected in the
!                                 los.f file.  See discussion above.
!     08-17-2022      "     "     Critical density was set in the wrong place.
!                                 Allow for reading the line of sight data as a
!                                 point-order Tecplot file, los.dat.  Allow for
!                                 variable function numbers for p, T, eND.
!     08-19-2022      "     "     Introduced the electrons-only method shown in
!                                 slides by summer intern Trevor Hedges, and
!                                 invoked with atmosphere = 'any'.  This meant
!                                 moving a bunch of constants up to the top
!                                 level so that variations of the attenuation
!                                 formulation can share them.
!     08-xx-2022      "     "     Initial testing completed.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_header_structure       ! Contained within Tecplot_io.f90
   use grid_block_structure        !     "       "       "      "
   use tecplot_io_module           ! Structured Tecplot I/O utilities
   use xyzq_io_module              ! PLOT3D-type file I/O

   implicit none

!  Constants:

   integer, parameter:: &
      lung    = 1,      &          ! los.g (Plot3D format)
      lunf    = 2,      &          ! los.f or los.dat
      lunctl  = 3,      &          ! Namelist attenuation controls
      lunplot = 4,      &          ! integrands.dat (Tecplotable columns)
      luncrt  = 6                  ! Screen diagnostics

   logical, parameter :: &
       false  = .false., &
       true   = .true.

!  More constants moved from the original p, T, eND attenuation routine so that
!  further variations can reuse them:

   real, parameter :: &
      half    = 0.5,                &
      one     = 1.,                 &
      zero    = 0.,                 &
      c       = 2.99792458e+8,      &  ! Speed of light, m/s
      cr      = one/c,              &
      echarge = 1.60217662e-19,     &  ! Electron charge, coulombs
      emass   = 9.109e-31,          &  ! Electron mass, kg
      eps0    = 8.854188e-12,       &  ! Permittivity of free space, farads/m
      eterm   = echarge**2 / (eps0*emass), &
      factor  = 8.980,              &  ! Plasma freq. = 8980 sqrt (eND) Hz
      twopi   = 6.2831853071795865     ! if eND is # /cm^3  [Reference?]

!  Variables:

   integer :: &
      ie, ios, ip, iT, k, los_form, nblocks, nf, nk

   real :: &
      atten, bp(5), eNDcritical, eNDmax, frequency, time

   real :: &  ! Moved here for sharing
      alphap, freq, k0, Kisq, Kr, nusq, Ne, omegasq, omegapsq, p, T, term

   real, allocatable :: &
      f(:), integrands(:,:)

   logical :: &
      cell_centered

   character (8) :: &
      atmosphere, receiver

   character (80) :: &
      descriptor(7)

   type (grid_header) :: &
      los_header
   type (grid_type), dimension(:), pointer :: &
      los

    namelist /attenuation/ &
       los_form, atmosphere, bp, frequency, receiver, time, ip, iT, ie

!  Execution:

!  Read the control inputs and the line of sight data:

   call read_inputs ()
   if (ios /= 0) go to 99

!  Compute the integrated attenuation:

   call compute_attenuation ()

!  Save plottable integrand data:

   call save_results ()

99 continue

   contains  ! Local procedures for program los_attenuation

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_inputs ()  ! Namelist controls & los.g, los.f or los.dat
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      open  (lunctl, file='attenuation.inp', status='old')

!     Control defaults:

      los_form   = 1
      atmosphere = 'Mars'
      bp         = [0., 0., 0., 0., 0.]
      frequency  = 400.
      receiver   = 'MRO'
      time       = 999.
      ip         = 1  ! Function numbers
      iT         = 2
      ie         = 3

      read  (lunctl, nml=attenuation)
      close (lunctl)
      write (luncrt, nml=attenuation)

      select case (los_form)

         case (1)  ! Plot3D formatted

            open (lung, file='los.g', status='old', iostat=ios)
            if (ios /= 0) then
               write (luncrt, '(a)') ' *** Unable to open los.g.'
               go to 99
            end if

            open (lunf, file='los.f', status='old', iostat=ios)
            if (ios /= 0) then
               write (luncrt, '(a)') ' *** Unable to open los.f.'
               go to 99
            end if
 
            call xyzq_read (lung, lunf, true, nblocks, nf, cell_centered, los, &
                            ios)
            if (ios /= 0) then
               write (luncrt, '(a)') ' *** Trouble reading LOS data.'
!!!            go to 99
            end if

         case (2)  ! Tecplot ASCII

            los_header%filename  = 'los.dat'
            los_header%formatted = true
            los_header%ndim      = 3

            call Tecplot_read (lunf, los_header, los, ios)
            if (ios /= 0) then
               write (luncrt, '(a)') ' *** Trouble reading LOS data.'
!!!            go to 99
            end if

         case default

            write (luncrt, '(/, (a, i9))') &
               ' Bad target format #: ', los_form, &
               ' 1 = Plot3D; 2 = Tecplot file.'
!!!         go to 99

      end select

 99   return

      end subroutine read_inputs

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine compute_attenuation ()
!
!     Reference:  See program header.  Rick Thompson determined that the plasma
!     resonant frequency, omegap, should be in radians/sec, and hence so should
!     the radio frequencies.  This original formulation for Earth has had to be
!     remodularized so that other formulations can share common constants..
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants (mostly moved to the top level now):

!!!   real, parameter :: &
!!!      one     = 1.,                 &
!!!      c       = 2.99792458e+8,      &  ! Speed of light, m/s
!!!      cr      = one/c,              &
!!!      echarge = 1.60217662e-19,     &  ! Electron charge, coulombs
!!!      emass   = 9.109e-31,          &  ! Electron mass, kg
!!!      eps0    = 8.854188e-12,       &  ! Permittivity of free space, farads/m
!!!      eterm   = echarge**2 / (eps0*emass), &
!!!      factor  = 8.980,              &  ! Plasma freq. = 8980 sqrt (eND) Hz
!!!      twopi   = 6.2831853071795865, &  ! if eND is # /cm^3  [Reference?]
!!!      half    = 0.5,                &
!!!      zero    = 0.

      character (1), parameter :: &
         method = 'M'                     ! Monotonic local spline quadra.

!     Local variables:

!!!   real :: alphap, freq, k0, Kisq, Kr, nusq, Ne, omegasq, omegapsq, &
!!!           p, T, term, total
      real :: total
      real, allocatable :: s(:), x(:), y(:), z(:)

      character (8) :: atmos

!     Execution:

!     Unit adjustments:

      freq = frequency*1.e+6            ! MHz -> Hz
      eNDcritical = (freq / factor)**2  ! m^-3; reference?
      freq = freq*twopi                 ! per sec -> radians/sec

      nk = los(1)%nk  ! Los is a (1,1,nk) block or a single zone

!     Set up the function to be integrated w.r.t. arc length:

      allocate (f(nk))

      atmos = atmosphere;  call upcase (atmos)

!     Potential awkwardnesses will have to be resolved as they arise.
!     One suggestion:  Resort to 'Earth2', say, to distinguished variations.

      select case (atmos)

         case ('ANY')
            call electrons_only ()
         case ('VENUS  ')
         case ('EARTH  ')
            call pTeND_formulation ()
         case ('MARS   ')
         case ('JUPITER')
         case ('SATURN ')
         case ('TITAN  ')
         case ('URANUS ')
         case ('NEPTUNE')
         case default

      end select

!     Calculate arc lengths:

      allocate (x(nk), y(nk), z(nk), s(nk))

      do k = 1, nk
         x(k) = los(1)%x(1,1,k)
         y(k) = los(1)%y(1,1,k)
         z(k) = los(1)%z(1,1,k)
      end do

      call chords3d (nk, x, y, z, false, total, s) 

!     Local cubic spline quadrature of f vs. s:

      call lcsquad (nk, s, f, zero, total, method, atten)

!     Plottable outputs:

      allocate (integrands(nk,0:4))

      integrands(:,0) = s(:)               ! arc length
      integrands(:,1) = los(1)%q(ip,1,1,:)
      integrands(:,2) = los(1)%q(iT,1,1,:)
      integrands(:,3) = los(1)%q(ie,1,1,:)
      integrands(:,4) = f(:)               ! integrand

      eNDmax = maxval (los(1)%q(ie,1,1,:))

      end subroutine compute_attenuation

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine electrons_only ()  ! Spencer, 1964
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real :: const, loge10term

!!!   log10(e^fun) = fun/ln(10)

      loge10term = 20./log (10.)
!!!   k0 = twopi*freq*cr
      k0 = freq*cr  ! 2 pi has already been applied to frequency
      const = loge10term*k0

      do k = 1, nk
         Ne   = los(1)%q(ie,1,1,k)
!!!      f(k) = const * sqrt (max (one, Ne*eterm/(twopi*freq)**2) - one)
         f(k) = const * sqrt (max (one, Ne*eterm/freq**2) - one)
      end do

      end subroutine electrons_only

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine pTeND_formulation ()
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      k0      = freq*cr
      omegasq = freq**2
      do k = 1, nk
         p        = los(1)%q(ip,1,1,k)
         T        = los(1)%q(iT,1,1,k)
         Ne       = los(1)%q(ie,1,1,k)
         omegapsq = Ne*eterm
         nusq     = (5.738e+7 * p)**2 / T
         term     = omegapsq / (omegasq + nusq)
         Kr       = one - term
         Kisq     = term**2 * (nusq/omegasq)
         alphap   = k0 * sqrt (half*(sqrt (Kr**2 + Kisq) - Kr))
         f(k)     = alphap
      end do

      end subroutine pTeND_formulation

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine save_results ()  ! Tecplotable output
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Execution:

      open (lunplot, file='integrands.dat', status='unknown')

      write (descriptor(1), '(a, f8.3)')   '# time (s):            ', time
      write (descriptor(2), '(a, 5f10.4)') '# bp x y z cone clock: ', bp(:)
      write (descriptor(3), '(a, a)')      '# atmosphere:          ', atmosphere
      write (descriptor(4), '(a, f9.2)')   '# frequency (MHz):     ', frequency
      write (descriptor(5), '(a, f12.3)')  '# attenuation (dB):    ', atten
      write (descriptor(6), '(a, es15.6)') '# eNDmax (per m^3):    ', eNDmax
      write (descriptor(7), '(a, es15.6)') '# eNDcrit (per m^3):   ',eNDcritical
      write (lunplot, '(a)') &
                         'variables=arc[m],p[Pa],T[K],eND[m^-3],integrand[dB/m]'
      write (lunplot, '(a)') (trim (descriptor(k)), k = 1, 7)
      write (lunplot, '(5es22.14)') (integrands(k,:), k = 1, nk)
      close (lunplot)

      end subroutine save_results

   end program los_attenuation
