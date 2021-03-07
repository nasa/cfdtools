!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      program inertial_and_relative

!     Convert a trajectory time history between inertial and rotational
!     coordinate systems, for Earth.
!
!     Input/output data format is that of 'traj_opt.ascent', but the latitude
!     must be geocentric as in 'traj.plt', and the heading must be between 0=N
!     and 90=E.
!
!     Title
!     Time  Vel, km/s  Gamma, deg   Heading, deg  Alt, km  Lat, geocntr  Long
!      x      xxx.xxx    xxx.xxx     xxx.xxx       xxx.xx     xxx.xx   xxx.xx
!      x      xxx.xxx    xxx.xxx     xxx.xxx       xxx.xx     xxx.xx   xxx.xx
!      :         :          :           :             :          :        :
!
!     10/06/03  D.A.Saunders  Adaptation of Traj formulation; Earth only;
!                             ellipsoid only
!     04/12/07    "     "     Safeguarded some arc sines.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Constants:

      integer, parameter ::
     >   luncrt = 6, lunkbd = 5, lundat = 1

      real, parameter ::
     >   dyn_flat      = 0.00335281,    ! Ellipsoid flattening factor
     >   r_equator     = 6378137.0,     ! Equatorial radius, m
     >   rotation_rate = 7.2895804E-5   ! Earth, radians/s

!     Variables:

      integer
     >   i, ios, n

      real
     >   cos_gam, cos_lat, deg_to_rad, e_2, e_cos_lat_2,
     >   gamma, pi, psi, r_planet, radius, sin_gam, sin_psi,
     >   v_atm, v_i, v_r 

      real, allocatable, dimension (:) ::
     >   time, velocity, fpa, heading, altitude, latitude, longitude

      logical
     >   i_to_r

      character
     >   answer * 1, filename * 64, headers * 80, title * 80

!     Execution:

      ios = 1
      do while (ios /= 0)
         write (luncrt, '(/, a)', advance='no') ' Input file name: '
         read  (lunkbd, '(a)') filename
         open  (lundat, file=filename, status='old', iostat=ios)
      end do

      write (luncrt, '(a)', advance='no') ' i to r [i] or r to i [r]?: '
      read  (lunkbd, '(a)') answer

      i_to_r = answer == 'i' .or. answer == 'I'

!     Count the lines:

      read (lundat, '(a)') ! Skip title
      read (lundat, '(a)') ! Skip headings

      n = 0
      do ! Until EOF
         read (lundat, '(a)', iostat=ios)
         if (ios /= 0) exit
         n = n + 1
      end do

      allocate (time(n), velocity(n), fpa(n), heading(n),
     >          altitude(n), latitude(n), longitude(n))

      rewind (lundat)
      read (lundat, '(a)') title
      read (lundat, '(a)') headers

      do i = 1, n
         read (lundat, *, iostat=ios)
     >      time(i), velocity(i), fpa(i), heading(i),
     >      altitude(i), latitude(i), longitude(i)
         if (ios /= 0) then
            write (luncrt, '(/, a, i5)')
     >         ' Trouble reading input file.  Line #: ', i
            go to 999
         end if
      end do

      close (lundat)

      write (luncrt, '(/, a)', advance='no') ' Output file name: '
      read  (lunkbd, '(a)') filename
      open  (lundat, file=filename, status='unknown', iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(/, 2a)') ' Unable to open ', filename
         go to 999
      end if

!     Transform the data:

      pi         = 2. * asin (1.)
      deg_to_rad = pi / 180.
      e_2        = (2.0 - dyn_flat) * dyn_flat

      if (i_to_r) then ! Inertial to rotational

         do i = 1, n
            sin_gam     = sin (fpa(i)       * deg_to_rad)
            cos_gam     = cos (fpa(i)       * deg_to_rad)
            cos_lat     = cos (latitude(i)  * deg_to_rad)
            e_cos_lat_2 = e_2 * cos_lat * cos_lat

            if (abs(1.0 - e_cos_lat_2) < 1.e-7) then
               r_planet = r_equator
            else
               r_planet = r_equator * sqrt((1.0 - e_2) /
     >                                     (1.0 - e_cos_lat_2))
            end if

            radius  = r_planet + altitude(i) * 1000.
            v_atm   = rotation_rate * radius * cos_lat

            v_i     = velocity(i) * 1000.
            v_r     = v_i * (v_i - 2. * v_atm * cos_gam * 
     >                       sin (heading(i) * deg_to_rad)) +
     >                v_atm ** 2 !heading = 90 is due East
            v_r     = sqrt (v_r)

            gamma   = min (1., max (-1., v_i * sin_gam / v_r))
            gamma   = asin (gamma)
            sin_psi = (v_i ** 2 - v_atm ** 2 - v_r ** 2) /
     >                (2. * v_atm * v_r * cos (gamma))
            sin_psi = min (1., max (-1., sin_psi))
            psi     = asin (sin_psi) / deg_to_rad

            if (heading(i) > 90.) then
               psi =  180. - psi
            else if (heading(i) < -90.) then
               psi = -180. - psi
            end if

            heading(i)  = psi ! sin_psi = 0. = due East
            fpa(i)      = gamma / deg_to_rad
            velocity(i) = v_r * 0.001
         end do

      else ! Rotating to inertial

         do i = 1, n
            sin_gam     = sin (fpa(i)       * deg_to_rad)
            cos_gam     = cos (fpa(i)       * deg_to_rad)
            cos_lat     = cos (latitude(i)  * deg_to_rad)
            e_cos_lat_2 = e_2 * cos_lat * cos_lat

            if (abs(1.0 - e_cos_lat_2) < 1.e-7) then
               r_planet = r_equator
            else
               r_planet = r_equator * sqrt((1.0 - e_2) /
     >                                     (1.0 - e_cos_lat_2))
            end if

            radius  = r_planet + altitude(i) * 1000.
            v_atm   = rotation_rate * radius * cos_lat

            v_r     = velocity(i) * 1000.
            v_i     = v_r * (v_r + 2. * v_atm * cos_gam *
     >                       sin (heading(i) * deg_to_rad)) +
     >                v_atm ** 2 !heading = 90 is due East
            v_i     = sqrt (v_i)

            gamma   = min (1., max (-1., v_r * sin_gam / v_i))
            gamma   = asin (gamma)
            sin_psi = (v_i ** 2 + v_atm ** 2 - v_r ** 2) /
     >                (2. * v_atm * v_i * cos (gamma))
            sin_psi = min (1., max (-1., sin_psi))
            psi     = asin (sin_psi) / deg_to_rad

            if (heading(i) > 90.) then
               psi =  180. - psi
            else if (heading(i) < -90.) then
               psi = -180. - psi
            end if

            heading(i)  = psi ! sin_psi = 0. = due East
            fpa(i)      = gamma / deg_to_rad
            velocity(i) = v_i * 0.001
         end do

      end if

      write (lundat, '(a)')
     >   title(1:len_trim(title)), headers(1:len_trim(headers))

      write (lundat,
     >   '(f8.2, f10.6, 1p, e16.7, 0p, f11.5, f10.5, 2f11.5)')
     >   (time(i), velocity(i), fpa(i), heading(i),
     >    altitude(i), latitude(i), longitude(i), i = 1, n)

  999 close (lundat)

      end program inertial_and_relative
