!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program running_statistics
!
!  This utility applies recursive calculations of mean and standard deviation
!  to body point time histories from DPLR, with a view to incorporating output
!  of running averages, etc., directly.
!
!  Surface body point dataset format:
!
!  Requested Point:   x    y    z
!    Nearest Match:   X    Y    Z
!  Variables = iter  time  dt  p  T  qw  tau_x  tau_y  tau_z  tau
!    nnnnn  tmtmtm  dtdtdt  pppp  TTTT  qwqwqw  txtxtx  tytyty  tztztz  tttt
!      :       :      :      :     :      :       :       :       :      :
!
!  The output format is similar, with lines 1 and 2 commented out for Tecplot
!  usage, and variable names also adjusted for Tecplot sub- & superscripting.
!
!  03/08/12  Initial implementation.
!  04/30/14  Output line 3 should be "variables=iter,t,..." line for Tecplot.
!  05/07/14  Ryan McDaniel asked for DPLR-specific variable names output.
!
!  Author:   David Saunders, ERC, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   integer, parameter :: lundat = 1, lunout = 2, lunkbd = 5, luncrt = 6
   integer, parameter :: nskip = 3   ! # dataset header lines to ignore
   integer, parameter :: numf = 9    ! Itn. # + numf reals expected each line
   real,    parameter :: zero = 0.

   integer   :: i, ios, j, n
   real      :: mean(numf), variance(numf), t1, tprev, fprev(numf)
   character :: buffer*256, filename*80

   integer, allocatable, dimension (:)   :: iter         ! Itn. read w/ raw data
   real,    allocatable, dimension (:,:) :: f            ! Function data

!  Execution:

   write (luncrt, '(a)', advance='no') 'Input time history file name: '
   read  (lunkbd, *) filename
   open  (lundat, file=filename, status='old')

   write (luncrt, '(a)', advance='no') 'Output statistics file name:  '
   read  (lunkbd, *) filename
   open  (lunout, file=filename, status='unknown')

   do n = 1, nskip
      read (lundat, '(a)') buffer
      if (n < nskip) write (lunout, '(2a)') '# ', trim (buffer)
   end do
   write (lunout, '(13a)') &
     'Variables =  iter  time  dt ', &
     'p<sub>mean</sub>(N/m<sup>2</sup>) T<sub>mean</sub>(K) ', &
     'qw<sub>mean</sub>(W/m<sup>2</sup>) ', &
     'taux<sub>mean</sub>(N/m<sup>2</sup>) ', &
     'tauy<sub>mean</sub>(N/m<sup>2</sup>) ', &
     'tauz<sub>mean</sub>(N/m<sup>2</sup>) ', &
     '|tau|<sub>mean</sub> ', &
     'p<sub>sigma</sub>(N/m<sup>2</sup>) T<sub>sigma</sub>(K) ', &
     'qw<sub>sigma</sub>(W/m<sup>2</sup>) ', &
     'taux<sub>sigma</sub>(N/m<sup>2</sup>) ', &
     'tauy<sub>sigma</sub>(N/m<sup>2</sup>) ', &
     'tauz<sub>sigma</sub>(N/m<sup>2</sup>) ', &
     '|tau|<sub>sigma</sub> '

!!!  'variables=iter,t,dt,m1,m2,m3,m4,m5,m6,m7,s1,s2,s3,s4,s5,s6,s7'

   n = 1
   do ! Until EOF  ! Count the data lines
      read (lundat, *, iostat=ios)
      if (ios < 0) exit
      n = n + 1
   end do
   n = n - 1
   write (luncrt, '(/, a, i6)') '# data points read:', n

   allocate (iter(n), f(numf,n))

   rewind (lundat)
   do i = 1, nskip
      read (lundat, '(a)')
   end do

   read (lundat, *) (iter(i), f(:,i), i = 1, n)
   close (lundat)

   t1 = f(1,1)  ! Initial time for all functions 3:numf

   tprev = zero;  fprev(:) = zero  ! Avoids uninitialized args. when i = 1

   do i = 1, n
      do j = 3, numf
         call trapezoid_recursion (i, t1, tprev, fprev(j), f(1,i), f(j,i), &
                                   mean(j), variance(j))
         fprev(j) = f(j,i)
         variance(j) = max (variance(j), zero)
      end do
      tprev = f(1,i)
      write (lunout, '(i8, 1p, 2e14.6, 2x, 7e14.6, 2x, 7e13.5)') &
         iter(i), f(1:2,i), mean(3:numf), sqrt (variance(3:numf))
   end do

   deallocate (iter, f)

   end program running_statistics
