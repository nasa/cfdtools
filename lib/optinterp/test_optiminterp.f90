!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program test_optiminterp

! Example of using the n-dimensional optimal interpolation Fortran module.
! Copyright (C) 2005 Alexander Barth <abarth@marine.usf.edu>
!  
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!  
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!  
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, 
! MA 02110-1301, USA.
!  
! Author: Alexander Barth <abarth@marine.usf.edu>
!
! 11/14/06   D.Saunders  Arranged to vary the problem parameters;
!                        extended to multiple functions per observation.
! 11/15/06    "     "    3-D variant with unnormalized data.
!                        Correlation lengths have to be normalized too!
!                        So do the observation variances.
! 02/09/07    "     "    Correlation function choice requires a new argument.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use optimal_interpolation

   implicit none

!  Constants:

   integer, parameter :: &
      m  = 30,           &  ! # nearest observations influencing interpolations
      n  = 100,          &  ! An n x n 2D grid is being interpolated to
      nd = 3,            &  ! # coordinate dimensions
      nf = 2,            &  ! # functions at each observation
      on = 200              ! # observations

   real, parameter :: &
      correlation_length   = 0.1,  &  ! Needs explanation
      observation_variance = 0.01, & 
      f_offset             = 0.,   &  ! Check automated normalization
      x_offset             = 6.,   &
      y_offset             = -8.,  &
      z_offset             = 0.5

   character, parameter :: &
      correlation_function * 1 = 'G'

!  We use a simple random generator instead of Fortran 90's random_number
!  in the hope that the results will be platform-independent:

   integer(8), parameter :: &
      A = 1664525_8, B = 1013904223_8, Mo = 4294967296_8
   integer(8) :: next = 0_8

!  Variables:

   integer :: &
      i,j,l

   real :: &
      fl,           &
      x(nd,on),     &       ! Coordinates of the observations
      f(nf,on),     &       ! 1 or more functions per observation
      var(on),      &       ! Variances of the observations
      cor_len(nd),  &       ! Correlation lengths in each dimension
      xi(nd,n*n),   &       ! Coordinates of the grid being interpolated to
      fi(nf,n*n),   &       ! Interpolated field
      vari(n*n)             ! Interpolated error variance

   character :: &
      fun(nd)*1

!  Execution:

   fun(:) = correlation_function

!  Create a regular 2D n x n grid:

   l = 1
   do j = 1, n
      do i = 1, n
         xi(1,l) = (i-1.)/(n-1.)
         xi(2,l) = (j-1.)/(n-1.)
         l = l + 1
      end do
   end do

   cor_len(:) = correlation_length
   var(:)     = observation_variance

!  Random locations of the observations:

   do j = 1, on
      do i = 1, nd
         next   = mod (A*next + B, Mo)     ! simple random generator
         x(i,j) = real (next, 8) / real (Mo, 8)
      end do
   end do

!  Analytic underlying function to interpolate:

   do l = 1, on
      fl = sin (x(1,l)*6.) * cos (x(2,l)*6.)
      f(1,l) =       fl + f_offset
      x(1,l) = x(1,l)   + x_offset
      x(2,l) = x(2,l)   + y_offset
      x(3,l) = 0.1 * fl + z_offset         ! Mildly wavy surface
   end do

!  Offset the observation (x,y)s and functions to see if normalization
!  appears necessary (as certain documentation suggests):

   if (nf > 1) f(2,:) = -f(1,:)

!  Offset the target (x,y,z)s to match the offset data:

   do l = 1, n * n
      xi(3,l) = 0.1 * sin (xi(1,l)*6.) * cos (xi(2,l)*6.) + z_offset
      xi(1,l) = xi(1,l) + x_offset
      xi(2,l) = xi(2,l) + y_offset
   end do

!  Perform the interpolations:

   call optiminterp (x, f, var, cor_len, fun, m, xi, fi, vari)

!  Display some control values and save results for plotting:

   write (6,'(a,f24.16)') 'Expected fi(1), m = 30:  ', 2.2205936104348591E-02
   write (6,'(a,i3,a,f24.16)') 'Obtained fi(1), m =', m, ':  ',fi(1,1)
   write (6,'(a,f24.16)') 'Obtained vari(1):        ', vari(1)

   open  (1, file='interpolations.g', status='unknown')
   write (1, '(i1)') 2
   write (1, '(3i5)') on, 1, 1, n, n, 1
   write (1, '(1p, 6e19.11)') x(1,:), x(2,:), x(3,:)
   write (1, '(1p, 6e19.11)') xi(1,:), xi(2,:), xi(3,:)
   close (1)

   open  (1, file='interpolations.f', status='unknown')
   write (1, '(i1)') 2
   write (1, '(4i5)') on, 1, 1, nf, n, n, 1, nf
   write (1, '(1p, 6e19.11)')  (f(i,:), i = 1, nf)
   write (1, '(1p, 6e19.11)') (fi(i,:), i = 1, nf)
   close (1)

   open  (1, file='variances.f', status='unknown')
   write (1, '(i1)') 2
   write (1, '(4i5)') on, 1, 1, 1, n, n, 1, 1
   write (1, '(1p, 6e19.11)') var
   write (1, '(1p, 6e19.11)') vari
   close (1)

   end program test_optiminterp
