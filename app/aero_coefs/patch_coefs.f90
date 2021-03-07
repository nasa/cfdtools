!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine patch_coefs (ni, nj, nf, ip, x, y, z, f, xyzm, Lref, Sref, coefs)
!
!  For one structured surface patch, calculate aerodynamic force and moment
!  coefficients.  Totals over all patches can be summed at the higher level.
!
!  06/01/2018  D.A.Saunders  Lower level routine for program aero_coefs.
!  06/05/2018     "    "     Dinesh found the Sref was missing from the moments.
!  02/28/2019     "    "     Surface normals point outward from a convex
!                            right-handed surface, so they are now negated
!                            in order for the x component of the forces to
!                            be positive as one would normally expect.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Reserach Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: ni, nj                     ! Patch dimensions
   integer, intent (in)  :: nf                         ! For dimensioning f(:,:)
   integer, intent (in)  :: ip                         ! Pressure index in f
   real,    intent (in), dimension (ni,nj) :: x, y, z  ! Patch coordinates
   real,    intent (in)  :: f(nf,ni,nj)                ! Vertex-centered fns.
   real,    intent (in)  :: xyzm(3)                    ! Moment center
   real,    intent (in)  :: Lref, Sref                 ! Reference length, area
   real,    intent (out) :: coefs(6)                   ! Force & moment coefs.
                                                       ! for this patch
!  Local constants:

   real, parameter :: fourth = 1./4., half = 0.5, zero = 0.

!  Local variables:

   integer :: i, j
   real :: area, fcenter, u(3), xyzc(3)

!  Procedures:

   real, external :: area4  ! Quadrilateral area

!  Execution:

   coefs(:) = zero

   do j = 1, nj - 1  ! Work with cell centroids and area
      do i = 1, ni - 1

!        The cell centroid isn't needed for this surface normal utility,
!        which is compatible with the ADT surface searching utility.

         call surface_normal (ni, nj, x, y, z, i, j, half, half, u)

         u(:) = -u(:)  ! We want the unit normal to point into the body

         area = area4 (x(i,j),     y(i,j),     z(i,j),     &
                       x(i+1,j),   y(i+1,j),   z(i+1,j),   &
                       x(i+1,j+1), y(i+1,j+1), z(i+1,j+1), &
                       x(i,j+1),   y(i,j+1),   z(i,j+1))

         xyzc(1) = (x(i,j) + x(i+1,j) + x(i,j+1) + x(i+1,j+1))*fourth
         xyzc(2) = (y(i,j) + y(i+1,j) + y(i,j+1) + y(i+1,j+1))*fourth
         xyzc(3) = (z(i,j) + z(i+1,j) + z(i,j+1) + z(i+1,j+1))*fourth
         fcenter = (f(ip,i,j)+f(ip,i+1,j)+f(ip,i,j+1)+f(ip,i+1,j+1))*fourth

         coefs(1:3) = coefs(1:3) + u(1:3)*fcenter*area
         coefs(4)   = coefs(4) - u(2)*area*(xyzc(3) - xyzm(3)) &
                               + u(3)*area*(xyzc(2) - xyzm(2))
         coefs(5)   = coefs(5) - u(3)*area*(xyzc(1) - xyzm(1)) &
                               + u(1)*area*(xyzc(3) - xyzm(3))
         coefs(6)   = coefs(6) - u(1)*area*(xyzc(2) - xyzm(2)) &
                               + u(2)*area*(xyzc(1) - xyzm(1))
      end do
   end do

   coefs(1:3) = coefs(1:3)/Sref
   coefs(4:6) = coefs(4:6)/(Sref*Lref)

   end subroutine patch_coefs
