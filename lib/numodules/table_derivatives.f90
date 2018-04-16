!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calculate 1st and 2nd derivatives w.r.t. each coordinate for nf functions
!  defined on a rectangular table with ndim dimensions, using 3-point finite
!  differencing.  A curvature-like quantity is also returned at each point.
!
!  The present approach takes advantage of multiple subscripts, but it's hard to
!  generalize.  Therefore, a different routine is provided for each number of
!  dimensions desired.  Add more dimensions as necessary.
!
!  The language supports implied passing of planes (and higher-rank sub-arrays)
!  as arguments.
!
!  All argument arrays are assumed to be sized exactly.  The coordinates and
!  functions at one table entry are contiguous (first index).
!
!  03/14/07  Initial adaptation of table_coordinates.f90.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine derivs_1d (nx, nf, x, f, f1, f2, fk)

!  This is a variant of FD12K for the case of more than one function.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: nx                  ! # data points
   integer, intent (in)  :: nf                  ! # functions at each point
   real,    intent (in)  :: x(nx)               ! Coordinates
   real,    intent (in)  :: f(nf,nx)            ! 1 or more functions per point
   real,    intent (out), dimension (nf,nx) :: f1, f2, fk  ! 1st & 2nd derivs.
                                                           ! and curvature fn.
!  Local variables:

   integer :: n

!  Execution:

   do n = 1, nf

      call fd12k (nx, x, f(n,:), f1(n,:), f2(n,:), fk(n,:))

   end do

   end subroutine derivs_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine derivs_2d (ndim, dim1, dim2, nf, x, f, f1, f2, fk)

   implicit none

!  Arguments:

   integer, intent (in)  :: ndim                ! Dimensionality of the table
   integer, intent (in)  :: dim1, dim2          ! That number of dimensions
   integer, intent (in)  :: nf                  ! # functions at each point
   real,    intent (in)  :: &                   ! Expanded form of coordinates
      x(ndim,dim1,dim2)
   real,    intent (in)  :: &                   ! Expanded form of functions
      f(nf,dim1,dim2)
   real,    intent (out), dimension (nf,dim1,dim2,ndim) :: &  ! Derivatives and
      f1, f2, fk                                              ! curvature-like f

!  Local variables:

   integer :: i, j, n

!  Execution:

!  Derivatives in the i direction:

   do j = 1, dim2
      do n = 1, nf
         call fd12k (dim1, x(1,:,1), f(n,:,j), f1(n,:,j,1), f2(n,:,j,1),       &
                     fk(n,:,j,1))
      end do
   end do

!  j direction:

   do i = 1, dim1
      do n = 1, nf
         call fd12k (dim2, x(2,1,:), f(n,i,:), f1(n,i,:,2), f2(n,i,:,2),       &
                     fk(n,i,:,2))
      end do
   end do

   end subroutine derivs_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine derivs_3d (ndim, dim1, dim2, dim3, nf, x, f, f1, f2, fk)

   implicit none

!  Arguments:

   integer, intent (in)  :: ndim                ! Dimensionality of the table
   integer, intent (in)  :: dim1, dim2, dim3    ! That number of dimensions
   integer, intent (in)  :: nf                  ! # functions at each point
   real,    intent (in)  :: &                   ! Expanded form of coordinates
      x(ndim,dim1,dim2,dim3)
   real,    intent (in)  :: &                   ! Expanded form of functions
      f(nf,dim1,dim2,dim3)
   real,    intent (out), dimension (nf,dim1,dim2,dim3,ndim) :: &
      f1, f2, fk                        ! Derivatives and the curvature-like fn.

!  Local variables:

   integer :: i, j, k, n

!  Execution:

!  i direction:

   do k = 1, dim3
      do j = 1, dim2
         do n = 1, nf
            call fd12k (dim1, x(1,:,1,1), f(n,:,j,k), f1(n,:,j,k,1),           &
                        f2(n,:,j,k,1), fk(n,:,j,k,1))
         end do
      end do
   end do

!  j direction:

   do k = 1, dim3
      do i = 1, dim1
         do n = 1, nf
            call fd12k (dim2, x(2,1,:,1), f(n,i,:,k), f1(n,i,:,k,2),           &
                        f2(n,i,:,k,2), fk(n,i,:,k,2))
         end do
      end do
   end do

!  k direction:

   do j = 1, dim2
      do i = 1, dim1
         do n = 1, nf
            call fd12k (dim3, x(3,1,1,:), f(n,i,j,:), f1(n,i,j,:,3),           &
                        f2(n,i,j,:,3), fk(n,i,j,:,3))
         end do
      end do
   end do

   end subroutine derivs_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine derivs_4d (ndim, dim1, dim2, dim3, dim4, nf, x, f, f1, f2, fk)

   implicit none

!  Arguments:

   integer, intent (in)  :: ndim                ! Dimensionality of the table
   integer, intent (in)  :: dim1, dim2, dim3, dim4  ! That number of dimensions
   integer, intent (in)  :: nf                  ! # functions at each point
   real,    intent (in)  :: &                   ! Expanded form of coordinates
      x(ndim,dim1,dim2,dim3,dim4)
   real,    intent (in)  :: &                   ! Expanded form of functions
      f(nf,dim1,dim2,dim3,dim4)
   real,    intent (out), dimension (nf,dim1,dim2,dim3,dim4,ndim) :: &
      f1, f2, fk                        ! Derivatives and the curvature-like fn.

!  Local variables:

   integer :: i, j, k, l, n

!  Execution:

!  i direction:

   do l = 1, dim4
      do k = 1, dim3
         do j = 1, dim2
            do n = 1, nf
               call fd12k (dim1, x(1,:,1,1,1), f(n,:,j,k,l), f1(n,:,j,k,l,1),  &
                           f2(n,:,j,k,l,1), fk(n,:,j,k,l,1))
            end do
         end do
      end do
   end do

!  j direction:

   do l = 1, dim4
      do k = 1, dim3
         do i = 1, dim1
            do n = 1, nf
               call fd12k (dim2, x(2,1,:,1,1), f(n,i,:,k,l), f1(n,i,:,k,l,2),  &
                           f2(n,i,:,k,l,2), fk(n,i,:,k,l,2))
            end do
         end do
      end do
   end do

!  k direction:

   do l = 1, dim4
      do j = 1, dim2
         do i = 1, dim1
            do n = 1, nf
               call fd12k (dim3, x(3,1,1,:,1), f(n,i,j,:,l), f1(n,i,j,:,l,3),  &
                           f2(n,i,j,:,l,3), fk(n,i,j,:,l,3))
            end do
         end do
      end do
   end do

!  l direction:

   do k = 1, dim3
      do j = 1, dim2
         do i = 1, dim1
            do n = 1, nf
               call fd12k (dim4, x(4,1,1,1,:), f(n,i,j,k,:), f1(n,i,j,k,:,4),  &
                           f2(n,i,j,k,:,4), fk(n,i,j,k,:,4))
            end do
         end do
      end do
   end do

   end subroutine derivs_4d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine derivs_5d (ndim, dim1, dim2, dim3, dim4, dim5, nf, x, f, f1, f2, &
                         fk)
   implicit none

!  Arguments:

   integer, intent (in)  :: ndim                ! Dimensionality of the table
   integer, intent (in)  :: dim1, dim2, dim3, & ! That number of dimensions
                            dim4, dim5
   integer, intent (in)  :: nf                  ! # functions at each point
   real,    intent (in)  :: &                   ! Expanded form of coordinates
      x(ndim,dim1,dim2,dim3,dim4,dim5)
   real,    intent (in)  :: &                   ! Expanded form of functions
      f(nf,dim1,dim2,dim3,dim4,dim5)
   real,    intent (out), dimension (nf,dim1,dim2,dim3,dim4,dim5,ndim) :: &
      f1, f2, fk                        ! Derivatives and the curvature-like fn.


!  Local variables:

   integer :: i, j, k, l, m, n

!  Execution:

!  i direction:

   do m = 1, dim5
      do l = 1, dim4
         do k = 1, dim3
            do j = 1, dim2
               do n = 1, nf
                  call fd12k (dim1, x(1,:,1,1,1,1), f(n,:,j,k,l,m),            &
                              f1(n,:,j,k,l,m,1), f2(n,:,j,k,l,m,1),            &
                              fk(n,:,j,k,l,m,1))
               end do
            end do
         end do
      end do
   end do

!  j direction:

   do m = 1, dim5
      do l = 1, dim4
         do k = 1, dim3
            do i = 1, dim1
               do n = 1, nf
                  call fd12k (dim2, x(2,1,:,1,1,1), f(n,i,:,k,l,m),            &
                              f1(n,i,:,k,l,m,2), f2(n,i,:,k,l,m,2),            &
                              fk(n,i,:,k,l,m,2))
               end do
            end do
         end do
      end do
   end do

!  k direction:

   do m = 1, dim5
      do l = 1, dim4
         do j = 1, dim2
            do i = 1, dim1
               do n = 1, nf
                  call fd12k (dim3, x(3,1,1,:,1,1), f(n,i,j,:,l,m),            &
                              f1(n,i,j,:,l,m,3), f2(n,i,j,:,l,m,3),            &
                              fk(n,i,j,:,l,m,3))
               end do
            end do
         end do
      end do
   end do

!  l direction:

   do m = 1, dim5
      do k = 1, dim3
         do j = 1, dim2
            do i = 1, dim1
               do n = 1, nf
                  call fd12k (dim4, x(4,1,1,1,:,1), f(n,i,j,k,:,m),            &
                              f1(n,i,j,k,:,m,4), f2(n,i,j,k,:,m,4),            &
                              fk(n,i,j,k,:,m,4))
               end do
            end do
         end do
      end do
   end do

!  m direction:

   do l = 1, dim4
      do k = 1, dim3
         do j = 1, dim2
            do i = 1, dim1
               do n = 1, nf
                  call fd12k (dim5, x(5,1,1,1,1,:), f(n,i,j,k,l,:),            &
                              f1(n,i,j,k,l,:,5), f2(n,i,j,k,l,:,5),            &
                              fk(n,i,j,k,l,:,5))
               end do
            end do
         end do
      end do
   end do

   end subroutine derivs_5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! This approach requires 8 subscripts which exceeds the Fortran 90 limit.
!! Deal with it later if we ever have to.

!! subroutine derivs_6d (ndim, dim1, dim2, dim3, dim4, dim5, dim6, nf, x, f,   &
!!                       f1, f2, fk)
!! implicit none

!  Arguments:

!! integer, intent (in)  :: ndim                ! Dimensionality of the table
!! integer, intent (in)  :: dim1, dim2, dim3, & ! That number of dimensions
!!                          dim4, dim5, dim6
!! integer, intent (in)  :: nf                  ! # functions at each point
!! real,    intent (in)  :: &                   ! Expanded form of coordinates
!!    x(ndim,dim1,dim2,dim3,dim4,dim5,dim6)
!! real,    intent (in)  :: &                   ! Expanded form of functions
!!    f(nf,dim1,dim2,dim3,dim4,dim5,dim6)
!! real,    intent (out), dimension (nf,dim1,dim2,dim3,dim4,dim5,dim6,ndim) :: &
!!    f1, f2, fk                        ! Derivatives and the curvature-like fn.


!  Local variables:

!! integer :: i, j, k, l, m, n, if

!  Execution:

!  i direction:

!! do n = 1, dim6
!!    do m = 1, dim5
!!       do l = 1, dim4
!!          do k = 1, dim3
!!             do j = 1, dim2
!!                do if = 1, nf
!!                   call fd12k (dim1, x(1,:,1,1,1,1,1), f(if,:,j,k,l,m,n),    &
!!                               f1(if,:,j,k,l,m,n,1), f2(if,:,j,k,l,m,n,1),   &
!!                               fk(if,:,j,k,l,m,n,1))
!!                end do
!!             end do
!!          end do
!!       end do
!!    end do
!! end do

!  j direction:

!! do n = 1, dim6
!!    do m = 1, dim5
!!       do l = 1, dim4
!!          do k = 1, dim3
!!             do i = 1, dim1
!!                do if = 1, nf
!!                   call fd12k (dim2, x(2,1,:,1,1,1,1), f(if,i,:,k,l,m,n),    &
!!                               f1(if,i,:,k,l,m,n,2), f2(if,i,:,k,l,m,n,2),   &
!!                               fk(if,i,:,k,l,m,n,2))
!!                end do
!!             end do
!!          end do
!!       end do
!!    end do
!! end do

!  k direction:

!! do n = 1, dim6
!!    do m = 1, dim5
!!       do l = 1, dim4
!!          do j = 1, dim2
!!             do i = 1, dim1
!!                do if = 1, nf
!!                   call fd12k (dim3, x(3,1,1,:,1,1,1), f(if,i,j,:,l,m,n),    &
!!                               f1(if,i,j,:,l,m,n,3), f2(if,i,j,:,l,m,n,3),   &
!!                               fk(if,i,j,:,l,m,n,3))
!!                end do
!!             end do
!!          end do
!!       end do
!!    end do
!! end do

!  l direction:

!! do n = 1, dim6
!!    do m = 1, dim5
!!       do k = 1, dim3
!!          do j = 1, dim2
!!             do i = 1, dim1
!!                do if = 1, nf
!!                   call fd12k (dim4, x(4,1,1,1,:,1,1), f(if,i,j,k,:,m,n),    &
!!                               f1(if,i,j,k,:,m,n,4), f2(if,i,j,k,:,m,n,4),   &
!!                               fk(if,i,j,k,:,m,n,4))
!!                end do
!!             end do
!!          end do
!!       end do
!!    end do
!! end do

!  m direction:

!! do n = 1, dim6
!!    do l = 1, dim4
!!       do k = 1, dim3
!!          do j = 1, dim2
!!             do i = 1, dim1
!!                do if = 1, nf
!!                   call fd12k (dim5, x(5,1,1,1,1,:,1), f(if,i,j,k,l,:,n),    &
!!                               f1(if,i,j,k,l,:,n,5), f2(if,i,j,k,l,:,n,5),   &
!!                               fk(if,i,j,k,l,:,n,5))
!!                end do
!!             end do
!!          end do
!!       end do
!!    end do
!! end do

!  n direction:

!! do m = 1, dim5
!!    do l = 1, dim4
!!       do k = 1, dim3
!!          do j = 1, dim2
!!             do i = 1, dim1
!!                do if = 1, nf
!!                   call fd12k (dim6, x(6,1,1,1,1,1,:), f(if,i,j,k,l,m,:),    &
!!                               f1(if,i,j,k,l,m,:,6), f2(if,i,j,k,l,m,:,6),   &
!!                               fk(if,i,j,k,l,m,:,6))
!!                end do
!!             end do
!!          end do
!!       end do
!!    end do
!! end do

!! end subroutine derivs_6d
