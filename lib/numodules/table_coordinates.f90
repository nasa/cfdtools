!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Expand the coordinates of a rectangular table from the compact form (a vector
!  of coordinates for each dimension, packed into one vector).
!
!  If anyone can see how to do this without itemizing the cases, please let me
!  know!  It requires a variable number of nested loops.  Recursion may be the
!  answer, but the weak way out is taken here.
!
!  The present approach takes advantage of multiple subscripts, but it's hard to
!  generalize.  Therefore, a different routine is provided for each number of
!  dimensions desired.  Add more dimensions as necessary.
!
!  All argument arrays are assumed to be sized exactly.  The coordinates at one
!  table entry are contiguous (first index).
!
!  01/22/07  Initial implementation.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine expand_2d (ndim, dim1, dim2, packed_coords, expanded_coords)

   implicit none

!  Arguments:

   integer, intent (in)  :: ndim                ! Dimensionality of the table
   integer, intent (in)  :: dim1, dim2          ! That number of dimensions
   real,    intent (in)  :: &                   ! Compact form of coordinates
      packed_coords(dim1+dim2)
   real,    intent (out) :: &                   ! Expanded form
      expanded_coords(ndim,dim1,dim2)

!  Local variables:

   integer :: i, j
   integer ::    jo  ! Offset

!  Execution:

   jo = dim1

   do j = 1, dim2
      do i = 1, dim1
         expanded_coords(1,i,j) = packed_coords(i)
         expanded_coords(2,i,j) = packed_coords(jo+j)
      end do
   end do

   end subroutine expand_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine expand_3d (ndim, dim1, dim2, dim3, packed_coords, expanded_coords)

   implicit none

!  Arguments:

   integer, intent (in)  :: ndim                ! Dimensionality of the table
   integer, intent (in)  :: dim1, dim2, dim3    ! That number of dimensions
   real,    intent (in)  :: &                   ! Compact form of coordinates
      packed_coords(dim1+dim2+dim3)
   real,    intent (out) :: &                   ! Expanded form
      expanded_coords(ndim,dim1,dim2,dim3)

!  Local variables:

   integer :: i, j, k
   integer ::    jo,ko  ! Offsets

!  Execution:

   jo = dim1;  ko = jo + dim2

   do k = 1, dim3
      do j = 1, dim2
         do i = 1, dim1
            expanded_coords(1,i,j,k) = packed_coords(i)
            expanded_coords(2,i,j,k) = packed_coords(jo+j)
            expanded_coords(3,i,j,k) = packed_coords(ko+k)
         end do
      end do
   end do

   end subroutine expand_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine expand_4d (ndim, dim1, dim2, dim3, dim4, packed_coords,          &
                         expanded_coords)

   implicit none

!  Arguments:

   integer, intent (in)  :: ndim                ! Dimensionality of the table
   integer, intent (in)  :: &                   ! That number of dimensions
      dim1, dim2, dim3, dim4
   real,    intent (in)  :: &                   ! Compact form of coordinates
      packed_coords(dim1+dim2+dim3+dim4)
   real,    intent (out) :: &                   ! Expanded form
      expanded_coords(ndim,dim1,dim2,dim3,dim4)

!  Local variables:

   integer :: i, j, k, l
   integer ::    jo,ko,lo  ! Offsets

!  Execution:

   jo = dim1;  ko = jo + dim2;  lo = ko + dim3

   do l = 1, dim4
      do k = 1, dim3
         do j = 1, dim2
            do i = 1, dim1
               expanded_coords(1,i,j,k,l) = packed_coords(i)
               expanded_coords(2,i,j,k,l) = packed_coords(jo+j)
               expanded_coords(3,i,j,k,l) = packed_coords(ko+k)
               expanded_coords(4,i,j,k,l) = packed_coords(lo+l)
            end do
         end do
      end do
   end do

   end subroutine expand_4d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine expand_5d (ndim, dim1, dim2, dim3, dim4, dim5, packed_coords,    &
                         expanded_coords)

   implicit none

!  Arguments:

   integer, intent (in)  :: ndim                ! Dimensionality of the table
   integer, intent (in)  :: &                   ! That number of dimensions
      dim1, dim2, dim3, dim4, dim5
   real,    intent (in)  :: &                   ! Compact form of coordinates
      packed_coords(dim1+dim2+dim3+dim4+dim5)
   real,    intent (out) :: &                   ! Expanded form
      expanded_coords(ndim,dim1,dim2,dim3,dim4,dim5)

!  Local variables:

   integer :: i, j, k, l, m
   integer ::    jo,ko,lo,mo  ! Offsets

!  Execution:

   jo = dim1;  ko = jo + dim2;  lo = ko + dim3;  mo = lo + dim4

   do m = 1, dim5
      do l = 1, dim4
         do k = 1, dim3
            do j = 1, dim2
               do i = 1, dim1
                  expanded_coords(1,i,j,k,l,m) = packed_coords(i)
                  expanded_coords(2,i,j,k,l,m) = packed_coords(jo+j)
                  expanded_coords(3,i,j,k,l,m) = packed_coords(ko+k)
                  expanded_coords(4,i,j,k,l,m) = packed_coords(lo+l)
                  expanded_coords(5,i,j,k,l,m) = packed_coords(mo+m)
               end do
            end do
         end do
      end do
   end do

   end subroutine expand_5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine expand_6d (ndim,  dim1, dim2, dim3, dim4, dim5, dim6,            &
                         packed_coords, expanded_coords)

   implicit none

!  Arguments:

   integer, intent (in)  :: ndim                ! Dimensionality of the table
   integer, intent (in)  :: &                   ! That number of dimensions
      dim1, dim2, dim3, dim4, dim5, dim6
   real,    intent (in)  :: &                   ! Compact form of coordinates
      packed_coords(dim1+dim2+dim3+dim4+dim5+dim6)
   real,    intent (out) :: &                   ! Expanded form
      expanded_coords(ndim,dim1,dim2,dim3,dim4,dim5,dim6)

!  Local variables:

   integer :: i, j, k, l, m, n
   integer ::    jo,ko,lo,mo,no  ! Offsets

!  Execution:

   jo = dim1;  ko = jo + dim2;  lo = ko + dim3;  mo = lo + dim4;  no = lo + dim5

   do n = 1, dim6
      do m = 1, dim5
         do l = 1, dim4
            do k = 1, dim3
               do j = 1, dim2
                  do i = 1, dim1
                     expanded_coords(1,i,j,k,l,m,n) = packed_coords(i)
                     expanded_coords(2,i,j,k,l,m,n) = packed_coords(jo+j)
                     expanded_coords(3,i,j,k,l,m,n) = packed_coords(ko+k)
                     expanded_coords(4,i,j,k,l,m,n) = packed_coords(lo+l)
                     expanded_coords(5,i,j,k,l,m,n) = packed_coords(mo+m)
                     expanded_coords(6,i,j,k,l,m,n) = packed_coords(no+n)
                  end do
               end do
            end do
         end do
      end do
   end do

   end subroutine expand_6d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
