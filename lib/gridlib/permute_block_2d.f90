!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine permute_block_2d (mode, nf, mi, mj, num_edge_on_body, x, y, f,   &
                                ni, nj)
!
!  Permute the indices of a grid block to make j = 1 the indicated edge, or
!  reverse the permutation.  The permuted block is returned in place.  An
!  associated function array is similarly permuted.
!
!  N.B.:  Packed arrays are assumed throughout.
!
!  11.05/05  D. A. Saunders  2-D version of permute_block.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in) :: &
      mode,                &  ! 1 = initial permutation; 2 = reverse permutation
      nf,                  &  ! Number of functions in the associated flow array
      mi, mj,              &  ! Input block dimensions
      num_edge_on_body        ! Input edge number to change to j = 1 if mode = 1
                              ! or to change from j = 1 if mode = 2;
                              ! 1 means the i = 1 line, 2 means the i = i(max)
                              ! line, 3 means the j = 1 line, and so on
   real, intent (inout),   &
      dimension (mi,mj) :: &
      x, y                    ! Packed block coordinates; input in the order
                              ! implied by mi, mj;    output with permuted
                              ! indices in the order implied by ni, nj
   real, intent (inout) :: &
      f(nf,mi,mj)             ! Associated function data to be permuted if
                              ! mode = 1; add it for mode = 2 if needed
   integer, intent (out) :: &
      ni, nj                  ! Dimensions of permuted block; when mode = 2,
                              ! ni/j swap roles with mi/j in the call

!  Local variables:

   integer :: &
      i, j, i1, i2, ii, j1, j2, jj, jnc, n, npt

   real, allocatable :: &
      xt(:,:), ft(:,:,:)

!  Execution:

   npt = mi*mj

   if (mode == 1) then ! Permute indices to achieve the desired edge at j = 1

      select case (num_edge_on_body)

      case (1, 2) ! The i = 1 or i = mi line becomes the j = 1 line

         ni = mj;  nj = mi

         if (num_edge_on_body == 1) then
            j1 = 1;  j2 = nj;  jnc =  1
         else
            j1 = nj;  j2 = 1;  jnc = -1
         end if

         allocate (xt(ni,nj))

!        Duplicate loops for x, y, f to reduce memory usage:

         ii = 0
         do j = j1, j2, jnc
            ii = ii + 1
            do i = 1, ni
               xt(i,j) = x(ii,i)
            end do
         end do

         call copy_as_vector_2d (npt, xt, x)

         ii = 0
         do j = j1, j2, jnc
            ii = ii + 1
            do i = 1, ni
               xt(i,j) = y(ii,i)
            end do
         end do

         call copy_as_vector_2d (npt, xt, y)

         allocate (ft(nf,ni,nj))

         ii = 0
         do j = j1, j2, jnc
            ii = ii + 1
            do i = 1, ni
               do n = 1, nf
                  ft(n,i,j) = f(n,ii,i)
               end do
            end do
         end do

         call copy_as_vector_2d (npt*nf, ft, f)

         deallocate (ft)

      case (3, 4) ! The j = 1 line is already edge 3, or j = mj becomes j = 1

         ni = mi;  nj = mj

         if (num_edge_on_body == 3) then
!           Do nothing
         else

            allocate (xt(ni,nj))

            jj = mj
            do j = 1, nj
               do i = 1, ni
                  xt(i,j) = x(i,jj)
               end do
               jj = jj - 1
            end do

            call copy_as_vector_2d (npt, xt, x)

            jj = mj
            do j = 1, nj
               do i = 1, ni
                  xt(i,j) = y(i,jj)
               end do
               jj = jj - 1
            end do

            call copy_as_vector_2d (npt, xt, y)

            deallocate (xt)

            allocate (ft(nf,ni,nj))

            jj = mj
            do j = 1, nj
               do i = 1, ni
                  do n = 1, nf
                     ft(n,i,j) = f(n,i,jj)
                  end do
               end do
               jj = jj - 1
            end do

            call copy_as_vector_2d (npt*nf, ft, f)

            deallocate (ft)

         end if

      end select ! Mode = 1 options (function array included)

   else ! Reverse mode = 2 (no need for permuting the function data)

     select case (num_edge_on_body)

      case (1, 2) ! The j = 1 line becomes the i = 1 or i = ni line

         ni = mj;  nj = mi

         if (num_edge_on_body == 1) then
            j1 = 1;  j2 = mj;  jnc =  1
         else
            j1 = mj;  j2 = 1;  jnc = -1
         end if

         allocate (xt(ni,nj))

         ii = 0
         do j = j1, j2, jnc
            ii = ii + 1
            do i = 1, mi
               xt(ii,i) = x(i,j)
            end do
         end do

         call copy_as_vector_2d (npt, xt, x)

         ii = 0
         do j = j1, j2, jnc
            ii = ii + 1
            do i = 1, mi
               xt(ii,i) = y(i,j)
            end do
         end do

         call copy_as_vector_2d (npt, xt, y)

         deallocate (xt)

      case (3, 4) ! The j = 1 line stays that way or becomes the j = nj line

         ni = mi;  nj = mj

         if (num_edge_on_body == 3) then
!           Do nothing
         else

            allocate (xt(ni,nj))

            jj = mj
            do j = 1, mj
               do i = 1, mi
                  xt(i,jj) = x(i,j)
               end do
               jj = jj - 1
            end do

            call copy_as_vector_2d (npt, xt, x)

            jj = mj
            do j = 1, mj 
               do i = 1, mi
                  xt(i,jj) = y(i,j)
               end do
               jj = jj - 1
            end do

            call copy_as_vector_2d (npt, xt, y)

            deallocate (xt)

         end if

      end select

   end if

   end subroutine permute_block_2d

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine copy_as_vector_2d (npt, xin, xout)

!  Copy a packed multidimensional array as a vector.
!
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Arguments:

   integer, intent (in)  :: npt
   real,    intent (in)  :: xin(npt)
   real,    intent (out) :: xout(npt)

!  Execution:

   xout = xin

   end subroutine copy_as_vector_2d
