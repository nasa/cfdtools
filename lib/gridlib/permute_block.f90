!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine permute_block (mode, nf, mi, mj, mk, num_face_on_body, x, y, z,  &
                             f, ni, nj, nk)
!
!  Permute the indices of a grid block to make k = 1 the indicated face, or
!  reverse the permutation.  The permuted block is returned in place.  An
!  associated function array is similarly permuted.
!
!  N.B.:  Packed arrays are assumed throughout.
!
!  03/03/04  David Saunders  Initial implementation (in-place permutation) for
!                            GRID_MORPH_1 application.
!  05/31/05    "      "      Added associated function data (mode = 1 only) for
!                            outer grid boundary adaptation purposes.
!  08/17/05    "      "      Revised to assume DPLR's (n,k,j,i) convention.
!  08/19/05    "      "      Reduce memory usage by doing one array at a time.
!  08/22/05    "      "      Conventional (i,j,k) version.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in) :: &
      mode,                &  ! 1 = initial permutation; 2 = reverse permutation
      nf,                  &  ! Number of functions in the associated flow array
      mi, mj, mk,          &  ! Input block dimensions
      num_face_on_body        ! Input face number to change to k = 1 if mode = 1
                              ! or to change from k = 1 if mode = 2;
                              ! 1 means the i = 1 plane, 2 means the i = i(max)
                              ! plane, 3 means the j = 1 plane, and so on
   real, intent (inout),   &
      dimension (mi,mj,mk) :: &
      x, y, z                 ! Packed block coordinates; input in the order
                              ! implied by mi, mj, mk;    output with permuted
                              ! indices in the order implied by ni, nj, nk
   real, intent (inout) :: &
      f(nf,mi,mj,mk)          ! Associated function data to be permuted if
                              ! mode = 1; add it for mode = 2 if needed
   integer, intent (out) :: &
      ni, nj, nk              ! Dimensions of permuted block; when mode = 2,
                              ! ni/j/k swap roles with mi/j/k in the call

!  Local variables:

   integer :: &
      i, j, k, i1, i2, ii, j1, j2, jj, k1, k2, kk, knc, n, npt

   real, allocatable :: &
      xt(:,:,:), ft(:,:,:,:)

!  Execution:

   npt = mi*mj*mk

   if (mode == 1) then ! Permute indices to achieve the desired face at k = 1

      select case (num_face_on_body)

      case (1, 2) ! The i = 1 or i = mi plane becomes the k = 1 plane

         ni = mk;  nj = mj;  nk = mi

         if (num_face_on_body == 1) then
            k1 = 1;  k2 = nk;  knc =  1
         else
            k1 = nk;  k2 = 1;  knc = -1
         end if

         allocate (xt(ni,nj,nk))

!        Duplicate loops for x, y, z, f to reduce memory usage:

         ii = 0
         do k = k1, k2, knc
            ii = ii + 1
            do j = 1, nj
               do i = 1, ni
                  xt(i,j,k) = x(ii,j,i)
               end do
            end do
         end do

         call copy_as_vector (npt, xt, x)

         ii = 0
         do k = k1, k2, knc
            ii = ii + 1
            do j = 1, nj
               do i = 1, ni
                  xt(i,j,k) = y(ii,j,i)
               end do
            end do
         end do

         call copy_as_vector (npt, xt, y)

         ii = 0
         do k = k1, k2, knc
            ii = ii + 1
            do j = 1, nj
               do i = 1, ni
                  xt(i,j,k) = z(ii,j,i)
               end do
            end do
         end do

         call copy_as_vector (npt, xt, z)

         deallocate (xt)

         allocate (ft(nf,ni,nj,nk))

         ii = 0
         do k = k1, k2, knc
            ii = ii + 1
            do j = 1, nj
               do i = 1, ni
                  do n = 1, nf
                     ft(n,i,j,k) = f(n,ii,j,i)
                  end do
               end do
            end do
         end do

         call copy_as_vector (npt*nf, ft, f)

         deallocate (ft)

      case (3, 4) ! The j = 1 or j = mj plane becomes the k = 1 plane

         ni = mi;  nj = mk;  nk = mj

         if (num_face_on_body == 3) then
            k1 = 1;  k2 = nk;  knc =  1
         else
            k1 = nk;  k2 = 1;  knc = -1
         end if

         allocate (xt(ni,nj,nk))

         jj = 0
         do k = k1, k2, knc
            jj = jj + 1
            do j = 1, nj
               do i = 1, ni
                  xt(i,j,k) = x(i,jj,j)
               end do
            end do
         end do

         call copy_as_vector (npt, xt, x)

         jj = 0 
         do k = k1, k2, knc
            jj = jj + 1
            do j = 1, nj
               do i = 1, ni
                  xt(i,j,k) = y(i,jj,j)
               end do 
            end do
         end do

         call copy_as_vector (npt, xt, y)

         jj = 0 
         do k = k1, k2, knc
            jj = jj + 1
            do j = 1, nj
               do i = 1, ni
                  xt(i,j,k) = z(i,jj,j)
               end do 
            end do
         end do

         call copy_as_vector (npt, xt, z)

         deallocate (xt)

         allocate (ft(nf,ni,nj,nk))

         jj = 0
         do k = k1, k2, knc
            jj = jj + 1
            do j = 1, nj
               do i = 1, ni
                  do n = 1, nf
                     ft(n,i,j,k) = f(n,i,jj,j)
                  end do
               end do
            end do
         end do

         call copy_as_vector (npt*nf, ft, f)

         deallocate (ft)

      case (5, 6) ! The k = 1 plane is already face 5, or k = mk becomes k = 1

         ni = mi;  nj = mj;  nk = mk

         if (num_face_on_body == 5) then
!           Do nothing
         else

            allocate (xt(ni,nj,nk))

            kk = mk
            do k = 1, nk
               do j = 1, nj
                  do i = 1, ni
                     xt(i,j,k) = x(i,j,kk)
                  end do
               end do
               kk = kk - 1
            end do

            call copy_as_vector (npt, xt, x)

            kk = mk
            do k = 1, nk 
               do j = 1, nj
                  do i = 1, ni
                     xt(i,j,k) = y(i,j,kk)
                  end do
               end do
               kk = kk - 1
            end do

            call copy_as_vector (npt, xt, y)

            kk = mk
            do k = 1, nk 
               do j = 1, nj
                  do i = 1, ni
                     xt(i,j,k) = z(i,j,kk)
                  end do
               end do
               kk = kk - 1
            end do

            call copy_as_vector (npt, xt, z)

            deallocate (xt)

            allocate (ft(nf,ni,nj,nk))

            kk = mk
            do k = 1, nk
               do j = 1, nj
                  do i = 1, ni
                     do n = 1, nf
                        ft(n,i,j,k) = f(n,i,j,kk)
                     end do
                  end do
               end do
               kk = kk - 1
            end do

            call copy_as_vector (npt*nf, ft, f)

            deallocate (ft)

         end if

      end select ! Mode = 1 options (function array included)

   else ! Reverse mode = 2 (no need for permuting the function data)

     select case (num_face_on_body)

      case (1, 2) ! The k = 1 plane becomes the i = 1 or i = ni plane

         ni = mk;  nj = mj;  nk = mi

         if (num_face_on_body == 1) then
            k1 = 1;  k2 = mk;  knc =  1
         else
            k1 = mk;  k2 = 1;  knc = -1
         end if

         allocate (xt(ni,nj,nk))

         ii = 0
         do k = k1, k2, knc
            ii = ii + 1
            do j = 1, mj
               do i = 1, mi
                  xt(ii,j,i) = x(i,j,k)
               end do
            end do
         end do

         call copy_as_vector (npt, xt, x)

         ii = 0
         do k = k1, k2, knc
            ii = ii + 1
            do j = 1, mj
               do i = 1, mi
                  xt(ii,j,i) = y(i,j,k)
               end do
            end do
         end do

         call copy_as_vector (npt, xt, y)

         ii = 0
         do k = k1, k2, knc
            ii = ii + 1
            do j = 1, mj
               do i = 1, mi
                  xt(ii,j,i) = z(i,j,k)
               end do
            end do
         end do

         call copy_as_vector (npt, xt, z)

         deallocate (xt)

      case (3, 4) ! The k = 1 plane becomes the j = 1 or j = nj plane

         ni = mi;  nj = mk;  nk = mj

         if (num_face_on_body == 3) then
            k1 = 1;  k2 = mk;  knc =  1
         else
            k1 = mk;  k2 = 1;  knc = -1
         end if

         allocate (xt(ni,nj,nk))

         jj = 0
         do k = k1, k2, knc
            jj = jj + 1
            do j = 1, mj
               do i = 1, mi
                  xt(i,jj,j) = x(i,j,k)
               end do
            end do
         end do

         call copy_as_vector (npt, xt, x)

         jj = 0
         do k = k1, k2, knc
            jj = jj + 1
            do j = 1, mj
               do i = 1, mi
                  xt(i,jj,j) = y(i,j,k) 
               end do
            end do
         end do

         call copy_as_vector (npt, xt, y)

         jj = 0
         do k = k1, k2, knc
            jj = jj + 1
            do j = 1, mj
               do i = 1, mi
                  xt(i,jj,j) = z(i,j,k) 
               end do
            end do
         end do

         call copy_as_vector (npt, xt, z)

         deallocate (xt)

      case (5, 6) ! The k = 1 plane stays that way or becomes the k = nk plane

         ni = mi;  nj = mj;  nk = mk

         if (num_face_on_body == 5) then
!           Do nothing
         else

            allocate (xt(ni,nj,nk))

            kk = mk
            do k = 1, mk
               do j = 1, mj
                  do i = 1, mi
                     xt(i,j,kk) = x(i,j,k)
                  end do
               end do
               kk = kk - 1
            end do

            call copy_as_vector (npt, xt, x)

            kk = mk
            do k = 1, mk 
               do j = 1, mj
                  do i = 1, mi
                     xt(i,j,kk) = y(i,j,k)
                  end do
               end do
               kk = kk - 1
            end do

            call copy_as_vector (npt, xt, y)

            kk = mk
            do k = 1, mk 
               do j = 1, mj
                  do i = 1, mi
                     xt(i,j,kk) = z(i,j,k)
                  end do
               end do
               kk = kk - 1
            end do

            call copy_as_vector (npt, xt, z)

            deallocate (xt)

         end if

      end select

   end if

   end subroutine permute_block

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine copy_as_vector (npt, xin, xout)

!  Copy a packed multidimensional array as a vector.
!
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Arguments:

   integer, intent (in)  :: npt
   real,    intent (in)  :: xin(npt)
   real,    intent (out) :: xout(npt)

!  Execution:

   xout = xin

   end subroutine copy_as_vector
