!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine swap_ijk_conventions (mode, nf, ni, nj, nk, x, y, z, f)
!
!  OUTBOUND needs this step to allow the key steps to work with DPLR's choice
!  of (n,k,j,i) indexing.
!
!  Note that x, y, z, f are reindexed in-place, so this routine should be
!  compiled WITHOUT array bounds checking, because the array dimensions are
!  changed by the routine.
!
!  08/17/05  David Saunders  Preferable to having two versions of PERMUTE_BLOCK.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in) :: &
      mode,                &  ! 1 = i,j,k -> k,j,i;  2 = k,j,i -> i,j,k
      nf,                  &  ! Number of functions in the flow array
      ni, nj, nk              ! Block dimensions in one order or the other

   real, pointer, dimension (:,:,:) :: &
      x, y, z                 ! Packed block coordinates

   real, pointer :: &
      f(:,:,:,:)              ! Associated function data to be permuted if
                              ! mode = 1; add it for mode = 2 if needed
!  Local variables:

   integer :: &
      i, j, k, n
   real, allocatable, dimension (:,:,:) :: &
      xt, yt, zt
   real, allocatable, dimension (:,:,:,:) :: &
      ft

!  Execution:


   select case (mode)

      case (1) ! (i,j,k) -> (k,j,i)

         allocate (xt(nk,nj,ni), yt(nk,nj,ni), zt(nk,nj,ni), ft(nf,nk,nj,ni))

!        Fill the temporary arrays in the DPLR order:

         do k = 1, nk
            do j = 1, nj
               do i = 1, ni
                  xt(k,j,i) = x(i,j,k)
                  yt(k,j,i) = y(i,j,k)
                  zt(k,j,i) = z(i,j,k)
                  do n = 1, nf
                     ft(n,k,j,i) = f(n,i,j,k)
                  end do
               end do
            end do
         end do

!        Deallocating and redimensioning seems to be the only way for pointers:

         deallocate (x, y, z, f)
         allocate   (x(nk,nj,ni), y(nk,nj,ni), z(nk,nj,ni), f(nf,nk,nj,ni))

         do i = 1, ni
            do j = 1, nj
               do k = 1, nk
                  x(k,j,i) = xt(k,j,i)
                  y(k,j,i) = yt(k,j,i)
                  z(k,j,i) = zt(k,j,i)
                  do n = 1, nf
                     f(n,k,j,i) = ft(n,k,j,i)
                  end do
               end do
            end do
         end do

         deallocate (xt, yt, zt, ft)

      case (2) ! (k,j,i) -> (i,j,k); no need to treat f()

         allocate (xt(ni,nj,nk), yt(ni,nj,nk), zt(ni,nj,nk))

         do i = 1, ni 
            do j = 1, nj
               do k = 1, nk
                  xt(i,j,k) = x(k,j,i)
                  yt(i,j,k) = y(k,j,i)
                  zt(i,j,k) = z(k,j,i)
               end do
            end do
         end do

         deallocate (x, y, z)
         allocate   (x(ni,nj,nk), y(ni,nj,nk), z(ni,nj,nk))

         do k = 1, nk
            do j = 1, nj
               do i = 1, ni
                  x(i,j,k) = xt(i,j,k)
                  y(i,j,k) = yt(i,j,k)
                  z(i,j,k) = zt(i,j,k)
               end do
            end do
         end do

         deallocate (xt, yt, zt)

      case default

   end select

   end subroutine swap_ijk_conventions
