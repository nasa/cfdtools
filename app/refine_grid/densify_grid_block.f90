!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine densify_grid_block (ib, nf, method_i, method_j, method_k,        &
                                  fraction_i,  fraction_j,  fraction_k,        &
                                  ni, nj, nk,  xin,   yin,  zin,   qin,        &
                                  mi, mj, mk,  xout, yout,  zout, qout, ier)
!
!  Adjust the cell counts of the indicated grid block via 1-D cubic spline
!  interpolation of the coordinates in the i, j, k directions in that order,
!  each direction employing the new data from the previous direction.
!
!  The adjusted block dimensions must be known at the higher level for output
!  header info. to be written to disk.
!
!  The index directions may be treated differently, depending on method_*.
!  Either the input relative spacings in a given index direction are preserved
!  or the existing initial increments are preserved as would be desirable for
!  maintaining cell Reynolds numbers at a wall boundary.
!
!  Optional function data may also be interpolated similarly (nf > 0).
!  Separating this from the x/y/z interpolations would be inefficient.
!
!  04/24/06  D. A. Saunders  OUTBOUND's subroutine adjust_point_counts
!                            (k direction only).
!  11/19/09     "      "     Adaptation for REFINE_GRID (all three directions).
!  11/21/09     "      "     Introduced the cspline_module and applied it
!                            to the grid coordinates (only) to improve
!                            curvature continuity on geometric surfaces.
!                            Block boundary effects are still likely.
!  07/28/20     "      "     The algebraic method 1 does not preserve smoothly
!                            varying cell spacings.  Introduce the technique of
!                            ARBDIS in an attempt to do so.  It seems to work
!                            well with the algebraic results as a starting soln.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use cspline_module        ! Conventional cubic spline, repackaged

   implicit none

!  Arguments:

   integer, intent (in) :: &
      ib,                  &                   ! Block number, for diagnostics
      nf,                  &                   ! nf > 0 means interp. fn. data
      method_i, method_j, method_k             ! 1 = same relative spacing;
                                               ! 2 = same absolute ds1
   real, intent (in) :: &
      fraction_i, fraction_j, fraction_k       ! If method_* = 2, this fraction
                                               ! is applied to the last spacing
                                               ! from a preliminary one-sided
                                               ! stretching with ds1 preserved,
                                               ! defining a two-sided stretching
   integer, intent (in) :: &
      ni, nj, nk                               ! Input block dimensions

   real, intent (in),  dimension (ni,nj,nk) :: &
      xin, yin, zin                            ! Input grid block coordinates

   real, intent (in), dimension (nf,ni,nj,nk) :: &
      qin                                      ! Input function data if nf > 0

   integer, intent (in) :: &
      mi, mj, mk                               ! Output block dimensions

   real, intent (out), dimension (mi,mj,mk) :: &
      xout, yout, zout                         ! Redistributed coordinates

   real, intent (out), dimension (nf,mi,mj,mk) :: &
      qout                                     ! Interpolated fn. data if nf > 0

   integer, intent (out) :: &
      ier                                      ! ier /= 0 means an allocation
                                               ! error occurred

!  Local constants:

   integer,   parameter :: lunerr = 6
   integer,   parameter :: ngeometric = 2      ! No option for partial geometric
   real,      parameter :: growth_rate = 1.01  ! Not used
   real,      parameter :: half = 0.5, one = 1., zero = 0.
   logical,   parameter :: false  = .false., true = .true.
   character, parameter :: routine * 20 = 'densify_grid_block: '
   character, parameter :: method * 1 = 'M'    ! Monotonic spline fits safeguard
                                               ! the symmetry plane
!  Local variables:

   integer :: i, ieval, ip, j, jp, jeval, k, keval, kp, n
   real    :: d1, d2, p, r, rp, s_total
   real    :: derivs(3)
   logical :: flowdata, new
   real, allocatable, dimension (:)       :: sold, xold, yold, zold, qold
   real, allocatable, dimension (:)       :: snew, xnew, ynew, znew, qnew, fp
   real, allocatable, dimension (:)       :: smid_in, spacing_in, work
   real, allocatable, dimension (:,:,:)   :: xi, yi, zi, xij, yij, zij
   real, allocatable, dimension (:,:,:,:) :: qi, qij

!  Execution:

   flowdata = nf > 0

!  Hard-code the spline end-condition choices for now (for the grid lines):

   iendl = 0;  derivl = zero  ! 3rd deriv. matches that of cubic thru pts. 1-4
   iendr = 0;  derivr = zero  ! (We really need 4th-degree polynomials ...)

!  *** Phase 1: interpolations in the i direction. ***

   allocate (xi(mi,nj,nk), yi(mi,nj,nk), zi(mi,nj,nk), stat=ier)  ! Interim xyzs

   if (ier /= 0) then
      write (lunerr, '(/, 2a, i5)') &
         routine, 'Trouble allocating xi, yi, zi; block #:', ib
      write (lunerr, '(a, 3i5)') 'Dimensions: ', mi, nj, nk
      go to 99
   end if

   if (flowdata) then
      allocate (qi(nf,mi,nj,nk), stat=ier)
      if (ier /= 0) then
         write (lunerr, '(/, 2a, i5)') &
            routine, 'Trouble allocating qi; block #:', ib
         write (lunerr, '(a, i2, 3i5)') 'Dimensions: ', nf, mi, nj, nk
         go to 99
      end if
   end if

   if (ni == mi) then  ! Just transcribe original data as intermediate data

      xi(:,:,:) = xin(:,:,:)
      yi(:,:,:) = yin(:,:,:)
      zi(:,:,:) = zin(:,:,:)

      if (flowdata) qi(:,:,:,:) = qin(:,:,:,:)

   else

      allocate (xold(ni), yold(ni), zold(ni), sold(ni), &
                xnew(mi), ynew(mi), znew(mi), snew(mi))

      if (flowdata) allocate (qold(ni), qnew(mi), fp(mi))  ! fp(:) = 1st derivs.

      if (method_i == 1) then  ! Same relative spacing/different arc lengths
         r = real (ni - 1) / real (mi - 1)
         allocate (smid_in(ni+1), spacing_in(ni+1), work(5*mi))
      else
         ! Method 2: retain input d1s; adjust d2s appropriately
      end if

      do k = 1, nk

         do j = 1, nj

            xold(:) = xin(:,j,k)
            yold(:) = yin(:,j,k)
            zold(:) = zin(:,j,k)

            call chords3d (ni, xold, yold, zold, false, s_total, sold)

            select case (method_i)  ! For new arc length distribution

               case (1) ! Same relative distribution

                  snew(1) = zero
                  do i = 2, mi - 1
                     rp = one + r * real (i - 1)
                     ip = int (rp)
                     p  = rp - real (ip)
                     snew(i) = (one - p) * sold(ip) + p * sold(ip+1)
                  end do
                  snew(mi) = s_total

!                 The above does not guarantee preservation of the smoothly
!                 varying cell growth rates/cell spacings that are presumed
!                 to exist in the input grid.  Use the input cell spacings as
!                 the shape function ARBDIS employs to impose a desired type
!                 of point distribution (as in CURVDIS).  Moreover, use the
!                 above algebraic result as a starting guess:

                  do i = 2, ni
                     smid_in(i)    = (sold(i) + sold(i-1))*half
                     spacing_in(i) =  sold(i) - sold(i-1)
                  end do

!                 Avoid extrapolation off the ends of the shape function:

                  smid_in(1) = zero;  smid_in(ni+1) = s_total
                  call linear_interp (smid_in(2),    spacing_in(2),    &
                                      smid_in(3),    spacing_in(3),    &
                                      smid_in(1),    spacing_in(1))
                  call linear_interp (smid_in(ni-1), spacing_in(ni-1), &
                                      smid_in(ni),   spacing_in(ni),   &
                                      smid_in(ni+1), spacing_in(ni+1))

                  call arbdis (mi, zero, s_total, ni+1, smid_in, &
                              spacing_in, 'I', -lunerr, work, snew, ier)
                  if (ier /= 0) then
                     write (lunerr, '(a, i2, /, a, 3i5)') &
                        ' *** densify_grid_block: Bad ier from arbdis:', ier, &
                        '     i direction; block, j, k:', ib, j, k
                     go to 99
                  end if

               case (2) ! Existing d1 everywhere

                  d1 = sold(2) - sold(1)

!                 1-sided Vinokur stretching (geometric if 2nd ds < ds1):

                  call expdis5 (1, zero, s_total, d1, mi, snew, -lunerr)

                  d2 = (s_total - snew(mi-1)) * fraction_i

                  call blgrid (mi, d1, d2, ngeometric, growth_rate, snew,      &
                               lunerr, ier)

               case default ! Can't happen

            end select

!           Interpolate the (x,y,z)s as functions of the new arc lengths:

            new = true

            call cspline (new, ni, sold, xold, mi, snew, xi(1,j,k), ier)
            call cspline (new, ni, sold, yold, mi, snew, yi(1,j,k), ier)
            call cspline (new, ni, sold, zold, mi, snew, zi(1,j,k), ier)

!!          ieval = 1;  new = true
!!
!!          do i = 1, mi
!!
!!             call plscrv3d (ni, xold, yold, zold, sold, method, new, false,  &
!!                            snew(i), ieval, xnew(i), ynew(i), znew(i), derivs)
!!             new = false
!!             xi(i,j,k) = xnew(i)
!!             yi(i,j,k) = ynew(i)
!!             zi(i,j,k) = znew(i)
!!
!!          end do

            if (flowdata) then  ! Equivalent interpolations of function data

               new = true

               do n = 1, nf
                  qold(:) = qin(n,:,j,k)

                  call lcsfit (ni, sold, qold, new, method, mi, snew, qnew, fp)

                  qi(n,:,j,k) = qnew(:)
               end do

            end if

         end do ! Next j

      end do ! Next k

      deallocate (xold, yold, zold, sold, xnew, ynew, znew, snew)
      if (method_i == 1) deallocate (smid_in, spacing_in, work)

      if (flowdata) deallocate (qold, qnew, fp)

   end if  ! End of i direction interpolations


!  *** Phase 2: interpolations in the j direction. ***

   allocate (xij(mi,mj,nk), yij(mi,mj,nk), zij(mi,mj,nk), stat=ier)  ! Interim

   if (ier /= 0) then
      write (lunerr, '(/, 2a, i5)') &
         routine, 'Trouble allocating xij, yij, zij; block #:', ib
      write (lunerr, '(a, 3i5)') 'Dimensions: ', mi, mj, nk
      go to 99
   end if

   if (flowdata) then
      allocate (qij(nf,mi,mj,nk), stat=ier)
      if (ier /= 0) then
         write (lunerr, '(/, 2a, i5)') &
            routine, 'Trouble allocating qi; block #:', ib
         write (lunerr, '(a, i2, 3i5)') 'Dimensions: ', nf, mi, mj, nk
         go to 99
      end if
   end if

   if (nj == mj) then  ! Just transcribe phase 1 data as intermediate data

      xij(:,:,:) = xi(:,:,:)
      yij(:,:,:) = yi(:,:,:)
      zij(:,:,:) = zi(:,:,:)

      if (flowdata) qij(:,:,:,:) = qi(:,:,:,:)

   else

      allocate (xold(nj), yold(nj), zold(nj), sold(nj), &
                xnew(mj), ynew(mj), znew(mj), snew(mj))

      if (flowdata) allocate (qold(nj), qnew(mj), fp(mj))  ! fp(:) = 1st derivs.

      if (method_j == 1) then  ! Same relative spacing/different arc lengths
         r = real (nj - 1) / real (mj - 1)
         allocate (smid_in(nj+1), spacing_in(nj+1), work(5*mj))
      else
         ! Method 2: retain input d1s; adjust d2s appropriately
      end if

      do k = 1, nk

         do i = 1, mi

            xold(:) = xi(i,:,k)
            yold(:) = yi(i,:,k)
            zold(:) = zi(i,:,k)

            call chords3d (nj, xold, yold, zold, false, s_total, sold)

            select case (method_j)  ! For new arc length distribution

               case (1) ! Same relative distribution

                  snew(1) = zero
                  do j = 2, mj - 1
                     rp = one + r * real (j - 1)
                     jp = int (rp)
                     p  = rp - real (jp)
                     snew(j) = (one - p) * sold(jp) + p * sold(jp+1)
                  end do
                  snew(mj) = s_total

!                 The above does not guarantee preservation of the smoothly
!                 varying cell growth rates/cell spacings that are presumed
!                 to exist in the input grid.  Use the input cell spacings as
!                 the shape function ARBDIS employs to impose a desired type
!                 of point distribution (as in CURVDIS).  Moreover, use the
!                 above algebraic result as a starting guess:

                  do j = 2, nj
                     smid_in(j)    = (sold(j) + sold(j-1))*half
                     spacing_in(j) =  sold(j) - sold(j-1)
                  end do

!                 Avoid extrapolation off the ends of the shape function:

                  smid_in(1) = zero;  smid_in(nj+1) = s_total
                  call linear_interp (smid_in(2),    spacing_in(2),    &
                                      smid_in(3),    spacing_in(3),    &
                                      smid_in(1),    spacing_in(1))
                  call linear_interp (smid_in(nj-1), spacing_in(nj-1), &
                                      smid_in(nj),   spacing_in(nj),   &
                                      smid_in(nj+1), spacing_in(nj+1))

                  call arbdis (mj, zero, s_total, nj+1, smid_in, &
                              spacing_in, 'I', -lunerr, work, snew, ier)
                  if (ier /= 0) then
                     write (lunerr, '(a, i2, /, a, 3i5)') &
                        ' *** densify_grid_block: Bad ier from arbdis:', ier, &
                        '     j direction; block, k, i:', ib, k, i
                     go to 99
                  end if

               case (2) ! Existing d1 everywhere

                  d1 = sold(2) - sold(1)

!                 1-sided Vinokur stretching (geometric if 2nd ds < ds1):

                  call expdis5 (1, zero, s_total, d1, mj, snew, -lunerr)

                  d2 = (s_total - snew(mj-1)) * fraction_j

                  call blgrid (mj, d1, d2, ngeometric, growth_rate, snew,      &
                               lunerr, ier)

               case default ! Can't happen

            end select

!           Interpolate the (x,y,z)s as functions of the new arc lengths:

            new = true

            call cspline (new, nj, sold, xold, mj, snew, xnew, ier)
            xij(i,:,k) = xnew(:)

            call cspline (new, nj, sold, yold, mj, snew, ynew, ier)
            yij(i,:,k) = ynew(:)

            call cspline (new, nj, sold, zold, mj, snew, znew, ier)
            zij(i,:,k) = znew(:)


!!          jeval = 1;  new = true
!!
!!          do j = 1, mj
!!
!!             call plscrv3d (nj, xold, yold, zold, sold, method, new, false,  &
!!                            snew(j), jeval, xnew(j), ynew(j), znew(j), derivs)
!!             new = false
!!             xij(i,j,k) = xnew(j)
!!             yij(i,j,k) = ynew(j)
!!             zij(i,j,k) = znew(j)
!!
!!          end do

            if (flowdata) then  ! Equivalent interpolations of function data

               new = true

               do n = 1, nf
                  qold(:) = qi(n,i,:,k)

                  call lcsfit (nj, sold, qold, new, method, mj, snew, qnew, fp)

                  qij(n,i,:,k) = qnew(:)
               end do

            end if

         end do ! Next i

      end do ! Next k

      deallocate (xold, yold, zold, sold, xnew, ynew, znew, snew)
      if (method_j == 1) deallocate (smid_in, spacing_in, work)

      if (flowdata) deallocate (qold, qnew, fp)

   end if  ! End of j direction interpolations

   deallocate (xi, yi, zi)

   if (flowdata) deallocate (qi)


!  *** Phase 3: interpolations in the k direction. ***

   if (nk == 1 .or. nk == mk) then  ! Transcribe phase 2 as the final result

      xout(:,:,:) = xij(:,:,:)
      yout(:,:,:) = yij(:,:,:)
      zout(:,:,:) = zij(:,:,:)

      if (flowdata) qout(:,:,:,:) = qij(:,:,:,:)

   else

      allocate (xold(nk), yold(nk), zold(nk), sold(nk), &
                xnew(mk), ynew(mk), znew(mk), snew(mk))

      if (flowdata) allocate (qold(nk), qnew(mk), fp(mk))  ! fp(:) = 1st derivs.

      if (method_k == 1) then  ! Same relative spacing/different arc lengths
         r = real (nk - 1) / real (mk - 1)
         allocate (smid_in(nk+1), spacing_in(nk+1), work(5*mk))
      else
         ! Method 2: retain input d1s; adjust d2s appropriately
      end if

      do j = 1, mj

         do i = 1, mi

            xold(:) = xij(i,j,:)
            yold(:) = yij(i,j,:)
            zold(:) = zij(i,j,:)

            call chords3d (nk, xold, yold, zold, false, s_total, sold)

            select case (method_k)  ! For new arc length distribution

               case (1) ! Same relative distribution

                  snew(1) = zero
                  do k = 2, mk - 1
                     rp = one + r * real (k - 1)
                     kp = int (rp)
                     p  = rp - real (kp)
                     snew(k) = (one - p) * sold(kp) + p * sold(kp+1)
                  end do
                  snew(mk) = s_total

!                 The above does not guarantee preservation of the smoothly
!                 varying cell growth rates/cell spacings that are presumed
!                 to exist in the input grid.  Use the input cell spacings as
!                 the shape function ARBDIS employs to impose a desired type
!                 of point distribution (as in CURVDIS).  Moreover, use the
!                 above algebraic result as a starting guess:

                  do k = 2, nk
                     smid_in(k)    = (sold(k) + sold(k-1))*half
                     spacing_in(k) =  sold(k) - sold(k-1)
                  end do

!                 Avoid extrapolation off the ends of the shape function:

                  smid_in(1) = zero;  smid_in(nk+1) = s_total
                  call linear_interp (smid_in(2),    spacing_in(2),    &
                                      smid_in(3),    spacing_in(3),    &
                                      smid_in(1),    spacing_in(1))
                  call linear_interp (smid_in(nk-1), spacing_in(nk-1), &
                                      smid_in(nk),   spacing_in(nk),   &
                                      smid_in(nk+1), spacing_in(nk+1))
                  
                  call arbdis (mk, zero, s_total, nk+1, smid_in, &
                              spacing_in, 'I', -lunerr, work, snew, ier)
                  if (ier /= 0) then
                     write (lunerr, '(a, i2, /, a, 3i5)') &
                        ' *** densify_grid_block: Bad ier from arbdis:', ier, &
                        '     k direction; block, i, j:', ib, i, j
                     go to 99
                  end if

               case (2) ! Existing d1 everywhere

                  d1 = sold(2) - sold(1)

!                 1-sided Vinokur stretching (geometric if 2nd ds < ds1):

                  call expdis5 (1, zero, s_total, d1, mk, snew, -lunerr)

                  d2 = (s_total - snew(mk-1)) * fraction_k

                  call blgrid (mk, d1, d2, ngeometric, growth_rate, snew,      &
                               lunerr, ier)

               case default ! Can't happen

            end select

!           Interpolate the (x,y,z)s as functions of the new arc lengths:

            new = true

            call cspline (new, nk, sold, xold, mk, snew, xnew, ier)
            xout(i,j,:) = xnew(:)

            call cspline (new, nk, sold, yold, mk, snew, ynew, ier)
            yout(i,j,:) = ynew(:)

            call cspline (new, nk, sold, zold, mk, snew, znew, ier)
            zout(i,j,:) = znew(:)

!!          keval = 1;  new = true
!!
!!          do k = 1, mk
!!
!!             call plscrv3d (nk, xold, yold, zold, sold, method, new, false,  &
!!                            snew(k), keval, xnew(k), ynew(k), znew(k), derivs)
!!             new = false
!!             xout(i,j,k) = xnew(k)
!!             yout(i,j,k) = ynew(k)
!!             zout(i,j,k) = znew(k)
!!
!!          end do

            if (flowdata) then  ! Equivalent interpolations of function data

               new = true

               do n = 1, nf
                  qold(:) = qij(n,i,j,:)

                  call lcsfit (nk, sold, qold, new, method, mk, snew, qnew, fp)

                  qout(n,i,j,:) = qnew(:)
               end do

            end if

         end do ! Next i

      end do ! Next j

      deallocate (xold, yold, zold, sold, xnew, ynew, znew, snew)
      if (method_k == 1) deallocate (smid_in, spacing_in, work)

      if (flowdata) deallocate (qold, qnew, fp)

   end if  ! End of k direction interpolations

   deallocate (xij, yij, zij)

   if (flowdata) deallocate (qij)

99 return

!  Internal procedure:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine linear_interp (x1, y1, x2, y2, xinterp, yinterp)  ! Obvious ...
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      real, intent (in)  :: x1, y1, x2, y2, xinterp
      real, intent (out) :: yinterp

!     Local variables:

      real :: slope

!     Execution:

      slope   = (y2 - y1)/(x2 - x1)
      yinterp = y1 + slope*(xinterp - x1)

      end subroutine linear_interp

   end subroutine densify_grid_block
