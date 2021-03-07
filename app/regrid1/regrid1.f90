!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program regrid1

!  Description:
!
!     REGRID1 is a specialized gridding application intended to adjust a 2-D
!     grid representing an arc-jet test article at zero angle of attack.  Just
!     the upper half of the cross-section is treated.  Given a 2-space curve
!     representing the test surface after some ablation, the program performs
!     regridding in a way that ensures orthogonality at the boundaries as far
!     as possible.  [The new wall may be omitted, in which case the input grid
!     is processed as for the ablation case.]
!
!     In order to use elliptic smoothing, we work initially with Euler-type
!     point distributions before reimposing the original type of stretching
!     in the radial direction.
!
!  Assumptions:
!
!     o  Single block, dimensions ni x nj, with i = 1 at the stagnation point
!        and j = 1 at the wall.
!     o  The ablation is moderate enough that the outer boundary does not need
!        to move.
!
!  Outline:
!
!     o  Read the existing grid (formatted, PLOT3D, 2-D, 1 block).
!     o  Read the ablated surface data (2-column ASCII, to EOF), if any.
!     o  Intersect the centerline and outflow boundaries with the new wall.
!     o  Transfer the old wall distribution to the [relevant part of] new wall.
!     o  Impose a uniform distribution along the new centerline boundary.
!     o  For the new outflow boundary, construct a curve from the new (ni,1)
!        point to the old (in,ni) point, where "in" is input to facilitate
!        achieving orthogonality at the outer boundary.  The curve should be
!        normal to both j boundaries.  It is used to impose another uniform
!        distribution for the new outflow boundary.
!     o  Redistribute the outer boundary points from (1,ni) to (in,ni) (if in
!        < ni) by using a Vinokur distribution with end increments from the old
!        outer boundary.
!     o  Fill the interim interior points via TFI2D.
!     o  Smooth them via ELLIP2D.
!     o  Redistribute all radial lines via Vinokur and the original end deltas.
!
!  History:
!
!     10/03/05  D. Saunders  Initial implementation (uniform distributions).
!     10/05/05   "     "     Radial redistribution options.  Lots of elliptic
!                            smoothing iterations are required to overcome the
!                            relatively poor TFI starting guess.
!     10/06/05   "     "     Allow the new wall to be optional, so that a grid
!                            prior to ablation can be processed the same way.
!     12/21/05   "     "     Avoid possible Zs in the new wall file.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      lungridin  = 1,    &  ! Input grid
      lungridout = 2,    &  ! Output grid
      lunwall    = 3,    &  ! XY curve defining new wall
      lunkbd     = 5,    &  ! Keyboard entries
      luncrt     = 6,    &  ! Screen
      ncubic     = 100      ! For interim definition of new outflow boundary

   real, parameter :: &
      half = 0.5, one = 1., zero = 0.

   logical, parameter :: &
      false = .false., true = .true.

   character, parameter :: &
      method * 1 = 'M'      ! Monotonic spline fits

!  Variables:

   integer :: &
      i, i1, i2, ier, in, ios, itmax, j, j1, j2, jb, nb, ni, nj, nwall

   real :: &
      dmax, dt, dt1, dt2, dx, dy, r, tint1, tint2, tol, total, tscale, tshift, &
      xint, yint, ymax, ymin

   real, dimension (ncubic) ::   &
      tcubic, xcubic, ycubic

   logical :: &
      show_smoothing_iters

   character :: &
      filename * 64

!  ELLIP2D variables other than itmax:

   character, parameter :: &
      bgmode * 4 = 'YYYY', &
      fgmode * 4 = 'YNYY', &
      spmode * 4 = 'NNNN'

   integer, parameter ::   &
      itfloat  = 0,        &
      itfreeze = 2000,     &
      jlayer   = 5

   real, parameter ::      &
      convg    = 7.0,      &
      convmin  = 0.0,      &
      dmaxrel  = 1.e-8,    &
      omega    = 1.3,      &
      poweri   = 1.0,      &
      powerj   = 1.0,      &
      expi1    = 0.2,      &
      expi2    = 0.2,      &
      expj1    = 0.2,      &
      expj2    = 0.2,      &
      fglimit  = 1.0,      &
      urfg     = 0.1,      &
      urfloat  = 0.

!  Dynamic arrays:

   real, allocatable, dimension (:) :: &
      derivsi, derivsj, tidir, tidir2, tjdir, tjdir2, twall, xjdir, yjdir,     &
      xjdir2, yjdir2, xwall, ywall

   real, allocatable, dimension (:,:) :: &
      xold, yold, xnew, ynew

   real, allocatable, dimension (:,:,:) :: &
      arcs

!  Execution:

!  Read the existing grid (formatted, PLOT3D, 2-D, 1 block).

   write (luncrt, '(/, a)', advance='no') ' Input grid name:  '
   read  (lunkbd, *) filename

   open  (lungridin, file=filename, status='old')

   read  (lungridin, *) nb  ! 1 for now
   read  (lungridin, *) ni, nj
   allocate (xold(ni,nj), yold(ni,nj))
   read  (lungridin, *) xold, yold
   close (lungridin)

   ymax = maxval (yold)
   ymin = minval (yold)
   dmax = dmaxrel * (ymax - ymin) ! ELLIP2D convergence tolerance on dx, dy
!! write (luncrt, *) 'dmax:', dmax

   write (luncrt, '(a, 2i5, a)', advance='no') &
      ' Dimensions found:', ni, nj, '.  I index of new outflow corner point: '
   read  (lunkbd, *) in

   write (luncrt, '(a)', advance='no') &
      ' Number of elliptic smoothing iterations (< 0 shows iterations): '
   read  (lunkbd, *) itmax

   show_smoothing_iters = itmax < 0
   if (show_smoothing_iters) itmax = abs (itmax)

!  Read the ablated surface data (2-column ASCII), if any:

   write (luncrt, '(a)', advance='no') &
      ' New wall curve file name ("none" means same wall): '
   read  (lunkbd, *) filename

   if (filename(1:4) == 'none') then

      nwall = ni
      allocate (xwall(nwall), ywall(nwall), twall(nwall))
      do i = 1, nwall
         xwall(i) = xold(i,1)
         ywall(i) = yold(i,1)
      end do

   else

      open  (lunwall, file=filename, status='old')

      nwall = 0
      do ! Until EOF
         read (lunwall, '(a)', iostat=ios)
         if (ios < 0) exit
         nwall = nwall + 1
      end do

      allocate (xwall(nwall), ywall(nwall), twall(nwall))
      rewind (lunwall)
      do i = 1, nwall
         read (lunwall, *) xwall(i), ywall(i) ! Avoid possible Z coordinate
      end do
      close  (lunwall)

!     Ensure that the new wall matches the assumption for the old wall indexing:

      if (ywall(nwall) < ywall(1)) then
         call rverse (nwall, xwall, xwall)
         call rverse (nwall, ywall, ywall)
      end if

   end if

   write (luncrt, '(a)', advance='no') ' Output grid name: '
   read  (lunkbd, *) filename

   open  (lungridout, file=filename, status='unknown')

!  Intersect the centerline and outflow boundaries with the new wall.

   tol = 100.* epsilon (tol)

   allocate (xjdir(nj),  yjdir(nj),  tjdir(nj),  &
             xjdir2(nj), yjdir2(nj), tjdir2(nj), tidir2(ni))

   write (luncrt, '(a)') ' Wall intersection diagnostics (possibly redundant):'

   do i = 1, ni, ni - 1

      do j = 1, nj
         xjdir(j) = xold(i,j)
         yjdir(j) = yold(i,j)
      end do

      if (i == 1) then ! Center line
         i1 = 1
      else             ! Outflow boundary
         i1 = nwall
      end if

      j1 = 1 ! These are starting guesses for the intersection indices

      call intsec2t (nwall, xwall, ywall, twall, i1, true,                     &
                     nj,    xjdir, yjdir, tjdir, j1, true,                     &
                     tol,   xint,  yint,  tint1, tint2, luncrt, ier)

      if (ier /= 0) then
         write (luncrt, '(/, a, i4)') &
            ' Intersection trouble.  i:', i, &
            ' IER from INTSEC2T:', ier
         go to 99
      end if

      tidir2(i) = tint1

      write (luncrt, '(a, 1p, 4e19.11)') &
         ' tint1, tint2, xint, yint: ', tint1, tint2, xint, yint

   end do

!  Transfer the old wall distribution to the relevant part of the new wall.
!  I.e., transform twall in [0, T] to [tint1(i=1), tint1(i=ni)] and evaluate
!  x/ywall at these ts to give the new wall:

   tscale = (tidir2(ni)   - tidir2(1)) / twall(nwall)
   tshift = (twall(nwall) * tidir2(1)) / twall(nwall)

   tidir2(:) = tscale * twall(:) + tshift

   allocate (xnew(ni,nj), ynew(ni,nj), derivsi(ni), derivsj(nj))

   call lcsfit (ni, twall, xwall, true, method, ni, tidir2, xnew(1,1), derivsi)
   call lcsfit (ni, twall, ywall, true, method, ni, tidir2, ynew(1,1), derivsi)

!  Impose a nonuniform distribution on the new centerline boundary:

   tjdir(1)  = zero
   tjdir(nj) = xnew(1,1) - xold(1,nj)
   dx = half * tjdir(nj) / real (nj - 1)   ! Half the uniform interval
   dt1 = dx  ! For reuse at outflow boundary

   call vinokur (1, nj, dx, dx, tjdir, luncrt, ier)

   do j = 2, nj
      xnew(1,j) = xnew(1,1) - tjdir(j)
      ynew(1,j) = zero
   end do

!  For the new outflow boundary, construct a curve from the new (ni,1)
!  point to the old (in,ni) point, where "in" is input to facilitate
!  achieving orthogonality at the outer boundary.  The curve should be
!  normal to both j boundaries.  A cubic should suffice.  It is used
!  to impose another interim distribution for the new outflow boundary.

   if (in < ni) then ! Strive for orthogonality at the outflow outer boundary:

!     Evaluate a cubic curve   x = ay**3 + by**2 + cy + d at uniform y.
!     CUBINTRP is intended for y = ax**3 + bx**2 + cx + d.

      xcubic(1) = xnew(ni,1);  xcubic(ncubic) = xold(in,nj)
      ycubic(1) = ynew(ni,1);  ycubic(ncubic) = yold(in,nj)

      dy = (ycubic(ncubic) - ycubic(1)) / real (ncubic - 1)

      do j = 2, ncubic - 1
         ycubic(j) = ycubic(1) + real (j-1) * dy
      end do
 
      dx = xcubic(ncubic) - xold(in-1,nj) ! Outer bndry. slope at in is dy/dx
      dy = ycubic(ncubic) - yold(in-1,nj) ! so normal slope is -dx/dy

      call cubintrp (ncubic-2, ycubic(2), ycubic(1), xcubic(1), zero,          &
                     ycubic(ncubic), xcubic(ncubic), -dy/dx, xcubic(2))

      call chords2d (ncubic, xcubic, ycubic, false, total, tcubic)

!!    write (6, '(a)') ' Cubic t, x, y:'
!!    write (6, '(1p, 3e19.11)') (tcubic(j), xcubic(j), ycubic(j), j = 1, ncubic)

!     Impose a non-uniform arc distribution on this cubic outflow boundary.
!     Keep the increment at the wall constant; relax the outer increment.

      tjdir(nj) = total
      dt2 = half * total / real (nj - 1)

      call vinokur (1, nj, dt1, dt2, tjdir, luncrt, ier)

      call lcsfit (ncubic, tcubic, xcubic, true, method, nj, tjdir, xjdir,      &
                   derivsj)
      call lcsfit (ncubic, tcubic, ycubic, true, method, nj, tjdir, yjdir,      &
                   derivsj)

!!    write (6, '(a)') ' New interim outflow boundary:'
!!    write (6, '(1p, 3e19.11)') (tjdir(j), xjdir(j), yjdir(j), j = 1, nj)
!!    write (6, '(a)') ' In loop:'
!!    call flush (6)

      do j = 1, nj
         xnew(ni,j) = xjdir(j)
         ynew(ni,j) = yjdir(j)
!!       write (6, '(1p, i4, 2e19.11)') j, xnew(ni,j), ynew(ni,j)
      end do

!     Redistribute the outer boundary points from (1,ni) to (in,ni).
!     --------------------------------------------------------------

!     Preserving the full original relative spacing can give poor results for
!     untailored grids with vertical outflow boundaries, so instead we blend
!     a uniform distribution with the wall distribution.

      allocate (tidir(ni))

      call chords2d (ni, xold(1,nj), yold(1,nj), false, total, tidir)

!!    write (6, '(a)') ' Chords 1:ni on outer boundary:'
!!    write (6, '(i4, 1p, e19.11)') (i, tidir(i), i = 1, ni)

      dt = tidir(in) / real (ni - 1)  ! Uniform new increment

      do i = 1, ni
         tidir2(i) = dt * real (i - 1)
      end do

!     Since this loses some correspondence between wall and outer boundary,
!     blend it moderately with the relative wall distribution:

      tscale = tidir2(ni) / twall(ni)
      twall(:) = tscale * twall(:)

!!    write (6, '(a)') ' i, tidir2, scaled twall'
!!    write (6, '(i4, 1p, 2e19.11)') (i, tidir2(i), twall(i), i = 1, ni)
!!    call flush (6)

      r = 0.1
      twall(:) = (one - r) * tidir2(:) + r * twall(:)

!!    write (6, '(a)') ' i, tidir2, updated twall'
!!    write (6, '(i4, 1p, 2e19.11)') (i, tidir2(i), twall(i), i = 1, ni)
!!    call flush (6)

      call lcsfit (in, tidir, xold(1,nj), true, method, ni, twall, xnew(1,nj), &
                   derivsi)
      call lcsfit (in, tidir, yold(1,nj), true, method, ni, twall, ynew(1,nj), &
                   derivsi)

!!    write (6, '(a)') ' t, x, y, yold along new outer boundary:'
!!    write (6, '(i4, 1p, 4e19.11)') &
!!       (i, twall(i), xnew(i,nj), ynew(i,nj), yold(i,nj), i = 1, ni)
!!    call flush (6)

   else ! in = ni: retain (most of) the old outflow boundary and outer boundary

!     Impose a non-uniform arc distribution on the old outflow boundary.
!     Keep the increment at the wall the same as at the centerline;
!     relax the outer increment.

      xjdir2(1) = xnew(ni,1)
      yjdir2(1) = ynew(ni,1)

      do j = 2, nj
         xjdir2(j) = xold(ni,j)
         yjdir2(j) = yold(ni,j)
      end do

      call chords2d (nj, xjdir2, yjdir2, false, total, tjdir2)

      tjdir(nj) = total
      dt2 = half * total / real (nj - 1)

      call vinokur (1, nj, dt1, dt2, tjdir, luncrt, ier) 

      call lcsfit (nj, tjdir2, xjdir2, true, method, nj, tjdir, xjdir, derivsj)
      call lcsfit (nj, tjdir2, yjdir2, true, method, nj, tjdir, yjdir, derivsj)

!!    write (6, '(a)') ' New interim outflow boundary:'
!!    write (6, '(1p, 3e19.11)') (tjdir(j), xjdir(j), yjdir(j), j = 1, nj)
!!    call flush (6)

      do j = 2, nj
         xnew(ni,j) = xjdir(j)
         ynew(ni,j) = yjdir(j)
      end do

!     Retain the outer boundary:

      do i = 2, ni - 1
         xnew(i,nj) = xold(i,nj)
         ynew(i,nj) = yold(i,nj)
      end do

   end if

!  Fill the interim interior points via TFI2D.

   allocate (arcs(ni,nj,2)) ! This big for ELLIP2D; just edges for TFI2D

   call tfi2d (ni, 1, ni, 1, nj, xnew, ynew, arcs)

!  Elliptic smoothing of interior points:

   call ellip2d (ni, nj, 1, ni, 1, nj, xnew, ynew, arcs, show_smoothing_iters, &
                 bgmode, fgmode, spmode, itmax, itfloat, itfreeze, jlayer,     &
                 convg, convmin, dmax, omega, poweri, powerj,                  &
                 expi1, expi2, expj1, expj2, fglimit, urfg, urfloat)


!  Redistribute all radial lines.  Several options may be desirable here.

!  End increments dt1 & dt2 at stag. point, in case we use them everywhere:

   call chordsrow (ni, nj, 1,    1,  2, xold, yold, false, dt1, tjdir)
   call chordsrow (ni, nj, 1, nj-1, nj, xold, yold, false, dt2, tjdir)

   tjdir2(1) = zero

   do i = 1, ni

      if (in == ni) then ! Cater to the original grid?  (Debatable.) 

!        Old grid's outermost increment dt2 for this i:

         call chordsrow (ni, nj, i, nj-1, nj, xold, yold, false, dt2, tjdir)

      else ! The old grid's dt2 can be too grossly different near i = ni

      end if

!     Gather the interim radial line:

      do j = 1, nj
         xjdir(j) = xnew(i,j)
         yjdir(j) = ynew(i,j)
      end do

      call chords2d (nj, xjdir, yjdir, false, total, tjdir)

!     New radial distribution:

      tjdir2(nj) = total

      call vinokur (1, nj, dt1, dt2, tjdir2, luncrt, ier) 

      call lcsfit (nj, tjdir, xjdir, true, method, nj, tjdir2, xjdir2, derivsj)
      call lcsfit (nj, tjdir, yjdir, true, method, nj, tjdir2, yjdir2, derivsj)

      do j = 2, nj - 1
         xnew(i,j) = xjdir2(j)
         ynew(i,j) = yjdir2(j)
      end do

   end do

   write (lungridout, '(i2, /, 2i5)') nb, ni, nj
   write (lungridout, '(1p, 6e19.11)') xnew
   write (lungridout, '(1p, 6e19.11)') ynew
   close (lungridout)

99 continue

   end program regrid1
