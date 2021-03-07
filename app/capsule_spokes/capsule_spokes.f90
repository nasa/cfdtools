!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program capsule_spokes

!  Description:
!
!     This utility is intended to complement the POLAR_INTERP scheme by defining
!     body points along capsule heat shield (or aft body) spokes at which the
!     NEQAIR radiation solver is intended to be run.  POLAR_INTERP turns the
!     radiative heat fluxes at those body points into a structured surface data-
!     set by interpolating nonlinearly across spokes in the azimuthal direction
!     after regularizing the spoke point counts.  This in turn may be interp-
!     olated to a desired mesh (SURFACE_INTERP) or contour-plotted directly.
!     [The same scheme is also applicable to mass loss data gathered with
!     material ablation solver FIAT in place of NEQAIR.]
!
!  Assumptions:
!
!     >  The capsule geometry smooth outer mold line is axisymmetric.
!     >  The analysis is to be at nonzero angle of attack (NOT axisymmetric).
!     >  The supplied portion of a generatrix consists of the points desired
!        along the vertical symmetry plane slice.  It is expected to be confined
!        to either the forebody or to the aft body--not both, in order to avoid
!        complications involving duplicate end points common to all spokes.
!     >  Further spokes are defined simply by rotating this generatrix.
!     >  The first spoke (input generatrix fragment) should be at the 12 o'clock
!        position, with the centerline point first        (as for POLAR_INTERP).
!     >  The spokes may proceed either clockwise or anticlockwise ("    "    ").
!        Angles 0 through +180 degrees in the one-line file mentioned below
!        produce the right-hand half of a whole body or forebody surface.
!     >  If the generatrix represents the aft body (only), the outputs will
!        start with the centerline point.  I.e., the input order along the
!        generatrix may be reversed if necessary.  Positive angles produce the
!        right hand half as for the forebody.
!     >  The spokes are traditionally NOT uniformly spaced.  Therefore, an
!        angular distribution defining the further spokes is read as input.
!     >  If spokes are really desired on the left of a body, enter negative
!        angles.
!
!  Inputs:
!
!     A few prompts suffice for file names and padded mesh dimensions.
!
!     The generatrix should look something like this (read to EOF):
!
!        x1  z1  [0.]  or  x1  0.  z1  ! 2 or 3 columns; 
!        x2  z2  [0.]      x2  0.  z2  ! if 3 columns, the column with smallest
!        x3  z3  [0.]      x3  0.  z3  ! data range is suppressed and the other
!         :   :  [ :]       :   :   :  ! two are treated as x and y.
!
!     If the first point is not at y = 0, it is reversed in place to make the
!     first point the centerline point.  If the resulting first point is still
!     not at y = 0, this is a fatal error because the spokes produced here are
!     supposed to have a common center point (that is suppressed in all but the
!     first spoke to avoid unnecessary NEQAIR runs; POLAR_INTERP handles the
!     suppressed points).
!
!     The angular spoke coordinates should appear in a one-line file something
!     like this (read as reals):
!
!        0 5 10 15 25 35 ... 155 165 170 175 180
!
!     See assumptions above for the meaning of positive angles here, in degrees.
!
!  Output Format:
!
!     Achieving compatibility with the PREPARE_NEQAIR_DATA scheme that will
!     probably be employed prior to the POLAR_INTERP step is a little awkward.
!     The solution is to number the points of all spokes from 1 through N where
!     N is the total number of points from all spokes.  The main output of
!     'x  y  z' coordinates produced here then serves as the body point list
!     expected by PREPARE_NEQAIR_DATA.  Then NEQAIR results will appear in
!     directories /LINE-1, /LINE-2, ..., /LINE-N.
!
!     Using the command "grep Total LINE-*/neqair.out > total.txt" (say) will
!     produce a total.txt file that starts like this:
!
!        LINE-1/neqair.out:Total radiative heating ... =     2.468401 W/cm2
!        LINE-11/neqair.out:Total radiative heating ... =     3.170634 W/cm2
!        LINE-12/...
!        ...
!        LINE-19/...
!        LINE-2/...
!        LINE-20/...
!
!     In order to obtain the desired order, LINE-1, LINE-2, ..., another utility
!     by the present author, named SORT_ROWS, can be used with 'LINE-' as the
!     prefix prompted for.  This will order the heat fluxes to match the list
!     of body point coordinates (for all spokes) written here.  Then some other
!     utility such as COLUMNEDIT (same author) can append the heat fluxes to the
!     table of coordinates.
!
!     This program also writes the POLAR_INTERP-type header that needs to be
!     inserted at the top of the combined "x y z f" table after the function
!     data (NEQAIR or maybe FIAT) have been appended to the body point coords.
!     generated here.  The file names are hard-coded as 'header.dat' and
!     'body-points.dat'.
!
!     The gyrations between running NEQAIR and running POLAR_INTERP may easily
!     be scripted.
!
!     Advantage is taken of POLAR_INTERP's option to avoid repeating NEQAIR runs
!     at the origin of all spokes.
!
!  History:
!
!     08/09/2017  D.A.Saunders  Initial design.
!     08/18/2017    "      "    Initial coding and testing after a hiatus.
!     08/22/2017    "      "    Diagnostics on input data are now written.
!                               Subroutine RVERSE has 3 arguments, not 2!
!     04/04/2019    "      "    Application to a forebody with a 3-column
!                               generatrix uncovered glitches that should
!                               have been found with more thorough testing.
!                               For the forebody/right-half case, positive
!                               input angles actually produced spokes on the
!                               left half.  Positive input angles need to be
!                               negated to produce the spokes on the right
!                               half for both forebody and aft body.  If
!                               negative angles are input, they are still
!                               negated and should produce left-half spokes.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use table_type_module    ! Both are in table_io.f90
   use table_io

   implicit none

!  Constants:

   integer, parameter :: &
      lunin    = 1,      &  ! For the input generatrix
      lunout   = 2,      &  ! For the generated spokes and for header data
      lunkbd   = 5,      &  ! Keyboard inputs
      luncrt   = 6,      &  ! Screen outputs
      mxangles = 91         ! Most azimuthal angles permitted (overkill)

!  Variables:

   integer :: &
      ier, ncol, nangles, ngen
   real :: &
      angles(mxangles)
   real, allocatable, dimension (:) :: &
      xgen, ygen, zgen
   real, allocatable, dimension (:,:) :: &
      xspoke, yspoke, zspoke
   type (table_type) :: &
      generatrix

!  Execution:

!  Read the [partial] generatrix as 2 or 3 columns starting and/or ending at a
!  point on Ox, and prompt for other controls:

   call read_generatrix_etc ()
   if (ier /= 0) go to 99

!  Generate the indicated spokes:

   call construct_spokes ()
   if (ier /= 0) go to 99

!  Save results suited to POLAR_INTERP:

   call save_results ()

99 continue

   contains  ! Local procedures for program Radio_Blackout

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_generatrix_etc ()  ! Read 2|3 columns & other controls.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, parameter :: eps = 1.e-5  ! For ignoring a spanwise coord. column
      integer :: i, icol
      real    :: vmax, vmin
      logical :: cr, eof

      cr = .true.
      do while (cr)
         call reads (luncrt, 'Input generatrix (whole|partial; 2|3 columns: ', &
                     lunkbd, generatrix%filename, cr, eof)
      end do
      if (eof) go to 99

      call table_io_read_real (lunin, generatrix, ier)
      if (ier /= 0) go to 99

      ngen = generatrix%nrows
      ncol = generatrix%ncols
      icol = 2

      write (luncrt, '(a, 2i5)') '# generatrix rows and columns:', ngen, ncol

      if (ncol == 3) then  ! Identify 2-space x & y columns
         vmax = maxval (generatrix%values(icol,:))
         vmin = minval (generatrix%values(icol,:))
         if (vmax - vmin > eps) icol = 3  ! Presumably the column to ignore
         write (luncrt, '(a, i2)') 'Ignoring column', icol
         icol = 5 - icol  ! The column to use
      else if (ncol /= 2) then
         write (luncrt, '(a, i4)') &
            '*** Two or three columns are expected; number found:', ncol
         ier = 1
         go to 99
      end if

!     Working with conventional 1-D arrays is more convenient from now on:

      allocate (xgen(ngen), ygen(ngen), zgen(ngen))

      xgen(:) = generatrix%values(1,:)
      ygen(:) = 0.
      zgen(:) = generatrix%values(icol,:)

!     Read the list of azimuthal coordinates (0 - 180, probably):

      nangles = mxangles
      call rdreals (luncrt, '$Enter azimuthal angles in 0:180 on this line: ', &
                    lunkbd, nangles, angles)
      ier = 0
      if (nangles <= 0) then
         write (luncrt, '(a, i4)') 'Bad return from rdreals. nangles:', nangles
         ier = 1
      else
         write (luncrt, '(a, i4)') '# angles read:', nangles
         write (luncrt, '(i4, f8.2)') (i, angles(i), i = 1, nangles)
      end if

99    return

      end subroutine read_generatrix_etc

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine construct_spokes ()  ! ... by rotating the generatrix.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, parameter :: one = 1.0, zero = 0.0
      integer :: j

!     Deal with duplicate points on the axis later, when we write the results,
!     but reverse the generatrix order if the last point is on the axis, as it
!     presumably represents an aft body.

      ier = 0
      if (zgen(1) /= zero) then
         if (zgen(ngen) == zero) then
            write (luncrt, '(a)') &
               'Reversing generatrix order for compatibility with POLAR_INTERP.'
            call rverse (ngen, xgen, xgen)
            call rverse (ngen, ygen, ygen)
            call rverse (ngen, zgen, zgen)
         else
            write (luncrt, '(a)') &
               'Neither end point of the generatrix is at y = z = 0.', &
               'This violates a POLAR_INTERP assumption.'
            ier = 1
            go to 99
         end if
      end if

      if (angles(2) > zero) then
         write (luncrt, '(a)') &
            'Positive angles assumed to mean right-side spokes desired.'
      else
         write (luncrt, '(a)') &
            'Negative angles assumed to mean left-side spokes desired.'
      end if

      angles(1:nangles) = -angles(1:nangles) ! Because of ROTATE3D

      allocate (xspoke(ngen,nangles), yspoke(ngen,nangles), &
                zspoke(ngen,nangles))

      do j = 1, nangles
         xspoke(:,j) = xgen(:)
         yspoke(:,j) = ygen(:)
         zspoke(:,j) = zgen(:)
         call rotate3d (ngen, xspoke(:,j), yspoke(:,j), zspoke(:,j), &
                        angles(j), zero, zero, zero, one, zero, zero)
      end do

99    return

      end subroutine construct_spokes

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine save_results ()  ! I.e., x y z columns suited to POLAR_INTERP.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, j

      open  (lunout, file='header.dat', status='unknown')
      write (lunout, '(2i4)') nangles, 1   ! # spokes & # fns. for POLAR_INTERP
      write (lunout, '(i3)') &
         ngen, &        ! Spoke 1 has the centerline point and is at 12 o'clock
         (ngen - 1, j = 2, nangles)
      close (lunout)

!     The POLAR_INTERP format just lists all points of all spokes contiguously.

      open  (lunout, file='body-points.dat', status='unknown')
      write (lunout, '(3es16.8)') &
          (xspoke(i,1), yspoke(i,1), zspoke(i,1), i = 1, ngen), &
         ((xspoke(i,j), yspoke(i,j), zspoke(i,j), i = 2, ngen), j = 2, nangles)
      close (lunout)

      write (luncrt, '(a)') 'Spoked body points are here in body-points.dat.'

      end subroutine save_results

   end program capsule_spokes
