!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine curvdis (ndat, xdat, ydat, n, power, ismooth, lunout, &
                          arcs_only, x, y, ier)
!
!     Acronym:  CURVature-based DIStribution (2-space)
!               ----            ---
!     Purpose:
!
!        CURVDIS redistributes points along a curve in 2-space which is
!     represented by a discrete dataset, not necessarily monotonic, such
!     that the local spacing is inversely proportional to the local
!     curvature (with adjustments to handle the low curvature regions).
!     CURVDIS is intended for geometric applications where the coord-
!     inates have comparable units, although it may be used for non-
!     geometric purposes as well (e.g. plotting applications).  In both
!     cases, if the units of the curve coordinates are large (producing
!     small values of curvatures), the coordinates should be normalized
!     to the range [0,1].  For geometric curves, only one of the coord-
!     inate axes should be scaled to [0,1] and the other axis should be
!     scaled such that the true geometric shape is retained.
!
!     Method:
!
!        Given the availability of subroutine ARBDIS for translating
!     arbitrary "shape" functions into 1-D distributions whose spacing
!     at any abscissa is proportional to the shape function at that
!     abscissa, the problem then is to set up a shape function related
!     to curvature; ARBDIS does the rest.  (Abscissas in this case are
!     measures of arc length along the curve, and the shape function is
!     defined at each of the original data points.)
!        The subroutine is passed NDAT coordinates (XDAT, YDAT), of the
!     curve, for which it computes the cumulative arc lengths (or more
!     precisely, the cumulative chord lengths), s(n).  Next, FD12K is
!     called twice to compute the 2nd derivatives (d2x/ds2 and d2y/ds2)
!     in order to approximate the local curvatures' magnitudes, k(s),
!     which are given by:
!
!        (1)   k(s) = sqrt((d2x/ds2)**2 + (d2y/ds2)**2)
!
!        A direct inverse proportionality would yield a shape function
!     with smaller spacings in regions of larger curvature.  However,
!     for straight segments of the curve where k(s) = 0., a simple
!     inverse proportionality such as shape(s) = 1./k(s) would produce
!     an infinite local spacing.  Therefore, some arbitrary constant
!     (say 1.) should be added to k(s) in the denominator.
!        Clustering may be controlled by raising this fraction to some
!     exponent in the range [0,1].  0 would produce shape(s) = 1., or
!     uniform spacing; 1 would yield "direct" inverse proportionality
!     spacing; and an exponent in between (0,1) would produce spacings
!     between these two extremes.  Thus, our spacing function becomes:
!
!         (2)   shape(s) = [1./(k(s) + 1.)]**exponent
!
!        After the relative spacings (defining the shape function) are
!     calculated from eqn. (2), an option is provided to smooth the
!     shape function to reduce abrupt changes observed in sphere/cone
!     capsule defining sections.  The ARBDIS utility is then called to
!     compute redistributed arc lengths for the specified number N of
!     output points.  A further option is provided to smooth those
!     arc-lengths (with index-based abscissas) as also found desirable
!     for sphere/cone applications.  LCSFIT then serves to evaluate the
!     new X-Y coordinates corresponding to the new arc lengths.
!
!        This version also has the option to constrain cell growth rates
!     by smoothing them until they lie within [0.8, 1.2].  In order not
!     to change the calling sequence, ISMOOTH < 0 invokes that option,
!     as it also does for the vertex-capturing option.
!
!        NO: Constraining growth rates is a bad idea, defeating the
!     purpose of curvature-based spacing.  ISMOOTH < 0 is now used
!     instead to specify special treatment of sharp corners.
!
!        LATER:  ISMOOTH < 0 also controls broadening of high-curvature
!     regions onto flat/low curvature segments as found desirable in
!     flow calculations for entry capsules at the shoulder.
!
!     Work-space Usage:
!
!        (1)  WORK(1:NDAT)                 Input curve chord lengths
!        (2)  WORK(NDAT+1:2*NDAT)          FP values for FD12K calls
!        (3)  WORK(2*NDAT+1:3*NDAT)        d2x/ds2 values
!        (4)  WORK(3*NDAT+1:4*NDAT)        d2y/ds2 values
!        (5)  WORK(NDAT+1:2*NDAT)          xdat-ydat relative spacings
!                                          defining the shape function
!        (6)  WORK(2*NDAT+1:2*NDAT+N)      New curve chord lengths
!        (7)  WORK(2*NDAT+N+1:5*NDAT+6*N)  Work area used by ARBDIS
!        (8)  WORK(2*NDAT+N+1:2*NDAT+2*N)  Reused for index-based abscissas
!                                           if ARBDIS results are smoothed
!     Procedures:
!
!        ARBDIS           Distributes points in an interval according to an
!                         arbitrary shape function
!        FD12K            Calculates derivatives by finite differencing
!        CHORDS2D         Computes the cumulative chord lengths of a curve
!        DETECT_VERTICES  Identifies curvature spikes likely due to sharp
!                         vertices where some clustering is desirable
!        LCSFIT           Interpolates x and y to the redistributed arc
!                         lengths along the parameterized curve
!        SMOOTH1D_NITERS  Explicit smoothing utility with iteration control
!        VERTEX_CURVATURE Broadens and moderates the height of curvature spikes
!                         to enable some clustering towards sharp vertices
!     Error Handling:
!
!        See IER description below.  It is the user's responsibility
!     to check IER on return from CURVDIS.  Error messages are sent to |LUNOUT|.
!
!     History:
!        08/12/88    D.A.Saunders   Initial design.
!        09/09/88    B.A.Nishida    Initial implementation.
!        06/05/89    DAS            Refined the description a little.
!        02/15/96     "             Cubic interpolation rather than linear.
!        12/01/10    DAS, ERC, Inc. Inserted explicit smoothing of the shape
!                                   function, as suggested by the defining
!                                   section of a sphere/cone capsule forebody.
!        12/03/10     "    "        Index-based smoothing of the redistributed
!                                   arc-lengths helps for cases with extreme
!                                   jumps in curvature, such as a sphere/cone.
!                                   New argument ISMOOTH allows suppression of
!                                   either or both smoothings.
!                                   The work-space is no longer an argument, as
!                                   an automatic array is now preferable.
!        12/06/10     "    "        Option to return only the redistributed arc-
!                                   length distribution (as X(1:N)), as is more
!                                   convenient for the HEAT_SHIELD program.
!        08/01/12     "    "        Very large grids can mean the original limit
!                                   of 30 smoothing iterations may not be enough
!                                   yet n/8 may still be too many.  Compromise
!                                   by using n/10 with a limit of 100.
!        10/23/13     "    "        Now that CURVDIS2 has been introduced on top
!                                   of CURVDIS to normalize the data before any
!                                   curvature calculations, less smoothing seems
!                                   to be desirable.  (A loop in CURVDIS2 that
!                                   lowers the exponent by 0.1 if CURVDIS does
!                                   not converge also helps adhering to the
!                                   curvature-based shape function as much as
!                                   possible without failing.)
!        11/01/13     "    "        Option to try and resolve apparent sharp
!                                   vertices in the curve (ISMOOTH < 0).
!        11/04/13     "    "        NV = 6 needs to scale with NDAT.
!        11/07/13     "    "        Option to smooth cell growth rates until
!                                   they're within [0.8, 1.2] (ISMOOTH < 0 too).
!        11/08/13     "    "        NO!  Such smoothing can lose the connection
!                                   between the final arc length distribution
!                                   and the original data curvature distribu-
!                                   tion!  Suppress this misguided addition.
!        04/18/16    DAS, AMA, Inc. Poor resolution of a small shoulder radius
!                                   by CAPSULE_GRID prompted reverting to more
!                                   smoothing iterations that should allow a
!                                   higher exponent to converge.
!        06/27/16     "    "        The 04/18/16 change turned out to be overly
!                                   influenced by a generatrix from a CAD system
!                                   that had unusually dense resolution on the
!                                   shoulder.  There is no way of safeguarding
!                                   any smoothing algorithm from wildly varying
!                                   numbers of data points!
!        06/29/16     "    "        The latest attempt to improve capsule grid
!                                   point distributions in the presence of flat
!                                   segments (sphere/cone) seeks to spread the
!                                   redistributed points further on to the low
!                                   curvature regions as found desirable for
!                                   good resolution of spikes in surface heat
!                                   flux at the shoulder.  New utilities
!                                   DETECT_FLATNESS and HANDLE_FLATNESS are
!                                   activated as for the earlier attempt at
!                                   treating sharp corners (ISMOOTH < 0).
!        05/11/18     "    "        Converted to Fortran 90 for inclusion in
!                                   CAPSULE_GRID source code to allow a verbose
!                                   mode. Handling of geometry vertices and flat
!                                   regions requires enough heuristics that we
!                                   need to see quantities like the normalized
!                                   curvature being worked with here to be able
!                                   to tell if those heuristics need adjusting.
!        09/25/19     "    "        Puzzling results observed with CAPSULE_GRID
!                                   application to a capsule on a sting case
!                                   (analytic, with corner rounding) were traced
!                                   to the earlier DETECT_VERTICES/VERTEX_CURV-
!                                   ATURE pair, which are now suppressed.  The
!                                   later blending on to low-curvature regions
!                                   should achieve the intended effect at sharp
!                                   corners as well.
!
!     Authors:  David Saunders, Brian Nishida, Sterling Software
!               NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:
!     ----------

      integer, intent (in)  :: ndat       ! Number of input curve points
      real,    intent (in)  :: xdat(ndat) ! X coordinates of input curve
      real,    intent (in)  :: ydat(ndat) ! Y coordinates of input curve
      integer, intent (in)  :: n          ! Number of output curve points
      real,    intent (in)  :: power      ! Exponent for spacing function;
                                          ! controls the point clustering.
                                          ! POWER must be in [0., 1.];
                                          ! 0.0 produces uniform spacing;
                                          ! 1.0 maximizes the curvature effect;
                                          ! 0.5 is suggested, but start with 1.0
                                          ! if CURVDIS2 is being employed, as is
                                          ! now recommended
      integer, intent (in)  :: ismooth    ! 0 => no smoothing;
                                          ! 1 => smooth crv-shape function only;
                                          ! 2 => smooth redistributed arcs only;
                                          ! 3 => perform both smoothings;
                                          ! use ISMOOTH = -1, -2, or -3 to allow
                                          ! automated treatment of apparently
                                          ! sharp vertices which otherwise are
                                          ! missed (1-pt. curvature spikes are
                                          ! essentially invisible to the scheme)
                                          ! (NO: this can interfere with the
                                          ! later blending onto low-curvature
                                          ! regions, which should handle sharp
                                          ! corners too;)
                                          ! and also to constrain cell growth
                                          ! rates via further smoothing of them.
                                          ! (NO: the latter is a BAD idea; see
                                          ! the History above).
      integer, intent (in)  :: lunout     ! Logical unit for showing convergence
                                          ! history; LUNOUT < 0 suppresses it
      logical, intent (in)  :: arcs_only  ! T => skip calculations of X & Y(1:N)
                                          ! but return revised arcs as X(1:N)
      real,    intent (out) :: x(n), y(n) ! X & Y coordinates of output curve,
                                          ! unless ARCS_ONLY = T, in which case
                                          ! just the arc-lengths are returned as
                                          ! X(1:N)
      integer, intent (out) :: ier        ! 0 => no errors;
                                          ! 1 => POWER is out of range;
                                          ! 2 => failure in ARBDIS

!------------------------------------------------------------------------------

!     Local Constants:
!     ----------------

      integer,   parameter :: nflmax   = 12,    &! Limit on # flat regions
                              nvertmax = 12,    &! Limit on # vertices handled
                              nvmin    = 5,     &! Limits on # vertex neighbors
                              nvmax    = 15      ! each side of curvature spikes
      real,      parameter :: alpha    = 1.e+0, &! In case curvature = 0.
                              exponent = 3.,    &! Controls shape fun. widths
                              fraction = 0.10,  &! Applied to spike heights
                              fractn   = 0.035, &! Applied to NDAT for NV
                              one      = 1.e+0, &
                              half     = 0.5+0, &
                              zero     = 0.e+0
      logical,   parameter :: false    = .false., &
                              true     = .true.
      character, parameter :: linear*1 = 'L',   &
                              tight*1  = 'M',   &
                              name*9   = 'CURVDIS: '

!     Local Variables:
!     ----------------

      integer :: i, ismth, iv, lunerr, nfl, niters, nv, nvertex
      integer :: ifl(nflmax), ifr(nflmax), ivertex(nvertmax)
      real    :: dsq, dsqmin, oneovern, stotal, xd, yd
      real    :: work(5*ndat + 6*n)
      logical :: constrain_growth_rates, treat_vertices, verbose

!     Procedures:
!     -----------

      external :: arbdis, chords2d, fd12k, lcsfit, smooth1d_niters
      external :: detect_vertices, vertex_curvature
      external :: detect_flatness, handle_flatness

!     Execution:
!     ----------

      treat_vertices = ismooth < 0             ! Retrofited option, now
                                               ! invoking better blending of
                                               ! low-curvature regions likewise
!!!   constrain_growth_rates = treat_vertices  ! Kludge for now
      constrain_growth_rates = .false.         ! See History above
      nvertex = 0                              ! See LCSFIT test near the end

      ismth   = abs (ismooth)
      lunerr  = abs (lunout)
      verbose = lunout > 0

      if (power < zero .or. power > one) then
         write (lunerr, '(/, 2a, es12.4)') &
            name, 'Bad POWER input outside [0, 1]:', power
         ier = 1
         go to 999
      end if

!     Compute the cumulative chord lengths of the input data:

      call chords2d (ndat, xdat, ydat, false, stotal, work(1))

!     Calculate the 2nd derivatives, d2x/ds2 and d2y/ds2:

      call fd12k (ndat, work(1), xdat, work(ndat+1), work(2*ndat+1), &
                  work(2*ndat+1))
      call fd12k (ndat, work(1), ydat, work(ndat+1), work(3*ndat+1), &
                  work(3*ndat+1))

!     Local curvature magnitudes:

      do i = 1, ndat
         work(ndat+i) = sqrt (work(2*ndat+i)**2 + work(3*ndat+i)**2)
      end do

      if (verbose) then
         write (lunout, '(a)') ' Arc length and unadjusted curvature:'
         write (lunout, '(i4, 2es14.6)') (i, work(i), work(ndat+i), i = 1, ndat)
      end if

!     Attempt to cluster towards any vertices that will otherwise be missed?
!     Broadening of high curvature onto flat/low-curvature regions is also
!     automated here now, as needed for atmospheric entry capsule grids.

      if (treat_vertices) then

         call detect_vertices (lunout, ndat, work(ndat+1), nvertmax, ivertex, &
                               nvertex)
         write (lunerr, '(/, (a))') &
            '*** Suppressing handling of sharp corners in favor of', &
            '    flatness blending only, to avoid likely interference. ***'
         nvertex = 0

         if (nvertex > 0) then  ! Adjust curvature in-place; x/y data untouched
            nv = min (nvmax, max (nvmin, nint (fractn * real (ndat))))
            if (verbose) write (lunout, '(a, i3, a, i3)') &
               '# vertices found:', nvertex, &
               '   # broadening neighbors either side:', nv

            do iv = 1, nvertex
               call vertex_curvature (ndat, work(1), work(ndat+1), &
                                      ivertex(iv), nv, exponent, fraction)
            end do
         end if

!        Broaden abrupt curvature changes onto low-curvature segments:

         call detect_flatness (lunout, ndat, work(ndat+1), &
                               nflmax, ifl, ifr, nfl)
         if (nfl > 0) then
            if (verbose) then
               write (lunout, '(a, i3)') &
                  '# flat regions detected:', nfl, 'ifl, ifr:'
               write (lunout, '(2i6)') (ifl(i), ifr(i), i = 1, nfl)
            end if

            call handle_flatness (lunout, ndat, work(1), work(ndat+1), &
                                  nfl, ifl, ifr)
         end if
      end if

!     Moderate the (possibly adjusted) curvatures:

      if (power == half) then  ! Avoid logarithms
         do i = 1, ndat
            work(ndat+i) = one / sqrt (work(ndat+i) + alpha)
         end do
      else
         do i = 1, ndat
            work(ndat+i) = (work(ndat+i) + alpha)**(-power)
         end do
      end if

      if (verbose) then
         write (lunout, '(a)') ' Unsmoothed curvature-related shape function'
         write (lunout, '(2es14.6)') (work(i), work(ndat+i), i = 1, ndat)
      end if

      if (ismth == 1 .or. ismth == 3) then

!        Smooth the shape function further with an explicit method:

!!!!     niters = max (5, min (ndat/10, 30))
!!!!     niters = max (4, min (ndat/40, 15))
         niters = max (9, min (ndat/20, 20))

         call smooth1d_niters (niters, 1, ndat, work(1), work(ndat+1))

         if (verbose) then
            write (lunout, '(a, i5)') ' NITERS for shape function:', niters
            write (lunout, '(a)') ' Smoothed curv.-based shape function'
            write (lunout, '(2es14.6)') (work(i), work(ndat+i), i = 1, ndat)
         end if

      end if

!     Redistribute the arc lengths based on the shape function:

      stotal = work(ndat)

      call arbdis (n, zero, stotal, ndat, work(1), work(ndat+1), 'a', &
                   lunout, work(2*ndat+n+1), work(2*ndat+1), ier)
      if (ier /= 0) then
         write (lunerr, '(/, 2a, i4)') name, 'ARBDIS ier:', ier
         ier = 2
         go to 999
      end if

      if (ismth >= 2) then  ! Smooth the redistributed arc-lengths

         oneovern = one / real (n)  ! Index-based abscissas in ARBDIS workspace
         do i = 1, n
            work(2*ndat+n+i) = real (i) * oneovern
         end do

!!!!     niters = max (5, min (n/10, 30))
         niters = max (5, min (n/20, 20))
         if (verbose) then
            write (lunout, '(a, i5)') &
               'NITERS for redistributed arc lengths:', niters
            write (lunout, '(a)') ' Unsmoothed redistributed arc lengths'
            write (lunout, '(2es14.6)') &
               (work(2*ndat+n+i), work(2*ndat+i), i = 1, n)
         end if

         call smooth1d_niters (niters, 1, n, work(2*ndat+n+1), work(2*ndat+1))

         if (verbose) then
            write (lunout, '(a)') ' Smoothed redistributed arc lengths'
            write (lunout, '(2es14.6)') &
               (work(2*ndat+n+i), work(2*ndat+i), i = 1, n)
         end if

      end if

!     Redistribute the curve coordinates, unless only revised arcs are wanted:

      if (arcs_only) then
         x(1:n) = work(2*ndat+1:2*ndat+n)
      else

         if (constrain_growth_rates) then   ! NOT A GOOD IDEA -- suppressed now
            call adjust_growths (n, work(2*ndat+1))         ! Internal procedure
         end if

         call lcsfit (ndat, work(1), xdat, true, tight, n, &
                      work(2*ndat+1), x, work(2*ndat+1+n))  ! Unused derivatives

         call lcsfit (ndat, work(1), ydat, true, tight, n, &
                      work(2*ndat+1), y, work(2*ndat+1+n))  !    "     "

         if (nvertex > 0) then  ! Even monotonic cubics are dubious at a vertex

            do iv = 1, nvertex  ! Messy!  Find grid pt. nearest to data vertex
               i = ivertex(iv)  ! Data index
               xd = xdat(i)
               yd = ydat(i)
               dsqmin = 1.e+6

               do i = 1, n
                  dsq = (x(i) - xd)**2 + (y(i) - yd)**2
                  if (dsq < dsqmin) then
                      dsqmin = dsq
                      ivertex(iv) = i  ! Now a redistributed index
                  end if
               end do

               i = ivertex(iv) - 1  ! Reinterpolate X, Y linearly at I-1, I, I+1

               call lcsfit (ndat, work(1), xdat, true, linear, 3, &
                            work(2*ndat+i), x(i), work(2*ndat+1+n))
               call lcsfit (ndat, work(1), ydat, true, linear, 3, &
                            work(2*ndat+i), y(i), work(2*ndat+1+n))
            end do
         end if

      end if

!     Error handling:
!     ---------------

  999 return

!     Internal procedure for CURVDIS if ISMOOTH < 0:

      contains  ! No longer recommended; retained for posterity

!        -----------------------------------------------------------------------
!
         subroutine adjust_growths (n, s)  ! Arguments simplify nomenclature
!
!        Evaluate cell growth rates 2:N-1.  If any lie outside [GMIN, GMAX],
!        apply one smoothing iteration to them and calculate the corresponding
!        adjusted arc lengths.  Repeat as often as needed, up to a limit.
!
!        Later:  This is a bad idea because it can lose the connection between
!        the redistributed arcs and the original data curvature distribution.
!        It is suppressed now, but left in-line just in case it's of interest.
!        -----------------------------------------------------------------------

!        Arguments:

         integer, intent (in)    :: n      ! # points in arc-length distribution
         real,    intent (inout) :: s(n)   ! Arc lengths adjusted in-place

!        Local constants:

         integer, parameter :: itmax = 10  ! Limit on adjustment iterations
         real, parameter    :: gmax = 1.2,&! Cell growth rate limits
                               gmin = 0.8

!        Local variables:

         integer :: i, iter, nbad
         real    :: b(n), c(n), g(n), r(n)

!        Execution:

         do iter = 1, itmax  ! Until all growth rates are acceptable

            nbad = 0
            g(1) = one
            do i = 2, n - 1
               g(i) = (s(i+1) - s(i)) / (s(i) - s(i-1))
               if (g(i) < gmin .or. g(i) > gmax) nbad = nbad + 1
            end do
            g(n) = one
!!!!        write (lunout, '(a, 2i4)') ' Itn. & # bad growth rates:', iter, nbad

!           Note that if growth rates are recovered later from the
!           redistributed (X, Y)s, they won't be identical to what we see here.

            if (nbad == 0) exit

            call smooth1d_niters (1, 1, n, s, g)

!           Back out the arc lengths that correspond to the adjusted growths:

            do i = 2, n - 1
!!!            a(i) =  g(i)          ! Subdiagonal is just G(2:N-1)
               b(i) = -(one + g(i))  ! Diagonal
               c(i) =  one           ! Superdiagonal
               r(i) =  zero          ! RHS, except for first & last elements
            end do
            r(2)   = -s(1)*g(2)
            r(n-1) = -s(n)

            call trdiag (g(2), b(2), c(2), r(2), s(2), n-2)

         end do

         end subroutine adjust_growths

      end subroutine curvdis
