!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program polar_interp
!
!  Description:
!
!     From surface data representing "spokes" of a symmetric body ordered by
!     azimuthal angle, interpolate a single-patch surface grid (and optional
!     function data) with specified dimensions.  The method can apply to a body
!     with off-center nose as well.  The interpolated grid has uniform spacing
!     (at least in the initial implementation).
!
!     The intent is to pad coarse input data via nonlinear interpolation in the
!     two surface directions using a 1-D spline method, by taking advantage of
!     the spoke structure, so that methods appropriate for high-resolution data
!     can then be applied.  The initial application is to (expensive) radiation
!     calculations on a coarse polar grid that need to be interpolated to a
!     dense CFD mesh.  Direct use of SURFACE_INTERP produces faceting because
!     the data cells are so big and the method is a one-cell method.
!
!     The azimuthal interpolations are performed versus clock angle, and the
!     angles of the data spokes are deduced from the given (x,y,z) coordinates.
!
!     Originally, only half the body was treated but this version has the option
!     to treat a full body. In the full body case, the first spoke is replicated
!     automatically, and periodic interpolations are employed in the azimuthal
!     direction.  In the half body case, the second and second-last spokes are
!     reflected automatically so that the interpolations match what would be
!     obtained if the missing half body were present (by reflection) explicitly.
!     If the first and/or last spokes do not appear to be on the center line,
!     the dataset is still interpolated to a structured surface mesh with no
!     symmetry assumptions.
!
!     A few prompts control the run (or a short control file on standard input).
!
!     This version provides for replicating the center point as makes sense for
!     expensive radiative heating data: one spoke may be indicated as containing
!     the center location as its first point; unless 0 is entered to suppress
!     this option, each other spoke then has that point inserted by this program
!     as its new first point.  Standard input entries:
!
!        xxx.dat             ! Input dataset name
!        ni nj               ! Output surface grid dimensions
!        yyy.dat             ! Output Tecplot ASCII dataset
!        1                   ! Spoke containing center point, or 0
!
!  Coordinate System:
!
!     Ox is assumed to be parallel to the symmetry axis, and right handedness
!     is assumed.  Oz is "up".
!
!  Input Data Format:
!
!     For historical reasons (NEQAIR practices), the input format is PLOT3D-like
!     (header records at the top), but in "point" rather than "block" order,
!     with at least one spoke beginning with the center point that is assumed
!     to be common to all spokes (see above option to replicate this point).
!     The first spoke should be at the 12 o'clock position (with the nose point
!     first).  The following spokes may proceed either clockwise or anticlock-
!     wise.
!
!     11 3                 [Number of spokes and number of functions per spoke]
!     15                   [Number of points on spoke 1]
!     15                   [  "   "   "   "   "   "   2]
!     :
!     :
!     9                    [  "   "   "   "   "   "  11]
!     x y z f1 f2 f3       [Coordinates and functions for point 1 of spoke 1]
!     x y z f1 f2 f3
!     : : : : : : :
!     : : : : : : :
!     x y z f1 f2 f3       [  "   "   "   "   "   "   "   "   " 9  "   "  11]
!
!  Output Data Format:
!
!     Single-zone Tecplot ASCII dataset, BLOCK order, suited to SURFACE_INTERP.
!
!  Option To Help Generate Target Spokes:
!
!     Constructing spokes on an "ellipsled" slender-body vehicle prompted this
!     option to read surface grid patch indices in the spoked form above, with
!     x and y replaced by i and j, and no z or function data.  The output is
!     then in the above spoked form (x, y, z only) as needed to define body
!     points for radiative heating calculations.  No interpolations are done.
!
!     Alternative Input Data Format (to help prepare the normal input data):
!
!        13  ! # spokes
!        18  ! # points on spoke 1:nspokes
!        17
!        17
!        :
!        :
!        8
!        1  1 ! (i,j) defining apex point 1 of spoke 1
!        8  1
!        16 1
!        :  :
!        :  :
!        8  3
!        16 3
!        :  :
!        :  :
!        8  65
!        17 65
!        :  :
!        :  :
!        137 65
!
!     If the first token on every line of the dataset prompted for initially
!     is an integer (no decimal point), this retrofitted option is assumed,
!     and a new prompt for the single-patch surface grid to which the indices
!     refer is issued so that indices can be converted to (x,y,z) coordinates.
!     This surface grid should be in PLOT3D form, ASCII or unformatted.
!
!  History:
!
!     08/29/08  D. Saunders  Initial implementation as part of an alternative
!                            to giving NEQAIR a (coarse) surface triangulation.
!                            Preserving some structure and padding should allow
!                            fewer coarse surface points than the unstructured
!                            approach to putting points where they're needed.
!
!     10/06/08   "      "    Added the option to replicate the first (center)
!                            point of one spoke as the first point for all other
!                            spokes, as makes sense for radiative heating data.
!
!     05/19/09   "      "    Retrofitted an option to read indices defining the
!                            spoked data points needed to set up body point
!                            coordinates for radiative heating calculations.
!                            A PLOT3D-type surface grid must accompany the input
!                            data file if it apparently contains (i,j) indices
!                            that define the spokes rather than (x,y,z,f) data.
!                            Only a single-patch grid appears necessary, as
!                            produced by HEAT_SHIELD (but with y and z swapped
!                            to be compatible with DPLR).
!
!     06/15/12   "      "    Jeff Brown had wind tunnel data on both halves of
!                            an inflatable decelerator, prompting the option to
!                            handle more than just half-body datasets.  The
!                            spoke order may be clockwise or anticlockwise but
!                            spoke 1 should be at 12 o'clock for half and whole
!                            body cases.  (Parts of bodies may still be OK no
!                            matter where they start, with no symmetry assump-
!                            tions made.)
!
!     06/18/12   "      "    The full-body case can use LCSFIT's cyclic option,
!                            at least for the geometry.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!                           now ERC, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      lunkbd  = 5,       &
      luncrt  = 6,       &
      lun_in  = 1,       &
      lun_out = 2

   logical, parameter ::     &
      diagnostics = .false., &
      new         = .true.,  &
      normalize   = .false.

   character, parameter :: &
      cyclic*1 = 'C',      &
      loose *1 = 'B',      &
      tight *1 = 'M'

!  Variables:

   integer :: &
      i, i1, ios, is, ispoke, j, maxnidata, mi, n, nf, ni, nj, npts, nspokes

   integer, allocatable, dimension (:) :: &
      nidata

   real :: &
      ds, dtheta, pi, r, stotal, xcenter, ycenter, zcenter

   real, allocatable, dimension (:) :: &
      derivs, odata, ogrid, sdata, sgrid, tdata, tgrid

   real, allocatable, dimension (:,:) :: &
      xdata, ydata, zdata, xgrid, ygrid, zgrid, xreg, yreg, zreg

   real, allocatable, dimension (:,:,:) :: &
      fdata, fgrid, freg

   logical :: &
      half_body, whole_body

   character :: &
      filename_in*80, filename_out*80, method*1

!  Execution:

!  Read the spoke header info, allocate storage, then read the spokes.
!  These may now be defined by indices into a spoked surface grid, though.

   ios = 1
   do while (ios /= 0)
      write (luncrt, '(a)', advance='no') '   Input dataset name: '
      read  (lunkbd, *) filename_in
      open  (lun_in, file=filename_in, status='old', iostat=ios)
   end do

!  Check for the retrofitted option to read indices.  If we return with
!  ios = 0, we assume that option wasn't detected, and proceed as originally,
!  else we're done.

   call look_for_indices_option ()  ! Internal procedure below

   if (ios /= 0) go to 99

   read (lun_in, *, iostat=ios) nspokes, nf
   if (ios /= 0) then
      write (luncrt, '(a)') '   Trouble reading nspokes and nf on line 1.'
      go to 99
   end if

   allocate (nidata(nspokes))

   maxnidata = 0
   do j = 1, nspokes
      read (lun_in, *, iostat=ios) nidata(j)  ! Avoid possible trailing comments
      if (ios /= 0) then
         write (luncrt, '(a)') '   Trouble reading # points for each spoke.'
         go to 99
      end if
      maxnidata = max (maxnidata, nidata(j))
   end do

!  Read the actual spoke data after it is known whether the center point needs
!  to be replicated or not.

   ios = 1
   do while (ios /= 0)
      write (luncrt, '(a)', advance='no') '   Output grid dimensions: '
      read  (lunkbd, *, iostat=ios) ni, nj
   end do

   ios = 1
   do while (ios /= 0)
      write (luncrt, '(a)', advance='no') '   Output dataset name: '
      read  (lunkbd, *) filename_out
      open  (lun_out, file=filename_out, status='unknown', iostat=ios)
   end do

   ios = 1;  ispoke = -1
   do while (ios /= 0 .or. ispoke < 0 .or. ispoke > nspokes)
      write (luncrt, '(a)', advance='no') &
         '   Spoke containing center point [0 => all do]: '
      read  (lunkbd, *, iostat=ios) ispoke
   end do

   if (ispoke /= 0) then
      maxnidata = maxnidata + 1  ! Be sure of room for inserting center point
      i1 = 2
      do j = 1, nspokes
         if (j /= ispoke) nidata(j) = nidata(j) + 1
      end do
   else
      i1 = 1
   end if

   allocate (xdata(maxnidata,nspokes), ydata(maxnidata,nspokes), &
             zdata(maxnidata,nspokes), fdata(nf,maxnidata,nspokes), &
             sdata(maxnidata))

   do j = 1, nspokes
      is = i1;  if (j == ispoke) is = 1
      read (lun_in, *, iostat=ios) &
         (xdata(i,j), ydata(i,j), zdata(i,j), fdata(:,i,j), i = is, nidata(j))
      if (ios /= 0) then
         write (luncrt, '(a, i4)') &
            '   Trouble reading data for spoke #', j
         go to 99
      end if
   end do

   close (lun_in)

   if (ispoke /= 0) then
      do j = 1, nspokes
         if (j /= ispoke) then
            xdata(1,j)   = xdata(1,ispoke)
            ydata(1,j)   = ydata(1,ispoke)
            zdata(1,j)   = zdata(1,ispoke)
            fdata(:,1,j) = fdata(:,1,ispoke)
         end if
      end do
   end if

   write (lun_out, '(3a)') 'TITLE = "', trim (filename_in), ', padded"'
   write (lun_out, '( a)') 'VARIABLES = "x" "y" "z"'
   write (lun_out, '(a, i3, a)') ('"function', n, '"', n = 1, nf)
   write (lun_out, '(a)') 'ZONE T = "Gridded data"'
   write (lun_out, '(a, i3, a, i3, a)')  &
      'I =', ni, ', J =', nj, ', K=1, ZONETYPE=ORDERED, DATAPACKING=BLOCK'

   xcenter = xdata(1,1)
   ycenter = ydata(1,1)
   zcenter = zdata(1,1)

!  Regularize the input spokes:

   allocate (sgrid(ni), derivs(ni), &
             xreg(ni,nspokes), yreg(ni,nspokes), zreg(ni,nspokes), &
             odata(maxnidata), ogrid(ni), freg(nf,ni,nspokes))

   do j = 1, nspokes

      mi = nidata(j)

      call chords3d (mi, xdata(1:mi,j), ydata(1:mi,j), zdata(1:mi,j), &
                     normalize, stotal, sdata(1:mi))

      ds = stotal / real (ni - 1)
      do i = 1, ni - 1
         sgrid(i) = ds * real (i - 1)
      end do
      sgrid(ni) = stotal

      call lcsfit (mi, sdata(1:mi), xdata(1:mi,j), new, loose, ni, sgrid, &
                   xreg(1,j), derivs)
      call lcsfit (mi, sdata(1:mi), ydata(1:mi,j), new, loose, ni, sgrid, &
                   yreg(1,j), derivs)
      call lcsfit (mi, sdata(1:mi), zdata(1:mi,j), new, loose, ni, sgrid, &
                   zreg(1,j), derivs)

      do n = 1, nf
         odata(1:mi) = fdata(n,1:mi,j)
         call lcsfit (mi, sdata(1:mi), odata(1:mi), new, tight, ni, sgrid, &
                      ogrid, derivs)
         freg(n,:,j) = ogrid(:)
      end do

   end do

!! write (8, '(i1)') 1
!! write (8, '(3i4)') ni, nspokes, 1
!! write (8, '(5es16.7)') xreg, yreg, zreg
!! close (8)

!! write (9, '(i1)') 1
!! write (9, '(4i4)') ni, nspokes, 1, nf
!! do n = 1, nf
!!    write (9, '(5es16.7)') freg(n,:,:)
!! end do
!! close (9)

!  Determine the azimuthal angles represented by the spokes.  They may be
!  adjusted later so as to be monotonically increasing or decreasing.
!  Recall that atan2 (y, x) is in [0, pi] for quadrants 1 & 2, and in [-pi, 0]
!  for quadrants 3 & 4.  This proves awkward unless we arrange for the angle
!  to be 0. at 6 o'clock and +/-pi at 12 o'clock.  Hence the odd minus signs.

   pi = atan2 (0., -1.)

   allocate (tdata(nspokes))

   do j = 1, nspokes
      mi = nidata(j)
      tdata(j) = atan2 (-(ydata(mi,j) - ycenter), -(zdata(mi,j) - zcenter))
   end do

   if (diagnostics) then
      write (6, '(a)') 'Unadjusted azimuthal angles (radians):'
      write (6, '(i3, es16.8)') (j, tdata(j), j = 1, nspokes)
   end if

!  Ascertain the symmetry to assume (if any), and adjust the data accordingly:

   call symmetry_adjustments ()

!  Interpolate the regularized spokes in the azimuthal direction:

   allocate (tgrid(nj))

   dtheta = (tdata(nspokes) - tdata(1)) / real (nj - 1)
   do j = 1, nj - 1
      tgrid(j) = tdata(1) + dtheta * real (j - 1)
   end do
   tgrid(nj) = tdata(nspokes)
!! write (6, '(i3, es16.8)') (j, tgrid(j), j = 1, nj)

   deallocate (odata, ogrid, derivs)

   if (half_body) then
      allocate (odata(0:nspokes+1))
      nspokes = nspokes + 2
   else
      allocate (odata(nspokes))
   end if

   allocate (ogrid(nj), derivs(nj), &
             xgrid(ni,nj), ygrid(ni,nj), zgrid(ni,nj), fgrid(nf,ni,nj))

!  Singular point at center:

   xgrid(1,:) = xdata(1,1)
   ygrid(1,:) = ydata(1,1)
   zgrid(1,:) = zdata(1,1)

   do n = 1, nf
      fgrid(n,1,:) = fdata(n,1,1)
   end do

   do i = 2, ni

      odata(:) = xreg(i,:)
      call lcsfit (nspokes, tdata, odata, new, method, nj, tgrid, ogrid, derivs)
      xgrid(i,:) = ogrid(:)

      odata(:) = yreg(i,:)
      call lcsfit (nspokes, tdata, odata, new, method, nj, tgrid, ogrid, derivs)
      ygrid(i,:) = ogrid(:)

      odata(:) = zreg(i,:)
      call lcsfit (nspokes, tdata, odata, new, method, nj, tgrid, ogrid, derivs)
      zgrid(i,:) = ogrid(:)

      do n = 1, nf
         odata(:) = freg(n,i,:)
         call lcsfit (nspokes, tdata, odata, new, tight, nj, tgrid, ogrid, &
                      derivs)
         fgrid(n,i,:) = ogrid(:)
      end do

   end do

!  Save the padded data:

   write (lun_out, '(5es16.7)') xgrid
   write (lun_out, '(5es16.7)') ygrid
   write (lun_out, '(5es16.7)') zgrid

   do n = 1, nf
      write (lun_out, '(5es16.7)') fgrid(n,:,:)
   end do

   close (lun_out)

99 continue

! *** stop ! Avoid system dependencies.

!  Internal procedures for program polar_interp:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine look_for_indices_option ()  ! See main program header

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      character, parameter :: format * 11 = 'unformatted'

!     Local variables:

      integer   :: if1, ii, jj, ni, nj
      logical   :: formatted, not_integer
      character :: token * 16
      integer, allocatable, dimension (:,:) :: idata, jdata

!     Execution:

      not_integer = .false.
      do ! Until EOF
         read (lun_in, *, iostat=ios) token
         if (ios < 0) exit
         if (index (token, '.') /= 0) not_integer = .true.
      end do

      rewind (lun_in)

      if (not_integer) then
         ios = 0
         go to 99
      end if

!     Read indices much as for the original spoked data:

      read (lun_in, *) nspokes
      write (luncrt, '(a, i4)') &
         'Spokes appear to be defined by indices. # spokes:', nspokes

      allocate (nidata(nspokes))

      maxnidata = 0
      do j = 1, nspokes
         read (lun_in, *, iostat=ios) nidata(j)  ! Avoid any trailing comments
         if (ios /= 0) then
            write (luncrt, '(a)') '   Trouble reading # points for each spoke.'
            go to 99
         end if
         maxnidata = max (maxnidata, nidata(j))
      end do

      allocate (idata(maxnidata,nspokes), jdata(maxnidata,nspokes), &
                xdata(maxnidata,nspokes), ydata(maxnidata,nspokes), &
                zdata(maxnidata,nspokes))

      do j = 1, nspokes
         read (lun_in, *, iostat=ios) (idata(i,j), jdata(i,j), i = 1, nidata(j))
         if (ios /= 0) then
            write (luncrt, '(a, i4)') 'Trouble reading indices for spoke #', j
            go to 99
         end if
      end do

      close (lun_in)

!     Get the surface grid to which the indices refer:

      write (luncrt, '(a)', advance='no') 'Surface grid for these indices: '
      read  (lunkbd, *) filename_in

      call determine_grid_form (filename_in, lun_in, formatted, ios)
      if (ios /= 0) go to 99

      if1 = 1;  if (formatted) if1 = 3
      open  (lun_in, file=filename_in, status='old', form=format(if1:11))

      if (formatted) then
         read (lun_in, *) ! nb is assumed to be 1
         read (lun_in, *) ni, nj
         allocate (xgrid(ni,nj), ygrid(ni,nj), zgrid(ni,nj))
         read (lun_in, *) xgrid, ygrid, zgrid
      else
         read (lun_in) ! nb is assumed to be 1
         read (lun_in) ni, nj
         allocate (xgrid(ni,nj), ygrid(ni,nj), zgrid(ni,nj))
         read (lun_in) xgrid, ygrid, zgrid
      end if
      close (lun_in)

      do j = 1, nspokes
         do i = 1, nidata(j)
            ii = idata(i,j)
            jj = jdata(i,j)
            xdata(i,j) = xgrid(ii,jj)
            ydata(i,j) = ygrid(ii,jj)
            zdata(i,j) = zgrid(ii,jj)
         end do
      end do

      write (luncrt, '(a)', advance='no') 'Output dataset name: '
      read  (lunkbd, *) filename_out
      open  (lun_out, file=filename_out, status='unknown')

      write (lun_out, '(2i3)') nspokes, 3
      write (lun_out, '(i3)') nidata(:)
      write (lun_out, '(3es15.7)') ((xdata(i,j), ydata(i,j), zdata(i,j), &
                                     i = 1, nidata(j)), j = 1, nspokes)
      close (lun_out)

      deallocate (xgrid, ygrid, zgrid, xdata, ydata, zdata, idata, jdata)

      ios = 1  ! Done

99    return

      end subroutine look_for_indices_option

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine symmetry_adjustments ()

!     Determine whether we appear to have half a body, the whole body, or just
!     some piece of a body, and reorganize the data if necessary so that the
!     interpolations reflect any apparent symmetry.  First, the azimuthal angles
!     may have to be made monotonically increasing or decreasing.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      real :: one_degree
      real, allocatable :: xreg_temp(:,:), yreg_temp(:,:), zreg_temp(:,:), &
                           freg_temp(:,:,:), tdat_temp(:)

!     Execution:

      one_degree = pi / 180.
      method = loose  ! Else periodic if whole body

!     We assume spoke 1 is at 12 o'clock.

      half_body  = abs (tdata(1) - pi) < one_degree .or. &
                   abs (tdata(1) + pi) < one_degree  ! ~ 12 o'clock at low end
      whole_body = half_body .and. &
                   sign (1., tdata(2)) /= sign (1., tdata(nspokes-1))
      if (whole_body) then
         half_body = .false.
         method = cyclic
      else  ! 6 o'clock <-> 0 degrees for left and right halves
         half_body = half_body .and. abs (tdata(nspokes)) < one_degree
         if (half_body) tdata(nspokes) = 0. ! exactly
      end if

      if (diagnostics) &
         write (6, '(a, 2l2)') 'half_body, whole_body:', half_body, whole_body

      if (half_body .or. whole_body) then  ! Ensure sign of tdata(1) = that of 2
         if (tdata(2) > 0.) then
             tdata(1) =  pi
         else
             tdata(1) = -pi
         end if
      end if

      if (diagnostics) write (6, '(a)') 'Adjusted azimuthal angles (radians):'

      if (half_body) then  ! Reflect the 2nd and 2nd-last spokes

         allocate (xreg_temp(ni,0:nspokes+1), yreg_temp(ni,0:nspokes+1),    &
                   zreg_temp(ni,0:nspokes+1), freg_temp(nf,ni,0:nspokes+1), &
                   tdat_temp(0:nspokes+1))
         xreg_temp(  :,0) =  xreg(  :,2)
         yreg_temp(  :,0) = -yreg(  :,2)
         zreg_temp(  :,0) =  zreg(  :,2)
         freg_temp(:,:,0) =  freg(:,:,2)
         tdat_temp(    0) =  tdata(1)*2. - tdata(2)
         xreg_temp(  :,1:nspokes) =  xreg(  :,1:nspokes)
         yreg_temp(  :,1:nspokes) =  yreg(  :,1:nspokes)
         zreg_temp(  :,1:nspokes) =  zreg(  :,1:nspokes)
         freg_temp(:,:,1:nspokes) =  freg(:,:,1:nspokes)
         tdat_temp(    1:nspokes) =  tdata(   1:nspokes)
         xreg_temp(  :,nspokes+1) =  xreg(  :,nspokes-1)
         yreg_temp(  :,nspokes+1) = -yreg(  :,nspokes-1)
         zreg_temp(  :,nspokes+1) =  zreg(  :,nspokes-1)
         freg_temp(:,:,nspokes+1) =  freg(:,:,nspokes-1)
         tdat_temp(    nspokes+1) =  tdata(nspokes)*2. - tdata(nspokes-1)

         deallocate (xreg, yreg, zreg, freg, tdata)
         allocate (xreg(ni,0:nspokes+1), yreg(ni,0:nspokes+1),    &
                   zreg(ni,0:nspokes+1), freg(nf,ni,0:nspokes+1), &
                   tdata(  0:nspokes+1))
         xreg(:,:)   = xreg_temp(:,:)
         yreg(:,:)   = yreg_temp(:,:)
         zreg(:,:)   = zreg_temp(:,:)
         freg(:,:,:) = freg_temp(:,:,:)
         tdata(:)    = tdat_temp(:)

         deallocate (xreg_temp, yreg_temp, zreg_temp, freg_temp, tdat_temp)

         if (diagnostics) &
            write (6, '(i3, es16.8)') (j, tdata(j), j = 0, nspokes+1)

      else if (whole_body) then  ! Copy spoke 1 to nspoke+1; use periodic method

         allocate (xreg_temp(ni,nspokes+1), yreg_temp(ni,nspokes+1),    &
                   zreg_temp(ni,nspokes+1), freg_temp(nf,ni,nspokes+1), &
                   tdat_temp(nspokes+1))
         xreg_temp(  :,1:nspokes) = xreg(  :,1:nspokes)
         yreg_temp(  :,1:nspokes) = yreg(  :,1:nspokes)
         zreg_temp(  :,1:nspokes) = zreg(  :,1:nspokes)
         freg_temp(:,:,1:nspokes) = freg(:,:,1:nspokes)
         tdat_temp(    1:nspokes) = tdata(   1:nspokes)
         xreg_temp(  :,nspokes+1) = xreg(  :,1)
         yreg_temp(  :,nspokes+1) = yreg(  :,1)
         zreg_temp(  :,nspokes+1) = zreg(  :,1)
         freg_temp(:,:,nspokes+1) = freg(:,:,1)
         tdat_temp(    nspokes+1) = -tdata(1)  ! +/- pi

         deallocate (xreg, yreg, zreg, freg, tdata)

         nspokes = nspokes + 1
         allocate (xreg(ni,nspokes), yreg(ni,nspokes), &
                   zreg(ni,nspokes), freg(nf,ni,nspokes), tdata(nspokes))
         xreg(:,:)   = xreg_temp(:,:)
         yreg(:,:)   = yreg_temp(:,:)
         zreg(:,:)   = zreg_temp(:,:)
         freg(:,:,:) = freg_temp(:,:,:)
         tdata(:)    = tdat_temp(:)

         deallocate (xreg_temp, yreg_temp, zreg_temp, freg_temp, tdat_temp)

         if (diagnostics) &
            write (6, '(i3, es16.8)') (j, tdata(j), j = 1, nspokes)

      end if

      end subroutine symmetry_adjustments

   end program polar_interp
