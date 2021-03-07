!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program ensemble
!
!  Description:
!
!     This specialized application interpolates CFD surface solutions in and
!  around idealized cavities (representing thermal protection system pieces
!  missing from a space vehicle) on to standardized surface patches with point
!  distributions derived from the cavity dimensions.  This allows ensemble
!  averaging of solutions for multiple cavities with different dimensions.
!
!     In retrospect, treating the general case from the start may have led to
!  a cleaner program able to handle all cases.  But nobody mentioned off-center
!  cavities till 04/13/05.
!
!  Original Assumptions (Topologies 1 - 4):
!
!     >  The cavities are rectangular and aligned with OX-OY-OZ.
!        Cavities with rounded corners are not handled initially.
!     >  X is downstream; Y is spanwise; Z is up.
!     >  Cavities at the center of a flat plate are treated as half cavities;
!        the half with Y < 0 is assumed.
!     >  For a given grid topology, the surface patch numbers need not be
!        consistent; input lists can point to patches in a consistent order;
!        the pointers group patches as follows:
!
!           cavity floor, downstream wall, upstream wall, side wall, smooth OML
!
!  More Generally (Topology 5):
!
!     >  No symmetry means treating the full cavity, not half.
!     >  Discrepancies resulting from side walls that aren't quite uniform in
!        height should be negligible.
!     >  As input, we can get by with the patch numbers of one zone per face
!        of the cavity, as long as each patch is completely on a single face.
!        (Some LaRC grids have patches that are on more than one side wall.)
!
!  Topologies Handled:
!
!     >  Topology 1 is the H-H topology from Bill Wood.
!        Topology 2 is the Ames O-H topology.
!        Topology 3 is the Ames Shuttle centerline cavity topology: reflect Y, Z
!        Topology 4 is the LARC Shuttle centerline morphed cavity:  reflect Y, Z
!        Topology 5 is the general case (arbitrary orientation, LAURA XYZ)
!        Topology 6 is the general case (arbitrary orientation, DPLR XYZ)
!
!  Control File ('ensemble.inp'):
!
!     Ensemble control file
!     ensemble.dat         ! Output file name (Tecplot format)
!     T                    ! T = formatted; F = unformatted
!     1+ line pairs:  topology, name, in/out BF var #s, ref. val; 5 patch lists
!     1  r150BW-HH-GEP.dat     3  7  1.
!     14  9    12 15     7 10    13  8 16 11    6 5 4 3 1 2
!     1  r155BW-HH-GEP.dat     3  7  1.
!     14  9    12 15     7 10    13  8 16 11    6 5 4 3 1 2
!      :  :     :  :     :  :     :  :  :  :    : : : : : :
!      :  :     :  :     :  :     :  :  :  :    : : : : : :
!     2  r150JJR-OH-GEP.dat
!     20 15    18 21    13 16    19 14 22 17     2  4  3  1
!     2  r155CT-OH-GEP.dat     3  7  1.
!     22 19    20 15    17 13    21 18 16 14     9 10 12 11
!      :  :     :  :     :  :     :  :  :  :     :  :  :  :
!
!     The above integers apply to topologies 1-4.
!
!     For topology 5, enter single patch numbers fully on each cavity face or
!     the OML, in the following order:
!
!        floor,  back wall,  front wall,  outboard wall,  inboard wall,  OML
!
!  Original Input Surface Solution Format (Tecplot ASCII, structured):
!
!        TITLE     = ""
!        VARIABLES = "x"
!        "y"
!        "z"
!        "T"
!        "p"
!        "qw"
!        "h"
!        "?"
!        "?"
!        "?"
!        ZONE T="G1"
!         I=66, J=66, K=1, ZONETYPE=Ordered
!         DATAPACKING=BLOCK
!         DT=(DOUBLE DOUBLE ... (variable # of them) ... DOUBLE DOUBLE )
!         6.647333145E+00 6.575638294E+00 6.489704609E+00 6.390774727E+00 6....
!          :               :               :               :               :
!
!  Alternative Input Solution Format:
!
!        As above except there may be different flow quantities (possibly
!        including or only a heating "bump" factor).
!
!        POINT data packing is also permitted.
!
!  Output Results (Tecplot binary or ASCII file with NORMALIZED target grid):
!
!     This is similar to the input format but with averaged flow quantities and
!     a second set of variables representing standard deviations.  If only one
!     file is found, the ensemble averaging results are suppressed.
!
!  History:
!
!     04/09/05  D. Saunders  Initial implementation (wind tunnel cavities).
!     04/11/05       "       Added Shuttle flight Ames cavity topology 3.
!     04/12/05       "       Added Shuttle flight morphed cavity topology 4.
!     04/13/05       "       Added the general Shuttle case (LAURA xyz).
!     04/15/05       "       LAURA wind tunnel cases have just 2 functions.
!     04/16/05       "       Full path names failed.  i_topology > 0 suppresses
!                            most of the diagnostic output; < 0 activates it;
!                            allowed scaling of a specified function as a way
!                            of converting it to "bump factor."
!     04/18/05       "       Topology 6 is as for toplogy 5, but it suppresses
!                            the initial conversion from LAURA to DPLR xyz.
!     05/13/05       "       Tecplot_read's datapacking argument is now an
!                            output, not an input.  Either POINT or BLOCK can
!                            be handled without knowing which here. Tecplot_read
!                            now returns variable names as it should have.
!
!  Author:  David Saunders, ELORET/NASA Ames, Moffett Field, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use adt_utilities
   use grid_block_structure    ! Must be the Tecplot_io version
   use grid_header_structure
   use surface_patch_utilities
   use tecplot_io_module
   use trigd

   implicit none

!  Constants:

   integer, parameter ::        &
      lunctl     = 1,           &  ! Control file
      lun_in     = 2,           &  ! Input CFD solution (Tecplot ASCII)
      lun_out    = 3,           &  ! Output normalized grid & soln. (Tecplot)
      lunsrf     = 4,           &  ! Actual interpolations for real space target
      luncrt     = 6,           &
      max_patch_per_region = 6, &  ! Most patches for any region in any topology
      name_limit = 32,          &  ! Must match the value in Tecplot_io_module
      nb_half    = 7,           &  ! # blocks in standardized half-cavity grid
      n_topology = 6               ! # topologies handled so far

   real, parameter :: &
      one  = 1.0,     &
      zero = 0.0

   character, parameter :: &
      blank  * 1  = ' ',   &
      format * 11 = 'unformatted'

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

!  Variables:

   integer :: &
      i, i_bump_var_in, i_bump_var_out, i_topology, ib, ic,       &
      ifile, ios, iquad, j, jc, n, nb_out, nblocks, nfiles, ni, ninside, nj,   &
      noutside, nquad, numf_in, numf_out

   integer :: &
      kf, kd, ku, ko, ki, kt  ! Patch numbers for the general case

   integer :: &
      first, last, mark       ! Parsing variables

   integer, dimension (5, n_topology) :: &
      n_patch_per_region                 ! See data statement below

   integer, dimension (max_patch_per_region) :: &
      list_floor, list_back, list_front, list_side, list_oml

   integer, allocatable, dimension (:,:) :: &
      conn               ! For patch & (i,j) of surface quads. being searched

   real :: &
      D, L, W, XB, ZT, dmax, dmean, dsqmin, dtolsq, p, pm1, q, qm1, term,      &
      bump_factor_reference_value, bump_scale

   real, dimension (3) :: &
      interp_xyz, target_xyz

   logical :: &
      convert_to_BF, diagnostics, formatted, whole

   character :: &
      buffer * 132, &
      filename_in * 80, filename_out * 80, title_out * 80

   character, pointer, dimension (:) :: &
      names_in  * (name_limit),         &
      names_out * (name_limit)

!  Composite data types:

   type (grid_type), pointer, dimension (:) :: &
      xyzq_in,     &  ! one CFD surface solution
      xyzq_out,    &  ! means and standard deviations, normalized standard grid
      xyzq_target, &  ! standardized patches on current CFD surface grid
      interp_surf     ! actual interpolated surface points found by ADT searches
   type (grid_header) :: header_in, header_out

!  Data:

   data n_patch_per_region &
     /2, 2, 2, 4, 6,       &  ! Topology 1 (H-H, WT, Bill Wood)
      2, 2, 2, 4, 4,       &  !     "    2 (O-H, WT, Ames)
      2, 1, 1, 2, 4,       &  !     "    3 (O-H, Flight, Ames)
      6, 1, 1, 2, 4,       &  !     "    4 (O-H, Flight, LaRC)
      1, 1, 1, 1, 2,       &  !     "    5 (General, Flight, LaRC xyz)
      1, 1, 1, 1, 2/          !     "    6 (General, Flight,  ARC xyz)

!  Execution:
!  !!!!!!!!!!

   open (lunctl, file='ensemble.inp', status='old', iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') &
         ' Unable to open ensemble.inp control file.'
      go to 99
   end if

!  Count the number of input data files, so we can suppress rms data if there
!  is only one of them.

   do i = 1, 999
      read (lunctl, '(a)', iostat=ios) buffer
      if (ios < 0) exit

      j = len_trim (buffer) ! Catch possible trailing blank lines
      if (j > 0) then
         write (luncrt, '(a)') buffer(1:j)
         nfiles = i
      end if
   end do

   nfiles = (nfiles - 4) / 2
   write (luncrt, '(/, a, i4)') ' Number of files to process:', nfiles

   rewind (lunctl)

   read (lunctl, *)                  ! Case description
   read (lunctl, *) filename_out     ! Output Tecplot file name
   read (lunctl, *) formatted
   read (lunctl, *)                  ! Description

!  Process one CFD solution at a time.  Finish initializing after reading # 1.
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do ifile = 1, 999 ! Until end of control file

!     Belated control inputs after the file name force parsing the input line:

      read (lunctl, '(a)', iostat=ios) buffer
      if (ios < 0) exit  ! EOF

      first = 1;  last = len_trim (buffer)
      call scan2 (buffer, blank, first, last, mark)
      if (buffer(first:first) == '-') then
         read (buffer(first:mark), '(i2)') i_topology
      else
         read (buffer(first:mark), '(i1)') i_topology
      end if

      first = mark + 2
      call scan2 (buffer, blank, first, last, mark)

      filename_in = buffer(first:mark)
      j = mark - first + 1

      first = mark + 2
      read (buffer(first:last), *, iostat=ios) &
         i_bump_var_in, i_bump_var_out, bump_factor_reference_value

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') &
            ' Trouble reading bump factor controls from this string:',         &
            buffer(first:last)
         go to 99
      end if

      diagnostics = i_topology < 0
      i_topology  = abs (i_topology)
      bump_scale  = one / bump_factor_reference_value
      convert_to_BF   =   bump_factor_reference_value /= one

      write (luncrt, '(/, 2a, /)') ' Reading flow solution ', filename_in(1:j)

      ! Use new API
      !call Tecplot_read (lun_in, filename_in, true, datapacking, title_in,     &
      !                   nblocks, numf_in, names_in, xyzq_in, ios)
      call Tecplot_read (lun_in, header_in, xyzq_in, ios)
      nblocks  = header_in%nblocks
      numf_in  = header_in%numq
      names_in = header_in%varname

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') ' Trouble reading Tecplot file ', &
            trim (header_in%filename)
         go to 99
      end if

      if (convert_to_BF) then
         if (i_bump_var_in < 0 .or. i_bump_var_in > numf_in) then
            write (luncrt, '(/, a, i3)') ' Bad i_bump_var_in:', i_bump_var_in
            go to 99
         end if
         if (i_bump_var_out < 0 .or. i_bump_var_out > numf_in) then
            write (luncrt, '(/, a, i3)') ' Bad i_bump_var_out:', i_bump_var_out
            go to 99
         end if
         do ib = 1, nblocks
            xyzq_in(ib)%q(i_bump_var_out,:,:,1) = bump_scale * &
            xyzq_in(ib)%q(i_bump_var_in, :,:,1)
         end do
      end if

!     The following patch numbers help determine actual cavity location & size:

      read (lunctl, *) &
         list_floor(1:n_patch_per_region(1,i_topology)), &
         list_back (1:n_patch_per_region(2,i_topology)), &
         list_front(1:n_patch_per_region(3,i_topology)), &
         list_side (1:n_patch_per_region(4,i_topology)), &
         list_oml  (1:n_patch_per_region(5,i_topology))  ! Inboard side and OML
                                                         ! for the general case
      if (i_topology <= 4) then ! Symmetric data
         nb_out = nb_half
         whole  = false
      else ! Whole cavity
         nb_out = nb_half + nb_half
         whole  = true
      end if

      if (ifile == 1) then

         if (nfiles == 1) then
            numf_out = numf_in
         else
            numf_out = numf_in * 2 ! Add rms values
         end if

!        Generate a "unit" standardized [half-]cavity grid:

!        The following is an internal procedure below.  It allocates too.

         call standard_cavity_grid (xyzq_out, numf_out, true)

         if (nfiles > 1) then

!           Set up for accumulating solution means & standard deviations:

            do ib = 1, nb_out
               xyzq_out(ib)%q(:,:,:,1) = zero
            end do

         end if

!        Save [the last] interpolated real-space grid for checking.

         open (lunsrf, file='interp_surf.p3da', status='unknown')

         allocate (interp_surf(nb_out))

         do i = 1, nb_out
            interp_surf(i)%ni = xyzq_out(i)%ni
            interp_surf(i)%nj = xyzq_out(i)%nj
            interp_surf(i)%nk = 1
         end do

         write (lunsrf, '(i2)') nb_out
         write (lunsrf, '(3i4)') &
            (interp_surf(i)%ni, interp_surf(i)%nj, 1, i = 1, nb_out)

      end if

!     Generate a standard 7-patch half-cavity surface grid in real space,
!     after removing any rotations off-center shifts.
!     For the general case, it is a 14-patch full cavity, including surrounds.
!     Adjust the xyz convention, rotate, etc., if necessary.

      call find_cavity_dimensions (nblocks, xyzq_in) ! Sets L, W, D, XB, ZT

      call standard_cavity_grid (xyzq_target, numf_in, false)

      dtolsq = (L * 0.001) ** 2  ! Tolerance for ADT search diagnostics

!     Interpolate the current grid and function values onto this target grid.
!     First, build a search tree for the current surface solution.
!     It doesn't matter much if there are some redundant patches in the soln.
!     The ADT technique is barely affected.

      nquad = 0
      do ib = 1, nblocks
         nquad = (xyzq_in(ib)%ni - 1) * (xyzq_in(ib)%nj - 1) + nquad
      end do

      allocate (conn(3,nquad)) ! For patch # and (i,j)

      call build_adt (nblocks, xyzq_in, nquad, conn)

      do ib = 1, nb_out

         ni = xyzq_target(ib)%ni
         nj = xyzq_target(ib)%nj

         if (ifile == 1) then
            allocate (interp_surf(ib)%x(ni,nj,1), interp_surf(ib)%y(ni,nj,1),  &
                      interp_surf(ib)%z(ni,nj,1))
         end if

         ninside  = 0;  dmax  = zero
         noutside = 0;  dmean = zero

         do j = 1, nj

            do i = 1, ni

               target_xyz(1) = xyzq_target(ib)%x(i,j,1)
               target_xyz(2) = xyzq_target(ib)%y(i,j,1)
               target_xyz(3) = xyzq_target(ib)%z(i,j,1)

               call search_adt (target_xyz, iquad, p, q, dsqmin, true,         &
                                nblocks, xyzq_in, nquad, conn, interp_xyz)

               if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
                  ninside  = ninside + 1
               else
                  noutside = noutside + 1
               end if

               dmax  = max (dmax, dsqmin)
               dmean = dmean + dsqmin

               if (ifile == 1) then
                  interp_surf(ib)%x(i,j,1) = interp_xyz(1)
                  interp_surf(ib)%y(i,j,1) = interp_xyz(2)
                  interp_surf(ib)%z(i,j,1) = interp_xyz(3)
               end if

!              Interpolate the surface flow at this target point:

               n   = conn(1,iquad) ! Block #
               ic  = conn(2,iquad) ! Lower left quad. indices
               jc  = conn(3,iquad)
               pm1 = one - p
               qm1 = one - q

               if (nfiles == 1) then

                  xyzq_out(ib)%q(:,i,j,1) = &
                     qm1 * (pm1 * xyzq_in(n)%q(:,ic,  jc,  1)  + &
                              p * xyzq_in(n)%q(:,ic+1,jc,  1)) + &
                       q * (pm1 * xyzq_in(n)%q(:,ic,  jc+1,1)  + &
                              p * xyzq_in(n)%q(:,ic+1,jc+1,1))
               else

                  xyzq_target(ib)%q(:,i,j,1) = &
                     qm1 * (pm1 * xyzq_in(n)%q(:,ic,  jc,  1)  + &
                              p * xyzq_in(n)%q(:,ic+1,jc,  1)) + &
                       q * (pm1 * xyzq_in(n)%q(:,ic,  jc+1,1)  + &
                              p * xyzq_in(n)%q(:,ic+1,jc+1,1))

!                 Update the mean flow values and their mean squared deviations:

                  do n = 1, numf_in
                     term = (xyzq_target(ib)%q(n,i,j,1) - &
                             xyzq_out   (ib)%q(n,i,j,1)) / real (ifile)
                     xyzq_out(ib)%q(n,i,j,1) = xyzq_out(ib)%q(n,i,j,1) + term
                     xyzq_out(ib)%q(n+numf_in,i,j,1) =  &
                     xyzq_out(ib)%q(n+numf_in,i,j,1) +  &
                        real (ifile * (ifile - 1)) * term ** 2
                  end do

               end if

            end do ! Next i for this target block

         end do ! Next j for this target block

         write (luncrt, '(a, i3, a, 2i7, a, 1p, 2e12.5)')                      &
            '   Standardized patch', ib,                                       &
            ':  # pts. in/outside tolerance:', ninside, noutside,              &
            ';  max/mean devn.:', sqrt (dmax), sqrt (dmean) / real (ni * nj)

         if (ifile == 1) then
            write (lunsrf, '(1p, 6e19.11)') &
               interp_surf(ib)%x, interp_surf(ib)%y, interp_surf(ib)%z

               deallocate (interp_surf(ib)%x, interp_surf(ib)%y,               &
                           interp_surf(ib)%z)
         end if

      end do ! Next target block

      deallocate (conn)

      do ib = 1, nblocks
         deallocate (xyzq_in(ib)%x, xyzq_in(ib)%y, xyzq_in(ib)%z,              &
                     xyzq_in(ib)%q)
      end do

      do ib = 1, nb_out
         deallocate (xyzq_target(ib)%x, xyzq_target(ib)%y, xyzq_target(ib)%z,  &
                     xyzq_target(ib)%q)
      end do

      if (ifile == 1) close (lunsrf)

   end do ! Next file to process


   if (nfiles > 1) then

      term = one / real (nfiles)

!     Turn the mean squared deviations into rms values:

      do ib = 1, nb_out
         do j = 1, xyzq_out(ib)%nj
            do i = 1, xyzq_out(ib)%ni
               do n = numf_in + 1, numf_out
                 xyzq_out(ib)%q(n,i,j,1) = sqrt (xyzq_out(ib)%q(n,i,j,1)) * term
               end do
            end do
         end do
      end do

   end if

   write (luncrt, '(/, a, i3)') ' # files processed:', nfiles

   allocate (names_out(3+numf_out))

!! Kludge until Tecplot_read handles variable names properly:

!! allocate (names_in(3+numf_in))

!! names_in(1) = 'x'
!! names_in(2) = 'y'
!! names_in(3) = 'z'

   names_out(1:3) = names_in(1:3)

!! if (numf_in == 7) then ! DPLR cases
!!    names_in(4) = 'T'
!!    names_in(5) = 'p'
!!    names_in(6) = 'qw'
!!    names_in(7) = 'h'
!!    names_in(8) = 'C_p'
!!    names_in(9) = 'Cf'
!!    names_in(10)= 'Ct'
!! else if (numf_in == 2) then ! LAURA WT / topology 1 after halving
!!    names_in(4) = 'qw'
!!    names_in(5) = 'BF'
!! else if (numf_in == 1) then
!!    names_in(4) = 'BF'
!! else
!!    write (luncrt, '(/, a, i4)') ' Unexpected # functions:', numf_in
!!    names_in(4:numf_in)(1:1) = 'f'
!!    do i = 4, min (numf_in, 9)
!!       if (i < 13) write (names_in(i)(2:2), '(i1)') i - 3
!!       if (i > 12) write (names_in(i)(2:3), '(i2)') i - 3
!!    end do
!! end if

   if (convert_to_BF) names_in(3 + i_bump_var_out) = 'BF'

   if (nfiles == 1) then ! Suppress 'mean' in the name, and the zero rms values

      do n = 1, numf_in  ! numf doesn't count x, y, z
         names_out(n+3) = trim (names_in(n+3))
      end do

   else ! We're ensemble averaging

      do n = 1, numf_in  ! numf doesn't count x, y, z
         names_out(n+3)         = trim (names_in(n+3)) // ' mean'
         names_out(n+3+numf_in) = trim (names_in(n+3)) // ' rms'
      end do

   end if

   title_out = trim (filename_out)
   write (luncrt, '(a, a)') ' Writing Tecplot file ', title_out

   ! Use new API
   !call Tecplot_write (lun_out, filename_out, formatted, 2, title_out,         &
   !                    nb_out, numf_out, names_out, xyzq_out, ios)
   header_out%filename  = filename_out
   header_out%formatted = formatted
   header_out%ndim      = 2
   header_out%nblocks   = nb_out
   header_out%numq      = numf_out
   header_out%title     = title_out
   header_out%varname   = names_out
   call Tecplot_write(lun_out, header_out, xyzq_out, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble writing the output Tecplot file.'
      go to 99
   end if

99 continue

!  Internal procedures for program ensemble:

   contains

!     -----------------------------------------
      subroutine find_cavity_dimensions (nb, g)
!     -----------------------------------------

!     Specialized routine for determining the size of a rectangular cavity
!     at the symmetry plane (or not):
!
!        L = length   W = half-width   D = depth   XB = back wall x   ZT = top z

!     Arguments:

      integer, intent (in) :: nb                            ! # blocks in g
      type (grid_type), dimension (nb), intent (inout) :: g ! Saves typing
                                                            ! xyzq_in; xmin out
!     Local constants:

      real, parameter :: big = 1.e+32, half = 0.5

!     Local variables:

      integer :: i, ib, j, ll, ni, nj, npts
      real    :: ci, si, theta
      real    :: dx, dy, dz, xd, px, py, pz, qx, qy, qz
      real    :: xf, xmax, xmin, ymax, ymin, zd, zf, zmax, zmin
      logical :: symmetric

!     Execution:

!     The wind tunnel cases originally treated have the Y < 0 half of the
!     cavity and Z pointing up from the floor.  It is too confusing to change
!     this, so convert other topologies to the initial convention:

      if (i_topology == 3 .or. i_topology == 4) then

         do ib = 1, nb
            g(ib)%y = -g(ib)%y  ! Switch to WT convention
            g(ib)%z = -g(ib)%z
         end do

!        Arrange to rotate the grid about the center of the cavity on the OML:

         ib = list_back(1)
         ni = g(ib)%ni;  nj = g(ib)%nj
         zd = -big

         do j = 1, nj, nj - 1
            do i = 1, ni, ni - 1
               if (zd < g(ib)%z(i,j,1)) then
                   zd = g(ib)%z(i,j,1)
                   xd = g(ib)%x(i,j,1)
               end if
            end do
         end do

         ib = list_front(1)
         ni = g(ib)%ni;  nj = g(ib)%nj
         zf = -big

         do j = 1, nj, nj - 1
            do i = 1, ni, ni - 1
               if (zf < g(ib)%z(i,j,1)) then
                   zf = g(ib)%z(i,j,1)
                   xf = g(ib)%x(i,j,1)
               end if
            end do
         end do

         write (luncrt, '(a, 1p, 4e15.7)') &
            ' Back & front (xd,zd) and (xf,zf): ', xd, zd, xf, zf

         theta = -atand ((zd - zf) / (xd - xf)) ! Anticlockwise in (x,z) space
         xd = half * (xd + xf)                  ! Center of rotation
         zd = half * (zd + zf)
         ci = cosd (theta)
         si = sind (theta)

         write (luncrt, '(a, 1p, e19.11, a, 1p, 2e19.11)') &
            ' Rotating CFD grid', theta, ' deg. about (xd,zd) =', xd, zd

         do ib = 1, nb
            npts = g(ib)%ni * g(ib)%nj

            call rotate2d (npts, g(ib)%x, g(ib)%z, theta, xd, zd)

         end do

      else if (i_topology == 5 .or. & ! General cavity in LAURA coordinates
               i_topology == 6) then  !    "      "     "  DPLR     "

         if (i_topology == 5) then    ! Flip LAURA x <- -z, y <- -y, & z <- -x.

            do ib = 1, nb
               call swap_coordinates (g(ib), 1, 3)
               do i = 1, 3
                  call reflect_x_y_or_z (g(ib), i)
               end do
            end do

            if (diagnostics) then
               open  (50, file='initial_transform.p3da', status='unknown')
               write (50, '(i3)') nb
               write (50, '(3i4)') (g(i)%ni, g(i)%nj, 1, i = 1, nb)
               do i = 1, nb
                  write (50, '(1p, 6e19.11)') g(i)%x, g(i)%y, g(i)%z
               end do
               close (50)
            end if

         end if

!        Take out any rotations in the arbitrary orientation.
!        We assume the given cavity wall patches are ESSENTIALLY rectangular.

         kf = list_floor(1)
         kd = list_back(1)
         ku = list_front(1)
         ko = list_side(1)
         ki = list_OML(1)
         kt = list_OML(2)

!        Find the angle in the back wall to rotate about a side edge.

         ni = g(kd)%ni
         dy = g(kd)%y(ni,1,1) - g(kd)%y(1,1,1)
         dz = g(kd)%z(ni,1,1) - g(kd)%z(1,1,1)
         theta = atand (dz / dy)

         px = g(ko)%x(1,1,1)
         py = g(ko)%y(1,1,1)
         pz = g(ko)%z(1,1,1)

         ni = g(ko)%ni
         qx = g(ko)%x(ni,1,1)
         qy = g(ko)%y(ni,1,1)
         qz = g(ko)%z(ni,1,1)

         if (px < qx) theta = -theta ! RH rule looking upstream from p to q
         write (6, '(/, a, 1p, 3e19.11)') 'dy,dz,theta1: ', dy,dz,theta
         write (6, '(a, 1p, 3e19.11)') 'px,py,pz: ', px,py,pz
         write (6, '(a, 1p, 3e19.11)') 'qx,qy,qz: ', qx,qy,qz

         symmetric = abs (g(ko)%y(1,1,1) + g(ki)%y(1,1,1)) * half  <           &
                     abs (g(ko)%y(1,1,1) - g(ki)%y(1,1,1)) * 0.01

         if (symmetric) then
            write (luncrt, '(a)') ' Treating cavity as being on the centerline.'
         else
            do ib = 1, nb
               npts = g(ib)%ni * g(ib)%nj

               call rotate3d (npts, g(ib)%x, g(ib)%y, g(ib)%z, theta, &
                              px, py, pz, qx, qy, qz)
            end do

            if (diagnostics) then
               open  (50, file='rotation1.p3da', status='unknown')
               write (50, '(i3)') nb
               write (50, '(3i4)') (g(i)%ni, g(i)%nj, 1, i = 1, nb)
               do i = 1, nb
                  write (50, '(1p, 6e19.11)') g(i)%x, g(i)%y, g(i)%z
               end do
               close (50)
            end if

         end if

!        Find angle in the new outer side wall to rotate about new back edge:

         dx = g(ko)%x(ni,1,1) - g(ko)%x(1,1,1)
         dz = g(ko)%z(ni,1,1) - g(ko)%z(1,1,1)
         theta = atand (dz / dx)

         px = g(kd)%x(1,1,1)
         py = g(kd)%y(1,1,1)
         pz = g(kd)%z(1,1,1)

         ni = g(kd)%ni
         qx = g(kd)%x(ni,1,1)
         qy = g(kd)%y(ni,1,1)
         qz = g(kd)%z(ni,1,1)

         if (py > qy) theta = -theta
         write (6, '(/, a, 1p, 3e19.11)') 'dx,dz,theta2: ', dx,dz,theta
         write (6, '(a, 1p, 3e19.11)') 'px,py,pz: ', px,py,pz
         write (6, '(a, 1p, 3e19.11)') 'qx,qy,qz: ', qx,qy,qz

         do ib = 1, nb
            npts = g(ib)%ni * g(ib)%nj

            call rotate3d (npts, g(ib)%x, g(ib)%y, g(ib)%z, theta, &
                           px, py, pz, qx, qy, qz)
         end do

         if (diagnostics) then
            open  (50, file='rotation2.p3da', status='unknown')
            write (50, '(i3)') nb
            write (50, '(3i4)') (g(i)%ni, g(i)%nj, 1, i = 1, nb)
            do i = 1, nb
               write (50, '(1p, 6e19.11)') g(i)%x, g(i)%y, g(i)%z
            end do
            close (50)
         end if

!        Finally, rotate in the new (x,y) plane:

         ni = g(ko)%ni
         dx = g(ko)%x(ni,1,1) - g(ko)%x(1,1,1)
         dy = g(ko)%y(ni,1,1) - g(ko)%y(1,1,1)
         theta = atand (dy / dx)

         px = g(kd)%x(1,1,1)
         py = g(kd)%y(1,1,1)

         if (dy < zero) theta = -theta
         write (6, '(/, a, 1p, 3e19.11)') 'dx,dy,theta3: ', dx,dy,theta
         write (6, '(a, 1p, 2e19.11, /)') 'px,py: ', px,py

         do ib = 1, nb
            npts = g(ib)%ni * g(ib)%nj

            call rotate2d (npts, g(ib)%x, g(ib)%y, theta, px, py)

         end do

!        Shift the grid spanwise to the center line:

         dy = half * (g(ko)%y(1,1,1) + g(ki)%y(1,1,1))
         do ib = 1, nb
            do j = 1, g(ib)%nj
               do i = 1, g(ib)%ni
                  g(ib)%y(i,j,1) = g(ib)%y(i,j,1) - dy
               end do
            end do
         end do

      end if ! Reflections and rotations to real space wind tunnel-like axes

      if (i_topology > 2) then
         open  (50, file='reflected_rotated.p3da', status='unknown')
         write (50, '(i3)') nb
         write (50, '(3i4)') (g(i)%ni, g(i)%nj, 1, i = 1, nb)
         do i = 1, nb
            write (50, '(1p, 6e19.11)') g(i)%x, g(i)%y, g(i)%z
         end do
         close (50)
      end if

!     Establish the (possibly adjusted) input grid data ranges:

      do ib = 1, nb
         call patch_data_range (g(ib))
         if (diagnostics) then
            write (6, '(a, i4, 6f15.10)') '  Range of [adjusted?] patch', ib,  &
               g(ib)%xmin, g(ib)%xmax, g(ib)%ymin, g(ib)%ymax, g(ib)%zmin,     &
               g(ib)%zmax
         end if
      end do

!     Establish cavity dimensions in (possibly adjusted) real space:

      if (.not. whole) then ! Original method

         xmin = big;  xmax = -big
         ymin = big;  ymax = -big
         zmin = big;  zmax = -big

!        Cavity floor:

         do ll = 1, n_patch_per_region(1,i_topology)
            ib = list_floor(ll)
            if (diagnostics) write (6, '(a, 2i4)') '   floor ll, ib: ', ll, ib
            xmax = max (xmax, g(ib)%xmax)
            xmin = min (xmin, g(ib)%xmin)
            ymax = max (ymax, g(ib)%ymax)
            ymin = min (ymin, g(ib)%ymin)
         end do

         L = xmax - xmin;  XB = xmax
         W = ymax - ymin

!        Cavity depth:

         do ll = 1, n_patch_per_region(2,i_topology)
            ib = list_back(ll)
            if (diagnostics) write (6, '(a, 2i4)') '   back wall ll, ib: ',ll,ib
            zmax = max (zmax, g(ib)%zmax)
            zmin = min (zmin, g(ib)%zmin)
         end do

         D = zmax - zmin;  ZT = zmax

      else ! Whole cavity case (later strategy)

         ZT = half * (g(kt)%zmax + g(kt)%zmin)      ! Top
         D  = ZT - half * (g(kf)%zmax + g(kf)%zmin) ! Top - floor
         XB = half * (g(kd)%xmax + g(kd)%xmin)      ! Back/downstream wall
         L  = XB - half * (g(ku)%xmax + g(ku)%xmin) ! Downstream - upstream
         W  = half * ((g(ki)%ymax + g(ki)%ymin) - & ! Full width
                      (g(ko)%ymax + g(ko)%ymin))
         W  = half * W                              ! Half width

      end if

      end subroutine find_cavity_dimensions

!     --------------------------------------------------
      subroutine standard_cavity_grid (g, nf, normalize)
!     --------------------------------------------------

!     Construct a standard surface grid for the current rectangular cavity.
!     For the whole cavity case, generate half a grid then reflect it.
!     The half grid is in the y <= 0 domain for historical reasons to do
!     with the CFD results for wind tunnel cavity cases.
!     The blocks and associated q arrays are necessarily allocated here.
!
!     Normalized:  L = 1       W = 0.5      D = 1      XB = 1       ZT = 1
!     Real space:  L = length  W = width/2  D = depth  XB = back x  ZT = top z
!
!     Patch    Location                  Size      Dimensions
!
!       1      Floor                     L x W     m x n  where n ~ 65
!       2      Back wall  (downstream)   W x D     n x m  and   m = 2n - 1
!       3      Front wall (upstream)     W x D     n x m
!       4      Side wall  (Y < 0)        L x D     m x m
!       5      OML, adjacent to cavity   L x 2W    m x n
!       6      Upstream of cavity        L x 3W    n x m
!       7      Downstream of cavity     2L x 3W    m x m
!     [ 8      Floor reflected ...
!       9      ...
!       :
!      14      Downstream patch reflected ]
!
!     Strategy:  Construct the edges with Vinokur distributions.
!                Filling the interiors doesn't really need TFIQ3D.
!
!     Arguments:

      type (grid_type), pointer :: g(:)   ! xyzq_target or xyzq_out
      integer, intent (in) :: nf          ! numf_in     or numf_out
      logical, intent (in) :: normalize   ! False       or true
                                          ! Tecplot_write requires pointer

!     Local constants:

      integer, parameter :: n = 65, m = 2*n - 1

!     Local variables:

      integer :: i, ier, j, j1, j2, lunerr, ni, nj
      integer, save :: lunout ! For saving standardized grids
      real    :: d1, d2, t(m)
      real, dimension (nb_half) :: d1i, d2i, d1j, d2j

!     Data:

      data d1i /0.001, 0.002, 0.002, 0.001, 0.001, 0.050, 0.001/
      data d2i /0.001, 0.005, 0.005, 0.001, 0.001, 0.001, 0.100/
      data d1j /0.002, 0.001, 0.001, 0.001, 0.050, 999.0, 999.0/
      data d2j /0.005, 0.001, 0.001, 0.001, 0.002, 999.0, 999.0/

!     The above have to be applied to L, W, or D appropriately.

!     Execution:

      lunerr = luncrt
      if (diagnostics) lunerr = -luncrt

!     We construct the y <= 0 half grid first, then reflect if necessary.

      if (normalize) then ! Assume this is true on the first call
         lunout = 10
         L = one;  W = 0.5;  D = one;  XB = one;  ZT = one
      else
         lunout = lunout + 1
      end if

      allocate (g(nb_out))

      g(1)%ni = m;  g(1)%nj = n
      g(2)%ni = n;  g(2)%nj = m
      g(3)%ni = n;  g(3)%nj = m
      g(4)%ni = m;  g(4)%nj = m
      g(5)%ni = m;  g(5)%nj = n
      g(6)%ni = n;  g(6)%nj = m
      g(7)%ni = m;  g(7)%nj = m

      if (.not. whole) then ! Half cavity

         g(1)%zone_title = 'Floor'
         g(2)%zone_title = 'Back wall'
         g(3)%zone_title = 'Front wall'
         g(4)%zone_title = 'Side wall'
         g(5)%zone_title = 'OML alongside'
         g(6)%zone_title = 'OML upstream'
         g(7)%zone_title = 'OML downstream'

      else ! Whole cavity

         do i = 1, nb_half
            g(i+nb_half)%ni = g(i)%ni
            g(i+nb_half)%nj = g(i)%nj
         end do

         g(1)%zone_title = 'L_Floor'
         g(2)%zone_title = 'L_Back wall'
         g(3)%zone_title = 'L_Front wall'
         g(4)%zone_title = 'L_Side wall'
         g(5)%zone_title = 'L_OML alongside'
         g(6)%zone_title = 'L_OML upstream'
         g(7)%zone_title = 'L_OML downstream'

         g( 8)%zone_title = 'R_Floor'
         g( 9)%zone_title = 'R_Back wall'
         g(10)%zone_title = 'R_Front wall'
         g(11)%zone_title = 'R_Side wall'
         g(12)%zone_title = 'R_OML alongside'
         g(13)%zone_title = 'R_OML upstream'
         g(14)%zone_title = 'R_OML downstream'

      end if

      do i = 1, nb_out
         g(i)%nk = 1;  call Tec_block_allocate (g(i), 3, nf, ier)
      end do

!     Half-cavity grid point distributions:

!     SUBROUTINE VINOKUR (I1, I2, D1, D2, X, LUNOUT, IER)
!     Vinokur expects the first and last X to be in place.

!     Floor Xs & Zs:

      ni = g(1)%ni;  t(1) = XB - L;  t(ni) = XB

      if (diagnostics) then
         write (6, '(a)') ' call Vinokur (1)'
         write (6, '(a, 5f13.8)') 'L, W, D, XB, ZT: ', L, W, D, XB, ZT
         write (6, '(a, 2f13.8)') 't1 & t(ni): ', t(1), t(ni)
      end if

      call vinokur (1, ni, d1i(1)*L, d2i(1)*L, t, lunerr, ier)

      do j = 1, g(1)%nj
         g(1)%x(1:ni,j,1) = t(1:ni)
         g(1)%z(1:ni,j,1) = (ZT - D)
      end do

!     Side wall Xs & Ys:

      do j = 1, g(4)%nj
         g(4)%x(1:ni,j,1) = t(1:ni)
         g(4)%y(1:ni,j,1) = -W
      end do

!     Adjacent top Xs & Zs:

      do j = 1, g(5)%nj
         g(5)%x(1:ni,j,1) = t(1:ni)
         g(5)%z(1:ni,j,1) = ZT
      end do

!     Floor Ys:

      nj = g(1)%nj;  t(1) = -W;  t(nj) = zero

      if (diagnostics) write (6, '(a)') ' call Vinokur (2)'
      call vinokur (1, nj, d1j(1)*W, d2j(1)*W, t, lunerr, ier)

      do i = 1, ni
         g(1)%y(i,1:nj,1) = t(1:nj)
      end do

!     Front and back wall Xs & Ys:

      do j = 1, g(2)%nj
         g(2)%x(1:nj,j,1) = XB
         g(2)%y(1:nj,j,1) = t(1:nj)
         g(3)%x(1:nj,j,1) = (XB - L)
         g(3)%y(1:nj,j,1) = t(1:nj)
      end do

!     OML Ys & Xs for subpatches fore and aft of cavity:

      j2 = g(6)%nj;  j1 = j2 - nj + 1

      do i = 1, g(6)%ni
         g(6)%y(i,j1:j2,1) = t(1:nj)
         g(6)%z(i,j1:j2,1) = ZT
      end do

      do i = 1, g(7)%ni
         g(7)%y(i,j1:j2,1) = t(1:nj)
         g(7)%z(i,j1:j2,1) = ZT
      end do

!     Side wall Zs:

      nj = g(4)%nj;  t(1) = ZT - D;  t(nj) = ZT

      if (diagnostics) write (6, '(a)') ' call Vinokur (3)'
      call vinokur (1, nj, d1j(4)*D, d2j(4)*D, t, lunerr, ier)

      do i = 1, g(4)%ni
         g(4)%z(i,1:nj,1) = t(1:nj)
      end do

!     Front and back wall Zs:

      do i = 1, g(2)%ni
         g(2)%z(i,1:nj,1) = t(1:nj)
         g(3)%z(i,1:nj,1) = t(1:nj)
      end do

!     Adjacent patch Ys:

      nj = g(5)%nj;  t(1) = -3.*W;  t(nj) = -W

      if (diagnostics) write (6, '(a)') ' call Vinokur (4)'
      call vinokur (1, nj, d1j(5)*W, d2j(5)*W, t, lunerr, ier)

      do i = 1, g(5)%ni
         g(5)%y(i,1:nj,1) = t(1:nj)
      end do

!     Ys for subpatches fore and aft of patch adjacent to cavity:

      do i = 1, g(6)%ni
         g(6)%y(i,1:nj,1) = t(1:nj)
      end do

      do i = 1, g(7)%ni
         g(7)%y(i,1:nj,1) = t(1:nj)
      end do

!     Upstream Xs and Zs:

      ni = g(6)%ni;  t(1) = XB - 2.*L;  t(ni) = XB - L

      if (diagnostics) write (6, '(a)') ' call Vinokur (5)'
      call vinokur (1, ni, d1i(6)*L, d2i(6)*L, t, lunerr, ier)

      do j = 1, g(6)%nj
         g(6)%x(1:ni,j,1) = t(1:ni)
         g(6)%z(1:ni,j,1) = ZT
      end do

!     Downstream Xs and Zs:

      ni = g(7)%ni;  t(1) = XB;  t(ni) = XB + 2.*L

      if (diagnostics) write (6, '(a)') ' call Vinokur (6)'
      call vinokur (1, ni, d1i(7)*L, d2i(7)*L, t, lunerr, ier)

      do j = 1, g(7)%nj
         g(7)%x(1:ni,j,1) = t(1:ni)
         g(7)%z(1:ni,j,1) = ZT
      end do

!     Reflect this half to make a whole?

      if (whole) then
         do i = 1, nb_half
            j = i + nb_half
            g(j)%x =  g(i)%x
            g(j)%y = -g(i)%y
            g(j)%z =  g(i)%z
         end do
      end if

      if (diagnostics) then

         if (normalize) then
            open (lunout, file='normalized_cavity.p3da', status='unknown')
         else
            j = index (trim (filename_in), '.') - 1
            open (lunout, file=filename_in(1:j) // '_standard.p3da', &
                  status='unknown')
         end if

         write (lunout, '(i1)') nb_out
         write (lunout, '(3i4)') (g(i)%ni, g(i)%nj, 1, i = 1, nb_out)
         do i = 1, nb_out
            write (lunout, '(1p, 6e19.11)') g(i)%x, g(i)%y, g(i)%z
         end do
         close (lunout)

      end if

      end subroutine standard_cavity_grid

   end program ensemble
