!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program tri_diffs
!
!  Description:
!
!     This a surface triangulation version of the earlier SURFACE_DIFFS.
!
!     Map each point of a first ("current") surface solution to a second
!     ("previous") surface (possibly the same grid, but not necessarily),
!     and output the first surface grid with the differenced solution for
!     the specified variable (one variable per run).
!
!     The surface solutions are expected to be vertex-centered triangulated
!     Tecplot datasets.  See related utility TRI_TO_TRI if handling of
!     cell-centered triangulations needs to be added here.
!
!     Input files must be ASCII; outputs are ASCII (.dat).
!
!     The difference distribution is written as a 2-variable dataset:
!
!        (1)  first - second                    (signed differences)
!        (2)  100 * (first - second) / |first|  (signed percentage differences)
!
!     Tabulations of minimum and maximum differences are also provided by zone.
!
!  Possible uses:
!
!     >  convergence testing for solutions n iterations apart on the same grid
!     >  grid resolution comparisons
!     >  other kinds of perturbed-surface comparisons
!
!     The variable treated is likely to be surface pressure, temperature, or
!     heat flux, or a radio signal attenuation measure, which is what prompted
!     this variation of SURFACE_DIFFS.
!
!  Input Tecplot data format (vertex-centered):
!
!     VARIABLES = "X", "Y", "Z", "TEMP", "PRESS", "CP", "MACHE", "ASOUNDE", ...
!     ZONE N=96000, E=32000, F=FEPOINT, ET=TRIANGLE
!     0.000000  0.000000 0.000000 2838.386719 51330.925781 1.552663 0.609412 ...
!     0.000883 -0.007150 0.113643 2838.386719 51330.925781 1.552663 0.609412 ...
!     0.000883  0.000000 0.113868 2838.386719 51330.925781 1.552663 0.609412 ...
!     ::::::::::::::::
!     ::::::::::::::::
!     4.882953  0.000000 0.011285 950.867676 16.506409 -0.001166 5.062649 ...
!     1 2 3
!     4 5 6
!     7 8 9
!     10 11 12
!     ::::::::
!     ::::::::
!     95992 95993 95994
!     95995 95996 95997
!     95998 95999 96000
!
!  Sample control file (standard input):
!
!     TRI_DIFFS control file
!     f1 attenuation, dB   ! Function name for plottable results
!     1           ! Function in current (first) dataset to treat
!     1           ! Function in previous (second) dataset to compare with
!     0           ! 1 = flip y of current data
!     1.0000      ! Scale applied to current function
!                 ! Current file first, previous file 2nd, output 3rd
!     ant1_13.attenuation_v6.20_h54.00.dat
!     ant1-regular.dat
!     ant1-freq1-diff.dat  ! Use *.plt to get binary output
!
!  08/08/13  D. A. Saunders  SURFACE_DIFFS version adapted here as TRI_DIFFS.
!  04/21/17-   "   "   "     Initial TRI_DIFFS adaptation, to enable comparison
!  04/24/17                  of radio attenuation datasets from a laminar flow
!                            solution with similar (perhaps turbulent) post-
!                            processing results from the scripts employed for
!                            quality-checking Orion/MPCV aerothermal database
!                            calculations.  Multiple zones work only if the
!                            two meshes match exactly.  Hooks for multizone
!                            form of ADT build and search when available.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure   ! Needed by adt_utilities
   use tri_header_structure   ! Part of triangulation_io.f90
   use tri_zone_structure     ! Likewise
   use triangulation_io       ! Unstructured data file I/O
   use adt_utilities          ! All ADT variants

   implicit none

!  Constants:

   integer, parameter :: &
      lun_current  =  1, &    ! First/current surface CFD soln. (Tecplot ASCII)
      lun_previous =  2, &    ! Second/previous CFD solution (Tecplot ASCII)
      lun_out      =  3, &    ! Output solution of differences
      lunctl       =  5, &    ! Control file (standard input)
      luncrt       =  6       ! For diagnostics

   real, parameter ::    &
      eps_sq = 1.e-14,   &    ! Grids match if |dxyz| < 1.e-7
      one  = 1.0,        &
      zero = 0.0

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

!  Variables:

   integer :: &
      if_current, if_previous, iflip, ios, itri, iv, iz, len_fname, nf, nnode, &
      ntri, numf, nz_current, nz_previous

   integer, allocatable, dimension (:,:) :: &
      conn               ! For 3 vertices and zone # of elements being searched

   real :: &
      current_scale, &   ! To handle W/m^2 vs. W/cm^2, say
      dxyz(3)

   real, allocatable, dimension (:,:) :: extrema

   logical :: &
      formatted, same_grids

   character (80) :: &
      fname

!  Derived data types:

   type (tri_header_type) :: &
      current_header, diffs_header, previous_header

   type (tri_type), pointer, dimension (:) :: &
      current, diffs, previous  ! (x,y,z,f) for the three surface solutions

!  Execution:

   read (lunctl, *)
   read (lunctl, '(a)') fname;  call find_fname_length ()  ! Local procedure

   read (lunctl, *) if_current       ! Function numbers
   read (lunctl, *) if_previous
   read (lunctl, *) iflip
   read (lunctl, *) current_scale
   read (lunctl, *)
   read (lunctl, *) current_header%filename

   current_header%fileform   = 1     ! Vertex-centered
   current_header%formatted  = true  ! Can't read binaries
   current_header%nvertices  = 3     ! Triangles, not tetrahedra

   ios = 1  ! Turn on verbose mode

   call tri_read (lun_current, current_header, current, ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Trouble reading current Tecplot file ', &
         trim (current_header%filename)
      go to 99
   end if

   nz_current = current_header%nzones
   numf       = current_header%numf

   if (if_current > numf) then
      write (luncrt, '(/, a, 2i5)') &
         ' Function # specified exceeds current dataset no. of functions:', &
         if_current, numf
      go to 99
   end if

   if (iflip /= 0) then
      do iz = 1, nz_current
         current(iz)%xyz(2,:) = -current(iz)%xyz(2,:)
      end do
   end if

   if (current_scale /= one) then
      do iz = 1, nz_current
         current(iz)%f(if_current,:) = current_scale*current(iz)%f(if_current,:)
      end do
   end if

   read (lunctl, *) previous_header%filename

   previous_header%formatted = true
   previous_header%nvertices = 3

   call tri_read (lun_previous, previous_header, previous, ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Trouble reading previous Tecplot file ',    &
         trim (previous_header%filename)
      go to 99
   end if

   nz_previous = previous_header%nzones
   numf        = previous_header%numf

   if (if_previous > numf) then
      write (luncrt, '(/, a, 2i5)') &
         ' Function # specified exceeds previous dataset no. of functions:',  &
         if_previous, numf
      go to 99
   end if

   read  (lunctl, *) diffs_header%filename
   close (lunctl)

!  Set up the output dataset:

   nf = 2                         ! Absolute and % diffs for one function
   diffs_header%fileform    = 1   ! Vertex-centered output
   diffs_header%formatted   = true
   diffs_header%nvertices   = 3
   diffs_header%numf        = nf
   diffs_header%nzones      = nz_current
   diffs_header%datapacking = 1  ! BLOCK order
   diffs_header%title = 'Differences for ' // trim (current_header%filename)  &
                        // ' - ' // trim (previous_header%filename)
   allocate (diffs_header%varname(3+2))
   diffs_header%varname(1) = 'x, m'
   diffs_header%varname(2) = 'y, m'
   diffs_header%varname(3) = 'z, m'
   diffs_header%varname(4) = fname(1:len_fname) // ' difference'
   diffs_header%varname(5) = '% ' // fname(1:len_fname) // ' difference'

   allocate (diffs(nz_current))

   do iz = 1, nz_current

      diffs(iz)%zone_title   = current(iz)%zone_title
      diffs(iz)%nnodes       = current(iz)%nnodes
      diffs(iz)%nelements    = current(iz)%nelements
      diffs(iz)%element_type = current(iz)%element_type

      allocate (diffs(iz)%conn(3,diffs(iz)%nelements), &
                diffs(iz)%xyz (3,diffs(iz)%nnodes),    &
                diffs(iz)%f   (2,diffs(iz)%nnodes))

      diffs(iz)%xyz  = current(iz)%xyz
      diffs(iz)%conn = current(iz)%conn

   end do

!  Check for matching surface meshes:

   same_grids = true
   if (nz_current == nz_previous) then
      do iz = 1, nz_current
         if (current(iz)%nnodes == previous(iz)%nnodes) then
            if (current(iz)%nelements == previous(iz)%nelements) then
               do iv = 1, current(iz)%nnodes
                  dxyz(:) = current(iz)%xyz(:,iv) - previous(iz)%xyz(:,iv)
                  if (dot_product (dxyz, dxyz) > eps_sq) same_grids = false
               end do
            else
               same_grids = false
            end if
         else
            same_grids = false
         end if
         if (.not. same_grids) exit
      end do
   else
      same_grids = false
   end if

!  For each point of each current zone, compute differences and % differences:

   write (luncrt, '(/, a, /)') ' Surface search statistics:'

   allocate (extrema(4,nz_current))

   if (same_grids) then  ! More likely than not

      call simple_case ()

   else  !  Handle the case of different surface triangulations

      call map_to_previous ()

   end if

   write (luncrt, '(/, 3a, //, a, /)') &
      ' Summary of ', fname(1:len_fname),                                   &
      ' differences (first - second)[ * 100 / |first|]:',                   &
      ' Zone         min abs        max abs          min %          max %'

   write (luncrt, '(i6, 4es15.6)') (iz, extrema(:,iz), iz = 1, nz_current)

   deallocate (extrema)

!  Save the difference results:

   call tri_write (lun_out, diffs_header, diffs, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble writing the difference results.'
   end if

99 continue

!  Local procedures for program tri_diffs

   contains

!     --------------------------------------------------------------------------

      subroutine find_fname_length ()

!     ! Isolate the function name from a possible trailing comment.
!     ! Allow embedded blanks in the name by searching backwards.

!     --------------------------------------------------------------------------

      integer :: i, l

      len_fname = len_trim (fname)
      l = index (fname(1:len_fname), '!')

      if (l /= 0) then
         do i = l - 1, 1, -1
            if (fname(i:i) /= ' ') then
               len_fname = i
               exit
            end if
         end do
      end if

      end subroutine find_fname_length

!     --------------------------------------------------------------------------

      subroutine simple_case ()  ! Both files have the same surface mesh

!     --------------------------------------------------------------------------

!     Local constants:

      real, parameter :: big = 1.e+30, &
                         eps = 1.e-12  ! Avoid zero divides in % diffs.
!     Local variables:

      integer :: i, iz
      real    :: current_f, diff, diff_pc, previous_f

!     Execution:

      do iz = 1, current_header%nzones
         extrema(1,iz) =  big  ! For smallest absolute diff
         extrema(2,iz) = -big  !  "  biggest ...
         extrema(3,iz) =  big  !  "  smallest % diff
         extrema(4,iz) = -big  !  "  biggest ...
         do i = 1, current(iz)%nnodes
            diffs(iz)%xyz(:,i) = current(iz)%xyz(:,i)
            diff = current(iz)%f(if_current,i) - previous(iz)%f(if_previous,i)
            diffs(iz)%f(1,i) = diff
            extrema(1,iz) = min (diff, extrema(1,iz))
            extrema(2,iz) = max (diff, extrema(2,iz))

            diff_pc = diff * 100. / max (abs (current(iz)%f(if_current,i)), eps)
            diffs(iz)%f(2,i) = diff_pc
            extrema(3,iz) = min (diff_pc, extrema(3,iz))
            extrema(4,iz) = max (diff_pc, extrema(4,iz))
         end do
      end do

      end subroutine simple_case

!     --------------------------------------------------------------------------

      subroutine map_to_previous ()

!     Map the current grid to the previous grid and do the differencing.

!     --------------------------------------------------------------------------

!     Local constants:

      real, parameter :: big = 1.e+30, &
                         eps = 1.e-12  ! Avoid zero divides in % diffs.
!     Local variables:

      integer :: i, i1, i2, i3, ib, itri, iz, ninside, noutside, ntri
      real    :: current_f, diff, diff_pc, dmax, dmean, dsqmin, dtolsq, &
                 p, q, r, previous_f, interp_xyz(3), target_xyz(3)

!     Execution:

      nnode = 0
      ntri  = 0
      do iz = 1, nz_previous  ! Needs to be 1 zone for now
         nnode = previous(iz)%nnodes    + nnode
         ntri  = previous(iz)%nelements + ntri
      end do

      if (nz_previous /= 1) then
         write (luncrt, '(/, a)') ' Second file must have just 1 zone for now.'
         stop
      end if

!!    allocate (conn(4,ntri))  ! 3 pointers & zone # for all elements, all zones

!     TO BE COMPLETED:  We need multizone versions of BUILD and SEARCH  :(

      call build_adt (nz_previous, ntri, previous(1)%conn, previous(1)%xyz)

!     Tolerance for search diagnostics (should be refined):

      dtolsq = (0.001) ** 2

      do iz = 1, nz_current

         extrema(1,iz) =  big  ! For smallest absolute diff
         extrema(2,iz) = -big  !  "  biggest ...
         extrema(3,iz) =  big  !  "  smallest % diff
         extrema(4,iz) = -big  !  "  biggest ...

         ninside  = 0;  dmax  = zero
         noutside = 0;  dmean = zero

         do i = 1, current(iz)%nnodes

            target_xyz(:) = current(iz)%xyz(:,i)

            call search_adt (target_xyz, itri, p, q, r, dsqmin, true, &
                             nnode, ntri, &
                             previous(1)%conn, previous(1)%xyz, interp_xyz)

            if (dsqmin < dtolsq) then ! The nearest cell was within tolerance
               ninside  = ninside + 1
            else
               noutside = noutside + 1
            end if

            dmax  = max (dmax, dsqmin)
            dmean = dmean + dsqmin

!           Interpolate the surface flow at this target point:

            i1  = previous(1)%conn(1,itri)
            i2  = previous(1)%conn(2,itri)
            i3  = previous(1)%conn(3,itri)
!!          ib  = conn(4,itri) ! 
               
            previous_f = p * previous(iz)%f(if_previous,i1) + &
                         q * previous(iz)%f(if_previous,i2) + &
                         r * previous(iz)%f(if_previous,i3)

            current_f = current(iz)%f(if_current,i)

            diff = current_f - previous_f
            diffs(iz)%f(1,i) = diff
            extrema(1,iz) = min (diff, extrema(1,iz))
            extrema(2,iz) = max (diff, extrema(2,iz))

            diff_pc = diff * 100. / max (abs (current_f), eps)
            diffs(iz)%f(2,i) = diff_pc
            extrema(3,iz) = min (diff_pc, extrema(3,iz))
            extrema(4,iz) = max (diff_pc, extrema(4,iz))
                                      
!!!         if (iz == 1) then
!!!            if (i == 1) then
!!!               write (6, '(a, 3f20.8)') ' target xyz: ', target_xyz
!!!               write (6, '(a, 3f20.8)') ' interp xyz: ', interp_xyz
!!!               write (6, '(a, 3i5)') ' zone,tri: ib, itri: ', ib, itri
!!!               write (6, '(3f20.8)') previous(ib)%(:,itri)
!!!               write (6, '(a, 4f20.8)') ' p, q, rm finterp: ', p,q,r, &
!!!                  previous_f
!!!               write (6, '(a, f20.8)')  ' Current f:   ', &
!!!                  current(iz)%f(if_current,i)
!!!               write (6, '(a, 2f20.8)') ' Differences: ', &
!!!                  diffs(iz)%f(:,i)
!!!            end if
!!!         end if

         end do ! Next vertex for this target zone

         write (luncrt, '(a, i3, a, 2i7, a, 2es12.5)')         &
            '   current zone', iz,                             &
            ':  # points in/outside tol.:', ninside, noutside, &
            ';  max/mean devn.:', sqrt (dmax), sqrt (dmean) /  &
            current(iz)%nnodes

      end do ! Next target zone

!!!   deallocate (conn)

      end subroutine map_to_previous

   end program tri_diffs
