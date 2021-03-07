!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program surface_diffs
!
!  Description:
!
!     Map each point of a first ("current") surface solution to a second
!     ("previous") surface (possibly the same grid, but not necessarily),
!     and output the first surface grid with the differenced solution for
!     the specified variable (one variable per run).
!
!     The surface solutions are expected to be structured Tecplot datasets.
!     Input files must be ASCII; outputs may be ASCII (.dat) or binary (.plt).
!
!     The difference distribution is written as a 2-variable dataset:
!
!        (1)  first - second                    (signed differences)
!        (2)  100 * (first - second) / |first|  (signed percentage differences)
!
!     Tabulations of minimum and maximum differences are also provided by block.
!
!  Possible uses:
!
!     >  convergence testing for solutions n iterations apart on the same grid
!     >  grid resolution comparisons
!     >  other kinds of perturbed-surface comparisons
!
!     The variable treated is likely to be surface pressure, temperature, or
!     heat flux.  See also the earlier BUMP_FACTORS program if ratios rather
!     than differences are preferred.  SURFACE_INTERP may also be relevant.
!
!  Sample control file (standard input):
!
!     SURFACE_DIFFS control file
!     qw, W/m^2   ! Function name for plottable results
!     3           ! Function in current (first) dataset to treat
!     3           ! Function in previous (second) dataset to compare with
!     0           ! 1 = flip y of current data
!     1.0000      ! Scale applied to current function
!                 ! Current file first, previous file 2nd, output 3rd
!     CFD06_10000.dat
!     CFD06_9500.dat
!     CFD06_10000-9500.dat  ! Use *.plt to get binary output
!
!  09/13/07  David Saunders  Initial adaptation of BUMP_FACTORS at Todd White's
!                            suggestion.
!  05/07/08    "      "      Todd's question about identifying the function in
!                            printable tabulation prompted showing its name.
!  08/08/13    "      "      All ADT variants have been merged into one module
!            now ERC, Inc.   with generic build_adt and search_adt interfaces.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_header_structure  ! Part of Tecplot_io.f90
   use grid_block_structure   !   "   "   "   "   "
   use tecplot_io_module      !   "   "   "   "   "
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
      one  = 1.0,        &
      zero = 0.0

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

!  Variables:

   integer :: &
      ib, iflip, ios, iq_current, iq_previous, len, len_fname, naux, &
      nb_current, nb_previous, nf, nquad, numq

   real :: &
      current_scale      ! To handle W/m^2 vs. W/cm^2, say

   integer, allocatable, dimension (:,:) :: &
      conn               ! For patch & (i,j) of surface quads. being searched

   logical :: &
      formatted

   character :: &
      fname * 64

!  Composite data types:

   type (grid_header) :: &
      current_header, diffs_header, previous_header

   type (grid_type), pointer, dimension (:) :: &
      current, diffs, previous  ! (x,y,z,f) for the three surface solutions

!  Execution:

   read (lunctl, *)
   read (lunctl, '(a)') fname;  call find_fname_length ()  ! Local procedure

   read (lunctl, *) iq_current
   read (lunctl, *) iq_previous
   read (lunctl, *) iflip
   read (lunctl, *) current_scale
   read (lunctl, *)
   read (lunctl, *) current_header%filename

   current_header%formatted = true
   current_header%ndim      = 3

   ios = 1  ! Turn on verbose mode

   call Tecplot_read (lun_current, current_header, current, ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Trouble reading current Tecplot file ', &
         trim (current_header%filename)
      go to 99
   end if

   nb_current = current_header%nblocks
   numq       = current_header%numq

   if (iq_current > numq) then
      write (luncrt, '(/, a, 2i5)') &
         ' Function # specified exceeds current dataset no. of functions:', &
         iq_current, numq
      go to 99
   end if

   if (iflip /= 0) then
      do ib = 1, nb_current
         current(ib)%y = -current(ib)%y
      end do
   end if

   if (current_scale /= one) then
      do ib = 1, nb_current
         current(ib)%q(iq_current,:,:,:) = current_scale * &
         current(ib)%q(iq_current,:,:,:)
      end do
   end if

   read (lunctl, *) previous_header%filename

   previous_header%formatted = true
   previous_header%ndim      = 3

   call Tecplot_read (lun_previous, previous_header, previous, ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Trouble reading previous Tecplot file ',    &
         trim (previous_header%filename)
      go to 99
   end if

   nb_previous = previous_header%nblocks
   numq        = previous_header%numq

   if (iq_previous > numq) then
      write (luncrt, '(/, a, 2i5)') &
         ' Function # specified exceeds previous dataset no. of functions:',  &
         iq_previous, numq
      go to 99
   end if

   read  (lunctl, *) diffs_header%filename
   close (lunctl)

!  Set up the output dataset:

   len = len_trim (diffs_header%filename)

   diffs_header%formatted = diffs_header%filename(len-2:len) /= 'plt'

   nf = 2

   call clone_header (current_header, nf, 1, diffs_header)  ! 1 = BLOCK order

   diffs_header%title = 'Differences for ' // trim (current_header%filename)  &
                        // ' - ' // trim (previous_header%filename)
   diffs_header%varname(4) = fname(1:len_fname) // ' difference'
   diffs_header%varname(5) = '% ' // fname(1:len_fname) // ' difference'

   allocate (diffs(nb_current))

   do ib = 1, nb_current

      call clone_zone (current(ib), 3, nf, diffs(ib), ios)  ! 3-D, 2 functions

      diffs(ib)%x = current(ib)%x
      diffs(ib)%y = current(ib)%y
      diffs(ib)%z = current(ib)%z

   end do

!  Construct a search tree from all previous (second) surface patches:

   nquad = 0
   do ib = 1, nb_previous
      nquad = (previous(ib)%ni - 1) * (previous(ib)%nj - 1) + nquad
   end do

   allocate (conn(3,nquad)) ! For patch # and (i,j)

   call build_adt (nb_previous, previous, nquad, conn)

!  Perform the mapping and conversion to differences and % differences:

   call map_to_previous ()

   deallocate (conn)

   call deallocate_header (current_header);  numq = current_header%numq
   call deallocate_blocks (1, nb_current, 3, numq, current, ios)

   call deallocate_header (previous_header);  numq = previous_header%numq
   call deallocate_blocks (1, nb_previous, 3, numq, previous, ios)

!  Save the difference results:

   call Tecplot_write (lun_out, diffs_header, diffs, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble writing the difference results.'
   end if

   call deallocate_header (diffs_header)

   call deallocate_blocks (1, nb_current, 3, nf, diffs, ios)

99 continue

!  Local procedures for program surface_diffs

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

      subroutine map_to_previous ()

!     Map the current grid to the previous grid and do the differencing.

!     --------------------------------------------------------------------------

!     Local constants:

      real, parameter :: big = 1.e+30, &
                         eps = 1.e-12  ! Avoid zero divides in % diffs.
!     Local variables:

      integer :: i, ib, ib1, ic, ier, iquad, j, jc, m, n, ni, ninside, nj,     &
                 noutside
      real    :: current_f, diff, diff_pc, dmax, dmean, dsqmin, dtolsq, fmin,  &
                 p, pm1, q, qm1, previous_f, interp_xyz(3), target_xyz(3)

      real, allocatable, dimension (:,:) :: extrema

!     Execution:

      allocate (extrema(4,nb_current))

!     Tolerance for search diagnostics (should be refined):

      write (luncrt, '(/, a, /)') ' Surface search statistics:'

      dtolsq = (0.001) ** 2

      do ib = 1, nb_current

         extrema(1,ib) =  big  ! For smallest absolute diff
         extrema(2,ib) = -big  !  "  biggest ...
         extrema(3,ib) =  big  !  "  smallest % diff
         extrema(4,ib) = -big  !  "  biggest ...

         ninside  = 0;  dmax  = zero
         noutside = 0;  dmean = zero

         ni = current(ib)%ni
         nj = current(ib)%nj

         do j = 1, nj

            do i = 1, ni

               target_xyz(1) = current(ib)%x(i,j,1)
               target_xyz(2) = current(ib)%y(i,j,1)
               target_xyz(3) = current(ib)%z(i,j,1)

               call search_adt (target_xyz, iquad, p, q, dsqmin, true,         &
                                nb_previous, previous, nquad, conn, interp_xyz)

               if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
                  ninside  = ninside + 1
               else
                  noutside = noutside + 1
               end if

               dmax  = max (dmax, dsqmin)
               dmean = dmean + dsqmin

!              Interpolate the surface flow at this target point:

               n   = conn(1,iquad) ! Block #
               ic  = conn(2,iquad) ! Lower left quad. indices
               jc  = conn(3,iquad)
               pm1 = one - p
               qm1 = one - q

               previous_f = &
                  qm1 * (pm1 * previous(n)%q(iq_previous,ic,  jc,  1)  + &
                           p * previous(n)%q(iq_previous,ic+1,jc,  1)) + &
                    q * (pm1 * previous(n)%q(iq_previous,ic,  jc+1,1)  + &
                           p * previous(n)%q(iq_previous,ic+1,jc+1,1))

               current_f = current(ib)%q(iq_current,i,j,1)

               diff = current_f - previous_f
               diffs(ib)%q(1,i,j,1) = diff
               extrema(1,ib) = min (diff, extrema(1,ib))
               extrema(2,ib) = max (diff, extrema(2,ib))

               diff_pc = diff * 100. / max (abs (current_f), eps)
               diffs(ib)%q(2,i,j,1) = diff_pc
               extrema(3,ib) = min (diff_pc, extrema(3,ib))
               extrema(4,ib) = max (diff_pc, extrema(4,ib))
                                      
!!!            if (ib == 1) then
!!!               if (i == 2 .and. j == 2) then
!!!                  write (6, '(a, 3f20.8)') ' target xyz: ', target_xyz
!!!                  write (6, '(a, 3f20.8)') ' interp xyz: ', interp_xyz
!!!                  write (6, '(a, 3i5)') ' b1,i2,j2: n, ic, jc: ', n, ic, jc
!!!                  write (6, '(4f20.8)') &
!!!                     previous(n)%q(iq_previous,ic:ic+1,jc:jc+1,1)
!!!                  write (6, '(a, 3f20.8)') ' p, q, finterp: ', p,q,previous_f
!!!                  write (6, '(a, f20.8)')  ' Current f:   ', &
!!!                     current(ib)%q(iq_current,i,j,1)
!!!                  write (6, '(a, 2f20.8)') ' Differences: ', &
!!!                     diffs(ib)%q(:,i,j,1)
!!!               end if
!!!            end if

            end do ! Next i for this target block

         end do ! Next j for this target block

         write (luncrt, '(a, i3, a, 2i7, a, 1p, 2e12.5)')                      &
            '   current patch', ib,                                            &
            ':  # points in/outside tol.:', ninside, noutside,                 &
            ';  max/mean devn.:', sqrt (dmax), sqrt (dmean) / real (ni * nj)

      end do ! Next target block

      write (luncrt, '(/, 3a, //, a, /)') &
         ' Summary of ', fname(1:len_fname),                                   &
         ' differences (first - second)[ * 100 / |first|]:',                   &
         ' Block        min abs        max abs          min %          max %'

      write (luncrt, '(i6, 1p, 4e15.6)') (ib, extrema(:,ib), ib = 1, nb_current)

      deallocate (extrema)

      end subroutine map_to_previous

   end program surface_diffs
