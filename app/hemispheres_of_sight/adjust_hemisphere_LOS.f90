!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program adjust_hemisphere_LOS

!  Description:
!
!     This is a specialized utility prompted by the need to deal with possible
!  concavities in the aft bodies of atmospheric entry vehicles such as the
!  Mars Science Laboratory aeroshell, which has nonmonotonic aft cone angles,
!  when full angular integration of radiative heat flux is being performed.
!  It reads the results of interpolating flow field data onto the initial
!  lines of sight produced by HEMISPHERE_LINES_OF_SIGHT, tests whether any of
!  those lines appear to encounter the body surface before their intersection
!  with the outer grid boundary, and suppresses the flow data at and outboard
!  of any such encounter.
!
!     The strategy is to march along each line of sight from the originating
!  body point end, testing whether each discretized line point appears to be
!  inside the body or not.  If such an inside point is found, it and all grid
!  points on that line further outboard have their flow data diminished as
!  opposed to truncating the line of sight, which would affect the angular
!  integration step by distorting even further the already distorted hemisphere
!  surface formed by HEMISPHERE_LINES_OF_SIGHT's intersection of all the lines
!  with the outer grid boundary.
!
!     A refinement would be to redistribute the unblocked portion of the line
!  in a way that captures the boundary layer at the point of blockage.  But
!  initially, we work with just the interpolated flow from the FLOW_INTERP step,
!  which has lost some of that resolution.  If the body is essentially convex,
!  only the lines at shallow angles are prone to blockage, and their radiation
!  contributions normal to the body at the body point are small because of the
!  sine term applied.
!
!     If the aft body is truly concave, as seen in some sample return capsules,
!  it remains to be seen if the above strategy suffices and if the suggested
!  refinement is essential.
!
!     In anticipation of likely complications, the underlying body surface is
!  obtained from the full volume data that may well be needed.  As with many
!  related utilities, this is taken to be the k = 1 surface of all blocks.
!
!     We assume that the HEMISPHERE_LINES_OF_SIGHT step successfully intersected
!  the outer boundary (k = kmax) for all lines, producing some lines that pass
!  through the body in places.
!
!     Since extreme cases of lines 90 degrees from the body normal (i.e.,
!  tangential to the body) may be completely clipped, NEQAIR will have nothing
!  to work with if we just zero out the flow data.  The angular integration
!  step needs a result at every hemisphere line too, even if it is only zero.
!  Therefore, for clipped portions of lines, we set the four temperature
!  values above the clipping limit (500 K) used in the NEQAIR_DATA step that
!  converts PLOT3D format to LOS.dat form, but set all the number densities
!  extremely low so that the NEQAIR result is essentially zero.
!
!  History:
!
!     02/18/15  D.A.Saunders  Initial implementation.
!     04/20/15    "     "     Record the line numbers that are affected.
!                             Aaron Brandis wants to include NEQAIR's black
!                             body option on such lines.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, California.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! For one structured grid block
   use adt_utilities         ! ADT search package (all variants)
   use xyzq_io_module        ! For PLOT3D files

   implicit none

!  Constants:

   integer, parameter :: lunvolg  = 1  ! Input flow volume grid and flow data
   integer, parameter :: lunvolf  = 2  ! used by HEMISPHERES... and FLOW_INTERP
   integer, parameter :: lunhemig = 3  ! Intersected line-of-sight data grid
   integer, parameter :: lunhemif = 4  ! Corresponding flow data from FLOW_INT.
   integer, parameter :: lunoutf  = 7  ! Adjusted line of sight function data
   integer, parameter :: lunadj   = 8  ! List of line numbers affected
   integer, parameter :: mxadjust = 20 ! Limits output diagnostics
   real,    parameter :: Tlow  = 600.  ! Must be above NEQAIR_DATA's low limit
   real,    parameter :: Dlow  = 1000. ! Low number density
   logical, parameter :: false = .false.
   logical, parameter :: true  = .true.

!  Variables:

   integer :: ib, ios, k, last, l, nadjust, nblocks_hemi, nblocks_vol, nf, nquad
   real    :: dsq
   logical :: cell_centered, formatted_hemi, formatted_vol, inside

   integer,          allocatable :: conn(:,:)
   type (grid_type), allocatable :: bodysurface(:)
   type (grid_type), pointer     :: flowvoldata(:), hemilosdata(:)

!  Execution:

!  Read the volume data:

   call file_prompt (lunvolg, -1, 'HEMISPHERES_OF_SIGHT volume grid', &
                     'old', true, formatted_vol, ios)
   if (ios /= 0) go to 99

   call file_prompt (lunvolf, -1, 'associated FLOW_INTERP volume data', &
                     'old', false, formatted_vol, ios)
   if (ios /= 0) go to 99

   call xyzq_read (lunvolg, lunvolf, formatted_vol, nblocks_vol, nf, &
                   cell_centered, flowvoldata, ios)
   if (ios /= 0) go to 99

!  Read the hemisphere LOS dataset that needs to be adjusted:

   call file_prompt (lunhemig, -1, 'hemisphere LOS grid', &
                     'old', true, formatted_hemi, ios)
   if (ios /= 0) go to 99

   call file_prompt (lunhemif, -1, 'associated LOS flow data', &
                     'old', false, formatted_hemi, ios)
   if (ios /= 0) go to 99

   call xyzq_read (lunhemig, lunhemif, formatted_hemi, nblocks_hemi, nf, &
                   cell_centered, hemilosdata, ios)
   if (ios /= 0) go to 99

!  Set up the body surface grid for inside/outside calculations:

   allocate (bodysurface(nblocks_vol))

   nquad = 0
   do ib = 1, nblocks_vol
      bodysurface(ib)%ni = flowvoldata(ib)%ni
      bodysurface(ib)%nj = flowvoldata(ib)%nj
      bodysurface(ib)%nk = 1
      nquad = (bodysurface(ib)%ni - 1)*(bodysurface(ib)%nj - 1) + nquad

      call xyz_allocate (bodysurface(ib), ios)

      bodysurface(ib)%x(:,:,1) = flowvoldata(ib)%x(:,:,1)
      bodysurface(ib)%y(:,:,1) = flowvoldata(ib)%y(:,:,1)
      bodysurface(ib)%z(:,:,1) = flowvoldata(ib)%z(:,:,1)
   end do

!  Construct a search tree for rapid determination of the nearest surface
!  point to each LOS point:

   allocate (conn(3,nquad))  ! For surface patch # and cell i, j of each quad

   call build_adt (nblocks_vol, bodysurface, nquad, conn)

!  Check each hemisphere line of sight for apparent intersection with the body:

   open (lunadj, file='adjusted_lines.txt', status='unknown')
   nadjust = 0

   do l = 1, nblocks_hemi  ! One LOS per block
      do k = 2, hemilosdata(l)%nk
         call inside_outside (nblocks_vol, bodysurface, nquad, conn, &
                              hemilosdata(l)%x(1,1,k), &
                              hemilosdata(l)%y(1,1,k), &
                              hemilosdata(l)%z(1,1,k), dsq, inside)
         if (inside) then
            hemilosdata(l)%q(1:4,1,1,k:) = Tlow
            hemilosdata(l)%q(5: ,1,1,k:) = Dlow
            nadjust = nadjust + 1
            if (nadjust <= mxadjust) then
               write (*, '(a, i6, a, i4)') '   Clipping LOS', l, ' at k =', k
            else if (nadjust == mxadjust + 1) then
               write (*, '(a)') '   Suppressing further commentary.'
            end if
            write (lunadj, '(2i5)') l, k
            exit  ! Go to the next LOS
         end if
      end do
   end do

   write (*, '(a, i6)') '   Total number of lines adjusted:', nadjust

   if (nadjust > 0) then  ! Save the adjusted function data
      open (lunoutf, file='adjusted_LOS.f', status='unknown')
      call q_write (lunoutf, true, nblocks_hemi, nf, hemilosdata, ios)
   end if

99 continue

   end program adjust_hemisphere_LOS
