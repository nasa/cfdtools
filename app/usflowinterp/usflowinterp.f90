!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program usflowinterp

!  Purpose:
!
!     Interpolate flow field data from a US3D flow solution on to target points
!     defined in the specified form such as a Plot3D multiblock grid.  This was
!     prompted by the need to set up line-of-sight data for radiation calcula-
!     tions, analogous to flow_interp which is restricted to structured multi-
!     block flow solutions.
!
!  Method:
!
!     The flow grid may contain mixed cell types defined in the Fluent format.
!     An ADT (alternating digital tree) search tree is built from all of the
!     flow field cells.  The search for each target point always returns the
!     cell closest to the target, and the point not outside that cell that is
!     closest, along with corresponding interpolation coefficients for the cell
!     which can be used to interpolate the flow variables from the cell vertex
!     values.  This implies use of vertex-centered flow data.  If the input data
!     are cell-centered, they will be converted to the vertices internally here.
!
!     Potential Issue:  Indications elsewhere are that all zones of an un-
!                       structured dataset contain the same kind of cells.
!                       Cross that bridge if we come to it.
!
!  Control File (Standard Input):
!
!     Usflowinterp control file
!     cylinder.dat             ! Input grid + flow field data file
!     1                        ! 1|2|? = Tecplot|x|? input volume data format
!     T                        ! T|F = formatted | unformatted input volume
!     T                        ! T|F = cell-centered | vertex-centered
!     los.g                    ! Target point file name
!     1                        ! 1|2 = Plot3D/structured | Tecplot/unstructured
!     los.f                    ! Interpolated flow file name (same format)
!
!  History:
!
!     03/08/2022     D.A.Saunders  Initial design.
!     03/28-31/2022    "     "     Initial testing; seems to be working.
!     04/04/2022       "     "     Combining of zones and moving centroid data
!                                  to the vertices are now done in the
!                                  triangulation_io package.
!     04/20/2022       "     "     So is determining %nvertices (via new public
!                                  subroutine get_element_type).
!     06/22-23/2022    "     "     Handled the unstructured target case that
!                                  had been left as a stub.
!     06/24/2022       "     "     Change zone%element_type = 'FE*' to '*' for
!                                  output.  Missing DATAPACKING keyword in the
!                                  doesn't seem to bother Tecplot.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use tri_header_structure  ! Part of triangulation_io.f90, surface or volume
   use tri_zone_structure    ! Likewise
   use triangulation_io      ! Unstructured surface/volume file I/O
   use grid_block_structure  ! One Plot3D block
   use xyzq_io_module        ! Plot3D I/O utilities
   use adt_utilities         ! ADT search package (all variants)

   implicit none

!  Local constants:

   integer, parameter ::    &
      lunin   = 1,          & ! Input volume grid & flow field data
      luntarg = 2,          & ! Input target grid
      lunout  = 3,          & ! Output interpolated flow data
      lunctl  = 5,          & ! Control file (standard input)
      luncrt  = 6             ! Screen diagnostics

   logical, parameter ::    &
      false = .false.,      &
      true  = .true.

   character, parameter :: &
      format*11 = 'unformatted'

!  Local variables:

   integer :: &
      i1, ib, ier, ios, iz, kind_target, nblocks_target, ncells, nf, nfutarg, &
      ni, nj, nk, nnodes, nnodes_target, npts, ntarget, nv, nvertices, nzones

   real :: &
      dtol, time_end, time_start

   logical :: &
      cell_centered, formatted_target

   character(132) :: &
      filename_out, filename_target

!  Derived data types:

   type (tri_header_type) :: &
      header_flow, header_interp, header_target
   type (tri_type), pointer, dimension (:) :: &
      xyz_flow, xyz_target
   type (grid_type), pointer, dimension(:) :: &
      target

!  Execution:

!  Read the control file:

   call read_controls ();                           if (ios /= 0) go to 99

!  Read the flow field grid and solution:

   call read_flow_data ();                          if (ios /= 0) go to 99

!  Read the target coordinates:

   call read_target ();                             if (ios /= 0) go to 99

!  Build the search tree:

   call build_search_tree ();                       if (ier /= 0) go to 99

!  Perform the flow field interpolations, and save results by block or zone:

   call interpflow ();                              if (ier /= 0) go to 99

99 continue

!  Internal procedures for program usflowinterp:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_controls ()  ! On standard input

!     Usflowinterp control file
!     cylinder.dat             ! Input grid + flow field data file
!     1                        ! 1|2|? = Tecplot|x|? input volume data format
!     T                        ! T|F = formatted | unformatted input volume
!     T                        ! T|F = cell-centered | vertex-centered
!     los.g                    ! Target point file name
!     1                        ! 1|2 = Plot3D/structured | Tecplot/unstructured
!     los.f                    ! Interpolated flow file name (same format)

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer :: line

!     Execution:

      line = 1
      read (lunctl, *, iostat=ios)  ! Skip title
      if (ios /= 0) go to 90

      line = line + 1
      read (lunctl, *, iostat=ios) header_flow%filename  ! Input flow data file
      if (ios /= 0) go to 90

      line = line + 1
      read (lunctl, *, iostat=ios) header_flow%fileform  ! 1|2|? = Tecplot|x|?
      if (ios /= 0) go to 90                             ! input flow data form

      line = line + 1
      read (lunctl, *, iostat=ios) header_flow%formatted  ! T|F = ASCII|binary
      if (ios /= 0) go to 90

      line = line + 1
      read (lunctl, *, iostat=ios) cell_centered  ! T|F
      if (ios /= 0) go to 90

      if (header_flow%fileform == 1) then
         if (cell_centered) header_flow%fileform = 2  ! Else Tec/vertex-centered
      else
         write (luncrt, '(a, i4)') &
            'Unhandled volume data file type:', header_flow%fileform
         ios = 1
         go to 99
      end if
         
      line = line + 1
      read (lunctl, *, iostat=ios) filename_target
      if (ios /= 0) go to 90

      line = line + 1
      read (lunctl, *, iostat=ios) kind_target  ! Target file type
      if (ios /= 0) go to 90

      if (kind_target /= 1 .and. kind_target /= 2) then
         write (luncrt, '(a, i4)') 'Unhandled target file type:', kind_target
         go to 99
      end if

      line = line + 1
      read (lunctl, *, iostat=ios) filename_out  ! Output file name
      if (ios == 0) go to 99

   90 write (luncrt, '(a, i2)') &
         '*** Read_controls: Trouble reading controls on line', line

   99 continue

      end subroutine read_controls

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_flow_data ()
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

!     Execution:

!     Determine the unstructured cell type (same for all zones):

      call get_element_type (lunin, header_flow, ios)  ! Gets %nvertices
      if (ios /= 0) go to 99

!     Read the input volume dataset:

      header_flow%combine_zones         = true
      header_flow%centroids_to_vertices = true
      ios = 1                           ! Verbose option

      call vol_read (lunin, header_flow, xyz_flow, ios)
      if (ios /= 0) go to 99

      nzones = header_flow%nzones

!!!   call combine_zones (header_flow, xyz_flow)

      nvertices = header_flow%nvertices
      nnodes    = header_flow%nnodes
      ncells    = header_flow%nelements
      nf        = header_flow%numf

   99 continue

      end subroutine read_flow_data

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_target ()  ! Or just its header; set up output file too
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

!     Execution:

      select case (kind_target)

         case (1)  ! Plot3D/structured

            call determine_grid_form (filename_target, luntarg, &
                                      formatted_target, ios)
            if (ios /= 0) go to 99

            i1 = 1;  if (formatted_target) i1 = 3
            open (luntarg, file=filename_target, status='old', &
                  form=format(i1:11))

            if (ios /= 0) then
               write (luncrt, '(/, 2a)') ' Unable to open target grid: ', &
                  filename_target(1:len_trim(filename_target))
               go to 99
            end if

!           There is no need to store the entire target grid.  Read its header:

            call xyz_header_io (1, luntarg, formatted_target, nblocks_target, &
                                target, ios)
            if (ios /= 0) then
               write (luncrt, '(/, 2a)') ' Trouble allocating the target ', &
                  'grid blocks or reading the header records.'
               go to 99
            end if

!           Set up the output file of interpolated results:

            open (lunout, file=filename_out, form=format(i1:11), &
                  status='unknown', iostat=ios)
            if (ios /= 0) then
               write (luncrt, '(2a)') &
                  '*** Read_target: Unable to open output file ', &
                  len_trim (filename_out)
               go to 99
            end if

            do ib = 1, nblocks_target
               target(ib)%mi = target(ib)%ni
               target(ib)%mj = target(ib)%nj
               target(ib)%mk = target(ib)%nk
            end do

            call q_header_io (2, lunout, formatted_target, nblocks_target, &
                              nf, target, ios)
            if (ios /= 0) go to 99

         case (2)  ! Unstructured target grid

!           Formatted Tecplot unstructured surface or volume assumed.

            header_target%filename  = filename_target
            header_target%formatted = true

!           Determine the unstructured cell type (same for all zones):

            call get_element_type (luntarg, header_target, ios)  ! %nvertices
            if (ios /= 0) go to 99

            header_target%combine_zones         = false
            header_target%centroids_to_vertices = false  ! Ignore any fn. data
            ios = 1                                      ! Verbose option

!           Read the unstructured target file header:

            call vol_header_read (luntarg, header_target, xyz_target, ios)
            if (ios /= 0) go to 99

            nfutarg = header_target%numf

!           Set up for vertex-centered interpolated output, zone by zone:

            header_interp%filename    = filename_out
            header_interp%fileform    = 1    ! Vertex-centered when written
            header_interp%formatted   = true
            header_interp%nvertices   = header_target%nvertices
            header_interp%numf        = nf   ! header_flow%numf upon output
            header_interp%nzones      = header_target%nzones
            header_interp%datapacking = header_target%datapacking
            header_interp%title       = header_target%title
            header_interp%zone_type   = header_target%zone_type

            allocate (header_interp%varname(3+nf))
            header_interp%varname(:) = header_flow%varname(:)

            call tri_header_write (lunout, header_interp, xyz_target, ios)
!!!         if (ios /= 0) go to 99                      ! xyz_target isn't used
                                                        ! in this routine
      end select

   99 continue

      end subroutine read_target

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine build_search_tree ()
!
!     A volume mesh of mixed cell types is planned for, although an assumption
!     that all zones contain the same cell type appears in triangulation_io.f90
!     and hence here.
!
!     The cell type convention below appears in adt_utilities.f90.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer :: i, icell, itype, nvert_per_cell(7)
      integer, allocatable :: icontemp(:,:)

!     Storage:

       data nvert_per_cell &  ! These follow Fluent convention
         /3,               &  ! 1 = triangle
          4,               &  ! 2 = tetrahedron
          4,               &  ! 3 = quadrilateral
          8,               &  ! 4 = hexahedron
          5,               &  ! 5 = pyramid with quadrilateral base
          6,               &  ! 6 = prism with triangular cross-section
          2/                  ! 7 = line segement (not a Fluent type)

!     Execution:

      call cpu_time (time_start)

!     Adjust the connectivity data to include the cell type with each cell:

      allocate (icontemp(nvertices,ncells))

      icontemp(:,:) = header_flow%conn(:,:)
      deallocate (header_flow%conn)
!!!   allocate (header_flow%conn(0:nvertices,ncells))  ! Search_mixed_cell_adt
      allocate (header_flow%conn(0:8,ncells))          ! dimensions it 0:8

      do i = 1, 7
         if (nvert_per_cell(i) == nvertices) then
            itype = i;  exit
         end if
      end do

      do icell = 1, ncells
         header_flow%conn(0,icell) = itype
         header_flow%conn(1:nvertices,icell) = icontemp(:,icell)
      end do

      deallocate (icontemp)

!     Build the search tree from all mesh cells:

      call build_adt (nnodes, ncells, header_flow%conn, header_flow%xyz, ier)
      if (ier /= 0) go to 99

      call cpu_time (time_end)
      write (luncrt, '(/, a, 2i10, f8.2)') &
         ' Search tree # points, # cells, and tree build time, s:', &
         nnodes, ncells, time_end - time_start

   99 continue

      end subroutine build_search_tree

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine interpflow ()
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer :: i, icell, j, k, n

      real :: coefs(8), dsqmin, fv(8), xyz_interp(3), xyz_targ(3)

!     Execution:

      nv = nvertices  ! Shorter name
      call cpu_time (time_start)

      ntarget = 0

      select case (kind_target)

         case (1)  ! Plot3D target points

            do ib = 1, nblocks_target

               ni = target(ib)%ni
               nj = target(ib)%nj
               nk = target(ib)%nk

               call xyz_allocate (target(ib), ios)
               if (ios /= 0) then
                  write (luncrt, '(/, a, i4, a, i4)') &
                     ' Trouble allocating target grid.  Block #:', ib, &
                     '  I/O status:', ios
                  write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
                  go to 99
               end if

!              Read one target block at a time:

               npts = ni*nj*nk
               ntarget = ntarget + npts

               call xyz_block_io (1, luntarg, formatted_target, npts, &
                                  target(ib)%x, target(ib)%y, target(ib)%z, ios)
               if (ios /= 0) then
                  write (luncrt, '(/, a, i4, a, i4)') &
                     ' Trouble reading target grid.  Block #:', ib, &
                     '  I/O status:', ios
                  write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
                  go to 99
               end if

!              Allocate space for the corresponding interpolated flow block:

               call q_allocate (target(ib), nf, ios)
               if (ios /= 0) then
                  write (luncrt, '(/, a, i4, a, i4)') &
                     ' Trouble allocating interpolated flow space.  Block #:', &
                     ib, '  I/O status:', ios
                  write (luncrt, '(a, 4i5)') ' ni, nj, nk, nf ', ni, nj, nk, nf
                  go to 99
               end if

               do k = 1, nk
                  do j = 1, nj
                     do i = 1, ni

                        xyz_targ(1) = target(ib)%x(i,j,k)
                        xyz_targ(2) = target(ib)%y(i,j,k)
                        xyz_targ(3) = target(ib)%z(i,j,k)

                        call search_adt (xyz_targ, icell, coefs, dsqmin, &
                                         true, nnodes, ncells, &
                                         header_flow%conn, header_flow%xyz, &
                                         xyz_interp, ier)
                        if (ier /= 0) then
                           write (luncrt, '(a, 7i8)') &
                  '** Search trouble. block, i, j, k, nnodes, ncells, icell:', &
                              ib, i, j, k, nnodes, ncells, icell
                           go to 99
                        end if

                        do n = 1, nf
                           fv(1:nv) = &
                              header_flow%f(n,header_flow%conn(1:nv,icell))
                           target(ib)%q(n,i,j,k) = &
                              dot_product (coefs(1:nv), fv(1:nv))
                        end do

                     end do
                  end do
               end do

!              Save the current block of interpolations:

               call q_block_io (2, lunout, formatted_target, nf, ni, nj, nk, &
                                target(ib)%q, ios)
               if (ios /= 0) then
                  write (luncrt, '(/, a, i4, a, i4)') &
                     ' Trouble writing interpolated flow solution.  Block #:', &
                     ib, '  I/O status:', ios
                  write (luncrt, '(a, 4i5)') ' Dimensions and # functions: ', &
                     ni, nj, nk, nf
                  go to 99
               end if

               deallocate (target(ib)%x, target(ib)%y, target(ib)%z, &
                           target(ib)%q)

            end do  ! Next target block

         case (2)  ! Unstructured Tecplot surface or volume target points

            call interp_at_unstructured_targets (ntarget, ios)
            if (ios /= 0) go to 99

         case default

            write (luncrt, '(a, i3)') &
               'Unhandled target coordinate type:', kind_target
            go to 99

      end select

      close (lunout)

      call cpu_time (time_end)
      write (luncrt, '(/, a, i10, a, f8.2)') &
         ' Number of target points:', ntarget, &
         '.  Search and interpolation time, s:', time_end - time_start

   99 continue

      end subroutine interpflow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine interp_at_unstructured_targets (ntarget, ios)

!     Modularized separately because it was belatedly added as an option.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (out) :: ntarget  ! Total mpoint count, all zones
      integer, intent (out) :: ios      ! Nonzero means trouble

!     Local variables:

      integer :: i, icell, n

      real :: coefs(8), dsqmin, fv(8), xyz_interp(3), xyz_targ(3)

!     Execution:

      ntarget = 0

      do iz = 1, header_target%nzones

         nnodes_target = xyz_target(iz)%nnodes
         ntarget = ntarget + nnodes_target

         call vol_zone_allocate (header_target, xyz_target(iz), ios)
         if (ios /= 0) then
            write (luncrt, '(/, a, i4, a, i4)') &
               ' Trouble allocating target grid.  Zone #:', iz, &
               '  I/O status:', ios
            write (luncrt, '(a, i5)') ' # nodes: ', nnodes_target
            go to 99
         end if

         call vol_zone_read (luntarg, header_target, xyz_target(iz), ios)
         if (ios /= 0) then
            write (luncrt, '(/, a, i4, a, i4)') &
               ' Trouble reading target grid.  Zone #:', iz, &
               '  I/O status:', ios
            write (luncrt, '(a, i5)') ' # nodes: ', nnodes_target
            go to 99
         end if

         allocate (xyz_target(iz)%f(nf,nnodes_target))

         do i = 1, nnodes_target

            xyz_targ(1) = xyz_target(iz)%xyz(1,i)
            xyz_targ(2) = xyz_target(iz)%xyz(2,i)
            xyz_targ(3) = xyz_target(iz)%xyz(3,i)

            call search_adt (xyz_targ, icell, coefs, dsqmin, &
                             true, nnodes, ncells, &
                             header_flow%conn, header_flow%xyz, &
                             xyz_interp, ier)
            if (ier /= 0) then
               write (luncrt, '(a, 5i8)') &
                  '** Search trouble. Zone, i, nnodes, ncells, icell:', &
                  iz, i, nnodes, ncells, icell
               go to 99
            end if

            do n = 1, nf
               fv(1:nv) = &
                  header_flow%f(n,header_flow%conn(1:nv,icell))
               xyz_target(iz)%f(n,i) = dot_product (coefs(1:nv), fv(1:nv))
            end do

         end do

!        Kludge for vertex-centered output:

         if (xyz_target(iz)%element_type(1:2) == 'FE') then
             xyz_target(iz)%element_type(:) = xyz_target(iz)%element_type(3:)
         end if

!        Save the current zone of interpolations:

         call vol_zone_write (lunout, header_interp, xyz_target(iz), ios)
         if (ios /= 0) then
            write (luncrt, '(/, a, i4, a, i4)') &
               ' Trouble writing interpolated flow solution.  Zone #:', &
               iz, '  I/O status:', ios
            write (luncrt, '(a, 2i5)') ' # nodes and # functions: ', &
               nnodes_target, nf
            go to 99
         end if

         deallocate (xyz_target(iz)%xyz, xyz_target(iz)%f)

      end do  ! Next target zone

 99   return

      end subroutine interp_at_unstructured_targets

   end program usflowinterp
