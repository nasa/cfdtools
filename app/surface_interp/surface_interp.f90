!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program surface_interp
!
!  Description:
!
!     Original:
!        Read a structured surface dataset (Tecplot ASCII) and a second target
!        dataset of one or more zones (possibly unstructured points in a simpler
!        format).  Interpolate the first dataset at each target point and save
!        results as a Tecplot file and (for the simpler target format) as a
!        tabulation for possible pasting into a spreadsheet.  The output file(s)
!        may contain just a subset of the input functions.
!
!     Current:
!        Interpolate within a surface dataset at the points of a second dataset.
!        Both datasets may be one of these choices:
!           Tecplot/structured, Tecplot unstructured (tri or quad cells), Plot3D
!        The target points may also be a simple list of (x,y,z) coordinates.
!        The interpolated results match the format of the target dataset, with
!        any subset of the input functions being output.
!
!  Strategy:
!
!     The ADT (Alternating Digital Tree) search technique employed always finds
!     the point not outside the surface grid being searched that is nearest to
!     the target point, along with corresponding interpolation coefficients for
!     use within the relevant cell.  Thus, extrapolation is not an option, so
!     the target points are assumed to be on or very close to the surface being
!     searched.
!
!  Control file format (standard input):
!
!     SURFACE_INTERP control file
!     SRL-surf.gu      ! File to interp. in (*.g[u] if Plot3D; *.f[u] assumed)
!     1                ! 1|2|3 = Plot3D|Tecplot structured|Tecplot unstructured
!     target-surf.dat  ! Target file to interpolate to
!     2                ! 1|2|3|4 = Plot3D|Tec struct.|Tec unstruct.|(x,y,z) list
!     Subset of functions in output file(s):
!     1 2 3
!     SRL.interp.dat   ! Plottable results; *.f[u] if Plot3D target grid
!     SRL.interp.txt   ! Tabulated results (target mode 4 only)
!     F                ! Interpolation details at initial target point?
!     2                ! 1|2 = BLOCK|POINT order output
!     0.779347  0.     ! X & F to use for x > X
!
!  Random target points format (mode 4):
!
!     x1   y1   z1   Body point 1       ! Any string may follow the coordinates
!     x2   y2   z2   Body point 2       ! else x is inserted as a tag
!     x3   y3   z3   Wing point 1
!     :    :    :    :
!
!  Tabulated output format (mode 1):
!
!     x1   y1   z1   f11  f41  f21  Body point 1 (or x1)  Min. distance
!     x2   y2   z2   f12  f42  f22  Body point 2 (or x2)    "      "
!     x3   y3   z3   f13  f43  f23  Wing point 1 (or x3)    "      "
!     :    :    :    :
!
!  History:
!
!     07/31/05  D.A.Saunders  Initial implementation, adapted from BUMP_FACTORS
!                             in a hurry during Discovery's return to flight.
!     08/12/05    "    "      Random target points may now be read from a simple
!                             list rather than requiring Tecplot format.  The
!                             smooth/damage nomenclature has been changed to
!                             surface/target.
!     12/04/05    "    "      Added optional X and F to deal with radiative
!                             heating that goes to zero beyond some X for CEV.
!     10/04/06    "    "      Installed Tecplot 360 version of the I/O package.
!     10/30/06    "    "      Turned on verbose Tecplot I/O after a mystery
!                             showed up reading output from BLAYER.
!     11/01/06    "    "      Avoid blank zone titles for unlabeled target pts.
!     11/03/06    "    "      For the random points input case, write the
!                             plottable results as a single (n,1,1) zone rather
!                             than one zone per point.
!     03/26/07    "    "      Todd White asked for minimum distances to be added
!                             to the tabulation for the list-of-points case.
!     05/02/07    "    "      The mean squared minimum distance didn't have npts
!                             inside the square root.
!     02/05/14    "    "      All ADT variants are now in a single module with
!                             generic build_adt and search_adt interfaces.
!     05/20/22    "    "      Redesigned the control file to handle all of the
!                             likely options for input surface and target data.
!     05/31/22    "    "      Finished initial rewrite.  Many permutations to
!                             test.
!     06/02/22    "    "      All options seem to be working.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!           (Later with ERC, Inc. then with AMA, Inc. at NASA ARC.)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use tecplot_io_module      ! For structured Tecplot zones
   use grid_header_structure  ! See Tecplot_io.f90 for these two modules
   use grid_block_structure
   use xyzq_io_module         ! Plot3D-type file I/O
   use tri_header_structure   ! Data structures used by triangulation_io
   use tri_zone_structure
   use triangulation_io       ! For unstructured surfaces and volumes
   use adt_utilities          ! All ADT variants (generic build_adt, search_adt)

   implicit none

!  Constants:

   integer, parameter :: &
      lunsurface =  1,   &  ! Input surface solution [grid if Plot3D)
      lun_f      =  2,   &  ! Input and output function data if Plot3D
      luntarget  =  3,   &  ! Input target points, surface or (x,y,z) list
      lunout     =  4,   &  ! Output results in format of the target points
      lunctl     =  5,   &  ! Control file (standard input)
      luncrt     =  6,   &  ! For diagnostics
      luntab     =  7,   &  ! Output tabulation (mode 4)
      lunscratch =  8,   &  ! For RDLIST kludge (it can't read from a buffer)
      ndim       =  3,   &
      nf_limit   = 64       ! Limit on subset of functions to interpolate

   real, parameter ::    &
      one  = 1.0,        &
      zero = 0.0

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

   character (11), parameter :: &
      fileform = 'unformatted'

!  Variables:

   integer :: &
      ib, il, iorder, ios, iz, l, l1, len_buffer, m, mode_surface, &
      mode_target, ncells, nf_out, nf_surface, nf_target, ni, nj, nnodes, &
      npts, nz_surface_s, nz_surface_u, nz_target_s, nz_target_u

   integer, allocatable, dimension (:) :: &
      list_of_functions

   integer, allocatable, dimension (:,:) :: &
      conn               ! For zone & (i,j) of structured surface to search

   real :: &
      f_beyond_xub, x_upper_bound

   logical :: &
      cell_centered, formatted, interp_details, structured_surface, &
      structured_target, tabulation

   character (38) :: &
      tab_format
   character (80) :: &  ! Filename length limit matches I/O package limits
      buffer, filename_f, filename_interp, filename_surface, filename_target, &
      filename_tab

!  Derived data types:

!  (a) Structured datasets (_s), Plot3D or Tecplot format; also for target list:

   type (grid_header) :: &
      header_surface_s, header_target_s

   type (grid_type), pointer, dimension (:) :: &
      one_zone, surface_s, target_s

!  (b) Unstructured datasets (_u), triangulations or quad surfaces:

   type (tri_header_type) :: &
      header_surface_u, header_target_u

   type (tri_type), pointer, dimension (:) :: &
      surface_u, target_u

!  Execution:
!  !!!!!!!!!!

   call read_controls ()
   if (ios /= 0) go to 99

   call read_surface_data ()
   if (ios /= 0) go to 99

   call read_target_data ()
   if (ios /= 0) go to 99

   call build_search_tree ()
   if (ios /= 0) go to 99

   call interpolate_at_targets ()
   if (ios /= 0) go to 99

   call save_results ()
   if (ios /= 0) go to 99

99 continue  ! Avoid system-dependent STOP issues; let F90 deallocate everything


!  Local procedures for program SURFACE_INTERP in ~alphabetical order:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   contains

!     --------------------------------------------------------------------------
      subroutine mixed_conn (header)  ! Turn conn(1:nv,:) into conn(0:8,:)
                                      ! with cell type in conn(0,:)
!     --------------------------------------------------------------------------

!     Argument:

      type (tri_header_type), intent (inout) :: header  ! One or more zones have
                                                        ! been merged here
!     Local variables:

      integer              :: i, itype, ne, nv
      integer, allocatable :: conn_temp(:,:)

!     Execution:

      nv = header%nvertices
      ne = header%nelements

      allocate (conn_temp(nv,ne))

      conn_temp(:,:) = header%conn(:,:)

      deallocate (header%conn)
      allocate   (header%conn(0:8,ne)) ! Not 0:nv because the mixed cell form of
                                       ! ADT searching expects up to 8 vertices

      if (nv == 3) then  ! See subroutine search_mixed_cell_adt in adt_utilities
          itype = 1  ! Triangular cells
      else
          itype = 3  ! Quad cells
      end if

      header%conn(0,:) = itype
      do i = 1, ne
         header%conn(1:nv,i) = conn_temp(:,i)
      end do
      deallocate (conn_temp)

      end subroutine mixed_conn

!     --------------------------------------------------------------------------
      subroutine read_controls ()

!     SURFACE_INTERP control file
!     SRL-surf.gu      ! File to interp. in (*.g[u] if Plot3D; *.f[u] assumed)
!     1                ! 1|2|3 = Plot3D|Tecplot structured|Tecplot unstructured
!     target-surf.dat  ! Target file to interpolate to
!     2                ! 1|2|3|4 = Plot3D|Tec struct.|Tec unstruct.|(x,y,z) list
!     Subset of functions in output file(s):
!     1 2 3
!     SRL.interp.dat   ! Plottable results
!     SRL.interp.txt   ! Tabulated results (target mode 4 only)
!     F                ! Interpolation details at initial target point?
!     2                ! 1|2 = BLOCK|POINT order output
!     0.779347  0.     ! X & F to use for x > X
!     --------------------------------------------------------------------------

      integer :: line

      line = 1
      read (lunctl, *, iostat=ios)
      if (ios /= 0) go to 90

      line = line + 1
      read (lunctl, *, iostat=ios) filename_surface
      if (ios /= 0) go to 90
      if (len_trim (filename_surface) == len (filename_surface)) go to 80

      line = line + 1
      read (lunctl, *, iostat=ios) mode_surface
      if (ios /= 0) go to 90

      structured_surface = mode_surface <= 2

      line = line + 1
      read (lunctl, *, iostat=ios) filename_target
      if (ios /= 0) go to 90
      if (len_trim (filename_target) == len (filename_target)) go to 80

      line = line + 1
      read (lunctl, *, iostat=ios) mode_target
      if (ios /= 0) go to 90

      structured_target = mode_target /= 3
      tabulation        = mode_target == 4

!     Allow for outputting a subset of the interpolated surface functions.

      line = line + 1
      read (lunctl, *, iostat=ios)
      if (ios /= 0) go to 90

      line = line + 1
      read (lunctl, '(a)', iostat=ios) buffer
      if (ios /= 0) go to 90

!     Kludge because RDLIST can't read from a buffer:

      open   (lunscratch, status='scratch')
      write  (lunscratch, '(a)') trim (buffer)
      rewind (lunscratch)

      nf_out = nf_limit  ! Because the input surface hasn't been read yet
      allocate (list_of_functions(nf_out))

      call rdlist (luncrt, &
                   'Reading list of output function numbers. ', &
                   lunscratch, nf_out, list_of_functions)
      close  (lunscratch)

      write (luncrt, '(a, 64i3)') 'Functions to output:', &
         list_of_functions(1:nf_out)

      line = line + 1
      read (lunctl, *, iostat=ios) filename_interp
      if (ios /= 0) go to 90
      if (len_trim (filename_interp) == len (filename_interp)) go to 80

      line = line + 1
      read (lunctl, *, iostat=ios) filename_tab
      if (ios /= 0) go to 90
      if (len_trim (filename_tab) == len (filename_tab)) go to 80

      line = line + 1
      read (lunctl, *, iostat=ios) interp_details  ! T|F
      if (ios /= 0) go to 90

      line = line + 1
      read (lunctl, *, iostat=ios) iorder  ! 1|2 = BLOCK|POINT order output
      if (ios /= 0) go to 90

      line = line + 1
      read (lunctl, *, iostat=ios) x_upper_bound, f_beyond_xub
      if (ios /= 0) go to 90

      close (lunctl)
      go to 99

   80 write (luncrt, '(a, i4)')  '** File name limit encountered:', &
         len (filename_surface), '   Control input line:         ', line
      ios = 1;  go to 99

   90 write (luncrt, '(a, i3)') '** Trouble reading controls.  Line #:', line

   99 return

      end subroutine read_controls

!     --------------------------------------------------------------------------
      subroutine read_surface_data ()  ! The data to be interpolated within
!     --------------------------------------------------------------------------

      l = len_trim (filename_surface)

      select case (mode_surface)

         case (1)  ! Plot3D

            filename_f(1: ) = filename_surface(1:l)

            if (filename_surface(l:l) == 'g') then
                l1 = 3
                filename_f(l:l) = 'f'
            else if (filename_surface(l-1:l) == 'gu') then
                l1 = 1
                filename_f(l-1:l-1) = 'f'
            else
                write (luncrt, '(2a)') '** Unknown surface grid (g[u]): ', &
                    filename_surface(l:l)
                    ios = 1
                    go to 99
            end if

            open (lunsurface, file=filename_surface(1:l), &
                  form=fileform(l1:11), status='old', iostat=ios)
            if (ios /= 0) then
               write (luncrt, '(2a)') &
                  '** Cannot open surface grid file ', filename_surface(1:l)
                  ios = 1
                  go to 99
            end if

            open (lun_f, file=filename_f(1:l), &
                  form=fileform(l1:11), status='old', iostat=ios)
            if (ios /= 0) then
               write (luncrt, '(2a)') &
                  '** Cannot open surface function file ', filename_f(1:l)
                  ios = 1
                  go to 99
            end if

            formatted = l1 == 3

            call xyzq_read (lunsurface, lun_f, formatted, nz_surface_s, &
                            nf_surface, cell_centered, surface_s, ios)
            if (ios /= 0) then
               write (luncrt, '(a)') &
                  '** Trouble reading surface data in Plot3D form:', &
                  trim (filename_surface), trim (filename_f)
               go to 99
            end if

            header_surface_s%numq = nf_surface

         case (2)  ! Tecplot structured

            header_surface_s%filename = filename_surface
            header_surface_s%formatted = true
            header_surface_s%ndim      = ndim
            ios                        = 1    ! Verbose input mode

            call Tecplot_read (lunsurface, header_surface_s, surface_s, ios)
            if (ios /= 0) then
               write (luncrt, '(/, 2a)') &
                  ' Trouble reading Tecplot structured surface data file ', &
                  trim (filename_surface)
               go to 99
            end if

            nz_surface_s = header_surface_s%nblocks
            nf_surface   = header_surface_s%numq

         case (3)  ! Tecplot unstructured

!           Determine the element type without a prompt:

            header_surface_u%filename  = filename_surface
            header_surface_u%formatted = true

            call get_element_type (lunsurface, header_surface_u, ios)
            if (ios /= 0) go to 99

            write (luncrt, '(a, i3)') &
               ' # vertices per element found for surface zone 1:', &
               header_surface_u%nvertices, &
               ' Any further zones are assumed to be of the same element type.'

            header_surface_u%combine_zones         = true
            header_surface_u%centroids_to_vertices = true

            ios = 1  ! Verbose mode

            if (header_surface_u%nvertices == 3) then  ! Triangulation

               call tri_read (lunsurface, header_surface_u, surface_u, ios)
               if (ios /= 0) then
                  write (luncrt, '(2a)') &
                     '*** Trouble reading triangulated surface ', &
                     trim (filename_surface)
                  go to 99
               end if

            else  ! Quad surface, read as for a volume dataset

              call vol_read (lunsurface, header_surface_u, surface_u, ios)
               if (ios /= 0) then
                  write (luncrt, '(2a)') &
                     '*** Trouble reading unstructured surface ', &
                        trim (filename_surface)
                  go to 99
               end if

            end if

            nz_surface_u = header_surface_u%nzones
            nf_surface   = header_surface_u%numf

!           Use the mixed cell form of surface searching for both tri & quad:

            call mixed_conn (header_surface_u)  ! conn(1:nv,:) --> conn(0:nv,:)

         case default

            write (luncrt, '(a, i10)') &
               '** Bad surface dataset type (1|2|3):', mode_surface
            ios = 1
            go to 99

      end select

      write (luncrt, '(a, i3)') ' # functions in surface solution:', nf_surface

   99 return

      end subroutine read_surface_data

!     --------------------------------------------------------------------------
      subroutine read_target_data ()  ! (x,y,z)s, structured/unstructured/list
!     --------------------------------------------------------------------------

!     We read the target points to interpolate at and allocate the interpolated
!     fn. values so that the output has the same format as the target input.

      l = len_trim (filename_target)

      select case (mode_target)

         case (1)  ! Plot3D

           if (filename_target(l:l) == 'g') then
                l1 = 3
                filename_f(l:l) = 'f'
            else if (filename_target(l-1:l) == 'gu') then
                l1 = 1
                filename_f(l-1:l-1) = 'f'
            else
                write (luncrt, '(2a)') '** Unknown target grid (g[u]): ', &
                    filename_target(l:l)
                    ios = 1
                    go to 99
            end if

            open (luntarget, file=filename_target(1:l), &
                  form=fileform(l1:11), status='old', iostat=ios)
            if (ios /= 0) then
               write (luncrt, '(2a)') &
                  '** Cannot open target grid file ', filename_target(1:l)
                  ios = 1
                  go to 99
            end if

            formatted = l1 == 3

            call xyzq_read (luntarget, -luntarget, formatted, nz_target_s, &
                            nf_target, cell_centered, target_s, ios)
            if (ios /= 0) then
               write (luncrt, '(a)') &
                  '** Trouble reading target data in Plot3D form:', &
                  trim (filename_target)
               go to 99
            end if

!           Allocate space for the interpolated function outputs:

            do iz = 1, nz_target_s
               target_s(iz)%mi = target_s(iz)%ni
               target_s(iz)%mj = target_s(iz)%nj
               target_s(iz)%mk = target_s(iz)%nk
               call q_allocate (target_s(iz), nf_out, ios)
               if (ios /= 0) then
                  write (luncrt, '(a,i6)') &
                     '** Trouble allocating interpolated fns. Block #:', iz
                  go to 99
               end if
            end do

         case (2)  ! Tecplot structured

            header_target_s%filename  = filename_target
            header_target_s%formatted = true
            header_target_s%ndim      = ndim

            call Tecplot_read (luntarget, header_target_s, target_s, ios)
            if (ios /= 0) then
               write (luncrt, '(2a)') &
                  '*** Trouble reading target Tecplot structured surface ', &
                  trim (filename_target)
               go to 99
            end if

            deallocate (header_target_s%varname)
            allocate (header_target_s%varname(ndim+nf_out))
            if (structured_surface) then
              header_target_s%varname(1:ndim) = header_surface_s%varname(1:ndim)
            else
              header_target_s%varname(1:ndim) = header_surface_u%varname(1:ndim)
            end if

            do m = 1, nf_out
               header_target_s%varname(ndim+m) = &
                  header_surface_s%varname(ndim+list_of_functions(m))
            end do

            nz_target_s = header_target_s%nblocks

            do iz = 1, nz_target_s
               ni = target_s(iz)%ni
               nj = target_s(iz)%nj
               if (header_target_s%numq > 0) deallocate (target_s(iz)%q)
               header_target_s%numq = nf_out
               allocate (target_s(iz)%q(nf_out,ni,nj,1))
            end do

         case (3)  ! Tecplot unstructured

!           Determine the element type without a prompt:

            header_target_u%filename  = filename_target
            header_target_u%formatted = true

            call get_element_type (luntarget, header_target_u, ios)
            if (ios /= 0) go to 99

            write (luncrt, '(a, i3)') &
               ' # vertices per element found for target surface zone 1:', &
               header_target_u%nvertices, &
               ' Any further zones are assumed to be of the same element type.'

            header_target_u%combine_zones = false

            ios = 1  ! Verbose mode

            if (header_target_u%nvertices == 3) then  ! Triangulation

               call tri_read (luntarget, header_target_u, target_u, ios)
               if (ios /= 0) then
                  write (luncrt, '(2a)') &
                     '*** Trouble reading target triangulated surface ', &
                     trim (filename_target)
                  go to 99
               end if

            else  ! Quad. surface, read as a volume

               call vol_read (luntarget, header_target_u, target_u, ios)
               if (ios /= 0) then
                  write (luncrt, '(2a)') &
                     '*** Trouble reading unstructured target surface ', &
                        trim (filename_target)
                  go to 99
               end if

            end if

            nz_target_u = header_target_u%nzones

            do iz = 1, nz_target_u
               ni = target_u(iz)%nnodes
               if (header_target_u%numf > 0) deallocate (target_u(iz)%f)
               header_target_u%numf = nf_out
               allocate (target_u(iz)%f(nf_out,ni))
            end do

            header_target_u%fileform = 1  ! Vertex-centered

         case (4)  ! List of (x,y,z)s

            call read_target_point_list ()

         case default

            write (luncrt, '(a, i10)') &
               '** Bad target dataset type (1|2|3|4):', mode_target
            ios = 1

      end select

   99 return

      end subroutine read_target_data

!     --------------------------------------------------------------------------
      subroutine read_target_point_list ()

!     Count and read surface points listed one per line with an optional label.
!     Treat them as one or more (1,1,1) structured Tecplotable surface zones.
!     --------------------------------------------------------------------------

!     Local constants:

      character (2), parameter :: delimiters = ' ,'

!     Local variables:

      integer       :: first, i, iz, last, mark, nchar
      real          :: xyz(ndim)
      character (7) :: format_string

!     Execution:

      open (luntarget, file=filename_target, status='old', iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(/, 2a)') &
            ' Unable to open target file ', trim (filename_target)
         go to 99
      end if

      nz_target_s = 0
      do while (ios == 0)
         nz_target_s = nz_target_s + 1
         read (luntarget, '(a)', iostat=ios)
      end do
      nz_target_s = nz_target_s - 1

      rewind (luntarget)

      header_target_s%ndim        = ndim
      header_target_s%numq        = 0
      header_target_s%nblocks     = nz_target_s
      header_target_s%ndatasetaux = 0

      allocate (target_s(nz_target_s))

      len_buffer = len (buffer);  format_string = '(fnn.0)'

      do iz = 1, nz_target_s

         target_s(iz)%ni = 1;  target_s(iz)%nj = 1;  target_s(iz)%nk = 1

         call Tec_block_allocate (target_s(iz), ndim, header_target_s%numq, ios)

         read (luntarget, '(a)') buffer

!        Parse the line in order to handle descriptors with embedded blanks:

         first = 1;  last = len_buffer

         do i = 1, ndim

            call scan2 (buffer, delimiters, first, last, mark)

            nchar = mark - first + 1
            if (nchar < 10) then
               write (format_string(2:4), '(a2, i1)') ' f', nchar
            else
               write (format_string(2:4), '(a1, i2)')  'f', nchar
            end if

            read (buffer(first:mark), format_string, iostat=ios) xyz(i)

            if (ios /= 0) then
               write (luncrt, '(/, a, i5)') &
                  ' Trouble reading target (x,y,z).  Line #:', iz
               go to 99
            end if

!           In case the points aren't labeled, enter x as the zone title:

            if (i == 1) target_s(iz)%zone_title = buffer(first:mark)

            first = mark + 2
         end do

         target_s(iz)%x = xyz(1)
         target_s(iz)%y = xyz(2)
         target_s(iz)%z = xyz(3)

         if (first <= last) target_s(iz)%zone_title = buffer(first:last)

         target_s(iz)%nzoneaux     = 0
         target_s(iz)%solutiontime = -999.  ! Undefined

      end do

 99   continue

      close (luntarget)

      end subroutine read_target_point_list

!     --------------------------------------------------------------------------
      subroutine build_search_tree ()
!     --------------------------------------------------------------------------

      select case (mode_surface)

         case (1, 2)  ! Plot3D or structured Tecplot zone(s)

!           Count the quad cells of all surface zones:

            ncells = 0
            do iz = 1, nz_surface_s
               ncells = (surface_s(iz)%ni - 1) * (surface_s(iz)%nj - 1) + ncells
            end do

!           Make use of the unstructured surface connectivity array that
!           won't be used otherwise in this run:

            allocate (header_surface_u%conn(3,ncells)) ! For patch # and (i,j)

            call build_adt (nz_surface_s, surface_s, ncells, &
                            header_surface_u%conn)

         case (3)  ! Tecplot triangulation or unstructured quad surface

            nnodes = header_surface_u%nnodes
            ncells = header_surface_u%nelements

            call build_adt (nnodes, ncells, header_surface_u%conn, &
                            header_surface_u%xyz, ios)
            if (ios /= 0) then
               write (luncrt, '(a)') &
                  '** Trouble building mixed cell search tree.'
               go to 99
            end if

         case default

            write (luncrt, '(a, i10)') &
               '** Build_search_tree: Bad surface dataset type (1|2|3):', &
               mode_surface
            ios = 1
            go to 99

      end select

      ios = 0

   99 return

      end subroutine build_search_tree

!     --------------------------------------------------------------------------
      subroutine interpolate_at_targets ()

!     Locate the target surface points in the surface solution grid and
!     interpolate all the flow quantities.  Return the specified subset.
!     The surface may be structured or unstructured; likewise for the target.
!     --------------------------------------------------------------------------

!     Local variables:

      integer :: i, ib, ib1, ic, ier, icell, iz, j, jc, n, ninside, noutside, &
                 nv, nz_target
      real    :: coefs(8), dmax, dmean, dmin, dsqmin, dtolsq, fv(4), p, pm1,  &
                 q, qm1, interp_xyz(3), target_xyz(3)
      real, allocatable :: interp_f(:)

!     Execution:

      allocate (interp_f(nf_surface))
      dtolsq = (0.001) ** 2  ! Tolerance for search diagnostics (refine later)

      if (.not. structured_surface) nv = header_surface_u%nvertices  ! 3|4
      if (structured_target) then
         nz_target = nz_target_s
      else
         nz_target = nz_target_u
      end if

!!!   write (luncrt, '(a, i6)') 'nz_target:', nz_target

      do iz = 1, nz_target

         ninside  = 0;  dmax  = zero
         noutside = 0;  dmean = zero
         npts     = 0

         if (structured_target) then
            ni = target_s(iz)%ni
            nj = target_s(iz)%nj
            allocate (target_s(iz)%q(nf_out,ni,nj,1))
         else
            ni = target_u(iz)%nnodes
            nj = 1
            allocate (target_u(iz)%f(nf_out,ni))
         end if

!! ??    results(iz)%xmin = zero  ! Avoid possible undefined minimum distance

!!!      write (luncrt, '(a, 3i6)') 'iz,ni,nj:', iz,ni,nj

         do j = 1, nj

            do i = 1, ni

!!!            write (luncrt, '(a, 3i6)') 'iz,i,j:', iz,i,j

               if (structured_target) then
                  target_xyz(1) = target_s(iz)%x(i,j,1)
                  target_xyz(2) = target_s(iz)%y(i,j,1)
                  target_xyz(3) = target_s(iz)%z(i,j,1)
               else
                  target_xyz(:) = target_u(iz)%xyz(:,i)
               end if

!              The following kludge was needed for CEV radiative heating on
!              aft-body patches of a full CFD surface grid:

               if (target_xyz(1) > x_upper_bound) then

                  if (structured_target) then
                     target_s(iz)%q(:,i,j,1) = f_beyond_xub
                  else
                     target_u(iz)%f(:,i) = f_beyond_xub
                  end if
                  cycle

               end if

               if (structured_surface) then

                  call search_adt (target_xyz, icell, p, q, dsqmin, true, &
                                   nz_surface_s, surface_s, ncells, &
                                   header_surface_u%conn, interp_xyz)
               else

                  call search_adt (target_xyz, icell, coefs, dsqmin, true,  &
                                   nnodes, ncells, header_surface_u%conn, &
                                   header_surface_u%xyz, interp_xyz, ios)
                  if (ios /= 0) then
                     write (luncrt, '(a, 3i16)') &
                        '** Search trouble.  Target zone, i, j:', iz, i, j
                     write (luncrt, '(a, 3es16.8)') '   Target xyz:', target_xyz
                     go to 99
                  end if

               end if

               if (dsqmin < dtolsq) then ! The nearest cell was within tolerance
                  ninside  = ninside + 1
               else
                  noutside = noutside + 1
               end if
               npts = npts + 1

               dmax  = max (dmax, dsqmin)
               dmean = dmean + dsqmin

!              Interpolate the surface flow at this target point:

               if (structured_surface) then  ! We're using the unstruct. conn
                  ib  = header_surface_u%conn(1,icell) ! Zone or block #
                  ic  = header_surface_u%conn(2,icell) ! Lower left cell indices
                  jc  = header_surface_u%conn(3,icell)
                  pm1 = one - p
                  qm1 = one - q

                  interp_f(:) = &
                     qm1 * (pm1 * surface_s(ib)%q(:,ic,  jc,  1)  + &
                              p * surface_s(ib)%q(:,ic+1,jc,  1)) + &
                       q * (pm1 * surface_s(ib)%q(:,ic,  jc+1,1)  + &
                              p * surface_s(ib)%q(:,ic+1,jc+1,1))
               else
                  do n = 1, nf_surface
                     fv(1:nv) = &
                        header_surface_u%f(n,header_surface_u%conn(1:nv,icell))
                     interp_f(n) = dot_product (coefs(1:nv), fv(1:nv))
                  end do
               end if

!              Debug output at first target point?

               if (interp_details) then
                  if (iz == 1 .and. i == 1 .and. j == 1) then
                     write (luncrt, '(a, 3f20.8)') ' Target xyz: ', target_xyz
                     write (luncrt, '(a, 3f20.8)') ' Interp xyz: ', interp_xyz
                     write (luncrt, '(a)') ' Interpolated values:'
                     write (luncrt, '(i3, es16.8)') &
                        (n, interp_f(n), n = 1, nf_surface)
                     dmin = sqrt (dsqmin)
                     write (luncrt, '(a, f20.10)') ' Shortest distance: ', dmin
                     if (structured_surface) then
                        write (luncrt, '(a, 3i5)') &
                           ' Search result: block, ic, jc: ', ib, ic, jc
                        write (luncrt, '(a, 2f16.10)') ' p, q: ', p, q
                        write (luncrt, '(a)') ' Vertex f(:,i:i+1,j:j+1):'
                        do n = 1, nf_surface
                           write (luncrt, '(i3, 4es16.8)') &
                              n, surface_s(ib)%q(n,ic:ic+1,jc:jc+1,1)
                        end do
                     else
                        write (luncrt, '(a, i10)') ' Best cell:', icell
                        write (luncrt, '(a, 4f16.10)') ' Coefs:', coefs(1:nv)
                        write (luncrt, '(a)') ' Vertex f(1:nv):'
                        do n = 1, nf_surface
                           write (luncrt, '(i3, 4es16.8)') n, &
                         header_surface_u%f(n,header_surface_u%conn(1:nv,icell))
                        end do
                     end if
                  end if
               end if

!              Store the indicated subset of functions:

               if (structured_target) then
                  do n = 1, nf_out
                     target_s(iz)%q(n,i,j,1) = interp_f(list_of_functions(n))
                  end do

!                 Enhance the tabulation from the list-of-points case?

                  if (tabulation) then                  ! Only one pt. per block
                     target_s(iz)%xmin = sqrt (dsqmin)  ! so use this field
                  end if
               else
                  do n = 1, nf_out
                     target_u(iz)%f(n,i) = interp_f(list_of_functions(n))
                  end do
               end if

            end do ! Next i for this target zone

         end do ! Next j for this target zone

         npts = max (npts, 1) ! In case entire patch was beyond x_upper_bound

         write (luncrt, '(a, i5, a, 2i7, a, 2es12.5)')                  &
            '  Target zone', iz,                                        &
            ':  # points in/outside tolerance:', ninside, noutside,     &
            ';  max/mean devn.:', sqrt (dmax), sqrt (dmean / real (npts))

      end do ! Next target zone

  99  return

      end subroutine interpolate_at_targets

!     --------------------------------------------------------------------------
      subroutine save_results ()  ! In the same format as the target dataset
!     --------------------------------------------------------------------------

      integer :: n

      if (tabulation) then  ! Simple tabulation; %xmin = min. distance

         open (luntab, file=filename_tab, status='unknown')

         tab_format = '(3es13.5, 2x, nnes14.6, 2x, a, es14.6)'
         write (tab_format(15:16), '(i2)') nf_out

         write (luntab, tab_format) &
         (target_s(ib)%x(1,1,1), target_s(ib)%y(1,1,1), target_s(ib)%z(1,1,1), &
          target_s(ib)%q(:,1,1,1), trim (target_s(ib)%zone_title),             &
          target_s(ib)%xmin, ib = 1, nz_target_s)

         close (luntab)

!        Save the tabulation in plottable form too:

         allocate (header_target_s%varname(ndim+nf_out))

         if (mode_surface == 1) then  ! Input *.f[u] has no function names
            do n = 1, header_target_s%numq  ! Assume nf < 100
               header_target_s%varname(ndim+n) = 'f'  ! Blanks the rest
               write (header_target_s%varname(ndim+n)(2:3), '(i2)') n
            end do
         end if

         header_target_s%filename    = filename_interp
         header_target_s%formatted   = true
         header_target_s%ndim        = ndim
         header_target_s%numq        = nf_out
         header_target_s%nblocks     = nz_target_s
         header_target_s%datapacking = 0  ! Point order
         header_target_s%title       = filename_interp
         header_target_s%ndatasetaux = 0
         header_target_s%varname(1)  = 'x'
         header_target_s%varname(2)  = 'y'
         header_target_s%varname(3)  = 'z'

         do n = 1, nf_out
            if (mode_surface <= 2) then
               header_target_s%varname(ndim+n) = &
                  header_surface_s%varname(ndim+(list_of_functions(n)))
            else
               header_target_s%varname(ndim+n) = &
                  header_surface_u%varname(ndim+(list_of_functions(n)))
            end if
         end do

!        One point per zone isn't much good if there are lots of them.
!        Merge them all into one zone:

         allocate (one_zone(1))

         one_zone(1)%zone_title   = 'Target points'
         one_zone(1)%nzoneaux     = 0
         one_zone(1)%solutiontime = -999.
         one_zone(1)%ni           = nz_target_s
         one_zone(1)%nj           = 1
         one_zone(1)%nk           = 1
         one_zone(1)%mi           = nz_target_s
         one_zone(1)%mj           = 1
         one_zone(1)%mk           = 1

         call Tec_block_allocate (one_zone(1), ndim, nf_out, ios)

         do iz = 1, nz_target_s
            one_zone(1)%x(iz,1,1)   = target_s(iz)%x(1,1,1)
            one_zone(1)%y(iz,1,1)   = target_s(iz)%y(1,1,1)
            one_zone(1)%z(iz,1,1)   = target_s(iz)%z(1,1,1)
            one_zone(1)%q(:,iz,1,1) = target_s(iz)%q(:,1,1,1)
         end do

         header_target_s%nblocks = 1

         call Tecplot_write (lunout, header_target_s, one_zone, ios)
         if (ios /= 0) then
            write (luncrt, '(2a)') &
              '** Trouble writing tabulated results in Tecplotable form):', &
              trim (filename_interp), trim (filename_interp)
         end if

      end if

      select case (mode_target)

         case (1)  ! Plot3D

            l = len_trim (filename_interp)
            filename_f(1: ) = filename_interp(1:l)

            if (filename_interp(l:l) == 'f') then
                l1 = 3
            else if (filename_interp(l-1:l) == 'fu') then
                l1 = 1
            else
                write (luncrt, '(2a)') &
                   '** Invalid interpolated function file name (f[u]): ', &
                   filename_interp(l:l)
                   ios = 1
                   go to 99
            end if

            open (lunout, file=filename_interp(1:l), &
                  form=fileform(l1:11), status='unknown', iostat=ios)
            if (ios /= 0) then
               write (luncrt, '(2a)') &
                  '** Cannot open output function file ', filename_interp(1:l)
                  ios = 1
                  go to 99
            end if

            formatted = l1 == 3

            call q_write (lunout, formatted, nz_target_s, nf_out, target_s, ios)
            if (ios /= 0) then
               write (luncrt, '(a)') &
                  '** Trouble writing interpolated Plot3D function file:', &
                  trim (filename_interp), filename_interp(1:l)
            end if

         case (2)  ! Tecplot structured

            header_target_s%filename = filename_interp

            call Tecplot_write (lunout, header_target_s, target_s, ios)
            if (ios /= 0) then
               write (luncrt, '(2a)') &
                 '** Trouble writing interpolated results file (structured):', &
                 trim (filename_interp), trim (filename_interp)
            end if

         case (3)  ! Tecplot unstructured

            header_target_u%filename = filename_interp

            if (header_target_u%nvertices == 3) then  ! Triangulation

               call tri_write (lunout, header_target_u, target_u, ios)
               if (ios /= 0) then
                  write (luncrt, '(2a)') &
                     '*** Trouble writing interpolated triangulation ', &
                     trim (filename_interp)
               end if

            else  ! Quad. surface, written as a volume

               call vol_write (lunout, header_target_u, target_u, ios)
               if (ios /= 0) then
                  write (luncrt, '(2a)') &
                     '*** Trouble writing interpolated unstructured results ', &
                        trim (filename_interp)
                  go to 99
               end if

            end if

         end select

   99 return

      end subroutine save_results

   end program surface_interp
