!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program surface_interp_2d
!
!  Description:
!
!     This is the 2D analogue of SURFACE_INTERP.
!
!     Read a multizone 2D surface dataset (Tecplot ASCII) and a second target
!  dataset of one or more zones (or possibly unstructured points in a simpler
!  format).  Interpolate the first dataset at each target point and save the
!  results as a Tecplotable file and (for the simpler input format) as a tab-
!  ulation for possible pasting into a spreadsheet.  The output file(s) may
!  contain just a subset of the input functions.
!
!  Control file format (standard input):
!
!     SURFACE_INTERP_2D control file
!     body_points.inp        ! Target file name
!     1                      ! Target format: 1 = list of points; 2 = Tecplot
!     Subset of functions in output file(s):
!     2 6
!     cssr-fullbody-wall.dat ! File to be interpolated
!     body_point_pw_qw.dat   ! Single-zone plottable results
!     body_point_pw_qw.txt   ! Tabulated results
!     1 1                    ! Target (zone, i) for debug output; 0s suppress
!    [BLOCK                  ! BLOCK | POINT overrides default of POINT output]
!    [7.890     0.           ! Optional X & F; if x > X, set finterp(:) = F]
!    [1. 1.e-4               ! Optional scale factors for output functions]
!
!  Random target points format (mode 1):
!
!     x1   y1   [Stag point 1]       ! Any string may follow the coordinates
!     x2   y2   [Shoulder point 2]   ! else x is inserted as a tag
!     x3   y3   [Aft stag point 3]
!     :    :    [:    :]
!
!  Tabulated output format (mode 1):
!
!     x1   y1   f21  f61  Stag point 1 (or x1)     Min. distance
!     x2   y2   f22  f62  Shoulder point 2 (or x2)    "      "
!     x3   y3   f23  f63  Aft stag point 3 (or x3)    "      "
!     :    :    :    :
!
!  07/31/05  David Saunders  Initial implementation of 3-space SURFACE_INTERP.
!  02/05/14    "      "      Variant of SURFACE_INTERP from which the 2-space
!                            version is being adapted.  (All ADT variants had
!                            been placed in a single module with generic
!                            BUILD_ADT and SEARCH_ADT interfaces.)
!  03/15/15    "      "      SURFACE_INTERP_2D adapted from SURFACE_INTERP.
!                            (Two-space multiblock curve utilities have just
!                            been added to the ADT package.)
!                            W/m^2 data suggested an option to scale results.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_header_structure  ! See Tecplot_io.f90 for these modules
   use grid_block_structure
   use tecplot_io_module
   use adt_utilities          ! All ADT variants (generic build_adt, search_adt)

   implicit none

!  Constants:

   integer, parameter :: &
      luntarget  =  1,   &  ! Input target pts.|2D surface  (list|Tecplot ASCII)
      lunsurface =  2,   &  ! Input structured 2D surface file   (Tecplot ASCII)
      lunout     =  3,   &  ! Output plottable results           (Tecplot ASCII)
      luntab     =  4,   &  ! Output tabulation (mode 1)
      lunctl     =  5,   &  ! Control file (standard input)
      luncrt     =  6,   &  ! For diagnostics
      lunscratch =  7,   &  ! For RDLIST kludge (it can't read from a buffer)
      ndim       =  2

   real, parameter ::    &
      one  = 1.0,        &
      zero = 0.0

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

!  Variables:

   integer :: &
      ib, idebug, ios, len_buffer, m, nb_surface, nb_target, ndebug,   &
      nf_out, nf_surface, npts, ncell, target_type

   integer, allocatable, dimension (:) :: &
      list_of_functions

   integer, allocatable, dimension (:,:) :: &
      conn               ! For zone & i index of each 2-pt. surface cell

   real :: &
      f_beyond_xub, x_upper_bound

   real, allocatable, dimension (:) :: &
      scale_factor

   logical :: &
      unstructured_targets

   character :: &
      buffer * 80, filename_tab * 80, output_packing * 5, tab_format * 38

!  Derived data types:

   type (grid_header) :: &
      header_one_zone, header_target, header_surface, header_results

   type (grid_type), pointer, dimension (:) :: &
      target, surface, results, &   ! (x,y,f) for the input surfaces & results
      one_zone                      ! Plottable results repacked for the
                                    ! unstructured input targets case

!  Execution:
!  !!!!!!!!!!

   read (lunctl, *)
   read (lunctl, *) header_target%filename
   read (lunctl, *) target_type

   len_buffer = len (buffer)        ! Used here and in read_target_points

   select case (target_type)

      case (1) ! Random list of points

         call read_target_points ()  ! Local procedure below

      case (2) ! Tecplot structured surface patch(es):

         header_target%formatted = true
         header_target%ndim      = ndim

         call Tecplot_read (luntarget, header_target, target, ios)

         nb_target = header_target%nblocks

      case default

         write (luncrt, '(/, a, i9)') &
            ' Bad target format #: ', target_type, &
            ' 1 = simple list of points; 2 = Tecplot file.'
         go to 99

   end select

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Trouble reading target points file ',        &
         trim (header_target%filename)
      go to 99
   end if

!  Allow for outputting a subset of the interpolated surface functions.
!  Buffer the function list so we can allocate the right size list later.

   read (lunctl, *)
   read (lunctl, '(a)') buffer

!  Read the surface solution to be interpolated:

   read (lunctl, *) header_surface%filename
   header_surface%formatted = true
   header_surface%ndim      = ndim
   ios                      = 1    ! Verbose input mode

   call Tecplot_read (lunsurface, header_surface, surface, ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Trouble reading surface data Tecplot file ', &
         trim (header_surface%filename)
      go to 99
   end if

   nb_surface = header_surface%nblocks
   nf_surface = header_surface%numq

   write (luncrt, '(a, i3)') ' # functions in surface solution:', nf_surface

   read (lunctl, *) header_results%filename
   read (lunctl, *) filename_tab
   read (lunctl, *) ndebug, idebug  ! ndebug <= 0 suppresses debug output

!  Optional controls:

   x_upper_bound = 1.e+30;  f_beyond_xub = zero   ! This f will not be used
   header_results%datapacking = 0                 ! POINT is the default

   read (lunctl, *, iostat=ios) output_packing    ! POINT or BLOCK

   allocate (scale_factor(nf_surface))            ! Probably more than enough

   if (ios == 0) then
      if (output_packing == 'BLOCK' .or.  &
          output_packing == 'block') header_results%datapacking = 1
  
      read (lunctl, *, iostat=ios) x_upper_bound, f_beyond_xub
   end if

!  Don't try to read scale factors till we know how many functions are output.

!! close (lunctl)

!  Kludge because RDLIST can't read from a buffer:

   open   (lunscratch, status='scratch')
   write  (lunscratch, '(a)') trim (buffer)
   rewind (lunscratch)

   nf_out = nf_surface
   allocate (list_of_functions(nf_out))

   call rdlist (luncrt, &
                'Reading list of output function numbers. ', &
                lunscratch, nf_out, list_of_functions)
   close  (lunscratch)

   write (luncrt, '(/, a, i4)') ' Number of output functions: ', nf_out

   scale_factor(:) = one
   if (ios == 0) read (lunctl, *, iostat=ios) scale_factor(1:nf_out)

   close (lunctl)

!  Special tabulation if target is really a list of body points (unstructured):

   unstructured_targets = true
   do ib = 1, nb_target
      if (target(ib)%ni > 1) unstructured_targets = false
   end do

!  Allocate interpolated results as though they are structured (may be 1 x 1):

   allocate (results(nb_target))

   do ib = 1, nb_target

      results(ib)%ni = target(ib)%ni
      results(ib)%nj = 1
      results(ib)%nk = 1
      results(ib)%mi = target(ib)%ni
      results(ib)%mj = 1
      results(ib)%mk = 1
      results(ib)%zone_title = target(ib)%zone_title
      results(ib)%nzoneaux   = target(ib)%nzoneaux

      call Tec_block_allocate (results(ib), ndim, nf_out, ios)

      results(ib)%x  = target(ib)%x
      results(ib)%y  = target(ib)%y

      if (results(ib)%nzoneaux > 0) then
          results(ib)%zoneauxname(:)  = target(ib)%zoneauxname(:)
          results(ib)%zoneauxvalue(:) = target(ib)%zoneauxvalue(:)
      end if

   end do

!  Construct a search tree from all patches of the surface data:

   ncell = 0
   do ib = 1, nb_surface
      ncell = (surface(ib)%ni - 1) + ncell
      allocate (surface(ib)%z(surface(ib)%ni,1,1))  ! For 3-sp. bounding boxes
   end do

   allocate (conn(2,ncell)) ! For patch # and (i,j)

   call build_adt (nb_surface, surface, ncell, 1, conn)


!  Perform the interpolations:
!  ---------------------------

   call interpolate_surface ()


!  Save results:

   if (unstructured_targets) then  ! Simple tabulation; %xmin = min. distance

      open (luntab, file=filename_tab, status='unknown')

      tab_format = '(2es13.5, 2x, nnes14.6, 2x, a, es14.6)'
      write (tab_format(15:16), '(i2)') nf_out

      write (luntab, tab_format) &
         (results(ib)%x(1,1,1), results(ib)%y(1,1,1), results(ib)%q(:,1,1,1),  &
          trim (results(ib)%zone_title), results(ib)%xmin, ib = 1, nb_target)

      close (luntab)

   end if

   allocate (header_results%varname(ndim + nf_out))

   header_results%varname(1:ndim) = header_surface%varname(1:ndim)

   do m = 1, nf_out
      header_results%varname(ndim + m) = &
      header_surface%varname(ndim + list_of_functions(m))
   end do

!  Save results in specified order, formatted:

   header_results%formatted   = true
   header_results%ndim        = ndim
   header_results%numq        = nf_out
   header_results%nblocks     = header_target%nblocks
   header_results%title       = 'Interpolated surface values'
   m                          = header_target%ndatasetaux
   header_results%ndatasetaux = m

   if (m > 0) then
      allocate (header_results%datasetauxname(m))
      header_results%datasetauxname(:)  = header_target%datasetauxname(:)
      allocate (header_results%datasetauxvalue(m))
      header_results%datasetauxvalue(:) = header_target%datasetauxvalue(:)
   end if

!  Repack results as a single zone?

   if (unstructured_targets) then

      header_one_zone%filename    = header_results%filename
      header_one_zone%formatted   = true
      header_one_zone%ndim        = ndim
      header_one_zone%numq        = nf_out
      header_one_zone%nblocks     = 1
      header_one_zone%datapacking = 0
      header_one_zone%title       = header_results%title
      header_one_zone%ndatasetaux = m

      if (m > 0) then
         allocate (header_one_zone%datasetauxname(m))
         header_one_zone%datasetauxname(:)  = header_target%datasetauxname(:)
         allocate (header_one_zone%datasetauxvalue(m))
         header_one_zone%datasetauxvalue(:) = header_target%datasetauxvalue(:)
      end if

      allocate (header_one_zone%varname(ndim + nf_out))

      do m = 1, ndim + nf_out
         header_one_zone%varname(m) = header_results%varname(m)
      end do

      allocate (one_zone(1))

      one_zone(1)%zone_title   = 'Target points'
      one_zone(1)%nzoneaux     = 0
      one_zone(1)%solutiontime = -999.
      one_zone(1)%ni           = nb_target
      one_zone(1)%nj           = 1
      one_zone(1)%nk           = 1
      one_zone(1)%mi           = nb_target
      one_zone(1)%mj           = 1
      one_zone(1)%mk           = 1

      call Tec_block_allocate (one_zone(1), ndim, nf_out, ios)

      do ib = 1, nb_target
         one_zone(1)%x(ib,1,1)   = results(ib)%x(1,1,1)
         one_zone(1)%y(ib,1,1)   = results(ib)%y(1,1,1)
         one_zone(1)%q(:,ib,1,1) = results(ib)%q(:,1,1,1)
      end do

      call Tecplot_write (lunout, header_one_zone, one_zone, ios)

   else

      call Tecplot_write (lunout, header_results, results, ios)

   end if

   if (ios /= 0) then
      write (luncrt, '(/, a)') &
         ' Trouble writing the interpolated results in Tecplot form.'
      stop
   end if

   close (lunout)

99 continue  ! Avoid system-dependent STOP issues; let F90 deallocate everything


!  Local procedures for program SURFACE_INTERP_2D
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   contains

!     --------------------------------------------------------------------------
      subroutine interpolate_surface ()

!     Locate the target surface points in the surface solution grid and
!     interpolate all the flow quantities.  Return the specified subset.
!     --------------------------------------------------------------------------

!     Local variables:

      integer :: i, ib, ib1, ic, ier, icell, m, n, ninside, noutside
      real    :: dmax, dmean, dmin, dsqmin, dtolsq, p, q, &
                 interp_xy(2), target_xy(2)
      real, allocatable :: interp_f(:)

!     Execution:

!     Tolerance for search diagnostics (refine later):

      dtolsq = (0.001) ** 2

      allocate (interp_f(nf_surface))

      do ib = 1, nb_target

         ninside  = 0;  dmax  = zero
         noutside = 0;  dmean = zero
         npts     = 0

         results(ib)%xmin = zero  ! Avoid possible undefined minimum distance

         do i = 1, target(ib)%ni

            target_xy(1) = target(ib)%x(i,1,1)
            target_xy(2) = target(ib)%y(i,1,1)

!           The following kludge was needed for CEV radiative heating on
!           aft-body patches of a full CFD surface grid:

            if (target_xy(1) > x_upper_bound) then

               results(ib)%q(:,i,1,1) = f_beyond_xub
               cycle

            end if

            call search_adt (target_xy, icell, p, q, dsqmin, true,         &
                             nb_surface, surface, ncell, 1, conn, interp_xy)

            if (dsqmin < dtolsq) then ! The nearest cell was within tolerance
               ninside  = ninside + 1
            else
               noutside = noutside + 1
            end if
            npts = npts + 1

            dmax  = max (dmax, dsqmin)
            dmean = dmean + dsqmin

!           Linearly interpolate the surface data at this target point:

            n  = conn(1,icell) ! Block #
            ic = conn(2,icell) ! Lower left cell indices

            interp_f(:) = p * surface(n)%q(:,ic,  1,1)  + &
                          q * surface(n)%q(:,ic+1,1,1)

!           Debug output:

            if (ib == ndebug) then
               if (i == idebug) then
                  write (6, '(a, 2f20.8)') ' target xy: ', target_xy
                  write (6, '(a, 2f20.8)') ' interp xy: ', interp_xy
                  write (6, '(a, 2i5)') &
                     ' Search result: block, ic: ', n, ic
                  write (6, '(a, 2f20.10)') ' p, q: ', p, q
                  dmin = sqrt (dsqmin)
                  write (6, '(a, f20.10)') ' shortest distance: ', dmin
                  write (6, '(/, a)') ' Cell vertex values:', &
                     ' fn #                ic                ic+1'
                  write (6, '(i3, 2f20.8)') &
                    (m, surface(n)%q(m,ic:ic+1,1,1), m = 1, nf_surface)
                  write (6, '(/, a)') ' Interpolated values:'
                  write (6, '(i5, f20.8)') (m, interp_f(m), m = 1, nf_surface)
               end if
            end if

!           Store the indicated subset of functions:

            do m = 1, nf_out
               results(ib)%q(m,i,1,1) = interp_f(list_of_functions(m)) * &
                                        scale_factor(m)
            end do

!           Kludge for enhancing the tabulation from the list-of-points case:

            if (unstructured_targets) then       ! Only one point per block,
               results(ib)%xmin = sqrt (dsqmin)  ! so make use of this field
            end if

         end do ! Next i for this target block

         npts = max (npts, 1) ! In case entire patch was beyond x_upper_bound

         write (luncrt, '(a, i5, a, 2i5, a, 2es12.5)') &
            ' Target zone or point', ib, &
            ':  # points in/outside tolerance:', ninside, noutside, &
            ';  max/mean devn.:', sqrt (dmax), sqrt (dmean / real (npts))

      end do ! Next target block

      end subroutine interpolate_surface

!     --------------------------------------------------------------------------
      subroutine read_target_points ()

!     Count and read the surface points listed one per line with optional label.
!     --------------------------------------------------------------------------

!     Local constants:

      character (2), parameter :: delimiters = ' ,'

!     Local variables:

      integer       :: first, i, ib, last, mark, nchar
      real          :: xy(2)
      character (7) :: format_string

!     Execution:

      open (luntarget, file=header_target%filename, status='old', iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') &
            ' Unable to open target file ', trim (header_target%filename)
         go to 99
      end if

      nb_target = 0
      do while (ios == 0)
         nb_target = nb_target + 1
         read (luntarget, '(a)', iostat=ios)
      end do
      nb_target = nb_target - 1

      rewind (luntarget)

      header_target%ndim        = ndim
      header_target%numq        = 0
      header_target%nblocks     = nb_target
      header_target%ndatasetaux = 0

      allocate (target(nb_target))

      len_buffer = len (buffer);  format_string = '(fnn.0)'

      do ib = 1, nb_target

         target(ib)%ni = 1;  target(ib)%nj = 1;  target(ib)%nk = 1

         call Tec_block_allocate (target(ib), ndim, header_target%numq, ios)

         read (luntarget, '(a)') buffer

!        Parse the line in order to handle descriptors with embedded blanks:

         first = 1;  last = len_buffer

         do i = 1, 2

            call scan2 (buffer, delimiters, first, last, mark)

            nchar = mark - first + 1
            if (nchar < 10) then
               write (format_string(2:4), '(a2, i1)') ' f', nchar
            else
               write (format_string(2:4), '(a1, i2)')  'f', nchar
            end if

            read (buffer(first:mark), format_string, iostat=ios) xy(i)

            if (ios /= 0) then
               write (luncrt, '(/, a, i5)') &
                  ' Trouble reading target (x,y).  Line #:', ib
               go to 99
            end if

!           In case the points aren't labeled, enter x as the zone title:

            if (i == 1) target(ib)%zone_title = buffer(first:mark)

            first = mark + 2
         end do

         target(ib)%x = xy(1);  target(ib)%y = xy(2)

         if (first <= last) target(ib)%zone_title = buffer(first:last)

         target(ib)%nzoneaux     = 0
         target(ib)%solutiontime = -999.  ! Undefined

      end do

 99   continue

      close (luntarget)

      end subroutine read_target_points

   end program surface_interp_2d
