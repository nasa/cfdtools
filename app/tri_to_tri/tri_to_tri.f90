!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program tri_to_tri

!  Description:
!
!     Interpolate one unstructured surface flow solution to another in 3-space.
!  This is an adaptation of the earlier TRI_TO_QUAD, replacing the structured
!  target surface grid with a surface triangulation.  The original single-zone
!  input dataset assumption has now been overcome via an extension of the tri-
!  angulation_io module.  The target (and hence output) triangulation may also
!  contain one or more zones as originally.
!
!  Input Tecplot data format (original, vertex-centered):
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
!     :::::::::::::::::
!
!  Further zones are now an option.
!
!  Alternative cell-centered Tecplot input format (DPLR overset grid tools):
!
!     TITLE     = ""
!     VARIABLES = "x"
!     "y"
!     "z"
!     "p"
!     "Chm"
!     "h"
!     "qw"
!     "Re_c"
!     ZONE T="ZONE 001"
!      STRANDID=0, SOLUTIONTIME=0
!      Nodes=9367, Elements=18438, ZONETYPE=FETriangle
!      DATAPACKING=BLOCK
!      VARLOCATION=([4-8]=CELLCENTERED)
!      FACENEIGHBORCONNECTIONS=55014
!      FACENEIGHBORMODE=LOCALONETOONE
!      FEFACENEIGHBORSCOMPLETE=YES
!      DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )
!      3.756540120E-01 3.601163924E-01 3.451932967E-01 3.309260905E-01  ...
!     ::::::::::::::::
!
!  Control file format (standard input):
!
!     TRI_TO_TRI control file
!     xxxx.xxx        ! Surface triangulation and flow solution (Tecplot format)
!     1               ! 1 => vertex-centered/Tecplot; 2 => cell-centered/Tecplot
!     xxxx.xxx        ! Target grid (Tecplot multizone surface triangulation)
!     T               ! T|F = formatted|unformatted
!     xxxx.xxx        ! Output file (Tecplot multizone surface triangulation)
!     T               ! T|F = formatted|unformatted
!     0.0001          ! Distance tolerance for target pt. inside solution grid
!    [7.890     0.    ! Optional X & F; if x > X, set finterp(:) = F]
!
!  Method:
!
!     >  The surface triangulation is read as a single zone and converted to an
!        ADT (Alternating Digital Tree) for search purposes.  If the function
!        values are cell-centered, they are interpolated to the cell vertices
!        first.  No unique best interpolation method exists.  The area-weighted
!        averaging used here suffices for typical triangulations.
!     >  Target zones are processed in order (read and written as needed).
!     >  For each target point, the ADT search locates the nearest point on the
!        surface triangulation (never outside it).  The function values can then
!        be interpolated with the indicated trilinear coefficients.
!     >  Optionally, target points beyond the indicated X will be set to the
!        indicated F (probably zero, as for radiation on an aft body).
!
!  History:
!
!     03/07/05  DAS  Initial TRI_TO_QUAD for vertex-centered function data.
!     08/29/08   "   Added optional X and F to deal with radiative heating
!                    that goes to zero beyond some X for CEV.
!     03/15/10   "   Provided for DPLR/Overset-related triangulated input data
!                    (type 2, cell-centered Tecplot format).
!     08/08/13   "   All ADT variants are in one module now with generic
!                    build_adt and search_adt interfaces.
!     07/03/14   "   TRI_TO_TRI adapted from TRI_TO_QUAD.
!     04/21/17   "   The tri_area(:) array name clashed with a callable utility
!                    added subsequently to triangulation_io.f90.  It has been
!                    renamed as triang_area(:).
!     07/20/18   "   Hemispherical integration of radiance data involves multi-
!                    zone triangulations.  Interpolation of such surface data
!                    to a different triangulation is now an option following
!                    an extension within triangulation_io that concatenates
!                    all zones as a single zone within new header fields, for
!                    rapid searching via the ADT scheme, which builds all cells
!                    of all zones into its search tree rather than treating one
!                    zone at a time.
!     05/03/22   "   One glitch following triangulation_io extensions.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA
!                  Now with AMA, Inc. at NASA ARC.
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use tri_header_structure  ! Part of triangulation_io.f90
   use tri_zone_structure    ! Likewise
   use triangulation_io      ! Unstructured data file I/O
   use adt_utilities         ! All ADT variants

   implicit none

!  Constants:

   integer, parameter :: &
      lunin     = 1,     &   ! Input triangulated surface data
      luntarget = 2,     &   ! Input target grid (multizone triangulation)
      lunout    = 3,     &   ! Output Tecplot multizone triangulation
      lunctl    = 5,     &   ! Control file on standard input (command line)
      luncrt    = 6          ! For diagnostics

   real, parameter :: &
      zero = 0.

!  Variables:

   integer :: &
      fileform, i, i1, i2, i3, ios, itri, iz, nelements, nf, nftarg, &
      ninside, nnode, nnodes, nout, ntri, nzone
   real :: &
      davg, dist, dmax, dsq, dtol, f_beyond_xub, p, q, r, x_upper_bound
   real, dimension (3) :: &
      interp_xyz
   logical :: &
      formatted_target
   character :: &
      filename * 80
   real, allocatable, dimension (:) :: &
      area_total, triang_area
   real, allocatable, dimension (:,:) :: &
      fnode

!  Derived data types:

   type (tri_header_type) :: &
      tri_header_1, tri_header_2, tri_header_3
   type (tri_type), pointer, dimension (:) :: &
      tri_xyzf_1, tri_xyzf_2

!  Execution:
!  !!!!!!!!!!

   read (lunctl, '(a)', iostat=ios) ! Skip title

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Usage:  tri_to_tri < xxx.xxx'
      go to 99
   end if

   tri_header_1%formatted = .true.  ! Unformatted has not been handled
   tri_header_1%nvertices = 3       ! Triangles, not tetrahedra

   read (lunctl, *) tri_header_1%filename  ! Input triangulated dataset
   read (lunctl, *, iostat=ios) fileform   ! 1|2 <-> vertex-|cell-centered

   if (fileform <= 2) then
      tri_header_1%fileform = fileform
   else
      write (luncrt, '(/, a, i6)') &
         ' Input dataset form must be 1 or 2:', fileform
      go to 99
   end if

   read (lunctl, *)         filename   ! Target triangulation
   tri_header_2%filename  = filename

   read (lunctl, *)         formatted_target
   tri_header_2%formatted = formatted_target
   tri_header_2%nvertices = 3

   read (lunctl, *)         filename   ! Output triangulation
   tri_header_3%filename  = filename
   tri_header_3%fileform  = 1          ! Vertex-centered

   read (lunctl, *) tri_header_3%formatted
   read (lunctl, *) dtol               ! Tolerance for diagnostic purposes only

!  Optional control:

   x_upper_bound = 1.e+30;  f_beyond_xub = zero   ! This f will not be used

   read  (lunctl, *, iostat=ios) x_upper_bound, f_beyond_xub
   close (lunctl)

!  Open and read the input surface dataset (including work-space allocation):
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ios = 1  ! Verbose mode
   tri_header_1%combine_zones         = .true.  ! Option prompted by TRI_TO_TRI
   tri_header_1%centroids_to_vertices = .true.  ! If input is cell-centered

   call tri_read (lunin, tri_header_1, tri_xyzf_1, ios)  ! Read all zones as one

   if (ios /= 0) then
      write (luncrt, '(/, a)') &
         'Trouble reading input triangulated dataset: aborting.'
      go to 99
   end if

   nf = tri_header_1%numf

!  The interpolation scheme requires vertex-centered data.
!  If necessary, area-average the input functions to the vertices, in-place.
!  This is now done within triangulation_io.

!! if (fileform == 2) then
!!
!!    do iz = 1, tri_header_1%nzones
!!
!!       nnode = tri_xyzf_1(iz)%nnodes
!!       ntri  = tri_xyzf_1(iz)%nelements
!!
!!       allocate (triang_area(ntri), area_total(nnode))
!!
!!       call tri_areas (nnode, ntri, tri_xyzf_1(iz)%xyz, tri_xyzf_1(iz)%conn, &
!!                       triang_area, area_total)
!!
!!       allocate (fnode(nf,nnode))
!!
!!       call tri_centers_to_vertices (nnode, ntri, nf, triang_area, &
!!                                     area_total, tri_xyzf_1(iz)%conn, &
!!                                     tri_xyzf_1(iz)%f, fnode)
!!
!!       deallocate (triang_area, area_total)
!!
!!       deallocate (tri_xyzf_1(iz)%f);  allocate (tri_xyzf_1(iz)%f(nf,nnode))
!!
!!       tri_xyzf_1(iz)%f(:,:) = fnode(:,:);  deallocate (fnode)
!!
!!    end do
!!
!! end if

!  Set up the search tree:

   nnode = tri_header_1%nnodes     ! Multiple zones have been merged now
   ntri  = tri_header_1%nelements

   call build_adt (nnode, ntri, tri_header_1%conn, tri_header_1%xyz)


!  Open and read the target triangulation header:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ios = 1  ! Verbose mode

   call tri_header_read (luntarget, tri_header_2, tri_xyzf_2, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') &
         'Trouble reading the target triangulation header.'
      go to 99
   end if

   nftarg = tri_header_2%numf  ! >= 0
   nzone  = tri_header_2%nzones


!  Set up the output file:
!  !!!!!!!!!!!!!!!!!!!!!!!

   tri_header_3%nvertices   = 3
   tri_header_3%numf        = nf
   tri_header_3%nzones      = nzone
   tri_header_3%datapacking = 0  ! Point order
   tri_header_3%title       = tri_header_2%title

   if (tri_header_2%title == " ") tri_header_3%title = tri_header_2%filename

   allocate (tri_header_3%varname(nf + 3))
   tri_header_3%varname(:)  = tri_header_1%varname(:)

   call tri_header_write (lunout, tri_header_3, tri_xyzf_2, ios)
   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble writing output file header.'
      go to 99
   end if


!  Process one target zone at a time:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   write (luncrt, '(a)')

   do iz = 1, nzone

      call tri_zone_allocate (tri_header_2, tri_xyzf_2(iz), ios)
      if (ios /= 0) go to 99

      call tri_zone_read (luntarget, tri_header_2, tri_xyzf_2(iz), ios)
      if (ios /= 0) then
         write (luncrt, '(/, a, i5)') &
            ' Trouble reading target zone #', iz
         go to 99
      end if

!     Allocate the corresponding output block of function values:

      deallocate (tri_xyzf_2(iz)%f)
      nnodes    = tri_xyzf_2(iz)%nnodes
      nelements = tri_xyzf_2(iz)%nelements
      allocate   (tri_xyzf_2(iz)%f(nf,nnodes))
      tri_xyzf_2(iz)%element_type = 'TRIANGLE'

!     Process all the points of this target zone:
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ninside = 0
      nout    = 0
      dmax    = zero
      davg    = zero

      do i = 1, nnodes

!        The following kludge was needed in TRI_TO_QUAD for CEV radiative
!        heating on aft-body patches of a full CFD surface grid:

         if (tri_xyzf_2(iz)%xyz(1,i) > x_upper_bound) then
             tri_xyzf_2(iz)%f(:,i) = f_beyond_xub
             cycle
         end if

         call search_adt (tri_xyzf_2(iz)%xyz(:,i), itri, p, q, r, dsq, &
                          .true., nnode, ntri,  &
                          tri_header_1%conn, tri_header_1%xyz, interp_xyz)
!        !!!!!!!!!!!!!!!

         dist = sqrt (dsq)
         dmax = max (dmax, dist)
         davg = davg + dist

         if (dist < dtol) then     ! The best triangle found was within the
            ninside = ninside + 1  ! distance tolerance 
         else
            nout = nout + 1
         end if

         i1 = tri_header_1%conn(1,itri)
         i2 = tri_header_1%conn(2,itri)
         i3 = tri_header_1%conn(3,itri)

         tri_xyzf_2(iz)%f(:,i) = p * tri_header_1%f(:,i1) + &
                                 q * tri_header_1%f(:,i2) + &
                                 r * tri_header_1%f(:,i3)
      end do ! Next node

      davg = davg / real (nnodes)

      write (luncrt, '(a, i4, a, i7, a, 2i7, a, 1p, 2e15.5)') &
        ' Target zone:', iz, '    # points:', nnodes, '    # in, out:', &
        ninside, nout, '    average & max. distance error:', davg, dmax

!     Save this output zone:

      call tri_zone_write (lunout, tri_header_3, tri_xyzf_2(iz), ios)

      if (ios /= 0) then
         write (luncrt, '(/, a, i4)') ' Trouble writing zone #', iz
         go to 99
      end if

      deallocate (tri_xyzf_2(iz)%xyz, tri_xyzf_2(iz)%f, tri_xyzf_2(iz)%conn)

   end do ! Next target zone

99 continue

   end program tri_to_tri
