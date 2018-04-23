!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program tri_to_quad

!  Description:
!
!     Interpolate an unstructured surface flow solution to a structured mesh.
!  Initially, the input is a single-zone Tecplot surface triangulation with a
!  variable number of quantities defined at each vertex, and the target mesh is
!  a multiblock PLOT3D-type grid.  Output is in either Tecplot or PLOT3D form.
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
!     TRI_TO_QUAD control file
!     xxxx.xxx        ! Surface triangulation and flow solution (Tecplot format)
!     1               ! 1 => vertex-centered/Tecplot; 2 => cell-centered/Tecplot
!     xxxx.xxx        ! Target grid (PLOT3D multiblock)
!     T               ! T|F = formatted|unformatted
!     xxxx.xxx        ! Output file 1 (grid file if PLOT3D format)
!     xxxx.xxx        ! Output file 2 (function file if PLOT3D, else 'none')
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
!     >  Target blocks are processed in order (read and written as needed).
!     >  For each target point, the ADT search locates the nearest point on the
!        surface triangulation (never outside it).  The function values can then
!        be interpolated with the indicated trilinear coefficients.
!     >  Optionally, target points beyond the indicated X will be set to the
!        indicated F (probably zero, as for radiation on an aft body).
!
!  Note:
!
!     This program still uses the original Tecplot I/O package.  It has not
!     been updated to use the I/O package version prompted by Tecplot 360.
!
!  History:
!
!     03/04/05  DAS  Initial implementation for vertex-centered function data.
!     03/07/05   "   Changed name from SURFACE_FLOW_INTERP to TRI_TO_QUAD in
!                    anticipation of a QUAD_TO_TRI analog.
!     08/29/08   "   Added optional X and F to deal with radiative heating
!                    that goes to zero beyond some X for CEV.
!     03/15/10   "   Provided for DPLR/Overset-related triangulated input data
!                    (type 2, cell-centered Tecplot format).
!     08/08/13   "   All ADT variants are in one module now with generic
!                    build_adt and search_adt interfaces.
!     07/02/14   "   The tri_* modules come from triangulation_io.f90, not
!                    tecplot_io_module (comments were wrong).
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!                           Now with ERC Incorporated/NASA ARC.
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! Must be the tecplot_io_module version
   use tecplot_io_module     ! Tecplot file I/O (structured data)
   use tri_header_structure  ! Part of triangulation_io.f90
   use tri_zone_structure    ! Likewise
   use triangulation_io      ! Unstructured data file I/O
   use xyzq_io_module        ! PLOT3D file I/O
   use adt_utilities         ! All ADT variants

   implicit none

!  Constants:

   integer, parameter :: &
      lunin     = 1,     &   ! Input triangulated surface data
      luntarget = 2,     &   ! Input target grid (PLOT3D multiblock)
      lunout1   = 3,     &   ! Output Tecplot file or PLOT3D grid
      lunout2   = 4,     &   ! Output PLOT3D function file if specified
      lunctl    = 5,     &   ! Control file on standard input (command line)
      luncrt    = 6          ! For diagnostics

   real, parameter :: &
      zero = 0.

!  Variables:

   integer :: &
      fileform, i, i1, i2, i3, ib, ios, itri, iz, j, nblocks, nf, ni, ninside, &
      nj, nnode, nout, npts, ntri, nzone

   real :: &
      davg, dist, dmax, dsq, dtol, f_beyond_xub, p, q, r, x_upper_bound

   real, dimension (3) :: &
      interp_xyz, surface_xyz

   logical :: &
      formatted_target, formatted_out, PLOT3D_out, Tecplot_out

   character :: &
      filename * 80, filename2 * 80, title_out * 80

   real, allocatable, dimension (:) :: &
      area_total, tri_area

   real, allocatable, dimension (:,:) :: &
      fnode

!  Derived data types:

   type (grid_type), pointer, dimension (:) :: &
      xyzf_interp

   type (tri_header_type) :: &
      tri_header

   type (tri_type), pointer, dimension (:) :: &
      tri_xyzf

!  Execution:
!  !!!!!!!!!!

   read (lunctl, '(a)', iostat=ios) ! Skip title

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Usage:  tri_to_quad < xxx.xxx'
      go to 99
   end if

   tri_header%formatted = .true.  ! Unformatted tri. reads are not handled
   tri_header%nvertices = 3       ! Triangles, not tetrahedra

   read (lunctl, *) tri_header%filename
   read (lunctl, *) fileform   ! 1|2 <-> vertex-|cell-centered

   if (fileform <= 2) then
      tri_header%fileform = fileform
   else
      write (luncrt, '(/, a, i6)') &
         ' Input dataset form must be 1 or 2:', fileform
      go to 99
   end if

   read (lunctl, *) filename      ! Target coordinates
   title_out = trim (filename)

   read (lunctl, *) formatted_target
   if (formatted_target) then
      open (luntarget, file=filename, status='old', form='formatted',          &
            iostat=ios)
   else
      open (luntarget, file=filename, status='old', form='unformatted',        &
            iostat=ios)
   end if
   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Unable to open target grid file:', filename
      go to 99
   end if

   read (lunctl, *) filename
   open (lunout1, file=filename, status='unknown', iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Unable to open output file:', filename
      go to 99
   end if

   read (lunctl, *) filename2
   read (lunctl, *) formatted_out

   Tecplot_out = filename2 (1:4) == 'none'
   PLOT3D_out  = .not. Tecplot_out

   if (PLOT3D_out) then
      open (lunout2, file=filename2, status='unknown', iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Unable to open output function file:',     &
            filename2
         go to 99
      end if
   end if

   read (lunctl, *) dtol ! Distance tolerance for diagnostic purposes only

!  Optional control:

   x_upper_bound = 1.e+30;  f_beyond_xub = zero   ! This f will not be used

   read  (lunctl, *, iostat=ios) x_upper_bound, f_beyond_xub
   close (lunctl)

!  Read the input surface data (including work-space allocation):
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ios = 1  ! Verbose mode

   call tri_read (lunin, tri_header, tri_xyzf, ios)  ! Read entire dataset

   if (ios /= 0) then
      write (luncrt, '(/, a)') 'Trouble reading triangulated dataset: aborting.'
      go to 99
   end if

   nzone = tri_header%nzones
   nf    = tri_header%numf

!  The interpolation scheme requires vertex-centered data.
!  If necessary, area-average the input functions to the vertices, in-place.

   if (fileform == 2) then

      do iz = 1, nzone

         nnode = tri_xyzf(iz)%nnodes
         ntri  = tri_xyzf(iz)%nelements

         allocate (tri_area(ntri), area_total(nnode))

         call tri_areas (nnode, ntri, tri_xyzf(iz)%xyz, tri_xyzf(iz)%conn,     &
                         tri_area, area_total)

         allocate (fnode(nf,nnode))

         call tri_centers_to_vertices (nnode, ntri, nf, tri_area, area_total,  &
                                       tri_xyzf(iz)%conn, tri_xyzf(iz)%f, fnode)

         deallocate (tri_area, area_total)

         deallocate (tri_xyzf(iz)%f);  allocate (tri_xyzf(iz)%f(nf,nnode))

         tri_xyzf(iz)%f(:,:) = fnode(:,:);  deallocate (fnode)

      end do

   end if

!  Set up the search tree:

   nnode = tri_xyzf(1)%nnodes     ! One zone is assumed for now
   ntri  = tri_xyzf(1)%nelements

   call build_adt (nnode, ntri, tri_xyzf(1)%conn, tri_xyzf(1)%xyz)


!  Read the target grid header and allocate the array of data structures:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call xyz_header_io (1, luntarget, formatted_target, nblocks, xyzf_interp,   &
                       ios)
   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble reading the target grid file header.'
      go to 99
   end if


!  Set up the output file(s):
!  !!!!!!!!!!!!!!!!!!!!!!!!!!

   do ib = 1, nblocks
      xyzf_interp(ib)%mi = xyzf_interp(ib)%ni
      xyzf_interp(ib)%mj = xyzf_interp(ib)%nj
      xyzf_interp(ib)%mk = 1
   end do

   if (Tecplot_out) then

      call Tec_header_write (lunout1, filename, formatted_out, 2, title_out,   &
                             nf, tri_header%varname, ios)
      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble writing Tecplot file header.'
         go to 99
      end if

   else ! PLOT3D-type output

      call xyz_header_io (2, lunout1, formatted_out, nblocks, xyzf_interp, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble writing the grid file header.'
         go to 99
      end if

      call q_header_io (2, lunout2, formatted_out, nblocks, nf, xyzf_interp,   &
                        ios)
      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble writing the function file header.'
         go to 99
      end if

   end if


!  Process one target block at a time:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   write (luncrt, '(a)')

   do ib = 1, nblocks

      ni = xyzf_interp(ib)%ni
      nj = xyzf_interp(ib)%nj

!     Read an input target block:

      call xyz_allocate (xyzf_interp(ib), ios)

      if (ios /= 0) then
         write (luncrt, '(/, a, i5)') &
            ' Allocation trouble for target grid.  Block #:', ib
         go to 99
      end if

      npts = ni * nj

      call xyz_block_io (1, luntarget, formatted_target, npts,                 &
                         xyzf_interp(ib)%x, xyzf_interp(ib)%y,                 &
                         xyzf_interp(ib)%z, ios)
      if (ios /= 0) then
         write (luncrt, '(/, a, i5)') &
            ' Trouble reading target grid.  Block #:', ib
         go to 99
      end if

!     Allocate the corresponding output block of function values:

      call q_allocate (xyzf_interp(ib), nf, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a, i5)') &
            ' Allocation trouble for function file output.  Block #:', ib
         go to 99
      end if

!     Process all the points of this target block:
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ninside = 0
      nout    = 0
      dmax    = zero
      davg    = zero

      do j = 1, nj

         do i = 1, ni

            surface_xyz(1) = xyzf_interp(ib)%x(i,j,1)
            surface_xyz(2) = xyzf_interp(ib)%y(i,j,1)
            surface_xyz(3) = xyzf_interp(ib)%z(i,j,1)

!           The following kludge was needed for CEV radiative heating on
!           aft-body patches of a full CFD surface grid:

            if (surface_xyz(1) > x_upper_bound) then

               xyzf_interp(ib)%q(:,i,j,1) = f_beyond_xub
               cycle

            end if

            call search_adt (surface_xyz, itri, p, q, r, dsq, .true., nnode,   &
                            ntri, tri_xyzf(1)%conn, tri_xyzf(1)%xyz, interp_xyz)
!           !!!!!!!!!!!!!!!

            dist = sqrt (dsq)
            dmax = max (dmax, dist)
            davg = davg + dist

            if (dist < dtol) then     ! The best triangle found was within the
               ninside = ninside + 1  ! distance tolerance 
            else
               nout = nout + 1
            end if

            i1 = tri_xyzf(1)%conn(1,itri)
            i2 = tri_xyzf(1)%conn(2,itri)
            i3 = tri_xyzf(1)%conn(3,itri)

            xyzf_interp(ib)%q(:,i,j,1) = p * tri_xyzf(1)%f(:,i1) +             &
                                         q * tri_xyzf(1)%f(:,i2) +             &
                                         r * tri_xyzf(1)%f(:,i3)
         end do ! Next i

      end do ! Next j

      davg = davg / real (npts)

      write (luncrt, '(a, i4, a, i7, a, 2i7, a, 1p, 2e15.5)') &
        ' Target block:', ib, '    # points:', npts, '    # in, out:',         &
        ninside, nout, '    average & max. distance error:', davg, dmax

!     Save this output block:

      if (Tecplot_out) then

         call Tec_block_write (lunout1, formatted_out, 2, title_out, ni, nj, 1,&
                               nf, xyzf_interp(ib)%x, xyzf_interp(ib)%y,       &
                               xyzf_interp(ib)%z, xyzf_interp(ib)%q, ios)
         if (ios /= 0) then
            write (luncrt, '(/, a, i4)') ' Trouble writing Tecplot block #', ib
            go to 99
         end if

      else ! PLOT3D-type output

         call xyz_block_io (2, lunout1, formatted_out, ni*nj,                  &
                            xyzf_interp(ib)%x, xyzf_interp(ib)%y,              &
                            xyzf_interp(ib)%z, ios)
         if (ios /= 0) then
            write (luncrt, '(/, a, i4)') ' Trouble writing grid block #', ib
            go to 99
         end if

         call q_block_io (2, lunout2, formatted_out, nf, ni, nj, 1,            &
                          xyzf_interp(ib)%q, ios)
         if (ios /= 0) then
            write (luncrt, '(/, a, i4)') ' Trouble writing function block #', ib
            go to 99
         end if

      end if

      deallocate (xyzf_interp(ib)%x, xyzf_interp(ib)%y, xyzf_interp(ib)%z,     &
                  xyzf_interp(ib)%q)

   end do ! Next target block

99 continue

   end program tri_to_quad
