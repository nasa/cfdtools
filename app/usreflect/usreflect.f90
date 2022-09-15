!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program usreflect
!
!  Description:
!
!     This is the unstructured-grid analogue of REFLECT_BLOCKS.
!
!     Provision is made for adding file types (probably HDF5 at some point),
!     but initially, Tecplottable files that US3D can write are handled.
!
!     Like REFLECT_BLOCKS, USREFLECT reflects the zones of a multizone grid
!     in the plane specified.  Optional function data may also be reflected.
!     Either both halves or just the reflected half may be output.  If function
!     data are present, one or more quantities may have the sign changed as part
!     of the reflection.  Normally this would be just one velocity component.
!
!     The option to clean up a symmetry plane is not viable here, for lack of a
!     way of identifying the points intended to be on the symmetry plane.
!
!     A handful of prompts suffice (no control file).
!
!  Strategy:
!
!     If both halves are to be output, all zones can be read and written as is.
!     The reflected half is always written, and these zones can be derived from
!     the zones already read and stored.  So we don't need to carry around an
!     extra set of zones.
!
!  Procedures:
!
!     triangulation_io   I/O utilities for multizone unstructured datasets
!     rdlist             Utility for reading an indefinite list of integers
!
!  History:
!
!     08/30/05  D.A.Saunders  Initial implementation of REFLECT_BLOCKS.
!     02/23/06    "     "     Belated realization that reflecting changes the
!                             handedness.  Restore the original handedness.
!     04/20/22    "     "     Adaptation as USREFLECT.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use tri_header_structure        ! Data structures used
   use tri_zone_structure          ! by triangulation_io
   use triangulation_io

   implicit none

!  Constants:

   integer, parameter :: &
      lunkbd = 5,        &
      luncrt = 6,        &
      lunin  = 7,        &
      lunout = 8

!  Variables:

   integer :: &
      i, ifile_type, ios, ir, ixyz, iz, l, m, n, nzones, nchange_sign, nfsets, &
      numf, nvertices

   integer, allocatable, dimension (:) :: &
      i_change_sign

   real :: &
      tol

   logical :: &
      both_halves, cell_centered, cr, clean_up_symmetry_plane, eof, &
      formatted_in, formatted_out, trisurf

   character (1) :: &
      answer

!  Derived data types:

   type (tri_header_type) :: &
      header

   type (tri_type), pointer, dimension (:) :: &
      xyzf

!  Execution:

   call reads (luncrt, 'Input dataset?  ', lunkbd, header%filename, cr, eof)
   if (cr .or. eof) go to 99

   open (lunin, file=header%filename, status='old', iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(2a)') 'Cannot open ', trim (header%filename)
      go to 99
   end if

   close (lunin)  ! tri_read or vol_read opens it

   l = len_trim (header%filename)
   if (header%filename(l-2:l) == 'dat') then
      ifile_type = 1  ! Tecplot ASCII
      header%formatted = .true.
   else if (header%filename(l-1:l) == 'h5') then
      ifile_type = 2  ! HDF5, to be completed
   end if

20 write (luncrt, '(a)', advance='no') &
      ' What coordinate is to be reflected?  [1|2|3 => x|y|z]: '
   read  (lunkbd, *, iostat=ios) ixyz
   if (ios /= 0) go to 20
   if (ixyz < 1 .or. ixyz > 3) go to 20

   write (luncrt, '(a)', advance='no') &
      ' Save both halves [y]?  [n = just save reflected half]: '
   read  (lunkbd, *) answer
   both_halves = answer == 'y' .or. answer == 'Y'

!! The following isn't viable because we can't be sure what cells are intended
!! to be on the symmetry plane.

!! write (luncrt, '(a)', advance='no') &
!!    ' Clean up the symmetry plane? [y|n]: '
!! read  (lunkbd, *) answer
!! clean_up_symmetry_plane = answer == 'y' .or. answer == 'Y'

!! if (clean_up_symmetry_plane) then
!!    write (luncrt, '(a)', advance='no') &
!!       ' Tolerance for average absolute value of appropriate coordinate: '
!!    read  (lunkbd, *) tol
!! end if

!  Read the input file.  The both halves case requires storing all the zones,
!  so there's no point in processing one zone at a time.

   if (ifile_type == 1) then
      call get_element_type (lunin, header, ios)
      if (ios /= 0) go to 99

      nvertices = header%nvertices
      header%combine_zones         = .false.
      header%centroids_to_vertices = .false.
      trisurf = nvertices == 3

      if (trisurf) then
         call tri_read (lunin, header, xyzf, ios)
      else
         call vol_read (lunin, header, xyzf, ios)
      end if
      if (ios /= 0) go to 99

      cell_centered = header%fileform == 2

   else if (ifile_type == 2) then
!     TBD
   end if

   numf   = header%numf
   nzones = header%nzones

   write (luncrt, '(a, i6)') ' # vertices per cell:', nvertices, &
                             ' # functions found:  ', numf

!  Display the zone details to reassure the user:

   write (luncrt, '(a)') '  Zone   # nodes  # elements'
   write (luncrt, '(i6, 2i10)') &
      (iz, xyzf(iz)%nnodes, xyzf(iz)%nelements, iz = 1, nzones)

!  Will functions like v need sign changes?

   if (numf > 0) then
      allocate (i_change_sign(numf))

      nchange_sign = numf

      call rdlist (luncrt, 'Function number(s) needing sign change: ',      &
                   lunkbd, nchange_sign, i_change_sign)
   end if

!  Set up the output file.  Reuse the input file header:

   call reads (luncrt, 'Output dataset? ', lunkbd, header%filename, cr, eof)
   if (cr .or. eof) go to 99

!  Transcribe the input zones?  (Avoid closing the file.)

   if (both_halves) then

      if (ifile_type == 1) then

         if (trisurf) then
            call tri_header_write (lunout, header, xyzf, ios)
            if (ios /= 0) go to 99

            do iz = 1, nzones
               call tri_zone_write (lunout, header, xyzf(iz), ios)
               if (ios /= 0) go to 99
            end do
         else
            call vol_header_write (lunout, header, xyzf, ios)
            if (ios /= 0) go to 99

            do iz = 1, nzones
               call vol_zone_write (lunout, header, xyzf(iz), ios)
               if (ios /= 0) go to 99
            end do
         end if

      else if (ifile_type == 2) then
!        TBD
      end if

   else ! Just the reflected half is written; write the header:

      if (ifile_type == 1) then

         if (trisurf) then
            call tri_header_write (lunout, header, xyzf, ios)
         else
            call vol_header_write (lunout, header, xyzf, ios)
         end if
         if (ios /= 0) go to 99

      else if (ifile_type == 2) then
!        TBD
      end if

   end if

!  Reflect the zones one at a time and write them:

   do iz = 1, nzones
      xyzf(iz)%xyz(ixyz,:) = -xyzf(iz)%xyz(ixyz,:)

      call fix_handedness (nvertices, xyzf(iz))

      if (nchange_sign > 0) then
         do m = 1, nchange_sign
            n = i_change_sign(m)
            xyzf(iz)%f(n,:) = -xyzf(iz)%f(n,:)
         end do
      end if

      if (ifile_type == 1) then
         if (trisurf) then
            call tri_zone_write (lunout, header, xyzf(iz), ios)
         else
            call vol_zone_write (lunout, header, xyzf(iz), ios)
         end if
         if (ios /= 0) go to 99

      else if (ifile_type == 2) then
!        TBD
      end if

      deallocate (xyzf(iz)%xyz, xyzf(iz)%f, xyzf(iz)%conn, stat=ios)
      if (ios /= 0) then
         write (luncrt, '(a, 2i5)') &
            ' Trouble deallocating grid zone #', iz, ios
         go to 99
      end if

   end do ! Next reflected zone to write

   close (lunout)

99 continue

! *** stop ! Avoid system dependencies.

   end program usreflect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine fix_handedness (nvertices, zone)

!  Change the handedness of a grid zone by reversing the node indices in place.
!  Does this work for all cell types?
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use tri_zone_structure  ! External module

   implicit none

!  Arguments:

   integer,         intent (in)    :: nvertices  ! Per cell
   type (tri_type), intent (inout) :: zone       ! %conn(:,:) are reversed here

!  Local variables:

   integer :: i, i1, ir, n, nvby2
   real    :: ivert(nvertices)

!  Execution:

   nvby2 = nvertices/2

   do n = 1, zone%nelements
      ivert(:) = zone%conn(:,n)
      do i = 1, nvby2
         i1        = ivert(i)
         ir        = nvertices + 1 - i
         ivert(i)  = ivert(ir)
         ivert(ir) = i1
      end do
      zone%conn(:,n) = ivert(:)
   end do

   end subroutine fix_handedness
