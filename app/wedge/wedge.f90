!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program wedge
!
!  Description:
!
!     For a given grid and optional flow solution in Plot2D form (as from the
!     DPLR2D solver), generate a 3D wedge form of the file(s) as needed by some
!     other flow solver such as LAURA and US3D.  Input files are prompted for
!     (no control file), as is the total wedge angle.  (The wedge is symmetric
!     about the xz plane.)  Any functions are necessarily vertex-centered.
!     Output file(s) are in Plot3D form.
!
!     The initial singular line along Ox may optionally be tweaked as two
!     separate lines either side of Ox, 2*epsilon apart, via simple y shifts.
!
!     Earlier utility DECONSTRUCT by the present author can convert the wedge
!     file(s) into unstructured form with no repeated cell vertices.  At the
!     time of writing, Tecplot form with hexahedral cells is the only option,
!     but that may be supplemented with an HDF5 output option.
!
!  History:
!
!     11/23/2022  D.A.Saunders  Initial implementation following a conversation
!                               with Khalil Bensassi.
!     12/03/2022    "      "    Included an option to avoid the singular line
!                               along Ox, by simple y shifts.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure ! Or equivalent to define "grid_type"
   use xyq_io_module        ! Plot2D I/O utilities
   use xyzq_io_module       ! Plot3D I/O utilities

   implicit none

!  Constants:

   integer, parameter :: &
      lung2d = 1,        &
      lunf2d = 2,        &
      lung3d = 3,        &
      lunf3d = 4,        &
      lunkbd = 5,        &
      luncrt = 6         ! Self-descriptive logical units

!  Variables:

   integer :: &
      lunf, nblocks, nf, ns, ios

   real :: &
      angle, eps

   logical :: &
      cell_centered, formatted_in, formatted_out, function_data, velocities

   character (1) :: &
      answer

   character (132) :: &
      filenamef_in, filenamef_out, filenameg_in, filenameg_out

!  Derived data types:

   type (grid_type), pointer, dimension (:) :: &
      xyf, xyzf

!  Execution:

   call read_controls ()
   if (ios /= 0) go to 99

!  Read the input file(s):

   lunf = -lunf2d
   if (function_data) lunf = lunf2d

   call xyq_read (lung2d, lunf, formatted_in, nblocks, nf, cell_centered, &
                  xyf, ios)
   if (ios /= 0) go to 99

!  Construct the wedge:

   call construct_wedge ()
   if (ios /= 0) go to 99

!  Save the results:

   call xyz_write (lung3d, formatted_out, nblocks, xyzf, ios)
   if (ios /= 0) go to 99

   if (function_data) then
      if (velocities) nf = nf + 1
      call q_write (lunf3d, formatted_out, nblocks, nf, xyzf, ios)
      if (ios /= 0) go to 99
   end if

!  Let Fortran do remaining deallocations.

99 continue

!  Internal procedures for program wedge:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_controls ()
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i1
      character (11), parameter :: &
         format = 'unformatted'

      write (luncrt, '(/, a)', advance='no') 'Input grid file (Plot2D):      '
      read  (lunkbd, '(a)') filenameg_in

      write (luncrt, '(a)', advance='no') 'Input flow file (Plot2D|none): '
      read  (lunkbd, '(a)') filenamef_in
      function_data = filenamef_in /= 'none'

      write (luncrt, '(a)', advance='no') &
         'Formatted|unformatted input files? [f|u]: '
      read  (lunkbd, '(a)') answer
      formatted_in = answer == 'f' .or. answer == 'F'
      i1 = 1;  if (formatted_in) i1 = 3

      open  (lung2d, file=filenameg_in, form=format(i1:), status='old', &
             iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(2a)') '*** Unable to open input grid file ', &
            trim (filenameg_in)
         go to 99
      end if

      if (function_data) then
         open (lunf2d, file=filenamef_in, form=format(i1:), status='old', &
               iostat=ios)
         if (ios /= 0) then
            write (luncrt, '(2a)') '*** Unable to open input function file ', &
               trim (filenamef_in)
            go to 99
         end if
      endif

      write (luncrt, '(a)', advance='no') 'Output grid file (Plot3D): '
      read  (lunkbd, '(a)') filenameg_out

      if (function_data) then
         write (luncrt, '(a)', advance='no') 'Output flow file (Plot3D): '
         read  (lunkbd, '(a)') filenamef_out
      end if

      write (luncrt, '(a)', advance='no') &
         'Formatted|unformatted output? [f|u]: '
      read  (lunkbd, '(a)') answer
      formatted_out = answer == 'f' .or. answer == 'F'
      i1 = 1;  if (formatted_out) i1 = 3

      open (lung3d, file=filenameg_out, form=format(i1:), status='unknown')
      if (function_data) then
         open (lunf3d, file=filenamef_out, form=format(i1:), status='unknown')
      end if

      write (luncrt, '(a)', advance='no') '(Total) wedge angle: '
      read  (lunkbd, *) angle
      angle = angle/2.  ! Either side of y = 0.

      write (luncrt, '(a)', advance='no') &
         'Singular line epsilon shift? [0.|y --> y +/- eps] '
      read  (lunkbd, *) eps

      if (function_data) then
         write (luncrt, '(a)', advance='no') &
     ' # species preceding velocities [ns > 0, or 0 if (u,v) are not present]: '
         read  (lunkbd, *) ns
         velocities = ns > 0
      end if

      ios = 0

 99   return

      end subroutine read_controls

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine construct_wedge ()
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, parameter :: zero = 0.

      integer :: i, ib, j, k, ni, nj, nk, npts
      real    :: px, py, pz, qx, qy, qz
      real, allocatable, dimension (:,:) :: vy, vz

      allocate (xyzf(nblocks))

!     The only rotations are 2D rotations about the origin.

      do ib = 1, nblocks
         ni = xyf(ib)%ni
         nj = xyf(ib)%nj
         npts = ni*nj

         xyzf(ib)%ni = ni
         xyzf(ib)%nj = 2
         xyzf(ib)%nk = nj

         if (function_data) then
            xyzf(ib)%mi = ni
            xyzf(ib)%mj = 2
            xyzf(ib)%mk = nj
         end if

         call xyz_allocate (xyzf(ib), ios)
         if (ios /= 0) go to 99

         allocate (vy(ni,nj), vz(ni,nj))

         do j = 1, 2
            xyzf(ib)%x(:,j,:) = xyf(ib)%x(:,:,1)
            vy(:,:)           = zero
            vz(:,:)           = xyf(ib)%y(:,:,1)

            call rotate2D (npts, vy, vz, angle, zero, zero)

            angle = -angle
            xyzf(ib)%y(:,j,:) = vy(:,:)
            xyzf(ib)%z(:,j,:) = vz(:,:)
         end do

         deallocate (vy, vz)

!        Fudge the singular edge?  If so, shift each side of the wedge sideways:

         if (eps /= zero) then
            eps = -eps
            do j = 1, 2
               xyzf(ib)%y(:,j,:) = xyzf(ib)%y(:,j,:) + eps
               eps = -eps
            end do
            eps = -eps
         end if

      end do  ! Next block

      if (function_data) then

         do ib = 1, nblocks

            call q_allocate (xyzf(ib), nf+1, ios)
            if (ios /= 0) go to 99

            ni = xyf(ib)%mi
            nj = xyf(ib)%mj
            nk = nj
            npts = ni*nj

            if (velocities) then

               allocate (vy(ni,nj), vz(ni,nj))

               do j = 1, 2
                  do k = 1, nk
                     xyzf(ib)%q(   1:ns+1,:,j,k) = xyf(ib)%q( 1:ns+1,:,k,1)
                     xyzf(ib)%q(ns+4:nf+1,:,j,k) = xyf(ib)%q(ns+3:nf,:,k,1)
                  end do

                  vy(:,:) = zero
                  vz(:,:) = xyf(ib)%q(ns+2,:,:,1)
                  angle   = -angle

                  call rotate2d (npts, vy, vz, angle, zero, zero)

                  xyzf(ib)%q(ns+2,:,j,:) = vy(:,:)
                  xyzf(ib)%q(ns+3,:,j,:) = vz(:,:)
               end do

               deallocate (vy, vz)

            else  ! No velocity components

               do k = 1, nk
                  do j = 1, 2
                     xyzf(ib)%q(:,:,j,k) = xyf(ib)%q(:,:,k,1)
                  end do
               end do

            end if

         end do

      end if

 99   return

      end subroutine construct_wedge

   end program wedge
