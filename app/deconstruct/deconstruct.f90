!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program deconstruct
!
!     Description:
!
!        Convert a multiblock structured surface or volume grid to a surface
!     triangulation or tetrahedral volume mesh.  Quadrilaterals become pairs
!     of triangles; hex cells become five tetrahedra each.
!
!     Outputs are in Tecplot format, as a single zone.
!
!     Procedures:
!
!        XYZQ_IO package  I/O utilities for PLOT3D grid and function files
!
!     History:
!
!        08/16/04  D. Saunders  Initial adaptation of EXTRACT_BLOCKS, to
!                               provide data for testing the tetrahedral mesh
!                               variant of the ADT search package.
!
!     Author:  David Saunders, ELORET/NASA Ames, Moffett Field, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Modules:

      use grid_block_structure
      use xyzq_io_module

      implicit none

!     Constants:

      integer, parameter :: &
         lunkbd  = 5,       &
         luncrt  = 6,       &
         lun_in  = 1,       &
         lun_out = 2

!     Variables:

      integer :: &
         i, j, k, i1, i2, i3, i4, ib, ie, ios, &
         nblocks, nelements, ni, nij, nj, nk, npoints, npts, offset

      integer, allocatable, dimension (:,:) :: &
         conn

      logical :: &
         cell_centered, formatted_in, formatted_out

!     Composite data types:

      type (grid_type), pointer, dimension (:) :: &
         xyz

!     Execution:

!     Prompt for and open the input grid:

      call file_prompt (lun_in, -luncrt, 'input grid', 'old', .true., &
                        formatted_in, ios)
      if (ios /= 0) go to 999

!     Read the input grid header:

      call xyz_header_io (1, lun_in, formatted_in, nblocks, xyz, ios)

      if (ios /= 0) go to 999

!     Display the block dimensions to reassure the user:

      write (luncrt, '(a)')
      write (luncrt, '(4(i6, 2X, 3i5))') &
         (ib, xyz(ib)%ni, xyz(ib)%nj, xyz(ib)%nk, ib = 1, nblocks)

!     Prompt for and open the output file:

      formatted_out = .true.  ! Avoid writing Tecplot binaries for now

      call file_prompt (lun_out, -luncrt, 'output mesh', 'unknown', .false., &
                        formatted_out, ios)
      if (ios /= 0) go to 999

      nk = xyz(1)%nk
      npoints = 0
      nelements = 0

      if (nk == 1) then ! Triangulate the surface grid

         do ib = 1, nblocks
            npoints   = npoints   +  xyz(ib)%ni * xyz(ib)%nj 
            nelements = nelements + (xyz(ib)%ni - 1) * (xyz(ib)%nj - 1) 
         end do
         nelements = nelements * 2

         write (lun_out, '(a)')                &
            'TITLE = "Surface triangulation"', &
            'VARIABLES = "x" "y" "z"',         &
            'ZONE T = ""'
         write (lun_out, '(a, i8, a, i8, a)')  &
            'N =', npoints, ', E =', nelements, ', ZONETYPE=FETRIANGLE' 
         write (lun_out, '(a)')                &
            'DATAPACKING=POINT',               &
            'DT=(SINGLE SINGLE SINGLE)'

         do ib = 1, nblocks

            call xyz_allocate (xyz(ib), ios)
            if (ios /= 0) go to 999

            npts = xyz(ib)%ni * xyz(ib)%nj

            call xyz_block_io (1, lun_in, formatted_in, npts,          &
                               xyz(ib)%x, xyz(ib)%y, xyz(ib)%z, ios)
            if (ios /= 0) go to 999

            write (lun_out, '(1p, 3e15.7)')                            &
               ((xyz(ib)%x(i,j,1), xyz(ib)%y(i,j,1), xyz(ib)%z(i,j,1), &
               i = 1, xyz(ib)%ni), j = 1, xyz(ib)%nj)

            deallocate (xyz(ib)%x, xyz(ib)%y, xyz(ib)%z, stat=ios)

            if (ios /= 0) then
               write (luncrt, '(/, a, 2i5)') &
                  ' Trouble deallocating input grid block #', ib, ios
               go to 999
            end if

         end do

         allocate (conn(3,nelements))

         offset = 0 ! For current block
         ie = 0

         do ib = 1, nblocks
            ni = xyz(ib)%ni
            nj = xyz(ib)%nj
            do j = 1, nj - 1
               do i = 1, ni - 1
                  ie = ie + 1
                  i1 = i + (j - 1) * ni + offset ! (i,j)
                  i2 = i1 + 1                    ! (i+1,j)
                  i3 = i1 + ni                   ! (i,j+1)
                  conn(1,ie) = i1
                  conn(2,ie) = i2 
                  conn(3,ie) = i3 
                  ie = ie + 1
                  conn(1,ie) = i2
                  conn(2,ie) = i2 + ni   ! (i+1,j+1)
                  conn(3,ie) = i3
               end do
            end do
            offset = offset + ni * nj
         end do

         if (ie /= nelements) then
            write (luncrt, '(/, a, 2i9)') & 
               ' ???:  ie, nelements =', ie, nelements
         else
            write (lun_out, '(1p, i8, 2i9)') conn
         end if

      else ! Generate tetrahedra from the hex cells

         do ib = 1, nblocks
            npoints   = npoints   +  xyz(ib)%ni * xyz(ib)%nj * xyz(ib)%nk
            nelements = nelements + (xyz(ib)%ni - 1) * (xyz(ib)%nj - 1) * &
                                    (xyz(ib)%nk - 1)
         end do
         nelements = nelements * 5

         write (lun_out, '(a)')                &
            'TITLE = "Tetrahedral mesh"',      &
            'VARIABLES = "x" "y" "z"',         &
            'ZONE T = ""'
         write (lun_out, '(a, i8, a, i8, a)')  &
            'N =', npoints, ', E =', nelements, ', ZONETYPE=FETETRAHEDRON'
         write (lun_out, '(a)')                &
            'DATAPACKING=POINT',               &
            'DT=(SINGLE SINGLE SINGLE SINGLE)'

         do ib = 1, nblocks

            call xyz_allocate (xyz(ib), ios)
            if (ios /= 0) go to 999

            npts = xyz(ib)%ni * xyz(ib)%nj * xyz(ib)%nk

            call xyz_block_io (1, lun_in, formatted_in, npts,           &
                               xyz(ib)%x, xyz(ib)%y, xyz(ib)%z, ios)
            if (ios /= 0) go to 999

            write (lun_out, '(1p, 3e15.7)')                             &
               (((xyz(ib)%x(i,j,k), xyz(ib)%y(i,j,k), xyz(ib)%z(i,j,k), &
               i = 1, xyz(ib)%ni), j = 1, xyz(ib)%nj), k = 1, xyz(ib)%nk)

            deallocate (xyz(ib)%x, xyz(ib)%y, xyz(ib)%z, stat=ios)

            if (ios /= 0) then
               write (luncrt, '(/, a, 2i5)') &
                  ' Trouble deallocating input grid block #', ib, ios
               go to 999
            end if

         end do

         allocate (conn(4,nelements))

         offset = 0 ! For current block
         ie = 0

         do ib = 1, nblocks
            ni  = xyz(ib)%ni
            nj  = xyz(ib)%nj
            nk  = xyz(ib)%nk
            nij = ni * nj
            do k = 1, nk - 1
               do j = 1, nj - 1
                  do i = 1, ni - 1
                     ie = ie + 1
                     i1 = i + (j - 1) * ni + (k - 1) * nij + offset  ! (i,j,k)
                     i2 = i1 + 1                                     ! (i+1,j,k)
                     i3 = i1 + ni                                    ! (i,j+1,k)
                     i4 = i1 + nij                                   ! (i,j,k+1)
                     conn(1,ie) = i1           ! (i,j,k)
                     conn(2,ie) = i2           ! (i+1,j,k)
                     conn(3,ie) = i3           ! (i,j+1,k)
                     conn(4,ie) = i4           ! (i,j,k+1)
                     ie = ie + 1
                     conn(1,ie) = i3 + 1       ! (i+1,j+1,k)
                     conn(2,ie) = i3           ! (i,j+1,k)
                     conn(3,ie) = i2           ! (i+1,j,k)
                     conn(4,ie) = i4 + 1 + ni  ! (i+1,j+1,k+1)
                     ie = ie + 1
                     conn(1,ie) = i2 + nij     ! (i+1,j,k+1)
                     conn(2,ie) = i2           ! (i+1,j,k)
                     conn(3,ie) = i4           ! (i,j,k+1)
                     conn(4,ie) = i4 + 1 + ni  ! (i+1,j+1,k+1)
                     ie = ie + 1
                     conn(1,ie) = i3 + nij     ! (i,j+1,k+1)
                     conn(2,ie) = i4           ! (i,j,k+1)
                     conn(3,ie) = i3           ! (i,j+1,k)
                     conn(4,ie) = i4 + 1 + ni  ! (i+1,j+1,k+1)
                     ie = ie + 1
                     conn(1,ie) = i2           ! (i+1,j,k)
                     conn(2,ie) = i3           ! (i,j+1,k)
                     conn(3,ie) = i4           ! (i,j,k+1)
                     conn(4,ie) = i4 + 1 + ni  ! (i+1,j+1,k+1)
                  end do
               end do
            end do
            offset = offset + nij * nk
         end do

         if (ie /= nelements) then
            write (luncrt, '(/, a, 2i9)') &
               ' ???:  ie, nelements =', ie, nelements
         else
            write (lun_out, '(1p, i8, 3i9)') conn
         end if

      end if

  999 continue

! *** stop ! Avoid system dependencies.

      end program deconstruct
