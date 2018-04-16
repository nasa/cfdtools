!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program test_block_utilities
!
!  Driver for the options in grid_block_utilities.f90.
!
!  07/26/2010  D.A.Saunders  Testing the cell volume option is the main concern.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! From grid_block_utilities or equivalent
   use grid_block_utilities  ! Structured grid block utilities
   use xyzq_io_module        ! I/O package for PLOT3D files

   implicit none

!  Constants:

   integer, parameter :: lunin = 1, lunout = 2, lunf = 3, lunkbd = 5, luncrt = 6

!  Variables:

   integer   :: ib, icase, ios, lun, m, n, nb, nf
   real      :: bbox(6), dx, dy, dz, scale
   logical   :: cell_centered, formatted
   character :: filein*80, fileout*80
   type (grid_type), pointer :: grid(:)
   type (grid_type), pointer :: cg_vol(:)

!  Execution:

   write (luncrt, '(/, (a))') 'Options:',                                   &
      '   1. block_data_range (block)',                                     &
      '   2. centroids_volumes (bin, bout)',                                &
      '   3. clone_structured_grid (nb, nf, grid1, grid2)',                 &
      '   4. get_bounding_boxes (lun, nblocks, grid, formatted, bbox)',     &
      '   5. reflect_block_xyz (block, n)',                                 &
      '   6. reverse_block_i (block, nf)',                                  &
      '   7. reverse_block_j (block, nf)',                                  &
      '   8. reverse_block_k (block, nf)',                                  &
      '   9. scale_shift_block (block, scale, dx, dy, dz, block_out)',      &
      '  10. swap_block_xyz (block, m, n)',                                 &
      '  11. update_block (block_in, nf, release_in, release_out, block_out)'

10 continue

   write (luncrt, '(a)', advance='no') 'Pick an option: '
   read  (lunkbd, *) icase

   formatted = .true.  ! Adequate here

   write (luncrt, '(a)', advance='no') 'Input grid:  '
   read  (lunkbd, *) filein
   open  (lunin, file=filein, form='formatted', status='old')
   call xyzq_read (lunin, -lunin, formatted, nb, nf, cell_centered, grid, ios)
   if (ios /= 0) go to 99

   write (luncrt, '(a)', advance='no') 'Output grid: '
   read  (lunkbd, *) fileout
   open  (lunout, file=fileout, form='formatted', status='new')

   if (icase == 2) then
      write (luncrt, '(a)', advance='no') 'Output function file: '
      read  (lunkbd, *) fileout
      open  (lunf, file=fileout, form='formatted', status='new')
   end if

   select case (icase)

   case (1)

      do ib = 1, nb

         call block_data_range (grid(ib))

         write (luncrt, '(i4, 1p, 6e15.6)') ib, &
            grid(ib)%xmin, grid(ib)%xmax, &
            grid(ib)%ymin, grid(ib)%ymax, &
            grid(ib)%zmin, grid(ib)%zmax
      end do

   case (2)

      allocate (cg_vol(nb))

      nf = 1

      do ib = 1, nb
         call centroids_volumes (grid(ib), cg_vol(ib))
      end do

      call xyz_write (lunout, formatted, nb, cg_vol, ios)
      call   q_write (lunf,   formatted, nb, nf, cg_vol, ios)

   case (3)

      call clone_structured_grid (nb, nf, grid, cg_vol)

   case (4)

      lun = -1
      call get_bounding_boxes (lun, nb, grid, formatted, bbox)

      write (luncrt, '(1p, 6e15.6)') bbox

   case (5)

      do ib = 1, nb
         do n = 1, 3
            call reflect_block_xyz (grid(ib), n)
         end do
      end do

   case (6)

      do ib = 1, nb
         call reverse_block_i (grid(ib), nf)
      end do

   case (7)

      do ib = 1, nb
         call reverse_block_j (grid(ib), nf)
      end do

   case (8)

      do ib = 1, nb
         call reverse_block_k (grid(ib), nf)
      end do

   case (9)

      dx = 2.;  dy = 4.;  dz = 6.;  scale = 10.

      do ib = 1, nb
         call scale_shift_block (grid(ib), scale, dx, dy, dz, grid(ib)) 
      end do

   case (10)

      do ib = 1, nb
         call swap_block_xyz (grid(ib), 1, 2)
         call swap_block_xyz (grid(ib), 2, 3)
         call swap_block_xyz (grid(ib), 3, 1)
      end do

   case (11)

      if (nb > 1) then
         call update_block (grid(1), nf, .false., .true., grid(2))
      end if

   case default

      write (luncrt, '(a)') 'Valid choices are 1 - 11.'
      go to 10

   end select

   if (icase >= 5) call xyz_write (lunout, formatted, nb, grid, ios)

99 continue

   end program test_block_utilities
