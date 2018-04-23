!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine read_surface_patch (lunsrf, surface_patch, ios)
!
!  Prompt for and read a structured surface patch file, assumed to be formatted.
!
!  03/16/04  David Saunders  Modularize part of the I/O for use with other
!                            utilities for morphing grids to damaged surfaces.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Data structure definition:

   use grid_block_structure

   implicit none

!  Arguments:

   integer, intent (in) :: lunsrf                  ! Logical unit, opened and
                                                   ! closed here
   type (grid_type), intent (out) :: surface_patch ! Unallocated on input;
                                                   ! (x,y,z)s defined on output
   integer, intent (out) :: ios                    ! Status code:
                                                   ! 0 = success, else failure
!  Local constants:

   integer, parameter :: lunkbd = 5, luncrt = 6    ! Keyboard & screen units
   logical, parameter :: formatted = .TRUE.        ! Unformatted is unlikely

!  Local variables:

   integer :: n
   character :: filename * 80

!  Execution:

   write (luncrt, '(/, a)', advance='no') ' Damage/repair surface file name: '

   do ! Until success or EOF (^D)
      read (lunkbd, '(a)', iostat=ios) filename
      if (ios < 0) go to 99            ! Single return

      open (lunsrf, file=filename, status='old', iostat=ios)
      if (ios == 0) exit
      write (luncrt, '(a)') ' Cannot open that file.  Try again.'
   end do

   read (lunsrf, *) n          ! Should be 1 (multiblock form)

   if (n /= 1) rewind (lunsrf) ! Try single-block form

   read (lunsrf, *) surface_patch%ni, surface_patch%nj

   surface_patch%nk = 1

   call allocate_block (surface_patch)

   n = surface_patch%ni * surface_patch%nj

   call block_io (1, lunsrf, formatted, n, &
                  surface_patch%x, surface_patch%y, surface_patch%z, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a, i4)') ' Trouble reading surface patch.  ios:', ios
      go to 99
   end if

   close (lunsrf)

99 return

   end subroutine read_surface_patch
