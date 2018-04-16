!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module f_io_module

!  This module supplements the earlier xyzq_io package by performing PLOT3D
!  function file I/O at maximum efficiency rather than storing all the functions
!  at a point together in memory, which is often preferrable other than for I/O.
!  It was prompted by the EXTRACT_FUNCTIONS utility, which performs nothing but
!  I/O and doesn't even need an associated grid file.
!
!  Header I/O can still be performed by xyzq_io utilities.
!
!  03/05/2016   D. A. Saunders   Initial implementation, to treat function files
!                                as efficiently as grid files have long been
!                                read and written.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure ! Or equivalent to define "grid_type"

   implicit none

   private

   public :: f_allocate     ! Allocates one flow block as f(mi,mj,mk,nf)
   public :: f_block_io     ! Reads or writes one such flow block

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine f_allocate (block, nf, ios)
!
!     Allocate a single block of a PLOT3D-type flow solution block in the order
!     that maximizes efficient I/O.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent(inout) :: block
      integer, intent (in)            :: nf
      integer, intent (out)           :: ios

!     Execution:

      allocate (block%q(block%mi, block%mj, block%mk, nf), stat=ios)

      end subroutine f_allocate

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine f_block_io (mode, lun, formatted, nf, mi, mj, mk, f, ios)
!
!     Read or write one block of a 3-space PLOT3D function file, with memory
!     storage order matching that of the disk file.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!     Arguments:

      integer, intent (in)    :: mode        ! 1 = read, 2 = write
      integer, intent (in)    :: lun         ! Logical unit number
      logical, intent (in)    :: formatted   ! T|F
      integer, intent (in)    :: nf          ! # flow variables at each grid pt.
      integer, intent (in)    :: mi, mj, mk  ! # grid points in each direction
      real,    intent (inout) :: f(mi,mj,mk,nf)  ! Flow function values
      integer, intent (out)   :: ios         ! 0 = no error

!     Execution:

      select case (mode)
         case (1)
            if (formatted) then
               read (lun, *,  iostat=ios) f
            else
               read (lun,     iostat=ios) f
            end if
         case (2)
           if (formatted) then
               write (lun, '(6es18.10)', iostat=ios) f
            else
               write (lun,               iostat=ios) f
            end if
         case default
            ios = 999
      end select

      end subroutine f_block_io

   end module f_io_module
