!     Dummy system utilities for porting Traj_opt to non-SGI systems.
!
      function iargc () ! Count the command line arguments.

      integer iargc

      iargc = 0 ! No command line arguments.

      end function iargc


      subroutine getarg (i, arg) ! Pick off ith command line switch

      character * (*) arg

      end subroutine getarg


      subroutine flush (lun) ! Flush output buffer

      end subroutine flush


!     subroutine second (cputime) ! CPU secs. since start of run

!     cputime = 0.

!     end subroutine second
