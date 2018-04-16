c----- This routine pushes two arguments onto a stack

      subroutine dtsze2 (x, y, ptr, stack)

      integer ptr
      double precision x, y, stack(*)

c----- Push the stuff on

      ptr = ptr + 1
      stack(2 * ptr - 1) = x
      stack(2 * ptr) = y
      end

