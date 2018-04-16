c----- This routine pops two arguments off a stack

      subroutine dtsze3 (x, y, ptr, stack)

      integer ptr
      double precision x, y, stack(*)

c----- Pop the stuff off

      if (ptr .gt. 0) then
        x = stack(2 * ptr - 1)
        y = stack(2 * ptr)
        endif
      ptr = ptr - 1
      end

