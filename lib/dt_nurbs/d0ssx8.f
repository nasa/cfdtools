      subroutine D0SSX8 (upts, vpts, spts, tpts, icptr, ix, work)

c----- Reverse the lists

      double precision upts(*), vpts(*), spts(*), tpts(*), work(*)
      integer icptr(*), ix, size
      size = icptr(ix+1) - icptr(ix)
      call dcopy (size, upts(icptr(ix)), -1, work, 1)
      call dcopy (size, work, 1, upts(icptr(ix)), 1)
      call dcopy (size, vpts(icptr(ix)), -1, work, 1)
      call dcopy (size, work, 1, vpts(icptr(ix)), 1)
      call dcopy (size, spts(icptr(ix)), -1, work, 1)
      call dcopy (size, work, 1, spts(icptr(ix)), 1)
      call dcopy (size, tpts(icptr(ix)), -1, work, 1)
      call dcopy (size, work, 1, tpts(icptr(ix)), 1)
      end
