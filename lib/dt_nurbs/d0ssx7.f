      subroutine D0SSX7 (upts, vpts, spts, tpts, icptr, lsize,
     +                   from, to, fpt, gpt)

c----- Insert material onto a list at a given point

      double precision upts(*), vpts(*), spts(*), tpts(*),
     +                 fpt(2), gpt(2)
      integer icptr(*), lsize, from, to, mcopy, ix, this

      if (from .gt. 0) then
        mcopy = icptr(from+1) - icptr(from) + 1
      else
        mcopy = 1
        endif
      do 10 ix = lsize,to,-1
        this = icptr(ix+1) - icptr(ix)
        call d0xcpy (this, upts, icptr(ix), -1, icptr(ix)+mcopy, -1)
        call d0xcpy (this, vpts, icptr(ix), -1, icptr(ix)+mcopy, -1)
        call d0xcpy (this, spts, icptr(ix), -1, icptr(ix)+mcopy, -1)
        call d0xcpy (this, tpts, icptr(ix), -1, icptr(ix)+mcopy, -1)
        icptr(ix+2) = icptr(ix+1) + mcopy
 10     continue
      icptr(to+1) = icptr(to) + mcopy
      if (from .gt. 0) then
        call d0xcpy (mcopy, upts, icptr(from), 1, icptr(to), 1)
        call d0xcpy (mcopy, vpts, icptr(from), 1, icptr(to), 1)
        call d0xcpy (mcopy, spts, icptr(from), 1, icptr(to), 1)
        call d0xcpy (mcopy, tpts, icptr(from), 1, icptr(to), 1)
        upts(icptr(to)+mcopy-1) = fpt(1)
        vpts(icptr(to)+mcopy-1) = fpt(2)
        spts(icptr(to)+mcopy-1) = gpt(1)
        tpts(icptr(to)+mcopy-1) = gpt(2)
      else
        upts(icptr(to)) = fpt(1)
        vpts(icptr(to)) = fpt(2)
        spts(icptr(to)) = gpt(1)
        tpts(icptr(to)) = gpt(2)
        endif
      if (from .eq. 0) return
      do 20 ix = from+1,lsize+1
        this = icptr(ix+1) - icptr(ix)
        call d0xcpy (this, upts, icptr(ix), 1, icptr(ix-1), 1)
        call d0xcpy (this, vpts, icptr(ix), 1, icptr(ix-1), 1)
        call d0xcpy (this, spts, icptr(ix), 1, icptr(ix-1), 1)
        call d0xcpy (this, tpts, icptr(ix), 1, icptr(ix-1), 1)
        icptr(ix) = icptr(ix-1) + this
 20     continue
      end
