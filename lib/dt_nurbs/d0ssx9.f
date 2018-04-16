      subroutine D0SSX9 (t, c, fval, fjac, ier)

c----- The next subroutine tracks the contours

      integer nwork
      parameter (nwork = 200)
      double precision t(*), c(*), fval(*), fjac(*), work(nwork)
      double precision temp(3)
      integer ier, st2, ider(3),i
      integer ier1, ier2, ier3, ier4, ier5, ier6
      data ider /1, 0, 1/

      do 10 i = 1, 4
        if ( t(i) .lt. 0.d0 ) t(i) = 0.d0
        if ( t(i) .gt. 1.d0 ) t(i) = 1.d0
   10 continue

      st2 = c(1)
      call dtnpvl (t, 1, c(2), work, nwork, fval, ier1)
      call dtnpvl (t(3), 1, c(st2), work, nwork, temp, ier2)
      call daxpy (3, -1.0d0, temp, 1, fval, 1)
      call dtnpdr (t, 1, ider(1), c(2), work, nwork, fjac(1), ier3)
      call dtnpdr (t, 1, ider(2), c(2), work, nwork, fjac(4), ier4)
      call dtnpdr (t(3), 1, ider(1), c(st2), work, nwork, fjac(7), ier5)
      call dscal (3, -1.0d0, fjac(7), 1)
      call dtnpdr (t(3), 1, ider(2), c(st2), work, nwork,
     +             fjac(10), ier6)
      call dscal (3, -1.0d0, fjac(10), 1)

      if ( ier1 + ier2 + ier3 + ier4 + ier5 + ier6 .lt. 0 ) ier = -200

      end
