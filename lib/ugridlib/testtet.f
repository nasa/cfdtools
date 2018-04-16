      program testtet

c  Test isintet and calcphi.
c
c  02-JUL-1999  Scott D. Thomas, Raytheon ITSS Corporation
c               Contract NAS 2-98080/172, NASA Ames Research Center.
c
c  11-AUG-2000  David Saunders  Removed the mixed precision for use with
c                               one-precision ISINTET and CALCPHI routines.
c                               (Use -r4 or -r8 compiler switch.)

      implicit none

      integer ndc(4)
      integer m,n1,n2,n3,n4,nmin,m1,m2,m3
      real    xyz(3,4),xyzd(3)
      real    phi0,phi(3)
      real    x,y,z,epsilon
      real    a(3),b(3),c(3),d(3),h
      logical isintet

      data xyz / 0.0, 0.0, 0.0,
     .           1.1, 0.0, 0.0,
     .           0.0, 0.9, 0.0,
     .           0.0, 0.0, 1.2 /
      data ndc / 1, 2, 3, 4 /

      n1 = 1
      n2 = 2
      n3 = 3
      n4 = 4
      epsilon = 1.e-6
      write(6,*) 'TESTTET:  EPSILON =', epsilon
c
c (4) n3,n4,n1 face.
c
      nmin = min(n3,n4,n1)
      if(n3.eq.nmin) then
        m1 = n3
        m2 = n4
        m3 = n1
      else if(n4.eq.nmin) then
        m1 = n4
        m2 = n1
        m3 = n3
      else
        m1 = n1
        m2 = n3
        m3 = n4
      endif
      do m=1,3
        a(m) = xyz(m,m2) - xyz(m,m1)
        b(m) = xyz(m,m3) - xyz(m,m1)
      end do
      c(1) = a(2)*b(3)-a(3)*b(2)
      c(2) = a(3)*b(1)-a(1)*b(3)
      c(3) = a(1)*b(2)-a(2)*b(1)
      do m=1,3
        d(m) = xyz(m,n2)-xyz(m,m1)
      end do
      h    = c(1)*d(1)+c(2)*d(2)+c(3)*d(3)

      write(6,*) 'n1,n2,n3,n4 =',n1,n2,n3,n4
      write(6,*) 'xyz(1..3,n1)=',xyz(1:3,n1)
      write(6,*) 'xyz(1..3,n2)=',xyz(1:3,n2)
      write(6,*) 'xyz(1..3,n3)=',xyz(1:3,n3)
      write(6,*) 'xyz(1..3,n4)=',xyz(1:3,n4)
      write(6,*) 'Volume #4   =',h/6.0
c
   10 continue
      write(6,*)
      write(6,*) 'Enter x,y,z'
      read(*,*,err=20,end=20) x,y,z
      write(6,*) '(x,y,z)     =',x,y,z

      if (isintet(x,y,z,xyz,ndc,epsilon)) then

        call calcphi (x,y,z,xyz,ndc,phi,epsilon)
c
c  Check to see whether it worked to within epsilon.
c
        phi0 = 1.0 - (phi(1)+phi(2)+phi(3))
        do m=1,3
          xyzd(m) = phi0   * xyz(m,n1) +
     .              phi(1) * xyz(m,n2) +
     .              phi(2) * xyz(m,n3) +
     .              phi(3) * xyz(m,n4)
        end do
        if(abs(x-xyzd(1)).gt.epsilon.or.
     .     abs(y-xyzd(2)).gt.epsilon.or.
     .     abs(z-xyzd(3)).gt.epsilon) then
          write(6,*) 'TESTTET:  Application of phi to vertices failed.'
          write(6,*) 'n1,n2,n3,n4 =',n1,n2,n3,n4
          write(6,*) '(x,y,z)     =',x,y,z
          write(6,*) 'xyzd(1..3)  =',xyzd
          write(6,*) 'epsilon     =',epsilon
          write(6,*) 'xyz(1..3,n1)=',xyz(1:3,n1)
          write(6,*) 'xyz(1..3,n2)=',xyz(1:3,n2)
          write(6,*) 'xyz(1..3,n3)=',xyz(1:3,n3)
          write(6,*) 'xyz(1..3,n4)=',xyz(1:3,n4)
          write(6,*) 'Volume #4   =',h/6.0
          write(6,*) 'phi(0..3)   =',phi0,phi
          stop 'ERROR'
        else
          write(6,*) 'TESTTET:  Application of phi to vertices worked.'
          write(6,*) '(x,y,z)     =',x,y,z
          write(6,*) 'xyzd(1..3)  =',xyzd
          write(6,*) 'epsilon     =',epsilon
          write(6,*) 'phi(0..3)   =',phi0,phi
        endif
      else
        write(6,*) 'TESTTET:  That point is outside the cell.'
      endif
      goto 10
   20 continue

      end program testtet
