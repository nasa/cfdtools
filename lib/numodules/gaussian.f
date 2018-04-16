c-------------------------------------------------------------------------------
c
      function gauss (iseed)
c
c --- Function to return Gaussian distributed points.
c     As received, this generates the same number every call after the first.
c     Why was it so poorly documented??
c     The functionality has been clarified as follows:
c
c     Set iseed to a negative value to initialize or reinitialize a sequence.
c
c     The output values belong to a normal distribution with mean zero and unit
c     standard deviation, meaning they lie in [-3.xxx, +3.xxx] with Gaussian-
c     type bunching towards zero, as might be handy for applying noise to an
c     analytic dataset.

      implicit none

      integer, intent (inout) :: iseed
      real*8       :: two, radius, theta, gauss1, gauss2
      real*8, save :: twopi
      real*8       :: gauss, ran3
      external     :: ran3
c
      if (iseed < 0) then  ! Avoid unnecessarily recalculating this
         twopi=4.d0*asin(1.d0)
      end if

      two=-2.d0*log(1.d0 - ran3(iseed))
      radius=sqrt(two)
      theta=twopi*ran3(iseed)
      gauss1=radius*cos(theta)
      gauss2=radius*sin(theta)

      if (ran3(iseed) < 0.5d0) then
         gauss=gauss1
      else
         gauss=gauss2
      end if

      end function gauss

c-------------------------------------------------------------------------------
c
      FUNCTION RAN3 (ISEED)
c
c     ran3.f from Numerical Recipes
      IMPLICIT REAL*8(A-H,M,O-Z)
      integer, intent (inout) :: ISEED
      real*8, parameter :: MBIG=4000000.d0,MSEED=1618033.d0,MZ=0.d0
      real*8, parameter :: FAC=2.5d-7
      real*8 MA(55)
      real*8 RAN3
      save
      DATA IFF /0/
c
      IF(ISEED.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        MJ=MSEED-IABS(ISEED)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1
        DO 11 I=1,54
          II=MOD(21*I,55)
          MA(II)=MK
          MK=MJ-MK
          IF(MK.LT.MZ)MK=MK+MBIG
          MJ=MA(II)
11      CONTINUE
        DO 13 K=1,4
          DO 12 I=1,55
            MA(I)=MA(I)-MA(1+MOD(I+30,55))
            IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
12        CONTINUE
13      CONTINUE
        INEXT=0
        INEXTP=31
        ISEED=1
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC

      END FUNCTION RAN3
