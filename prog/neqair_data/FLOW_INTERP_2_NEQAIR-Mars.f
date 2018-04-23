      PROGRAM FLOWINTERP2NEQAIR

c.... Program to convert FLOW_INTERP LOS data (16-sp Mars model) to a format
c     usable by NEQAIR.
c
c.... Dinesh Prabhu, December 3, 2005

      PARAMETER (il=201, izn=999)

      CHARACTER*1  dig1,dig2,dig3
      CHARACTER*3  lbl
      CHARACTER*28 ofile

      INTEGER      imax(izn), jmax(izn), kmax(izn), nmax(izn)

      REAL         x(il,izn), y(il,izn), z(il,izn), xshk(il)
      REAL         T_in(il,izn), Tv_in(il,izn)
      REAL         Ns_in(il,izn,16)

      REAL         x_out(il)
      REAL         T_out(il), Tv_out(il)
      REAL         Tr_out(il), Te_out(il)
      REAL         Ns_out(il,16), Ns_dum(il,16)
      REAL         Nstot_out(il)

      ns = 16

c.... Read in the LOS data
c     Assumed to be in PLOT3D (Function) ASCII format
c     All quantities are in SI units

      OPEN(UNIT=1,FILE='target.g',FORM='FORMATTED',STATUS='OLD')
         READ(1,*) ng
         READ(1,*) (imax(n),jmax(n),kmax(n),n=1,ng)
         DO n = 1,ng
            READ(1,*) (x(k,n),k=1,kmax(n)),
     &                (y(k,n),k=1,kmax(n)),
     &                (z(k,n),k=1,kmax(n))
         ENDDO
      CLOSE(UNIT=1)

      OPEN(UNIT=1,FILE='target.f',FORM='FORMATTED',STATUS='OLD')
         READ(1,*) ng
         READ(1,*) (imax(n),jmax(n),kmax(n),nmax(n),n=1,ng)
         DO n = 1,ng
            READ(1,*) (T_in(k,n),k=1,kmax(n)),
     &                (Tv_in(k,n),k=1,kmax(n)),
     &                ((Ns_in(k,n,m),k=1,kmax(n)),m=1,ns)
         ENDDO
      CLOSE(UNIT=1)

      np = kmax(1)

c *********************************************************************
c ****** process LOS data
c *********************************************************************
c
c ****** 16-species Mars/Venus model:
c
c....   input species order:   CO2  CO   CO+  C2   N2   O2   NO   NO+  CN   C
c....                         ( 1) ( 2) ( 3) ( 4) ( 5) ( 6) ( 7) ( 8) (9) (10)
c....                          C+   N    O+  O+    Ar   e-
c....                         (11) (12) (13) (14) (15) (16)
c....   output species order:  Ar   CO   N2   O2   NO   NO+  CN   C2   C   C+
c....                         ( 1) ( 2) ( 3) ( 4) ( 5) ( 6) ( 7) ( 8) (9) (10)
c....                           N    N+   O   O+   Ar+  e-
c....                         (11) (12) (13) (14) (15) (16)

      DO n = 1,ng

         IF (n .LT. 10) THEN
            dig3 = "0"
            dig2 = "0"
            dig1 = CHAR(48+n)
         ELSEIF (n .LT. 100) THEN
            dig3 = "0"
            dig2 = CHAR(48+n/10)
            dig1 = CHAR(48+n-(ICHAR(dig2)-48)*10)
         ELSEIF (n .LT. 1000) THEN
            dig3 = CHAR(48+n/100)
            dig2 = CHAR(48+(n-(ICHAR(dig3)-48)*100)/10)
            dig1 = CHAR(48+n-(ICHAR(dig2)-48)*10-(ICHAR(dig3)-48)*100)
         ENDIF
         lbl  = dig3//dig2//dig1
         ofile = "./LOSDATA-NEQAIR/line"//lbl//".los"

c.... Conversion
c     Compute body normal for 3D

         xshk = 0.d0
         DO i = 2,np
            xshk(i) = xshk(i-1) - SQRT((x(i,n)-x(i-1,n))**2 +
     &                                 (y(i,n)-y(i-1,n))**2 +
     &                                 (z(i,n)-z(i-1,n))**2)
         ENDDO

c.... Conversion
c     Distance coordinate converted from m to cm
c     Number densities converted from per m^3 to per cm^3

         DO i = 1,np
            xshk(i) = xshk(i)*100.0
            DO m = 1,ns
               Ns_in(i,n,m) = Ns_in(i,n,m)*1.0E-06
               IF (Ns_in(i,n,m) .LT. 1.0E+03) Ns_in(i,n,m) = 1.0E+03
            ENDDO
         ENDDO

c....      Index reversal 
c      LOS begins in the freestream and ends at the body surface

         DO i = 1,np
            x_out(i)  = xshk(np+1-i)-xshk(np)
            T_out(i)  = T_in(np+1-i,n)
            Tv_out(i) = Tv_in(np+1-i,n)
            DO m = 1,ns
               Ns_dum(i,m) = Ns_in(np+1-i,n,m)
            ENDDO
         ENDDO

c....   Tv limiting -- do not let Tv > T

c          DO i = 1,np
c           Tv_out(i) = min(Tv_out(i),T_out(i))
c         ENDDO

c.... fill in the number densities correctly

         DO i = 1,np
            Tr_out(i) = T_out(i)
            Te_out(i) = Tv_out(i)
            nsout = 16

            Ns_out(i,1)  = Ns_dum(i,1) + Ns_dum(i,15)
            Ns_out(i,2)  = Ns_dum(i,2)
            Ns_out(i,3)  = Ns_dum(i,5)
            Ns_out(i,4)  = Ns_dum(i,6)
            Ns_out(i,5)  = Ns_dum(i,7)
            Ns_out(i,6)  = Ns_dum(i,8)
            Ns_out(i,7)  = Ns_dum(i,9)
            Ns_out(i,8)  = Ns_dum(i,4)
            Ns_out(i,9)  = Ns_dum(i,10)
            Ns_out(i,10) = Ns_dum(i,11)
            Ns_out(i,11) = Ns_dum(i,12)
            Ns_out(i,12) = 1.0d3
            Ns_out(i,13) = Ns_dum(i,13)
            Ns_out(i,14) = Ns_dum(i,14)
            Ns_out(i,15) = Ns_dum(i,3)
            Ns_out(i,16) = Ns_dum(i,16)

            sum = 0.0
            DO m = 1,nsout
               sum = sum + Ns_out(i,m)
            ENDDO
            Nstot_out(i) = sum
         ENDDO

c....      Determine which point marks the beginning of the shock wave

         ipeak = 1
          DO i = 2,np
             IF (ipeak .EQ. 1) THEN
                IF (Tv_out(i) .GT. 5.0d2) ipeak = i
             ENDIF
          ENDDO
          xpeak = x_out(ipeak)


c....      Formatted output
c      This formatted output file has to be pasted into the NEQAIR input file

         OPEN(UNIT=1,FILE=ofile,FORM='FORMATTED',STATUS='UNKNOWN')
            DO i = ipeak,np
               WRITE(1,2000) i-ipeak+1, x_out(i)-xpeak, Nstot_out(i), 
     &                       T_out(i), Tr_out(i), Tv_out(i), Te_out(i)
               WRITE(1,2010) (Ns_out(i,m),m=1,nsout)
               WRITE(1,2020) 
            ENDDO
         CLOSE(UNIT=1)

      ENDDO

 2000 FORMAT(I5,2E15.7,4F10.1)
 2010 FORMAT((6X,4(1PE15.7)))
 2020 FORMAT(" ")

      END
