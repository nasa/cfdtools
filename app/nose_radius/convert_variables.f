********************************************************************************
*
      program convert_variables
*
*     Convert NOSE_RADIUS shape function results to the forms required for
*     input to NOSE_RADIUS and AEROSURF.
*
*     INPUT FILE (1):  'convert_variables.in' (probably from 'nose_radius.out')
*
*     DESIGN VARIABLES AT DESIGN ITERATION  11
*       10
*       1  4.443512173477E+00    2  1.587512222995E+00    3  8.345607033491E-01
*       :   :                    :   :                    :   :
*       <any number of pairs on any number of lines, read list-directed>
*
*     INPUT FILE (2):  'nose_radius.inp' corresponding to 'nose_radius.out'
*
*     OUTPUT FORMAT (1):  'nose_radius.in'
*
*     #   VTYPE   UP  LO  POWER  CENTER  XMIN  XMAX     V  VSCALE   H    BL  BU
*     1  'SIN1 '  0.  1.    4.0   0.200   0.0   1.0 4.444   0.001 1.E-5 -10. 10.
*     2  'SIN1 '  0.  1.    4.0   0.250   0.0   1.0 1.588   0.001 1.E-5 -10. 10.
*     3  'SIN1 '  0.  1.    4.0   0.300   0.0   1.0 8.346   0.001 1.E-5 -10. 10.
*     :     :     :
*
*     OUTPUT FORMAT (2):  'aerosurf.in'
*
*     #  VTYPE   UTYPE UP LO XPWR  XCEN XMIN XMAX LOFT KCEN KMIN KMAX    V  VSCALE  H  COMP
*     1 'EXP  ' 'COSR' 0. 1.  30. 0.002 0.00 1.00   3.  60.  45.  60. -0.34 0.001 1.E-3  1
*     2 'EXP  ' 'COSR' 0. 1.  30. 0.004 0.00 1.00   3.  60.  45.  60.  0.31 0.001 1.E-3  1
*     3 'SIN2 ' 'COSR' 0. 1.   2. 0.003 0.00 1.00   3.  60.  45.  60.  0.12 0.001 1.E-3  1
*     :     :     :
*
*     Some columns of 'aerosurf.in' will need editing to loft the results of
*     NOSE_RADIUS across a wing.
*
*     11/05/03  DAS  Initial implementation.
*
*     Author:  David Saunders, ELORET/NASA Ames Research Cntr, Moffett Field, CA
*
********************************************************************************

      implicit none

*     Constants:

      integer, parameter ::
     >   lunin1 = 1, lunin2 = 2, lunout = 3

*     Variables:

      integer
     >   i, icomp, n, nvar

      real
     >   bl, bu, h, lo, up, vscale, vtemp, xcen, xmax, xmin, xpower,
     >   zcen, zmax, zmin, zpower

      real, allocatable ::
     >   v(:)

      character
     >   utype * 4, vtype5 * 5, vtype6 * 6

*     Execution:

      open (unit=lunin1, file='convert_variables.in', status='old')
      open (unit=lunin2, file='nose_radius.inp',      status='old')
      open (unit=lunout, file='nose_radius.in',       status='unknown')

      read (lunin1, '(a)') ! Skip title
      read (lunin1, *) nvar

      allocate (v(nvar))

      read (lunin1, *) (n, v(i), i = 1, nvar)

      close (lunin1)

*     Plug the variables back into a NOSE_RADIUS input fragment:

      do i = 1, 11
         read (lunin2, '(a)') ! Skip the first 11 lines of 'nose_radius.inp'
      end do

      write (lunout, '(2a)')
     >   ' #   VTYPE   UP  LO  POWER  CENTER  XMIN  XMAX         V',
     >   '             VSCALE    H       BL    BU'

      do i = 1, nvar
         read (lunin2, *) n, vtype5, up, lo, xpower, xcen, xmin, xmax,
     >                    vtemp, vscale, h, bl, bu
         write (lunout,
     >      '(i2, 3a, 2f4.0, f7.1, f8.3, 2f6.3, 1p, e21.13, 0p, f8.4,
     >        1p, e8.1, 0p, 2f6.1)')
     >      i, '  ''', vtype5, '''', up, lo, xpower, xcen, xmin, xmax,
     >                    v(i),  vscale, h, bl, bu
      end do

      close  (lunout)
      rewind (lunin2)

*     Repeat for AEROSURF-type input:

      open (unit=lunout, file='aerosurf.in',       status='unknown')

      do i = 1, 11
         read (lunin2, '(a)')
      end do

      write (lunout, '(3a)')
     >   ' #  VTYPE    UTYPE  UP  LO XPWR  XCEN  XMIN  XMAX  LOFT',
     >   '  ZCEN  ZMIN  ZMAX         V             VSCALE    H   COMP'

      zpower = 3.
      icomp  = 1

      do i = 1, nvar
         read (lunin2, *) n, vtype6, up, lo, xpower, xcen, xmin, xmax,
     >                    vtemp, vscale, h, bl, bu
         write (lunout,
     >      '(i2, 3a, 2f4.0, f5.1, 3f6.3, f6.2, 3f6.3, 1p, e21.13, 0p,
     >        f8.4, 1p, e8.1, i4)')
     >      i, ' ''', vtype6, ''' ''RCOS''', up, lo, xpower, xcen, xmin,
     >      xmax, zpower, zcen, zmin, zmax, v(i),  vscale, h, icomp 
      end do

      close (lunout)
      close (lunin2)

      end program convert_variables
