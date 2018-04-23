      program convert2d
c
c     Split a single-block 64-bit (no blanking) 2D PLOT3D grid file
c     into a 32-bit grid file (unformatted in both cases).
c
c     07/01/97  D.Saunders  Adapted from convert_grid (3D), which was
c                           adapted from Mark Rimlinger's
c                           multiblock --> single-block converter.

      implicit   real*8 (a-h,o-z)
      parameter (mpts=50000,mgrid=1)
      dimension  x(mpts,2),ni(mgrid),nj(mgrid)
      logical    toobig
      character  file(mgrid)*40, infile*40

      write(*,'(a, $)') ' Input 64-bit PLOT3D grid filename? '
      read (*, '(a)', end=99) infile

      open (unit=3,file=infile,status='old',form='unformatted')
      ngrid = 1
      read(3) (ni(ig), nj(ig), ig = 1, ngrid)
      if (ngrid.gt.mgrid) then
         write (*, 2) ngrid, mgrid
    2    format('File has too many grids: ',i7,', dimensioned ',i7)
         goto 99
      end if

      toobig = .false.
      do ig= 1,ngrid
         write (*, 4) ig,ni(ig),nj(ig)
    4    format ('    Grid ', i2, ' is dimensioned ', 2i7)
         npts  = ni(ig)*nj(ig)
         if (npts.gt.mpts) then
            write (*, 6) npts, mpts
    6       format ('Grid has too many points: ',i7,', dimensioned ',i7)
            toobig = .true.
         end if
      end do

      if (toobig) go to 99

      do ig = 1, ngrid
         write (*, '(a, i2, a, $)') ' Output grid filename #', ig, ': '
         read (*, '(a)', end=99) file(ig)
      end do

c     Generate single grid files.  If the filename is blank, skip this grid.

      do 30 ig = 1, ngrid
         npts = ni(ig)*nj(ig)
         read (3) ((x(i,nx), i=1,npts), nx=1,2)
         if (file(ig).eq.' ') go to 30

         open(unit=2,file=file(ig),status='unknown',form='unformatted')
         write (2) ni(ig),nj(ig)
         write (*, '(a, i2)') ' Writing output grid # ', ig
         write(2)  ((sngl (x(i,nx)), i=1,npts), nx=1,2)
         close (unit=2)
   30 continue

      close (unit=3)

   99 continue
      stop
      end
