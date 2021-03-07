c-----------------------------------------------------------------------
c
      program template
c
c  Original Description:
c  ---------------------
c
c  Make a FLO107MB connection file template from a PLOT3D grid file.
c  The name 'template' is used because most but not all of the block
c  face connectivity/boundary condition information can be deduced
c  from just the grid, although some educated guesses are applied to
c  the faces that seem not to match other block faces.
c
c  The grid file may be single or mgrid, real*4 or real*8, though it
c  only makes sense to process multiblock grids.
c
c  Use "template -v" for verbose template.con file (last column will
c  indicate what fraction of quads of discovered abutting faces match
c  within epsilon distance at all four corners.
c
c  (This version turns on "-v" explicitly.)
c
c  There are nb blocks, the nth block has dimensions nj(n),nk(n),nl(n).
c  Each block becomes a network of six patches, one for each side.  The
c  sides are numbered 1 to 6 as:  jmin, jmax, kmin, kmax, lmin, lmax.
c  Each quad face of each side of each block is stored as a four-tuple
c  of integers that point into xyz.  Component numbers are the sides
c  1-6 for the first block, then 7-12 for the second block, and so on
c  to 6*nb.  An index array iaf into ndf points to the nth component's
c  quads which range from iaf(n) to iaf(n+1)-1.  The xtrema array is
c  used to store the xyz minmax values for each component (for each
c  side of each block).  The component xyz extrema are used to match
c  sides within epsilon.  ihand(n) is the handedness of the nth block,
c  right=1, left=-1; it is found roughly by adding some tet volumes and
c  checking the sign of the result.
c
c  Scott D. Thomas, Sterling Software, Inc.
c  Contract NAS 2-13210/6/5, NASA Ames Research Center.
c  25-JUN-1998..30-JUN-1998.
c
c  Scott D. Thomas, Raytheon STX Corporation
c  Contract NAS 2-98080/53, NASA Ames Research Center.
c  1-JUL-1998..9-JUL-1998.
c
c  Adaptation within the Space Technology Division at Ames:
c  --------------------------------------------------------
c
c  DPLR Input and Interface Files:
c
c     This version writes a GASP control file and two DPLR control files
c     as well as the 'template.con' file used by FLO107MB.  It also handles
c     formatted grids as well as the original real*4 or real*8 unformatted.
c
c  Ancillary control file 'template.inp.2' usage:
c
c     No 'template.inp' control file is in the picture here, but the name
c     'template.inp.2' is used for this optional control file to be consistent
c     with other grid manipulation utilities.
c
c     The leading keyword on the first of each optional pair of input lines is
c     converted to upper case here, and only the first 4 characters are used.
c
c     Pairs of control lines are read until EOF, so the order of the pairs
c     does not matter.
c
c  Optional usage (1) (obsolete with DPLR V4.03, since ntx/y/z are gone now):
c
c     Control the type of sequencing (grid coarsening) specified in the
c     'dplr.inputs.2' output.  Originally, 2 2 1 coarsening was hard-coded.
c     Now, a pair of lines as shown can produce whatever is really needed.
c     The default coarsening remains 2 2 1.  Example:
c
c     Sequencing controls
c     4 4 2              ! Coarsening factors in the i, j, k directions
c
c  Optional usage (2):
c
c     Force boundary condition changes in the two DPLR control files as needed
c     for proper analysis of cavity or plug repair cases.
c
c     Cavity block(s)    ! Any BC 1 or 2 will be changed to BC 26 (wall)
c     6                  ! Any reasonable list of block numbers
c
c  or
c
c     Plug blocks        ! BC 2 at jmin will be changed to BC 26 (wall)
c     10:13              ! Any reasonable list of block numbers
c
c  History beyond Scott's original implementation:
c  -----------------------------------------------
c
c  [Note that Scott uses a pre-Fortran 90 form of dynamic memory allocation
c   that appears more trouble than it is worth to replace with the modern
c   forms of allocate and deallocate.]
c
c  13-MAY-2004  David Saunders  Command line arguments don't work with
c               ELORET Corp./   Intel's ifc compiler;
c               NASA Ames       the grid may be formatted now.
c  14-MAY-2004    "             Added output of two DPLR-related files,
c                               derived from the saved FLO107MB data.
c  17-MAY-2004    "             Added output of GASP format as well.
c  23-OCT-2004    "             Dinesh asked for dplr.inputs.2 with the
c                               i and j block dimensions halved.
c  23-DEC-2004    "             DPLR-related BCs for a grid with more
c                               than one layer of blocks were being
c                               affected by the one-layer case, but
c                               we shouldn't make assumptions about
c                               k = 1 or nk faces unless matching faces
c                               have not been found.
c  13-JAN-2005    "             Changed "init1" to "initi" so that the
c                               Perl script for updating DPLR control
c                               files doesn't misbehave.
c  26-MAY-2005    "             Option to output complete dplr.inputs:
c                               if 'sample.inputs' is present, its
c                               header and trailer are transcribed.
c  14-NOV-2005    "             Distinguish the symmetry planes in the
c                               dplr.inputs files (not always y = 0).
c  18-NOV-2005    "             With DPLR V3.04 or later now the usual
c                               choice, use BC 2 where BC 3 was used
c                               (specified inflow or supersonic outflow,
c                               rather than just supersonic outflow).
c  03-JAN-2006    "             The symmetry plane test should not use
c                               too loose a tolerance. Make it <= 0.001.
c  11-APR-2006    "             Introduced ancillary control file
c                               'template.inp.2' to allow automation of
c                               BC changes in 'dplr.inputs' for cavity
c                               or plug repair cases.
c  13-APR-2006    "             Chun Tang asked for a way to control the
c                               type of grid sequencing in 'dplr.inputs.2',
c                               so that has been added to the options
c                               enabled by 'template.inp.2' (which still
c                               need not be present).
c  10-MAY-2006    "             DPLR V3.05's new grid tailoring controls
c                               required changing how the header lines
c                               from 'sample.inputs' are transcribed.
c  27-MAY-2009    "             DPLR V4 expects "ibadpt" on the end of the
c                               first line of BC info.  DPLR V3 will not be
c                               affected by its presence.
c  27-MAR-2010    "             Long chemistry file names were being
c                               truncated beyond column 80.
c  11-FEB-2014  DAS, ERC Inc./  The dplr.inputs output file now suits
c               NASA Ames       DPLR V4.03: the grid block dimensions of
c                               V4.02 and earlier are no longer present.
c                               dplr.inputs.2 is therefore redundant.
c-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

      parameter (MSTRINGS=576) ! # different ways to position two faces;
                               ! see subroutine diddle.
      character string*20, strings(MSTRINGS)*20
      character argument*80, inmeshfn*80, line*132, outconfn*80
      character aiformat*7, keyword*4
      logical verbose, formatted, header_trailer, ismgrid, matches
      logical cavity, plug, cr, eof
      integer ib, ibc(6), icon(4,6), ios, izone1(8), izone2(8)
      integer iseq, jseq, kseq
      integer, allocatable, dimension (:) :: iblist
      integer, allocatable, dimension (:,:,:) :: iconsave

      data machisz,machrsz,machr4sz/4,8,4/

c  Convention for dynamic array naming (using mempack):
c  If array name is x then its pointer is px and its handle is kx.

      dimension   nj(1),   nk(1),   nl(1),   ihand(1)
      pointer   (pnj,nj),(pnk,nk),(pnl,nl),(pihand,ihand)
      dimension   xyz(3,1),  xtrema(6,1)
      pointer   (pxyz,xyz),(pxtrema,xtrema)
      real*4      xyz4(3,1)
      pointer   (pxyz4,xyz4)
      real*8      xyz8(3,1)
      pointer   (pxyz8,xyz8)
      dimension   ndf(4,1),  iaf(1),    ipoint(1)
      pointer   (pndf,ndf),(piaf,iaf),(pipoint,ipoint)
      dimension   ndc(4,1)
      pointer   (pndc,ndc)

      indx(j,k,l,n) = j+nj(n)*((k-1)+nk(n)*(l-1))

c  Execution:

c  Get file names and epsilon (used to compare xyz coordinates).

      write (*,'(a)')
     .  ' TEMPLATE:  Make a FLO107-MB connection file template.',
     .  '            This version also writes GASP & DPLR controls.'

c     verbose = .false.
c     if(iargc().eq.1) then
c       call getarg(1,argument)
c       if(argument.eq.'-v') verbose = .true.
c     endif

c     write(*,'(a)',advance='no') ' Verbose mode? [y|n]: '
c     read (*,*) string(1:1)
c     verbose = string(1:1) == 'y'

      verbose = .true.

      if(verbose) then
        write(*,*) 'Verbose option:  show quad hit fractions,'//
     .             ' in last column of connectivity file.'
      else
        write(*,*) 'Use "template -v" for a verbose template.con '//
     .             'file with quad hit fractions.'
      endif

      write(*,'(/, a)',advance='no') ' Input grid name: '
      read(*,1010) inmeshfn
 1010 format(a)
      if(lenstr(inmeshfn).eq.0) inmeshfn='xyz.unf'

ccc   write(*,1010,advance='no') ' Name of output connectivity file: '
ccc   read(*,1010) outconfn
      outconfn = ' '
      if(lenstr(outconfn).eq.0) outconfn='template.con' ! FLO107-MB version

      epsilon = 1.e-4
      call readr (6,
     >   'Distance tolerance for face matching: [<cr> = 0.0001] ',
     >   5, epsilon, cr, eof)

    1 continue

c  Initialize the strings array which maps first (n) and second (m)
c  block abutting face numbers (nf and mf, respectively), face
c  connection info (ia,ib) and first block handedness (ihand(n)), or
c  (nf,mf,ia,ib,ihand(n)), to the first (n) block face number (nf),
c  block connection info, and first block handedness (ihand(n)), or
c  (nf,(icon(i,nf),i=2,4),ihand(n)).  See diddle source code for
c  detailed documentation.

      call diddle(strings,nstrings,MSTRINGS)

c  Check and read the mesh, allocating memory as you go.

      n = lenstr(inmeshfn)
      write(*,'(3a)') ' Opening file ', inmeshfn(1:n), ' for read ...'
      formatted = .true.
      ismgrid   = .true.  ! Multiblock
      nrpgrid   = 8       ! Double precision
      open (7,file=inmeshfn(1:n),form='FORMATTED',status='OLD',
     .      iostat=ios)
      if (ios .ne. 0) then
        write (*,'(2a)') ' No such file: ', inmeshfn(1:n)
        stop
      end if
      read (7, *, iostat=ios) nb
      if (ios .ne. 0 .or. nb .le. 0 .or. nb .gt. 500) then
        write (*,*) ' Formatted read of # blocks gives nb & iostat = ',
     .    nb, ios
        write (*,*) ' Try again as unformatted.'
        formatted = .false.
        close (7)
        open (7,file=inmeshfn(1:n),form='UNFORMATTED')
        call gridprop(7,ismgrid,nrpgrid)
        if(ismgrid) then
          read (7) nb
        else
          nb = 1
        endif
      end if
      write(*,*) '   NB: ',nb
      if(nb.le.1) stop ' nblocks <= 1'

c     Now that we know the number of blocks, set up for the optional
c     ancillary control file:

      iseq   = 2
      jseq   = 2
      kseq   = 1
      cavity = .false.
      plug   = .false.

      open (20, file='template.inp.2', status='old', iostat=ios)

      do while (ios == 0)

        read (20, *, iostat=ios) keyword

        if (ios /= 0) then
          close (20)
          exit
        end if

        call upcase (keyword)

        select case (keyword)

        case ('SEQU')  ! Grid sequencing controls
          read (20, *) iseq, jseq, kseq

        case ('CAVI')  ! Damage case with cavity blocks
          cavity = .true.

        case ('PLUG')  ! Plug repair case with thin blocks around the plug
          plug = .true.

        end select

        if (cavity .or. plug) then
          allocate (iblist(nb))
          iblist(:) = 0
          niblist   = nb ! Upper limit

          if (cavity) then
            call rdlist
     >        (6, 'Reading template.inp.2 for cavity blocks. ',
     >         20, niblist, iblist)
          else if (plug) then
            call rdlist
     >        (6, 'Reading template.inp.2 for blocks around plug. ',
     >         20, niblist, iblist)
          end if

          do i = 1, niblist
            ib = iblist(i)
            if (ib < 1 .or. ib > nb) then
              write (6, '(a, i6)') ' Bad block # in template.inp.2:', ib
              cavity = .false.
              plug   = .false.
            end if
          end do
        end if

      end do  ! Next pair of control lines if any

      numr = max(1,(nb+1)/(machrsz/machisz))
      call ema(pnj,numr,knj,'nj',1)
      call ema(pnk,numr,knk,'nk',2)
      call ema(pnl,numr,knl,'nl',3)
      if (formatted) then
        read (7, *) (nj(n),nk(n),nl(n),n=1,nb)
      else
        read (7) (nj(n),nk(n),nl(n),n=1,nb)
      end if

      npts  = 0
      nnode = 0
      nquad = 0
      do n=1,nb
        if(nj(n).lt.1.or.nk(n).lt.1.or.nl(n).lt.1) then
          write(*,*) ' n,nj,nk,nl=',n,nj(n),nk(n),nl(n)
          stop ' NJ,NK,NL'
        endif
        njkl = nj(n)*nk(n)*nl(n)
        npts = npts + njkl
        do l=1,nl(n)
          do k=1,nk(n)
            do j=1,nj(n)
              if(j.eq.1.or.j.eq.nj(n).or.
     .           k.eq.1.or.k.eq.nk(n).or.
     .           l.eq.1.or.l.eq.nl(n)) nnode = nnode + 1
            end do
          end do
        end do
        nquad = nquad + 2*(((nk(n)-1)*(nl(n)-1)) +
     .                     ((nj(n)-1)*(nl(n)-1)) +
     .                     ((nj(n)-1)*(nk(n)-1)))
        if(n.eq.1) then
          njklmax = njkl
        else
          njklmax = max(njklmax,njkl)
        endif
      end do

      ncomp = 6*nb
      write(*,*) ' NPTS: ',npts
      write(*,*) 'NNODE: ',nnode
      write(*,*) 'NCOMP: ',ncomp
      write(*,*) 'NQUAD: ',nquad
c
c  Allocate storage for all points (xyz), quads (ndf), for the index
c  array (iaf) that points component-wise into the quad list.
c  And for the corner index (ndc) array.
c
      numr = max(1,(3*nnode+1)/(machrsz/machrsz))
      call ema(pxyz,numr,kxyz,'xyz',4)
      numr = max(1,(4*nquad+1)/(machrsz/machisz))
      call ema(pndf,numr,kndf,'ndf',5)
      numr = max(1,(ncomp+1+1)/(machrsz/machisz))
      call ema(piaf,numr,kiaf,'iaf',6)
      numr = max(1,(4*ncomp+1)/(machrsz/machisz))
      call ema(pndc,numr,kndc,'ndc',7)
c
c  Allocate space for the handedness array (ihand).
c
      numr = max(1,(nb+1)/(machrsz/machisz))
      call ema(pihand,numr,kihand,'ihand',8)
c
c  Allocate space for reading each block (xyz4 or xyz8).
c
      if(nrpgrid.eq.4) then
        numr = max(1,(3*njklmax+1)/(machrsz/machr4sz))
        call ema(pxyz4,numr,kpxyz4,'pxyz4',9)
      else
        numr = max(1,(3*njklmax+1)/(machrsz/machrsz))
        call ema(pxyz8,numr,kpxyz8,'pxyz8',10)
      endif
c
c  Work array (ipoint) for determining storage locations of points.
c
      numr = max(1,(njklmax+1)/(machrsz/machisz))
      call ema(pipoint,numr,kipoint,'ipoint',11)
c
c  Read each block and store face points (xyz) and quad pointers (ndf).
c
      nn = 0
      nc = 0
      nq = 0
      iaf(1) = 1

      do n=1,nb

        njkl = nj(n)*nk(n)*nl(n)
        if (formatted) then
          read (7, *) ((xyz8(m,i),i=1,njkl),m=1,3)
          ihand(n) = ihand8(xyz8,nj(n),nk(n),nl(n),n)
        else
          if(nrpgrid.eq.4) then
            read(7) ((xyz4(m,i),i=1,njkl),m=1,3)
            ihand(n) = ihand4(xyz4,nj(n),nk(n),nl(n),n)
          else
            read(7) ((xyz8(m,i),i=1,njkl),m=1,3)
            ihand(n) = ihand8(xyz8,nj(n),nk(n),nl(n),n)
          endif
        end if
c
c  Store the face points (xyz) and pointers to points (ipoint).
c
        do l=1,nl(n)
          do k=1,nk(n)
            do j=1,nj(n)
              if(j.eq.1.or.j.eq.nj(n).or.
     .           k.eq.1.or.k.eq.nk(n).or.
     .           l.eq.1.or.l.eq.nl(n)) then
                nn = nn + 1
                ndx = indx(j,k,l,n)
                if(nrpgrid.eq.4) then
                  do m=1,3
                    xyz(m,nn) = xyz4(m,ndx)
                  end do
                else
                  do m=1,3
                    xyz(m,nn) = xyz8(m,ndx)
                  end do
                endif
                ipoint(ndx) = nn
              endif
            end do
          end do
        end do
c
c  Store the faces (ndf) using the pointers to points (ipoint) in the
c  face numbering order:  jmin,jmax,kmin,kmax,lmin,lmax.
c
        do j=1,nj(n),nj(n)-1
          do l=1,nl(n)-1
            do k=1,nk(n)-1
              nq = nq + 1
              ndf(1,nq) = ipoint(indx(j,k  ,l  ,n))
              ndf(2,nq) = ipoint(indx(j,k+1,l  ,n))
              ndf(3,nq) = ipoint(indx(j,k+1,l+1,n))
              ndf(4,nq) = ipoint(indx(j,k  ,l+1,n))
            end do
          end do
          nc = nc + 1
          iaf(nc+1) = nq + 1
          ndc(1,nc) = ipoint(indx(j,1    ,1    ,n))
          ndc(2,nc) = ipoint(indx(j,nk(n),1    ,n))
          ndc(3,nc) = ipoint(indx(j,nk(n),nl(n),n))
          ndc(4,nc) = ipoint(indx(j,1    ,nl(n),n))
        end do
c
        do k=1,nk(n),nk(n)-1
          do l=1,nl(n)-1
            do j=1,nj(n)-1
              nq = nq + 1
              ndf(1,nq) = ipoint(indx(j  ,k,l  ,n))
              ndf(2,nq) = ipoint(indx(j+1,k,l  ,n))
              ndf(3,nq) = ipoint(indx(j+1,k,l+1,n))
              ndf(4,nq) = ipoint(indx(j  ,k,l+1,n))
            end do
          end do
          nc = nc + 1
          iaf(nc+1) = nq + 1
          ndc(1,nc) = ipoint(indx(1    ,k,1    ,n))
          ndc(2,nc) = ipoint(indx(nj(n),k,1    ,n))
          ndc(3,nc) = ipoint(indx(nj(n),k,nl(n),n))
          ndc(4,nc) = ipoint(indx(1    ,k,nl(n),n))
        end do
c
        do l=1,nl(n),nl(n)-1
          do k=1,nk(n)-1
            do j=1,nj(n)-1
              nq = nq + 1
              ndf(1,nq) = ipoint(indx(j  ,k  ,l,n))
              ndf(2,nq) = ipoint(indx(j+1,k  ,l,n))
              ndf(3,nq) = ipoint(indx(j+1,k+1,l,n))
              ndf(4,nq) = ipoint(indx(j  ,k+1,l,n))
            end do
          end do
          nc = nc + 1
          iaf(nc+1) = nq + 1
          ndc(1,nc) = ipoint(indx(1    ,1    ,l,n))
          ndc(2,nc) = ipoint(indx(nj(n),1    ,l,n))
          ndc(3,nc) = ipoint(indx(nj(n),nk(n),l,n))
          ndc(4,nc) = ipoint(indx(1    ,nk(n),l,n))
        end do

      end do ! Next block

      close (7)
c
c  A message about memory usage.
c
      call memuser(numhu,numru,0)
      numhumax = numhu
      numrumax = numru
      call memsizr(numhs,numrs,0)
      nbytesu = numru*machrsz
      nbytess = numrs*machrsz
      pctused = 100.0 * float(nbytesu) / float(nbytesu+nbytess)
      write(string,1020) pctused
 1020 format(f6.2)
c
c  Check a few things..
c
      if(nn.ne.nnode) stop 'nn.ne.nnode'
      if(nc.ne.ncomp) stop 'nc.ne.ncomp'
      if(nq.ne.nquad) stop 'nq.ne.nquad'
      call emf(kipoint,'ipoint',12)
      if(nrpgrid.eq.4) then
        call emf(kpxyz4,'pxyz4',13)
      else
        call emf(kpxyz8,'pxyz8',14)
      endif
c
c  Allocate storage for xyz minmax extrema (xtrema) for each component,
c  then run through the components and find xyz extrema.  Elements
c  xtrema(1,n) to xtrema(6,n) are xmin,xmax,ymin,ymax,zmin,zmax for the
c  nth component.  The block number of the nth component is ((n-1)/6)+1.
c
      numr = max(1,(6*ncomp+1)/(machrsz/machrsz))
      call ema(pxtrema,numr,kxtrema,'xtrema',17)
      do n=1,ncomp
        xmin = xyz(1,ndf(1,iaf(n)))
        xmax = xmin
        ymin = xyz(2,ndf(1,iaf(n)))
        ymax = ymin
        zmin = xyz(3,ndf(1,iaf(n)))
        zmax = zmin
        do i=iaf(n),iaf(n+1)-1
          do j=1,4
            xmin = min(xmin,xyz(1,ndf(j,i)))
            xmax = max(xmax,xyz(1,ndf(j,i)))
            ymin = min(ymin,xyz(2,ndf(j,i)))
            ymax = max(ymax,xyz(2,ndf(j,i)))
            zmin = min(zmin,xyz(3,ndf(j,i)))
            zmax = max(zmax,xyz(3,ndf(j,i)))
          end do
        end do
        xtrema(1,n) = xmin
        xtrema(2,n) = xmax
        xtrema(3,n) = ymin
        xtrema(4,n) = ymax
        xtrema(5,n) = zmin
        xtrema(6,n) = zmax
      end do
c
c  For each block figure out which other blocks it touches.  Create
c  a connection file that looks like this (more or less):
c
c 72
c    1  0 0 0 0  13 1 2 3  61-1-2 3   2 1 2 3   0 0 0 0   4 1 2 3   1    1
c    2  0 0 0 0  14 1 2 3   1 1 2 3   3 1 2 3   0 0 0 0   5 1 2 3   1    1
c    3  0 0 0 0  15 1 2 3   2 1 2 3   0 0 0 0   0 0 0 0   6 1 2 3   1    1
c    4  0 0 0 0  16 1 2 3  64-1-2 3   5 1 2 3   1 1 2 3   7 1 2 3   1    1
c[...]
c   71 59 1 2 3  -5 1 2 3  70 1 2 3  72 1 2 3  68 1 2 3  -5 1 2 3   1    1
c   72 60 1 2 3  -5 1 2 3  71 1 2 3  -4 1 2 3  69 1 2 3  -5 1 2 3   1    1
cWAKES
c   0
c
c  The form of the connection file is first an integer number of blocks,
c  nb, and then nb lines consisting of six fields of four integers, then
c  two more integers for block handedness and undefined.  Each of the nb
c  lines gives block abutment information for the nth block.  It may be
c  assumed that the blocks are numbered in ascending order.  Each of the
c  six fields (one per hex face) lists which block it abuts and what the
c  relative basis orientations are.
c
c  If the block number is -4 it is a z = 0 symmetry plane.
c  [Later: use -4, -8, -9 to signal x = 0, y = 0, and z = 0 symmetry planes,
c  respectively for DPLR purposes.  These choices avoid other flags used by
c  FLO107MB, namely -5 for outer boundary, -2 for an Euler wall, and -7 for
c  a Navier-Stokes wall.]
c
c  If the block number is 0 then the face may need an outer or body boundary
c  condition (typically another negative integer) keyed in by the user with
c  a text editor.  Sometimes block abutments will not be found because of a
c  gap or other misalignment, in which case the user will have to fix the
c  mesh or assume responsibility for connecting blocks through "hyperspace".
c  The connection file ends with an empty "wakes" section.  The undefined
c  value at the end of the nth block's line may be assigned a verbose quality
c  number if the -v option is selected.  It will be the fraction of faces that
c  actually match on four corners within epsilon.

      n = lenstr(outconfn)
      write(*,'(3a)') ' Opening file ', outconfn(1:n), ' for write ...'
      open(8,file=outconfn(1:lenstr(outconfn)),form='FORMATTED')
      write(8,1040) nb
 1040 format(i3)
      nguesses = 0
      ninta = 0    ! DPLR's count of interior face matches
      allocate (iconsave(4,6,nb))
      symtol = min (10.* epsilon, 0.001) ! Avoid a nearly flat face at x ~ 0.

      do n=1,nb
        nquads = 0
        nhitquad = 0
        do 20 nf=1,6
          nc = nf+6*(n-1)
          do i=1,4
            icon(i,nf) = 0
          end do
          nmatches = 0
          do 10 m=1,nb
            if(m.eq.n) goto 10
            do mf=1,6
              mc = mf+6*(m-1)
              if(matches(xtrema,nc,mc,epsilon,nj,nk,nl,nb)) then
                nmatches = nmatches + 1
                ninta    = ninta    + 1
                if(nmatches.eq.1) then
                  call geticon(icon(1,nf),nc,mc,epsilon,xyz,ndf,ndc,
     .                         ihand,nguesses,strings,imatch,MSTRINGS)
                  if(verbose) then
                    call qhitter(icon(1,nf),nc,mc,epsilon,xyz,ndf,iaf,
     .                           nj,nk,nl,nhitquad,nquads,
     .                           strings(imatch))
                  endif
                else
                  write(*,1045) nmatches,n,nf,m,mf
 1045             format(' WARNING:  Repeat match number',i2,' found f',
     .                   'or block',i3,' face',i2,' at block',i3,' fac',
     .                   'e',i2,'.')
                endif
              endif
            end do
   10     continue

c  Check symmetry conditions for faces that do not appear to abut to
c  any other.  The tolerance need not be as tight as for matching faces.

          if (icon(1,nf).eq.0) then
            if (abs(xtrema(1,nc)).lt.symtol .and.
     .          abs(xtrema(2,nc)).lt.symtol) then
              icon(1,nf) = -4
              icon(2,nf) = 1
              icon(3,nf) = 2
              icon(4,nf) = 3
            else if (abs(xtrema(3,nc)).lt.symtol .and.
     .               abs(xtrema(4,nc)).lt.symtol) then
              icon(1,nf) = -8
              icon(2,nf) = 1
              icon(3,nf) = 2
              icon(4,nf) = 3
            else if (abs(xtrema(5,nc)).lt.symtol .and.
     .               abs(xtrema(6,nc)).lt.symtol) then
              icon(1,nf) = -9
              icon(2,nf) = 1
              icon(3,nf) = 2
              icon(4,nf) = 3
            endif
          endif

   20   continue ! Next face
c
c  hitfrac is a measure of fidelity, a figure of merit, a metric
c  that reports what fraction of the number of points in the face
c  actually come within epsilon of the points they're supposed to.
c
        if(nquads.eq.0) then
          hitfrac = 1.0
        else
          hitfrac = float(nhitquad)/float(nquads)
        endif
c
        if(abs(hitfrac-1.0).le.1.e-6) then
          write(8,1050) n,((icon(i,j),i=1,4),j=1,6),ihand(n)
 1050     format(i5,i3,3i2,5(1x,i3,3i2),2x,i2,3x,' 1')
        else
          write(8,1051) n,((icon(i,j),i=1,4),j=1,6),ihand(n),hitfrac
 1051     format(i5,i3,3i2,5(1x,i3,3i2),2x,i2,3x,f8.5)
        endif

        iconsave(:,:,n) = icon(:,:) ! For GASP & DPLR purposes

        call flush(8)

      end do ! Next block

      close (8)

      if(nguesses.gt.0) then
        write(*,1062) nguesses
 1062   format(' WARNING:  Had to guess face orientations',i5,' times.')
      endif

c     GASP- and DPLR-specific outputs for a single layer of blocks:

      write(*,'(a)')
     .  ' Opening GASP and DPLR control files for write ...'

      open (8, file='dplr.inputs',   status='unknown')
ccc   open (10,file='dplr.inputs.2', status='unknown') ! Coarser ni/nj/nk
      open (9, file='gasp.inputs',   status='unknown')

c     Option to add header and trailer to DPLR control files:

      open (11,file='sample.inputs', status='old', iostat=ios)
      header_trailer = ios == 0

      write (9, '(a)') 'GASP Input Deck:'
      write (9, 1063)
 1063 format ('-------------------------------------------------------')
      write (9, '(23x, a)') 'FILE INFO'
      write (9, 1063)
      write (9, '(a)') 'chemModPath    gridDir    solutionDir    bcDir',
     .  '''./db1.bin''     ''.''    ''.''    ''.''',
     .  'gridFileMode  solnFileMode  bcFilemode',
     .  '0             0             0'
      write (9, 1063)
      write (9, '(21x, a)') 'GENERAL INFO'
      write (9, 1063)
      write (9, '(a)')
     .   'iunits  rhoNDim   vNDim   tNDim   lNDim  pGauge',
     .   '  1      1.e+00  1.e+00  1.e+00  1.e+00  0.e+00',
     .   'irest  initByZone  memoryMode  residMode  cpuMode',
     .   '  1      0           2           1          0',
     .   'nZone nZonalBounds nPartStyle nPhysMod nBlock iBlkStt iBlkEnd'
      write (9, '(i4, i9, a)')
     .   nb, ninta/2, '       1         1        4      2       4'
      write (9, 1063)
      write (9, '(22x, a)') 'ZONE INFO'

      if (header_trailer) then ! Transcribe the sample DPLR file header

        nheader = 0
        do ! Until the start of the block-specific controls
          nheader = nheader + 1
          read (11, '(a)') line
          n = len_trim (line)
!!!       if (index (line(1:n), '==') > 0) exit  ! Fails for DPLR V3.05
          if (n > 0) then
            call upcase (line(1:n))
            if (index (line(1:n), 'BLOCK') > 0) then
              nheader = nheader - 2
              exit
            end if
          end if
        end do

        rewind (11)

        do l = 1, nheader
          read (11, '(a)') line
          n = len_trim (line)
          write (8, '(a)') line (1:n)
        end do

      end if

      aiformat = '(a, in)'

      do n = 1, nb

        write (8, '(a)') '=========================================='
        write (9, 1063)

        if (n < 10)  then
          aiformat(6:6) = '2'
        else if (n < 100)  then
          aiformat(6:6) = '3'
        else if (n < 1000)  then
          aiformat(6:6) = '4'
        else
          aiformat(6:6) = '5'
        end if

        write (8, aiformat) 'Block #', n
        write (9, aiformat) 'ZONE #', n

c       Most of the DPLR format:

        write (8, '(a)') '=========================================='
        write (8, '(/, a)') ' iconr isim ifree initi ibrad ibadpt isurf'
        write (8, '(7i6)') -1, 1, 1, 1, 1, 1, 1
        write (8, '(/, a)') '  iflx  iord  ilim idiss  epsi'
        write (8, '(a)')    '     4     3     1     1   0.3'
        write (8, '(/, a)') '  jflx  jord  jlim jdiss  epsj'
        write (8, '(a)')    '     4     3     1     1   0.3'
        write (8, '(/, a)') '  kflx  kord  klim kdiss  epsk'
        write (8, '(a)')    '     4     3     1     0   0.3'
        write (8, '(/, a)')
     .    'iextst  nrlx ildir  ibcu iblag   ilt ibdir   cflm'
        write (8, '(a)')
     .    '    -1     4     0     1    -1    -1     1  1.d20'
        write (8, '(/, a)') ' Boundary condition type [ibc]:'
        write (8,    '(a)') ' imin imax jmin jmax kmin kmax'

        do nf = 1, 6
          i = iconsave(1,nf,n)
          if      (i == -4) then
            ibc(nf) = 17       ! u = -u (x = 0 symmetry plane)
          else if (i == -8) then
            ibc(nf) = 18       ! v = -v (y = 0    "       "  )
          else if (i == -9) then
            ibc(nf) = 19       ! w = -w (z = 0    "       "  )
          else if (i == -5) then ! ? Where is -5 ever assigned?
            ibc(nf) = 3        ! Supersonic exit, first order extrapolation
          else if (i > 0) then
            ibc(nf) = 20       ! Zonal boundary in the volume
          else                 ! Educated guess for one-layer Shuttle-type grids
            if (nf == 5) then
              ibc(5) = 26      ! Catalytic radiative equilibrium wall (k = 1)
            else if (nf == 6) then
              ibc(6) = 1       ! Freestream (k = nk)
            else
              ibc(nf) = 2      ! Specified inflow, or supersonic outflow
!!!           ibc(nf) = 3      ! Supersonic exit, first order extrapolation
            end if
          end if
        end do

!       Look for known solid wall BC corrections from 'template.inp.2':

        if (cavity) then
          do i = 1, niblist
            if (iblist(i) == n) then
              do nf = 1, 6
                if (ibc(nf) == 1 .or. ibc(nf) == 2) ibc(nf) = 26
              end do
            end if
          end do
        else if (plug) then
          do i = 1, niblist
            if (iblist(i) == n) ibc(3) = 26  ! jmin face
          end do
        end if

        write (8, '(6i5, /)') ibc

c       Rest of the GASP BC format:

        do nf = 1, 6
          if (ibc(nf) == 26) then
            ibc(nf) = 35
          else
            ibc(nf) = -ibc(nf)
          end if
        end do

        write (9, '(a)')
     .    'Initial Conditions',
     .    '   initPhysMod',
     .    '       1',
     .    'Boundary Conditions',
     .    'nSurfTypes    i0   id   j0   jd   k0   kd'
        write (9, '(a, 6i5)') '     6     ', ibc
        write (9, '(a)')
     .    '   twall   pback',
     .    '  288.00   101325',
     .    'Mesh Sequencing',
     .    '   nSeqLevel',
     .    '       1',
     .    '    seqNum   initSolFile   intrpQ   ilev   jlev   klev',
     .    '       1         0           0       1      1      1'

      end do ! Next block

      if (header_trailer) then ! Add the rest of the sample DPLR inputs

        write (8, '(a)') ' ========================================='

        do ! Until 'Freestream' string
          read (11, '(a)') line
          n = len_trim (line)
          if (n > 0) then
            if (index (line(1:n), 'Free') > 0) exit
          end if
        end do

        do ! Until EOF
          write (8, '(a)') line(1:n)
          read (11, '(a)', iostat=ios) line
          if (ios /= 0) exit
          n = len_trim (line)
        end do

        close (11)

      end if

      close (8)

      write (9, 1063)
      write (9, '(18x, a)') 'ZONAL BOUNDARY INFO'

      write(*,'(a)') ' Opening file dplr.interfaces for write ...'

      open (8, file='dplr.interfaces', status='unknown')

      write (8, '(a)') ' ZONAL BOUNDARY INFORMATION',
     .                 ' Cell matching - without dummy cells'
      write (8, 1065)
 1065 format (' -----------------------------------------')

      write (8, '(a)') '   zvers izdum',
     .                 '   3.00    0'
      write (8, '(/, a)') ' nblk ninta nintc'
      write (8, '(3i5)') nb, ninta/2, 0

c     List the interior block-to-block boundary pairs by scanning the saved
c     connectivity info for positive block numbers.  Avoid listing them twice.

      nzb = 0

      do n = 1, nb - 1 ! The last block can't match a face of itself

        do nf = 1, 6

          m = iconsave(1,nf,n)

          if (m > 0) then ! m is the block number of a matching face

            izone1(1) = n
            izone2(1) = m
            izone1(2) = nf

            do mf = 1, 6 ! Locate the matching face number of block m
               if (iconsave(1,mf,m) == n) then
                 iconsave(1,mf,m) = -n ! Suppress it from remaining searches 
                 izone2(2) = mf
                 exit
               end if
            end do

            select case (nf) ! Block n info tells us all the rest we need

            case (1, 2)            ! An i face of block n

              izone1(3) = 2        ! j direction
              j = iconsave(3,nf,n) ! Retain cyclic order in block n info
              izone2(3) = abs (j)
              numj = nk(n) - 1     ! Scott's use of j,k,l is unfortunate
              if (j > 0) then
                 izone1(4) = 1;  izone1(5) = numj
              else
                 izone1(4) = numj;  izone1(5) = 1
              end if
              izone2(4) = 1;  izone2(5) = numj

              izone1(6) = 3        ! k direction
              k = iconsave(4,nf,n)
              izone2(6) = abs (k)
              numk = nl(n) - 1
              if (k > 0) then
                 izone1(7) = 1;  izone1(8) = numk
              else
                 izone1(7) = numk;  izone1(8) = 1
              end if
              izone2(7) = 1;  izone2(8) = numk

            case (3, 4)            ! A j face of block n

              izone1(3) = 3        ! k direction
              k = iconsave(4,nf,n) ! Retain cyclic order
              izone2(3) = abs (k) 
              numk = nl(n) - 1
              if (k > 0) then
                 izone1(4) = 1;  izone1(5) = numk
              else
                 izone1(4) = numk;  izone1(5) = 1
              end if
              izone2(4) = 1;  izone2(5) = numk

              izone1(6) = 1        ! i direction
              i = iconsave(2,nf,n)
              izone2(6) = abs (i)
              numi = nj(n) - 1
              if (i > 0) then
                 izone1(7) = 1;  izone1(8) = numi
              else
                 izone1(7) = numi;  izone1(8) = 1
              end if
              izone2(7) = 1;  izone2(8) = numi

            case (5, 6)            ! A  k face of block n

              izone1(3) = 1        ! i direction
              i = iconsave(2,nf,n)
              izone2(3) = abs (i)
              numi = nj(n) - 1
              if (i > 0) then
                 izone1(4) = 1;  izone1(5) = numi
              else
                 izone1(4) = numi;  izone1(5) = 1
              end if
              izone2(4) = 1;  izone2(5) = numi

              izone1(6) = 2        ! j direction
              j = iconsave(3,nf,n)
              izone2(6) = abs (j)
              numj = nk(n) - 1
              if (j > 0) then
                 izone1(7) = 1;  izone1(8) = numj
              else
                 izone1(7) = numj;  izone1(8) = 1
              end if
              izone2(7) = 1;  izone2(8) = numj

            end select

            write (8, '(a)')
            write (8, 1065)
            write (9, 1063)

            nzb = nzb + 1
            if (nzb < 10) then
              write (8, '(a, i2)') ' Zonal Boundary #', nzb
              write (9, '(a, i2)') ' Zonal Boundary #', nzb
            else if (nzb < 100) then
              write (8, '(a, i3)') ' Zonal Boundary #', nzb
              write (9, '(a, i3)') ' Zonal Boundary #', nzb
            else
              write (8, '(a, i4)') ' Zonal Boundary #', nzb
              write (9, '(a, i4)') ' Zonal Boundary #', nzb
            end if

            write (9, '(a)')
     .        'izbpass  zbType  zbFluxCrct',
     .        '  0        0       0'

            do l = 8, 9
              write (l, '(a)')
     .          '   nz nface ndr1 nst1 nen1 ndr2 nst2 nen2'
            end do

            write (8, '(8i5)') izone1, izone2

            write (9, '(8i5)') izone1
            write (9, '(a)')
     .         'neqn    zbmap(1:neqn)',
     .         '  9     1  2  3  4  5  6  7  8  9',
     .         '   nz nface ndr1 nst1 nen1 ndr2 nst2 nen2'
            write (9, '(8i5)') izone2
            write (9, '(a)')
     .         'neqn    zbmap(1:neqn)',
     .         '  9     1  2  3  4  5  6  7  8  9'

          end if

        end do ! Next face of block n

      end do ! Next block n

      close (8)

      call gasp_tail_end (nb)    ! Add the rest of the GASP control data

      close (9)

c  Print out some peak memory stats.

      call memuser(numhu,numru,0)
      numhumax = max(numhumax,numhu)
      numrumax = max(numrumax,numru)
      nbytesumax = numrumax*machrsz
      pctused = 100.0 * float(nbytesumax) / float(nbytesu+nbytess)
      write(string,1020) pctused
      if(string(1:1).ne.' ') then
        ib = 1
      else
        ib = 2
      endif
      write(*,1070) numhumax,nbytesumax,
     .              string(ib:lenstr(string))
 1070 format(/,' There were max',i3,' arrays allocated using max',i10,
     .       ' bytes mem (',a,'%).')

      write(*,*) 'TEMPLATE:  Normal Termination.'

      end program template

c-----------------------------------------------------------------------

      subroutine gasp_tail_end (nb)

c     Add the rest of the GASP control data for the single-sequence case.

      implicit real*8 (a-h,o-z)

      integer ncycles (8)
      data ncycles          /250, 500, 250, 1000,
     .                       100, 500, 100, 1300/
      integer itmstep (8)
      data itmstep          /1, 1, 2, 2,
     .                       1, 1, 2, 2/
      character dt (8) * 24
      data dt               /'   -1.00e-02   -8.00e-01',
     .                       '   -8.00e-01   -8.00e-01',
     .                       '   -1.00e-02   -4.00e-01',
     .                       '   -4.00e-01   -4.00e-01',
     .                       '   -1.00e-02   -8.00e-01',
     .                       '   -8.00e-01   -8.00e-01',
     .                       '   -1.00e-02   -8.00e-01',
     .                       '   -8.00e-01   -8.00e-01'/

      write (9, 1063)
      write (9, '(20x, a)') 'PARTITION STYLES'
      write (9, 1063)
      write (9, '(a)')
     .  'PARTITION STYLE # 1',
     .  'nDir',
     .  ' 1',
     .  'dir  numPart',
     .  ' 1      1'
      write (9, 1064)
      write (9, '(17x, a)') 'PHYSICAL MODELING INFO'
      write (9, 1064)
      write (9, '(a)')
     .  'PHYSICAL MODEL # 1'
      write (9, 1063)
      write (9, '(15x, a)') 'CHEMISTRY & THERMODYNAMICS'
      write (9, 1063)
      write (9, '(a)')
     .  'nspec   nnev   prefDiss  vibRelax  itherm  chemmod  ieq',
     .  '  5      0        0         0        1    ''Park 1''    3'
      write (9, 1063)
      write (9, '(20x, a)') 'INVISCID FLUXES'
      write (9, 1063)
      write (9, '(a)')
     .  'imarch',
     .  '   0',
     .  'invflxi  invflxj  invflxk',
     .  '   2        2        2',
     .  '  rkapi    rkapj    rkapk    sdm2    sdm4',
     .  ' 0.3333   0.3333   0.3333  0.0000  0.0000',
     .  '   limi     limj     limk  rk_ven',
     .  '   2        2        2      1.000'
      write (9, 1063)
      write (9, '(21x, a)') 'VISCOUS FLUXES'
      write (9, 1063)
      write (9, '(a)')
     .  'isViscous',
     .  '   1',
     .  'visflxi  visflxj  visflxk',
     .  '   1        1        1',
     .  ' modlmu    modlk   imodld   ivac',
     .  '   1        1        1       2',
     .  '  prl      prt      scl     sct',
     .  ' 0.72     0.90     0.70    0.50',
     .  ' ikeps    kemin   wallfunc  igb',
     .  '   0        0        0       0'
      write (9, 1065)
      write (9, '(20x, a)') 'INITIAL CONDITIONS'
      write (9, 1065)
      write (9, '(a)')
     .  'icond',
     .  '   2',
     .  'Vel/Mach       cx    cy        cz  temp/press turbi tkelref',
     .  '17.88    0.776952  0.00  0.629560    238.48    0.05   0.001',
     .  'rho_spec',
     .  '2.00206e-04',
     .  '6.07846e-05',
     .  '0',
     .  '0',
     .  '0'
      write (9, 1065)
      write (9, '(18x, a)') 'REFERENCE STATE FOR BCs'
      write (9, 1065)
      write (9, '(a)')
     .  'icond',
     .  '   2',
     .  'Vel/Mach       cx    cy        cz  temp/press turbi tkelref',
     .  '17.88    0.776952  0.00  0.629560    238.48    0.05   0.001',
     .  'rho_spec',
     .  '2.00206e-04',
     .  '6.07846e-05',
     .  '0',
     .  '0',
     .  '0'
      write (9, 1065)
      write (9, '(22x, a)') 'BLOCK INFO'

      do n = 5, 8
        write (9, 1065)
        write (9, '(a, i2)')
     .    'BLOCK #', n
        write (9, '(a)')
     .    '  imcont   nswp   ncycle   nwres   mstage   rtolr     rtola'
        write (9, '(a, i5, a)')
     .    '     1        1    ', ncycles(n),
     .                            '      50     1   1.00e-05  1.00e-12'
        write (9, '(a)')
     .    '   mgstyle  nitfine  nitcrct  nitsmth     smthfac     dtfac',
     .    '     0        1        1       0         3.00e-01  1.00e+00',
     .    '     SWEEP # 1',
     .    '     minZone maxZone  iseq npm iswpdir swpAll iplstt iplend'
        write (9, '(a, i4, a)')
     .    '     1      ', nb,
     .                    '       1    1    1      1      1      16'
        write (9, '(a)')
     .    '     impl itmstep keslv ichemslv nplane  inner  mxin  tolin'
        write (9, '(a, i8, a)')
     .    '     2', itmstep(n),
     .                   '    1     1         1       0      5     0.01'
        write (9, '(a)')
     .    '     dtmin         dtmax  irelu   nremax    tolreu'
        write (9, '(a24, a)')
     .                     dt(n), '      0     0     1.00e-02'
      end do

 1063 format ('-------------------------------------------------------')
 1064 format ('*******************************************************')
 1065 format
     .   ('-----------------------------------------------------------')

      end subroutine gasp_tail_end

c-----------------------------------------------------------------------

      subroutine ema(px,numr,kx,serr,nerr)

      dimension  x(1)
      pointer  (px,x)
      character(*) serr*(*)
c
c  error memory allocation, an error handling wrapper for memallr.
c
c  Scott D. Thomas, Sterling Software, Inc.
c  Contract NAS 2-13210/6/5, NASA Ames Research Center.
c
      if(memptrr(px,memallr(numr,kx)).le.0) then
        write(*,1010) serr,nerr
 1010   format(' ERROR ema ',a,' nerr=',i4)
      endif
c
      return
      end

c-----------------------------------------------------------------------

      subroutine emf(kx,serr,nerr)

      character(*) serr*(*)
c
c  error memory free, an error handling wrapper for memfrer.
c
c  Scott D. Thomas, Sterling Software, Inc.
c  Contract NAS 2-13210/6/5, NASA Ames Research Center.
c
      if(memfrer(kx).lt.0) then
        write(*,1010) serr,nerr
 1010   format(' ERROR emf ',a,' nerr=',i10)
      endif
c
      return
      end

c-----------------------------------------------------------------------

      subroutine emp(px,kx,serr,nerr)

      dimension  x(1)
      pointer  (px,x)
      character(*) serr*(*)
c
c  error memory pointer, an error handling wrapper for memptrr.
c
c  Scott D. Thomas, Sterling Software, Inc.
c  Contract NAS 2-13210/6/5, NASA Ames Research Center.
c
      if(memptrr(px,kx).le.0) then
        write(*,1010) serr,nerr
 1010   format(' ERROR emp ',a,' nerr=',i4)
      endif
c
      return
      end

c-----------------------------------------------------------------------

      subroutine emr(numr,kx,serr,nerr)

      character(*) serr*(*)
c
c  error memory reallocation, an error handling wrapper for memrear.
c
c  Scott D. Thomas, Sterling Software, Inc.
c  Contract NAS 2-13210/6/5, NASA Ames Research Center.
c
      if(memrear(numr,kx).lt.0) then
        write(*,1010) serr,nerr
 1010   format(' ERROR emr ',a,' nerr=',i4)
        stop
      endif
c
      return
      end

c-----------------------------------------------------------------------

      function memallr(nreal,handle)

      implicit real*8 (a-h,o-z)
      integer  handle
      parameter (MXNUMARR=4000,MXNUMREAL=100000000)
      parameter (MXALLOCR=MXNUMREAL+2*MXNUMARR+1)
      common /meminte/ narrayr,iarrayr(MXNUMARR)
      common /memreal/ realmem(MXALLOCR)
c
c  Allocate memory for nreal real*8's from the realmem array.  The
c  value returned in handle is the integer array handle such that the
c  allocated real*8 array starts at realmem(n) where
c  n=iarrayr(handle) and ends at realmem(n+nreal-1).  When freeing
c  memory that is allocated here in memallr, use handle as the argument
c  to the memfrer routine.  The function return value memallr=handle.
c
c  The value of nreal is also stored in realmem at realmem(n-1) and at
c  realmem(n+nreal).  By convention, the value stored at
c  realmem(n+nreal+1) is zero for the last array or else it indicates
c  the size of the next array.
c
c  This works only if enough memory is available.  There has to be
c  another handle left as well as enough memory at the top of realmem.
c
c  If 0 integers are requested then no action will be taken and the
c  value of handle will be set to 0.  Values .lt. 0 indicate failure.
c
c  Scott D. Thomas, Sterling Software, Inc.
c  Contract NAS 2-13210/6/5, NASA Ames Research Center.
c
      if(nreal.eq.0) then
        handle = 0
        goto 20
      endif
      if(nreal.lt.0) then
        write(*,1010) 'may not allocate negative number:',nreal
        handle = -1
        goto 20
      endif
      if(narrayr.ge.MXNUMARR) then
        write(*,1010) 'no more array handles, narrayr:',narrayr,MXNUMARR
        handle = -2
        goto 20
      endif
      next = memchkr()
      if(next.le.0) then
        write(*,1010) 'memchkr returned bad value:',next
        handle = -3
        goto 20
      endif
      memsize = max(0,MXALLOCR-(2*(MXNUMARR-narrayr)+1)-next+1)
      if(nreal.gt.memsize) then
        write(*,1010) 'not enough memory:  nreal,memsize:',
     .                nreal,memsize
        handle = -4
        goto 20
      endif
      realmem(next)         = nreal
      realmem(next+nreal+1) = nreal
      realmem(next+nreal+2) = 0
      narrayr               = narrayr + 1
      if(narrayr.eq.1) then
        nexthand = 1
      else
c
c  Run down to the next available array handle location.
c
        numbhand = 0
        do i=1,MXNUMARR
          if(iarrayr(i).eq.0) then
            nexthand = i
            goto 10
          else
            numbhand = numbhand + 1
          endif
          if(numbhand.eq.narrayr-1) then
            nexthand = i+1
            goto 10
          endif
        end do
   10   continue
      endif
      if(nexthand.gt.MXNUMARR) then
        write(*,1010) 'nexthand.gt.MXNUMARR:',nexthand,MXNUMARR
        handle = -5
        goto 20
      endif
      iarrayr(nexthand) = next+1
      handle            = nexthand
   20 continue
      memallr = handle
 1010 format(' MEMALLR:  ERROR, ',a,4i10)
c
      return
      end

c-----------------------------------------------------------------------

      function memaszr(handle)

      implicit real*8 (a-h,o-z)
      integer handle,address
      parameter (MXNUMARR=4000,MXNUMREAL=100000000)
      parameter (MXALLOCR=MXNUMREAL+2*MXNUMARR+1)
      common /meminte/ narrayr,iarrayr(MXNUMARR)
      common /memreal/ realmem(MXALLOCR)
c
c  Scott D. Thomas, Sterling Software, Inc.
c  Contract NAS2-13210/6/5, NASA Ames Research Center.
c  21-NOV-1996..4-DEC-1996.
c
c  Return the current size of the real array with handle value equal
c  to handle.  This array starts at realmem(iarrayr(handle)).
c  Negative return values indicate an error.
c
      if(handle.lt.1) then
        memaszr = -1
        goto 10
      endif
      if(handle.gt.MXNUMARR) then
        memaszr = -2
        goto 10
      endif
      address = iarrayr(handle)
      if(address.lt.2) then
        memaszr = -3
        goto 10
      endif
      if(address.gt.MXALLOCR-2) then
        memaszr = -4
        goto 10
      endif
      memaszr = realmem(address-1)
   10 continue
c
      return
      end

c-----------------------------------------------------------------------

      function memchkr()

      implicit real*8 (a-h,o-z)
      parameter (MXNUMARR=4000,MXNUMREAL=100000000)
      parameter (MXALLOCR=MXNUMREAL+2*MXNUMARR+1)
      common /meminte/ narrayr,iarrayr(MXNUMARR)
      common /memreal/ realmem(MXALLOCR)
c
c  Scott D. Thomas, Sterling Software, Inc.
c  Contract NAS2-13210/6/5, NASA Ames Research Center.
c  21-NOV-1996..4-DEC-1996.
c
c  Checks the real memory arrays, iarrayr and realmem, to make sure
c  all the record information is consistent.  Return values:
c
c  memchkr .le. 0 means something is wrong.
c  memchkr .gt. 0 is the address where the next memallr may start.
c
      if(MXNUMARR.lt.1) then
        write(*,1010) 'MXNUMARR.lt.1:',MXNUMARR
        memchkr = -1
        goto 30
      endif
      if(MXALLOCR.lt.1) then
        write(*,1010) 'MXALLOCR.lt.1:',MXALLOCR
        memchkr = -2
        goto 30
      endif
      if(narrayr.lt.0) then
        write(*,1010) 'narrayr.lt.0:',narrayr
        memchkr = -3
        goto 30
      endif
      if(narrayr.gt.MXNUMARR) then
        write(*,1010) 'narrayr.gt.MXNUMARR:',narrayr,MXNUMARR
        memchkr = -4
        goto 30
      endif
c
c  Check the iarrayr array to make sure that every valid pointer into
c  the realmem array corresponds to consistent storage of data.
c
      if(narrayr.gt.0) then
        i      = 0
        numarr = 0
    5   continue
          i = i + 1
          if(i.gt.MXNUMARR) then
            write(*,1010) 'i.gt.MXNUMARR:',i,MXNUMARR
            memchkr = -5
            goto 30
          endif
          irec = iarrayr(i)
          if(irec.eq.0) goto 5
          numarr = numarr + 1
          if(irec.lt.2) then
            write(*,1010) 'iarrayr(i).lt.2 at i=',i,iarrayr(i)
            memchkr = -6
            goto 30
          endif
          if(irec.gt.MXALLOCR-2) then
            write(*,1010) 'iarrayr(i).gt.MXALLOCR-2 at i=',i,iarrayr(i)
            memchkr = -7
            goto 30
          endif
          if(i.gt.1) then
            do j=1,i-1
              if(iarrayr(i).eq.iarrayr(j)) then
                write(*,1010) 'iarrayr(i).eq.iarrayr(j),i,j:',
     .                         iarrayr(i),   iarrayr(j),i,j
                memchkr = -8
                goto 30
              endif
            end do
          endif
          irecsz = realmem(irec-1)
          irsend = realmem(irec+irecsz)
          if(irecsz.ne.irsend) then
            write(*,1010) 'irecsz.ne.irsend:',irecsz,irsend
            memchkr = -9
            goto 30
          endif
        if(numarr.lt.narrayr) goto 5
      endif
c
c  Run through the realmem array to make sure all allocated arrays
c  are stored consistently.  The number of arrays found should be
c  equal to narrayr.
c
      numarr = 0
      next   = 1
   10 continue
        if(realmem(next).eq.0) goto 20
        irecsz = realmem(next)
        irsend = realmem(next+irecsz+1)
        if(irecsz.ne.irsend) then
          write(*,1010) 'irecsz.ne.irsend:',irecsz,irsend
          memchkr = -10
          goto 30
        endif
        itmp = realmem(next)
        next = next + itmp + 2
        if(next.lt.0) then
          write(*,1010) 'next.lt.0:',next
          memchkr = -11
          goto 30
        else if(next.gt.MXALLOCR) then
          write(*,1010) 'next.gt.MXALLOCR:',next,MXALLOCR
          memchkr = -12
          goto 30
        endif
        numarr = numarr + 1
        goto 10
   20 continue
      if(numarr.ne.narrayr) then
        write(*,1010) 'numarr.ne.narrayr:',numarr,narrayr
        memchkr = -13
        goto 30
      endif
      memchkr = next
   30 continue
 1010 format(' MEMCHKR:  ERROR, ',a,4i10)
c
      return
      end

c-----------------------------------------------------------------------

      function memfrer(handle)

      implicit real*8 (a-h,o-z)
      integer handle
      parameter (MXNUMARR=4000,MXNUMREAL=100000000)
      parameter (MXALLOCR=MXNUMREAL+2*MXNUMARR+1)
      common /meminte/ narrayr,iarrayr(MXNUMARR)
      common /memreal/ realmem(MXALLOCR)
c
c  Scott D. Thomas, Sterling Software, Inc.
c  Contract NAS2-13210/6/5, NASA Ames Research Center.
c  21-NOV-1996..4-DEC-1996.
c
c  Free the space reserved for the array starting at iarrayr(handle)
c  in the realmem array.  If handle comes in .le. 0 then do nothing.
c  Freeing memory involves shifting higher elements down to fill the
c  gap, if any, and adjusting the array handles.  The value of handle
c  will be changed to 0 if memory is freed correctly, but on error it
c  will be set to a negative value.  The function return value is the
c  same as the value of handlei.
c
      if(handle.le.0) goto 40
      if(handle.gt.MXNUMARR) then
        write(*,1010) 'handle.gt.MXNUMARR:',handle,MXNUMARR
        handle = -1
        goto 40
      endif
      nreal = realmem(iarrayr(handle)-1)
      nsize = realmem(iarrayr(handle)+nreal)
      if(nreal.ne.nsize) then
        write(*,1010) 'nreal.ne.nsize:',nreal,nsize
        handle = -2
        goto 40
      endif
      next = iarrayr(handle) + nreal + 1
      if(next.lt.1.or.next.gt.MXALLOCR) then
        write(*,1010) 'next out of range:  next,MXALLOCR:',
     .                next,MXALLOCR
        handle = -3
        goto 40
      endif
      if(realmem(next).eq.0.0) then
c
c  If you free the last allocated array there is little to do.
c
        realmem(iarrayr(handle)-1) = 0.0
      else
        last = memchkr()
        if(next.ge.last) then
          write(*,1010) 'next.ge.last:',next,last
          handle = -4
          goto 40
        endif
c
c  Otherwise, adjust the handles and ...
c
        irecsz = realmem(iarrayr(handle)-1)
        locnxt = iarrayr(handle)+irecsz+2
   10   continue
          irecsz = realmem(locnxt-1)
          if(irecsz.eq.0) goto 30
          do i=1,MXNUMARR
            if(iarrayr(i).eq.locnxt) then
              iarrayr(i) = locnxt-(nreal+2)
              goto 20
            endif
          end do
   20     continue
          locnxt = locnxt+irecsz+2
        goto 10
   30   continue
c
c  ... shift the array elements down to fill the gap.
c
        do i=next,last
          realmem(i-(nreal+2)) = realmem(i)
        end do
      endif
      iarrayr(handle) = 0
      narrayr         = narrayr - 1
      handle          = 0
   40 continue
      memfrer         = handle
 1010 format(' MEMFRER:  ERROR, ',a,2i10)
c
      return
      end

c-----------------------------------------------------------------------

      subroutine meminit

      implicit real*8 (a-h,o-z)
      parameter (MXNUMARR=4000,MXNUMREAL=100000000)
      parameter (MXALLOCR=MXNUMREAL+2*MXNUMARR+1)
      common /meminte/ narrayr,iarrayr(MXNUMARR)
      common /memreal/ realmem(MXALLOCR)
c
c  Scott D. Thomas, Sterling Software, Inc.
c  Contract NAS 2-13210/6/5, NASA Ames Research Center.
c
      narrayr    = 0
      iarrayr(1) = 0
c
      return
      end

c-----------------------------------------------------------------------

      function memptrr(raptr,rahandle)

      implicit real*8 (a-h,o-z)
      dimension ra(1)
      pointer  (raptr,ra)
      integer   rahandle
      parameter (MXNUMARR=4000,MXNUMREAL=100000000)
      parameter (MXALLOCR=MXNUMREAL+2*MXNUMARR+1)
      common /meminte/ narrayr,iarrayr(MXNUMARR)
      common /memreal/ realmem(MXALLOCR)
c
c  Scott D. Thomas, Sterling Software, Inc.
c  Contract NAS2-13210/6/5, NASA Ames Research Center.
c  21-NOV-1996..4-DEC-1996.
c
c  Set the pointer so that the array it points to will coincide with
c  realmem(iarrayr(rahandle).  This permits the variables and memory
c  arrays defined in memparm.h to be hidden from the caller's view.
c
c  If no error is detected then the return value will be rahandle.
c  If rahandle comes in as zero then the pointer will remain
c  uninitialized.  Negative return values indicate an error.
c
c  Example:
c
c     implicit real*8 (a-h,o-z)
c     dimension ra(1)
c     pointer  (raptr,ra)
c     integer   rahandle
c     nreal = 100
c     if(memptrr(raptr,memallr(nreal,rahandle)).lt.0) stop 'memptrr'
c     do i=1,nreal
c       ra(i) = 0.0
c     end do
c
c  Note that if nreal is zero, then no array will be allocated.
c  The value of rahandle will not be changed in this routine.
c
c  Caller must be careful to avoid using a pointer which may go stale
c  after memfrer or memrear routines are called.  To be safe, call
c  memptrr again to refresh the pointer.
c
      if(rahandle.le.0) goto 10
      last = memchkr()
      if(last.le.0) then
        write(*,1010) 'memchkr.lt.1:',last
        memptrr = -1
        goto 10
      endif
      if(rahandle.gt.MXNUMARR) then
        write(*,1010) 'rahandle.gt.MXNUMARR:',rahandle,MXNUMARR
        memptrr = -2
        goto 10
      endif
      m = iarrayr(rahandle)
      if(m.lt.2.or.m.gt.last-2.or.m.gt.MXALLOCR-2) then
        write(*,1010) 'iarrayr(rahandle) out of range:',rahandle,m
        memptrr = -3
        goto 10
      endif
      raptr = loc(realmem(m))
   10 continue
      memptrr = rahandle
 1010 format(' MEMPTRR:  ERROR, ',a,4i10)
c
      return
      end

c-----------------------------------------------------------------------

      function memrear(nreal,handle)

      implicit real*8 (a-h,o-z)
      integer  nreal,handle,address,curarsz,curarsz2,shift
      parameter (MXNUMARR=4000,MXNUMREAL=100000000)
      parameter (MXALLOCR=MXNUMREAL+2*MXNUMARR+1)
      common /meminte/ narrayr,iarrayr(MXNUMARR)
      common /memreal/ realmem(MXALLOCR)
c
c  Scott D. Thomas, Sterling Software, Inc.
c  Contract NAS2-13210/6/5, NASA Ames Research Center.
c  21-NOV-1996..4-DEC-1996.
c
c  The memallr routine must be called first to get handle.
c
c  Reallocate memory for the handle array which starts at realmem(n)
c  where n=iarrayr(handle).  That is, perform an in-place adjustment
c  of the size of the array from what it is now to the new value,
c  nreal, which comes in through the argument list.  The array may
c  shrink to as little as one element or grow to as many elements as
c  can be allocated.
c
c  This involves adjusting array handles and shifting arrays in memory
c  that are located after the current array in order to avoid gaps in
c  the realmem array.  The function return value will be handle
c  unless there is an error, in which case the value will be .lt. 0.
c  The value of handle will not be changed.
c
c  If nreal .lt. current array size, higher elements will be lost.
c  If nreal .eq. current array size, no action will be taken.
c  If nreal .gt. current array size, added elements will be junk.
c  In any case, the original array elements will be unchanged.
c
      if(handle.lt.1) then
        write(*,1010) 'handle is less than 1:',handle
        memrear = -1
        goto 70
      endif
      if(handle.gt.MXNUMARR) then
        write(*,1010) 'handle is .gt. MXNUMARR:',handle,MXNUMARR
        memrear = -2
        goto 70
      endif
      address = iarrayr(handle)
      if(address.lt.2) then
        write(*,1010) 'address is less than 2:',address
        memrear = -3
        goto 70
      endif
      if(address.gt.MXALLOCR-2) then
        write(*,1010) 'address is .gt. MXALLOCR-2:',address,MXALLOCR-2
        memrear = -4
        goto 70
      endif
      curarsz = realmem(address-1)
      if(curarsz.lt.1) then
        write(*,1010) 'curarsz is less than 1:',curarsz
        memrear = -5
        goto 70
      endif
      curarsz2 = realmem(address+curarsz)
      if(curarsz.ne.curarsz2) then
        write(*,1010) 'curarsz .ne. curarsz2:',curarsz,curarsz2
        memrear = -6
        goto 70
      endif
c
c  The handle is valid and the data is probably not corrupted.
c  If nreal .eq. current array size, no action will be taken.
c
      if(nreal.eq.curarsz) then
        memrear = handle
        goto 70
      endif
c
c  This only works if enough memory remains.
c
      call memsizr(junk,memremr,0)
      if(nreal-curarsz.gt.memremr) then
        write(*,1010) 'not enough memory, nreal-curarsz,memremr:',
     .                                    nreal-curarsz,memremr
        memrear = -7
        goto 70
      endif
c
c  Finally, ready to reallocate.
c
      if(realmem(address+curarsz+1).eq.0.0) then
c
c  It's trivial if this is the last array.
c
        realmem(address-1)       = nreal
        realmem(address+nreal)   = nreal
        realmem(address+nreal+1) = 0
        memrear                  = handle
      else
        shift = nreal - curarsz
        if(shift.lt.0) then
c
c  The new size is smaller than the current size, so start at the array
c  after handle and shift arrays until the end.
c
          realmem(address-1)       = nreal
          realmem(address+nreal)   = nreal
          memrear                  = handle
          next                     = address + curarsz + 1
   10     continue
            if(next.lt.1.or.next.gt.MXALLOCR) then
              write(*,1010) 'next out of range:',next,MXALLOCR
              memrear = -8
              goto 70
            endif
            if(realmem(next).eq.0.0) goto 30
c
c  Adjust the handle.
c
            do i=1,MXNUMARR
              if(iarrayr(i).eq.next+1) then
                iarrayr(i) = iarrayr(i) + shift
                goto 20
              endif
            end do
   20       continue
c
c  Shift the array back.
c
            curarsz = realmem(next)
            do i=1,curarsz+2
              realmem(i+(next+shift-1)) = realmem(i+(next-1))
            end do
            next = next + curarsz + 2
            goto 10
   30     continue
          realmem(next+shift) = 0.0
        else
c
c  The new size is larger than the current size, so start at the last
c  array and shift arrays until the next one after handle.
c
          next = memchkr()
          if(next.le.1) then
            write(*,1010) 'memchkr returned bad value:',next
            memrear = -9
            goto 70
          endif
          if(next+shift.lt.1.or.next+shift.gt.MXALLOCR) then
            write(*,1010) 'next+shift out of range:',next+shift,MXALLOCR
            memrear = -10
            goto 70
          endif
          realmem(next+shift) = 0.0
          next = next - realmem(next-1) - 2
   40     continue
            if(next.lt.1.or.next.ge.MXALLOCR-2) then
              write(*,1010) 'next out of range:',next,MXALLOCR
              memrear = -11
              goto 70
            endif
            if(next+1.eq.address) goto 60
c
c  Adjust the handle.
c
            do i=1,MXNUMARR
              if(iarrayr(i).eq.next+1) then
                iarrayr(i) = iarrayr(i) + shift
                goto 50
              endif
            end do
   50       continue
c
c  Shift the array forward.
c
            curarsz = realmem(next)
            do i=curarsz+2,1,-1
              realmem(i+(next+shift-1)) = realmem(i+(next-1))
            end do
            next = next - realmem(next-1) - 2
            goto 40
   60     continue
          realmem(address-1)       = nreal
          realmem(address+nreal)   = nreal
          memrear                  = handle
        endif
      endif
   70 continue
 1010 format(' MEMREAR:  ERROR, ',a,4i10)
c
      return
      end

c-----------------------------------------------------------------------

      subroutine memsizr(handrem,memrem,iprint)

      implicit real*8 (a-h,o-z)
      integer handrem
      parameter (MXNUMARR=4000,MXNUMREAL=100000000)
      parameter (MXALLOCR=MXNUMREAL+2*MXNUMARR+1)
      common /meminte/ narrayr,iarrayr(MXNUMARR)
      common /memreal/ realmem(MXALLOCR)
c
c  Scott D. Thomas, Sterling Software, Inc.
c  Contract NAS2-13210/6/5, NASA Ames Research Center.
c  21-NOV-1996..4-DEC-1996.
c
c  Tell how much memory is left for array allocation.  handrem will be
c  set to the remaining number of array handles available.  The value of
c  memrem will be the current maximum number of real*8's that might be
c  allocated if all the handles were used.
c
c  If iprint is set to a nonzero value then a message will be written
c  to the standard output channel (*).
c
      handrem = MXNUMARR - narrayr
      memused = 0
      next    = 1
   10 continue
        if(realmem(next).eq.0) goto 20
        memused = memused + realmem(next)
        next    = next + realmem(next)+2
        goto 10
   20 continue
      memrem = max(0,MXALLOCR-(2*MXNUMARR+1)-memused)
      if(iprint.ne.0) write(*,1010) handrem,memrem
 1010 format(' MEMSIZR:',i5,'  real*8 array handles left, room for',
     .       i10,' real*8''s.')
c
      return
      end

c-----------------------------------------------------------------------

      subroutine memuser(handuse,memsum,iprint)

      implicit real*8 (a-h,o-z)
      integer handuse
      parameter (MXNUMARR=4000,MXNUMREAL=100000000)
      parameter (MXALLOCR=MXNUMREAL+2*MXNUMARR+1)
      common /meminte/ narrayr,iarrayr(MXNUMARR)
      common /memreal/ realmem(MXALLOCR)
c
c  Scott D. Thomas, Sterling Software, Inc.
c  Contract NAS2-13210/6/5, NASA Ames Research Center.
c  21-NOV-1996..4-DEC-1996.
c
c  Tell how much memory is used for array allocation.  handuse will be
c  set to the number of array handles currently being used.  The value
c  of memsum will be the sum of the sizes of currently allocated real*8
c  arrays.
c
c  If iprint is not zero, print a table of real memory use to the
c  standard output channel (*).
c
      if(iprint.ne.0) then
        if(narrayr.eq.0) then
          write(*,1010)
 1010     format(' MEMUSER: There are no  real*8 arrays currently',
     .           ' allocated.')
        else if(narrayr.eq.1) then
          write(*,1011)
 1011     format(' MEMUSER: There is 1  real*8 array currently',
     .           ' allocated.')
        else
          write(*,1012) narrayr
 1012     format(' MEMUSER: There are',i4,'  real*8 arrays currently',
     .           ' allocated.')
        endif
        if(narrayr.ge.1) write(*,1020)
 1020   format(' MEMUSER:   Rank  Handle   FirstValue   Address',
     .         '      Size   SizeSum')
      endif
      nrank    = 0
      handuse  = narrayr
      memsum   = 0
      next     = 1
   10 continue
        if(realmem(next).eq.0.0) goto 30
        nrank   = nrank + 1
        memsum  = memsum + realmem(next)
        if(iprint.ne.0) then
          nhandle = 0
          do i=1,MXNUMARR
            if(iarrayr(i).eq.next+1) then
              nhandle = i
              goto 20
            endif
          end do
   20     continue
          isize  = realmem(next)
          rfirst = realmem(next+1)
          write(*,1030) nrank,nhandle,rfirst,next+1,isize,memsum
 1030     format(' MEMUSER:',i7,i8,1pe13.5,3i10)
        endif
        next = next + realmem(next) + 2
        goto 10
   30 continue
c
      return
      end

c-----------------------------------------------------------------------

      function lenstr(string)

      character(*) string*(*)
c
c  Return length of string less trailing blanks.
c
c  Scott D. Thomas, Sterling Software, Inc.
c  Contract NAS 2-13210/6/5, NASA Ames Research Center.
c
      length = len(string)
      lenstr = length
      do 10 i=length,1,-1
        if(string(i:i).ne.' ') goto 20
        lenstr = i-1
   10 continue
   20 continue
c
      return
      end

c-----------------------------------------------------------------------

      subroutine gridprop(inunit,ismgrid,nrpgrid)
c
c  Read inunit to get grid properties:
c
c  single block or multiple blocks (ismgrid = .false. or .true.)
c  real*4 or real*8 precision (nrpgrid = 4 or 8)
c
c  Test read only the first block.  Rewind inunit before and after.
c
c  Scott D. Thomas, Sterling Software, Inc.
c  NASA Ames Research Center, Contract NAS 2-13210/6/3.
c
      logical ismgrid
      data machisz,machrsz,machr4sz/4,8,4/
      dimension   nj(1),   nk(1),   nl(1)
      pointer   (pnj,nj),(pnk,nk),(pnl,nl)
      real*4      xyz4(3,1)
      pointer   (pxyz4,xyz4)
      real*8      xyz8(3,1)
      pointer   (pxyz8,xyz8)
c
      ismgrid = .true.
      rewind(inunit)
      read(inunit,err=10) i1,i2,i3
      ismgrid = .false.
      nb = 1
   10 continue
      rewind(inunit)
      if(ismgrid) then
        read(inunit) nb
        if(nb.lt.1) then
          write(*,*) 'NB=',nb
          stop 'NB'
        endif
      endif
      numr = max(1,(nb+1)/(machrsz/machisz))
      call ema(pnj,numr,knj,'nj',1)
      call ema(pnk,numr,knk,'nk',2)
      call ema(pnl,numr,knl,'nl',3)
      read(inunit) (nj(n),nk(n),nl(n),n=1,nb)
      if(nj(1).lt.1.or.nk(1).lt.1.or.nl(1).lt.1) then
        write(*,*) 'NB,NJ,NK,NL=',NB,NJ,NK,NL
        stop 'NJ,NK,NL'
      endif
      njkl = nj(1)*nk(1)*nl(1)
      if(njkl.lt.1) stop 'njkl'
      numr = max(1,(3*njkl+1)/(machrsz/machrsz))
      call ema(pxyz8,numr,kpxyz8,'pxyz8',4)
      read(inunit,err=30) ((xyz8(m,i),i=1,njkl),m=1,3)
      nrpgrid = 8
      call emf(kpxyz8,'pxyz8',5)
      goto 40
   30 continue
      call emf(kpxyz8,'pxyz8',6)
      rewind(inunit)
      if(ismgrid) read(inunit) nb
      read(inunit) (nj(n),nk(n),nl(n),n=1,nb)
      numr = max(1,(3*njkl+1)/(machrsz/machr4sz))
      call ema(pxyz4,numr,kpxyz4,'pxyz4',7)
      read(inunit) ((xyz4(m,i),i=1,njkl),m=1,3)
      nrpgrid = 4
      call emf(kpxyz4,'pxyz4',8)
   40 continue
      call emf(knl,'nj',9)
      call emf(knk,'nj',10)
      call emf(knj,'nj',11)
      rewind(inunit)
c
      if(ismgrid.and.nrpgrid.eq.4) then
        write(*,*) 'It appears to be a real*4 multiple block grid file.'
      else if(ismgrid.and.nrpgrid.eq.8) then
        write(*,*) 'It appears to be a real*8 multiple block grid file.'
      else if(.not.ismgrid.and.nrpgrid.eq.4) then
        write(*,*) 'It appears to be a real*4 single grid file.'
      else if(.not.ismgrid.and.nrpgrid.eq.8) then
        write(*,*) 'It appears to be a real*8 single grid file.'
      else
        stop 'never'
      endif
c
      return
      end

c-----------------------------------------------------------------------

      logical function matches(xtrema,nc,mc,epsilon,nj,nk,nl,nb)

      implicit real*8 (a-h,o-z)
      dimension xtrema(6,nb*6)
      integer   nj(nb),nk(nb),nl(nb)

      matches = .false.
      if(abs(xtrema(1,nc)-xtrema(1,mc)).le.epsilon.and.
     .   abs(xtrema(2,nc)-xtrema(2,mc)).le.epsilon.and.
     .   abs(xtrema(3,nc)-xtrema(3,mc)).le.epsilon.and.
     .   abs(xtrema(4,nc)-xtrema(4,mc)).le.epsilon.and.
     .   abs(xtrema(5,nc)-xtrema(5,mc)).le.epsilon.and.
     .   abs(xtrema(6,nc)-xtrema(6,mc)).le.epsilon) then
        n  = (nc-1)/6+1
        nf = nc-6*(n-1)
        if(nf.le.2) then
          n1 = nk(n)
          n2 = nl(n)
        else if(nf.ge.5) then
          n1 = nj(n)
          n2 = nk(n)
        else
          n1 = nj(n)
          n2 = nl(n)
        endif
        m  = (mc-1)/6+1
        mf = mc-6*(m-1)
        if(mf.le.2) then
          m1 = nk(m)
          m2 = nl(m)
        else if(mf.ge.5) then
          m1 = nj(m)
          m2 = nk(m)
        else
          m1 = nj(m)
          m2 = nl(m)
        endif
        if((n1.eq.m1.and.n2.eq.m2).or.
     .     (n1.eq.m2.and.n2.eq.m1)) then
          matches = .true.
        else
          write(*,1010) n,nf,m,mf,n1,n2,m1,m2
 1010     format(' WARNING:  BLK',i3,' F',i2,' near BLK',i3,' F',i2,
     .           ' but dimensions are',i4,',',i4,' and',i4,',',i4,'.')
        endif
      endif
c
      return
      end

c-----------------------------------------------------------------------

      subroutine geticon(icon,nc,mc,epsilon,xyz,ndf,ndc,ihand,
     .                   nguesses,strings,imatch,MSTRINGS)
      implicit real*8 (a-h,o-z)
      integer   icon(4),ndf(4,*),ndc(4,*),ihand(*),ndq(4),mdq(4),it(2)
      dimension xyz(3,*)
      character string*20, strings(MSTRINGS)*20
c
c  Figure out the four connection integers for the nth block's nf face,
c  which has already been found to be close to the mth block's mf face.
c  where m,n range from 1 to nb and mf,nf range from 1 to 6.
c  The nth block's nf face is component nc=nf+6*(n-1) and
c  the mth block's mf face is component mc=mf+6*(m-1).
c  The ihand array records the handedness of the nth block.
c
      n  = (nc-1)/6+1
      m  = (mc-1)/6+1
      nf = nc-6*(n-1)
      mf = mc-6*(m-1)
c
c  Default values (show up in the output case of programming error).
c
      icon(1) = m
      icon(2) = 99
      icon(3) = 99
      icon(4) = 99
c
c  Put the four corners of nc into ndq, and
c  put the four corners of mc into mdq.
c
      do i=1,4
        ndq(i) = ndc(i,nc)
        mdq(i) = ndc(i,mc)
      end do
c
c  There are 8 ways for two faces to align, represented by (ia,ib) where
c  the first coordinate describes which of the abutting face's two
c  basis vectors align with it, and the second coordinate is for the
c  remaining basis vector.  The sign indicates direction.  They are:
c
c  (1,2),(2,-1),(-1,-2),(-2,1)
c  (2,1),(-1,2),(-2,-1),(1,-2)
c
      ia = 0
      if(dist(xyz(1,ndq(1)),xyz(1,mdq(1))).le.epsilon) then
        if(dist(xyz(1,ndq(2)),xyz(1,mdq(2))).le.epsilon.and.
     .     dist(xyz(1,ndq(3)),xyz(1,mdq(3))).le.epsilon.and.
     .     dist(xyz(1,ndq(4)),xyz(1,mdq(4))).le.epsilon) then
          ia = 1
          ib = 2
        else if(dist(xyz(1,ndq(2)),xyz(1,mdq(4))).le.epsilon.and.
     .          dist(xyz(1,ndq(3)),xyz(1,mdq(3))).le.epsilon.and.
     .          dist(xyz(1,ndq(4)),xyz(1,mdq(2))).le.epsilon) then
          ia = 2
          ib = 1
        endif
      else if(dist(xyz(1,ndq(1)),xyz(1,mdq(2))).le.epsilon) then
        if(dist(xyz(1,ndq(2)),xyz(1,mdq(3))).le.epsilon.and.
     .     dist(xyz(1,ndq(3)),xyz(1,mdq(4))).le.epsilon.and.
     .     dist(xyz(1,ndq(4)),xyz(1,mdq(1))).le.epsilon) then
          ia = 2
          ib = -1
        else if(dist(xyz(1,ndq(2)),xyz(1,mdq(1))).le.epsilon.and.
     .          dist(xyz(1,ndq(3)),xyz(1,mdq(4))).le.epsilon.and.
     .          dist(xyz(1,ndq(4)),xyz(1,mdq(3))).le.epsilon) then
          ia = -1
          ib = 2
        endif
      else if(dist(xyz(1,ndq(1)),xyz(1,mdq(3))).le.epsilon) then
        if(dist(xyz(1,ndq(2)),xyz(1,mdq(4))).le.epsilon.and.
     .     dist(xyz(1,ndq(3)),xyz(1,mdq(1))).le.epsilon.and.
     .     dist(xyz(1,ndq(4)),xyz(1,mdq(2))).le.epsilon) then
          ia = -1
          ib = -2
        else if(dist(xyz(1,ndq(2)),xyz(1,mdq(2))).le.epsilon.and.
     .          dist(xyz(1,ndq(3)),xyz(1,mdq(1))).le.epsilon.and.
     .          dist(xyz(1,ndq(4)),xyz(1,mdq(4))).le.epsilon) then
          ia = -2
          ib = -1
        endif
      else if(dist(xyz(1,ndq(1)),xyz(1,mdq(4))).le.epsilon) then
        if(dist(xyz(1,ndq(2)),xyz(1,mdq(1))).le.epsilon.and.
     .     dist(xyz(1,ndq(3)),xyz(1,mdq(2))).le.epsilon.and.
     .     dist(xyz(1,ndq(4)),xyz(1,mdq(3))).le.epsilon) then
          ia = -2
          ib = 1
        else if(dist(xyz(1,ndq(2)),xyz(1,mdq(3))).le.epsilon.and.
     .          dist(xyz(1,ndq(3)),xyz(1,mdq(2))).le.epsilon.and.
     .          dist(xyz(1,ndq(4)),xyz(1,mdq(1))).le.epsilon) then
          ia = 1
          ib = -2
        endif
      endif
c
c  If these attempts all missed then print a warning and guess.
c
      if(ia.eq.0) then
c
c  Make a guess based on proximity of points.
c
        d11 = dist(xyz(1,ndq(1)),xyz(1,mdq(1)))
        d12 = dist(xyz(1,ndq(1)),xyz(1,mdq(2)))
        d13 = dist(xyz(1,ndq(1)),xyz(1,mdq(3)))
        d14 = dist(xyz(1,ndq(1)),xyz(1,mdq(4)))
        e1  = min(d11,d12,d13,d14)
        d21 = dist(xyz(1,ndq(2)),xyz(1,mdq(1)))
        d22 = dist(xyz(1,ndq(2)),xyz(1,mdq(2)))
        d23 = dist(xyz(1,ndq(2)),xyz(1,mdq(3)))
        d24 = dist(xyz(1,ndq(2)),xyz(1,mdq(4)))
        e2  = min(d21,d22,d23,d24)
        d31 = dist(xyz(1,ndq(3)),xyz(1,mdq(1)))
        d32 = dist(xyz(1,ndq(3)),xyz(1,mdq(2)))
        d33 = dist(xyz(1,ndq(3)),xyz(1,mdq(3)))
        d34 = dist(xyz(1,ndq(3)),xyz(1,mdq(4)))
        e3  = min(d31,d32,d33,d34)
        d41 = dist(xyz(1,ndq(4)),xyz(1,mdq(1)))
        d42 = dist(xyz(1,ndq(4)),xyz(1,mdq(2)))
        d43 = dist(xyz(1,ndq(4)),xyz(1,mdq(3)))
        d44 = dist(xyz(1,ndq(4)),xyz(1,mdq(4)))
        e4  = min(d41,d42,d43,d44)
        eps = max(e1,e2,e3,e4)
c
c  (1,2),(2,-1),(-1,-2),(-2,1)
c  (2,1),(-1,2),(-2,-1),(1,-2)
c
        if(dist(xyz(1,ndq(1)),xyz(1,mdq(1))).le.e1) then
          if(dist(xyz(1,ndq(2)),xyz(1,mdq(2))).le.eps.and.
     .       dist(xyz(1,ndq(3)),xyz(1,mdq(3))).le.eps.and.
     .       dist(xyz(1,ndq(4)),xyz(1,mdq(4))).le.eps) then
            ia = 1
            ib = 2
          else if(dist(xyz(1,ndq(2)),xyz(1,mdq(4))).le.eps.and.
     .            dist(xyz(1,ndq(3)),xyz(1,mdq(3))).le.eps.and.
     .            dist(xyz(1,ndq(4)),xyz(1,mdq(2))).le.eps) then
            ia = 2
            ib = 1
          endif
        else if(dist(xyz(1,ndq(1)),xyz(1,mdq(2))).le.e1) then
          if(dist(xyz(1,ndq(2)),xyz(1,mdq(3))).le.eps.and.
     .       dist(xyz(1,ndq(3)),xyz(1,mdq(4))).le.eps.and.
     .       dist(xyz(1,ndq(4)),xyz(1,mdq(1))).le.eps) then
            ia = 2
            ib = -1
          else if(dist(xyz(1,ndq(2)),xyz(1,mdq(1))).le.eps.and.
     .            dist(xyz(1,ndq(3)),xyz(1,mdq(4))).le.eps.and.
     .            dist(xyz(1,ndq(4)),xyz(1,mdq(3))).le.eps) then
            ia = -1
            ib = 2
          endif
        else if(dist(xyz(1,ndq(1)),xyz(1,mdq(3))).le.e1) then
          if(dist(xyz(1,ndq(2)),xyz(1,mdq(4))).le.eps.and.
     .       dist(xyz(1,ndq(3)),xyz(1,mdq(1))).le.eps.and.
     .       dist(xyz(1,ndq(4)),xyz(1,mdq(2))).le.eps) then
            ia = -1
            ib = -2
          else if(dist(xyz(1,ndq(2)),xyz(1,mdq(2))).le.eps.and.
     .            dist(xyz(1,ndq(3)),xyz(1,mdq(1))).le.eps.and.
     .            dist(xyz(1,ndq(4)),xyz(1,mdq(4))).le.eps) then
            ia = -2
            ib = -1
          endif
        else if(dist(xyz(1,ndq(1)),xyz(1,mdq(4))).le.e1) then
          if(dist(xyz(1,ndq(2)),xyz(1,mdq(1))).le.eps.and.
     .       dist(xyz(1,ndq(3)),xyz(1,mdq(2))).le.eps.and.
     .       dist(xyz(1,ndq(4)),xyz(1,mdq(3))).le.eps) then
            ia = -2
            ib = 1
          else if(dist(xyz(1,ndq(2)),xyz(1,mdq(3))).le.eps.and.
     .            dist(xyz(1,ndq(3)),xyz(1,mdq(2))).le.eps.and.
     .            dist(xyz(1,ndq(4)),xyz(1,mdq(1))).le.eps) then
            ia = 1
            ib = -2
          endif
        endif
        write(*,1011) n,nf,m,mf,ia,ib
 1011   format(' WARNING:  BLK',i3,' F',i2,' <--> BLK',i3,' F',i2,
     .         ' corners don''t match, guessing ia,ib=',i2,',',i2)
c
c  There should be at least one guess, but just in case...
c
        if(ia.eq.0) then
          write(*,*) 'ERROR:  ia.eq.0'
          stop 'geticon'
        endif
        nguesses = nguesses + 1
      endif
c
c  Fill in icon(2..4) with the coded neighbor basis vectors.
c
c  strings looks like this:
c
c  nfmfiaibnhnfc1c2c3nh
c   2 1 1 2 1 2 1 2 3 1
c   2 1-2 1 1 2 1-3 2 1
c   2 1-1-2 1 2 1-2-3 1
c   2 1 2-1 1 2 1 3-2 1
c   [...]
c
      write(string,1030) nf,mf,ia,ib,ihand(n)
 1030 format(5i2)
c
c  Given the first 10 characters, find a match in strings(i)(1:10), then
c  the 6 characters in strings(i)(13:18) contain the icon values.
c  Kick out after first find, no need to go all the way down the list.
c
      nmatch = 0
      do i=1,MSTRINGS
        if(string(1:10).eq.strings(i)(1:10)) then
          nmatch = nmatch + 1
          imatch = i
          goto 10
        endif
      end do
   10 continue
      if(nmatch.eq.0) then
        write(*,1040) n,nf,m,mf,string(1:10)
 1040   format(' BLK',i4,' F',i2,' <--> BLK',i4,' F',i2,
     .         ' but no match was found to string: ',a)
        stop 'geticon'
      else if(nmatch.ne.1) then
        write(*,1050) n,nf,m,mf,nmatch,imatch,strings(i)
 1050   format(' BLK',i4,' F',i2,' <--> BLK',i4,' F',i2,
     .         ' nmatch=',i3,', strings(',i4,')=',a)
        write(*,*) 'ERROR:  found more than one match.'
        stop 'geticon'
      endif
c
c  All is well, get icon values from the rest of strings(imatch).
c
      read(strings(imatch)(13:18),1060) (icon(j),j=2,4)
 1060 format(3i2)
c
      return
      end

c-----------------------------------------------------------------------

      integer function ihand4(xyz4,nj,nk,nl,n)

      implicit real*8 (a-h,o-z)
      real*4 xyz4(3,nj,nk,nl)
c
c  Compute the volume of a hex cell using six tets.
c  Number the vertices as follows.
c
c  (1)=(j  ,k  ,l  )   (2)=(j+1,k  ,l  )
c  (3)=(j+1,k+1,l  )   (4)=(j  ,k+1,l  )
c  (5)=(j  ,k  ,l+1)   (6)=(j+1,k  ,l+1)
c  (7)=(j+1,k+1,l+1)   (8)=(j  ,k+1,l+1)
c
c  The six tets are 1245, 7368, 4528, 5628, 6328, 3428.
c  Order matters, these orientations give positive volumes for a right-
c  handed convex hex cell.
c
      volsum = 0.0
      do l=1,nl-1
        do k=1,nk-1
          do j=1,nj-1
            v1     = voltet4(xyz4(1,j  ,k  ,l  ),xyz4(1,j+1,k  ,l  ),
     .                       xyz4(1,j  ,k+1,l  ),xyz4(1,j  ,k  ,l+1))
            v2     = voltet4(xyz4(1,j+1,k+1,l+1),xyz4(1,j+1,k+1,l  ),
     .                       xyz4(1,j+1,k  ,l+1),xyz4(1,j  ,k+1,l+1))
            v3     = voltet4(xyz4(1,j  ,k+1,l  ),xyz4(1,j  ,k  ,l+1),
     .                       xyz4(1,j+1,k  ,l  ),xyz4(1,j  ,k+1,l+1))
            v4     = voltet4(xyz4(1,j  ,k  ,l+1),xyz4(1,j+1,k  ,l+1),
     .                       xyz4(1,j+1,k  ,l  ),xyz4(1,j  ,k+1,l+1))
            v5     = voltet4(xyz4(1,j+1,k  ,l+1),xyz4(1,j+1,k+1,l  ),
     .                       xyz4(1,j+1,k  ,l  ),xyz4(1,j  ,k+1,l+1))
            v6     = voltet4(xyz4(1,j+1,k+1,l  ),xyz4(1,j  ,k+1,l  ),
     .                       xyz4(1,j+1,k  ,l  ),xyz4(1,j  ,k+1,l+1))
            vol6 = v1 + v2 + v3 + v4 + v5 + v6
            if(volsum.gt.0.0.and.vol6.lt.0.0) then
              write(*,1010) n
 1010         format(' WARNING:  inverted +- volume in block',i4)
            else if(volsum.lt.0.0.and.vol6.gt.0.0) then
              write(*,1011) n
 1011         format(' WARNING:  inverted -+ volume in block',i4)
            endif
            volsum = volsum + vol6
          end do
        end do
      end do
      if(volsum.gt.0.0) then
        ihand4 = 1
      else if(volsum.eq.0.0) then
        write(*,1020) n
 1020   format(' WARNING:  volume is zero in block',i4)
        ihand4 = 0
      else
        ihand4 = -1
      endif
      return
      end

c-----------------------------------------------------------------------

      integer function ihand8(xyz8,nj,nk,nl,n)

      implicit real*8 (a-h,o-z)
      real*8 xyz8(3,nj,nk,nl)
c
c  Compute the volume of a hex cell using six tets.
c  Number the vertices as follows.
c
c  (1)=(j  ,k  ,l  )   (2)=(j+1,k  ,l  )
c  (3)=(j+1,k+1,l  )   (4)=(j  ,k+1,l  )
c  (5)=(j  ,k  ,l+1)   (6)=(j+1,k  ,l+1)
c  (7)=(j+1,k+1,l+1)   (8)=(j  ,k+1,l+1)
c
c  The six tets are 1245, 7368, 4528, 5628, 6328, 3428.
c  Order matters, these orientations give positive volumes for a right-
c  handed convex hex cell.
c
      volsum = 0.0
      do l=1,nl-1
        do k=1,nk-1
          do j=1,nj-1
            v1     = voltet8(xyz8(1,j  ,k  ,l  ),xyz8(1,j+1,k  ,l  ),
     .                       xyz8(1,j  ,k+1,l  ),xyz8(1,j  ,k  ,l+1))
            v2     = voltet8(xyz8(1,j+1,k+1,l+1),xyz8(1,j+1,k+1,l  ),
     .                       xyz8(1,j+1,k  ,l+1),xyz8(1,j  ,k+1,l+1))
            v3     = voltet8(xyz8(1,j  ,k+1,l  ),xyz8(1,j  ,k  ,l+1),
     .                       xyz8(1,j+1,k  ,l  ),xyz8(1,j  ,k+1,l+1))
            v4     = voltet8(xyz8(1,j  ,k  ,l+1),xyz8(1,j+1,k  ,l+1),
     .                       xyz8(1,j+1,k  ,l  ),xyz8(1,j  ,k+1,l+1))
            v5     = voltet8(xyz8(1,j+1,k  ,l+1),xyz8(1,j+1,k+1,l  ),
     .                       xyz8(1,j+1,k  ,l  ),xyz8(1,j  ,k+1,l+1))
            v6     = voltet8(xyz8(1,j+1,k+1,l  ),xyz8(1,j  ,k+1,l  ),
     .                       xyz8(1,j+1,k  ,l  ),xyz8(1,j  ,k+1,l+1))
            vol6 = v1 + v2 + v3 + v4 + v5 + v6
            if(volsum.gt.0.0.and.vol6.lt.0.0) then
              write(*,1010) n
 1010         format(' WARNING:  inverted +- volume in block',i4)
            else if(volsum.lt.0.0.and.vol6.gt.0.0) then
              write(*,1011) n
 1011         format(' WARNING:  inverted -+ volume in block',i4)
            endif
            volsum = volsum + vol6
          end do
        end do
      end do
      if(volsum.gt.0.0) then
        ihand8 = 1
      else if(volsum.eq.0.0) then
        write(*,1020) n
 1020   format(' WARNING:  volume is zero in block',i4)
        ihand8 = 0
      else
        ihand8 = -1
      endif
      return
      end

c-----------------------------------------------------------------------

      real*8 function voltet4(a,b,c,d)

      implicit real*8 (a-h,o-z)
      real*4 a(3),b(3),c(3),d(3)
      dimension u(3),v(3),w(3)
c
c  Compute volume of tet on points a,b,c,d.
c
      do m=1,3
        u(m) = b(m) - a(m)
        v(m) = c(m) - a(m)
        w(m) = d(m) - a(m)
      end do
      voltet4 = (u(1) * (v(2) * w(3) - v(3) * w(2)) -
     .           u(2) * (v(1) * w(3) - v(3) * w(1)) +
     .           u(3) * (v(1) * w(2) - v(2) * w(1))) / 6.0
      return
      end

c-----------------------------------------------------------------------

      real*8 function voltet8(a,b,c,d)

      implicit real*8 (a-h,o-z)
      real*8 a(3),b(3),c(3),d(3)
      dimension u(3),v(3),w(3)
c
c  Compute volume of tet on points a,b,c,d.
c
      do m=1,3
        u(m) = b(m) - a(m)
        v(m) = c(m) - a(m)
        w(m) = d(m) - a(m)
      end do
      voltet8 = (u(1) * (v(2) * w(3) - v(3) * w(2)) -
     .           u(2) * (v(1) * w(3) - v(3) * w(1)) +
     .           u(3) * (v(1) * w(2) - v(2) * w(1))) / 6.0
      return
      end

c-----------------------------------------------------------------------

      real*8 function dist(x,y)

      implicit real*8 (a-h,o-z)
      dimension x(3),y(3)
      dist = sqrt((x(1)-y(1))**2+(x(2)-y(2))**2+(x(3)-y(3))**2)
      return
      end

c-----------------------------------------------------------------------

      subroutine qhitter(icon,nc,mc,epsilon,xyz,ndf,iaf,
     .                   nj,nk,nl,nhitquad,nquads,string)
      implicit real*8 (a-h,o-z)
      integer   icon(4),ndf(4,*),nj(*),nk(*),nl(*),iaf(*)
      dimension xyz(3,*)
      character string*20
      logical qmatches
c
c  For each quad in the first component (nc) look to see whether it hits
c  a quad in the second component (mc) in the proper order but any
c  vertex orientation within epsilon distance for all four points.
c
c  The number of quads (nquads) will be incremented by the number of
c  candidate quads in the face (nc) only if the icon(1) value is
c  positive, indicating a block-to-block abutment between components nc
c  and mc.
c
c  The nhitquad value will be incremented by 1 for each quad in nc that
c  matches the quad in mc by epsilon at all four vertex points, but if
c  the number of quads in mc is different from the number in nc the
c  search will not be done and the nhitquad value will not be
c  incremented at all.  Again, the search is done only if icon(1).gt.0.
c
c  string.eq.strings(imatch) contains nf,mf,ia,ib,nh,nf,icon(3..4),nh,
c  as 10i2 format.
c
c  Program stops if the component sizes (number of quads) do not match.
c
c  Abutting blocks' sides are indicated by positive icon(1).
c
      if(icon(1).le.0) return
c
c  The total number of quads should be hit for a perfect match.
c
      nquads = nquads + iaf(nc+1)-iaf(nc)
c
c  Check whether quads match at their vertices within epsilon. 
c  Get block numbers, face numbers, and face dimensions.
c
      n  = (nc-1)/6+1
      m  = (mc-1)/6+1
      nf = nc-6*(n-1)
      mf = mc-6*(m-1)
      if(nf.eq.1.or.nf.eq.2) then
        n1 = nk(n)
        n2 = nl(n)
      else if(nf.eq.3.or.nf.eq.4) then
        n1 = nj(n)
        n2 = nl(n)
      else
        n1 = nj(n)
        n2 = nk(n)
      endif
      if(mf.eq.1.or.mf.eq.2) then
        m1 = nk(m)
        m2 = nl(m)
      else if(mf.eq.3.or.mf.eq.4) then
        m1 = nj(m)
        m2 = nl(m)
      else
        m1 = nj(m)
        m2 = nk(m)
      endif
c
c  The abutting face dimensions should match.
c
      if(n1*n2.ne.m1*m2) then
        write(*,*) 'ERROR:  n1*n2.ne.m1*m2'
        write(*,*) 'n1,n2,n1*n2=',n1,n2,n1*n2
        write(*,*) 'm1,m2,m1*m2=',m1,m2,m1*m2
        write(*,*) 'n,nf,nc    =',n,nf,nc
        write(*,*) 'm,mf,mc    =',m,mf,mc
        stop 'qhitter'
      endif
c
c  The number of quads in abutting faces should match.
c
      if(iaf(nc+1)-iaf(nc).ne.iaf(mc+1)-iaf(mc)) then
        write(*,*) 'ERROR:  iaf(nc+1)-iaf(nc).ne.iaf(mc+1)-iaf(mc)'
        write(*,*) 'iaf(nc+1),iaf(nc),iaf(nc+1)-iaf(nc)=',
     .              iaf(nc+1),iaf(nc),iaf(nc+1)-iaf(nc)
        write(*,*) 'iaf(mc+1),iaf(mc),iaf(mc+1)-iaf(mc)=',
     .              iaf(mc+1),iaf(mc),iaf(mc+1)-iaf(mc)
        write(*,*) 'n,nf,nc                            =',n,nf,nc
        write(*,*) 'm,mf,mc                            =',m,mf,mc
        stop 'qhitter'
      endif
c
c  Get ia and ib out of string to determine face-to-face orientation.
c  The nf and mf in string should match values computed above.
c
      read(string,1010) nfs,mfs,ia,ib
 1010 format(4i2)
      if(nfs.ne.nf.or.mfs.ne.mf) then
        write(*,*) 'ERROR:  nf,mf do not match string.'
        write(*,*) 'n,nf,nc=',n,nf,nc
        write(*,*) 'm,mf,mc=',m,mf,mc
        write(*,*) 'nfs,mfs=',nfs,mfs
        stop 'qhitter'
      endif
c
c  There are 8 ways for two faces to align, represented by (ia,ib) where
c  the first coordinate describes which of the abutting face's two
c  basis vectors align with it, and the second coordinate is for the
c  remaining basis vector.  The sign indicates direction.  They are:
c
c  (1,2),(2,-1),(-1,-2),(-2,1)
c  (2,1),(-1,2),(-2,-1),(1,-2)
c
c
      it = iaf(nc)-1
      if(ia.eq.1.and.ib.eq.2) then
        do j2=1,m2-1
          do j1=1,m1-1
            it = it + 1
            jt = iaf(mc)-1+j1+(m1-1)*(j2-1)
            if(qmatches(ndf(1,it),ndf(1,jt),xyz,epsilon)) then
              nhitquad = nhitquad + 1
            endif
          end do
        end do
      else if(ia.eq.2.and.ib.eq.1) then
        do j1=1,m1-1
          do j2=1,m2-1
            it = it + 1
            jt = iaf(mc)-1+j1+(m1-1)*(j2-1)
            if(qmatches(ndf(1,it),ndf(1,jt),xyz,epsilon)) then
              nhitquad = nhitquad + 1
            endif
          end do
        end do
      else if(ia.eq.2.and.ib.eq.-1) then
        do j1=m1-1,1,-1
          do j2=1,m2-1
            it = it + 1
            jt = iaf(mc)-1+j1+(m1-1)*(j2-1)
            if(qmatches(ndf(1,it),ndf(1,jt),xyz,epsilon)) then
              nhitquad = nhitquad + 1
            endif
          end do
        end do
      else if(ia.eq.-1.and.ib.eq.2) then
        do j2=1,m2-1
          do j1=m1-1,1,-1
            it = it + 1
            jt = iaf(mc)-1+j1+(m1-1)*(j2-1)
            if(qmatches(ndf(1,it),ndf(1,jt),xyz,epsilon)) then
              nhitquad = nhitquad + 1
            endif
          end do
        end do
      else if(ia.eq.-1.and.ib.eq.-2) then
        do j2=m2-1,1,-1
          do j1=m1-1,1,-1
            it = it + 1
            jt = iaf(mc)-1+j1+(m1-1)*(j2-1)
            if(qmatches(ndf(1,it),ndf(1,jt),xyz,epsilon)) then
              nhitquad = nhitquad + 1
            endif
          end do
        end do
      else if(ia.eq.-2.and.ib.eq.-1) then
        do j1=m1-1,1,-1
          do j2=m2-1,1,-1
            it = it + 1
            jt = iaf(mc)-1+j1+(m1-1)*(j2-1)
            if(qmatches(ndf(1,it),ndf(1,jt),xyz,epsilon)) then
              nhitquad = nhitquad + 1
            endif
          end do
        end do
      else if(ia.eq.-2.and.ib.eq.1) then
        it = iaf(nc)-1
        do j1=1,m1-1
          do j2=m2-1,1,-1
            it = it + 1
            jt = iaf(mc)-1+j1+(m1-1)*(j2-1)
            if(qmatches(ndf(1,it),ndf(1,jt),xyz,epsilon)) then
              nhitquad = nhitquad + 1
            endif
          end do
        end do
      else if(ia.eq.1.and.ib.eq.-2) then
        it = iaf(nc)-1
        do j2=m2-1,1,-1
          do j1=1,m1-1
            it = it + 1
            jt = iaf(mc)-1+j1+(m1-1)*(j2-1)
            if(qmatches(ndf(1,it),ndf(1,jt),xyz,epsilon)) then
              nhitquad = nhitquad + 1
            endif
          end do
        end do
      else
        write(*,*) 'ERROR:  ia,ib out of range:',ia,ib
        stop 'qhitter,ia,ib'
      endif
c
      return
      end

c-----------------------------------------------------------------------

      logical function qmatches(ndf1,ndf2,xyz,epsilon)

      implicit real*8 (a-h,o-z)
      integer   ndf1(4),ndf2(4)
      dimension xyz(3,*)
c
c  If the first quad matches the second quad in any orientation
c  within epsilon distance at all four points then return .true.,
c  else .false.
c
      qmatches = .false.
      if(dist(xyz(1,ndf1(1)),xyz(1,ndf2(1))).le.epsilon) then
        if(dist(xyz(1,ndf1(2)),xyz(1,ndf2(2))).le.epsilon.and.
     .     dist(xyz(1,ndf1(3)),xyz(1,ndf2(3))).le.epsilon.and.
     .     dist(xyz(1,ndf1(4)),xyz(1,ndf2(4))).le.epsilon) then
          qmatches = .true.
        else if(dist(xyz(1,ndf1(2)),xyz(1,ndf2(4))).le.epsilon.and.
     .          dist(xyz(1,ndf1(3)),xyz(1,ndf2(3))).le.epsilon.and.
     .          dist(xyz(1,ndf1(4)),xyz(1,ndf2(2))).le.epsilon) then
          qmatches = .true.
        endif
      else if(dist(xyz(1,ndf1(1)),xyz(1,ndf2(2))).le.epsilon) then
        if(dist(xyz(1,ndf1(2)),xyz(1,ndf2(3))).le.epsilon.and.
     .     dist(xyz(1,ndf1(3)),xyz(1,ndf2(4))).le.epsilon.and.
     .     dist(xyz(1,ndf1(4)),xyz(1,ndf2(1))).le.epsilon) then
          qmatches = .true.
        else if(dist(xyz(1,ndf1(2)),xyz(1,ndf2(1))).le.epsilon.and.
     .          dist(xyz(1,ndf1(3)),xyz(1,ndf2(4))).le.epsilon.and.
     .          dist(xyz(1,ndf1(4)),xyz(1,ndf2(3))).le.epsilon) then
          qmatches = .true.
        endif
      else if(dist(xyz(1,ndf1(1)),xyz(1,ndf2(3))).le.epsilon) then
        if(dist(xyz(1,ndf1(2)),xyz(1,ndf2(4))).le.epsilon.and.
     .     dist(xyz(1,ndf1(3)),xyz(1,ndf2(1))).le.epsilon.and.
     .     dist(xyz(1,ndf1(4)),xyz(1,ndf2(2))).le.epsilon) then
          qmatches = .true.
        else if(dist(xyz(1,ndf1(2)),xyz(1,ndf2(2))).le.epsilon.and.
     .          dist(xyz(1,ndf1(3)),xyz(1,ndf2(1))).le.epsilon.and.
     .          dist(xyz(1,ndf1(4)),xyz(1,ndf2(4))).le.epsilon) then
          qmatches = .true.
        endif
      else if(dist(xyz(1,ndf1(1)),xyz(1,ndf2(4))).le.epsilon) then
        if(dist(xyz(1,ndf1(2)),xyz(1,ndf2(1))).le.epsilon.and.
     .     dist(xyz(1,ndf1(3)),xyz(1,ndf2(2))).le.epsilon.and.
     .     dist(xyz(1,ndf1(4)),xyz(1,ndf2(3))).le.epsilon) then
          qmatches = .true.
        else if(dist(xyz(1,ndf1(2)),xyz(1,ndf2(3))).le.epsilon.and.
     .          dist(xyz(1,ndf1(3)),xyz(1,ndf2(2))).le.epsilon.and.
     .          dist(xyz(1,ndf1(4)),xyz(1,ndf2(1))).le.epsilon) then
          qmatches = .true.
        endif
      endif
c
      return
      end

c-----------------------------------------------------------------------

      subroutine diddle(strings,nstrings,MSTRINGS)
c
c  Count and map the ways that two blocks can abut.
c
c  Scott D. Thomas, Raytheon STX Corporation
c  Contract NAS 2-98080/53, NASA Ames Research Center.
c  1-JUL-1998..4-JUL-1998.
c
c  The connection file for FLO107MB records connection, orientation,
c  and handedness information for each block.  It looks like this:
c
c 72
c    1 -5 1 2 3  13 1 2 3  61-1-2 3   2 1 2 3 -21 1 2 3   4 1 2 3   1    1
c    2 -5 1 2 3  14 1 2 3   1 1 2 3   3 1 2 3 -21 1 2 3   5 1 2 3   1    1
c    3 -5 1 2 3  15 1 2 3   2 1 2 3  -4 1 2 3 -21 1 2 3   6 1 2 3   1    1
c    4 -5 1 2 3  16 1 2 3  64-1-2 3   5 1 2 3   1 1 2 3   7 1 2 3   1    1
c[...]
c   71 59 1 2 3  -5 1 2 3  70 1 2 3  72 1 2 3  68 1 2 3  -5 1 2 3   1    1
c   72 60 1 2 3  -5 1 2 3  71 1 2 3  -4 1 2 3  69 1 2 3  -5 1 2 3   1    1
cWAKES
c   0
c
c  The connection file consists of the block number followed by four
c  integers for each of six faces, where the first integer, if
c  positive, is the abutting block number, and the next three integers
c  record the orientation of the block with respect to the block it
c  abuts.  If first integer is not positive then the orientation
c  integers are meaningless; it is a boundary condition other than
c  abutting flow-through.  The last two integers are the block's
c  handedness (1 for right-handed or -1 for left-handed) followed by an
c  unused value (set here to 1).
c
c  There are 6*4=24 ways to rotate a right-handed block:  six faces and
c  four rotations on the axis through each face.  If reflection is
c  allowed this makes 24*2=48 permutations of the vertices.  There are
c  48*48=2304 ways to position two adjacent blocks.  Some relative
c  orientations are repeated (overcounted) by this method of counting.
c  In fact, we are only interested in 2304/4=576 abutments.  Divide by
c  4 because each time both blocks are rotated by 90 degrees about the
c  axis normal to their abutting faces the relative alignment remains
c  the same. 
c
c  This program assumes that block 1 abuts block 2 on one face and it
c  goes through all the permutations.  It creates a one-to-one map 
c  between the following:
c
c  . first block's abutting face number
c  . second block's abutting face number
c  . face-to-face orientation
c  . first block's handedness
c    (6*6*8*2=576 choices)
c
c  and
c
c  . first block's abutting face number
c  . first block's orientation with respect to the second block
c  . first block's handedness
c    (6*48*2=576 choices)
c
c  The face numbers are 1 to 6 for jmin,jmax,kmin,kmax,lmin,lmax.
c
c  Face-to-face orientaton has to do with the way the abutting faces
c  align.  There are 8 ways for two faces to align, represented by
c  (ia,ib) where the first coordinate describes which of the second
c  abutting face's two basis vectors aligns with the first abutting
c  face's two basis vectors, and the second coordinate is for the
c  second of the first abutting face's two basis vectors.  The sign
c  indicates direction.  They are (ia,ib)=
c
c  (1,2) (-2,1) (-1,-2) (2,-1) (2,1) (-1,2) (-2,-1) (1,-2)
c
c  Block-to-block orientation has to do with the way the abutting
c  blocks align.  There are 24 ways that two right-handed blocks can
c  align, represented by (c1,c2,c3) where the first coordinate
c  describes which of the second abutting block's three basis vectors
c  aligns with the first abutting block's first basis vector, and so
c  on.  The sign indicates direction.  They are (c1,c2,c3)=
c
c  ( 1, 2, 3) ( 1, 3,-2) ( 1,-2,-3) ( 1,-3, 2) ( 2, 1,-3) ( 2, 3, 1)
c  ( 2,-1, 3) ( 2,-3,-1) ( 3, 1, 2) ( 3, 2,-1) ( 3,-1,-2) ( 3,-2, 1)
c  (-1, 2,-3) (-1, 3, 2) (-1,-2, 3) (-1,-3,-2) (-2, 1, 3) (-2, 3,-1)
c  (-2,-1,-3) (-2,-3, 1) (-3, 1,-2) (-3, 2, 1) (-3,-1, 2) (-3,-2,-1)
c
c  If left-handed block orientation is allowed then the number of
c  orientations is 3*2*8=48.
c
      character string*20, strings(MSTRINGS)*20
      dimension a(3,2,2,2),b(3,2,2,2)
      dimension r1(3,3),r2(3,3),r3(3,3)
      integer   permute(6,24),c1,c2,c3
      data permute /0,0,0,0,0,0, 1,0,0,0,0,0, 1,1,0,0,0,0, 1,1,1,0,0,0,
     .              2,0,0,0,0,0, 2,3,0,0,0,0, 2,3,3,0,0,0, 2,3,3,3,0,0,
     .              2,2,0,0,0,0, 2,2,1,0,0,0, 2,2,1,1,0,0, 2,2,1,1,1,0,
     .              2,2,2,0,0,0, 2,2,2,3,0,0, 2,2,2,3,3,0, 2,2,2,3,3,3,
     .              3,0,0,0,0,0, 3,2,0,0,0,0, 3,2,2,0,0,0, 3,2,2,2,0,0,
     .              3,3,3,0,0,0, 3,3,3,2,0,0, 3,3,3,2,2,0, 3,3,3,2,2,2/
c
c  r1 is a matrix for 90 deg rotation about the x axis.
c
      data (r1(1,j),j=1,3) / 1.0, 0.0, 0.0 /
      data (r1(2,j),j=1,3) / 0.0, 0.0,-1.0 /
      data (r1(3,j),j=1,3) / 0.0, 1.0, 0.0 /
c
c  r2 is a matrix for 90 deg rotation about the y axis.
c
      data (r2(1,j),j=1,3) / 0.0, 0.0,-1.0 /
      data (r2(2,j),j=1,3) / 0.0, 1.0, 0.0 /
      data (r2(3,j),j=1,3) / 1.0, 0.0, 0.0 /
c
c  r3 is a matrix for 90 deg rotation about the z axis.
c
      data (r3(1,j),j=1,3) / 0.0,-1.0, 0.0 /
      data (r3(2,j),j=1,3) / 1.0, 0.0, 0.0 /
      data (r3(3,j),j=1,3) / 0.0, 0.0, 1.0 /
c
c  Initialize strings.
c
      do j=1,MSTRINGS
        do i=1,20
          strings(j)(i:i) = ' '
        end do
      end do
      nstrings = 0
c
c  Outer loop through all the rotation group elements.
c
      do ip=1,48
        if(ip.le.24) then
          ipt = ip
          nhand = 1
        else
          ipt = ip - 24
          nhand = -1
        endif
        call unit(a,nhand)
        do i=1,6
          if(     permute(i,ipt).eq.1) then
            call leftmult(r1,a)
          else if(permute(i,ipt).eq.2) then
            call leftmult(r2,a)
          else if(permute(i,ipt).eq.3) then
            call leftmult(r3,a)
          endif
        end do
c
c  Inner loop through all the rotation group elements.
c
        do iq=1,48
          if(iq.le.24) then
            iqt = iq
            mhand = 1
          else
            iqt = iq - 24
            mhand = -1
          endif
          call unit(b,mhand)
          do i=1,6
            if(     permute(i,iqt).eq.1) then
              call leftmult(r1,b)
            else if(permute(i,iqt).eq.2) then
              call leftmult(r2,b)
            else if(permute(i,iqt).eq.3) then
              call leftmult(r3,b)
            endif
          end do
c
c  Compare the two cubes.
c
          call getcon(a,b,nf,mf,ia,ib,c1,c2,c3)
c
c  Write info that shows the relationship between:
c
c  . first block's abutting face number
c  . second block's abutting face number
c  . face-to-face orientation
c  . first block's handedness
c    (6*6*8*2=576 choices)
c
c  . first block's abutting face number
c  . first block's orientation with respect to the second block
c  . first block's handedness
c    (6*48*2=576 choices)
c
          write(string,1010) nf,mf,ia,ib,nhand, nf,c1,c2,c3,nhand
 1010     format(10i2)
          call addstrin(string,strings,nstrings,MSTRINGS)
        end do
      end do
c
      return
      end

c-----------------------------------------------------------------------

      subroutine addstrin(string,strings,nstrings,MSTRINGS)

      character string*20, strings(MSTRINGS)*20
c
      if(nstrings.le.0) then
        nstrings = 1
        strings(nstrings) = string
      else
        do i=1,nstrings
          if(string.eq.strings(i)) return 
        end do
        if(nstrings.ge.MSTRINGS) stop 'MSTRINGS'
        nstrings = nstrings + 1
        strings(nstrings) = string
      endif
c
      return
      end

c-----------------------------------------------------------------------

      function dot(a,b)

      dimension a(3),b(3)
      dot=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      return
      end

c-----------------------------------------------------------------------

      subroutine getcon(a,b,nf,mf,ia,ib,c1,c2,c3)

      dimension a(3,2,2,2),b(3,2,2,2)
      integer   nf,mf,ia,ib,c1,c2,c3
      dimension a1(3),a2(3),a3(3)
      dimension b1(3),b2(3),b3(3)
      data eps / 1.e-4 /
c
      nf = maxface(a,1)
      mf = minface(b,1)
c
      do m=1,3
        a1(m) = a(m,2,1,1)-a(m,1,1,1)
        a2(m) = a(m,1,2,1)-a(m,1,1,1)
        a3(m) = a(m,1,1,2)-a(m,1,1,1)
        b1(m) = b(m,2,1,1)-b(m,1,1,1)
        b2(m) = b(m,1,2,1)-b(m,1,1,1)
        b3(m) = b(m,1,1,2)-b(m,1,1,1)
      end do
c
      ia = 0
      ib = 0
c
      if(nf.eq.1.or.nf.eq.2) then
c a2,a3
        if(mf.eq.1.or.mf.eq.2) then
c b2,b3
          if(     dot(a2,b2).gt. eps) then
            ia =  1
          else if(dot(a2,b2).lt.-eps) then
            ia = -1
          else if(dot(a2,b3).gt. eps) then
            ia =  2
          else if(dot(a2,b3).lt.-eps) then
            ia = -2
          endif
          if(     dot(a3,b2).gt. eps) then
            ib =  1
          else if(dot(a3,b2).lt.-eps) then
            ib = -1
          else if(dot(a3,b3).gt. eps) then
            ib =  2
          else if(dot(a3,b3).lt.-eps) then
            ib = -2
          endif
        else if(mf.eq.3.or.mf.eq.4) then
c b1,b3
          if(     dot(a2,b1).gt. eps) then
            ia =  1
          else if(dot(a2,b1).lt.-eps) then
            ia = -1
          else if(dot(a2,b3).gt. eps) then
            ia =  2
          else if(dot(a2,b3).lt.-eps) then
            ia = -2
          endif
          if(     dot(a3,b1).gt. eps) then
            ib =  1
          else if(dot(a3,b1).lt.-eps) then
            ib = -1
          else if(dot(a3,b3).gt. eps) then
            ib =  2
          else if(dot(a3,b3).lt.-eps) then
            ib = -2
          endif
        else
c b1,b2
          if(     dot(a2,b1).gt. eps) then
            ia =  1
          else if(dot(a2,b1).lt.-eps) then
            ia = -1
          else if(dot(a2,b2).gt. eps) then
            ia =  2
          else if(dot(a2,b2).lt.-eps) then
            ia = -2
          endif
          if(     dot(a3,b1).gt. eps) then
            ib =  1
          else if(dot(a3,b1).lt.-eps) then
            ib = -1
          else if(dot(a3,b2).gt. eps) then
            ib =  2
          else if(dot(a3,b2).lt.-eps) then
            ib = -2
          endif
        endif
      else if(nf.eq.3.or.nf.eq.4) then
c a1,a3
        if(mf.eq.1.or.mf.eq.2) then
c b2,b3
          if(     dot(a1,b2).gt. eps) then
            ia =  1
          else if(dot(a1,b2).lt.-eps) then
            ia = -1
          else if(dot(a1,b3).gt. eps) then
            ia =  2
          else if(dot(a1,b3).lt.-eps) then
            ia = -2
          endif
          if(     dot(a3,b2).gt. eps) then
            ib =  1
          else if(dot(a3,b2).lt.-eps) then
            ib = -1
          else if(dot(a3,b3).gt. eps) then
            ib =  2
          else if(dot(a3,b3).lt.-eps) then
            ib = -2
          endif
        else if(mf.eq.3.or.mf.eq.4) then
c b1,b3
          if(     dot(a1,b1).gt. eps) then
            ia =  1
          else if(dot(a1,b1).lt.-eps) then
            ia = -1
          else if(dot(a1,b3).gt. eps) then
            ia =  2
          else if(dot(a1,b3).lt.-eps) then
            ia = -2
          endif
          if(     dot(a3,b1).gt. eps) then
            ib =  1
          else if(dot(a3,b1).lt.-eps) then
            ib = -1
          else if(dot(a3,b3).gt. eps) then
            ib =  2
          else if(dot(a3,b3).lt.-eps) then
            ib = -2
          endif
        else
c b1,b2
          if(     dot(a1,b1).gt. eps) then
            ia =  1
          else if(dot(a1,b1).lt.-eps) then
            ia = -1
          else if(dot(a1,b2).gt. eps) then
            ia =  2
          else if(dot(a1,b2).lt.-eps) then
            ia = -2
          endif
          if(     dot(a3,b1).gt. eps) then
            ib =  1
          else if(dot(a3,b1).lt.-eps) then
            ib = -1
          else if(dot(a3,b2).gt. eps) then
            ib =  2
          else if(dot(a3,b2).lt.-eps) then
            ib = -2
          endif
        endif
      else
c a1,a2
        if(mf.eq.1.or.mf.eq.2) then
c b2,b3
          if(     dot(a1,b2).gt. eps) then
            ia =  1
          else if(dot(a1,b2).lt.-eps) then
            ia = -1
          else if(dot(a1,b3).gt. eps) then
            ia =  2
          else if(dot(a1,b3).lt.-eps) then
            ia = -2
          endif
          if(     dot(a2,b2).gt. eps) then
            ib =  1
          else if(dot(a2,b2).lt.-eps) then
            ib = -1
          else if(dot(a2,b3).gt. eps) then
            ib =  2
          else if(dot(a2,b3).lt.-eps) then
            ib = -2
          endif
        else if(mf.eq.3.or.mf.eq.4) then
c b1,b3
          if(     dot(a1,b1).gt. eps) then
            ia =  1
          else if(dot(a1,b1).lt.-eps) then
            ia = -1
          else if(dot(a1,b3).gt. eps) then
            ia =  2
          else if(dot(a1,b3).lt.-eps) then
            ia = -2
          endif
          if(     dot(a2,b1).gt. eps) then
            ib =  1
          else if(dot(a2,b1).lt.-eps) then
            ib = -1
          else if(dot(a2,b3).gt. eps) then
            ib =  2
          else if(dot(a2,b3).lt.-eps) then
            ib = -2
          endif
        else
c b1,b2
          if(     dot(a1,b1).gt. eps) then
            ia =  1
          else if(dot(a1,b1).lt.-eps) then
            ia = -1
          else if(dot(a1,b2).gt. eps) then
            ia =  2
          else if(dot(a1,b2).lt.-eps) then
            ia = -2
          endif
          if(     dot(a2,b1).gt. eps) then
            ib =  1
          else if(dot(a2,b1).lt.-eps) then
            ib = -1
          else if(dot(a2,b2).gt. eps) then
            ib =  2
          else if(dot(a2,b2).lt.-eps) then
            ib = -2
          endif
        endif
      endif
      if(ia.eq.0.or.ib.eq.0) stop 'ia.or.ib'
c
      c1 = 0
      c2 = 0
      c3 = 0
c
      if(     dot(a1,b1).gt. eps) then
        c1 =  1
      else if(dot(a1,b1).lt.-eps) then
        c1 = -1
      else if(dot(a1,b2).gt. eps) then
        c1 =  2
      else if(dot(a1,b2).lt.-eps) then
        c1 = -2
      else if(dot(a1,b3).gt. eps) then
        c1 =  3
      else if(dot(a1,b3).lt.-eps) then
        c1 = -3
      endif
      if(     dot(a2,b1).gt. eps) then
        c2 =  1
      else if(dot(a2,b1).lt.-eps) then
        c2 = -1
      else if(dot(a2,b2).gt. eps) then
        c2 =  2
      else if(dot(a2,b2).lt.-eps) then
        c2 = -2
      else if(dot(a2,b3).gt. eps) then
        c2 =  3
      else if(dot(a2,b3).lt.-eps) then
        c2 = -3
      endif
      if(     dot(a3,b1).gt. eps) then
        c3 =  1
      else if(dot(a3,b1).lt.-eps) then
        c3 = -1
      else if(dot(a3,b2).gt. eps) then
        c3 =  2
      else if(dot(a3,b2).lt.-eps) then
        c3 = -2
      else if(dot(a3,b3).gt. eps) then
        c3 =  3
      else if(dot(a3,b3).lt.-eps) then
        c3 = -3
      endif
      if(c1.eq.0.or.c2.eq.0.or.c3.eq.0) stop 'c[123]'
c
      return
      end

c-----------------------------------------------------------------------

      subroutine leftmult(r,a)

      dimension a(3,2,2,2),r(3,3),ar(3,2,2,2)
      do l=1,2
        do k=1,2
          do j=1,2
            do n=1,3
              sum = 0
              do m=1,3
                sum = sum + r(n,m)*a(m,j,k,l)
              end do
              ar(n,j,k,l) = sum
            end do
          end do
        end do
      end do
      do l=1,2
        do k=1,2
          do j=1,2
            do n=1,3
              a(n,j,k,l) = ar(n,j,k,l)
            end do
          end do
        end do
      end do
      return
      end

c-----------------------------------------------------------------------

      function maxface(a,m)

      dimension a(3,2,2,2)
      if(     a(m,1,1,1).gt.0.0.and.a(m,1,2,1).gt.0.0.and.
     .        a(m,1,2,2).gt.0.0.and.a(m,1,1,2).gt.0.0) then
        maxface = 1
      else if(a(m,2,1,1).gt.0.0.and.a(m,2,2,1).gt.0.0.and.
     .        a(m,2,2,2).gt.0.0.and.a(m,2,1,2).gt.0.0) then
        maxface = 2
      else if(a(m,1,1,1).gt.0.0.and.a(m,2,1,1).gt.0.0.and.
     .        a(m,2,1,2).gt.0.0.and.a(m,1,1,2).gt.0.0) then
        maxface = 3
      else if(a(m,1,2,1).gt.0.0.and.a(m,2,2,1).gt.0.0.and.
     .        a(m,2,2,2).gt.0.0.and.a(m,1,2,2).gt.0.0) then
        maxface = 4
      else if(a(m,1,1,1).gt.0.0.and.a(m,2,1,1).gt.0.0.and.
     .        a(m,2,2,1).gt.0.0.and.a(m,1,2,1).gt.0.0) then
        maxface = 5
      else if(a(m,1,1,2).gt.0.0.and.a(m,2,1,2).gt.0.0.and.
     .        a(m,2,2,2).gt.0.0.and.a(m,1,2,2).gt.0.0) then
        maxface = 6
      else
        maxface = 0
      endif
      return
      end

c-----------------------------------------------------------------------

      function minface(a,m)

      dimension a(3,2,2,2)
      if(     a(m,1,1,1).lt.0.0.and.a(m,1,2,1).lt.0.0.and.
     .        a(m,1,2,2).lt.0.0.and.a(m,1,1,2).lt.0.0) then
        minface = 1
      else if(a(m,2,1,1).lt.0.0.and.a(m,2,2,1).lt.0.0.and.
     .        a(m,2,2,2).lt.0.0.and.a(m,2,1,2).lt.0.0) then
        minface = 2
      else if(a(m,1,1,1).lt.0.0.and.a(m,2,1,1).lt.0.0.and.
     .        a(m,2,1,2).lt.0.0.and.a(m,1,1,2).lt.0.0) then
        minface = 3
      else if(a(m,1,2,1).lt.0.0.and.a(m,2,2,1).lt.0.0.and.
     .        a(m,2,2,2).lt.0.0.and.a(m,1,2,2).lt.0.0) then
        minface = 4
      else if(a(m,1,1,1).lt.0.0.and.a(m,2,1,1).lt.0.0.and.
     .        a(m,2,2,1).lt.0.0.and.a(m,1,2,1).lt.0.0) then
        minface = 5
      else if(a(m,1,1,2).lt.0.0.and.a(m,2,1,2).lt.0.0.and.
     .        a(m,2,2,2).lt.0.0.and.a(m,1,2,2).lt.0.0) then
        minface = 6
      else
        minface = 0
      endif
      return
      end

c-----------------------------------------------------------------------

      subroutine unit(a,ihand)

      dimension a(3,2,2,2)
      do l=1,2
        do k=1,2
          do j=1,2
            a(1,j,k,l) =  -1+2*(j-1)
            a(2,j,k,l) =  -1+2*(k-1)
            a(3,j,k,l) = (-1+2*(l-1))*isign(1,ihand)
          end do
        end do
      end do
      return
      end
