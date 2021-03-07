!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program gu

!  Grid utility to manipulate multiblock PLOT3D grids or function files.

!  05/01/13  Ryan McDaniel   Original release (mostly FORTRAN 77).
!  May 2013  David Saunders  Fortran 90 version.  Swap x/y/z option added.
!  09/09/13    "      "      Ryan found mirrored block handedness was forgotten.
!                            Also: the header record of unformatted blocks
!                            extracted as a chunk was wrong (not one record).
!  07/08/16    "      "      Ryan provided a revised subroutine scale with an
!                            option to vary the scaling with x.
!  12/12/17    "      "      Ryan found that the handedness wasn't maintained
!                            for the "swapijk" option.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Variables:

   integer        :: idone, ios, mmax, qmax
   logical        :: grid
   character (1)  :: dowhat, gorf, iform, iform1, iform2, oform
   character (64) :: name, name2

   integer, allocatable, dimension (:)       :: imax, jmax, kmax
   integer, allocatable, dimension (:)       :: imaxtemp, jmaxtemp, kmaxtemp
   real,    allocatable, dimension (:,:,:)   :: x, y, z
   real,    allocatable, dimension (:,:,:,:) :: q

!  Execution:

!  Prompt for the initial grid or function file and determine its format:

   call get_input_name_and_form &
      (1, 'PLOT3D multiblock grid or function file: ', name, iform, ios)
   if (ios /= 0) go to 99

   do ! Until valid
      write (*, '(a)', advance='no') 'Grid or function file? [g|f] '
      read  (*, *, iostat=ios) gorf
      if (ios /= 0) go to 99
      grid = gorf == 'g'
      if (grid .or. gorf == 'f') exit
   end do

   write (*, '(/, (3x, a))') &
      'reorder blocks (r)',                      &
      'extract blocks (e)',                      &
      'combine grids (b)',                       &
      'convert formatted <--> unformatted (c)',  &
      'sequence (s)',                            &
      'count points (n)',                        &
      'unit conversion (u)',                     &
      'add zeros (pad 2d for Gridgen) (a)',      &
      'remove zeros (for dplr2d) (z)',           &
      'extract surface k = 1 or kmax (x)',       &
      'swap i/j/k (i)',                          &
      'swap x/y/z (y)',                          &
      'fix a surface file to have kmax = 1 (f)', &
      'mirror grid (m)',                         &
      'probe for max and min (p)',               &
      'scale outer boundary (o)',                &
      'extrude 2d plane to create 3d volume (v)'

10 write (*, '(/, a)', advance='no') 'Choice: '
   read  (*, *, iostat=ios) dowhat
   if (ios /= 0) go to 99

!  Some options don't write a new file; trap a bad choice too:

   select case (dowhat)
      case ('a', 'b', 'e', 'f', 'i', 'm', 'o', 'r', 's', 'u', 'v', 'x', 'y', 'z')
         do ! Until valid
            write (*,'(a)', advance='no') 'Formatted|unformatted output? [f|u] '
            read  (*, *, iostat=ios) oform
            if (ios /= 0) go to 99  ! E.g., ^D means quit
            if (oform == 'f' .or. oform == 'u') exit
         end do
      case ('c', 'n', 'p')
         ! No need for output format prompt
      case default
         print *, 'Entry not recognized.  Try again.'
         go to 10
   end select

!  The combine-files case is the only one with a second input file:

   if (dowhat == 'b') then
      iform1 = iform
      call get_input_name_and_form &
         (1, '2nd file name: ', name2, iform2, ios)
      if (ios /= 0) go to 99
   end if

!  Allocate work-space for all cases:

   call makespace

!  Main branch:

   select case (dowhat)

      case ('a') ! Add z = 0
         call pad

      case ('b') ! Combine grids
         call combine

      case ('c') ! Convert format
         call convert

      case ('e') ! Extract 1 or more blocks, individually or as 1 file
         call extract

      case ('f') ! Fix a surface file so that kmax = 1
         call fix

      case ('i') ! Swap i/j/k
         call swapijk

      case ('m') ! Mirror
         call mirror

      case ('n') ! Count points
         call countnum

      case ('o') ! Scale outer boundary
         call scale

      case ('p') ! Probe for max/min
         call probe

      case ('r') ! Reorder blocks
         call reorder

      case ('s') ! Sequence
         call sequence

      case ('u') ! Unit conversion
         call unitconv

      case ('v') ! Extrude 2d --> 3d volume
         call extrude

      case ('x') ! Extract surface k = 1 or kmax
         call surf

      case ('y') ! Swap x/y/z
         call swapxyz

      case ('z') ! Remove z = 0
         call zeros

   end select

   if (dowhat /= 'p') write (*, '(/, 2a)') 'Wrote ', name

99 continue  ! STOP behavior can be system-dependent

!  Internal procedures for program gu:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_input_name_and_form (lun, prompt, name, iform, ios)

!     Determine an existing input file name and its format.  A complication is
!     that an open with a valid name can succeed with an invalid format, so we
!     try to read the first record, which is assumed to be the number of blocks
!     in a multiblock PLOT3D grid or function file.
!     Allow for repeated tries, and signal failure if no good file is found.
!     Any opened file is closed before return.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)        :: lun    ! Logical unit to use
      character (*), intent (in)  :: prompt ! Allows use on more than one file
      character (*), intent (out) :: name   ! File name found (valid if ier = 0)
      character (*), intent (out) :: iform  ! 'f' or 'u' (valid if ier = 0)
      integer, intent (out)       :: ios    ! 0 for success; 1 for failure

!     Local variables:

      integer :: nblocks

!     Execution:

      do ! Until success or quit
         write (*, '(a)', advance='no') prompt  ! E.g., 'Second input file: '
         read  (*, *, iostat=ios) name
         if (ios <  0) go to 99  ! ^D = bail out
         open (lun, file=name, status='old', iostat=ios)
         if (ios == 0) then  ! Opening alone is not sufficient proof of format
            read (lun, *, iostat=ios) nblocks
            if (ios == 0) then
               iform = 'f'
               exit
            else
               close (lun)
            end if
         end if
         open (lun, file=name, status='old', form='unformatted', iostat=ios)
         if (ios == 0) then
            read (lun, iostat=ios) nblocks
            if (ios == 0) then
               iform = 'u'
               exit
            else
               close (lun)
            end if
         end if
         write (*, '(2a)') 'Problem opening file ', trim (name)
      end do
      close (lun)

   99 return

      end subroutine get_input_name_and_form

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine makespace  ! Allocate work-space for all cases.

!     Once the dimensions of all blocks are known, one block of largest possible
!     size is allocated and used for all blocks.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer :: icount, idim, jdim, kdim, m, mf, mmax2
      character (64) :: nname

!     Execution:

      idim   = 0
      jdim   = 0
      kdim   = 0
      mmax2  = 0
      icount = 0
      nname  = name

321   continue

      if (grid) then
         if (iform == 'f') then
            open (1, file=nname, status='unknown')
            read (1,*) mmax
            mmax2 = max (mmax, mmax2)
            allocate (imaxtemp(mmax2), jmaxtemp(mmax2), kmaxtemp(mmax2))
            if (dowhat /= 'a') then
               read (1,*) (imaxtemp(m), jmaxtemp(m), kmaxtemp(m), m = 1, mmax)
            else  ! We'll be adding z = 0
               read (1,*) (imaxtemp(m), jmaxtemp(m), m = 1, mmax)
            end if
         else
            open (1, file=nname, status='unknown', form='unformatted')
            read (1) mmax
            mmax2 = max (mmax, mmax2)
            allocate (imaxtemp(mmax2),jmaxtemp(mmax2),kmaxtemp(mmax2))
            if (dowhat /= 'a') then
               read (1) (imaxtemp(m), jmaxtemp(m), kmaxtemp(m), m = 1, mmax)
            else
               read (1) (imaxtemp(m), jmaxtemp(m), m = 1, mmax)
            end if
         end if
         do m = 1, mmax
            if (dowhat == 'a') kmaxtemp(m) = 1
            if (imaxtemp(m) > idim) idim = imaxtemp(m)
            if (jmaxtemp(m) > jdim) jdim = jmaxtemp(m)
            if (kmaxtemp(m) > kdim) kdim = kmaxtemp(m)
         end do
         allocate (x(idim,jdim,kdim), y(idim,jdim,kdim), z(idim,jdim,kdim))
      else  ! Function file
         if (iform == 'f') then
            open (1, file=nname, status='unknown')
            read (1,*) mmax
            mmax2 = max (mmax, mmax2)
            allocate (imaxtemp(mmax2), jmaxtemp(mmax2), kmaxtemp(mmax2))
            read (1,*) (imaxtemp(m), jmaxtemp(m), kmaxtemp(m), qmax, m = 1,mmax)
         else
            open (1, file=nname, status='unknown', form='unformatted')
            read (1) mmax
            mmax2 = max (mmax, mmax2)
            allocate (imaxtemp(mmax2), jmaxtemp(mmax2), kmaxtemp(mmax2))
            read (1) (imaxtemp(m), jmaxtemp(m), kmaxtemp(m), qmax, m = 1, mmax)
         end if
         do m = 1, mmax
            if (dowhat == 'a') kmaxtemp(m) = 1
            if (imaxtemp(m) > idim) idim = imaxtemp(m)
            if (jmaxtemp(m) > jdim) jdim = jmaxtemp(m)
            if (kmaxtemp(m) > kdim) kdim = kmaxtemp(m)
         end do
         allocate (q(idim,jdim,kdim,qmax))
      end if

      close (1)

!     Combine's file 2 may have a bigger block than the biggest in file 1:

      icount = icount + 1
      if (dowhat == 'b') then
         if (icount == 1) then
            deallocate (imaxtemp, jmaxtemp, kmaxtemp)
            if (grid) then
               deallocate (x, y, z)
            else
               deallocate (q)
            end if
            nname = name2
            iform = iform2
            go to 321
         else  ! icount = 2
            iform = iform1
         end if
      end if
      deallocate (imaxtemp, jmaxtemp, kmaxtemp)

      mf = mmax2
      if (dowhat == 'm' .or. dowhat == 'b') mf = 2*mf
      allocate (imax(mf), jmax(mf), kmax(mf))

      end subroutine makespace

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine iohead (ioro)  ! Read or write a multiblock PLOT3D file header.

!     The file is opened here and read on unit 1 or written on unit 2.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, intent (in) :: ioro  ! 1 => read
      integer :: m

      if (ioro == 1) then
         if (iform == 'f') then
            open (1, file=name, status='unknown')
            read (1, *) mmax
            if (grid) then
               read (1, *) (imax(m), jmax(m), kmax(m), m = 1, mmax)
            else
               read (1, *) (imax(m), jmax(m), kmax(m), qmax, m = 1, mmax)
            end if
         else
            open (1, file=name, status='unknown', form='unformatted')
            read (1) mmax
            if (grid) then
               read (1) (imax(m), jmax(m), kmax(m), m = 1, mmax)
            else
               read (1) (imax(m), jmax(m), kmax(m), qmax, m = 1, mmax)
            end if
         end if
      else
         if (oform == 'f') then
            open  (2, file=name, status='unknown')
            write (2, '(i4)') mmax
            if (grid) then
               write (2, '(3i5)') (imax(m), jmax(m), kmax(m), m = 1, mmax)
            else
               write (2, '(4i5)') (imax(m), jmax(m), kmax(m), qmax, m = 1 ,mmax)
            end if
         else
            name = trim (name) // 'u'
            open  (2, file=name, status='unknown', form='unformatted')
            write (2) mmax
            if (grid) then
               write (2) (imax(m), jmax(m), kmax(m), m = 1, mmax)
            else
               write (2) (imax(m), jmax(m), kmax(m), qmax, m = 1, mmax)
            end if
         end if
      end if

      end subroutine iohead

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine iogrid (ioro, m)

!     Read or write block m of the indicated grid or function file.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, intent (in) :: ioro, m  ! ioro = 1 => read; m = block #
      integer :: i, j, k, l, ni, nj, nk

      ni = imax(m)
      nj = jmax(m)
      nk = kmax(m)

      if (ioro == 1) then
         write (*, '(a, i5)') '   Reading block', m
         if (iform == 'f') then
            if (grid) then
               read (1, *) (((x(i,j,k),i=1,ni),j=1,nj),k=1,nk), &
                           (((y(i,j,k),i=1,ni),j=1,nj),k=1,nk), &
                           (((z(i,j,k),i=1,ni),j=1,nj),k=1,nk)
            else
               read (1, *) ((((q(i,j,k,l),i=1,ni),j=1,nj),k=1,nk),l=1,qmax)
            end if
         else
            if (grid) then
               read (1) (((x(i,j,k),i=1,ni),j=1,nj),k=1,nk), &
                        (((y(i,j,k),i=1,ni),j=1,nj),k=1,nk), &
                        (((z(i,j,k),i=1,ni),j=1,nj),k=1,nk)
            else
               read (1) ((((q(i,j,k,l),i=1,ni),j=1,nj),k=1,nk),l=1,qmax)
            end if
         end if
      else
         write (*, '(a, i5)') '   Writing block', m
         if (oform == 'f') then
            if (grid) then
               write (2, '(4es19.11)') &
                  (((x(i,j,k),i=1,ni),j=1,nj),k=1,nk), &
                  (((y(i,j,k),i=1,ni),j=1,nj),k=1,nk), &
                  (((z(i,j,k),i=1,ni),j=1,nj),k=1,nk)
            else
               write (2, '(4es19.11)') &
                  ((((q(i,j,k,l),i=1,ni),j=1,nj),k=1,nk),l=1,qmax)
            end if
         else
            if (grid) then
               write (2) (((x(i,j,k),i=1,ni),j=1,nj),k=1,nk), &
                         (((y(i,j,k),i=1,ni),j=1,nj),k=1,nk), &
                         (((z(i,j,k),i=1,ni),j=1,nj),k=1,nk)
            else
               write (2) ((((q(i,j,k,l),i=1,ni),j=1,nj),k=1,nk),l=1,qmax)
            end if
         end if
      end if

      end subroutine iogrid

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine countnum  ! Count the numbers of points in each block.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: m, icount, itotal

      call iohead (1)
      name = 'head.out'
      open  (654, file=name, status='unknown')
      write (654, '(i4)') mmax
      icount = 0
      itotal = 0
      do m = 1, mmax
         icount = imax(m)*jmax(m)*kmax(m)
         itotal = itotal + icount
         if (grid) then
            write (*, '(i4, 3i5, a, i8)') m,imax(m),jmax(m),kmax(m),':', icount
            write (654, '(3i5)') imax(m),jmax(m),kmax(m)
         else
            write (*, '(i4, 3i5, i3, a, i8)') &
               m,imax(m),jmax(m),kmax(m),qmax,':', icount
            write (654, '(3i5, i3)') imax(m),jmax(m),kmax(m),qmax
         end if
      end do
      write (*, '(a, i9)') 'Total number of grid points:', itotal

      end subroutine countnum

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine convert  ! .. between formatted and unformatted.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: m

      call iohead (1)
      if (iform == 'f') then
         name  = 'unformatted.p3d'
         oform = 'u'
      else
         name  = 'formatted.p3d'
         oform = 'f'
      end if
      call iohead (2)

      do m = 1, mmax
         call iogrid (1, m)
         call iogrid (2, m)
      end do

      end subroutine convert

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine probe  ! ... for max. & min. coordinates or function values.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, parameter :: big = 9.e9
      integer :: i, j, k, l, m
      integer :: ixmax,ixmin,iymax,iymin,izmax,izmin,iqmax,iqmin
      integer :: jxmax,jxmin,jymax,jymin,jzmax,jzmin,jqmax,jqmin
      integer :: kxmax,kxmin,kymax,kymin,kzmax,kzmin,kqmax,kqmin
      real    ::  xmax, xmin, ymax, ymin, zmax, zmin
      real, allocatable, dimension (:) :: maxq, minq

      call iohead (1)

      if (.not. grid) allocate (maxq(qmax), minq(qmax))

      do m = 1, mmax
         if (grid) then
            xmin =  big
            ymin =  big
            zmin =  big
            xmax = -big
            ymax = -big
            zmax = -big
         else
            minq(1:qmax) =  big
            maxq(1:qmax) = -big
         end if
         call iogrid (1, m)
         print *
         if (grid) then
            do k=1,kmax(m)
               do j=1,jmax(m)
                  do i=1,imax(m)
                     if (x(i,j,k) > xmax) then
                        xmax=x(i,j,k)
                        ixmax=i
                        jxmax=j
                        kxmax=k
                     else if (x(i,j,k) < xmin) then
                        xmin=x(i,j,k)
                        ixmin=i
                        jxmin=j
                        kxmin=k
                     end if
                     if (y(i,j,k) > ymax) then
                        ymax=y(i,j,k)
                        iymax=i
                        jymax=j
                        kymax=k
                     else if (y(i,j,k) < ymin) then
                        ymin=y(i,j,k)
                        iymin=i
                        jymin=j
                        kymin=k
                     end if
                     if (z(i,j,k) > zmax) then
                        zmax=z(i,j,k)
                        izmax=i
                        jzmax=j
                        kzmax=k
                     else if (z(i,j,k) < zmin) then
                        zmin=z(i,j,k)
                        izmin=i
                        jzmin=j
                        kzmin=k
                     end if
                  end do
               end do
            end do
            write (*, '(a, es19.11, a, 3i6)') &
               'xmin =', xmin, ' @', ixmin,jxmin,kxmin, &
               'xmax =', xmax, ' @', ixmax,jxmax,kxmax, &
               'ymin =', ymin, ' @', iymin,jymin,kymin, &
               'ymax =', ymax, ' @', iymax,jymax,kymax, &
               'zmin =', zmin, ' @', izmin,jzmin,kzmin, &
               'zmax =', zmax, ' @', izmax,jzmax,kzmax
            print *
         else  ! Function file
            do l=1,qmax
               do k=1,kmax(m)
                  do j=1,jmax(m)
                     do i=1,imax(m)
                        if (q(i,j,k,l) > maxq(l)) then
                           maxq(l)=q(i,j,k,l)
                           iqmax=i
                           jqmax=j
                           kqmax=k
                        else if (q(i,j,k,l) < minq(l)) then
                           minq(l)=q(i,j,k,l)
                           iqmin=i
                           jqmin=j
                           kqmin=k
                        end if
                     end do
                  end do
               end do
               write (*, '(i4, 2x, a, es19.11, a, 3i6)') &
                  l, 'qmin =', minq(l), ' @', iqmin,jqmin,kqmin, &
                  l, 'qmax =', maxq(l), ' @', iqmax,jqmax,kqmax
            end do
         end if
      end do

      if (.not. grid) deallocate (maxq, minq)

      end subroutine probe

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mirror  ! ... in x|y|z (grid|function file; save both halves,
                         ! and preserve reflected-block handedness (reverse i).
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, j, k, l, m, n, ni, nip1, niby2, nj, nk, num
      real    :: t
      character (1)  :: dir
      character (64) :: name1, name2
      integer, allocatable, dimension (:) :: iqrev, s

      if (grid) then
         write (*, '(a)', advance='no') 'Mirror in x, y, or z direction? '
         read  (*, *) dir
      else
         allocate (iqrev(qmax))
         write (*, '(a, i3, a)', advance='no') &
            'How many of the', qmax, ' function variables need sign changes? '
         read  (*, *) num
         if (num > 0) then
            write (*, '(a, i2, a)', advance='no') &
               'Enter', num, ' variable numbers: '
            read  (*, *) iqrev(1:num)
         end if
      end if

      print *
      call iohead (1)

      do m = 1, mmax
         imax(mmax+m) = imax(m)
         jmax(mmax+m) = jmax(m)
         kmax(mmax+m) = kmax(m)
      end do

      mmax  = mmax*2
      name1 = name
      name2 = 'mirror.p3d'
      name  = name2

      call iohead (2)

      do m = 1, mmax/2
         call iogrid (1, m)
         call iogrid (2, m)
      end do
      close (1)

!     Rereading the file is preferable to storing all blocks for reflecting.
!     Also, reversing the i indices to preserve handedness is best done in place
!     to keep the writing efficient.

      if (.not. grid) allocate (s(qmax))
      name = name1
      call iohead (1)

      do m = 1, mmax
         call iogrid (1, m)
         ni = imax(m);  nip1 = ni + 1;  niby2 = nip1 / 2
         nj = jmax(m)
         nk = kmax(m)
         if (grid) then
            select case (dir)
               case ('x')
                  x(1:ni,1:nj,1:nk) = -x(1:ni,1:nj,1:nk)
               case ('y')
                  y(1:ni,1:nj,1:nk) = -y(1:ni,1:nj,1:nk)
               case ('z')
                  z(1:ni,1:nj,1:nk) = -z(1:ni,1:nj,1:nk)
            end select
            do k = 1, nk
               do j = 1, nj
                  do i = 1, niby2
                     l = nip1 - i
                     t = x(i,j,k)
                     x(i,j,k) = x(l,j,k)
                     x(l,j,k) = t
                     t = y(i,j,k)
                     y(i,j,k) = y(l,j,k)
                     y(l,j,k) = t
                     t = z(i,j,k)
                     z(i,j,k) = z(l,j,k)
                     z(l,j,k) = t
                  end do
               end do
            end do
         else
            do l = 1, num
               n = iqrev(l)
               q(1:ni,1:nj,1:nk,n) = -q(1:ni,1:nj,1:nk,n)
            end do
            do k = 1, nk
               do j = 1, nj
                  do i = 1, niby2
                     l = nip1 - i
                     do n = 1, qmax
                        t = q(i,j,k,n)
                        q(i,j,k,n) = q(l,j,k,n)
                        q(l,j,k,n) = t
                     end do
!!                   s(:) = q(i,j,k,:)        ! Why does this misbehave?
!!                   q(i,j,k,:) = q(l,j,k,:)
!!                   q(l,j,k,:) = s(:)
                  end do
               end do
            end do
         end if
         call iogrid (2, m)
      end do

      if (.not. grid) deallocate (iqrev, s)

      name = 'mirror.p3d'
      if (oform == 'u') name = 'mirror.p3du'

      end subroutine mirror

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine pad  ! ... 2d grid with with z = 0.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, j, k, m

      if (iform == 'f') then
          open (1, file=name, status='old')
          read (1, *) mmax
          read (1, *) (imax(m), jmax(m), m = 1, mmax)
      else
          open (1, file=name, status='old', form='unformatted')
          read (1) mmax
          read (1) (imax(m), jmax(m), m = 1, mmax)
      end if

      do m = 1, mmax
         kmax(m) = 1
         write (*, '(2i5)') imax(m), jmax(m)
      end do

      name = 'padded.p3d'
      call iohead (2)

      do m = 1, mmax
         if (iform == 'f') then
            read (1, *) &
               (((x(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)), &
               (((y(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m))
         else
            read (1) &
               (((x(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)), &
               (((y(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m))
         end if
         call iogrid (2, m)
      end do

      end subroutine pad

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine zeros  ! Remove z = 0 from what is really a 2-space grid.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, j, m, ni, nj

      call iohead (1)
      if (oform == 'f') then
         name='2d.p3d'
         open  (2, file=name, status='unknown')
         write (2, '(i4)') mmax
         write (2, '(2i5)') (imax(m), jmax(m), m = 1, mmax)
      else
         name='2d.p3du'
         open  (2, file=name, status='unknown', form='unformatted')
         write (2) mmax
         write (2) (imax(m), jmax(m), m = 1, mmax)
      end if

      do m = 1, mmax
         call iogrid (1, m)
         ni = imax(m)
         nj = jmax(m)
         write (*, '(a, i5)') '   Writing block', m
         if (oform == 'f') then
            write (2, '(4es19.11)') &
               ((x(i,j,1),i=1,ni),j=1,nj), ((y(i,j,1),i=1,ni),j=1,nj)
         else
            write (2) ((x(i,j,1),i=1,ni),j=1,nj), ((y(i,j,1),i=1,ni),j=1,nj)
         end if
      end do

      end subroutine zeros

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine combine  ! Combine two files into one.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: m, mmax1, mmax2, mmaxboth
      character (64) :: name1, name3
      integer, allocatable, dimension (:) :: imaxtemp, jmaxtemp, kmaxtemp

      name1 = name
      name  = name2
      iform = iform2
      call iohead (1)
      close (1)

      mmax2 = mmax
      name  = name1
      iform = iform1
      call iohead (1)
      close (1)

      mmax1 = mmax
      mmaxboth = mmax1 + mmax2
      allocate (imaxtemp(mmaxboth), jmaxtemp(mmaxboth), kmaxtemp(mmaxboth))

      do m = 1,mmax
         imaxtemp(m) = imax(m)
         jmaxtemp(m) = jmax(m)
         kmaxtemp(m) = kmax(m)
      end do

      name  = name2
      iform = iform2
      call iohead (1)
      close (1)

      do m = 1, mmax
         imaxtemp(m+mmax1) = imax(m)
         jmaxtemp(m+mmax1) = jmax(m)
         kmaxtemp(m+mmax1) = kmax(m)
      end do

      if (oform == 'f') then
         name = 'merged.p3d'
         open  (2, file=name, status='unknown')
         write (2, '(i4)') mmaxboth
         if (grid) then
            write (2, '(3i5)') &
               (imaxtemp(m), jmaxtemp(m), kmaxtemp(m), m = 1, mmaxboth)
         else
            write (2, '(3i5, i4)') &
               (imaxtemp(m), jmaxtemp(m), kmaxtemp(m), qmax, m = 1, mmaxboth)
         end if
      else
         name = 'merged.p3du'
         open  (2, file=name, form='unformatted', status='unknown')
         write (2) mmaxboth
         if (grid) then
            write (2) (imaxtemp(m), jmaxtemp(m), kmaxtemp(m), m = 1, mmaxboth)
         else
            write (2) (imaxtemp(m),jmaxtemp(m),kmaxtemp(m),qmax,m = 1, mmaxboth)
         end if
      end if

      write (*, '(2a)') 'Reading ', trim (name1)
      name3 = name
      name  = name1
      iform = iform1
      call iohead (1)

      do m = 1, mmax
         call iogrid (1, m)
         call iogrid (2, m)
      end do
      close (1)

      write (*, '(2a)') 'Reading ', trim (name2)
      name = name2
      iform = iform2
      call iohead (1)

      do m = 1, mmax
         call iogrid (1, m)
         call iogrid (2, m)
      end do
      close (1)

      name = name3

      end subroutine combine

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine unitconv  ! Scale x, y, z by some factor.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: m, ni, nj, nk
      real    :: xfac

      write (*, '(a)', advance='no') 'Conversion factor: '
      read  (*, *) xfac
      call iohead (1)
      name = 'scaled.p3d'
      call iohead (2)

      do m = 1, mmax
         call iogrid (1, m)
         ni = imax(m)
         nj = jmax(m)
         nk = kmax(m)
         x(1:ni,1:nj,1:nk) = xfac*x(1:ni,1:nj,1:nk)
         y(1:ni,1:nj,1:nk) = xfac*y(1:ni,1:nj,1:nk)
         z(1:ni,1:nj,1:nk) = xfac*z(1:ni,1:nj,1:nk)
         call iogrid (2, m)
      end do

      end subroutine unitconv

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine extract  ! Extract 1 or more blocks, individually or as 1 file.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, is, j, m, mextract, n, n1, n2, n3
      character (1) :: what
      character (3) :: digit
      integer, allocatable :: mz(:)

      call iohead (1)
      allocate (mz(mmax))
      write (*, '(i4, a)') mmax, ' zones'
      write (*, '(i4, 3i5)') (m, imax(m), jmax(m), kmax(m), m = 1, mmax)
      write (*, '(/, a)', advance='no') 'Extract how many zones? '
      read  (*, *) mextract
      write (*, '(a)', advance='no') &
         'Extract as individual zones (i) or one chunk of zones (c)? '
      read  (*, *) what
      write (*, '(a, i3, a)', advance='no') &
         'Enter',  mextract, ' zones, in order: '
      read  (*, *) mz(1:mextract)

      if (what == 'c') then  ! Write a subset of blocks to one file
         if (oform == 'f') then
            name = 'chunk.f';  if (grid) name = 'chunk.p3d'
            open  (2, file=name, status='unknown')
            write (2, '(i4)') mextract
            do m = 1, mextract
               n = mz(m)
               if (grid) then
                  write (2, '(3i5)') imax(n), jmax(n), kmax(n)
               else
                  write (2, '(3i5, i4)') imax(n), jmax(n), kmax(n), qmax
               end if
            end do
         else
            name = 'chunk.fu';  if (grid) name = 'chunk.p3du'
            open  (2, file=name, status='unknown', form='unformatted')
            write (2) mextract
            if (grid) then
               write (2) &
                  (imax(mz(m)), jmax(mz(m)), kmax(mz(m)), m = 1, mextract)
            else
               write (2) &
                  (imax(mz(m)), jmax(mz(m)), kmax(mz(m)), qmax, m = 1, mextract)
            end if
         end if
         do m = 1, mz(mextract)
            call iogrid (1, m)
            do j = 1, mextract
               if (mz(j) == m) call iogrid (2, m)
            end do
         end do
      else  ! Write each block to a separate file
         is = 1
         do m = 1, mextract
            n = mz(m)
            if (n <= 9) then
               digit = char (48+n)
            else if (n <= 99) then
               n1 = n / 10
               n2 = n - 10*n1
               digit = char (48+n1) // char (48+n2)
            else  ! Assume n <= 999
               n1 =  n / 100
               n2 = (n - 100*n1) / 10
               n3 =  n - 100*n1 -  10*n2
               digit = char (48+n1) // char (48+n2) // char(48+n3)
            end if
            if (grid) then
               name = 'z' // trim (digit) // '.p3d'
            else
               name = 'z' // trim (digit) // '.f'
            end if
            if (oform == 'f') then
               open (2, file=name, status='unknown')
               write (2, '(i1)') 1
               if (grid) then
                  write (2, '(3i5)') imax(n), jmax(n), kmax(n)
               else
                  write (2, '(3i5, i4)') imax(n), jmax(n), kmax(n), qmax
               end if
            else
               name = trim (name) // 'u'
               open (2, file=name, form='unformatted', status='unknown')
               write (2) 1
               if (grid) then
                  write (2) imax(n), jmax(n), kmax(n)
               else
                  write (2) imax(n), jmax(n), kmax(n), qmax
               end if
            end if
            do i = is, n
               call iogrid (1, i)
            end do
            call iogrid (2, n)
            is = n + 1
         end do
      end if

      deallocate (mz)

      end subroutine extract

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine reorder  ! ... the blocks.
!
!     The strategy is to avoid storing the whole grid by rereading as often as
!     necessary to output the blocks in the specified order, which is read from
!     a control file containing one or more block numbers per line (and nothing
!     else).  The default control file is 'reorder.inp' although another name
!     may be used.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, ios, j, k, l, m, mm
      character (64) :: name1, name2
      integer, allocatable, dimension (:) :: imaxtemp,jmaxtemp,kmaxtemp,neworder

      name1 = name
      call iohead (1)
      close (1)
      allocate (imaxtemp(mmax), jmaxtemp(mmax), kmaxtemp(mmax), neworder(mmax))

      name2 = 'reorder.inp'
      do  ! Until control file is found
         open (2, file=name2, status='old', iostat=ios)
         if (ios == 0) exit
         write (*, '(2a)', advance='no') &
            trim (name2), ' not found; enter a control file name: '
         read  (*, *, iostat=ios) name2
         if (ios < 0) stop  ! ^D will quit cleanly
      end do

      read (2, *, iostat=ios) neworder(:)
      if (ios /= 0) then
         write (*, '(a, /, a, i5)') 'Trouble reading the new block order.', &
            'Expected number of entries:', mmax
         stop
      end if

      close (2)
      write (*, '(a, /, (15i5))') 'New order:', neworder(:)

      do m = 1, mmax
         mm = neworder(m)
         imaxtemp(m) = imax(mm)
         jmaxtemp(m) = jmax(mm)
         kmaxtemp(m) = kmax(mm)
      end do

      do m = 1, mmax
         imax(m) = imaxtemp(m)
         jmax(m) = jmaxtemp(m)
         kmax(m) = kmaxtemp(m)
      end do

      name  = 'reorder.p3d'
      call iohead (2)
      name2 = name
      name  = name1

      do m = 1, mmax
         call iohead (1)
         do mm = 1, neworder(m)  ! Read to next output block;
            if (grid) then       ! avoid the commentary of iogrid (1, mm)
               if (iform == 'f') then
                  read (1, *) &
                     (((x(i,j,k),i=1,imax(mm)),j=1,jmax(mm)),k=1,kmax(mm)), &
                     (((y(i,j,k),i=1,imax(mm)),j=1,jmax(mm)),k=1,kmax(mm)), &
                     (((z(i,j,k),i=1,imax(mm)),j=1,jmax(mm)),k=1,kmax(mm))
               else
                  read (1) &
                     (((x(i,j,k),i=1,imax(mm)),j=1,jmax(mm)),k=1,kmax(mm)), &
                     (((y(i,j,k),i=1,imax(mm)),j=1,jmax(mm)),k=1,kmax(mm)), &
                     (((z(i,j,k),i=1,imax(mm)),j=1,jmax(mm)),k=1,kmax(mm))
               end if
            else
               if (iform == 'f') then
                  read (1, *) &
                     ((((q(i,j,k,l),i=1,imax(mm)),j=1,jmax(mm)),k=1,kmax(mm)), &
                         l=1,qmax)
               else
                  read (1) &
                     ((((q(i,j,k,l),i=1,imax(mm)),j=1,jmax(mm)),k=1,kmax(mm)), &
                         l=1,qmax)
               end if
            end if
         end do
         mm = neworder(m)
         write (*, '(a, i4)') 'writing block', mm
         if (grid) then
            if (oform == 'f') then
               write (2, '(4es19.11)') &
                  (((x(i,j,k),i=1,imax(mm)),j=1,jmax(mm)),k=1,kmax(mm)), &
                  (((y(i,j,k),i=1,imax(mm)),j=1,jmax(mm)),k=1,kmax(mm)), &
                  (((z(i,j,k),i=1,imax(mm)),j=1,jmax(mm)),k=1,kmax(mm))
            else
               write (2) &
                  (((x(i,j,k),i=1,imax(mm)),j=1,jmax(mm)),k=1,kmax(mm)), &
                  (((y(i,j,k),i=1,imax(mm)),j=1,jmax(mm)),k=1,kmax(mm)), &
                  (((z(i,j,k),i=1,imax(mm)),j=1,jmax(mm)),k=1,kmax(mm))
            end if
         else
            if (oform == 'f') then
              write (2, '(4es19.11)') &
              ((((q(i,j,k,l),i=1,imax(mm)),j=1,jmax(mm)),k=1,kmax(mm)),l=1,qmax)
            else
              write (2) &
              ((((q(i,j,k,l),i=1,imax(mm)),j=1,jmax(mm)),k=1,kmax(mm)),l=1,qmax)
            end if
         end if
         rewind (1)
      end do

      name = name2

      end subroutine reorder

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine surf  ! Extract the k = 1 or k = kmax surface.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, j, kk, l, m

      write (*, '(a, /, a)', advance='no') &
         'Extract the k = 1 or k = kmax surface from all blocks.', &
         'Enter k = 1 or k = 0: '
      read (*, *) kk

      call iohead (1)

      if (kk == 0) kk = kmax(1)

      if (oform == 'f') then
         name = 'surf.p3d'
         open  (2, file=name, status='unknown')
         write (2, '(i4)') mmax
         if (grid) then
            write (2, '(3i5)') (imax(m), jmax(m), 1, m = 1, mmax)
         else
            write (2, '(3i5, i4)') (imax(m), jmax(m), 1, qmax, m = 1, mmax)
         end if
      else
         name = 'surf.p3du'
         open  (2, file=name, status='unknown', form='unformatted')
         write (2) mmax
         if (grid) then
            write (2) (imax(m), jmax(m), 1, m = 1, mmax)
         else
            write (2) (imax(m), jmax(m), 1, qmax, m = 1, mmax)
         end if
      end if

      do m = 1, mmax
         call iogrid (1, m)
         write (*, '(a, i5)') '   Writing block', m
         if (oform == 'f') then
            if (grid) then
               write (2, '(4es19.11)') &
                  ((x(i,j,kk), i = 1, imax(m)), j = 1, jmax(m)), &
                  ((y(i,j,kk), i = 1, imax(m)), j = 1, jmax(m)), &
                  ((z(i,j,kk), i = 1, imax(m)), j = 1, jmax(m))
            else
               write (2, '(4es19.11)') &
                  (((q(i,j,kk,l), i = 1, imax(m)), j = 1, jmax(m)), l = 1, qmax)
            end if
         else
            if (grid) then
               write (2) &
                  ((x(i,j,kk), i = 1, imax(m)), j = 1, jmax(m)), &
                  ((y(i,j,kk), i = 1, imax(m)), j = 1, jmax(m)), &
                  ((z(i,j,kk), i = 1, imax(m)), j = 1, jmax(m))
            else
               write (2) &
                  (((q(i,j,kk,l), i = 1, imax(m)), j = 1, jmax(m)), l = 1, qmax)
            end if
         end if
      end do

      end subroutine surf

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sequence  ! ... a grid or a function file.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, j, k, l, m, mm
      integer, allocatable, dimension (:) :: mi, mj, mk, ni, nj, nk

      call iohead (1)

      if (oform == 'f') then
         name='sequenced.p3d'
         open  (2, file=name, status='unknown')
         write (2, '(i4)') mmax
      else
         name='sequenced.p3du'
         open  (2, file=name, status='unknown', form='unformatted')
         write (2) mmax
      end if

      allocate (ni(mmax), nj(mmax), nk(mmax), &
                mi(mmax), mj(mmax), mk(mmax))

      if (mmax == 1) then
         write (*, '(a)', advance='no') 'Sequence levels (e.g., 2 2 2): '
         read  (*, *) ni, nj, nk
      else
         write (*, '(a, i4, a, /, a)') &
            'Enter sequence levels for', mmax, ' zones (ni nj nk).', &
            '0 0 0 for block 2 means sequence all blocks as for block 1.'
         do m = 1, mmax
            read (*, *) ni(m), nj(m), nk(m)
            if (m == 2 .and. ni(m) == 0) then
               do mm = 2, mmax
                  ni(mm) = ni(1)
                  nj(mm) = nj(1)
                  nk(mm) = nk(1)
                  write (*, '(i4, 3i5)') ni(1), nj(1), nk(1)
               end do
               exit
            end if
         end do
      end if

      do m = 1, mmax  ! Sequenced dimensions
         mi(m) = (imax(m)-1)/ni(m) + 1
         mj(m) = (jmax(m)-1)/nj(m) + 1
         mk(m) = (kmax(m)-1)/nk(m) + 1
      end do

      write (*, '(a, /, (3i6, a, 3i6))') &
         '    Original block         Sequenced block', &
         (imax(m), jmax(m), kmax(m), '  --> ', mi(m), mj(m), mk(m), m = 1, mmax)

      if (oform == 'f') then
         if (grid) then
            write (2, '(3i5)') (mi(m), mj(m), mk(m), m = 1, mmax)
         else
            write (2, '(3i5, i4)') (mi(m), mj(m), mk(m), qmax, m = 1, mmax)
         end if
      else
         if (grid) then
            write (2) (mi(m), mj(m), mk(m), m = 1, mmax)
         else
            write (2) (mi(m), mj(m), mk(m), qmax, m = 1, mmax)
         end if
      end if

      do m = 1, mmax
         call iogrid (1, m)
         write (*, '(a, i5)') '   Writing block', m
         if (oform == 'f') then
            if (grid) then
               write (2, '(4es19.11)') &
                  (((x(i,j,k), i=1,imax(m),ni(m)), j=1,jmax(m),nj(m)), &
                               k=1,kmax(m),nk(m)), &
                  (((y(i,j,k), i=1,imax(m),ni(m)), j=1,jmax(m),nj(m)), &
                               k=1,kmax(m),nk(m)), &
                  (((z(i,j,k), i=1,imax(m),ni(m)), j=1,jmax(m),nj(m)), &
                               k=1,kmax(m),nk(m))
            else
               write (2, '(4es19.11)') &
                  ((((q(i,j,k,l), i=1,imax(m),ni(m)), j=1,jmax(m),nj(m)), &
                                  k=1,kmax(m),nk(m)), l=1,qmax)
            end if
         else  ! Unformatted output
            if (grid) then
               write (2) &
                  (((x(i,j,k), i=1,imax(m),ni(m)), j=1,jmax(m),nj(m)), &
                               k=1,kmax(m),nk(m)), &
                  (((y(i,j,k), i=1,imax(m),ni(m)), j=1,jmax(m),nj(m)), &
                               k=1,kmax(m),nk(m)), &
                  (((z(i,j,k), i=1,imax(m),ni(m)), j=1,jmax(m),nj(m)), &
                               k=1,kmax(m),nk(m))
            else
               write (2) &
                  ((((q(i,j,k,l), i=1,imax(m),ni(m)), j=1,jmax(m),nj(m)), &
                                  k=1,kmax(m),nk(m)), l=1,qmax)
            end if
         end if
      end do

      deallocate (ni, nj, nk, mi, mj, mk)

      end subroutine sequence

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fix  ! Fix a surface file to have kmax = 1 (not imax or jmax).
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: m

      call iohead (1)

      if (oform == 'f') then
         name = 'fixed.p3d'
         open (2, file=name, status='unknown')
      else
         name = 'fixed.p3du'
         open (2, file=name, status='unknown', form='unformatted')
      end if

      do m = 1, mmax
         if (imax(m) == 1) then
            imax(m) = kmax(m)   ! What about handedness?
            kmax(m) = 1
         end if
         if (jmax(m) == 1) then
            jmax(m) = kmax(m)
            kmax(m) = 1
         end if
         if (kmax(m) /= 1) then
            write (*, '(a4, 3i5)') m, imax(m), jmax(m), kmax(m)
            write (*, '(a)') 'This is not a surface file.  Aborting ...'
            stop
         end if
      end do

      call iohead (2)

      do m = 1, mmax
         call iogrid (1, m)  ! Mismatched dimensions handily permute the indices
         call iogrid (2, m)
      end do

      end subroutine fix

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine swapijk  ! Swap any pair of indices; handedness may suffer.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, parameter :: lun = 122  ! Why not 2?

      integer :: i, j, k, l, m, swapwhat

      write (*, '(a)') &
         'Swap i, j? (1)', &
         'Swap j, k? (2)', &
         'Swap k, i? (3)'
      read (*, *) swapwhat

      call iohead (1)

      name = 'swapped.p3d'
      if (.not. grid) name = 'swapped.f'

      if (oform == 'f') then
         open (lun, file=name, status='unknown')
      else
         name = trim (name) // 'u'
         open (lun, file=name, status='unknown', form='unformatted')
      end if

      if (oform == 'f') then

         write (lun, '(i4)') mmax
         if (grid) then
            select case (swapwhat)
               case (1)
                  write (lun, '(3i5)') (jmax(m), imax(m), kmax(m), m = 1, mmax)
               case (2)
                  write (lun, '(3i5)') (imax(m), kmax(m), jmax(m), m = 1, mmax)
               case (3)
                  write (lun, '(3i5)') (kmax(m), jmax(m), imax(m), m = 1, mmax)
            end select
            do m = 1, mmax
               call iogrid (1, m)
               write (*, '(a, i5)') '   Writing block', m
               select case (swapwhat)  ! Note preservation of handedness
                  case (1)
                     write (lun, '(4es19.11)') &
                        (((x(i,j,k),j=jmax(m),1,-1),i=1,imax(m)),k=1,kmax(m)), &
                        (((y(i,j,k),j=jmax(m),1,-1),i=1,imax(m)),k=1,kmax(m)), &
                        (((z(i,j,k),j=jmax(m),1,-1),i=1,imax(m)),k=1,kmax(m))
                  case (2)
                     write (lun, '(4es19.11)') &
                        (((x(i,j,k),i=imax(m),1,-1),k=1,kmax(m)),j=1,jmax(m)), &
                        (((y(i,j,k),i=imax(m),1,-1),k=1,kmax(m)),j=1,jmax(m)), &
                        (((z(i,j,k),i=imax(m),1,-1),k=1,kmax(m)),j=1,jmax(m))
                  case (3)
                     write (lun, '(4es19.11)') &
                        (((x(i,j,k),k=kmax(m),1,-1),j=1,jmax(m)),i=1,imax(m)), &
                        (((y(i,j,k),k=kmax(m),1,-1),j=1,jmax(m)),i=1,imax(m)), &
                        (((z(i,j,k),k=kmax(m),1,-1),j=1,jmax(m)),i=1,imax(m))
               end select
            end do
         else  ! Function file
            select case (swapwhat)
               case (1)
                  write (lun, '(3i5, i4)') &
                     (jmax(m), imax(m), kmax(m), qmax, m = 1, mmax)
               case (2)
                  write (lun, '(3i5, i4)') &
                     (imax(m), kmax(m), jmax(m), qmax, m = 1, mmax)
               case (3)
                  write (lun, '(3i5, i4)') &
                     (kmax(m), jmax(m), imax(m), qmax, m = 1, mmax)
            end select
            do m = 1, mmax
               call iogrid (1, m)
               write (*, '(a, i5)') '   Writing block', m
               select case (swapwhat)   ! Note preservation of handedness
                  case (1)
                     write (lun, '(4es19.11)') &
              ((((q(i,j,k,l),j=jmax(m),1,-1),i=1,imax(m)),k=1,kmax(m)),l=1,qmax)
                  case (2)
                     write (lun, '(4es19.11)') &
              ((((q(i,j,k,l),i=imax(m),1,-1),k=1,kmax(m)),j=1,jmax(m)),l=1,qmax)
                  case (3)
                     write (lun, '(4es19.11)') &
              ((((q(i,j,k,l),k=kmax(m),1,-1),j=1,jmax(m)),i=1,imax(m)),l=1,qmax)
               end select
            end do
         end if

      else  ! Unformatted

         write (lun) mmax
         if (grid) then
            select case (swapwhat)
               case (1)
                  write (lun) (jmax(m), imax(m), kmax(m), m = 1, mmax)
               case (2)
                  write (lun) (imax(m), kmax(m), jmax(m), m = 1, mmax)
               case (3)
                  write (lun) (kmax(m), jmax(m), imax(m), m = 1, mmax)
            end select
            do m = 1, mmax
               call iogrid (1, m)
               write (*, '(a, i5)') '   Writing block', m
               select case (swapwhat)
                  case (1)
                     write (lun) &
                        (((x(i,j,k),j=jmax(m),1,-1),i=1,imax(m)),k=1,kmax(m)), &
                        (((y(i,j,k),j=jmax(m),1,-1),i=1,imax(m)),k=1,kmax(m)), &
                        (((z(i,j,k),j=jmax(m),1,-1),i=1,imax(m)),k=1,kmax(m))
                  case (2)
                     write (lun) &
                        (((x(i,j,k),i=imax(m),1,-1),k=1,kmax(m)),j=1,jmax(m)), &
                        (((y(i,j,k),i=imax(m),1,-1),k=1,kmax(m)),j=1,jmax(m)), &
                        (((z(i,j,k),i=imax(m),1,-1),k=1,kmax(m)),j=1,jmax(m))
                  case (3)
                     write (lun) &
                        (((x(i,j,k),k=kmax(m),1,-1),j=1,jmax(m)),i=1,imax(m)), &
                        (((y(i,j,k),k=kmax(m),1,-1),j=1,jmax(m)),i=1,imax(m)), &
                        (((z(i,j,k),k=kmax(m),1,-1),j=1,jmax(m)),i=1,imax(m))
               end select
            end do
         else  ! Function file
            select case (swapwhat)
               case (1)
                  write (lun) &
                     (jmax(m), imax(m), kmax(m), qmax, m = 1, mmax)
               case (2)
                  write (lun) &
                     (imax(m), kmax(m), jmax(m), qmax, m = 1, mmax)
               case (3)
                  write (lun) &
                     (kmax(m), jmax(m), imax(m), qmax, m = 1, mmax)
            end select
            do m = 1, mmax
               call iogrid (1, m)
               write (*, '(a, i5)') '   Writing block', m
               select case (swapwhat)
                  case (1)
                     write (lun) &
              ((((q(i,j,k,l),j=jmax(m),1,-1),i=1,imax(m)),k=1,kmax(m)),l=1,qmax)
                  case (2)
                     write (lun) &
              ((((q(i,j,k,l),i=imax(m),1,-1),k=1,kmax(m)),j=1,jmax(m)),l=1,qmax)
                  case (3)
                     write (lun) &
              ((((q(i,j,k,l),k=kmax(m),1,-1),j=1,jmax(m)),i=1,imax(m)),l=1,qmax)
               end select
            end do
         end if

      end if

      close (lun)

      end subroutine swapijk

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine swapxyz  ! Swap any pair of coordinates; handedness may suffer.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, j, k, m, swapwhat
      logical :: formatted

      write (*, '(a)') &
         'Switch x, y? (1)', &
         'Switch y, z? (2)', &
         'Switch z, x? (3)'
      read (*, *) swapwhat

      call iohead (1)

      name = 'switched.p3d'
      formatted = oform == 'f'

      call iohead (2)

      do m = 1, mmax
         call iogrid (1, m)
         write (*, '(a, i5)') '   Writing block', m
         select case (swapwhat)
            case (1)
               if (formatted) then
                  write (2, '(4es19.11)') &
                     (((y(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)), &
                     (((x(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)), &
                     (((z(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m))
               else
                  write (2) &
                     (((y(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)), &
                     (((x(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)), &
                     (((z(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m))
               end if
            case (2)
               if (formatted) then
                  write (2, '(4es19.11)') &
                     (((x(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)), &
                     (((z(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)), &
                     (((y(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m))
               else
                  write (2) &
                     (((x(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)), &
                     (((z(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)), &
                     (((y(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m))
               end if
            case (3)
               if (formatted) then
                  write (2, '(4es19.11)') &
                     (((z(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)), &
                     (((y(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)), &
                     (((x(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m))
               else
                  write (2) &
                     (((z(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)), &
                     (((y(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)), &
                     (((x(i,j,k),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m))
               end if
         end select
      end do

      end subroutine swapxyz

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine scale  ! Scale volume coordinates; leave surface (k = 1) alone;
                        ! option to vary the scaling with x.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, j, k, m
      real    :: delx, dely, delz, r, x1, x2, xfac, xfac1, xfac2, xp, yp, zp

      write (*, '(a, /, a)', advance='no') &
         'Move the outer boundary, leaving the k = 1 surface alone.', &
         'Constant arc length scaling [1] or scaling varying with x [2]: '
      read (*, *) m
      if (m == 1) then
         write (*, '(a)', advance='no') &
            'Scale factor (> 1 to grow, < 1 to shrink): '
         read (*, *) xfac1
         r = 0.
      else
         write (*, '(a)', advance='no') &
            'Nose scale factor (> 1 to grow, < 1 to shrink): '
         read (*, *) xfac1
         write (*, '(a)', advance='no') &
            'Base scale factor (> 1 to grow, < 1 to shrink): '
         read (*, *) xfac2
         write (*, '(a)', advance='no') 'x nose: '
         read (*, *) x1
         write (*, '(a)', advance='no') 'x base: '
         read (*, *) x2
         r = (xfac2 - xfac1)/(x2 - x1)
      end if

      call iohead (1)
      name = 'newob.p3d'
      call iohead (2)

      do m = 1, mmax
         call iogrid (1, m)
         do j = 1, jmax(m)
            do i = 1, imax(m)
               xp = x(i,j,1)
               yp = y(i,j,1)
               zp = z(i,j,1)
               xfac = xfac1 + r*(xp - x1)
               do k = 2, kmax(m)
                  delx = x(i,j,k) - xp
                  dely = y(i,j,k) - yp
                  delz = z(i,j,k) - zp
                  xp = x(i,j,k)
                  yp = y(i,j,k)
                  zp = z(i,j,k)
                  x(i,j,k) = x(i,j,k-1) + xfac*delx
                  y(i,j,k) = y(i,j,k-1) + xfac*dely
                  z(i,j,k) = z(i,j,k-1) + xfac*delz
               end do
            end do
         end do
         call iogrid (2, m)
      end do

      end subroutine scale

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine extrude  ! Turn a 2d grid into 2 planes one unit apart.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, parameter :: one  =  1.0
      integer :: i, j, m, ni, nj
      character (1) :: dir
      logical :: formatted

      write (*, '(a)', advance='no') 'Extrude in x,y,or z? '
      read  (*, *) dir

      call iohead (1)

      name = 'extruded.p3d'
      formatted = oform == 'f'

      do m = 1, mmax
         kmax(m) = 2
      end do

      call iohead (2)

      do m = 1, mmax
         ni = imax(m)
         nj = jmax(m)
         kmax(m) = 1
         call iogrid (1, m)
         write (*, '(a, i5)') '   Writing block', m
         select case (dir)
            case ('x')
               if (formatted) then
                  write (2, '(4es19.11)') &
                     ((x(i,j,1),i=1,ni),j=1,nj),     &
                     ((x(i,j,1)+one,i=1,ni),j=1,nj), &
                     ((y(i,j,1),i=1,ni),j=1,nj),     &
                     ((y(i,j,1),i=1,ni),j=1,nj),     &
                     ((z(i,j,1),i=1,ni),j=1,nj),     &
                     ((z(i,j,1),i=1,ni),j=1,nj)
               else
                  write (2) &
                     ((x(i,j,1),i=1,ni),j=1,nj),     &
                     ((x(i,j,1)+one,i=1,ni),j=1,nj), &
                     ((y(i,j,1),i=1,ni),j=1,nj),     &
                     ((y(i,j,1),i=1,ni),j=1,nj),     &
                     ((z(i,j,1),i=1,ni),j=1,nj),     &
                     ((z(i,j,1),i=1,ni),j=1,nj)
               end if
            case ('y')
               if (formatted) then
                  write (2, '(4es19.11)') &
                     ((x(i,j,1),i=1,ni),j=1,nj),     &
                     ((x(i,j,1),i=1,ni),j=1,nj),     &
                     ((y(i,j,1),i=1,ni),j=1,nj),     &
                     ((y(i,j,1)+one,i=1,ni),j=1,nj), &
                     ((z(i,j,1),i=1,ni),j=1,nj),     &
                     ((z(i,j,1),i=1,ni),j=1,nj)
               else
                  write (2) &
                     ((x(i,j,1),i=1,ni),j=1,nj),     &
                     ((x(i,j,1),i=1,ni),j=1,nj),     &
                     ((y(i,j,1),i=1,ni),j=1,nj),     &
                     ((y(i,j,1)+one,i=1,ni),j=1,nj), &
                     ((z(i,j,1),i=1,ni),j=1,nj),     &
                     ((z(i,j,1),i=1,ni),j=1,nj)
               end if
            case ('z')
               if (formatted) then
                  write (2, '(4es19.11)') &
                     ((x(i,j,1),i=1,ni),j=1,nj),     &
                     ((x(i,j,1),i=1,ni),j=1,nj),     &
                     ((y(i,j,1),i=1,ni),j=1,nj),     &
                     ((y(i,j,1),i=1,ni),j=1,nj),     &
                     ((z(i,j,1),i=1,ni),j=1,nj),     &
                     ((z(i,j,1)+one,i=1,ni),j=1,nj)
               else
                  write (2) &
                     ((x(i,j,1),i=1,ni),j=1,nj),     &
                     ((x(i,j,1),i=1,ni),j=1,nj),     &
                     ((y(i,j,1),i=1,ni),j=1,nj),     &
                     ((y(i,j,1),i=1,ni),j=1,nj),     &
                     ((z(i,j,1),i=1,ni),j=1,nj),     &
                     ((z(i,j,1)+one,i=1,ni),j=1,nj)
               end if
         end select
      end do

      end subroutine extrude

   end program gu
