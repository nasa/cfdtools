
       Program zb

! *************************************************************************
! **** This program will fix a bug in POSTFLOW that causes slight ZB   ****
! **** mismatches in output plot3d function files at symmetry planes.  ****
! **** This bug primarily affects 3D SAGe adaptions for multiple block ****
! **** grid files. Run this code, specifying the function file and the ****
! **** master interface file as inputs. The code will make sure that   ****
! **** ZB data match exactly and generate an output function file.     ****
! **** --------------------------------------------------------------- ****
! **** Michael Wright                   last modified: 02/19/03        ****
! **** David Saunders  12/20/06  Unformatted I/O version.              ****
! ****                           Note that the algorithm is imperfect. ****
! ****                           Matching of one pair of faces can     ****
! ****                           be changed by matching another pair.  ****
! ****                           A fudged interface file may work      ****
! ****                           around this weakness by treating      ****
! ****                           blocks with common edges in the       ****
! ****                           right order.                          ****
! ****   "      "       05/04/09 Changed grid/flow indexing from       ****
! ****                           q(nb,ix,jx,kx,nv) to q(ix,jx,kx,nv,nb)****
! ****                           to make the I/O more efficient, in    ****
! ****                           preparation for dealing with the      ****
! ****                           edge and corner issue.                ****
! ****                           Installed automatic detection of the  ****
! ****                           input file formatting.  The output    ****
! ****                           file format matches the input file.   ****
! ****                           Thus we need only one version now.    ****
! ****   "      "       12/16/13 The 2-D case is working now.          ****
! *************************************************************************

       Implicit none

       character, parameter :: format * 11 = 'unformatted'

       Integer i,j,k,nb,ii,jj,kk,n,ix,jx,kx,nn,m,izdum,nbl2,nitot,inz,ier
       Integer if1, leni, leno
       Integer intyp,idim,i1,i2,i3,j1,j2,j3,ns,nd,intloc(8)
       Integer, dimension(3) :: mbs,mes,mss,mbd,med,msd,md,ms
       Integer, allocatable  :: iil(:,:),nv(:),intin(:,:,:)

       Real*8 zvers
       Real*8, allocatable :: q(:,:,:,:,:)

       Logical formatted

       Character*128 inflow,outflow,inint

! **** Read input information

       write (6,111) 'Enter input PLOT3D grid/function filename:   '
       read  (5,*) inflow
       leni = len_trim (inflow)

       call determine_grid_form (inflow(1:leni), 10, formatted, ier)
       if (ier /= 0) go to 99

       if1 = 1;  if (formatted) if1 = 3

       write (6,111) 'Enter input interface filename:              '
       read  (5,*) inint

       write (6,111) 'Is this a (1) grid or (2) function file?:    '
       read  (5,*) intyp

       write (6,111) 'Enter input file dimension (2=2D; 3=3D):     '
       read  (5,*) idim

       write (6,111) 'Enter output PLOT3D grid/function filename:  '
       read  (5,*) outflow
       leno = len_trim (outflow)

       write (6,*)

       open  (10, file=inflow(1:leni), form=format(if1:11), status='old')

       if (formatted) then
          read (10, *) nb
       else
          read (10)    nb
       end if

       write (6, *) '# blocks: ', nb

       allocate (iil(3,nb), nv(nb), stat=ier)

       if (intyp == 1) then
          nv(:) = idim
          if (formatted) then
             if (idim == 2) then
                read (10, *) (iil(1,n),iil(2,n), n=1,nb)
             else
                read (10, *) (iil(1,n),iil(2,n),iil(3,n), n=1,nb)
             end if
          else
             if (idim == 2) then
                read (10)    (iil(1,n),iil(2,n), n=1,nb)
             else
                read (10)    (iil(1,n),iil(2,n),iil(3,n), n=1,nb)
             end if
          end if
       else  ! intyp == 2)
          if (formatted) then
             if (idim == 2) then
                read (10, *) (iil(1,n),iil(2,n),nv(n), n=1,nb)
             else
                read (10, *) (iil(1,n),iil(2,n),iil(3,n),nv(n), n=1,nb)
             end if
          else
             if (idim == 2) then
                read (10)    (iil(1,n),iil(2,n),nv(n), n=1,nb)
             else
                read (10)    (iil(1,n),iil(2,n),iil(3,n),nv(n), n=1,nb)
             end if
          end if
       end if

       if (idim == 2) iil(3,:) = 1

       ix = maxval (iil(1,:))
       jx = maxval (iil(2,:))
       kx = maxval (iil(3,:))
       nn = maxval (nv)
       write (6,*) 'Max. dimensions:', ix,jx,kx,nn

       allocate (q(ix,jx,kx,nn,nb), stat=ier)

       if (ier /= 0) then
          write (6, '(a)') '*** Trouble allocating q(:,:,:,:,:).  Too big?'
          write (6, *)     '    Storage needed: ', ix*jx*kx*nn*nb
          go to 99
       end if

       do n = 1,nb
          write (6, '(a, i6)') ' Processing block', n
          ii = iil(1,n)
          jj = iil(2,n)
          kk = iil(3,n)

          if (formatted) then
             read (10,*) ((((q(i,j,k,m,n), i=1,ii), j=1,jj), k=1,kk), m=1,nv(n))
          else
             read (10)   ((((q(i,j,k,m,n), i=1,ii), j=1,jj), k=1,kk), m=1,nv(n))
          end if
       enddo

       close (10)

! **** Read the face interface information:

       open (17, file=trim(inint), status='old')

       read(17,*)
       read(17,*)
       read(17,*)
       read(17,*)
       read(17,*) zvers,izdum
       read(17,*)
       read(17,*)
       read(17,*) nbl2,nitot,inz

       allocate (intin(nitot,2,8), stat=ier)

       do n = 1,nitot
          read(17,*)
          read(17,*)
          read(17,*)
          read(17,*)
          read(17,*) (intin(n,1,m), m=1,8)
          read(17,*) (intin(n,2,m), m=1,8)
       enddo

       close (17)

       Do n = 1,nitot
          intloc = intin(n,1,:)
          ns     = intloc(1)
          call gzrange (iil(:,ns),3,intloc,mbs,mes,mss)

          i1 = intloc(3)
          i2 = intloc(6)
          i3 =(intloc(2)+1)/2

          intloc = intin(n,2,:)
          nd     = intloc(1)
          call gzrange (iil(:,nd),3,intloc,mbd,med,msd)

          j1 = intloc(3)
          j2 = intloc(6)
          j3 =(intloc(2)+1)/2

          md(j1) = mbd(j1) - msd(j1)

          do i = mbs(i1),mes(i1),mss(i1)
             ms(i1) = i
             md(j1) = md(j1)  + msd(j1)
             md(j2) = mbd(j2) - msd(j2)

             do j = mbs(i2),mes(i2),mss(i2)
                ms(i2) = j
                md(j2) = md(j2)  + msd(j2)
                md(j3) = mbd(j3) - msd(j3)

                do k = mbs(i3),mes(i3),mss(i3)
                   ms(i3) = k
                   md(j3) = md(j3) + msd(j3)

                   q(md(1),md(2),md(3),:,nd) = q(ms(1),ms(2),ms(3),:,ns)
                enddo
             enddo
          enddo
       Enddo

! **** Write the "fixed" file:

       open (20, file=outflow(1:leno), form=format(if1:11), status='unknown')

       if (formatted) then
          write (20, *) nb
       else
          write (20)    nb
       end if

       if (intyp == 1) then
          if (formatted) then
             if (idim == 2) then
                write (20, '(2i7)') (iil(1,n),iil(2,n), n=1,nb)
             else
                write (20, '(3i7)') (iil(1,n),iil(2,n),iil(3,n), n=1,nb)
             end if
             do n = 1,nb
                ii = iil(1,n)
                jj = iil(2,n)
                kk = iil(3,n)

                write (20, '(6es19.11)') &
                   ((((q(i,j,k,m,n), i=1,ii), j=1,jj), k=1,kk), m=1,nv(n))
             enddo
          else
             if (idim == 2) then
                write(20) (iil(1,n),iil(2,n), n=1,nb)
             else
                write(20) (iil(1,n),iil(2,n),iil(3,n), n=1,nb)
             end if
             do n = 1,nb
                ii = iil(1,n)
                jj = iil(2,n)
                kk = iil(3,n)

                write (20) &
                   ((((q(i,j,k,m,n), i=1,ii), j=1,jj), k=1,kk), m=1,nv(n))
             enddo
          end if
       else  ! intyp == 2
          if (formatted) then
             if (idim == 2) then
                write (20, '(2i7,i4)') (iil(1,n),iil(2,n),nv(n),n=1,nb)
             else
                write (20, '(3i7,i4)') (iil(1,n),iil(2,n),iil(3,n),nv(n),n=1,nb)
             end if
             do n = 1,nb
                ii = iil(1,n)
                jj = iil(2,n)
                kk = iil(3,n)

                write (20, '(6es19.11)') &
                   ((((q(i,j,k,m,n), i=1,ii), j=1,jj), k=1,kk), m=1,nv(n))
             enddo
          else
             if (idim == 2) then
                write (20) (iil(1,n),iil(2,n),nv(n), n=1,nb)
             else
                write (20) (iil(1,n),iil(2,n),iil(3,n),nv(n), n=1,nb)
             end if
             do n = 1,nb
                ii = iil(1,n)
                jj = iil(2,n)
                kk = iil(3,n)

                write (20) &
                   ((((q(i,j,k,m,n), i=1,ii), j=1,jj), k=1,kk), m=1,nv(n))
             enddo
          end if
       end if

       close (20)

       deallocate (q,iil,nv,intin)

   99  continue

  111  format (a,$)
  112  format (a)

       end program zb

       Subroutine gzrange (iil,idim,intloc,mbq,meq,msq)

! ******************************************************************
! **** compute ranges for grid dummy cells at zone boundaries   ****
! **** written by: Michael Wright (mjwright@mail.arc.nasa.gov)  ****
! ******************************************************************

       Implicit none

       Integer, Dimension(3) :: mbq,meq,msq,iil
       Integer  idim,iface,i1,i2,i3,ide,intloc(8)

! **** Set-up and default ranges:

       ide     = intloc(1)
       iface   = intloc(2)
       i1      = intloc(3)
       i2      = intloc(6)
       i3      =(intloc(2)+1)/2

       mbq(i1) = intloc(4)
       meq(i1) = intloc(5)
       msq(i1) = sign(1, intloc(5)-intloc(4))
       mbq(i2) = intloc(7)
       meq(i2) = intloc(8)
       msq(i2) = sign(1, intloc(8)-intloc(7))

       if (mod(intloc(2),2) == 1) then
          mbq(i3) = 1
          meq(i3) = 1
          msq(i3) = 1
       else
          mbq(i3) = iil(i3)
          meq(i3) = iil(i3)
          msq(i3) = 1
       endif

! **** Correction to ranges for points rather than cells:

       if (meq(i1)-mbq(i1) /= 0) then
          if (msq(i1) == -1) mbq(i1) = mbq(i1) + 1
          if (msq(i1) ==  1) meq(i1) = meq(i1) + 1
       endif

       if (meq(i2)-mbq(i2) /= 0) then
          if (msq(i2) == -1) mbq(i2) = mbq(i2) + 1
          if (msq(i2) ==  1) meq(i2) = meq(i2) + 1
       endif

       end subroutine gzrange

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine determine_grid_form (grid_file_name, lun, formatted, ios)
!
!  Open a multiblock grid or function file (PLOT3D-type) and attempt to read the
!  number of blocks, thus determining whether it seems to be unformatted or not.
!  Close the file upon return.  Nonexistence of the named file is also trapped.
!
!  Note that use of the INQUIRE statement with the FORM or FORMATTED keyword
!  does not solve the problem because unless the file is open already, the
!  result is 'UNDEFINED' or 'UNKNOWN'.  This seems to be a Catch-22 situation.
!
!  Using this utility in a while loop checking ios allows retries in case of bad
!  file names as shown below.  Here, "reads" is part of a prompting utility
!  "reader" available from the present author, "cr" means "carriage return" or
!  "default", and entering control-D is the standard "end-of-file" or "quit"
!  indicator under Unix):
!
!     integer,   parameter :: lunin = 1, lunkbd = 5, luncrt = 6
!     character, parameter :: format * 11 = 'unformatted'
!
!     integer   :: i1, ios
!     logical   :: cr, eof, formatted
!     character :: filename * 64
!
!     :::::::
!
!     ios = 1
!     do while (ios /= 0)  ! Check existence and format
!
!        filename = 'both-halves.cvertex.gu'
!        call reads (luncrt, &
!                    'Baseline grid [<CR> = both-halves.cvertex.gu]: ',  &
!                    lunkbd, filename, cr, eof)
!        if (eof) go to 99  ! Abort, or whatever
!
!        call determine_grid_form (filename, lunin, formatted, ios)
!
!     end do
!
!     i1 = 1;  if (formatted) i1 = 3
!
!     open (lunin, file=filename, form=format(i1:11), status='OLD')
!
!     :::::::::
!
!  History:
!
!     Nov. 2005  D.A.Saunders  Original "determine_form" internal procedure as
!                              part of prepare_rapid_analysis.f90.
!     01/23/2008    "     "    Reusable "determine_grid_form" utility.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character, intent (in)  :: grid_file_name * (*)
   integer,   intent (in)  :: lun
   logical,   intent (out) :: formatted
   integer,   intent (out) :: ios      ! ios /= 0 suggests a bad file name

!  Local variables:

   integer :: l, nbf, nbu
   logical :: doubtful

   l = len_trim (grid_file_name)

!  First, try opening the file as formatted (the default):

   open (lun, file=grid_file_name(1:l), status='old', iostat=ios)

   if (ios /= 0) then
      write (*, '(2a)') ' Unable to open ', grid_file_name(1:l)
      go to 99
   end if

   read (lun, *, iostat=ios) nbf

   doubtful = ios /= 0
   if (.not. doubtful) doubtful = nbf <= 0 .or. nbf > 2000

   if (doubtful) then
!!!   write (*,*) 'Formatted read gives nb & ios = ', nbf, ios
!!!   write (*,*) 'Try again as unformatted.'
      close (lun)
      open (lun, file=grid_file_name(1:l), status='old', form='unformatted')
      formatted = .false.
      read  (lun, iostat=ios) nbu
!!!   write (*, *) 'unformatted read gives nb & ios = ', nbu, ios
      if (ios == 0) then
         write (*, *) 'The file appears to be unformatted.'
      else
         write (*, *) &
         '# blocks found from formatted and unformatted reads: ', nbf, nbu, &
         'Bad file name?'
      end if
   else
      formatted = .true.
      write (*, *) 'The file appears to be formatted.'
!!!   write (*, *) 'nb, ios:', nbf, ios
   end if

   close (lun)

99 return

   end subroutine determine_grid_form
