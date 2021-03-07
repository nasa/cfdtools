!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program merge_files

!  Merge the contents of two Tecplot files which apply to the same structured
!  grid (2-D or 3-D, same number of zones and same zone sizes).  The coordinates
!  are assumed to match, so only those from the first file are carried over.
!  The inputs must be ASCII files; the merged output may be ASCII or binary.
!  File types are deduced from their names (*.dat or *.plt), which are prompted
!  for.  Zones are processed one at a time (no need to store more than one).
!
!  10/31/06  D.A.Saunders  Initial implementation.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_header_structure  ! See Tecplot_io.f90 for these modules
   use grid_block_structure
   use tecplot_io_module

   implicit none

!  Constants:

   integer, parameter :: &
      lunin1 = 1,        &    ! First input file  (Tecplot ASCII)
      lunin2 = 2,        &    ! Second input file (  "      "   )
      lunout = 3,        &    ! Output merged file (Tecplot ASCII or binary)
      lunkbd = 5,        &    ! Keyboard entries
      luncrt = 6              ! For prompts and diagnostics

   logical, parameter :: &
      false  = .false.,  &
      true   = .true.

!  Variables:

   integer :: &
      i, ib, ios, ipack, j, k, len, n1, n2, naux, nblocks, ndim, nf1, nf2,     &
      nfout, num

!  Composite data types:

   type (grid_header) :: &
      header_in1, header_in2, header_out

   type (grid_type), pointer, dimension (:) :: &
      zones_in1, zones_in2, zones_out

!  Execution:
!  !!!!!!!!!!

   write (luncrt, '(/, a)', advance='no') ' Input file # 1 (Tecplot ASCII):   '
   read  (lunkbd, *) header_in1%filename
   header_in1%formatted = true
   header_in1%ndim      = -1   ! Option to determine it from presence of 'z'

   write (luncrt, '(   a)', advance='no') ' Input file # 2 (Tecplot ASCII):   '
   read  (lunkbd, *) header_in2%filename
   header_in2%formatted = true
   header_in2%ndim      = -1

   write (luncrt, '(   a)', advance='no') ' Output merged file [*.dat|*.plt]: '
   read  (lunkbd, *) header_out%filename
   write (luncrt, '(a)')
   len  =  len_trim (header_out%filename)
   header_out%formatted = header_out%filename(len-2:len) == 'dat'

!  Open the input files and read their headers:

   ios = 1  ! Verbose mode
   call Tec_header_read (lunin1, header_in1, zones_in1, ios)
   if (ios /= 0) go to 99

   ios = 1
   call Tec_header_read (lunin2, header_in2, zones_in2, ios)
   if (ios /= 0) go to 99

!  Check their compatibility:

   call check_compatibility ()  ! Internal procedure sets ndim, nf1, nf2, nfout

   if (ios /= 0) go to 99

   header_out%ndim    = ndim
   header_out%numq    = nfout
   header_out%nblocks = nblocks

   write (luncrt, '(a)', advance='no') &
      ' Output datapacking [0 => POINT | 1 => BLOCK]: '
   read  (lunkbd, *) header_out%datapacking
   header_out%title = header_in1%title

   num = ndim + nf1
   allocate (header_out%varname(num + nf2))

   header_out%varname(1:num) = header_in1%varname(:)
   header_out%varname(num+1:num+nf2) = header_in2%varname(ndim+1:ndim+nf2)

   n1   = header_in1%ndatasetaux
   n2   = header_in2%ndatasetaux
   naux = n1 + n2

   if (naux > 0) then
      allocate (header_out%datasetauxname(naux), &
                header_out%datasetauxvalue(naux))
      header_out%datasetauxname(1:n1)      = header_in1%datasetauxname(:)
      header_out%datasetauxname(n1+1:naux) = header_in2%datasetauxname(:)
   end if

!  Write the output file header:

   ios = 1  ! Verbose mode

   call Tec_header_write (lunout, header_out, ios)
   if (ios /= 0) go to 99

   allocate (zones_out(nblocks))

!  Process one zone at a time:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!

   do ib = 1, nblocks

      call Tec_block_allocate (zones_in1(ib), ndim, nf1, ios)
      if (ios /= 0) go to 99

      call Tec_block_read (lunin1, header_in1, zones_in1(ib), ios)
      if (ios /= 0) then
         write (luncrt, '(/, a, i5)') ' Trouble reading file 1 block #', ib
         go to 99
      end if

      call Tec_block_allocate (zones_in2(ib), ndim, nf2, ios)
      if (ios /= 0) go to 99

      call Tec_block_read (lunin2, header_in2, zones_in2(ib), ios)
      if (ios /= 0) then
         write (luncrt, '(/, a, i5)') ' Trouble reading file 2 block #', ib
         go to 99
      end if

      zones_out(ib)%zone_title   = zones_in1(ib)%zone_title
      zones_out(ib)%nzoneaux     = 0   ! For now
      zones_out(ib)%solutiontime = zones_in1(ib)%solutiontime
      zones_out(ib)%ni           = zones_in1(ib)%ni
      zones_out(ib)%nj           = zones_in1(ib)%nj
      zones_out(ib)%nk           = zones_in1(ib)%nk
      zones_out(ib)%mi           = zones_in1(ib)%mi
      zones_out(ib)%mj           = zones_in1(ib)%mj
      zones_out(ib)%mk           = zones_in1(ib)%mk

      call Tec_block_allocate (zones_out(ib), ndim, nfout, ios)
      if (ios /= 0) go to 99

      zones_out(ib)%x = zones_in1(ib)%x
      zones_out(ib)%y = zones_in1(ib)%y
      zones_out(ib)%z = zones_in1(ib)%z

      do k = 1, zones_out(ib)%mk
         do j = 1, zones_out(ib)%mj
            do i = 1, zones_out(ib)%mi
               zones_out(ib)%q(1:nf1,i,j,k)       = zones_in1(ib)%q(1:nf1,i,j,k)
               zones_out(ib)%q(nf1+1:nfout,i,j,k) = zones_in2(ib)%q(1:nf2,i,j,k)
            end do
         end do
      end do

      call Tec_block_write (lunout, header_out, zones_out(ib), ios)
      if (ios /= 0) then
         write (luncrt, '(/, a, i5)') ' Trouble writing merged block #', ib
         go to 99
      end if

      call deallocate_blocks (ib, ib, ndim, nf1, zones_in1, ios)
      if (ios /= 0) then
         write (luncrt, '(/, a, i5)') ' Trouble deallocating file 1 block #', ib
         go to 99
      end if

      call deallocate_blocks (ib, ib, ndim, nf2, zones_in2, ios)
      if (ios /= 0) then
         write (luncrt, '(/, a, i5)') ' Trouble deallocating file 2 block #', ib
         go to 99
      end if

      call deallocate_blocks (ib, ib, ndim, nfout, zones_out, ios)
      if (ios /= 0) then
         write (luncrt, '(/, a, i5)') ' Trouble deallocating merged block #', ib
         go to 99
      end if

   end do ! Next zone

   call deallocate_header (header_in1);  close (lunin1)
   call deallocate_header (header_in2);  close (lunin2)
   call deallocate_header (header_out);  close (lunout)

99 continue  ! Avoid system-dependent STOP behavior

!  Internal procedures for program merge_files:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine check_compatibility ()

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: l1, l2, m, n

      ndim = header_in1%ndim

      if (ndim /= header_in2%ndim) then
         write (luncrt, '(/, a, 2i5)') &
            ' Different file dimensions:', ndim, header_in2%ndim
         ios = 1
         go to 99
      end if

      nblocks = header_in1%nblocks

      if (nblocks /= header_in2%nblocks) then
         write (luncrt, '(/, a, 2i5)') &
            ' Different zone counts:', nblocks, header_in2%nblocks
         ios = 1
         go to 99
      end if

      do ib = 1, nblocks

         if (zones_in1(ib)%ni /= zones_in2(ib)%ni) ios = 1
         if (zones_in1(ib)%nj /= zones_in2(ib)%nj) ios = 1
         if (header_in1%ndim  == 3) then
            if (zones_in1(ib)%nk /= zones_in2(ib)%nk) ios = 1
         else
             zones_in1(ib)%nk = 1;  zones_in2(ib)%nk = 1
         end if
 
         if (ios /= 0) then
            write (luncrt, '(/, a, i4, a, /, (a, 3i6))') &
               ' Different sizes for zone # ', ib, ':', &
               '    1:', zones_in1(ib)%ni, zones_in1(ib)%nj, zones_in1(ib)%nk, &
               '    2:', zones_in2(ib)%ni, zones_in2(ib)%nj, zones_in2(ib)%nk
            go to 99
         end if

      end do

      nf1   = header_in1%numq
      nf2   = header_in2%numq
      nfout = nf1 + nf2

      do m = 1, nf1
         l1 = len_trim (header_in1%varname(ndim + m))
         do n = 1, nf2
            l2 = len_trim (header_in2%varname(ndim + n))
            if (l1 == l2) then
               if (header_in1%varname(ndim+m)(1:l1) ==    &
                   header_in2%varname(ndim+n)(1:l1)) then
                  write (luncrt, '(/, a, /, (a, i4, 2x, a))') &
                     ' Non-unique function names:', &
                     ' File 1, fn. #', m, header_in1%varname(ndim+m)(1:l1),    &
                     ' File 2, fn. #', n, header_in2%varname(ndim+n)(1:l1)
                  ios = 1
                  go to 99
               end if
            end if
         end do
      end do

99    return

      end subroutine check_compatibility

   end program merge_files
