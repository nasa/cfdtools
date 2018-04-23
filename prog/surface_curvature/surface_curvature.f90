!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program surface_curvature
!
!  Description:
!
!     For a structured multiblock surface or volume grid, calculate four forms
!  of curvature and save the surface results in Tecplotable form.  If the grid
!  is a volume grid, determine the face of each block most likely to be a solid
!  wall (via the smallest average off-wall increment) and output results for
!  that face of the block.
!
!     This version has the option to add unit surface normals to the output.
!
!     The two main forms of curvature are Gaussian and mean curvature, from
!  which two principal curvatures can be derived.  Most of the calculations are
!  performed in subroutine "gaussian_curvature", q.v. for further details.
!
!     Note that the partial derivatives of x, y, z with respect to arc length
!  along the grid lines are not true partial derivatives if the surface grid
!  lines are not orthogonal.  Thus, in general, these curvature results are
!  approximate only, less so the more orthogonal the surface grid lines.  A
!  hemispherical test case gives very good results on the three main patches
!  but the end cap patches that avoid singular points have almost 180-degree
!  interior angles.  The method breaks down at those corners (though it is
!  safeguarded), and is noisy near them.
!
!     Another weakness is that there is no guarantee of continuity across
!  surface patches at common boundaries, since only one patch is treated at a
!  time, with one-sided differencing at grid line end points.
!
!  History:
!
!     09/07/05  D. A. Saunders  Initial implementation, for CEV capsules.
!     09/08/05     "     "      Added output of the principal curvatures.
!     09/09/05     "     "      Handled volume grids as well as surface grids.
!     09/12/05     "     "      More careful solution of the quadratic for the
!                               principal curvatures.
!     10/03/06     "     "      Installed Tecplot 360 version of the I/O pkg.
!     11/27/07     "     "      A question from Todd White about surface normals
!                               led to adding such an option here.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_header_structure  ! See Tecplot_io.f90 for this module
   use grid_block_structure   !  "     "     "
   use xyzq_io_module
   use tecplot_io_module      !  "     "     "

   implicit none

!  Constants:

   integer, parameter :: &
      lunkbd = 5,        &
      luncrt = 6,        &
      lunxyz = 2,        &
      lunplt = 3,        &
      ndim   = 3

   real, parameter :: &
      one = 1., zero = 0.

   character, parameter :: &
      blank * 1   = ' ',   &
      format * 11 = 'unformatted'

!  Variables:

   integer :: &
      i, i1, i2, ib, ic, id, ios, iwall, i_array(1), j, jc, len, nblocks, nf,  &
      ni, nj, nk, npts

   real :: &
      average_ds(6), H, K, k1, k2, p, q, term, unit_normal(3)

   real, allocatable, dimension (:,:) :: &  ! The grid_type (n,i,j,k) order is
      curv_gauss, curv_mean                 ! incompatible with the crv. utility

   logical :: &
      formatted_in, normals, surface

   character :: &
      answer * 1, filename * 80

!  Composite data types:

   type (grid_header) :: &
      header_surf

   type (grid_type), pointer, dimension (:) :: &
      grid, surf

!  Tecplot function:

   integer :: TecEnd110

!  Execution:

   write (luncrt, '(/, a)', advance='no') ' Input grid name: '
   read  (lunkbd, *) filename
   write (luncrt, '(a)', advance='no') ' Formatted?  [y|n]: '
   read  (lunkbd, *) answer
   formatted_in = answer == 'y' .or. answer == 'Y'
   i1 = 1;  if (formatted_in) i1 = 3

   open (lunxyz, file=filename, form=format(i1:11), status='old', iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') &
         ' Unable to open input grid: ', filename(1:len_trim(filename))
      go to 99
   end if

   nf = 4
   write (luncrt, '(a)', advance='no') &
      ' Unit normals as well as curvature?  [y|n]: '
   read  (lunkbd, *) answer
   normals = answer == 'y' .or. answer == 'Y'
   if (normals) nf = 7

!  Allocate grid coordinate arrays and read the input file header records:

   call xyz_header_io (1, lunxyz, formatted_in, nblocks, grid, ios)
   if (ios /= 0) go to 99

   surface = grid(1)%nk == 1

!  Set up the output surface curvature file header, although the patch
!  dimensions may change for some input volume grids:

   write (luncrt, '(a)', advance='no') &
      ' Output curvature file name (*.dat | *.plt => Tec ASCII | Tec binary): '

   read  (lunkbd, *) header_surf%filename
   len  =  len_trim (header_surf%filename)
   header_surf%formatted   = header_surf%filename(len-2:len) /= 'plt'
   header_surf%ndim        = ndim
   header_surf%numq        = nf
   header_surf%nblocks     = nblocks
   header_surf%datapacking = 0  ! 0 = POINT; 1 = BLOCK
   header_surf%title       = filename
   header_surf%ndatasetaux = 0

   allocate (header_surf%varname(3 + nf))

   header_surf%varname(1)  = 'x'
   header_surf%varname(2)  = 'y'
   header_surf%varname(3)  = 'z'
   header_surf%varname(4)  = 'Gaussian curvature'
   header_surf%varname(5)  = 'Mean curvature'
   header_surf%varname(6)  = 'Principal k1'
   header_surf%varname(7)  = 'Principal k2'

   if (normals) then
      header_surf%varname(8)  = 'ux'
      header_surf%varname(9)  = 'uy'
      header_surf%varname(10) = 'uz'
   end if

   allocate (surf(nblocks))

   do ib = 1, nblocks
      surf(ib)%zone_title   = blank
      write (surf(ib)%zone_title(1:8), '(a, i4)') 'Zone', ib
      surf(ib)%nzoneaux     = 0
      surf(ib)%solutiontime = -999.  ! Undefined
      surf(ib)%ni = grid(ib)%ni;  surf(ib)%mi = grid(ib)%ni
      surf(ib)%nj = grid(ib)%nj;  surf(ib)%mj = grid(ib)%nj
      surf(ib)%nk = 1;            surf(ib)%mk = 1
   end do

   call Tec_header_write (lunplt, header_surf, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble writing the output file header.'
      go to 99
   end if

!  Process one block at a time:

   do ib = 1, nblocks

      call xyz_allocate (grid(ib), ios)
      if (ios /= 0) go to 99

      ni = grid(ib)%ni;  nj = grid(ib)%nj;  nk = grid(ib)%nk
      npts = ni * nj * nk

      call xyz_block_io (1, lunxyz, formatted_in, npts, grid(ib)%x,            &
                         grid(ib)%y, grid(ib)%z, ios)
      if (ios /= 0) go to 99

      if (surface) then ! k = 1 = face 5 is the wall

         iwall = 5

      else ! Determine the most likely solid wall face

         call average_increments (ni, nj, nk, grid(ib)%x, grid(ib)%y,          &
                                  grid(ib)%z, average_ds)

         i_array = minloc (average_ds)
         iwall   = i_array(1)

         id = (iwall + 1) / 2 ! Convert face # to dimension 1, 2, or 3
         i1 = mod (id, 3) + 1 ! 1st dimension of face iwall
         i2 = mod (i1, 3) + 1 ! 2nd ....................

!        Variables i1 and i2 come out cyclic for id = 1, 2, 3.
!        I.e., they are 2,3 and 3,1 and 1,2 respectively, meaning
!        nj,nk and nk,ni and ni,nj respectively.

         select case (iwall)

            case (1, 2) ! i = 1 or i = ni face

               surf(ib)%ni = nj;  surf(ib)%mi = nj
               surf(ib)%nj = nk;  surf(ib)%mj = nk

            case (3, 4) ! j = 1 or j = nj face

               surf(ib)%ni = nk;  surf(ib)%mi = nk
               surf(ib)%nj = ni;  surf(ib)%mj = ni

            case (5, 6) ! k = 1 or k = nk

               ! Already set above for the input surface grid case

         end select

      end if

      call Tec_block_allocate (surf(ib), ndim, nf, ios)
      if (ios /= 0) go to 99

!     Transcribe the surface or volume face to the output surface patch:

      call copy_face (ni, nj, nk, grid(ib)%x, grid(ib)%y, grid(ib)%z, iwall,   &
                      surf(ib)%ni, surf(ib)%nj,                                &
                      surf(ib)%x,  surf(ib)%y,  surf(ib)%z)

      deallocate (grid(ib)%x, grid(ib)%y, grid(ib)%z)

      ni = surf(ib)%ni;  nj = surf(ib)%nj  ! From now on

!     We can't use surf(ib)%q because of its (n,i,j,k) indexing.

      allocate (curv_gauss(ni,nj), curv_mean(ni,nj))

!     Gaussian and mean curvature share the bulk of their computations:

      call gaussian_curvature (ni, nj, surf(ib)%x, surf(ib)%y, surf(ib)%z,     &
                               curv_gauss, curv_mean)

!     Transcribe to the Tecplot I/O package's data structure, deriving the
!     two principal curvatures as roots of a quadratic in the process:

      do j = 1, nj
         do i = 1, ni
            K = curv_gauss(i,j)
            H = curv_mean (i,j)
            term = sqrt (max (zero, H**2 - K))
            if (H > zero) then ! Avoid possible catastrophic cancellation
               k1 = H + term
               k2 = K / k1
            else if (H < zero) then
               k2 = H - term
               k1 = K / k2
            else
               k1 =  term
               k2 = -term
            end if
            surf(ib)%q(1,i,j,1) = K
            surf(ib)%q(2,i,j,1) = H
            surf(ib)%q(3,i,j,1) = k1
            surf(ib)%q(4,i,j,1) = k2
         end do
      end do

      deallocate (curv_gauss, curv_mean)

      if (normals) then

         q = zero

         do j = 1, nj

            if (j < nj) then ! Must point to lower left cell corner
               jc = j
            else
               jc = nj - 1;  q = one
            end if

            p = zero

            do i = 1, ni

               if (i < ni) then
                  ic = i
               else
                  ic = ni - 1;  p = one
               end if

               call surface_normal (ni, nj, surf(ib)%x, surf(ib)%y, surf(ib)%z,&
                                    ic, jc, p, q, unit_normal)

               surf(ib)%q(5,i,j,1) = unit_normal(1)
               surf(ib)%q(6,i,j,1) = unit_normal(2)
               surf(ib)%q(7,i,j,1) = unit_normal(3)

            end do

         end do

      end if

      call Tec_block_write (lunplt, header_surf, surf(ib), ios)
      if (ios /= 0) go to 99

      deallocate (surf(ib)%x, surf(ib)%y, surf(ib)%z, surf(ib)%q)

   end do ! Next block

   if (header_surf%formatted) then
      close (lunplt)
   else
      ios = TecEnd110 ()
   end if

99 continue

!! stop ! Avoid system dependencies

   end program surface_curvature
