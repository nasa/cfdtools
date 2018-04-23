!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program revolve_grid
!
!  REVOLVE_GRID is a revision of Todd White's original "hotdog" utility for
!  turning a 2-D grid and optional flow solution into a 3-space equivalent.
!  For axisymmetric vehicles, even at angle of attack, an initial 2-space grid
!  and flow solution can provide an effective start for 3-space cases once the
!  resulting singular grid lines are removed (see RADIAL_INTERP in conjunction
!  with CAPSULE_GRID procedures).
!
!  The right (starboard) half of the body of revolution is produced with a
!  180-degree rotation.
!
!  Input  Coordinate System:    X streamwise, Y up
!  Output Coordinate System:    X streamwise, Y positive to the right, Z up
!
!  The input dataset[s] are rotated about OX by the indicated angle (normally
!  180 or perhaps 360 degrees) to produce the specified number of grid points
!  uniformly spaced in the azimuthal direction.
!
!  The input 2-space grid may or may not include Z = 0.  Reading is attempted
!  first as multiblock 3-D (X/Y/Z), then retried as 2-D if necessary.
!  All file names and other control inputs are prompted for, along with the
!  number of species densities, ns, from which any velocity component indices
!  are deduced.
!
!  Any input flow field file should normally contain state variables in the
!  following order, with zero or more extra variables appended:
!
!  species densities (ns of them), vx & vy (no vz), temperature[s], any extras
!
!  However, this does not suit flow inputs to a radiation solver such as NEQAIR,
!  which expects four temperatures followed by species number densities.  There-
!  fore, provision is made for NO velocity components, meaning all of the flow
!  variables are simply replicated at each azimumthal station.
!
!  The resulting output file MAY have one extra flow variable (3rd v component)
!  with the y and z components adjusted appropriately at each angular station
!  for axial symmetry.
!
!  11/09/2011  David Saunders  Initial adaptation of Todd White's hotdog tool.
!  05/10/2012    "      "      Testing of the function file option had been
!                              overlooked, and there were indexing errors.
!                              Apologies to Jay Hyatt.
!  04/08/2014    "      "      Todd found that flow datasets for radiation
!                              calculations don't contain velocity components,
!                              so there's now an option to handle that.
!  04/09/204     "      "      Allow for Z = 0. or not (3-D or 2-D) in the
!                              input grid (automatically detected).
!
!  Authors:  Todd White/David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure   ! Derived data type for one grid block
   use xyzq_io_module         ! I/O package for PLOT3D files ...
   use xyq_io_module          ! ... and its 2-space analogue

   implicit none

!  Local constants:

   integer, parameter :: &
      luncrt  = 6,       &    ! Screen
      lunkbd  = 5,       &    ! Keyboard
      luning  = 1,       &    ! Input 2-D grid
      luninf  = 2,       &    ! Input 2-D flow field
      lunoutg = 3,       &    ! Output 3-D grid
      lunoutf = 4             ! Output 3-D flow field

   real, parameter :: &
      zero = 0.

!  Local variables:

   integer :: &
      i, ib, ios, j, k, nblocks, ndim, nfin, nfout, ni, nj, nk, npts, ns, numj

   real :: &
      angle, delta_angle, total_angle

   logical :: &
      formatted_in, formatted_out, qfile, velocities

   character :: &
      answer*1

   real, allocatable :: &
      vy(:,:), vz(:,:)

   type (grid_type), pointer, dimension (:) :: &
      xyzq_in, xyzq_out

!  Execution:

   call file_prompt (luning, -luncrt, 'input grid', 'old', .true., &
                     formatted_in, ios)
   if (ios /= 0) go to 99

   call determine_grid_dim ('none', luning, formatted_in, ndim, ios)
   if (ios /= 0) go to 99

   write (luncrt, '(a)', advance='no') ' Is there a function file too? [y|n]: '
   read  (lunkbd, *) answer
   qfile = answer == 'y' .or. answer == 'Y'

   if (qfile) then
      call file_prompt (luninf, -luncrt, 'input flow', 'old', .false., &
                        formatted_in, ios)
      if (ios /= 0) go to 99
   end if

!  Read the input file header(s).

   if (ndim == 3) then
      call xyz_header_io (1, luning, formatted_in, nblocks, xyzq_in, ios)
      if (ios /= 0) go to 99
   else
      call xy_header_io (1, luning, formatted_in, nblocks, xyzq_in, ios)
      if (ios /= 0) go to 99

      xyzq_in(:)%nk = 1
   end if

   if (qfile) then
      if (ndim == 3) then
         call q_header_io (1, luninf, formatted_in, nblocks, nfin, xyzq_in, ios)
      else
         call q_header_io_2d (1, luninf, formatted_in, nblocks, nfin, xyzq_in, &
                              ios)
         xyzq_in(:)%mk = 1
      end if
      if (ios /= 0) go to 99
   end if

   write (luncrt, '(a)', advance='no') ' Total angle of revolution: '
   read  (lunkbd, *) total_angle
   write (luncrt, '(a)', advance='no') ' Number of circumferential grid pts.: '
   read  (lunkbd, *) numj
   delta_angle = total_angle / real (numj - 1)

   if (qfile) then
      write (luncrt, '(a)', advance='no') &
     ' # species preceding velocities [ns > 0, or 0 if (u,v) are not present]: '
      read  (lunkbd, *) ns
      velocities = ns > 0
      if (velocities) then
         nfout = nfin + 1
      else
         nfout = nfin
      end if
   else
      nfin = 0;  nfout = 0
   end if

!  Set up the output file(s):

   call file_prompt (lunoutg, -luncrt, 'output grid', 'unknown', .true., &
                     formatted_out, ios)
   if (ios /= 0) go to 99

   if (qfile) then
      call file_prompt (lunoutf, -luncrt, 'output flow', 'unknown', .false., &
                        formatted_out, ios)
      if (ios /= 0) go to 99
   end if

   allocate (xyzq_out(nblocks))

!  Set up the output block dimensions:

   do ib = 1, nblocks
      xyzq_out(ib)%ni = xyzq_in(ib)%ni
      xyzq_out(ib)%nj = numj
      xyzq_out(ib)%nk = xyzq_in(ib)%nj
      if (qfile) then
         xyzq_out(ib)%mi = xyzq_in(ib)%mi
         xyzq_out(ib)%mj = numj
         xyzq_out(ib)%mk = xyzq_in(ib)%mj
      end if
   end do

   call xyz_header_io (2, lunoutg, formatted_out, nblocks, xyzq_out, ios)
   if (ios /= 0) go to 99

   if (qfile) then
      call q_header_io (2, lunoutf, formatted_out, nblocks, nfout, xyzq_out, &
                        ios)
      if (ios /= 0) go to 99
   end if

!  Process the blocks one at a time:

   do ib = 1, nblocks

      write (luncrt, '(a, i5)') ' Grid block #:', ib

      call xyz_allocate (xyzq_in(ib), ios)
      if (ios /= 0) go to 99

      ni = xyzq_in(ib)%ni;  nj = xyzq_in(ib)%nj
      npts = ni * nj

      if (ndim == 3) then
         call xyz_block_io (1, luning, formatted_in, npts, &
                            xyzq_in(ib)%x, xyzq_in(ib)%y, xyzq_in(ib)%z, ios)
         if (ios /= 0) go to 99
      else
         call xy_block_io (1, luning, formatted_in, npts, &
                           xyzq_in(ib)%x, xyzq_in(ib)%y, ios)
         if (ios /= 0) go to 99
         xyzq_in(ib)%z(:,:,:) = zero
      end if

      call xyz_allocate (xyzq_out(ib), ios)
      if (ios /= 0) go to 99

      allocate (vy(ni,nj), vz(ni,nj))

      do j = 1, numj
         angle = -delta_angle * real (j - 1)  ! Negative for right-handedness
         if (j == numj) angle = -total_angle  ! Exactly

         xyzq_out(ib)%x(:,j,:) = xyzq_in(ib)%x(:,:,1)
         vy(:,:)               = zero
         vz(:,:)               = xyzq_in(ib)%y(:,:,1)

         if (j > 1) call rotate2d (npts, vy, vz, angle, zero, zero)

         xyzq_out(ib)%y(:,j,:) = vy(:,:)
         xyzq_out(ib)%z(:,j,:) = vz(:,:)
      end do

      deallocate (xyzq_in(ib)%x, xyzq_in(ib)%y, xyzq_in(ib)%z, vy, vz)

      npts = npts * numj

      call xyz_block_io (2, lunoutg, formatted_out, npts, &
                         xyzq_out(ib)%x, xyzq_out(ib)%y, xyzq_out(ib)%z, ios)
      if (ios /= 0) go to 99

      deallocate (xyzq_out(ib)%x, xyzq_out(ib)%y, xyzq_out(ib)%z)

   end do  ! Next grid block

   close (lunoutg)

   if (qfile) then

      do ib = 1, nblocks

         write (luncrt, '(a, i5)') ' Flow block #:', ib

         call q_allocate (xyzq_in(ib), nfin, ios)
         if (ios /= 0) go to 99

         ni = xyzq_in(ib)%mi
         nj = xyzq_in(ib)%mj;  npts = ni * nj
         nk = nj

         if (ndim == 2) then
            call q_block_io_2d (1, luninf, formatted_in, nfin, ni, nj, &
                                xyzq_in(ib)%q, ios)
         else
            call q_block_io (1, luninf, formatted_in, nfin, ni, nj, 1, &
                             xyzq_in(ib)%q, ios)
         end if
         if (ios /= 0) go to 99

         call q_allocate (xyzq_out(ib), nfout, ios)

         if (velocities) then

            allocate (vy(ni,nj), vz(ni,nj))

            do j = 1, numj
               do k = 1, nk
                  xyzq_out(ib)%q(   1:ns+1, :,j,k) = &
                     xyzq_in(ib)%q(   1:ns+1,:,k,1)
                  xyzq_out(ib)%q(ns+4:nfout,:,j,k) = &
                     xyzq_in(ib)%q(ns+3:nfin,:,k,1)
               end do

               angle = -delta_angle * real (j - 1)  
               if (j == numj) angle = -total_angle  ! Exactly

               vy(:,:) = zero
               vz(:,:) = xyzq_in(ib)%q(ns+2,:,:,1)

               if (j > 1) call rotate2d (npts, vy, vz, angle, zero, zero)

               xyzq_out(ib)%q(ns+2,:,j,:) = vy(:,:)
               xyzq_out(ib)%q(ns+3,:,j,:) = vz(:,:)
            end do

            deallocate (vy, vz)

         else  ! No velocity components

            do k = 1, nk
               do j = 1, numj
                  xyzq_out(ib)%q(:,:,j,k) = xyzq_in(ib)%q(:,:,k,1)
               end do
            end do

         end if

         deallocate (xyzq_in(ib)%q)

         call q_block_io (2, lunoutf, formatted_out, nfout, ni, numj, nk, &
                          xyzq_out(ib)%q, ios)
         if (ios /= 0) go to 99

         deallocate (xyzq_out(ib)%q)

      end do  ! Next flow block

      close (lunoutf)

   end if

99 continue

   end program revolve_grid
