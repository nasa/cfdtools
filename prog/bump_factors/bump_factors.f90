!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program bump_factors
!
!  Description:
!
!  Map each point of a damaged surface flow solution to the smooth surface, and
!  convert the damaged flow quantities to ratios of damaged to undamaged values
!  by locating the smooth interpolated point nearest to each damaged point.
!
!  Sample control file (standard input):
!
!     BUMP_FACTORS control file
!     3      ! Function in damage data to be converted to bump factor
!     3      ! Function in smooth data to divide by
!     0      ! 1 = flip y of damage data
!     1.0000 ! Scale applied to damage function
!            ! Damage file first, smooth file 2nd, output 3rd
!     140-02_CASE6_surface-008.dat
!     case6-rtf_221_surface.dat
!     140-02_CASE6_surface-008_BF.dat  ! Use *.plt to get binary output
!
!  07/31/05  David Saunders  Initial implementation for Shuttle mission STS-114.
!  08/10/06     "     "      Extended as a test vehicle for Tecplot 360 I/O.
!  08/12/13     "     "      All ADT variants are now in one module with generic
!            Now ERC, Inc.   build_adt and search_adt interfaces.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_header_structure  ! The first 3 of these are part of Tecplot_io.f90
   use grid_block_structure
   use tecplot_io_module
   use adt_utilities          ! All variants of the ADT search package

   implicit none

!  Constants:

   integer, parameter :: &
      lun_damage =  1,   &  ! Input damaged surface CFD soln. (Tecplot ASCII)
      lun_smooth =  2,   &  ! Input smooth OML CFD solution (Tecplot ASCII)
      lun_out    =  3,   &  ! Output solution of bump factors
      lunctl     =  5,   &  ! Control file (standard input)
      luncrt     =  6       ! For diagnostics

   real, parameter ::    &
      one  = 1.0,        &
      zero = 0.0

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

!  Variables:

   integer :: &
      ib, iflip, ios, iq_damage, iq_smooth, len, naux, &
      nb_damage, nb_smooth, nquad, numq

   real :: &
      damage_scale       ! To handle W/m^2 vs. W/cm^2

   integer, allocatable, dimension (:,:) :: &
      conn               ! For patch & (i,j) of surface quads. being searched

   logical :: &
      formatted

!  Composite data types:

   type (grid_header) :: &
      bfactors_header, damage_header, smooth_header

   type (grid_type), pointer, dimension (:) :: &
      bfactors, damage, smooth  ! (x,y,z,f) for the two surface solutions

!  Execution:

   read (lunctl, *)
   read (lunctl, *) iq_damage
   read (lunctl, *) iq_smooth
   read (lunctl, *) iflip
   read (lunctl, *) damage_scale
   read (lunctl, *)
   read (lunctl, *) damage_header%filename

   damage_header%formatted = true
   damage_header%ndim      = 3

   ios = 1  ! Turn on verbose mode

   call Tecplot_read (lun_damage, damage_header, damage, ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Trouble reading damage Tecplot file ', &
         trim (damage_header%filename)
      go to 99
   end if

   nb_damage = damage_header%nblocks
   numq      = damage_header%numq

   if (iflip /= 0) then
      do ib = 1, nb_damage
         damage(ib)%y = -damage(ib)%y
      end do
   end if

   if (damage_scale /= one) then
      do ib = 1, nb_damage
         damage(ib)%q(iq_damage,:,:,:) = damage_scale * &
         damage(ib)%q(iq_damage,:,:,:)
      end do
   end if

   read (lunctl, *) smooth_header%filename

   smooth_header%formatted = true
   smooth_header%ndim      = 3

   call Tecplot_read (lun_smooth, smooth_header, smooth, ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Trouble reading smooth Tecplot file ', &
         trim (smooth_header%filename)
      go to 99
   end if

   nb_smooth = smooth_header%nblocks
   numq      = smooth_header%numq

   read  (lunctl, *) bfactors_header%filename
   close (lunctl)

!  Set up the output dataset:

   len = len_trim (bfactors_header%filename)

   bfactors_header%formatted = bfactors_header%filename(len-2:len) /= 'plt'

   call clone_header (damage_header, 1, 0, bfactors_header)  ! 1 function, POINT

   bfactors_header%title      = 'Heat flux bump factors for ' // &
                                 trim (damage_header%filename)
   bfactors_header%varname(4) = 'Bump_factor_(local_qw)'

   allocate (bfactors(nb_damage))

   do ib = 1, nb_damage

      call clone_zone (damage(ib), 3, 1, bfactors(ib), ios)

      bfactors(ib)%x = damage(ib)%x
      bfactors(ib)%y = damage(ib)%y
      bfactors(ib)%z = damage(ib)%z

   end do

!  Construct a search tree from all smooth OML surface patches:

   nquad = 0
   do ib = 1, nb_smooth
      nquad = (smooth(ib)%ni - 1) * (smooth(ib)%nj - 1) + nquad
   end do

   allocate (conn(3,nquad)) ! For patch # and (i,j)

   call build_adt (nb_smooth, smooth, nquad, conn)

!  Perform the mapping and conversion to bump factors vs. local smooth flow:

   call map_to_smooth ()

   deallocate (conn)

   call deallocate_header (damage_header);  numq = damage_header%numq
   call deallocate_blocks (1, nb_damage, 3, numq, damage, ios)

   call deallocate_header (smooth_header);  numq = smooth_header%numq
   call deallocate_blocks (1, nb_smooth, 3, numq, smooth, ios)

!  Save the bump factors:

   call Tecplot_write (lun_out, bfactors_header, bfactors, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble writing the bump factor results.'
   end if

   call deallocate_header (bfactors_header)

   call deallocate_blocks (1, nb_damage, 3, 1, bfactors, ios)

99 continue

!  Local procedures for program bump_factors

   contains

!     --------------------------------------------------------------------------

      subroutine map_to_smooth ()

!     Map the damage grid to the smooth grid and convert the damage flow to
!     bump factors by dividing by the local smooth flow values.

!     --------------------------------------------------------------------------

!     Local variables:

      integer :: i, ib, ib1, ic, ier, iquad, j, jc, m, n, ni, ninside, nj, &
                 noutside
      real    :: dmax, dmean, dsqmin, dtolsq, p, pm1, q, qm1, smooth_qw,   &
                 interp_xyz(3), target_xyz(3)

!     Execution:

!     Tolerance for search diagnostics (refine later):

      dtolsq = (0.001) ** 2

      do ib = 1, nb_damage

         ninside  = 0;  dmax  = zero
         noutside = 0;  dmean = zero

         ni = damage(ib)%ni
         nj = damage(ib)%nj

         do j = 1, nj

            do i = 1, ni

               target_xyz(1) = damage(ib)%x(i,j,1)
               target_xyz(2) = damage(ib)%y(i,j,1)
               target_xyz(3) = damage(ib)%z(i,j,1)

               call search_adt (target_xyz, iquad, p, q, dsqmin, true,         &
                                nb_smooth, smooth, nquad, conn, interp_xyz)

               if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
                  ninside  = ninside + 1
               else
                  noutside = noutside + 1
               end if

               dmax  = max (dmax, dsqmin)
               dmean = dmean + dsqmin

!              Interpolate the surface flow at this target point:

               n   = conn(1,iquad) ! Block #
               ic  = conn(2,iquad) ! Lower left quad. indices
               jc  = conn(3,iquad)
               pm1 = one - p
               qm1 = one - q

               smooth_qw = &
                  qm1 * (pm1 * smooth(n)%q(iq_smooth,ic,  jc,  1)  + &
                           p * smooth(n)%q(iq_smooth,ic+1,jc,  1)) + &
                    q * (pm1 * smooth(n)%q(iq_smooth,ic,  jc+1,1)  + &
                           p * smooth(n)%q(iq_smooth,ic+1,jc+1,1))

               bfactors(ib)%q(1,i,j,1) = damage(ib)%q(iq_damage,i,j,1) / &
                                         smooth_qw
               if (ib == 8) then
                  if (i ==17 .and. j == 1) then
                     write (6, '(a, 3f20.8)') ' target xyz: ', target_xyz
                     write (6, '(a, 3f20.8)') ' interp xyz: ', interp_xyz
                     write (6, '(a, 3i5)') ' b19,i7,j24: n, ic, jc: ', n, ic, jc
                     write (6, '(4f20.8)') &
                        smooth(n)%q(iq_smooth,ic:ic+1,jc:jc+1,1)
                     write (6, '(a, 3f20.8)') ' p, q, finterp: ', p,q, smooth_qw
           write (6, '(a, f20.8)') ' Damage qw: ', damage(ib)%q(iq_damage,i,j,1)
               write (6, '(a, f20.8)') ' Bump factor: ', bfactors(ib)%q(1,i,j,1)
                  end if
               end if

            end do ! Next i for this target block

         end do ! Next j for this target block

         write (luncrt, '(a, i3, a, 2i7, a, 1p, 2e12.5)')                      &
            '   damage patch', ib,                                             &
            ':  # points in/outside tolerance:', ninside, noutside,            &
            ';  max/mean devn.:', sqrt (dmax), sqrt (dmean) / real (ni * nj)

      end do ! Next target block

      end subroutine map_to_smooth

   end program bump_factors
