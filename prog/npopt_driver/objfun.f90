!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine objfun (mode, nvpar, vpar, objpar, grdpar, nstate)
!
!     OBJFUN acts as the interface between NPOPT and FUN.  See also CONFUN if
!  nonlinear constraints are present.  Its main task is to determine whether
!  NPOPT wants just the objective function or the objective function plus
!  gradient information, and then call FUN appropriately.
!
!  History:
!  1996      James Reuther/David S.  Initial implementation for SYN87-SB.
!  2000+       "      "      "       HEAT_SHIELD/Traj_opt reuse (much the same).
!  08/07/07  D.A.Saunders            Provided a central differencing option.
!  05/21/12    "      "              NPOPT_DRIVER form (essentially the same).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Global variables:

   use opt_mod

   implicit none

!  Arguments:

   integer, intent (inout) :: mode          ! 0 means return the objective;
                                            ! 1 means return the gradient;
                                            ! 2 means return both;
                                            !-1 output means the objective
                                            !   function calculation was aborted
   integer, intent (in)    :: nvpar         ! Number of optimization variables
   real,    intent (inout) :: vpar(nvpar)   ! Optimization variables
   real,    intent (out)   :: objpar        ! Objective function value
   real,    intent (out)   :: grdpar(nvpar) ! Gradient vector
   integer, intent (in)    :: nstate        ! 1 means this is the first call

!  Local variables:

   integer :: i, ncall, nhalf
   real    :: objfor, objbak, temph, tempv
   logical :: fail

!  Execution:

!  Treat the case in which nonlinear constraints have forced calls to FUN,
!  thus providing OBJPAR and GRDPAR.  See CONFUN.

   if (confun_does_obj) then

      objpar = obj_confun

      if (mode > 0) then

         grdpar(:) = ffwd(:)  ! CONFUN now prints the gradient

!!!      write (lunlog, '(/, 1x, a)') &
!!!         'OBJFUN: Objective gradient via CONFUN:', ' ', &
!!!         '   I            V (I)            G (I)            H (I)'
!!!      write (lunlog, '(i5, 3es17.8)') &
!!!         (i, vpar(i), grdpar(i), h(i), i = 1, nvpar)
      end if

   else

!     Case of no nonlinear constraints:

!!!   write (lunlog, '(/, 2a, 2i4)') 'OBJFUN: MODE and NSTATE: ', mode, nstate

!     Calculate the objective function:

      if (nstate == 1) then ! First OBJFUN call from NPOPT

         ncall = -3

         call fun (nvpar, vpar, objpar, ncall, fail)

         if (fail) then
            write (lunlog, '(/, a)') 'OBJFUN: Bad return from call to FUN.'
            go to 90
         end if

         if (major_iteration_limit == 0) then
            write (lunlog, '(/, a)') 'OBJFUN: Exiting from NPOPT (niter = 0).'
            go to 90
         end if

      else ! Not the first call to OBJFUN

         if (mode == 0) then ! Objective function only
            ncall =  0       ! Line search function evaluation
         else                ! Objective and gradient
            ncall = -2       ! Recover the best line search solution
         end if

         call fun (nvpar, vpar, objpar, ncall, fail)

         if (fail) then
            if (mode == 0) then
               write (lunlog, '(/, a)')'OBJFUN: NPOPT may shorten step & retry.'
               objpar = 9999.
               go to 99
            else
               write (lunlog, '(/, a)') 'OBJFUN: Aborting.'
               go to 90
            end if
         end if

      end if

      if (modE > 0) then ! Calculate the gradient too

         if (mode_gradients == 2) then ! Forward differencing

            do ncall = 1, nvpar
               tempv = vpar(ncall)
               temph = h(ncall)
               nhalf = 0

10             continue

               vpar(ncall) = tempv + temph

               call fun (nvpar, vpar, objfor, ncall, fail)

               if (fail) then
                  write (lunlog, '(/, a, i4)') &
                     'OBJFUN: Function failure for gradient element ', ncall

                  if (nhalf < 3) then
                     nhalf = nhalf + 1
                     temph = temph * 0.5
                     write (lunlog, '(a)') &
                        'Halve the stepsize for this variable.'
                     go to 10
                  end if

                  objfor = objpar
                  write (lunlog, '(/, 2a)') &
                     'Three halvings of gradient step size failed.', &
                     'Proceeding with zero for this component.'
               end if

               grdpar(ncall) = (objfor - objpar) / temph
               vpar(ncall)   = tempv
            end do

         else ! mode_gradients = 3 (central differencing)

            do ncall = 1, nvpar
               tempv = vpar(ncall)
               temph = h(ncall) * epsm6th  ! Since eps ** (1/3) is ~optimal
               vpar(ncall) = tempv + temph

               call fun (nvpar, vpar, objfor, ncall, fail)

               vpar(ncall) = tempv - temph

               call fun (nvpar, vpar, objbak, ncall, fail)

               grdpar(ncall) = (objfor - objbak) / (temph + temph)
               vpar(ncall)   = tempv
            end do

         end if

         write (lunlog, '(/, a)') 'OBJFUN: Gradient vector:', &
            '   I            V (I)            G (I)            H (I)'
         write (lunlog, '(i4, 3es17.8)') &
            (i, vpar(i), grdpar(i), h(i), i = 1, nvpar)

      end if

   end if

   go to 99

90 mode = -1 ! NPOPT should quit

99 return

   end subroutine objfun
