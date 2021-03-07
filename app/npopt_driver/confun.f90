!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine confun (mode, ncpar, nvpar, nrowjpar, needc, vpar, &
                         cpar, cjacpar, nstate)
!
!     CONFUN computes the nonlinear constraint functions C(*) and their
!  gradients.  Its calling sequence is that required by NPSOL/NPOPT.
!  Finite differencing uses the "h"s associated with each variable.
!  The fixed calling sequence requires these and a few other items to
!  be accessed from a module (opt_mod).
!
!     The "par" nomenclature allows certain arguments (parameters) to appear
!  in opt_mod as well (e.g., v(*) and vpar(*) may well be the same array).
!
!     This version ignores needc(*) -  it is more efficient to evaluate all
!  constraints in groups than to evaluate each one independently.
!
!     This version also allows CONFUN to evaluate the objective function and
!  its derivatives.  Otherwise, the optimization repeats the same calls to
!  FUN for the objective and for the nonlinear constraints.
!
!     Since OBJFUN is called following each CONFUN call, it is safe for
!  CONFUN to carry around the most recent objective and its gradient (but
!  OBJFUN and FUN still have to handle the no-nonlinear-constraint case).
!
!  History:
!
!  1996      James Reuther/David S.  Initial implementation for SYN87-SB.
!  2000+       "      "      "       HEAT_SHIELD/Traj_opt reuse (much the same).
!  08/07/07  D.A.Saunders            Provided a central differencing option.
!  05/21/12    "      "              NPOPT_DRIVER form (essentially the same). 
!  05/24/12    "      "              NLCON didn't have access to the variables.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Global variables:

   use opt_mod    ! Optimization variables common to any application

   implicit none

!  Arguments:

   integer, intent (inout) :: &
      mode                    ! mode = 2 means the elements of C(*) indicated
                              !          by NEEDC(*) > 0 need to be set (and
                              !          similarly for the available elements
                              !          of the rows of CJAC(*,*));
                              ! mode = 1 means the available elements of the
                              !          rows of CJAC(*,*) corresponding to
                              !          NEEDC(*) > 0 must be set; C(*) and
                              !          other rows of CJAC will be ignored;
                              ! mode = 0 means the elements of C(*) indicated
                              !          by NEEDC(*) > 0 must be set; CJAC and
                              !          other elements of C(*) are ignored;
                              ! mode =-1 means terminate the solution of the
                              !          current problem.
   integer, intent (in) :: &
      ncpar,               &  ! The number of nonlinear constraints.
      nvpar,               &  ! The number of variables.
      nrowjpar,            &  ! The maximum number of nonlinear constraints.
      needc(ncpar)            ! Indices of the elements of C(*) and CJAC(*,*)
                              ! which must be returned by CONFUN (ignored).
   real, intent (inout) :: &
      vpar(nvpar)             ! Variables at which to evaluate the constraints,
                              ! temporarily perturbed here (but restored).
   real, intent (out) ::   &
      cpar(ncpar)             ! Constraint values, returned if MODE = 0 or 2
                              ! for constraints indicated by NEEDC(*) > 0.
   real, intent (out) ::   &
      cjacpar(nrowjpar,nvpar) ! Constraint derivatives; see MODE.

   integer, intent(in) ::  &
      nstate                  ! nstate = 1 means this is the first call (which
                              ! precedes the first call to OBJFUN).

!  Local constants:

   real,      parameter :: half = 0.5
   logical,   parameter :: noprint = .false.

!  Local variables:

   integer :: i, j, ncall
   real    :: hcen, objfor, objbak, vtemp
   logical :: fail

!  Execution:

!! write (lunlog, '(/, a, 2i4)') ' confun: mode & nstate: ', mode, nstate

   if (nstate == 1) then ! First call to CONFUN (precedes OBJFUN)
      ncall = -3
   else ! Not the first call to CONFUN
      if (mode == 0) then ! Constraints only
         ncall =  0
      else                ! Constraints + derivatives, best line search result
         ncall = -2
      end if
   end if

!  Perform the analysis corresponding to the input variables:

   call fun (nvpar, vpar, objfor, ncall, fail)

   if (fail) then
      write (lunlog, '(/, 1x, a)') 'CONFUN:  Bad return from FUN.'
      go to 90
   end if

   obj_confun = objfor  ! Avoid repeat: OBJFUN uses it

!  Evaluate all nonlinear constraints, scaled:

   call nlcon (nvpar, vpar, noprint, cpar)

   if (major_iteration_limit == 0) then
      if (confun_does_obj) then
         write (lunlog, '(/, 1x, a)') &
            'CONFUN:  Forcing exit from NPOPT (niter = 0).'
         go to 90
      else ! No point in doing the constraint gradients
         cjacpar(:,:) = 999. ! Keep NPOPT happy
         go to 99
      end if
   end if

!  Constraint derivatives as well?

   if (mode > 0) then

      if (mode_gradients == 2) then ! Forward differencing

!        Perturb a given variable only once for all constraints:

         do j = 1, nvpar

            vtemp   = vpar(j)
            vpar(j) = vtemp + h(j)

            call fun   (nvpar, vpar, objfor, j, fail)     ! ncall = j here
            call nlcon (nvpar, vpar, noprint, cjacpar(1,j))

            vpar(j) = vtemp

            do i = 1, ncpar
               cjacpar(i,j) = (cjacpar(i,j) - cpar(i)) / h(j)
            end do

            ffwd(j) = (objfor - obj_confun) / h(j)

         end do

      else ! MODE_GRADIENTS = 3: central differencing

!        Perturb a given variable twice for all constraints:

         do j = 1, nvpar

            vtemp   = vpar(j)
            hcen    = h(j) * epsm6th  ! Since eps**(1/3) is ~optimal
            vpar(j) = vtemp + hcen

            call fun   (nvpar, vpar, objfor, j, fail)
            call nlcon (nvpar, vpar, noprint, cjacpar(1,j))

            vpar(j) = vtemp - hcen

            call fun   (nvpar, vpar, objbak, j, fail)
            call nlcon (nvpar, vpar, noprint, cjacbak)

            vpar(j) = vtemp

            hcen = half / hcen

            do i = 1, ncpar
               cjacpar(i,j) = (cjacpar(i,j) - cjacbak(i)) * hcen
            end do

            ffwd(j) = (objfor - objbak) * hcen

         end do

      end if

      write (lunlog, '(/, a)' ) &
         ' IVAR              VAR     dOBJ/dVAR   dC1/dVAR   dC2/dVAR ...'

      i = min (ncpar, 10)

      do j = 1, nvpar
         write (lunlog, '(i4, es18.9, es14.5, 10es11.3)') &
            j, vpar(j), ffwd(j), cjacpar(1:i,j)
      end do

   end if

   go to 99

90 mode = -1 ! NPOPT should quit

99 return

   end subroutine confun
