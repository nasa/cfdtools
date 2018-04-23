!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program gestimate  ! Gamma estimate
!
!  Description:
!
!     For the given set(s) of free stream density, speed, and temperature,
!     estimate the value(s) of gamma (ratio of specific heats) appropriate for
!     introducing real gas effects to a perfect gas solver such as CART3D.
!
!     Subroutine TGAS1 supplied by Dinesh Prabhu is used to solve two nonlinear
!     equations for e and rho, derived by Dinesh as follows:
!
!        rhoinf Rinf Tinf + rhoinf uinf**2 [1 -  rhoinf/rho]  =  p           (1)
!        G Rinf Tinf + 0.5 [1 - (rhoinf/rho)**2] uinf**2      =  e + p/rho   (2)
!
!     where  e = energy (J/kg)
!            p = pressure (Pa)
!          rho = density (kg/m^3),
!            R = gas constant (J/kg)
!            T = temperature (deg K)
!            G = ginf / (ginf-1); ginf = 1.4
!           Cv = specific heat at constant volume (J/kg)
!            u = speed (m/s)
!
!     Then the value of gamma behind the shock is given by
!
!        gamma = 1 + p / (rho e)
!
!     Starting guesses for e and p have also been supplied by Dinesh using
!     Mathematica.  Inputs and outputs may be in SI or English units.
!     N.B.: English units pressures are in slug/ft**3, not lb/ft**3.
!
!  Sample control file (standard input):
!
!     Control file for program GESTIMATE
!     SI Units?                ! T means kg, Joule, Pa,  m/s, deg K, sq.m.
!     T                        ! F means slug, BTU, psf, fps, deg R, sq.ft.
!     rhoinf         uinf      Tinf
!     3.31773E-4   5144.11     242.62
!     4.58023E-4   4743.08     249.47
!     6.21805E-4   4334.90     255.57
!      :               :          :
!
!  Sample results (standard output):
!
!       rhoinf       uinf   Tinf    Minf   Energy     Density  Press.   M  Gamma
!       kg/m**3       m/s      K            J/kg      kg/m**3     Pa
!     3.31773E-4  5144.11  242.62  16.47  1.123E+07  1.123E-03  1.123 15.1 1.123
!     4.58023E-4  4743.08  249.47  14.98  1.012E+07  2.345E-03  2.345 14.2 1.159
!     6.21805E-4  4334.90  255.57  13.53  9.876E+06  3.456E-03  4.789 13.3 1.234
!      :              :       :      :     :          :          :      :   :
!
!  TGAS1 Reference:  "Simplified curve fits for the thermodynamic properties
!                     of equilibrium air," Srinivasan, S., Tannehill, J.C
!
!  History:
!
!     05/12/06  DKP/DAS  Initial implementation, adapted from THROAT_CONDITIONS
!                        (central difference gradients and a safeguarded Newton
!                        iteration).  This appears to converge well even with
!                        unscaled SI units for the variables, so do without
!                        scaling for now.
!     05/15/06   "   "   The safeguarding against negative energy and density is
!                        surely redundant now. (Slight errors in the formulation
!                        prompted the lower bounds.)
!     05/18/06   "   "   The Cygwin compiler disallows negative logical units.
!                        Use a logical instead of a negative parameter constant.
!
!  Authors:
!
!     Analysis:  Dinesh Prabhu,  ELORET/NASA Ames Research Center, Moffett Field
!     Numerics:  David Saunders, ELORET/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Local constants:

   logical, parameter :: &
      details = .false.     ! .true. shows the iterations

   integer, parameter :: &
      luninp = 5,        &  ! Control file (standard input)
      luncrt = 6,        &
      mflag  = 1,        &  ! Tells TGAS1 to return p and a for given e and rho
      nv     = 2,        &  ! Number of variables being solved for
      itmax  = 40           ! Outer iteration limit
                            ! Halving the step is the inner iteration, using tol

   real, parameter ::    &
      slugft3_per_kgm3 = 515.3788,          &
      kg_per_lb        = 0.45359237,        &
      Joules_per_BTU   = 1054.8,            &
      Pa_per_psf       = 47.8801724,        &
      m_per_ft         = 0.3048,            &
      sqm_per_sqft     = 0.09290304,        &
      ginf             = 1.4,               &
      G                = ginf / (ginf-1.),  &
      Rinf             = 287.06,            & ! Gas constant for air (J/kg)
      Tscale           = 1.8,               &
      half             = 0.5,               &
      one              = 1.0,               &
      third            = one / 3.0,         &
      tol              = 1.E-10,            &
      two              = 2.0,               &
      zero             = 0.0

!  Local variables:

   integer :: &
      i, ier, ios, iter, ip(nv), j, k, line

   real ::    &
      a, e, einf, gamma, M, Minf, p, rho, rho_ratio, rhoinf, T, Tinf, uinf,    &
      cube_root_eps, alpha, fnorm, fnorm0, toler,                              &
      AJ(nv,nv), aitch(nv), dv(nv), f(nv), ffd(-1:1,nv,nv), v(nv),             &
      vfd(nv), vlast(nv), vmin(nv)

   logical :: &
      SI_units

!  Procedures.

   external lusolve  ! Solution of square system by LU-factorization
   external tgas1    ! Gas properties routine (Dinesh Prabhu)

!  Execution:
!  !!!!!!!!!!

   cube_root_eps = epsilon (cube_root_eps) ** third
   line = 0

   read (luninp, *) ! Title
   read (luninp, *)
   read (luninp, *, iostat=ios) SI_units
   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Enter T or F for the units switch.'
      go to 99
   end if
   read (luninp, *)

   do ! Until EOF

      read (luninp, *, iostat=ios) rhoinf, uinf, Tinf
      if (ios < 0) go to 99

      line = line + 1
      if (ios /= 0) then
         write (luncrt, '(/, a, i4)') &
            ' Trouble reading case variables.  Case number:', line
         go to 99
      end if

      if (.not. SI_units) then
         rhoinf = rhoinf * slugft3_per_kgm3
         uinf   = uinf   * m_per_ft
         Tinf   = Tinf   / Tscale
      end if

!     Starting guesses for energy e and density rho:

      Minf = uinf / sqrt (ginf * Rinf * Tinf)

      call starting_guess ()  ! Internal procedure below

!     Perform a Newton iteration for e and rho behind the shock.
!     The iteration involves solving  J dv = f, where J is the matrix
!     of partial derivatives df(i)/dv(j), updating v as v - dv, and
!     repeating until || f || becomes arbitrarily small.

!     Safeguarding involves halving each pure Newton step until || f ||
!     is reduced from its value for the previous iteration, meaning that
!     the update step is actually v <-- v - alpha dv, where alpha <= 1.
!
!     Tests for convergence are subject to the choice of tolerance, which
!     is dependent upon the data scaling.  The tolerance chosen (1.E-8 for
!     64-bit precision) is used relative to the norm of the starting guess
!     in SI units.

      ier      = 0
      iter     = 0
      v(1)     = e              ! Starting guesses for the variables
      vmin(1)  = e * 1.e-5
      v(2)     = rho
      vmin(2)  = rho * 1.e-3
      aitch(:) = (one + v(:)) * cube_root_eps
      dv(:)    = zero           ! For iteration 0 printout
      alpha    = zero           !  "      "     "     "
      fnorm0   = 1.E+30         ! Big to avoid iteration 0 test
      fnorm    = max (abs (v(1)), abs (v(2)))
      toler    = tol * fnorm

      if (line == 1 .or. details) then
         write (luncrt, '(/, a, 1p, 2e15.6, a)') &
            ' Initial energy & density:', v, '   (SI units)'
         write (luncrt, '(a, 1p, 2e15.6)') &
            ' Central difference steps:', aitch
         write (luncrt, '(a, 1p, e15.6)') &
            ' Convergence tolerance:   ', toler
      end if

30    continue

!        Evaluate the functions at the current energy and density.
!        Actually, to calculate central derivatives while we're at it,
!        we evaluate the functions for each v at v - h, v, v + h:

         do j = 1, nv

            vfd(:) = v(:)

            do k = -1, 1

               vfd(j) = v(j) + real (k) * aitch(j)

               vfd(j) = max (vfd(j), vmin(j)) ! Must not go negative

               call tgas1 (vfd(1), vfd(2), p, a, T, mflag)

               if (details) write (luncrt, '(a, 1p, 2e14.6, a, e14.6)') &
                  '   e, rho input to TGAS1:', vfd(:), '   p:', p

!              f(1) = rhoinf Rinf Tinf + rhoinf uinf**2 [1 -  rhoinf/rho]  -  p
!              f(2) = G Rinf Tinf  +  0.5 [1 - (rhoinf/rho)**2] uinf**2
!                     - e - p/rho

               rho_ratio  = rhoinf / vfd(2)
               ffd(k,1,j) = rhoinf * Rinf * Tinf + rhoinf * uinf**2 *          &
                            (one - rho_ratio) - p
               ffd(k,2,j) = G * Rinf * Tinf + half * (one - rho_ratio**2) *    &
                            uinf**2  -  vfd(1)  -  p / vfd(2)
            end do

!           Partial derivatives dfi / dvj:

            do i = 1, nv
              AJ(i,j) = (ffd(1,i,j) - ffd(-1,i,j)) / (two * aitch(j))
            end do

         end do

!        Set up the residual vector f of function values as the RHS of
!        the system J dv = f.  The norm of f should converge to ~0.

         f(1) = ffd(0,1,1)
         f(2) = ffd(0,2,2)

         fnorm = zero
         do i = 1, nv
            fnorm = max (fnorm, abs (f (i)))
         end do

!        Halve the step until || f || is reduced (except first time through).

         if (fnorm > fnorm0) then
            if (alpha > tol) then   ! Presume tol > 2 * machine epsilon
               alpha = half * alpha
               do i = 1, nv
                  v(i) = vlast(i) - alpha * dv(i)
               end do
               v(1) = max (v(1), vmin(1)) ! Must not go negative
               v(2) = max (v(2), vmin(2))
               go to 30
            end if
            go to 83
         end if

         if (details) then
            dv = alpha * dv
            write (luncrt, '(a, i3, a, 1p, e9.2, a, e9.2, 2(a, 2e16.8))') &
               ' Itn.', iter, '  ||f||:', fnorm, '  step:', alpha,        &
               '  dv:', dv, '  v:', v
         end if

         if (fnorm < toler) go to 90  ! Success


         iter = iter + 1
         if (iter == itmax) go to 81

         fnorm0 = fnorm

!        Set up the LHS matrix for the next iteration.  AJ(i,j) = dfi/dvj.

!!!      AJ(1, 1) = ...  ! Already done above during the finite differencing
!!!      AJ(2, 1) = ...
!!!      AJ(1, 2) = ...
!!!      AJ(2, 2) = ...

!        Solve  J dv = f.

         call lusolve (nv, nv, AJ, f, ier)      ! Solution overwrites RHS

         if (ier /= 0) go to 84

         do i = 1, nv
            dv(i)    = f(i)
            vlast(i) = v(i)
            v(i)     = v(i) - dv(i)
         end do

         v(1)  = max (v(1), vmin(1))     ! Must not go negative
         v(2)  = max (v(2), vmin(2))
         alpha = one

      go to 30                           ! Do another iteration


!     Error handling.

81    ier = 1
      write (luncrt, '(/, a)') ' Iteration limit reached.'
      if (fnorm < 10. * toler) go to 90

      ier = 2
!!!   go to 99   ! Bad data scaling can prevent the convergence test from
      go to 90   ! being satisfied even when the solution is still adequate.

83    ier = 3
      write (luncrt, '(/, a)') ' Step-halving iteration failed.'
      go to 90

84    ier = 4
      write (luncrt, '(/, a)') ' Singular matrix.'  ! Show solution anyway
!!!   go to 99

90    continue    ! Wrap up, having converged or almost converged.

!     Ensure the final iterate values are used:

      e   = v(1)
      rho = v(2)

      call tgas1 (e, rho, p, a, T, mflag)

      gamma = one + p / (e * rho)
      M     = (rhoinf / rho) * (uinf / a)

      if (line == 1 .or. details) then

         write (luncrt, '(/, 2a)') &
            '      Rhoinf          uinf           Tinf         Minf      ',    &
            ' Energy        Density       Pressure       M       Gamma'

         if (SI_units) then
            write (luncrt, '(2a)') &
               '      kg/m**3          m/s         Degrees K             ',    &
               '  Joules/kg       kg/m**3          Pa'
         else
            write (luncrt, '(2a)') &
               '    slug/ft**3        ft/s         Degrees R             ',    &
               '    BTU/lb      slug/ft**3         psf'
         end if

      end if

      if (.not. SI_units) then
         rhoinf = rhoinf / slugft3_per_kgm3
         uinf   = uinf   / m_per_ft
         Tinf   = Tinf   * Tscale
         v(1)   = v(1)   * (kg_per_lb / Joules_per_BTU)
         v(2)   = v(2)   / slugft3_per_kgm3
         p      = p      / Pa_per_psf
      end if

      write (luncrt, '(1p, 3e15.6, 0p, f9.4, 1p, 3e15.6, 0p, 2f9.4)') &
         rhoinf, uinf, Tinf, Minf, v(:), p, M, gamma

   end do ! Next case

99 continue

!  Internal procedures for program gestimate:

   contains

      subroutine starting_guess ()

!     Code fragment from Mathematica/Dinesh Prabhu, in SI units using gamma 1.2:

      real :: e_ratio, M2, M4, rho_ratio, term

      einf = Rinf * Tinf / (ginf - one)
      M2   = Minf**2
      M4   = Minf**4
      term = sqrt (0.29387755102040813 - 0.05714285714285694*M2 + 0.4*M4)

      e_ratio = &
         (-0.9333333333333333 + 6.346666666666667*M2 + 1.3066666666666667*M4)/ &
         (-0.4 + 0.6533333333333333*M4 - 0.7378647873726215*term +             &
            M2 * (6.066666666666666667 - 1.0330107023216701*term))

      e = e_ratio * einf

      rho_ratio = &
         (0.5421047417431506*Minf*sqrt (-0.4 + 2.8*M2) *                       &
          sqrt (2. + (1.361111111111111*(0.29387755102040813 +                 &
                0.34285714285714297* M2  + 0.48*M4 -                           &
                0.45175395145262554*(1.2 + 1.68*M2)*term))/(-0.4 + 2.8*M2))) / &
         (sqrt (2. + 0.4*M2) *                                                 &
          sqrt (0.29387755102040813 + 0.34285714285714297*M2 + 0.48*M4 -       &
                0.45175395145262554*(1.2 + 1.68*M2) *                          &
                sqrt (0.29387755102040813 - 0.05714285714285694*M2 + 0.4*M4)))

      rho = rho_ratio * rhoinf

      end subroutine starting_guess

   end program gestimate
