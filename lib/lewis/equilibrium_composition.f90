!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine equilibrium_composition (specified_state, state_1, state_2,      &
                                       nspecies, names, mole_fractions,        &
                                       pressure, density, temperature,         &
                                       enthalpy, entropy, speed_of_sound,      &
                                       status)

!  Calculate the equilibrium composition of a gas mixture.  Two of the mixture
!  state variables are given along with a list of species names and starting
!  guesses for their mole fractions.  The inputs should be consistent with the
!  C implementation by Gary Allen of the Gordon and McBride technique for
!  performing Gibbs free energy minimization.
!
!  All input and output values should be in SI units.
!
!  12/10/05  D. A. Saunders  Initial design and implementation, employing a
!                            system call for the C program that drives the
!                            "lewis" procedure via an input control file written
!                            here.  The resulting output file is parsed here.
!                            NOZZLE_THROAT_CONDITIONS is the first application.
!
!  10/08/08    "        "    Avoid the system call by passing the input data and
!                            output results as arguments to an interface routine
!                            (equilibrium_gas) written in C.  This version is
!                            intended for the [NOZZLE_]THROAT_CONDITIONS[_3D]
!                            pair and for the DPLR flow solver (potentially).
!
!  10/29/08    "        "    Made use of status argument to control printing.
!
!  09/06/11    "        "    Gary supplied an updated lewis.h, lewis.c, etc.;
!                            max_species = 54 now, not 48 (some insertions).
!
!  09/12/11    "        "    Leave any error message to the higher level, as
!                            opposed to supplying a logical unit number to
!                            write to (or writing to standard output), because
!                            of apparent clashes between use of standard output
!                            by both Fortran 90 and C in the same application.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!                      now: ERC, Inc./NASA ARC
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character, intent (in) :: specified_state * 5 ! Descriptor defining the two
                                                 ! state variables input:
                                                 ! Rho_T, Rho_H, Rho_S,
                                                 ! P_T, P_H, P_S are the options

   real,      intent (in) :: state_1, state_2    ! The two input state variables

   integer,   intent (in) :: nspecies            ! # gas species

   character, intent (in) :: names(nspecies) * 4 ! Species names consistent with
                                                 ! Gary's lewis.h header file

   real,   intent (inout) :: mole_fractions  &   ! Input with starting guesses;
                             (nspecies)          ! output with equilib. values

   real,     intent (out) :: pressure,       &   ! Output mixture properties
                             density,        &
                             temperature,    &
                             enthalpy,       &
                             entropy,        &
                             speed_of_sound      ! What else?

   integer, intent (inout):: status              ! On input:  0 => no printing;
                                                 !            1 => print;
                                                 ! On output, 0 => no problem
!  Local constants:

   integer, parameter :: &
      max_species = 54,  &                       ! Same as N_species in lewis.h
      max_states  = 6

!  Local variables:

   integer   :: &
      i, ier, ispecies(nspecies), istate, n
   real      :: &
      sum
   character, save :: &
      species(0:max_species-1)*4, states(0:max_states-1)*5

!  Procedures:

   integer  :: equilibrium_gas
   external :: equilibrium_gas  ! C interface to lewis.c

!  Storage:

   data species &  ! Must be consistent with #define Species_list in lewis.h
     /'N   ',   &
      'C   ',   &
      'Ar  ',   &
      'H   ',   &
      'He  ',   &
      'O   ',   &
      'El  ',   &
      'Cgr ',   &
      'H2  ',   &
      'NO  ',   &
      'N2  ',   &
      'O2  ',   &
      'H2O ',   &
      'Wtr ',   &
      'Ice ',   &
      'OH  ',   &
      'CO  ',   &
      'CO2 ',   &
      'CH4 ',   &
      'CN  ',   &
      'C2N2',   &
      'HCN ',   &
      'HNC ',   &
      'C2  ',   &
      'C3  ',   &
      'C4  ',   &
      'C5  ',   &
      'C2H ',   &
      'CH  ',   &
      'CH2 ',   &
      'CH3 ',   &
      'C2H2',   &
      'Vin ',   &
      'C2H3',   &
      'C3H3',   &
      'C4H2',   &
      'C6H2',   &
      'C2O ',   &
      'Np  ',   &
      'Op  ',   &
      'Cp  ',   &
      'Nn  ',   &
      'On  ',   &
      'Cn  ',   &
      'Arp ',   &
      'NOp ',   &
      'N2p ',   &
      'O2p ',   &
      'COp ',   &
      'Hp  ',   &
      'Hep ',   &
      'Npp ',   &
      'Opp ',   &
      'Cpp '/

   data states  &
     /'Rho_T',  &
      'P_H  ',  &
      'P_T  ',  &
      'Rho_H',  &
      'Rho_S',  &
      'P_S  '/

!  Execution:
!  !!!!!!!!!!

!  Avoid passing character data from Fortran to C (though there may be a way).
!  If a bad string is passed in, Gary's CEA code will catch it, so don't repeat
!  the user error traps here.  (The driving program should be trapping bad
!  character strings.)

   do n = 1, nspecies
      do i = 0, max_species - 1
         if (names(n) == species(i)) then
            ispecies(n) = i
            exit
         end if
      end do
   end do

   do i = 0, max_states - 1
      if (specified_state == states(i)) then
         istate = i
         exit
      end if
   end do

!  Calculate the equilibrium mixture via this C interface to lewis.c:

!! These diagnostics have been known to disappear when they were written to
!! unit 6.  In turn, they seemed to interfere with singular matrix messages
!! from lewis.c  Unit 7 matches the log file now used by the nozzle codes.

!! write (7, '(a, a)') 'specified_state: ', specified_state
!! write (7, '(a, 1p, 2e20.10)') 'state_1, state_2: ', state_1, state_2
!! write (7, '(a, i4)') 'istate:', istate
!! write (7, '(a, i4)') 'status input:', status
!! write (7, '(a, i4)') 'nspecies:', nspecies
!! write (7, '(i2, 2x, a, i3)') (n, names(n), ispecies(n), n = 1, nspecies)

   ier = equilibrium_gas (status,   istate,   state_1,  state_2,            &
                          nspecies, ispecies, mole_fractions,               &
                          pressure, density,  temperature,                  &
                          enthalpy, entropy,  speed_of_sound)
   status = ier
   if (status /= 0) go to 99

!  Normalize the mole fractions to avoid a check in Gary's CEA code:

   sum = 0.
   do n = 1, nspecies
      sum = sum + mole_fractions(n)
   end do
   mole_fractions(:) = mole_fractions(:) / sum

99 return

   end subroutine equilibrium_composition
