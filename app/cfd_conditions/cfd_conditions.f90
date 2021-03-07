!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program cfd_conditions
!
!  Description:
!
!     This utility converts between (V, h) and (M, Q) in either direction, and
!     tabulates corresponding densities, temperatures, and sound speeds.  Here,
!     V = freestream velocity, h = altitude, M = Mach number and Q = dynamic
!     pressure, and the atmosphere model is the U.S. 1976 Standard Atmosphere
!     model for air.  SI units are output; input V = m/s|km/s; input h = m|km.
!
!     Choosing additional aerothermal database points for Orion/MPCV prompted
!     the capability.
!
!     A related requirement is to choose such CFD points in normalized space
!     then convert those coordinates to unnormalized form.  This has been
!     retrofitted as a rather trivial separate capability.  The remaining
!     description applies to the original functionality.
!
!     Most recently, for normalizing V & h, in addition to using the input
!     dataset coordinates for the scales and shifts, the user may now enter
!     different scales and shifts, probably from a prior run of this utility.
!     The intent is first to normalize a list of points comprising a database,
!     then to apply those factors to a trajectory dataset so that the two may
!     be overlaid in a plot.  This enables plotting database points that are
!     being proposed before they are actually in hand, in normalized space.
!     The author's radio blackout software does its interpolations in normalized
!     space, so the points the interpolation will be working with can be shown
!     before the necessary flow calculations have been performed.
!
!     While M and Q can be derived from V and h explicitly, going the other way
!     involves solving a nonlinear equation in h as follows:
!
!            M = V/a                    ! a = speed of sound = sqrt (gamma Rs T)
!        =>  V = M a                    ! gamma = specific heat ratio, 1.4
!            Q = 0.5 rho V**2           ! Rs = specific gas constant
!              = 0.5 rho M**2 gamma Rs T
!        =>  rho T = 2Q / (M**2 Rs)
!
!     From the atmosphere model, rho = rho(h) and T = T(h), so we need to find
!     the altitude h* at which rho(h*) T(h*) = rho T for the given M, Q.  Non-
!     derivative zero finder ZERORC is appropriate.  The solution is assumed to
!     be unique in a sensible altitude range.
!
!        Solve  rho(h) T(h) - rho T = 0  for h = h*
!        Then   V* = M a = M sqrt (gamma Rs T(h*))
!
!     The input table may have any number of header lines (ignored) and any
!     number of columns.  The columns for V & h or M & Q are prompted for.  Any
!     other columns are ignored.  Meters & meters/sec are the working units.
!
!     A radio attenuation database requirement has prompted an option to include
!     normalized V & h in the output table, which looks something like this:
!
!        V,km/s h,km  rho,kg/m^3     T,K   p,Pa   a,m/s   Mach    Q,Pa  [Vn  hn]
!        8.2    68.0  1.09189e-4 225.065 xxx.xx  xxx.xx  xx.xx  xxx.xx  [xx  xx]
!         :       :    :            :       :       :      :       :    [ :  : ]
!
!     The author's radio blackout software uses a triangulation of normalized V
!     & h to identify target (V, h) points that are outside the database range.
!
!     In anticipation of possibly switching to use of (M, Q) space, the option
!     to normalize the coordinates now allows for normalized M & Q rather than
!     normalized V & h if desired.  A "table_delaunay" utility can produce the
!     associated triangulation.
!
!  History:
!
!     01/26/17  D.A.Saunders  Adaptation of Grant Palmer's 1976 U.S. Standard
!                             Atmosphere model as a resuable subroutine named
!                             earth_atmosphere (rho, T, p, a for given h).
!     01/27/17    "   "   "   Initial implementation of CFD_Conditions, for MPCV
!                             aerothermal & radio attenuation database purposes.
!     02/07/17    "   "   "   If output of normalized coordinates is requested,
!                             they can now be for either V & h or M & Q.
!     03/17/17    "   "   "   Retrofitted a way of denormalizing a list of
!                             normalized (V, h) (or (M, Q)) points.  This was
!                             prompted by some radio blackout software by the
!                             present author, which can now output the target
!                             trajectory coordinates in normalized form for
!                             overlaying on the normalized database coordinates
!                             that are used during interpolation in 2-space.
!     03/31/17    "   "   "   Provided for normalizing a trajectory using shifts
!                             and scales produced by a prior run on a set of
!                             [proposed?] radio attenuation database pts. (say).
!       
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use table_type_module    ! Both are in table_io.f90
   use table_io

   implicit none

!  Constants:

   integer, parameter :: &
      icol_V_out = 1,    &  ! Output column for velocity
      icol_h_out = 2,    &  !    "     "     "  altitude
      icol_rho   = 3,    &  !    "     "     "  density
      icol_T     = 4,    &  !    "     "     "  temperature
      icol_p     = 5,    &  !    "     "     "  pressure
      icol_a     = 6,    &  !    "     "     "  speed of sound
      icol_M_out = 7,    &  !    "     "     "  Mach number
      icol_Q_out = 8,    &  !    "     "     "  dynamic pressure
      icol_Vnorm = 9,    &  !    "     "     "  normalized velocity
      icol_hnorm = 10,   &  !    "     "     "     "     " altitude
      lunin      = 1,    &  ! Input table
      lunout     = 2,    &  ! Output table
      lunkbd     = 5,    &  ! Keyboard inputs
      luncrt     = 6        ! Screen outputs

   real, parameter ::    &
      gamma  = 1.4,      &  ! Gamma, ratio of specific heats for air
      half   = 0.5,      &
      Rs     = 287.058,  &  ! Specific gas constant, J/kg/K
      thous  = 1000.

!  Variables:

   integer :: &
      i, icol_h, icol_V, icol_M, icol_Q, ier, l, ncols_out, nrows

   real :: &
      a, h_km, M, p, Q, rho, T, V_mps, vmean, vsd

   logical :: &
      cr, dataset_mean_sd, denorm, eof, full_precision, h_in_km, norm, &
      V_in_kmps, Vh_input, Vh_norm

   type (table_type) :: &
      table_in, table_out

!  Execution:

   denorm = .false.  ! Retrofitted option
   call ready (luncrt, &
               'Simply denormalize given coordinates? (y); ' // &
               '(V, h) <--> (M, Q) or normalization (n); <cr>=n: ', &
               lunkbd, denorm, cr, eof)
   if (eof) go to 99

   if (denorm) then
      call denormalize ()  ! Keep it well separated at some duplication cost
      go to 99
   end if

   Vh_input = .true.
   call ready (luncrt, &
               'Are (velocity, altitude) to be input? (y); ' // &
               '(Mach, Q) input? (n); [y] ', &
               lunkbd, Vh_input, cr, eof)
   if (eof) go to 99

   call reads (luncrt, 'Input table (any # columns & header rows): ', &
               lunkbd, table_in%filename, cr, eof)
   if (eof) go to 99

   call table_io_read_real (lunin, table_in, ier)
   if (ier /= 0) go to 99

   if (Vh_input) then
      icol_V = 1
      call readi (luncrt, 'Input column for V: [1] ', lunkbd, icol_V, cr, eof)
      if (eof) go to 99

      icol_h = 2
      call readi (luncrt, 'Input column for h: [2] ', lunkbd, icol_h, cr, eof)
      if (eof) go to 99

      V_in_kmps = .true.
      call ready (luncrt, 'Input velocities in km/s? (y);  m/s (n): [y] ', &
                  lunkbd, V_in_kmps, cr, eof)
      if (eof) go to 99

      ! Ensure that the working velocities are m/s:

      if (V_in_kmps) table_in%values(icol_V,:) = table_in%values(icol_V,:)*thous

      h_in_km = .true.
      call ready (luncrt, 'Input altitudes in km? (y);  meters (n): [y] ', &
                  lunkbd, h_in_km, cr, eof)
      if (eof) go to 99

      if (.not. h_in_km) &  ! Ensure working altitudes are km:
         table_in%values(icol_h,:) = table_in%values(icol_h,:)/thous

   else ! M, Q are input
      icol_M = 1
      call readi (luncrt, 'Input column for M: [1] ', lunkbd, icol_M, cr, eof)
      if (eof) go to 99

      icol_Q = 2
      call readi (luncrt, 'Input column for Q: [2] ', lunkbd, icol_Q, cr, eof)
      if (eof) go to 99

   end if

!  Option to output normalized coordinates:

   norm = .true.
   call ready (luncrt, 'Add normalized coordinates to the output? [y] ', &
               lunkbd, norm, cr, eof)
   if (eof) go to 99

   if (norm) then
      Vh_norm = .true.
      call ready (luncrt, 'Normalized V & h? (y);  M & Q? (n): [y] ', &
                  lunkbd, Vh_norm, cr, eof)
      if (eof) go to 99
      if (Vh_norm) then
         dataset_mean_sd = .true.
         call ready (luncrt, &
                     'Use input dataset means & sigmas? n = enter them [y] ', &
                     lunkbd, dataset_mean_sd, cr, eof)
         if (eof) go to 99
      end if
   end if

10 call reads (luncrt, 'Output table with V, h, rho, T, p, a, Mach, Q: ', &
               lunkbd, table_out%filename, cr, eof)
   if (eof) go to 99
   if (cr)  go to 10

!  Set up the output table:

   nrows  = table_in%nrows
   table_out%nrows = nrows
   ncols_out = icol_Q_out
   if (norm) ncols_out = ncols_out + 2
   table_out%ncols = ncols_out
   table_out%column_type(:) = 2  ! All reals
   table_io_format = &
      '(f10.4, f11.4, es15.7, f10.4, es15.7, f11.4, f9.4, f12.4, 2f10.6)'
   table_out%nheader = 1
   allocate (table_out%header(1))
   table_out%header(1) = &
      '    V,km/s       h,km     rho,kg/m^3       T,K           p,Pa' // &
      '      a,m/s     Mach        Q,Pa'
   l = len_trim (table_out%header(1))

   if (norm) then
      if (Vh_norm) then
         table_out%header(1)(l+1:l+20) = '     Vnorm     hnorm'
      else
         table_out%header(1)(l+1:l+20) = '     Mnorm     Qnorm'
      end if
   end if

   allocate (table_out%values(ncols_out,nrows))

!  Do the calculations for each input condition:

   if (Vh_input) then
      do i = 1, nrows
         V_mps = table_in%values(icol_V,i)
         h_km  = table_in%values(icol_h,i)

         call earth_atmosphere (h_km, rho, T, p, a)

         table_out%values(icol_V_out,i) = V_mps/thous  ! Only km/s out
         table_out%values(icol_h_out,i) = h_km         ! Only km   out
         table_out%values(icol_rho,  i) = rho
         table_out%values(icol_T,    i) = T
         table_out%values(icol_p,    i) = p
         table_out%values(icol_a,    i) = a
         table_out%values(icol_M_out,i) = V_mps/a
         table_out%values(icol_Q_out,i) = half*rho*V_mps**2
      end do

   else  ! (Mach, Q) have been input

      do i = 1, nrows
         M = table_in%values(icol_M,i)
         Q = table_in%values(icol_Q,i)

         call vel_alt_from_Mach_Q (M, Q, V_mps, h_km, rho, T, p, a, ier)
         if (ier /= 0) then
            write (luncrt, '(a, i3, a, f9.4, es15.7)') &
               'ZERORC trouble: ier =', ier, ';  Mach, Q:', M, Q
            go to 99
         end if

         table_out%values(icol_V_out,i) = V_mps/thous  ! Only km/s out
         table_out%values(icol_h_out,i) = h_km         ! Only km   out
         table_out%values(icol_rho,  i) = rho
         table_out%values(icol_T,    i) = T
         table_out%values(icol_p,    i) = p
         table_out%values(icol_a,    i) = a
         table_out%values(icol_M_out,i) = M
         table_out%values(icol_Q_out,i) = Q
      end do
   end if

   if (norm) then
      if (Vh_norm) then  ! Normalize V and h:

         if (dataset_mean_sd) then
            call stdev1p (nrows, table_out%values(icol_V_out,:), 1, vmean, vsd)
            write (luncrt, '(/, a, 2es16.8)') &
               'Mean and standard deviation for V:', vmean, vsd
         else
            write (luncrt, '(a)', advance='no') ' Mean & sigma for "V": '
            read  (lunkbd, *) vmean, vsd
         end if
         table_out%values(icol_Vnorm,:) = &
            (table_out%values(icol_V_out,:) - vmean)/vsd

         if (dataset_mean_sd) then
            call stdev1p (nrows, table_out%values(icol_h_out,:), 1, vmean, vsd)
            write (luncrt, '(a, 2es16.8, /)') &
               'Mean and standard deviation for h:', vmean, vsd
         else
            write (luncrt, '(a)', advance='no') ' Mean & sigma for "h": '
            read  (lunkbd, *) vmean, vsd
         end if
         table_out%values(icol_hnorm,:) = &
            (table_out%values(icol_h_out,:) - vmean)/vsd

      else  ! Normalize M and Q:

         call stdev1p (nrows, table_out%values(icol_M_out,:), 1, vmean, vsd)
         write (luncrt, '(/, a, 2es16.8)') &
            'Mean and standard deviation for M:', vmean, vsd
         table_out%values(icol_Vnorm,:) = &
            (table_out%values(icol_M_out,:) - vmean)/vsd
         call stdev1p (nrows, table_out%values(icol_Q_out,:), 1, vmean, vsd)
         write (luncrt, '(a, 2es16.8, /)') &
            'Mean and standard deviation for Q:', vmean, vsd
         table_out%values(icol_hnorm,:) = &
            (table_out%values(icol_Q_out,:) - vmean)/vsd

      end if
   end if

!  Save the results:

   call table_io_write_real (lunout, table_out, table_io_format_sep, ier)
!!!if (ier /= 0) go to 99

   full_precision = .false.
   call ready (luncrt, &
              'Full precision output for reversibility check? [y|n; <cr>=n] ', &
               lunkbd, full_precision, cr, eof)
   if (full_precision) then
      do i = 1, nrows
         write (luncrt, '(10es23.15)') table_out%values(:,i)
      end do
   end if

99 continue

   contains

!     Internal routine for program cfd_conditions:

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine denormalize ()  ! Retrofitted option
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call reads (luncrt, 'Input table (any # columns & header rows): ', &
                  lunkbd, table_in%filename, cr, eof)
      if (eof) go to 99

      call table_io_read_real (lunin, table_in, ier)
      if (ier /= 0) go to 99

      icol_V = 1
      call readi (luncrt, 'Input column for "V": [1] ', lunkbd, icol_V, cr, eof)
      if (eof) go to 99

      write (luncrt, '(2a)') 'If only one column is to be denormalized, ', &
                             'enter 0 next.'
      icol_h = 2
      call readi (luncrt, 'Input column for "h": [2] ', lunkbd, icol_h, cr, eof)
      if (eof) go to 99

      table_out%nheader = 1
      table_out%nrows   = table_in%nrows
      table_out%ncols   = 1;  if (icol_h > 0) table_out%ncols = 2
      table_io_format   = '(2f11.6)'
      table_out%column_type(:) = 2  ! All reals

      allocate (table_out%header(1))
      allocate (table_out%values(table_out%ncols,table_out%nrows))

      write (luncrt, '(a)', advance='no') ' Mean & sigma for "V": '
      read  (lunkbd, *) vmean, vsd

      table_out%values(1,:) = vmean + vsd*table_in%values(icol_V,:)
      table_out%header(1) = '        "V"'

      if (icol_h > 0) then
         write (luncrt, '(a)', advance='no') ' Mean & sigma for "h": '
         read  (lunkbd, *) vmean, vsd
         table_out%values(2,:) = vmean + vsd*table_in%values(icol_h,:)
         table_out%header(1)(12:22) = '        "h"'
      end if

   10 call reads (luncrt, 'Output table with "V"[& "h"] denormalized: ', &
               lunkbd, table_out%filename, cr, eof)
      if (eof) go to 99
      if (cr)  go to 10

      call table_io_write_real (lunout, table_out, table_io_format_sep, ier)

   99 return

      end subroutine denormalize

   end program cfd_conditions
