!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program prepare_fiat_data
!
!  Description:
!
!     FIAT is a material response solver developed at NASA ARC to help determine
!  the thickness and composition of the thermal protection system (TPS) required
!  by an atmospheric entry vehicle for a given mission.  FIAT treats one body
!  point at a time.  Its normal inputs include a time history of the aerothermal
!  environment predicted by other solvers for that body point.  An available
!  document summarizes 13 steps required to prepare such an environment dataset.
!  These steps are typically required for a handful of body points when a new
!  entry probe mission is proposed.  They are cumbersome enough that automating
!  them to some extent is desirable.  Tying all of them together somehow is not
!  really practical.  Rather, several utilities have been developed to simplify
!  some of the steps.  In its initial form, PREPARE_FIAT_DATA performs just the
!  final step:  it gathers time histories from the curve-fitting step and writes
!  two input environment datasets for FIAT: one in English units as required by
!  FIAT v3, and one in CGS units as expected by Traj_FIAT.
!
!     This version distinguishes between ablating and non-ablating cases - see
!  below.
!
!  Usage Example:
!
!     prepare_fiat_data < prepare_fiat_shoulder_pt.inp
!
!  Control File Example (Ablating Case):
!
!     VASE MHL Trajectory (-11 deg EFPA)     ! Line 1 descriptor, FIAT V3
!     Shoulder Pt. Margined Data: HTC x 2.0; qrad x 1.5  ! Line 2 title
!
!     2.0                                    ! Factor applied here to HTC
!     1.5                                    ! Margin applied here to qrad
!     -8731292.                              ! Hformation (CO2, Venus, J/kg)
!     98.5                                   ! Time > 0 at which option 1 --> 3
!
!     free_g11.dat                           ! Trajectory freestream dataset
!     pw-shoulder-pt.eval.dat                ! Pressure, Pa from Traj_Fit
!     HTC-shoulder-pt.eval.dat               ! Heat Transfer Coef., kg/m^2.s
!     qrad-shoulder-pt.eval.dat              ! Radiative heat flux, W/cm^2
!
!     VASE-MHL-shoulder-pt.FIATv3.dat        ! Output in English units
!     VASE-MHL-shoulder-pt.Traj_FIAT.dat     ! Output in CGS units
!
!     2                                      ! Optional output thinning factor
!
!  Control File Notes (Ablating Case):
!
!     o  Trailing ! comments are ignored.
!     o  Blank lines improve readability but are optional.
!     o  Scale factors for heat transfer coefficient and qrad are provided
!        here, but specifying 1.0 may be the wise choice, as FIAT v3 can also
!        scale the data similarly and more than one such scaling may be of
!        interest.
!     o  N.B.:  If qconv data are in W/cm^2 instead of W/m^2, the scale factor
!        can be set to 1.e+4 to get unmargined W/cm^2 in the output files.
!     o  Traj_FIAT does not provide for such scaling.  It is intended for
!        stagnation point sizing only, and the two bounds provided by the max.
!        heat load trajectory (shallow entry) and the max. heat rate trajectory
!        (steep entry) are certainly not the same as what scaling some nominal
!        trajectory time history data would give.  Still, some users may
!        choose to make use of the scale factors provided here.
!     o  Recovery enthalpy is traditionally approximated by boundary layer edge
!        enthalpy.  Use of free-stream enthalpy has even been advocated as the
!        least dubious choice in separated flows.  But in such regions, the TPS
!        is most likely to be nonablating, so no enthalpy issue arises (see
!        the nonablating case below).
!     o  In the example, -8731292. is -8932880. (Hform) * 0.977433 (mass fr.);
!        then,        enthalpy = 0.5 V**2 + Hformation       is calculated here
!     o  [This Hform0 explanation is evidently about to be changed ...]
!     o  A freestream trajectory dataset from Traj is needed for V here, and
!        the time steps from such a trajectory should be what was used for the
!        evaluations of the other quantities.
!     o  Three curve-fit time histories from Traj_Fit (or equivalent) are
!        expected, in the units shown.
!     o  Two forms of the output file are written as shown below.
!     o  If a final thinning factor is present, the outputs will be thinned
!        accordingly; the default is 1 (no thinning).
!
!  Control File Example (Nonablating Case):
!
!     VASE MHL Trajectory (-11 deg EFPA)     ! Line 1 descriptor, FIAT V3
!     Aft Stag. Pt. Margined Data: qconv x 3.0; qrad x 1.7  ! Line 2 title
!
!     3.0                                    ! Factor applied here to HTC
!     1.7                                    ! Margin applied here to qrad
!     -8731292.                              ! Hformation (NOT USED)
!     -999.                                  ! Time <= 0 => nonablating case
!
!     free_g11.dat                           ! Trajectory freestream dataset
!     pw-aft-stag-pt.eval.dat                ! Pressure (Pa) time history
!     qconv-aft-stag-pt.eval.dat             ! Convective heat flux, W/m^2
!     qrad-aft-stag-pt.eval.dat              ! Radiative heat flux, W/cm^2
!
!     VASE-MHL-aft-stag-pt.FIATv3.dat        ! Output in English units
!     VASE-MHL-aft-stag-pt.Traj_FIAT.dat     ! Output in CGS units - see notes
!
!     2                                      ! Optional output thinning factor
!
!  Control File Notes (Nonablating Case):
!
!     o  Trailing ! comments, blank lines, and scale factors: as for ablating.
!     o  Freestream and pressure time histories: as for the ablating case.
!     o  Convective heat flux time history is entered in place of enthalpy.
!     o  N.B.:  If qconv data are in W/cm^2 instead of W/m^2, the scale factor
!        can be set to 1.e+4 to get unmargined W/cm^2 in the output files.
!     o  Two forms of the output file are written as shown below.
!     o  If a final thinning factor is present, the outputs will be thinned
!        accordingly; the default is 1 (no thinning).
!
!  Input Freestream Trajectory File Example:
!
!     # Time   Vel.  FS Prs.  FS Dens. FS Temp. Stg. Prs. Mach q Cnv. q Rad.
!     # sec.  km/sec pascals   kg/m^3   deg. K   pascals  no.  W/cm^2 W/cm^2
!     #------ ------ -------- -------- -------- --------- ---- ------ ------
!        0.00 11.938 2.80e-07 1.97e-12  258.68   2.80e-04 46.6    0.1    0.0
!        0.10 11.938 2.84e-07 2.00e-12  258.68   2.85e-04 46.6    0.1    0.0
!        0.20 11.938 2.88e-07 2.03e-12  258.67   2.89e-04 46.6    0.1    0.0
!         :     :     :        :           :      :         :      :      :
!
!     o  Any number of top non-numeric header lines
!     o  No embedded header lines (see FILTER_LINES utility)
!     o  Note the indicated units as output by Traj
!
!  Input Curve Fit Evaluation Data:
!
!     # t,s  f4(rho(t), V(t))
!       0.0  4.51205504E-10
!       0.1  4.56654810E-10
!       0.2  4.63894489E-10
!        :    :
!
!     o  As output by Traj_Fit
!     o  Same time discretization as the trajectory file is assumed
!
!  Output FIAT v3 Environment File Example, Ablating Case:
!  (squeezed to fit, tab-delimited)
!
!     1       VASE MHL Trajectory (-11 deg EFPA)
!     Shoulder Pt. Margined Data: HTC x 2.0; qrad x 1.5  ! Line 2 title
!     -1   1.00   1.00   1.00   0.00
!     Opt Time  Enthalpy    qrad        HTC         Pressure   Blow
!         sec   BTU/lb      BTU/ft2.s   lb/ft2.s    atm
!     1   0.0   2.6881E+04  4.6283E-08  3.0034E-15  3.2894E-12  0.5
!     1   0.1   2.6881E+04  4.6908E-08  3.0630E-15  3.3413E-12  0.5
!     1   0.2   2.6881E+04  4.7532E-08  3.1229E-15  3.3932E-12  0.5
!     :    :     :           :           :           :           :
!
!     o  Some minor editing may be needed before running FIAT v3.
!
!  Output Traj_FIAT Environment File Example, Ablating Case:
!
!     t,s  Enthalpy,J/kg  qrad,W/cm^2  HTC,kg/cm^2.s  Pressure,Pa  Blow
!     0.0  6.2527E+07     4.8589E-47   6.0657E-18     6.9625E-11   0.5
!     0.1  6.2527E+07     5.2445E-47   6.1794E-18     7.11593-11   0.5
!     0.2  6.2527E+07     5.6542E-47   6.2934E-18     7.27038-11   0.5
!      :    :              :            :              :            :
!
!  Output FIAT v3 Environment File Example, Nonablating Case:
!
!     1       VASE MHL Trajectory (-11 deg EFPA)
!     Aft Stag. Pt. Margined Data: HTC x 3.0; qrad x 1.7  ! Line 2 title
!     -1   1.00   1.00   1.00   0.00
!     Opt Time  qconv       qrad       Unused CT   Pressure     TF
!         sec   BTU/ft2.s   BTU/ft2.s              atm
!     0   0.0   2.6881E-06  4.6283E-09       0.0    3.2894E-13  0.
!     0   0.1   2.6881E-06  4.6908E-09       0.0    3.3413E-13  0.
!     0   0.2   2.6881E-06  4.7532E-09       0.0    3.3932E-13  0.
!     :    :     :           :                :      :           :
!
!     o  Some minor editing may be needed before running FIAT v3.
!
!  Output Traj_FIAT Environment File Example, Nonablating Case:
!
!     t,s  qconv,W/cm^2  qrad,W/cm^2  Unused CT  Pressure,Pa  TF
!     0.0  6.2527E-05     4.8589E-48        0.0   6.9625E-12  0.
!     0.1  6.2527E-05     5.2445E-48        0.0   7.11593-12  0.
!     0.2  6.2527E-05     5.6542E-48        0.0   7.27038-12  0.
!      :    :              :                 :     :           :
!
!  History:
!
!     09/25/2014  D.A.Saunders  Initial design.
!     10/02/2014    "     "     Completed initial implementation.
!     12/12/2014    "     "     Fixed one typo in the header.
!     12/16/2014    "     "     Use a new Mach 2 time input to switch the option
!                               output from 1 (ablating) to 3 (non-ablating with
!                               radiation).
!     01/14/2015    "     "     Distinguished between ablating and nonablating
!                               cases.  Why didn't we think of that long ago?
!                               Thinning the output tables is also an option.
!     01/20/2015    "     "     FIAT expects 7 input columns, not 8.  Where
!                               did the "laminar" column come from?  It has
!                               been eliminated.
!     06/26/2015    "     "     A CORSAIR trajectory had 7-digit times.  Since
!                               input data came from FLOW_INTERP_2D, not BLAYER,
!                               the assumption that input qconv is in W/m^2 can
!                               be worked around by entering 1.e+4 as the scale
!                               factor normally used for applying a margin.
!                               This would give unmargined qconv in W/cm^2 and
!                               BTU/ft2.s.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use table_type_module  ! Both are in table_io.f90
   use table_io

   implicit none

!  Constants:

   integer, parameter :: luntraj           = 1   ! Input trajectory time history
   integer, parameter :: luneval           = 2   ! Evaluated curve fits
   integer, parameter :: lunout            = 3   ! Output files for FIAT
   integer, parameter :: maxbuf            = 132 ! Limit on control file inputs
   integer, parameter :: traj_col_t        = 1   ! Unlikely to change
   integer, parameter :: traj_col_V        = 2   !    "     "     "
   real,    parameter :: J_kg_to_BTU_lb    = 0.000238846
   real,    parameter :: kg_cm2_to_lb_ft2  = 2048.16144
   real,    parameter :: Pa_to_atm         = 9.86923267e-6
   real,    parameter :: W_cm2_to_BTU_ft2  = 0.88055
   real,    parameter :: half              = 0.5
   real,    parameter :: zero              = 0.
   character (1), parameter :: comment     = '!' ! Control file comments
   character (maxbuf) :: buffer, descriptor, title

!  Variables:

   integer :: i, ios, ithin, len_descriptor, len_title, nrows
   real    :: factor_HTC, factor_qrad, Hformation, option, toption3
   logical :: ablating

   type (table_type) :: freestream_table, pressure_fit_table, &
                        HTC_fit_table, qrad_fit_table, FIAT_table, &
                        Traj_FIAT_table

!  Execution:

   call read_controls ()  ! This also reads the four input datasets as tables
   if (ios /= 0) go to 99

!  Option to scale HTC/qconv and qrad:

    HTC_fit_table%values(2,:) =  HTC_fit_table%values(2,:) * factor_HTC ! qconv?
   qrad_fit_table%values(2,:) = qrad_fit_table%values(2,:) * factor_qrad

!  Set up the Traj_FIAT table, ablating or nonablating case:

   Traj_FIAT_table%nheader = 1
   Traj_FIAT_table%nrows   = nrows
   Traj_FIAT_table%ncols   = 6

   allocate (Traj_FIAT_table%header(Traj_FIAT_table%nheader), &
             Traj_FIAT_table%values(Traj_FIAT_table%ncols,nrows))

   if (ablating) then
      Traj_FIAT_table%header(1) = &
  '#   t,s   Enthalpy,J/kg     qrad,W/cm^2   HTC,kg/cm^2.s     Pressure,Pa  Blow'
      Traj_FIAT_table%values(1,:)=freestream_table%values(traj_col_t,:) ! t, s
      Traj_FIAT_table%values(2,:)=freestream_table%values(traj_col_V,:)**2 * &
                                  0.5e+6 + Hformation           ! Enthlp.,J/kg
      Traj_FIAT_table%values(3,:)=     qrad_fit_table%values(2,:) ! qrad,W/cm^2
      Traj_FIAT_table%values(4,:)=1.e-4*HTC_fit_table%values(2,:) ! CH,kg/cm^2.s
      Traj_FIAT_table%values(5,:)= pressure_fit_table%values(2,:) ! pw,Pa
      Traj_FIAT_table%values(6,:)= half                           ! Blowing par.
   else
      Traj_FIAT_table%header(1) = &
        '#   t,s  qconv,W/cm^2    qrad,W/cm^2  Unused CT  Pressure,Pa  Unused TF'
      Traj_FIAT_table%values(1,:)=freestream_table%values(traj_col_t,:) ! t, s
      Traj_FIAT_table%values(2,:)=1.e-4*HTC_fit_table%values(2,:) ! qconv,W/cm^2
      Traj_FIAT_table%values(3,:)=     qrad_fit_table%values(2,:) ! qrad,W/cm^2
      Traj_FIAT_table%values(4,:)= zero                           ! Unused CT
      Traj_FIAT_table%values(5,:)= pressure_fit_table%values(2,:) ! pw,Pa
      Traj_FIAT_table%values(6,:)= zero                           ! Unused TF
   end if

!  Set up the FIAT v3 table, ablating or nonablating case:

   FIAT_table%nheader = 5
   FIAT_table%nrows   = nrows
   FIAT_table%ncols   = 7  ! Store column 0 as real; write integer

   allocate (FIAT_table%header(FIAT_table%nheader), &
             FIAT_table%values(0:FIAT_table%ncols-1,nrows))

   FIAT_table%header(1)      = '1    ' // descriptor(1:len_descriptor)
   FIAT_table%header(2)      = title(1:len_title)
   FIAT_table%header(3)      = '-1   1.00   1.00   1.00   0.00'

   if (ablating) then
      FIAT_table%header(4)   = 'Opt  Time  Enthalpy        qrad            ' // &
                               'HTC             Pressure      Blow'
      FIAT_table%header(5)   = '      sec  BTU/lb          BTU/ft2.s       ' // &
                               'lb/ft2.s        atm'
      FIAT_table%values(0,:) = 1.0  ! Ablating option; switches to 3[.0] below
      FIAT_table%values(1,:) = freestream_table%values(traj_col_t,:)  ! t, s
      FIAT_table%values(2,:) =  Traj_FIAT_table%values(2,:) * J_kg_to_BTU_lb
      FIAT_table%values(3,:) =  Traj_FIAT_table%values(3,:) * W_cm2_to_BTU_ft2
      FIAT_table%values(4,:) =  Traj_FIAT_table%values(4,:) * kg_cm2_to_lb_ft2
      FIAT_table%values(5,:) =  Traj_FIAT_table%values(5,:) * Pa_to_atm
      FIAT_table%values(6,:) =  half
!!!   FIAT_table%values(7,:) =  1.0  ! "Laminar" flag, now omitted
      do i = 1, nrows
         if (freestream_table%values(traj_col_t,i) >= toption3) then
            FIAT_table%values(0,i) = 3.0  ! Suppress ablation calculations
         end if
      end do
   else
      FIAT_table%header(4)   = 'Opt  Time  qconv           qrad        ' // &
                               'Unused CT  Pressure  Unused TF'
      FIAT_table%header(5)   = '      sec  BTU/ft2.s       BTU/ft2.s   ' // &
                               '           atm'
      FIAT_table%values(0,:) = 3.0
      FIAT_table%values(1,:) = freestream_table%values(traj_col_t,:)  ! t, s
      FIAT_table%values(2,:) =  Traj_FIAT_table%values(2,:) * W_cm2_to_BTU_ft2
      FIAT_table%values(3,:) =  Traj_FIAT_table%values(3,:) * W_cm2_to_BTU_ft2
      FIAT_table%values(4,:) =  zero
      FIAT_table%values(5,:) =  Traj_FIAT_table%values(5,:) * Pa_to_atm
      FIAT_table%values(6,:) =  zero
   end if

!  Option to thin the tables before writing them:

   if (ithin /= 1) then
      call table_io_thin (     FIAT_table, ithin)
      call table_io_thin (Traj_FIAT_table, ithin)
   end if

!  Write the two forms of the output table, Traj_FIAT first:

   if (ablating) then
      table_io_format = '(f7.1, 4es16.8, f5.1)'
   else
      table_io_format = '(f7.1, 2es16.8, f5.1, es17.8, f5.0)'
   end if

   call table_io_write_real (lunout, Traj_FIAT_table, table_io_format_sep, ios)

!  Now the FIAT table in English units:

   if (ablating) then
      table_io_format = '(i1, f8.1, 4es16.8, f4.1)'
   else
      table_io_format = '(i1, f8.1, 2es16.8, f5.1, es16.8, f5.0)'
   end if

   call table_io_write_real (lunout,      FIAT_table, table_io_format_sep, ios)

99 continue

   contains

!     Internal procedures for program PREPARE_FIAT_DATA:

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_controls () ! See the program header for the input format.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer :: inumber, l, len_name, len_number
      real    :: rnumber

!     Execution:

      read (*, '(a)') descriptor           ! Will be output for FIAT v3
      len_descriptor = len_trim (descriptor)
      l = index (descriptor(1:len_descriptor), comment) - 1   ! Ignore any
      if (l > 0) len_descriptor = len_trim (descriptor(1:l))  ! trailing comment

      read (*, '(a)') title                ! Likewise
      len_title = len_trim (title)
      l = index (title(1:len_title), comment) - 1
      if (l > 0) len_title = len_trim (title(1:l))

      do  ! Allow blank lines for readability, looking for two factors and
         read (*, '(a)') buffer              ! Hformation and toption3
         l = len_trim (buffer)
         if (l > 0) exit
      end do
      len_number = l
      l = index (buffer(1:len_number), comment) - 1
      if (l > 0) len_number = len_trim (buffer(1:l))

      call decode_number (buffer(1:len_number), inumber, rnumber, ios)

      if (ios == 0) then  ! Neither an integer nor a real
         write (*, '(2a)') '*** Bad margin factor for Heat Transfer Coef.: ', &
                           buffer(1:len_number)
         ios = 1
         go to 99
      end if

      factor_HTC = rnumber

      read (*, '(a)') buffer
      len_number = len_trim (buffer)
      l = index (buffer(1:len_number), comment) - 1
      if (l > 0) len_number = len_trim (buffer(1:l))

      call decode_number (buffer(1:len_number), inumber, rnumber, ios)

      if (ios == 0) then  ! Neither an integer nor a real
         write (*, '(2a)') '*** Bad margin factor for qrad: ', &
                           buffer(1:len_number)
         ios = 1
         go to 99
      end if

      factor_qrad = rnumber

      read (*, '(a)') buffer
      len_number = len_trim (buffer)
      l = index (buffer(1:len_number), comment) - 1
      if (l > 0) len_number = len_trim (buffer(1:l))

      call decode_number (buffer(1:len_number), inumber, rnumber, ios)

      if (ios == 0) then  ! Neither an integer nor a real
         write (*, '(2a)') '*** Bad Hformation input: ', buffer(1:len_number)
         ios = 1
         go to 99
      end if

      Hformation = rnumber

      read (*, '(a)') buffer
      len_number = len_trim (buffer)
      l = index (buffer(1:len_number), comment) - 1
      if (l > 0) len_number = len_trim (buffer(1:l))

      call decode_number (buffer(1:len_number), inumber, rnumber, ios)

      if (ios == 0) then  ! Neither an integer nor a real
         write (*, '(2a)') '*** Bad toption3 input: ', buffer(1:len_number)
         ios = 1
         go to 99
      end if

      toption3 = rnumber
      ablating = toption3 > 0

      do  ! Allow blank lines for readability, looking for a file name
         read (*, '(a)') buffer
         len_name = len_trim (buffer)
         if (len_name > 0) exit
      end do

      l = index (buffer(1:len_name), comment) - 1
      if (l > 0) len_name = len_trim (buffer(1:l))

      freestream_table%filename = buffer(1:len_name)

      call table_io_read_real (luntraj, freestream_table, ios)
      if (ios /= 0) go to 99

      nrows = freestream_table%nrows

      read (*, '(a)') buffer
      len_name = len_trim (buffer)
      l = index (buffer(1:len_name), comment) - 1
      if (l > 0) len_name = len_trim (buffer(1:l))

      pressure_fit_table%filename = buffer(1:len_name)

      call table_io_read_real (luneval, pressure_fit_table, ios)
      if (ios /= 0) go to 99

      if (pressure_fit_table%nrows /= nrows) then
         write (*, '(a, 2i6)') &
            '*** # pw times does not match # trajectory times:', &
            pressure_fit_table%nrows, nrows
         ios = 1
         go to 99
      end if

      read (*, '(a)') buffer
      len_name = len_trim (buffer)
      l = index (buffer(1:len_name), comment) - 1
      if (l > 0) len_name = len_trim (buffer(1:l))

      HTC_fit_table%filename = buffer(1:len_name)  ! qconv if nonablating

      call table_io_read_real (luneval, HTC_fit_table, ios)
      if (ios /= 0) go to 99

      if (HTC_fit_table%nrows /= nrows) then
         write (*, '(a, 2i6)') &
            '*** # HTC times does not match # trajectory times:', &
            HTC_fit_table%nrows, nrows
         ios = 1
         go to 99
      end if

      read (*, '(a)') buffer
      len_name = len_trim (buffer)
      l = index (buffer(1:len_name), comment) - 1
      if (l > 0) len_name = len_trim (buffer(1:l))

      qrad_fit_table%filename = buffer(1:len_name)

      call table_io_read_real (luneval, qrad_fit_table, ios)
      if (ios /= 0) go to 99

      if (qrad_fit_table%nrows /= nrows) then
         write (*, '(a, 2i6)') &
            '*** # qrad times does not match # trajectory times:', &
            qrad_fit_table%nrows, nrows
         ios = 1
         go to 99
      end if

      do  ! Allow blank lines for readability, looking for an output file name
         read (*, '(a)') buffer
         len_name = len_trim (buffer)
         if (len_name > 0) exit
      end do

      l = index (buffer(1:len_name), comment) - 1
      if (l > 0) len_name = len_trim (buffer(1:l))

      FIAT_table%filename = buffer(1:len_name)

      read (*, '(a)') buffer
      len_name = len_trim (buffer)
      l = index (buffer(1:len_name), comment) - 1
      if (l > 0) len_name = len_trim (buffer(1:l))

      Traj_FIAT_table%filename = buffer(1:len_name)

      ithin = 1  ! No output table thinning is the default

      read (*, *, iostat=ios) ithin  ! Need not be there

 99   return

      end subroutine read_controls

   end program prepare_fiat_data
