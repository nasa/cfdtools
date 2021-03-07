!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program free_stream_properties
!
!  Description:
!
!     For flow field calculations along a given trajectory, certain desirable
!     quantities are not necessarily provided with the trajectory time history
!     data.  The two likely items prompting this utility are dynamic pressure
!     and Knudsen number.  Derivation of others could be added readily.
!
!        Dynamic pressure Q is easy: 0.5*rho*V**2
!
!     Knudsen number (particle mean free path / characteristic length L)
!     requires certain atmosphere-dependent constants that are problematic
!     at this level:
!
!         Kn = (M/Re) sqrt (gamma pi/2) where
!         Re = rho V L / mu             Reynolds number
!         mu = (C1 T**1.5)/(T + C2)     dynamic viscosity
!
!     For Earth entry, we use C1 = 1.458E-6 and C2 = 110.4 K.
!
!     For other atmospheres, only Q is added at present.
!
!     Flow regimes defined in terms of Knudsen number are:
!
!        continuum             Kn < 0.001
!        slip flow             0.001 - 0.1
!        transition flow       0.1 - 10
!        free molecular flow   > 10
!
!  Implementation Details:
!
!     A trajectory dataset is read as columns of data and these columns are
!     extended with extra columns for (initially) Q and (maybe) Kn.
!
!     SI units are assumed.
!
!     Use of the table_io package allows any number of header lines to be
!     transcribed to the output dataset and use of a run-time output format.
!
!  Sample Control File (Standard Input):
!
!     traj.dat                   Trajectory dataset file name
!     atmosphere                 Earth, Mars, Venus, Titan; Kn for Earth only
!     1.4                        Ratio of specific heats
!     1                          Column number for time t, s
!     3                              "     "     "  relative velocity V, m/s
!     4                              "     "     "  density, kg/m**3
!     5                              "     "     "  temperature T, K
!     9                              "     "     "  Mach number
!     0.23                       Characteristic length L (e.g., nose radius)
!     1.458E-6                   Sutherland's law constant C1, kg/(m.s.K**0.5)
!     110.4                      Sutherland's law constant C2, K
!     traj.out.dat               Output dataset file name
!     '(f10.4, ......, 2f10.6)'  Output table format
!
!  Output formatting note:
!     The enclosing quotes allow suppression of a trailing comment.
!     Alternatively, enter a format in the form (......) without quotes.
!
!  References:
!
!     (1) https://www.cfd-online.com/Wiki
!     (2) Dinesh Prabhu advises use of L = capsule nose radius, not diameter.
!
!  History:
!
!     04/19/2019  D.A.Saunders  Initial design (structure for any atmosphere).
!     04/22/2019    "      "    Initial coding (good for Earth entry only).
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use table_type_module     ! Both are in table_io.f90
   use table_io

   implicit none

!  Constants:

   integer, parameter ::  &
      lunin   = 1,        &  ! Input trajectory dataset
      lunout  = 2,        &  ! Output     "        "
      lunctl  = 5,        &  ! Control file (standard input)
      luncrt  = 6            ! Diagnostics  (standard output)
   real, parameter ::     &
      half    = 0.5,      &
      pi      = 3.14159,  &
      C1_air  = 1.458E-6, &  ! Sutherland's law C1, kg/(m.s.K**0.5)
      C2_air  = 110.4        ! Sutherland's law C2, K
      

!  Variables:

   integer :: &
      i, index_time, index_V, index_rho, index_T, index_Mach, ios, j, len, &
      nheader, nrows, ncols
   real :: &
      C1, C2, gamma, Kn, L, Mach, mu, Q, Re, rho, T, V
   character (8) :: &
      atmosphere

   type (table_type) :: &
      table_in, table_out

!  Execution:

   call read_controls ()     ! Internal procedure below
   if (ios /= 0) go to 99    ! Single STOP phhilosophy

   call computations ()
   if (ios /= 0) go to 99

   call save_results ()

99 continue

!  Internal procedures for program free_stream_properties:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_controls ()  ! ... and the input trajectory dataset
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      read (lunctl, *, iostat=ios) table_in%filename
      if (ios /= 0) then
         write (luncrt, '(a)') &
            ' Trouble reading input trajectory file name on standard input.'
         go to 99
      end if

      call table_io_read_real (lunin, table_in, ios)
      if (ios /= 0) then
         write (luncrt, '(a)') &
            ' Trouble reading input trajectory dataset. File name:', &
            trim (table_in%filename)
         go to 99
      end if

      read (lunctl, *, iostat=ios) atmosphere
      if (ios /= 0) then
         write (luncrt, '(a)') &
            ' Trouble reading atmosphere name (gas definition).'
         go to 99
      end if
      read (lunctl, *, iostat=ios) gamma
      if (ios /= 0) then
         write (luncrt, '(a)') ' Trouble reading gamma (e.g., 1.4).'
         go to 99
      end if
      read (lunctl, *, iostat=ios) index_time
      if (ios /= 0) then
         write (luncrt, '(a)') ' Trouble reading column # for time.'
         go to 99
      end if
      read (lunctl, *, iostat=ios) index_V
      if (ios /= 0) then
         write (luncrt, '(a)') &
            ' Trouble reading column # for relative velocity.'
         go to 99
      end if
      read (lunctl, *, iostat=ios) index_rho
      if (ios /= 0) then
         write (luncrt, '(a)') ' Trouble reading column # for density.'
         go to 99
      end if
      read (lunctl, *, iostat=ios) index_T
      if (ios /= 0) then
         write (luncrt, '(a)') ' Trouble reading column # for temperature.'
         go to 99
      end if
      read (lunctl, *, iostat=ios) index_Mach
      if (ios /= 0) then
         write (luncrt, '(a)') ' Trouble reading column # for Mach number.'
         go to 99
      end if
      read (lunctl, *, iostat=ios) L
      if (ios /= 0) then
         write (luncrt, '(a)') ' Trouble reading characteristic length L.'
         go to 99
      end if
      read (lunctl, *, iostat=ios) C1
      if (ios /= 0) then
         write (luncrt, '(a)') ' Trouble reading viscosity constant C1.'
         go to 99
      end if
      read (lunctl, *, iostat=ios) C2
      if (ios /= 0) then
         write (luncrt, '(a)') ' Trouble reading viscosity constant C2.'
         go to 99
      end if

      read (lunctl, *, iostat=ios) table_out%filename
      if (ios /= 0) then
         write (luncrt, '(a)') ' Trouble reading output file name.'
         go to 99
      end if

      read (lunctl, '(a)', iostat=ios) table_io_format  ! Whole input line
      len = index (table_io_format, ')''')   ! Suppress any trailing comment
      if (len > 0) then                      ! and any enclosing apostrophes
         table_io_format(1:1)    = ' '        
         table_io_format(len+1:) = ' '
      end if

      close (lunctl)

99    return

      end subroutine read_controls

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine computations ()
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      logical :: Knudsen

!     Set up the output table as the input table plus columns for Q and Kn:

      nheader           = table_in%nheader
      nrows             = table_in%nrows
      ncols             = table_in%ncols
      table_out%nheader = nheader
      table_out%nrows   = nrows
      table_out%ncols   = ncols + 2
      table_out%column_type(:)  = table_in%column_type(:)
      table_out%column_width(:) = table_in%column_width(:)

      allocate (table_out%header(nheader), &
                table_out%values(ncols+2,nrows), stat=ios)
      if (ios /= 0) then
         write (luncrt, '(a)') ' Trouble allocating output table.'
         write (luncrt, '(a, 2i6)') ' # rows & columns:', nrows, ncols+2
         go to 99
      end if

!     Derive the new columns. Knudsen number is gas-specific.

      call upcase (atmosphere);  ios = 0

      select case (atmosphere)

         case ('EARTH   ')

            Knudsen = .true.
            do i = 1, nrows
               rho  = table_in%values(index_rho, i)
               V    = table_in%values(index_V,   i)
               T    = table_in%values(index_T,   i)
               Mach = table_in%values(index_Mach,i)
               Q    = 0.5*rho*V**2
               mu   = (C1 * T**1.5)/(T + C2)  ! Dynamic visc. (Sutherland, air)
               Re   = rho*V*L/mu
               Kn   = (Mach/Re)*sqrt (gamma*pi*0.5)
               table_out%values(1:ncols,i) = table_in%values(1:ncols,i)
               table_out%values(ncols+1,i) = Q
               table_out%values(ncols+2,i) = Kn
            end do

         case ('MARS    ')
            Knudsen = .false.
         case ('VENUS   ')
            Knudsen = .false.
         case ('TITAN   ')
            Knudsen = .false.

         case default
            write (luncrt, '(2a)') ' Unknown atmosphere spec.: ', atmosphere
            ios = -1

      end select

      table_out%header(:) = table_in%header(:)
      len = len_trim (table_in%header(nheader))

      if (Knudsen) then
         table_out%header(nheader)(len+1:) = '    "Q,Pa"          "Kn"'
      else if (ios /= -1) then
         table_out%header(nheader)(len+1:) = '    "Q,Pa"    "Kn unknown"'
         write (luncrt, '(2a)') ' Suppressing Knudsen number output.', &
               ' Formulation is uncertain for this atmosphere.'
         do i = 1, nrows
            rho  = table_in%values(index_rho, i)
            V    = table_in%values(index_V,   i)
            Q    = 0.5*rho*V**2
            Kn   = 999.
            table_out%values(1:ncols,i) = table_in%values(1:ncols,i)
            table_out%values(ncols+1,i) = Q
            table_out%values(ncols+2,i) = 999.
         end do
      end if

99    return

      end subroutine computations

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine save_results ()
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     The table_io_format_sep argument (a constant within the module
!     table_type_module) means use run-time formatting via table_io_format:

      call table_io_write_real (lunout, table_out, table_io_format_sep, ios)
      if (ios /= 0) then
         write (luncrt, '(a)') ' Trouble writing output table. File name:', &
            trim (table_out%filename), &
            ' Active run-time format:', &
            trim (table_io_format)
         write (luncrt, '(a, /, 30es15.8)') 'Line 1:', table_out%values(:,1)
      end if

      end subroutine save_results

   end program free_stream_properties
