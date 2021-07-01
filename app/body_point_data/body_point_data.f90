!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program body_point_data
!
!  Description:
!
!     This utility is intended to ease the pain of gathering time histories for
!     body points from aerothermal solution surface datasets.  One body point is
!     treated at a time using a control file that should be straightforward to
!     prepare.  The intended input files are surface datasets from DPLR's
!     POSTFLOW utility, or outputs from BLAYER, and (optionally) outputs
!     gathered from a sequence of NEQAIR runs.
!
!  Control file format (read on standard input):
!
!     Body point data control file
!     2                               ! 2D | 3D
!     0.1234  1.2345                  ! Body point (x,y) or (x,y,z)
!     T       T       T               ! Laminar?   Turbulent?   Radiative?
!     times.dat                       ! List of times common to dataset groups
!     9                               ! Number of output variables
!       x y p T tau qw                ! Input variables to extract (laminar)
!     t x y p T tau qlam qturb qrad   ! utput var. names; t,x,y[,z]  mandatory
!     s m m Pa K Pa W/cm2 W/cm2 W/cm2 ! Output units
!     1 1 1 1 1 1 0.0001 0.0001 1     ! Scale factors to apply before output
!     body_point_2.dat                ! Output table file name
!     (f5.1, 2f10.6, es14.6, f8.2, es14.6, 3f8.3)  ! Output format
!     ! Laminar cases in time order (Tecplot surface datasets)
!     t149-surf.dat
!     t171-surf.dat
!     t192-surf.dat
!     t211-surf.dat
!     t222-surf.dat
!     t234-surf.dat
!     t245-surf.dat
!     t255-surf.dat
!     t263-surf.dat
!     t275-surf.dat
!     ! Turbulent cases in time order (Tecplot surface datasets; only qw used)
!     t149-surf.SST.dat
!     t171-surf.SST.dat
!     t192-surf.SST.dat
!     t211-surf.SST.dat
!     t222-surf.SST.dat
!     t234-surf.SST.dat
!     t245-surf.SST.dat
!     t255-surf.SST.dat
!     t263-surf.SST.dat
!     t275-surf.SST.dat
!     ! (x,y,qrad)  or (x,y,z,qrad) datasets
!     t149-xyqrad.dat
!     t171-xyqrad.dat
!     t192-xyqrad.dat
!     t211-xyqrad.dat
!     t222-xyqrad.dat
!     t234-xyqrad.dat
!     t245-xyqrad.dat
!     t255-xyqrad.dat
!     t263-xyqrad.dat
!     t275-xyqrad.dat
!     ! Output table file name
!     nose-tangency-point.dat
!
!  Assumptions:
!
!     o  Control file comments should begin with '!'.
!        This may be changed to another character in one place if necessary.
!
!     o  The input file of times that will become column 1 in the output table
!        should contain one time per line.  The number of times found should
!        match the length of each list of datasets to be searched for the
!        indicated body point.
!
!     o  Input dataset file names may be symbolic links or full paths.
! 
!     o  For the indicated body point coordinates, the nearest data point is
!        chosen from each dataset.  Exact coordinates are not required, and no
!        interpolation is attempted.  The coordinates are assumed to be common
!        to each set of like files.
!
!     o  Any of the laminar, turbulent, and radiative file lists may be
!        omitted in the control file (F instead of T), but the header in the
!        control file should not be omitted.  (LATER:) Allowing the laminar
!        input group to be suppressed is more trouble than it's worth, so
!        Laminar = T is assumed but retained for aesthetic reasons.
!
!     o  If laminar and turbulent are both specified, only qw is extracted from
!        the turbulent datasets, and it is appended to the laminar variables.
!
!     o  If radiative datasets are present, qrad is assumed to be in the last
!        input column (no need for column headers) and it will form the last
!        output column.
!
!  History:
!
!     05/13/2021  D.A.Saunders  Initial design.  This is overdue, and should be
!                               straightforward to implement using existing
!                               Tecplot I/O utilities and the relatively recent
!                               table_io module.
!     05/18/2021   "   "   "    Extracting the time column from the file names
!                               as in the initial design proved more trouble
!                               than it was worth, so read the times from a
!                               file.  One time per line simplifies work-space
!                               allocation.
!     05/20/2021   "   "   "    Initial coding is complete (on a laptop), but
!                               our only computer is down for maintenance.
!                               This has proved to be a lot messier than one
!                               might imagine, so testing is expected to take
!                               considerable effort.
!     05/24/2021   "   "   "    Debugged the compilation.  Testing has to wait.
!     05/26/2021   "   "   "    Testing proved arduous--mostly due to confusion
!                               associated with the various lists of variable
!                               names.  The 80-character filename limit in
!                               tecplot_io.f90 had to be raised (to 160 chars.).
!     05/26/2021   "   "   "    Fixed run-time formatting glitch.  The debug
!                               output is helpful as a dialogue, so leave it in
!                               there for now.  It would be nice to right-
!                               justify the column header tokens, based on the
!                               run-time format string.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! The tecplot_io version
   use grid_header_structure ! Likewise
   use string_justify        ! String justification utilities
   use table_type_module     ! Table data structure
   use table_io              ! Table I/O utilities
   use Tecplot_io_module     ! Tecplot I/O utilities

   implicit none

!  Constants:

   integer, parameter :: &
      lunin   = 1,       &   ! Input surface datasetrs
      lunout  = 2,       &   ! Output table
      lunctl  = 5,       &   ! Control file (standard input)
      luncrt  = 6,       &   ! Output diagnostics
      maxbuf  = 128,     &   ! Character limit on input control lines
      maxname = 32           ! Character limit on variable names | other tokens
                             ! as for tecplot_io.f90's varname(:)

   character (1), parameter :: &
      comment = '!'

!  Derived data types:

   type (grid_header) :: &
      header
   type (grid_type), dimension (:), pointer :: &
      surface_dataset
   type (table_type) :: &
      table

!  Variables:

   integer :: &
      i, icase, icol, idim, ios, iqturb, irow, lformat, line, lname, ntimes, &
      nvars, nvout
   integer, pointer :: &
      width(:)    ! Apparent column widths
   logical :: &
      laminar, radiative, turbulent
   real :: &
      bp_coords(3)
   character (maxbuf) :: &
      buffer, filename, format_string, table_name
   real, allocatable, dimension (:) :: &
      scale, times
   character (maxname), allocatable, dimension (:) :: &
      varname_in, varname_out, var_units

!  Explicit interface to allow array allocation in a called procedure:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   interface
      subroutine column_widths (format_string, ncol, width, ier)

      implicit none
      character*(*), intent (in)  :: format_string  ! Fortran format string,
                                                    ! including leading and
                                                    ! trailing parentheses
      integer,          intent (out) :: ncol        ! Apparent # columns
      integer, pointer, intent (out) :: width(:)    ! Apparent column widths,
                                                    ! allocated here
      integer,          intent (out) :: ier         ! 0 => no error detected;
                                                    ! 1 => flawed format?
      end subroutine column_widths
   end interface

!  Execution:

   call read_main_controls ()
   if (ios /= 0) go to 99

!  Set up the output table:

   table%nheader              = 2     ! Variable names and units
   table%nrows                = ntimes
   table%ncols                = nvars
   table%column_type(1:nvars) = 2     ! Reals
   table%filename             = table_name(1:lname)
   table_io_format            = format_string

   allocate (table%header(2))

   table%header(:) = ' '  ! Single-space the header names initially
   i = 0
   do icol = 1, nvars
      lname = len_trim (varname_out(icol))
      table%header(1)(i+1:i+lname) = varname_out(icol)(1:lname)
      i = i+lname+1
   end do
   i = 0
   do icol = 1, nvars
      lname = len_trim (var_units(icol))
      table%header(2)(i+1:i+lname) = var_units(icol)(1:lname)
      i = i+lname+1
   end do

   allocate (table%values(nvars,ntimes))

!  Process each list of datasets:

   do icase = 1, 3  ! Laminar, turbulent, radiative
      read (lunctl, '(a)', iostat=ios)  ! Skip header line in the control file
      if (ios /= 0) then
         if (icase == 3) exit  ! Handle missing radiative header
         write (luncrt, '(a)') 'Unexpected end-of-control-file?'
         go to 99
      end if
      if (icase == 1) then  ! Laminar
         if (.not. laminar) then
            write (luncrt, '(a)') 'Omitting laminar cases is not an option.'
            ios = 1
            go to 99
         end if
         nvout = nvars - 1 - idim  ! Not counting t or x,y[,z]
         if (turbulent) nvout = nvout - 1
         if (radiative) nvout = nvout - 1
      else if (icase == 2) then  ! Turbulent?
         if (.not. turbulent) cycle
         nvout = 1
      else  ! icase = 3; Radiative?
         if (.not. radiative) exit
      end if

      call process_list ()  ! Read ntimes files & transfer indicated vars.

      if (ios /= 0) then
         write (luncrt, '(a, i2)') 'Trouble processing dataset group', icase
         go to 99
      end if
   end do

!  Save the assembled table:

   do irow = 1, ntimes
      table%values(:,irow) = table%values(:,irow)*scale(:)
   end do

!  Right-justify the column headers:

   call column_widths (format_string, nvars, width, ios)

   if (ios /= 0) then
      write (luncrt, '(a)') 'Trouble determining column widths.'
   else
      call right_justify_names (nvars, varname_out, width, table%header(1))
      call right_justify_names (nvars, var_units,   width, table%header(2))
   end if

   call table_io_write_real (lunout, table, table_io_format_sep, ios)
   if (ios /= 0) then
      write (luncrt, '(2a)') 'Trouble writing table ', table%filename(1:lname),&
      'format_string: ', format_string
   end if

99 continue

!  Internal procedures for program body_point_data:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_main_controls ()  ! Except file lists
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, l, nfound

      line = 1
      read (lunctl, '(a)', iostat=ios)  ! Skip title
      if (ios /= 0) go to 90

      line = line + 1
      read (lunctl, *, iostat=ios) idim
      if (ios /= 0) go to 90

      header%ndim      = idim
      header%formatted = .true.

      line = line + 1
      read (lunctl, *, iostat=ios) bp_coords(1:idim)
      if (ios /= 0) go to 90

      line = line + 1
      read (lunctl, *, iostat=ios) laminar, turbulent, radiative
      if (ios /= 0) go to 90

!     Look for a list of times, one time per line:

      line = line + 1
      read (lunctl, '(a)', iostat=ios) filename
      if (ios /= 0) go to 90

      l = index (filename, comment) - 1  ! Suppress possible trailing comment
      if (l < 0) then
          l = len_trim (filename)
      else
          l = len_trim (filename(1:l))
      end if

      open (lunin, file=filename(1:l), status='old', iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(2a)') 'Unable to open list of times.  File name: ', &
            trim (filename)
         go to 99
      end if

!     Count the number of times:

      call count_records (lunin, .true., ntimes, ios)  ! File is rewound
      if (ios /= 0) go to 90

      allocate (times(ntimes))
      read (lunin, *, iostat=ios) times(:)
      if (ios /= 0) then
         write (luncrt, '(a)') 'Error reading the list of times.'
         go to 99
      end if
      close (lunin)

      write (luncrt, '(3a, /, (f7.2))') 'Times found in ', filename(1:l), ':', &
         times(:) 

      line = line + 1
      read (lunctl, *, iostat=ios) nvars  ! Large enough for all name lists
      if (ios /= 0) go to 90

      write (luncrt, '(a, i3)') 'Number of output variables:', nvars

      allocate (varname_in(nvars), varname_out(nvars), var_units(nvars), &
                scale(nvars))

!     The input variable names may be simple tokens or quoted strings.
!     Take the presence of any " to signify that all names are quoted, with
!     embedded blanks a possibility.  Since the rdlist* utilities read from a
!     file while token4 parses a buffer, we have a rare case where back-spacing
!     is a reasonable choice:

      line = line + 1
      read (lunctl, '(a)') buffer
      l = index (buffer, comment) - 1  ! Suppress possible trailing comment
      if (l < 0) l = len_trim (buffer)
      i = index (buffer(1:l), '"')
      nfound = nvars  ! Input with maximum number to read

      if (i == 0) then  ! Simple tokens are assumed
         backspace (lunctl)
         call rdlistc (lunctl, ' ', 1, nvars, nfound, varname_in)

!        Nasty one: rdlistc doesn't handle trailing conmments. Handle one here.

         l = -1
         do i = 1, nfound
            if (varname_in(i)(1:1) == comment) then
               l = i - 1  ! Avoid redefining loop limit nfound
               exit
            end if
         end do
         if (l > 0) nfound = l

      else  ! Handle quated variable names
         call token4 (buffer(1:l), ' ', '"', nfound, varname_in)
      end if

      write (luncrt, '(/, a)') 'Variable names expected as inputs:'
      do i = 1, nfound
         l = len_trim (varname_in(i))
         write (luncrt, '(i4, 3x, a)') i, varname_in(i)(1:l)
      end do

!     Repeat for the output variable names:

      line = line + 1
      read (lunctl, '(a)') buffer
      l = index (buffer, comment) - 1  ! Suppress possible trailing comment
      if (l < 0) l = len_trim (buffer)
      l = index (buffer(1:l), '"')
      nfound = nvars  ! Input with maximum number to read

      if (l == 0) then  ! Simple tokens are assumed
         backspace (lunctl)
         call rdlistc (lunctl, ' ', 1, nvars, nfound, varname_out)
      else  ! Handle quated variable names
         call token4 (buffer(1:l), ' ', '"', nfound, varname_out)
      end if

      nvars = nfound  ! nvars = actual number of output variables now

      write (luncrt, '(/, a)') 'Variable names specified for output:'
      do i = 1, nvars
         l = len_trim (varname_out(i))
         write (luncrt, '(i4, 3x, a)') i, varname_out(i)(1:l)
      end do

!     Repeat for the output unit names:

      line = line + 1
      read (lunctl, '(a)') buffer
      l = index (buffer, comment) - 1  ! Suppress possible trailing comment
      if (l < 0) l = len_trim (buffer)
      l = index (buffer(1:l), '"')
      nfound = nvars  ! Input with maximum number to read

      if (l == 0) then  ! Simple tokens are assumed
         backspace (lunctl)
         call rdlistc (lunctl, ' ', 1, nvars, nfound, var_units)
      else  ! Handle quated variable units
         call token4 (buffer(1:l), ' ', '"', nfound, var_units)
      end if

      write (luncrt, '(/, a)') 'Variable units specified for output:'
      do i = 1, nfound
         l = len_trim (var_units(i))
         write (luncrt, '(i4, 3x, a)') i, var_units(i)(1:l)
      end do

      if (nfound /= nvars) then
         write (luncrt, '(a, i4)') &
            'Number of output variable names:', nvars, &
            'Number of output variable units:', nfound
         ios = 1
         go to 99
      end if

!     Allow for (say) converting W/m^2 to W/cm^2:

      line = line + 1
      read (lunctl, *, iostat=ios) scale(1:nvars)
      if (ios /= 0) go to 90

      line = line + 1
      read (lunctl, '(a)', iostat=ios) table_name
      if (ios /= 0) go to 90

      lname = index (table_name, comment)  ! Suppress any trailing comment
      if (lformat == 0) then
         lname = len_trim (table_name)
      else
         lname = len_trim (table_name(1:lname - 1))
      end if

!     Run-time formatting option:

      line = line + 1
      read (lunctl, '(a)', iostat=ios) format_string
      if (ios /= 0) go to 90

      lformat = index (format_string, comment)  ! Suppress any trailing comment
      if (lformat == 0) then
         lformat = len_trim (format_string)
      else
         lformat = len_trim (format_string(1:lformat - 1))
      end if
      format_string(lformat+1:) = ' '
      go to 99

 90   write (luncrt, '(a, i4)') &
         '*** Trouble reading control file. Line number:', line

 99   return

      end subroutine read_main_controls

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine process_list ()

!     Read the indicated datasets and transfer the body point data variable(s)
!     according to icase:
!        icase = 1 means transfer t, x, y[, z] and nvout of the laminar vars.
!        icase = 2 means add turbulent qconv as output variable nvout+1;
!        icase = 3 means add radiative heat flux as column nvars.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, ibp, iq, irow, it, iv, izone, j, l, nbp, numq
      integer, save :: jcol
      real :: coords(3)
      real, allocatable :: xyqrad(:,:)  ! For the 2D case

      select case (icase)

         case (1)  ! t, x, y[, z] + nvout variables 

            do it = 1, ntimes  ! For each dataset in the laminar group
               read (lunctl, '(a)', iostat=ios) header%filename
               if (ios /= 0) then
                  write (luncrt, '(a, i2, a, i4)') &
                     'Error reading dataset name.  Group:', icase, '; time:', it
                  go to 99
               end if
               l = len_trim (header%filename)

               call Tecplot_read (lunin, header, surface_dataset, ios)
               if (ios /= 0) then
                  write (luncrt, '(2a)') &
                     'Trouble reading dataset ', header%filename(1:l)
                  go to 99
               end if

               irow = it
               table%values(1,irow) = times(it)
               table%values(2:1+idim,irow) = bp_coords(1:idim)
               jcol = 1 + idim
               numq = header%numq  ! # flow variables not counting coordinates

               do iv = 1, nvout

!                 What is the index of the input variable to be transferred?

                  call find_var (varname_in(idim+iv), numq, &
                                 header%varname(idim+1), iq)
                  if (iq == 0) then
                     write (luncrt, '(4a)') &
                        'Variable ', trim (varname_in(1+idim+iv)), &
                        ' not found in dataset ', header%filename(1:l)
                     ios = 1
                     go to 99
                  end if

                  if (index (varname_in(idim+iv), 'qw') > 0) iqturb = iq

!                 Find the variable value in some zone at the closest (x,y[,z]):

                  call find_body_pt (izone, i, j)

                  jcol = jcol + 1
                  table%values(jcol,irow) = surface_dataset(izone)%q(iq,i,j,1)

               end do  ! Next variable

               call deallocate_Tec ()

            end do  ! Next time

         case (2)  ! One variable (qturb)

            jcol = jcol + 1
            do it = 1, ntimes  ! For each dataset in the turbulent group
               read (lunctl, '(a)', iostat=ios) header%filename
               if (ios /= 0) then
                  write (luncrt, '(a, i2, a, i4)') &
                     'Error reading dataset name.  Group:', icase, '; time:', it
                  go to 99
               end if
               l = len_trim (header%filename)

               call Tecplot_read (lunin, header, surface_dataset, ios)
               if (ios /= 0) then
                  write (luncrt, '(2a)') &
                     'Trouble reading dataset ', header%filename(1:l)
                  go to 99
               end if

               call find_body_pt (izone, i, j)

               irow = it
               table%values(jcol,irow) = surface_dataset(izone)%q(iqturb,i,j,1)

               call deallocate_Tec ()

            end do  ! Next time

         case (3)  ! Radiation dataset group

            jcol = jcol + 1  ! Should match nvars

            do it = 1, ntimes  ! For each dataset in the radiation group
               irow = it
               read (lunctl, '(a)', iostat=ios) filename
               if (ios /= 0) then
                  write (luncrt, '(a, i2, a, i4)') &
                     'Error reading dataset name.  Group:', icase, '; time:', it
                  go to 99
               end if
               write (luncrt, '(2a)') 'Radiation file: ', trim (filename)

               l = len_trim (filename)
               open (lunin, file=filename(1:l), iostat=ios)
               if (ios /= 0) then
                  write (luncrt, '(3a, i4)') &
                     'Error opening radiation file ', filename(1:l), &
                     '; time:', it
                  go to 99
               end if

               if (idim == 2) then  ! Expect one or more x y qrad lines
                  call count_records (lunin, .true., nbp, ios)
                  allocate (xyqrad(3,nbp))
                  read (lunin, *, iostat=ios) xyqrad(:,:)
                  if (ios /= 0) then
                     write (luncrt, '(2a)') &
                       'Trouble reading x,y,qrad data from file ', filename(1:l)
                     go to 99
                  end if
                  call find_xy (nbp, xyqrad, ibp)
                  table%values(jcol,irow) = xyqrad(3,ibp)
                  deallocate (xyqrad)
               else  ! 3D qrad: x,y,z,qrad or polar_interp/surface_interp Tec
                  allocate (xyqrad(4,1))  ! Try reading line 1 as columns
                  read (lunin, *, iostat=ios) xyqrad(:,:)
                  deallocate (xyqrad)
                  rewind (lunin)
                  if (ios == 0) then
                     allocate (xyqrad(4,1))  ! xyzqrad really
                     read (lunin, *, iostat=ios) xyqrad(:,:)
                     if (ios /= 0) then
                        write (luncrt, '(2a)') &
                           'Trouble reading x,y,z,qrad as columns from ', &
                           filename(1:l)
                        go to 99
                     end if
                     table%values(jcol,irow) = xyqrad(4,1)
                  else  ! Must be 3D radiation in Tecplot surface form
                     header%filename = filename(1:l)
                     call Tecplot_read (lunin, header, surface_dataset, ios)
                     if (ios /= 0) then
                        write (luncrt, '(2a)') &
                           'Trouble reading dataset ', header%filename(1:l)
                        go to 99
                     end if

                     call find_body_pt (izone, i, j)

                     table%values(jcol,irow) = &
                        surface_dataset(izone)%q(iq,i,j,1)
                     call deallocate_Tec ()
                  end if
               end if
            end do

      end select

99    return

      end subroutine process_list

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine find_var (var, nlist, varlist, j)  ! Find it in a list of names
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      character (maxname), intent (in)  :: var            ! Target variable name
      integer,             intent (in)  :: nlist          ! # vars. to search
      character (maxname), intent (in)  :: varlist(nlist) ! Names to search
      integer,             intent (out) :: j              ! Index of var in list
                                                          ! after x,y[,z];
                                                          ! j = 0 if var is not
                                                          ! found in the list
      integer :: i, l

      l = len_trim (var);  j = 0

      write (luncrt, '(2a)') 'Input var: ', var(1:l)
      do i = 1, nlist
         write (luncrt, '(2a)') &
            ' List var: ', varlist(i)(1:len_trim (varlist(i)))
         if (var(1:l) == varlist(i)(1:l)) then
            j = i
            exit
         end if
      end do

      end subroutine find_var

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine find_body_pt (izbp, ibp, jbp)  ! Nearest grid pt, 2D or 3D Tec.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, intent (out) :: izbp, ibp, jbp

      integer :: izone, i, j
      real    :: diff(3), dsq, dsqmin

      dsqmin = 1.e+10
      do izone = 1, header%nblocks
         do j = 1, surface_dataset(izone)%nj
            do i = 1, surface_dataset(izone)%ni
               diff(1) = bp_coords(1) - surface_dataset(izone)%x(i,j,1)
               diff(2) = bp_coords(2) - surface_dataset(izone)%y(i,j,1)
               if (idim == 3) &
               diff(3) = bp_coords(3) - surface_dataset(izone)%z(i,j,1)
               dsq = dot_product (diff(1:idim), diff(1:idim))
               if (dsq < dsqmin) then
                  dsqmin = dsq;  izbp = izone;  ibp = i;  jbp = j
               end if
            end do
         end do
      end do

      write (luncrt, '(3(a,i4), a, es14.6)') &
         'izbp:', izbp, ' ibp:', ibp, ' jbp:', jbp, &
         ' dsqmin:', dsqmin

      end subroutine find_body_pt

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine find_xy (nbp, xyqrad, ibp)  ! Nearest grid pt, 2D radn. data
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, intent (in)  :: nbp            ! # body points to search
      real,    intent (in)  :: xyqrad(3,nbp)  ! One or more x,y,qrad
      integer, intent (out) :: ibp            ! Index of nearest grid point

      integer :: i
      real    :: diff(2), dsq, dsqmin

      dsqmin = 1.e+10
      do i = 1, nbp
         diff(1) = bp_coords(1) - xyqrad(1,i)
         diff(2) = bp_coords(2) - xyqrad(2,i)
         dsq = dot_product (diff(:), diff(:))
         if (dsq < dsqmin) then
            dsqmin = dsq;  ibp = i
         end if
      end do

      end subroutine find_xy

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine deallocate_Tec ()  ! Deallocate a surface dataset's variables
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: izone

      deallocate (header%varname)
      do izone = 1, header%nblocks
         deallocate (surface_dataset(izone)%x)
         deallocate (surface_dataset(izone)%y)
         if (idim == 3) deallocate (surface_dataset(izone)%z)
         deallocate (surface_dataset(izone)%q)
      end do
      deallocate (surface_dataset)

      end subroutine deallocate_Tec

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine right_justify_names (nnames, name, width, line)
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)  :: nnames        ! # column headers to be justified
      character (*), intent (in), dimension (nnames) :: name  ! Column headers
      integer, intent (in), dimension (nnames) :: width       ! Column widths
      character (*), intent (out) :: line    ! Character line of column headers,
                                             ! each right justified

!     Local variables:

      integer :: i1, i2, icol, l

!     Execution:

      l = len (line)
      i2 = 0

!     Rudimentary safeguards avoid going off the end of the line:

      do icol = 1, nnames
         i1 = min (l, i2 + 1)
         i2 = min (l, i2 + width(icol))
         call right_justify (width(icol), name(icol), line(i1:i2))
      end do

      end subroutine right_justify_names

   end program body_point_data
