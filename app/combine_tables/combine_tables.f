********************************************************************************

      program combine_tables

*        Combine related sets of single-function CBAERO-type tables into one
*        set of multi-function tables in Traj_opt's distributed heating
*        "aerothermal.database" format.
*
*        CBAERO-type CL, CD & CM tables can also be combined into a single table
*        with this utility as part of preparing a Traj_opt aerodynamic database.
*
*     Control file format:               ! Standard input on the command line
*     --------------------
*
*        Control file for merging of multiple 1-function tables from CBAERO
*        3                               ! # independent vars. (all tables)
*        1   3   2                       ! Permutation of input vars. to output
*        Variable  Scale   Shift   Format! Name, factors, precision for outputs
*        Mach        1.      0.     F4.1 ! Independent output variable 1
*        Alpha       1.      0.     F6.1 ! ........................... 2
*        Q      101325.      0.     F9.1 ! ........................... 3
*        2                               ! # dependent vars. combined in output
*        Function  Scale   Shift   Format! Name, factors, precision for outputs
*        Temp        1.      0.     F9.2 ! Dependent variable 1
*        Qdot        1.      0.    F11.5 ! .................. 2
*        3                               ! # multi-function tables in the output
*        HL-20.aerothermal.database      ! Output table name
*        hl20.TemperaturePoint.1.tbl     ! Input table for function 1, point 1
*        hl20.HeatFluxPoint.1.tbl        ! ........................ 2, ..... 1
*        hl20.TemperaturePoint.2.tbl     ! ........................ 1, ..... 2
*        hl20.HeatFluxPoint.2.tbl        ! ........................ 2, ..... 2
*        hl20.TemperaturePoint.3.tbl     ! ........................ 1, ..... 3
*        hl20.HeatFluxPoint.3.tbl        ! ........................ 2, ..... 3
*
*     Input table format:       ! Same rectangular dimensions for all inputs
*     -------------------
*
*        Title
*        nx1  nx2  nx3          ! Numbers of independent vars., 3-space assumed
*        x1 range               ! Heading
*         0.3                   ! Values of x1 (e.g., Mach) in rectangular table
*         0.6
*          :
*        26.0
*        x2 range
*        0.005
*        0.050
*         :
*        1.000
*        x3 range
*         0.0
*         2.0
*         4.0
*          :
*        50.0
*        Function data          ! Opposite of preferred Fortran order
*          x1      x2    x3     f
*         0.3   0.005   0.0   217.00026
*         0.3   0.005   2.0   216.54283
*          :     :       :       :
*
*     Output table format:
*     --------------------
*
*        Title
*        # Mach numbers
*        # Alphas
*        # Qs
*        # function types (temperature, heat flux, ...)
*        # multi-function tables which are concatenated here
*        SURFACE POINT 1
*        Mach    Alpha    Q       Temp     [  Qdot ...]
*         0.3     0.0     506.6    217.00  [  0.02987 ]
*         0.6     0.0     506.6    217.00  [  0.02987 ]
*          :       :         :        :        :
*        26.0    50.0   47880.3   3070.48  [402.38443 ]
*        SURFACE POINT 2
*        Mach    Alpha    Q       Temp     [  Qdot ...]
*         0.3     0.0     526.6    222.00  [  0.03222 ]
*          :       :         :        :        :
*
*     History:
*     --------
*
*        08/05/02  DAS  Initial implementation.
*        03/11/05   "   Minor description changes to accommodate application
*                       to CL, CD & CM data in the ModelCenter environment.
*        09/22/10   "   Use standard input for the control file.
*
*     Author:
*     -------
*
*        David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
*
********************************************************************************

      implicit none

*     Constants:
*     ----------

      integer, parameter ::
     >   lunin  = 2,      ! Logical unit for input table[s]
     >   lunout = 3,      ! ... and for output table
     >   lunctl = 5,      ! Control file, not necessarily combine_tables.inp
     >   luncrt = 6       ! For diagnostics

*     Variables:
*     ----------

      integer
     >   i, i1, i2, i3, ios, j, k, l, m,
     >   ncoords, nindep(3), nf, n_mf_tables, nx

      real
     >   x1, x2, x3

      real, allocatable, dimension (:) ::
     >   scalef, scalex, ! For units conversions
     >   shiftf, shiftx,
     >   v1, v2, v3      ! Coordinates of all rectangular tables

      real, allocatable, dimension (:,:) ::
     >   v               ! Allows permuting input order of dependent vars.

      real, allocatable, dimension (:,:,:) ::
     >   fin             ! Dependent functions for one set of input tables

      real, allocatable, dimension (:,:,:,:) ::
     >   fout            ! Dependent functions for one output combined table

      character
     >   filename * 80,  ! For all table files
     >   formout  * 64   ! For output format string

      character, allocatable, dimension (:) ::
     >   formf * 8, formx * 8, namef * 8, namex * 8

*     Execution:
*     ----------

***   open (lunctl, file='combine_tables.inp', status='old', iostat=ios)

***   if (ios /= 0) then
***      write (luncrt, '(/, a)')
***  >      ' Cannot open control file "combine_tables.inp".'
***      go to 999
***   end if

      read (lunctl, *) ! Title
      read (lunctl, *) nx                         ! # independent variables

      allocate (namex(nx), scalex(nx), shiftx(nx), formx(nx))

      read (lunctl, *) i1, i2, i3                 ! permuting indices
      read (lunctl, *)

      formout = '('
      m = 2
      do i = 1, nx
         read (lunctl, *) namex(i), scalex(i), shiftx(i), formx(i)
         l = len_trim (formx(i))
         formout(m:m+l-1) = formx(i)(1:l)
         m = m+l
         formout(m:m+1) = ', ' 
         m = m + 2
      end do

      read (lunctl, *) nf
      read (lunctl, *)

      allocate (namef(nf), scalef(nf), shiftf(nf), formf(nf))

      do i = 1, nf
         read (lunctl, *) namef(i), scalef(i), shiftf(i), formf(i)
         l = len_trim (formf(i))
         formout(m:m+l-1) = formf(i)(1:l)
         m = m+l
         formout(m:m+1) = ', '
         m = m + 2
      end do
      formout(m:m) = ')'

      read (lunctl, *) n_mf_tables  ! # multi-f. tables in output (concatenated)
      read (lunctl, *) filename     ! Output file name

      open (lunout, file=filename, status='unknown')

      write (lunout, '(a)') filename(1:len_trim (filename))

      do l = 1, n_mf_tables  ! For each set of related input tables

         do m = 1, nf  ! For each dependent function per set of input tables

            read (lunctl, *) filename
            open (lunin, file=filename, status='old', iostat=ios)
            if (ios /= 0) then
               i = len_trim (filename)
               write (luncrt, '(/, 3a)')
     >            ' Unable to open input file "', filename(1:i), '".'
               go to 999
            end if

            read (lunin, *) ! Title
            read (lunin, *) nindep(1:nx) ! Should check for changes

            if (l == 1 .and. m == 1) then

*              Only now can we write the output file header:

               write (lunout, '(i3, 2x, 3a)')
     >         nindep(i1), '! # ', namex(1)(1:len_trim (namex(1))), 's',
     >         nindep(i2), '! # ', namex(2)(1:len_trim (namex(2))), 's',
     >         nindep(i3), '! # ', namex(3)(1:len_trim (namex(3))), 's',
     >         nf,         '! # ', 'function', 's',
     >         n_mf_tables,'! # ', 'multi-function', ' tables'

*              Set up the rectangular table coords. once:

               ncoords = 0
               do i = 1, nx
                  ncoords = max (ncoords, nindep(i))
               end do

               allocate (v(ncoords,nx))

               do i = 1, nx
                  read (lunin, *)
                  read (lunin, *) v(1:nindep(i),i)
               end do

               allocate (v1(nindep(i1)), v2(nindep(i2)),
     >                   v3(nindep(i3)))

               v1 = v(1:nindep(i1),i1) * scalex(1) + shiftx(1)
               v2 = v(1:nindep(i2),i2) * scalex(2) + shiftx(2)
               v3 = v(1:nindep(i3),i3) * scalex(3) + shiftx(3)

               allocate (fin(nindep(1),nindep(2),nindep(3)),
     >                   fout(nindep(i1),nindep(i2),nindep(i3),nf))

            else ! Skip the rectangular coordinates (check some day)

               do i = 1, nx
                  read (lunin, *)
                  read (lunin, *) v(1:nindep(i),i)
               end do

            end if

            read (lunin, *) ! C++ order is opposite that of Fortran
            read (lunin, *) !  x1  x2  x3  f
            read (lunin, *) (((x1, x2, x3, fin(i,j,k),
     >                         k = 1, nindep(3)),
     >                         j = 1, nindep(2)),
     >                         i = 1, nindep(1))

*           Multiple subscripting means itemizing the permutations:

            select case (i1)

               case (1)

                  if (i2 == 2) then ! i3 = 3
                     do k = 1, nindep(3)
                        do j = 1, nindep(2)
                           fout(1:nindep(1),j,k,m) =
     >                        fin(1:nindep(1),j,k) * scalef(m) +
     >                           shiftf(m)
                        end do
                     end do

                  else ! i2 = 3, i3 = 2

                     do k = 1, nindep(2)
                        do j = 1, nindep(3)
                           fout(1:nindep(1),j,k,m) =
     >                        fin(1:nindep(1),k,j) * scalef(m) +
     >                           shiftf(m)
                        end do
                     end do

                  end if

               case default

                  write (luncrt, '(/, a)')
     >               ' This permutation of variables is not handled.'
                  go to 999

            end select

            close (lunin)

         end do ! Next function m for table set l

*        Write combined table l to the output file:

         write (lunout, '(a, i3)') 'Multi-function table', l
         write (lunout, '(a, 10(2x, a))')
     >      (namex(i)(1:len_trim (namex(i))), i = 1, nx),
     >      (namef(i)(1:len_trim (namef(i))), i = 1, nf)
         write (lunout, formout)
     >      (((v1(i), v2(j), v3(k), fout(i,j,k,1:nf),
     >         i = 1, nindep(i1)),
     >         j = 1, nindep(i2)),
     >         k = 1, nindep(i3))

      end do ! Next set l of related input tables


  999 continue

      end program combine_tables
