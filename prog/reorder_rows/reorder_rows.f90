!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program reorder_rows
!
!  Purpose:
!
!     Reorder the (significant) rows in a tabulation so that the leading n
!     columns (presumed to be coordinates of a rectangular table for one or
!     more equal-sized "blocks") vary most rapidly in a different way.
!
!     The number of dimensions, n, may be 2 or 3 initially, but this could be
!     extended up to 7, which is the Fortran limit on number of subscripts.
!
!     This was prompted by a need to reorder an aerothermal database so that
!     the first 3 columns vary in the Fortran order (column 1 most rapidly,
!     column 2 next most rapidly, and so on).
!
!     File names, etc., are prompted for - there is no control file.
!
!  Strategy:
!
!     Header records tell how many dimensions, subtables, etc., but the number
!     of pure text lines heading each subtable is prompted for, along with the
!     desired output order.
!
!     All lines are read as character strings and stored in arrays for writing
!     in the specified order.
!
!     The column order is preserved; see "extract_columns" for reordering those.
!
!     Sample Input (aerothermal database):
!
!        2    # surface points  (1: stagnation, 2: shoulder)  (# subtables)
!        11   # alphas
!        25   # Machs
!        19   # dynamic pressures
!        // Point: 1 triangle: 5371     3.41388  0.0434033  -1.56086, id = 1010
!           Alpha       Mach       Qbar         Temp        Press        Ch
!           (deg)        (-)      (atm)      (deg.K)         (Pa)    (Kg/s-m2)
!             147        1.3      1e-06      121.856     0.221915  8.21526e-05
!             147        1.3      1e-05      155.482      2.21892  0.000264979
!             147        1.3     0.0001      205.681      22.1908  0.000864174
!              :          :         :           :           :           :
!             147        1.3      0.475      320.558       105399     0.141527 
!             147        2.0      1e-06      194.261     0.187776  8.98921e-05
!              :          :         :           :           :           :
!
!     The example shown has column 3 values varying most rapidly and column 1
!     values varying least rapidly.  Enter 321 at the appropriate prompt for
!     the existing order or permutation.
!
!     To change the row order so that the 2nd column values vary most rapidly
!     and the 3rd column values least rapidly, enter 213 as the desired output
!     permutation.
!
!  History:
!     
!        04/14/07  D.Saunders  Initial implementation for ndims = 2 and 3.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      lunin = 1, lunout = 2, lunkbd = 5, luncrt = 6, maxlen = 512

!  Variables:

   integer &
      iperm_in, iperm_out, itable, l, n, ndims, ntables, ntitle_lines

   integer, allocatable, dimension (:) :: &
      ncoords

   character :: &
      buffer * (maxlen), filein * 48, fileout * 48

   character, allocatable, dimension (:) :: &
      subtable * (maxlen)

!  Execution:

   write (luncrt, '(/, a)', advance='no') ' Input file name:  '
   read  (lunkbd, *) filein
   open  (lunin,  file=filein,  status='old')

   write (luncrt, '(a)', advance='no') ' Output file name: '
   read  (lunkbd, *) fileout
   open  (lunout, file=fileout, status='unknown')

   write (luncrt, '(a)', advance='no') ' # table dimensions (2 | 3): '
   read  (lunkbd, *) ndims

   allocate (ncoords(ndims))

!  Line 1 should contain the number of subtables in the file:

   read  (lunin, '(a)') buffer;  l = len_trim (buffer)
   read  (buffer(1:l), *) ntables
   write (lunout, '(a)') buffer(1:l)

   do n = 1, ndims
      read  (lunin, '(a)') buffer;  l = len_trim (buffer)
      read  (buffer(1:l), *) ncoords(n)
      write (lunout, '(a)') buffer(1:l)
   end do

   n = product (ncoords)

   allocate (subtable(n)) ! One character string for each line of a subtable

   write (luncrt, '(a)', advance='no') ' # column header lines (e.g., 1 or 2): '
   read  (lunkbd, *) ntitle_lines

   if (ndims > 2) then
      write (luncrt, '(/, 2a)') &
         ' A column permutation of 321 means column 3 is', &
         ' varying most rapidly and column 1 least rapidly.'
      write (luncrt, '(a)', advance='no') &
         ' Enter the column permutation existing in the input file: '
      read  (lunkbd, *) iperm_in
      write (luncrt, '(a)', advance='no') &
         ' Enter the desired column permutation (e.g., 123): '
      read  (lunkbd, *) iperm_out
   end if

!  For each subtable:

   do itable = 1, ntables

      do n = 1, ntitle_lines
         read  (lunin, '(a)') buffer
         l = len_trim (buffer)
         write (lunout, '(a)') buffer(1:l)
      end do

      select case (ndims)

         case (2)

            call reorder2  (ncoords(1), ncoords(2), subtable)  ! Local procedure

         case (3)

            call reorder3 (ncoords(1), ncoords(2), ncoords(3), subtable)

         case default

            write (luncrt, '(a)') ' Extend the utility for ndims > 3.'

      end select

   end do

   close (lunin)
   close (lunout)

   deallocate (ncoords, subtable)

!  Internal procedures for program reorder_rows:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine reorder2 (n1, n2, t)  ! Transpose a 2-dimensional table.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in) :: n1, n2
      character, intent (in), dimension (n1,n2) :: t * (*)

!     Local variables:

      integer :: i, j

!     Execution:

      write (lunout, '(a)') ((t(i,j)(1:len_trim(t(i,j))), j = 1, n2), i = 1, n1)

      end subroutine reorder2

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine reorder3 (n1, n2, n3, t)  ! Reorder the rows of a 3-D table.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in) :: n1, n2, n3
      character, intent (inout), dimension (n1,n2,n3) :: t * (*)

!     Local variables:

      integer :: i, j, k

!     Execution:

      select case (iperm_in)

         case (123)

            read (lunin, '(a)') (((t(i,j,k), &
                                i = 1, n1), j = 1, n2), k = 1, n3)

         case (132)

            read (lunin, '(a)') (((t(i,j,k), &
                                i = 1, n1), k = 1, n3), j = 1, n2)

         case (213)

            read (lunin, '(a)') (((t(i,j,k), &
                                j = 1, n2), i = 1, n1), k = 1, n3)

         case (231)

            read (lunin, '(a)') (((t(i,j,k), &
                                j = 1, n2), k = 1, n3), i = 1, n1)

         case (312)

            read (lunin, '(a)') (((t(i,j,k), &
                                k = 1, n3), i = 1, n1), j = 1, n2)

         case (321)

            read (lunin, '(a)') (((t(i,j,k), &
                                k = 1, n3), j = 1, n2), i = 1, n1)

         case default

            write (luncrt, '(/, a, i4)') ' Bad input permutation:', iperm_in

      end select

      select case (iperm_out)

         case (123)

            write (lunout, '(a)') (((t(i,j,k)(1:len_trim(t(i,j,k))), &
                                  i = 1, n1), j = 1, n2), k = 1, n3)

         case (132)

            write (lunout, '(a)') (((t(i,j,k)(1:len_trim(t(i,j,k))), &
                                  i = 1, n1), k = 1, n3), j = 1, n2)

         case (213)

            write (lunout, '(a)') (((t(i,j,k)(1:len_trim(t(i,j,k))), &
                                  j = 1, n2), i = 1, n1), k = 1, n3)

         case (231)

            write (lunout, '(a)') (((t(i,j,k)(1:len_trim(t(i,j,k))), &
                                  j = 1, n2), k = 1, n3), i = 1, n1)

         case (312)

            write (lunout, '(a)') (((t(i,j,k)(1:len_trim(t(i,j,k))), &
                                  k = 1, n3), i = 1, n1), j = 1, n2)

         case (321)

            write (lunout, '(a)') (((t(i,j,k)(1:len_trim(t(i,j,k))), &
                                  k = 1, n3), j = 1, n2), i = 1, n1)

         case default

            write (luncrt, '(/, a, i4)') ' Bad output permutation:', iperm_out

      end select

      end subroutine reorder3

   end program reorder_rows
