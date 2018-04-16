C+-----------------------------------------------------------------------------+
C
      subroutine HSORTIIND ( ilist, iperm, Length )
C
C Title: Sort integer list, non-destructive of list.
C
C Purpose:
C   Sort a list (ilist) of integer values, and rearrange a second integer list
C   so it matches the first.  If iperm is initialized to 1, 2, 3, etc. before
C   this module is called, consecutive elements of iperm will identify
C   monotonically increasing elements of ilist on return.
C   (i.e. ilist(iperm(i)) increases with i.)
C
C Method:
C   Use the HeapSort algorithm of J.W.J. Williams to sort ilist, and change
C   iperm whenever ilist would change.  Since this implementation of HeapSort
C   is a sort-by-reference algorithm, no scratch storage is needed, and the
C   original order of the ilist array is retained.
C
C References:
C   Press, et al. "Numerical Recipes", chapter 8.
C   Try Knuth, "Searching and Sorting" for a full discussion.
C
C Error Handling:
C   Length <=0 will probably cause Fortran array dimensions error - there's
C   nothing we can do about it here.  The calling routine should insure this
C   doesn't happen.
C
C Notes: This is a modification of the HSORTII subroutine to sort ilist by
C   reference so that ilist is not destroyed.  Original HSORTII note follows.
C
C        This is a simple modification of the HSORTRR subroutine to define
C   a sort permutation index array to be applied to any number of associated
C   vectors that are to be sorted.
C
C Standards violations: (any or all of the following)
C   >lower case in code
C   >long variable names
C   >underbar "_" used in variable names
C   >"!" used for end-of-line comments
C
C Coding conventions:
C   UPPERCASE names are external references (both code and data)
C   Capitalized names are constants (parameter, read-only arguments, etc.)
C   lowercase names are local variables, read-write arguments, Fortran keywords
C   ">"s prefix normal comments, ">" for high level comments, ">>" for details
C   { } is used for assertions (statements assumed to be true when they appear)
C   C#  comments for noteworthy code (modifications, machine dependent
C         code, Fortran extensions, debugging code, etc)
C   Labels below 999 are for normal code (999 is always for STOP or RETURN)
C   Labels from 1000 to 1999 are format statements
C   Labels from 2000 on      are for error handling
C   Syntax for code modifications:
C#mod0.0new# -- comment about modification
C     ...new code
C#mod0.0old#
C#    ...replaced code
C#end0.0#
C
C Environment:  FORTRAN 77 + extensions
C
C Author:  David Serafini, Sterling Software, Palo Alto, CA
C
C Modification history:
C   ver 0.0  Jul 87  D.B.Serafini  Initial design and coding
C            Jul 88  P.J.Trosin    Built "integer" version
C            Jul 91  D.B.Serafini  Added check for Length=1
C            Dec 91  D.A.Saunders  Composite IF changed to two IFs as
C                                  suggested by a Macintosh problem.
C            Oct 99  D.L.Hermstad  Changed length to Length in argument
C                                  list for case-sensitive Absoft compiler.
C            Nov 00  D.L.Hermstad  Changed all ilist(i) to ilist(iperm(i))
C                                  to achieve a non-destructive sort.  Changed
C                                  routine name from HSORTII to HSORTIIND.
C+---------------------------------------------------------------------+

C     >declare arguments
      integer Length           !length of ilist, iperm
      integer ilist(Length)    !list of integer values to sort, used as sort key
      integer iperm(Length)    !list of permutation values, not used for key

C-----------------------------------------------------------------------

C     >declare local variables
      integer lu ,ld           !pointers to the left and right boundaries of
     &       ,ru ,rd           ! the heap during sift-up and sift-down passes
      integer iltemp           !temporaries used to swap elements
      integer iptemp           ! in the lists


C     >begin execution

C     >handle special case of 1 element list
      if( Length .EQ. 1 )goto 999

C     >initialize
      lu = ( Length / 2 ) +1
      ru = Length

C     >loop until whole array is sorted, ie. both down and up sifts are complete
 10   continue
C         >if lu>1, still sifting up
          if( lu .GT. 1 )then
C             >>sift up phase: move the left boundary down, and set the temp
              lu = lu -1
              iltemp = ilist(iperm(lu))
              iptemp = iperm(lu)

          else
C             >>sift down: move the min in the heap to the low end of the
C               sorted array
              iltemp = ilist(iperm(ru))
              iptemp = iperm(ru)
              iperm(ru) = iperm(1)
              ru = ru -1
C             >>see if we're done
              if( ru .EQ. 1 )then
                  iperm(1) = iptemp
                  goto 999   !return
              endif

          endif

C         >>set up for sift-down
          ld = lu
          rd = lu + lu

C         >>loop while right end of sorted region is before end of heap
 20       if( rd .LE. ru )then

C             >>if not at the end of the heap, compare this element to the
C               next higher element, and expand the sorted region by 1 if it's
C               less
C****         if( rd .LT. ru  .AND.
C****&            ilist(iperm(rd)) .LT. ilist(iperm(rd+1)) ) rd = rd +1
              if( rd .LT. ru )then
                  if( ilist(iperm(rd)) .LT. ilist(iperm(rd+1)) )
     &               rd = rd +1
              endif

C             >>sift down until the right element is sorted
              if( ilist(iperm(rd)) .GE. iltemp )then
                  iperm(ld) = iperm(rd)
                  ld = rd
                  rd = rd + rd

              else
                  rd = ru +1    !indicates end of sift-down

              endif

C             >>endloop 20
              goto 20

          endif      !{rd .GT. ru}

C         >>save the temp values
          iperm(ld) = iptemp

C     >>endloop 10 - keep looping until the smallest value has sifted to the top
      goto 10


C     >terminate execution
 999  return

      end
