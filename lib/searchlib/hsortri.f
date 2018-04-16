C+-----------------------------------------------------------------------------+
C
      subroutine HSORTRI ( rlist, ilist, Length )
C
C Title: Sort Real list and corresponding integer list
C
C Purpose:
C   Sort a list (rlist) of real values, and rearrange a second integer list
C   so it matches the first.  If ilist is initialized to 1, 2, 3, etc. before
C   this module is called, consecutive elements of ilist will identify
C   monotonically increasing elements of rlist on return.  This is useful
C   for sorting one or more vectors with regard to another.
C
C Method:
C   Use the HeapSort algorithm of J.W.J. Williams to sort rlist, and change
C   ilist whenever rlist is changed.  Since HeapSort is a sort-in-place
C   algorithm, no scratch storage is needed, and the original order of the
C   arrays is lost.
C
C References:
C   Press, et al. "Numerical Recipes", chapter 8.
C   Try Knuth, "Searching and Sorting" for a full discussion.
C 
C Error Handling: 
C   Length <=0 will probably cause a Fortran array dimensions error, which 
C   cannot be caught here, so the calling application should make sure
C   this doesn't happen.
C
C Notes: This is a simple modification of the HSORTRR subroutine to define
C   a sort permutation index array to be applied to any number of associated
C   vectors that are to be sorted in terms of one vector.
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
C            Jul 88  P.J.Trosin    Built "ilist=integer" version
C            Jul 91  D.B.Serafini  Added check for Length=1
C            Dec 91  D.A.Saunders  Use on a Macintosh by Roy Hampton
C                                  encountered failure.  Changed a
C                                  composite IF test into two IFs.
C            Oct 99  D.L.Hermstad  Macintosh/Absoft version honors case
C                                  of variables.  Therefore, changed
C                                  argument length to Length.  (Prior
C                                  version should have failed to compile
C                                  "real    rlist(Length)"!)
C+---------------------------------------------------------------------+

C     >declare arguments
      integer Length            !length of rlist, ilist
      real    rlist(Length)     !list of real values to sort, used as sort key
      integer ilist(Length)     !list of (integer) values, not used for key

C-----------------------------------------------------------------------

C     >declare local variables
      integer lu ,ld            !pointers to the left and right boundaries of
     &       ,ru ,rd            ! the heap during sift-up and sift-down passes
      real    rtemp             !temporaries used to swap elements 
      integer itemp             !in the lists


C     >begin execution

C     >handle special case of 1 element list
      if( Length .EQ. 1 ) goto 999

C     >initialize
      lu = ( Length / 2 ) +1
      ru = Length

C     >loop until whole array is sorted, i.e. both down & up sifts are complete
 10   continue
C         >if lu>1, still sifting up
          if( lu .GT. 1 )then
C             >>sift up phase: move the left boundary down, and set the temp
              lu = lu -1
              rtemp = rlist(lu)
              itemp = ilist(lu)

          else
C             >>sift down: move the min in the heap to the low end of the
C               sorted array
              rtemp = rlist(ru)
              itemp = ilist(ru)
              rlist(ru) = rlist(1)
              ilist(ru) = ilist(1)
              ru = ru -1
C             >>see if we're done
              if( ru .EQ. 1 )then
                  rlist(1) = rtemp
                  ilist(1) = itemp
                  goto 999   !return
              endif

          endif

C         >>set up for sift-down
          ld = lu
          rd = lu + lu

C         >>loop while right end of sorted region is before end of heap
 20       if( rd .LE. ru )then

C             >>if not at the end of the heap, compare this element to the
C               next higher element, and expand the sorted region by 1 if its
C               less
C****         if( rd .LT. ru  .AND.
C****&            rlist(rd) .LT. rlist(rd+1) ) rd = rd +1
              if( rd .LT. ru )then
                  if( rlist(rd) .LT. rlist(rd+1) ) rd = rd +1
              endif

C             >>sift down until the right element is sorted
              if( rlist(rd) .GE. rtemp )then
                  rlist(ld) = rlist(rd)
                  ilist(ld) = ilist(rd)
                  ld = rd
                  rd = rd + rd

              else
                  rd = ru +1    !indicates end of sift-down

              endif

C             >>endloop 20
              goto 20

          endif      !{rd .GT. ru}

C         >>save the temp values
          rlist(ld) = rtemp
          ilist(ld) = itemp

C     >>endloop 10 - keep looping until the smallest value has sifted to the top
      goto 10


C     >terminate execution
 999  return

      end
