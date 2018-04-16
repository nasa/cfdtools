C+-----------------------------------------------------------------------------+
C
      subroutine HSORTRR( list1 ,list2 ,Length )
C
C Title: Sort 2 Real lists, using the first list as the key for the sort
C
C Purpose:
C   Sort a list (list1) of real values, and rearrange a second list so it
C   matches the first.  This is useful for sorting (X,Y) coordinates.
C
C Method:
C   Use the HeapSort algorithm of J.W.J. Williams to sort list1, and change
C   list2 whenever list1 is changed.  Since HeapSort is a sort-in-place
C   algorithm, no scratch storage is needed, and the original order of the
C   arrays is lost.
C
C References:
C   Press, et al. "Numerical Recipes", chapter 8.
C   Try Knuth, "Searching and Sorting" for a full discussion.
C 
C Error Handling:
C   Length <=0 will probably cause a Fortran array dimensions error, which
C   cannot be caught here, so the application should make sure it doesn't
C   happen.
C
C Notes: <none>
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
C            Jul 91  "             Added check for Length=1
C            Dec 91  D.A.Saunders  Composite IF changed to two IFs after
C                                  a Macintosh problem was reported.
C+---------------------------------------------------------------------+

C     >declare arguments
      integer Length            !length of list1, list2
      real    list1(Length)     !1st list of values to sort, used as sort key
     &       ,list2(Length)     !2nd list of values, not used for key

C-----------------------------------------------------------------------

C     >declare local variables
      integer lu ,ld            !pointers to the left and right boundaries of
     &       ,ru ,rd            ! the heap during sift-up and sift-down passes
      real    temp1 ,temp2      !temporaries used to swap elements in the lists


C     >begin execution

C     >handle special case of 1 element list
      if( Length .EQ. 1 ) goto 999

C     >initialize
      lu = ( Length / 2 ) +1
      ru = Length

C     >loop until whole array is sorted, ie. both down and up sifts are complete
 10   continue
C         >if lu>1, still sifting up
          if( lu .GT. 1 )then
C             >>sift up phase: move the left boundary down, and set the temp
              lu = lu -1
              temp1 = list1(lu)
              temp2 = list2(lu)

          else
C             >>sift down: move the min in the heap to the low end of the
C               sorted array
              temp1 = list1(ru)
              temp2 = list2(ru)
              list1(ru) = list1(1)
              list2(ru) = list2(1)
              ru = ru -1
C             >>see if we're done
              if( ru .EQ. 1 )then
                  list1(1) = temp1
                  list2(1) = temp2
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
C****&            list1(rd) .LT. list1(rd+1) ) rd = rd +1
              if( rd .LT. ru )then
                  if( list1(rd) .LT. list1(rd+1) ) rd = rd +1
              endif

C             >>sift down until the right element is sorted
              if( list1(rd) .GE. temp1 )then
                  list1(ld) = list1(rd)
                  list2(ld) = list2(rd)
                  ld = rd
                  rd = rd + rd

              else
                  rd = ru +1    !indicates end of sift-down

              endif

C             >>endloop 20
              goto 20

          endif      !{rd .GT. ru}

C         >>save the temp values
          list1(ld) = temp1
          list2(ld) = temp2

C     >>endloop 10 - keep looping until the smallest value has sifted to the top
      goto 10


C     >terminate execution
 999  return

      end
