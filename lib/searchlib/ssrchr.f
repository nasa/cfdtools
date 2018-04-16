C+----------------------------------------------------------------------
C
      logical function SSRCHR( Value ,List ,Length ,position )
C
C Acronym: Sequential SeaRCH for Real number
C
C Purpose:
C   Searches a real array (1D) looking for the given value.  If it 
C   finds it, the function is true and position contains the index
C   of the value in the array.  If not found, the function is false
C   and position contains the index of the last element +1.
C
C Method: 
C   Sequential search.
C
C Arguments:
C   Name       Dimens  Type  IOS  Description
C   ----       ------  ----  ---  -----------
C   SSRCHR            logical O   function value is true if found, else false
C   Value              real  I    number to search list for
C   List()     Length  real  I    values to search through
C   Length             int   I    number of elements to search in list
C   position           int    O   returned index of value, (see PURPOSE)
C
C Local constants:
C   Name     Type   Description(value)
C   ----     ----   ------------------
C   Epsilon  real   neighborhood of List() that value must be in to give a match
C                   (approx. 10*machine accuracy)
C
C Error Handling:
C   If length of list is zero, return false, with position = 1.
C
C Portability concerns:
C   Epsilon constant is machine dependent and set for VAX single precision
C
C Standards violations:
C   -  lower case in code
C   -  long variable names
C   -  underbar "_" used in variable names
C   -  "!" used for end-of-line comments
C
C Coding conventions:
C   UPPERCASE names are external references (both code and data)
C   Capitalized names are constants (parameter, read-only arguments, etc.)
C   lowercase names are local variables, read-write arguments, Fortran keywords
C   >>  notation is used for normal comments
C   { } is used for assertions (statements assumed to be true when they appear)
C   C#  comments for noteworthy code (modifications, machine dependent
C         code, Fortran extensions, debugging code, etc)
C   &   always indicates a continuation lines
C   Indentation is 4 spaces per level
C   Syntax for code modifications:
C#mod0.0new# -- comment about modification 
C     ...new code
C#mod0.0old# 
C#    ...replaced code
C#end0.0#
C
C Environment:  VAX/VMS FORTRAN 77 + extensions
C
C Author:  David Serafini, Sterling Software, Palo Alto, CA
C 
C Modification history:
C   ver 0.0   may86   D.B.Serafini  Initial design and coding
C 
C-----------------------------------------------------------------------

C     arguments
C     ---------
      real Value ,List(*)
      integer Length ,position

C     constants
C     ---------
      real Epsilon
      parameter( Epsilon = 1.E-6 )

C     execution
C     ---------

C     >>initialize
      SSRCHR = .False.
      position = 1

C     >>check for empty list
      if(  Length .LE. 0  )goto 999

C     >>loop through all elements in list, looking for a match
      SSRCHR = .True.
      do 10 position = 1, Length
           if( Value .GE. ( List(position)-Epsilon ) .AND. 
     &         Value .LE. ( List(position)+Epsilon ) )goto 999
 10   continue

C     >>if loop completes without a match, return false
C     { position eq Length+1 }
      SSRCHR = .False.


C     terminate execution
 999  return

      END
