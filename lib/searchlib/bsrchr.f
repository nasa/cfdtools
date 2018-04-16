C+-------------------------------------------------------------------------
C
      LOGICAL FUNCTION BSRCHR( VALUE
     &                        ,LIST
     &                        ,LENGTH
     &                        ,POSITION  )

C ACRONYM: Binary SeaRCH through a Real number list
C
C PURPOSE: Searches a real array (1D) looking for the given value.  If it 
C          finds it, the function is true and POSITION contains the index
C          of the value in the array.  If not found, the function is false
C          and POSITION contains the index where the value would go in the
C          list (i.e. the index after the greatest list element less than
C          the value).  POSITION is LENGTH+1 if the value is greater than
C          the end of the list. POSITION is 1 if less than the beginning.
C          If the list is not sorted in ascending order, the result will
C          be meaningless, but the routine will run without error.
C
C METHOD: Binary search.
C
C ARGUMENTS:  type  IOS  DESCRIPTION
C   BSRCHR   logical O   function value is true if found, else false
C   VALUE     real  I    number to search list for
C   LIST      real  I    values to search through. MUST be sorted ascending
C   LENGTH    int   I    number of elements to search in list
C   POSITION  int    O   returned index of value, (see PURPOSE)
C
C LOCAL VARIABLES:
C    name     type   description
C  LOWER      int    lower boundary of search partition
C  UPPER      int    upper boundary of search partition, returned in POSITION
C  EPSILON    real   maximum error between VALUE and LIST that will be 
C                      considered a match (Machine dependent)
C
C FILES USED: none
C
C EXTERNAL REFERENCES: none
C
C ERROR HANDLING:
C   Empty list: function returns .false. and POSITION returns 1
C
C NOTES: none
C
C DEVOLOPMENT HISTORY:
C   DATE  INITIALS  DESCRIPTION
C   9/85    dbs     Initial coding.
C   3/86    dbs     #mod1.1# - added EPSILON value
C
C AUTHOR: David Serafini, Informatics General, Palo Alto  CA.
C
C--------------------------------------------------------------------
C     >>declare arguments

      INTEGER  LENGTH ,POSITION
      REAL VALUE ,LIST(LENGTH)

C     >>declare locals

      INTEGER LOWER ,UPPER
C#    #mod1.1# -- machine dependent value roughly 10*machine accuracy
      REAL EPSILON
      PARAMETER  (EPSILON = 0.000001)

C ---------  begin execution  ------------------

C     >>initialize
      BSRCHR = .FALSE.
      POSITION = 1

C     >>check for empty list
      IF(  LENGTH .LE. 0  )  RETURN

C     >>initialize boundary of partition
      LOWER = 0
      UPPER = LENGTH + 1

C     >>check for degenerate cases, and if none, execute
      IF(  VALUE .GT. LIST(LENGTH)+EPSILON  )THEN
           POSITION = UPPER

      ELSEIF(  VALUE .LT. LIST(1)-EPSILON  )THEN
C          >>same as default 
C          { POSITION = 1 }

      ELSE
 10        IF(  LOWER+1 .LT. UPPER  )THEN
                POSITION = ( LOWER + UPPER ) /2
                IF (  VALUE .GT. LIST(POSITION)+EPSILON  )THEN
                      LOWER = POSITION
                ELSE
                      UPPER = POSITION
                ENDIF

                GOTO 10
           ENDIF

           POSITION = UPPER

           BSRCHR =      ( VALUE .GE. LIST(POSITION)-EPSILON ) 
     &             .AND. ( VALUE .LE. LIST(POSITION)+EPSILON )

      ENDIF

      RETURN
      END
