      SUBROUTINE D2LPUT (CMEM, IMEM, DMEM, LP)

C     Locate pointers in User-defined data types.  Provided to serve as
C     a default value for the SUBROU argument in D2WRIT.  It declares all
C     user-defined types to be "unknown" to it.  Thus, it is the correct
C     routine to use if the user has not defined any user data types.  It may
C     also be satisfactory if no user-defined types contain pointers to other
C     entities.  The cost is extraneous warning messages.  For a better template
C     for a user-developed version of this subroutine see D0LPLT.
C
C   DYNAMIC MEMORY INPUT-OUTPUT
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT-OUTPUT:
C     LP      Array of context data.
C             LP(1) = Header index for entity
C             LP(2) = MEM data type number of entity.
C             LP(3) = 0  Initial call for this entity.
C                   > 0  Beginning absolute index of previous pointer block.
C             LP(4) = Ending absolute index of previous pointer block.
C             LP(5) = Absolute index of end of entity's integer data space.
C             Upon return, LP(3) and LP(4) have been updated to locate the
C             beginning and end of the next pointer block in the entity's
C             integer data space, or to LP(5)+1 if no more pointers.
C
C   CALLS:    none
C
C   HISTORY:
C     20Sep93 P. Kraushar   Created
C*****
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER LP(14)

C****
      LP(3) = -1
      RETURN
      END
