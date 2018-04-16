      SUBROUTINE D2FEAC (CMEM, IMEM, DMEM, IDE, JLO, JHI, CA, IER)

C     PURPOSE:
C        Fetch subarray of character data from the data space of
C        entity IDE
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IDE      Pointer to entity to fetch
C        JLO      Relative index into the character space
C        JHI      Relative index into the character space
C
C     OUTPUT:
C        CA       Character array to put data into
C        IER      Error flag.  If negative, the subarray is unchanged
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = Null pointer (IDE = 0)
C                 -3  = Garbage value in pointer (points outside valid
C                       range in IMEM)
C                 -4  = Either a deleted entity or a garbage value in
C                       pointer
C                 -5  = Pointer refers to a deleted entity
C                 -8  = JLO <= 0
C                 -9  = JHI > LENC
C                 -10 = JLO > JHI
C
C     CALLS:
C        D0PTR
C        DTERR
C        D1FEAC
C
C     HISTORY:
C        10Mar92  D. Parsons   Created.
C        29Jun93  D. Parsons   Extracted guts to D1FEAC.
C
C     ------

C     Long name alias:
C        ENTRY D2_FETCH_CMEM_SUBARRAY (CMEM, IMEM, DMEM, IDE, JLO,
C    +         JHI, CA, IER)
      CHARACTER         CMEM*(*), CA*(*)
      INTEGER           IMEM(*), IDE, JLO, JHI
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IHDR, ITYP, IER

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2FEAC'/

C     ------

      IER = 0

C     Call the general pointer-check utility routine

      CALL D0PTR (CMEM, IMEM, DMEM, IDE, 0, 0, IHDR, ITYP, IER)
      IF (IER .EQ. 0) THEN
         CALL D1FEAC (CMEM, IMEM, DMEM, IHDR, JLO, JHI, CA, IER) 
      ENDIF
      
      IF (IER .NE. 0) THEN
         CALL DTERR (1, SUBNAM, IER, 0)
         CA = ' '
      ENDIF

 9999 CONTINUE

      RETURN
      END
