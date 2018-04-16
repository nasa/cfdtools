      SUBROUTINE D0AHDR (CMEM, IMEM, DMEM, IHDR, LENC, LENI, LEND,
     +    LHDR)
C
C   Adjust header blocks whose data spaces are above IHDR to reflect
C   elimination of LENC characters, LENI integers, and LEND doubles from
C   the data spaces of the entity with header at IHDR.  Also used by D2DEFX
C   with negative LENC, LENI, and/or LEND.  Does not adjust the data spaces
C   themselves, only the indexes in the header blocks.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IHDR    Index to header experiencing the actual change in its data
C             spaces.
C     LENC    Number of characters being removed from IHDR entity
C     LENI    Number of integers being removed from IHDR entity
C     LEND    Number of double precision values being removed from IHDR entity
C   OUTPUT:
C     LHDR    Index to header immmediately preceding IHDR in active entity
C             list, or zero, if IHDR is the head of the list.
C
C   CALLS:    none
C   HISTORY:
C     03Feb92 P Kraushar    Created.
C*****

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER IHDR, LENC, LENI, LEND, LHDR, I

      INTEGER KHDR
C****
      LHDR = 0
      KHDR = IMEM(17)
      DO 100 I = 1, 1+(IMEM(2)-IMEM(16))/8
        IF (KHDR .EQ. IHDR) GO TO 110
        IMEM(KHDR+1) = IMEM(KHDR+1) - LENC
        IMEM(KHDR+2) = IMEM(KHDR+2) - LENI
        IMEM(KHDR+3) = IMEM(KHDR+3) - LEND
        LHDR = KHDR
        KHDR = IMEM(KHDR+7)/256
  100 CONTINUE
C     The preceding loop should never terminate by completing.  The count
C     is one more than the number of header blocks currently allocated.
      STOP 'Dynamic memory fatally corrupted - D0AHDR'
  110 CONTINUE
C     Adjust master parameter block free space indexes
      IMEM(4) = IMEM(4) - LENC
      IMEM(5) = IMEM(5) - LENI
      IMEM(6) = IMEM(6) - LEND
C     Adjust lengths in IHDR block itself
      IMEM(IHDR+4) = IMEM(IHDR+4) - LENC
      IMEM(IHDR+5) = IMEM(IHDR+5) - LENI
      IMEM(IHDR+6) = IMEM(IHDR+6) - LEND

      RETURN
      END
