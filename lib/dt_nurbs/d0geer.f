      SUBROUTINE D0GEER (CMEM, IMEM, DMEM, IGEHDR, IGI, IER)

C     Erase an IGES Entity.  Do not attempt to locate and erase IGES
C     "physically dependent" entities.
C
C   DYNAMIC MEMORY INPUT-OUTPUT
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IGEHDR  MEM header index to IGES Entity
C     IGI     If non-zero, MEM pointer to IGES Index in which IGES Entity is
C             located.  Also indicates caller is deleting Index reference to
C             this IGES Entity.  If zero, this subroutine should delete the
C             Index reference.
C   OUTPUT:
C     IER     Returned error code.  Zero implies no errors.
C             -8  = Given IGES Index, IGI, is not the one the IGES Entity
C                   points to.  Since the internal pointer appears to be
C                   valid, it is assumed that the IGES Entity properly
C                   belongs to this other IGES Index and it is not deleted.
C
C   CALLS:
C     D0PTR
C     D2ERAS
C     DTERPT
C     DTERR
C
C   HISTORY:
C     08Sep93 P. Kraushar   Created
C*****
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER IGEHDR, IGI, IER

      INTEGER ENTGI
      PARAMETER (ENTGI=237)

      INTEGER IHDR, ITYP, IERX, ICQ, IXI, IXC, I, NPAR
      INTEGER IGDEX, IXGI, JGI

C****
      IER = 0
C     Extract critical data
      IXC = IMEM(IGEHDR+1) + 23
      IXI = IMEM(IGEHDR+2)
      JGI  = IMEM(IXI)
      ICQ  = IMEM(IXI+1)
      NPAR = IMEM(IXI+2)
      IGDEX  = (IMEM(IXI+13) + 1)/2
C     First, disconnect it from the Index when requested (IGI.eq.0).  Ignore
C     bad data other than being given the wrong IGES Index when a valid one
C     exists.  The first IF actually covers two logical conditions in one test.
      IF (IGI .NE. JGI) THEN
C        Check validity of internal back link to IGES Index
         CALL D0PTR (CMEM, IMEM, DMEM, JGI, ENTGI, 0, IHDR, ITYP,
     +               IERX)
         IF (IERX .EQ. 0) THEN
            IF (IGI .NE. 0) THEN
               IER = -8
               GOTO 9000
            ENDIF
            IXGI = IMEM(IHDR+2) + 1
            IF (IGDEX .GT. 0 .AND. IGDEX .LE. IMEM(IXGI)) THEN
C              Don't reference if IGDEX is not in range
               IF (IMEM(IXGI+IGDEX) .EQ. IMEM(IGEHDR)) THEN
                  IMEM(IXGI+IGDEX) = 0
                  JGI = 0
               ENDIF
            ENDIF
            IF (JGI .NE. 0) THEN
C              Check all used locations in the Index
               DO 100 I=IXGI+1,IXGI+IMEM(IXGI)
                  IF (IMEM(I) .EQ. IMEM(IGEHDR))  IMEM(I) = 0
  100          CONTINUE
            ENDIF
         ENDIF
      ENDIF
C     First scan for any TYPCOD 'T' or 'R' parameters and delete the
C     associated Character String Sequence entity.  Ignore any errors.
      CALL DTERPT (0)
      IXI = IXI + 18
      DO 200 I=1,NPAR
         IF (CMEM(IXC+I:IXC+I) .EQ. 'T') THEN
            IF (IMEM(IXI+I) .NE. 0) THEN
               CALL D2ERAS (CMEM, IMEM, DMEM, IMEM(IXI+I), IERX)
            ENDIF
         ENDIF
  200 CONTINUE
      IF (CMEM(IXI+NPAR:IXI+NPAR) .EQ. 'R') THEN
         IF (IMEM(IXI+NPAR) .NE. 0) THEN
            CALL D2ERAS (CMEM, IMEM, DMEM, IMEM(IXI+NPAR), IERX)
         ENDIF
      ENDIF
C     Erase the directly associated Character String Sequence entity
      IF (ICQ .NE. 0)  CALL D2ERAS (CMEM, IMEM, DMEM, ICQ, IERX)
C     Finally, erase the IGES Entity
      CALL D2ERAS (CMEM, IMEM, DMEM, IMEM(IGEHDR), IERX)
C     Resume normal error processing
      CALL DTERPT (1)
      RETURN

C     Error Handling

 9000 CONTINUE
      RETURN
      END
