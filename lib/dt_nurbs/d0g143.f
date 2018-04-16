      SUBROUTINE D0G143 (CMEM, IMEM, DMEM, IGEHDR, IGI, IER)

C     Erase an IGES type 143 (Bounded Surface) Entity and its assumed
C     "physically dependent" entities, namely the untrimmed surface and
C     the Boundary entities.
C
C   DYNAMIC MEMORY INPUT-OUTPUT
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IGEHDR  MEM header index to IGES Entity
C     IGI     If non-zero, MEM pointer to IGES Index in which IGES Entity is
C             located.  Also indicates caller is deleting Index reference to
C             this IGES Entity.  If zero, delete the Index reference.
C             (Actually, this task is delegated to D0GEER.)
C   OUTPUT:
C     IER     Returned error code.  Zero implies no errors.
C             +1  = Could not erase all subordinate entities.
C             -8  = Given IGES Index, IGI, is not the one the IGES Entity
C                   points to.  Since the internal pointer appears to be
C                   valid, it is assumed that the IGES Entity properly
C                   belongs to this other IGES Index and it is not deleted.
C
C   CALLS:
C     D0G141
C     D0GEER
C     D0PTR
C
C   HISTORY:
C     15Sep93 P. Kraushar   Created
C*****
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER IGEHDR, IGI, IER

      INTEGER ENTGE
      PARAMETER (ENTGE=238)
      INTEGER ENTGI
      PARAMETER (ENTGI=237)

      INTEGER IHDR, ITYP, IERX, IBGE, IXI, I, NBDRY, NDEX
      INTEGER IGDEX, IXGI, JGI

C****
      IER = 0
C     First get and verify the associated IGES Index
      IXI = IMEM(IGEHDR+2)
      JGI = IMEM(IXI)
      CALL D0PTR (CMEM, IMEM, DMEM, JGI, ENTGI, 0, IHDR, ITYP, IERX)
      IF (IERX .NE. 0) THEN
C        Given IGES Entity does not point to a good Index, check IGI
         JGI = IGI
         IF (JGI .NE. 0) THEN
            CALL D0PTR (CMEM, IMEM, DMEM, JGI, ENTGI, 0, IHDR, ITYP,
     +                  IERX)
            IF (IERX .NE. 0)  JGI = 0
         ENDIF
      ELSE IF (IGI .NE. JGI .AND. IGI .NE. 0) THEN
C        The Entity points to a good Index, but not the one given as IGI
         IER = -8
         GOTO 9000
      ENDIF
      IF (JGI .NE. 0) THEN
C        We have a good Index, locate the list of MEM pointers to Enitites
         IXGI = IMEM(IHDR+2) + 1
         NDEX = IMEM(IXGI)/2
C        Locate and erase the surface entity
         IGDEX = (IMEM(IXI+20) + 1)/2
         IF (IGDEX .GT. 0 .AND. IGDEX .LE. NDEX .AND.
     +       2*IGDEX .NE. IMEM(IXI+20)) THEN
C           The index into the Index looks okay, check the object
            IBGE = IMEM(IXGI+IGDEX)
            CALL D0PTR (CMEM, IMEM, DMEM, IBGE, ENTGE, 0, IHDR,
     +                  ITYP, IERX)
            IF (IERX .EQ. 0) THEN
C              The object is an IGES Entity
               ITYP = IMEM(IMEM(IHDR+2)+4)
               IF (ITYP .EQ. 128) THEN
C                 The IGES Entity is a B-spline Surface (type 128),
C                 disconnect and erase it.
                  IMEM(IXGI+IGDEX) = 0
                  CALL D0GEER (CMEM, IMEM, DMEM, IHDR, JGI, IERX)
                  IF (IERX .NE. 0)  IER = 1
C              ELSE IF (other types of known surface entities) THEN
C                 erase them
               ELSE
                  IER = 1
               ENDIF
            ENDIF
         ENDIF
C        Prepare to loop through the Boundary entities
         IXI = IXI + 21
         NBDRY = IMEM(IXI)
         DO 100 I=1,NBDRY
            IGDEX = (IMEM(IXI+I) + 1)/2
            IF (IGDEX .GT. 0 .AND. IGDEX .LE. NDEX .AND.
     +          2*IGDEX .NE. IMEM(IXI+I)) THEN
C              The index into the Index looks okay, check the object pointed to
               IBGE = IMEM(IXGI+IGDEX)
               CALL D0PTR (CMEM, IMEM, DMEM, IBGE, ENTGE, 0, IHDR,
     +                     ITYP, IERX)
               IF (IERX .EQ. 0) THEN
C                 The object is an IGES Entity
                  ITYP = IMEM(IMEM(IHDR+2)+4)
                  IF (ITYP .EQ. 141) THEN
C                    The IGES Entity is a Boundary (type 141),
C                    disconnect and erase it.
                     IMEM(IXGI+IGDEX) = 0
                     CALL D0G141 (CMEM, IMEM, DMEM, IHDR, JGI, IERX)
                     IF (IERX .NE. 0)  IER = 1
                  ELSE
                     IER = 1
                  ENDIF
               ENDIF
            ENDIF
  100    CONTINUE
      ELSE
         IER = 1
      ENDIF
C     Finally, erase the given entity
      CALL D0GEER (CMEM, IMEM, DMEM, IGEHDR, IGI, IERX)

C     Error Handling

 9000 CONTINUE
      RETURN
      END
