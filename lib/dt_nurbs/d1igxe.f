      SUBROUTINE D1IGXE (CMEM, IMEM, DMEM, IGEHDR, MPAR, MSTR, MSTRCH, 
     +                   MREAL, NEED, IER)

C     PURPOSE:
C        IGES Extend Entity.  Extend an IGES entity. Shadow routine of
C        D2IGXE.
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IGEHDR   Header index of IGES entity ID
C        MPAR     Amount by which to increase number of Parameters
C        MSTR     Amount by which to increase number of string 
C                 parameters to allocate
C        MSTRCH   Amount by which to increase total number of 
C                 characters in strings
C        MREAL    Amount by which to increase number of double 
C                 precision parameters
C
C     OUTPUT:
C        IER      Error flag.  If negative, the entity is not defined.
C                 -7  = Character space is locked.
C                 -8  = Integer space is locked.
C                 -9  = Double Precision space is locked.
C                 -10 = Insufficient space available for character data.
C                 -11 = Insufficient space available for integer data.
C                 -12 = Insufficient space available for double
C                       precision data.
C                 -17 = Entity does not point to a valid Character Sequence
C                       Entity
C                 -999= Unexpected error in a called subroutine.
C
C     CALLS:
C        D0PTR    pointer check
C        D1DEFX   entity definition extend
C        D1FEEI   fetch entity element
C        D1STEI   store entity element
C        DTERR    error handler
C
C
C     HISTORY:
C        21Jun93  D. Parsons   Created.
C        06Aug93  D. Parsons   Add support of DMEM "Free" stack.
C
C     ------

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

C     Get the Reserved Entity Type ID Numbers

C =====================================================================
C
C     INCLUDE FILE DITYID.INC
C
C     Reserved Entity Type ID Numbers
C
C     HISTORY
C        12May92     D. Parsons     created
C
C                          Character Sequence Entity
      INTEGER     ENTCQ
      PARAMETER  (ENTCQ = 249)
C =====================================================================

      INTEGER     MPAR, MSTR, MSTRCH, MREAL, IER
      INTEGER     ICQ, IERX, IGEHDR, ICQHDR, ITYP
      INTEGER     MOREC, MOREI, MORED, NEED, NPAR, NSTR
      INTEGER     ISTART, INEXT, I, LEND
      LOGICAL     INITIZ

      IER = 0
      IERX = 0
      
C     Fetch and check pointer to Character Sequence
      CALL D1FEEI (CMEM, IMEM, DMEM, IGEHDR, 2, ICQ, IERX)
      IF (IERX .LT. 0) THEN 
         IER = -17
         GOTO 9000
      ENDIF
      
      CALL D0PTR (CMEM, IMEM, DMEM, ICQ, ENTCQ, 0, ICQHDR, ITYP, IERX)
      IF (IERX .LT. 0) THEN 
         IER = -17
         GOTO 9000
      ENDIF
      
C     Extend the data in IGE

      MOREC = MPAR
      MOREI = MPAR
      MORED = MREAL
      INITIZ = .TRUE.

      IF ((MOREC .GT. 0) .OR. (MOREI .GT. 0) .OR. (MORED .GT. 0)) THEN
         CALL D1DEFX (CMEM, IMEM, DMEM, IGEHDR, MOREC, MOREI, MORED,
     +      INITIZ, NEED, IERX)
         IF (IERX .NE. 0) THEN
            IF ((IERX .LE. -8) .AND. (IERX .GE. -13)) THEN
               IER = IERX + 1
            ELSE
               IER = -999
            ENDIF
            GOTO 9000
         ENDIF
         
C        Add new REAL data to "free" stack.

         IF (MORED .GT. 0) THEN
            CALL D1FEEI (CMEM, IMEM, DMEM, IGEHDR, 4, ISTART, IERX)
            
            LEND = ABS(IMEM(IGEHDR+6))
            INEXT = LEND-MORED+1
            
            CALL D1STEI (CMEM, IMEM, DMEM, INEXT, 4, IGEHDR, IERX)
            
            DO 100 I = 1, MORED-1
               CALL D1STED (CMEM, IMEM, DMEM, DBLE(INEXT+1), INEXT, 
     +                      IGEHDR, IERX)
               INEXT = INEXT+1
  100       CONTINUE
  
            CALL D1STED (CMEM, IMEM, DMEM, DBLE(ISTART), LEND, IGEHDR, 
     +                  IERX)
         ENDIF  
      ENDIF
            
C     Extend the data in ICQ

      MOREC = MSTRCH
      MOREI = MSTR
      MORED = 0
      INITIZ = .TRUE.

      IF ((MOREC .GT. 0) .OR. (MOREI .GT. 0)) THEN
         CALL D1DEFX (CMEM, IMEM, DMEM, ICQHDR, MOREC, MOREI, MORED,
     +      INITIZ, NEED, IERX)
         IF (IERX .NE. 0) THEN
            IF ((IERX .LE. -8) .AND. (IERX .GT. -13)) THEN
               IER = IERX + 1
            ELSE
               IER = -999
            ENDIF
            GOTO 9000
         ENDIF
      ENDIF


C     Store the updated NPAR into IGE
      IF (MPAR .GT. 0) THEN
         CALL D1FEEI (CMEM, IMEM, DMEM, IGEHDR, 3, NPAR, IERX)
         CALL D1STEI (CMEM, IMEM, DMEM, NPAR+MPAR, 3, IGEHDR, IERX)
      ENDIF

C     Store the updated NSTR into ICQ
      IF (MSTR .GT. 0) THEN
         CALL D1FEEI (CMEM, IMEM, DMEM, ICQHDR, 1, NSTR, IERX)
         CALL D1STEI (CMEM, IMEM, DMEM, NSTR+MSTR, 1, ICQHDR, IERX)
      ENDIF

 9000 CONTINUE
      RETURN
      END
