C+**********************************************************************
C
      SUBROUTINE CONSCL(ICOPT,FMIN,FMAX,NCONT,AINC)
C
C ACRONYM: CONtour SCaLe
C          ---     -- -
C
C PURPOSE:
C   For a given data range CONSCL computes a "nice" scaling of about
C   NCONT values or else returns NCONT for a given a contour level
C   increment.
C
C ARGUMENTS:
C     ARG   TYPE I/O/S  DIM     DESCRIPTION
C    ICOPT    I    I     -      ICOPT = 0 means use input AINC to set 
C                               contour levels;
C                               ICOPT < 0 means get about NCONT nicely
C                               scaled values between FMIN and FMAX.
C    FMIN,    R   I/O    -      User supplied min and max on input;
C     FMAX                      on output are overwritten with updated
C                               min and max that totally enclose the
C                               original function range.
C    NCONT    I   I/O    -      User supplied number of contour levels
C                               for ICOPT<0; actual number of levels
C                               on output.
C    AINC     R   I/O    -      Input contour level increment when 
C                               ICOPT=0 (not overwritten on output).
C                               When ICOPT<0, AINC is returned with
C                               computed contour level increment.
C
C EXTERNAL REFERENCES: None
C
C KNOWN SYSTEM DEPENDENCIES: None
C
C STANDARDS VIOLATIONS: None.
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE   INITIALS   DESCRIPTION 
C    Jul. 86   RCL      Adapted for use with PLLEVS 
C
C AUTHOR: Pieter G. Buning, NASA Ames Research Center
C
C-**********************************************************************
C
      PARAMETER (NNICE=4)
      DIMENSION RNICE(4)
C
      DATA RNICE/.1,.2,.25,.5/
C
C   Compute the contour increment.
C
      IF (ICOPT.LT.0) THEN
C
C   As a first approximation, get the difference, its characteristic and
C   mantissa.
C
         DIFF= ABS(FMAX-FMIN)/(NCONT+1)
         IF (DIFF.LE.0.) THEN
C
C   Check for all values being the same -- just set up one contour level.
C
            NCONT= 1
            AINC = 1.
            GOTO 60
         ENDIF
C
         CHAR= ALOG10(DIFF)+1.
C
C   Round CHAR down and get the mantissa.
C
         IF (CHAR.GE.0.) THEN
            ICHAR = CHAR
         ELSE
            ICHAR = CHAR-1.
         ENDIF
         RMANT = DIFF*10.**(-ICHAR)
C
C   What's the next largest "nice" mantissa?
C
         DO 10 I = 1,NNICE
            IF (RMANT.LE.RNICE(I)) GOTO 20
   10    CONTINUE
         I = NNICE
C
C   Got a guess.  Check the number of levels.
C
   20    CONTINUE

         AINC = SIGN(RNICE(I),FMAX-FMIN)*10.**ICHAR
         IMIN = FMIN/AINC
         IF (FMIN.LT.IMIN*AINC) IMIN = IMIN-1
         IMAX = FMAX/AINC
         IF (FMAX.GT.IMAX*AINC) IMAX = IMAX+1
         NNEED = IMAX+1-IMIN
C
C   Are we under?
C
         IF (NNEED.GT.NCONT) THEN
C
C   Nope.  Try the next nice number.
C
            IF (I.LT.NNICE) THEN
               I = I+1
            ELSE
               ICHAR = ICHAR+1
               I = 1
            ENDIF
            GOTO 20
         ENDIF
C
C   Now just set up FMIN and FMAX and update NCONT.
C
         NCONT= NNEED
         FMIN = IMIN*AINC
         FMAX = IMAX*AINC
C
C   Compute the number of contour levels based on the increment.
C
      ELSE IF (ICOPT.EQ.0) THEN
         IMIN = FMIN/AINC
         IF (FMIN.LT.IMIN*AINC) IMIN= IMIN-1
         IMAX = FMAX/AINC
         IF (FMAX.GT.IMAX*AINC) IMAX= IMAX+1
         NCONT= IMAX-IMIN+1
         FMIN = IMIN*AINC
         FMAX = IMAX*AINC
      ENDIF
C
   60 CONTINUE
      RETURN
      END
