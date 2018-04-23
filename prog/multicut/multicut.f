C***********************************************************************
C
      PROGRAM MULTICUT
C
C     Description:
C
C     MULTICUT is an adaptation of MULTIPLOT.  Rather than producing
C     scalar plots of surface data from a pair of "XYZ" and "Q" files
C     containing one or more surface datasets in PLOT3D /mgrid format,
C     it treats just the "XYZ" file and saves results as sections in
C     columns which should be straightforward to adapt as the input
C     geometry for a typical flow solver.
C
C     The cuts are made at the indicated stations (constant X and/or Z).
C
C     Input data format (PLOT3D /mgrid /whole /[un]formatted):
C
C        NGRID
C        IDIM1 JDIM1 1 [IDIM2 JDIM2 1]
C        All Xs for surface # 1
C        All Ys  "     "      "
C        All Zs  "     "      "
C       [All Xs for surface # 2]
C         :   :   :   :   :   :
C         :   :   :   :   :   :
C         :   :   :   :   :   :
C
C     Control inputs:
C
C        multicut.inp      See a sample
C
C     Outputs:
C
C        xcuts.geom        Cuts at specified Xs
C        zcuts.geom          "        "      Zs
C
C     History:
C
C     08/12/97  D.Saunders  Initial adaptation of MULTIPLOT, which was
C                           design around SYN87-SB plot utilities for
C                           SYN87-MB purposes.  MULTICUT was prompted
C                           by the need to derive body sections from a
C                           PLOT3D /mgrid surface definition.
C     04/14/98      "       Switched from CONT2D to CUTSURF/CUTORDER,
C                           which (apart from overcoming possible skewed
C                           (u,v) problems) permits assembling the line
C                           string(s) across all of the blocks for a given
C                           cut before writing them, as opposed to the
C                           original approach of processing each block
C                           independently.
C     10/10/03      "       The call to CUTORDER had an obsolete extra
C                           argument.
C
C     Author:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C
C***********************************************************************

      IMPLICIT NONE

C     Constants:

      INTEGER
     >   LUNCRT, LUNCTL, LUNDAT, LUNCUT,
     >   MAXBLK, MAXCUT, MAXI, MAXJ, MAXSEG, NADDIM
      REAL
     >   DEFAULT, EPS, ONE, ZERO

      PARAMETER
     >  (EPS    = 2.E-7, ! Fraction of data range for CUTSURF
     >   ONE    = 1.,
     >   ZERO   = 0.,
     >   NADDIM = 20,    ! 1 + max. # line strings per cut over all blocks
     >   LUNCRT = 6,     ! Screen
     >   LUNCTL = 1,     ! Control inputs
     >   LUNDAT = 2,     ! XYZ input file in PLOT3D /mgrid format
     >   LUNCUT = 3,     ! Output file(s) of cuts at specified stations
     >   MAXBLK = 140,   ! Max. # blocks (surfaces) per set
     >   MAXCUT = 300,   ! Max. # cuts in either direction
     >   MAXSEG = 3000,  ! Max. # 2-point line segments handled per cut
     >   MAXI   = 128,   ! Max. sizes of any one surface
     >   MAXJ   = 128, 
     >   DEFAULT= 999.)  ! Inputs read as DEFAULT are reset to sensible values

C     Variables:

      INTEGER
     >   NAD(NADDIM), NI(MAXBLK),  NJ(MAXBLK),
     >   NDGCUT(2,MAXSEG),  ! Pointers used by CUTSURF & CUTORDER
     >   I, IL, J, JL, L, N, NBLK, NXCUTS, NZCUTS
      REAL
     >   X(MAXI,MAXJ,MAXBLK), Y(MAXI,MAXJ,MAXBLK), Z(MAXI,MAXJ,MAXBLK), 
     >   XCUTS(MAXCUT), ZCUTS(MAXCUT),
     >   XCRV(MAXSEG), YCRV(MAXSEG), ZCRV(MAXSEG), FCRV(MAXSEG),
     >   P(3), S(3), SEG(4,2*MAXSEG), DATAMAX, DATAMIN, DRANGE, EPSXYZ
      LOGICAL
     >   BLANK(MAXBLK), FORMATTED, DOXCUTS, DOZCUTS
      CHARACTER
     >   FILENAME*48, FORM*11, TITLE*80

C     Execution:

C     ********************** Read control inputs ***********************

C     Control inputs are read as they are needed, as opposed to reading
C     them all ahead of time.

      OPEN (UNIT=LUNCTL, FILE='multicut.inp', STATUS='OLD')

      READ (LUNCTL, '(A)') TITLE
      L = LEN_TRIM (TITLE)

      READ (LUNCTL, *) FORMATTED
      IF (FORMATTED) THEN
         FORM = 'FORMATTED'
      ELSE
         FORM = 'UNFORMATTED'
      END IF

      READ (LUNCTL, *) FILENAME  ! Multiblock XYZ file

      OPEN (UNIT=LUNDAT, FILE=FILENAME, STATUS='OLD', FORM=FORM)

      IF (FORMATTED) THEN
         READ (LUNDAT, *) NBLK
      ELSE
         READ (LUNDAT) NBLK
      END IF

      IF (NBLK .GT. MAXBLK) THEN
         WRITE (LUNCRT, '(A)') ' Too many blocks in XYZ file.'
         GO TO 999
      END IF

      IF (FORMATTED) THEN
         READ (LUNDAT, *) (NI(N), NJ(N), I, N = 1, NBLK)
      ELSE
         READ (LUNDAT) (NI(N), NJ(N), I, N = 1, NBLK)
      END IF

      WRITE (LUNCRT, '(/, A)') ' Surface Block Dimensions'

      DO N = 1, NBLK

         IL = NI(N)
         JL = NJ(N)
         WRITE (LUNCRT, '(I3, 2I5)') N, IL, JL

         IF (IL .GT. MAXI .OR. JL .GT. MAXJ) THEN
            WRITE (LUNCRT, '(A, I3)') ' Too many XYZs in block', N
            GO TO 999
         END IF

         IF (FORMATTED) THEN
            READ (LUNDAT, *)
     >         ((X(I,J,N), I = 1, IL), J = 1, JL),
     >         ((Y(I,J,N), I = 1, IL), J = 1, JL),
     >         ((Z(I,J,N), I = 1, IL), J = 1, JL)
         ELSE
            READ (LUNDAT)
     >         ((X(I,J,N), I = 1, IL), J = 1, JL),
     >         ((Y(I,J,N), I = 1, IL), J = 1, JL),
     >         ((Z(I,J,N), I = 1, IL), J = 1, JL)
         END IF
      END DO
            
      CLOSE (LUNDAT)

      DO N = 1, NBLK                 ! Suppress cutting of any surfaces?
         READ (LUNCTL, *) I, BLANK(N)
      END DO

      READ (LUNCTL, *) DOXCUTS      ! Allows retention of suppressed cuts
      READ (LUNCTL, *) DOZCUTS
      READ (LUNCTL, *) NXCUTS

      IF (NXCUTS .GT. MAXCUT) THEN
         WRITE (LUNCRT, '(A)') ' Too many X cuts.'
         GO TO 999  ! More trouble than it's worth skipping over the excess
      END IF

      DO N = 1, NXCUTS
         READ (LUNCTL, *) I, XCUTS(N)
      END DO

      READ (LUNCTL, *) NZCUTS

      IF (NZCUTS .GT. MAXCUT) THEN
         WRITE (LUNCRT, '(A)') ' Too many Z cuts.'
         GO TO 999
      END IF

      DO N = 1, NZCUTS
         READ (LUNCTL, *) I, ZCUTS(N)
      END DO

      CLOSE (LUNCTL)

C     ********************* End of control inputs **********************


C     Determine the data range to provide CUTSURF with a tolerance:

      DATAMAX = X(1,1,1)  ! Necessary initialization for BOUNDS
      DATAMIN = DATAMAX
      DRANGE  = ZERO

      DO N = 1, NBLK
         IL = NI(N)
         JL = NJ(N)

         CALL BOUNDS (IL, JL, MAXI, X(1,1,N), DATAMIN, DATAMAX)
         DRANGE = MAX (DRANGE, DATAMAX - DATAMIN)
      END DO

      WRITE (LUNCRT, '(/, A, 1P, 2E16.6)')
     >   ' X data range: ', DATAMIN, DATAMAX

      DATAMAX = Y(1,1,1)  ! Necessary initialization for BOUNDS
      DATAMIN = DATAMAX

      DO N = 1, NBLK
         IL = NI(N)
         JL = NJ(N)

         CALL BOUNDS (IL, JL, MAXI, Y(1,1,N), DATAMIN, DATAMAX)
         DRANGE = MAX (DRANGE, DATAMAX - DATAMIN)
      END DO

      WRITE (LUNCRT, '(A, 1P, 2E16.6)')
     >   ' Y data range: ', DATAMIN, DATAMAX

      DATAMAX = Z(1,1,1)  ! Necessary initialization for BOUNDS
      DATAMIN = DATAMAX

      DO N = 1, NBLK
         IL = NI(N)
         JL = NJ(N)

         CALL BOUNDS (IL, JL, MAXI, Z(1,1,N), DATAMIN, DATAMAX)
         DRANGE = MAX (DRANGE, DATAMAX - DATAMIN)
      END DO

      WRITE (LUNCRT, '(A, 1P, 2E16.6)')
     >   ' Z data range: ', DATAMIN, DATAMAX

      WRITE (LUNCRT, '(/, A, 1P, E16.6)')
     >   ' Largest data range over all directions: ', DRANGE

      EPSXYZ = EPS * DRANGE


C     ********************* Cuts at specified Xs ***********************


      IF (NXCUTS .GT. 0 .AND. DOXCUTS) THEN

         OPEN (UNIT=LUNCUT, FILE='xcuts.geom', STATUS='UNKNOWN')

         WRITE (LUNCUT, '(A)') TITLE(1:L)

         P(1) = ZERO  ! Pretend it's a Z cut by switching X and Z data
         P(2) = ZERO
         P(3) = ONE   ! Reset to X for each cut - P is a point on the cut

         S(1) = ZERO  ! Vector from the point to define the cutting plane
         S(2) = ZERO
         S(3) = ONE

         DO I = 1, NXCUTS

            P(3) = XCUTS(I)

C           Cut all blocks at this X and write results to LUNCUT:

            CALL CUT ('X', P(3), MAXI, MAXJ, Z, Y, X)

         END DO

         CLOSE (LUNCUT)
         WRITE (LUNCRT, '(/, A)')
     >      ' Cross-stream cuts completed: xcuts.geom'

      END IF


C     ********************* Cuts at specified Zs ***********************


      IF (NZCUTS .GT. 0 .AND. DOZCUTS) THEN

         OPEN (UNIT=LUNCUT, FILE='zcuts.geom', STATUS='UNKNOWN')

         WRITE (LUNCUT, '(A)') TITLE(1:L)

         P(1) = ZERO
         P(2) = ZERO
         P(3) = ONE   ! Reset to Z for each cut

         S(1) = ZERO  ! Vector from the point to define the cutting plane
         S(2) = ZERO
         S(3) = ONE

         DO J = 1, NZCUTS

            P(3) = ZCUTS(J)

C           Cut all blocks at this Z and write results to LUNCUT:

            CALL CUT ('Z', P(3), MAXI, MAXJ, X, Y, Z)

         END DO

         CLOSE (LUNCUT)
         WRITE (LUNCRT, '(/, A)')
     >      '   Streamwise cuts completed: zcuts.geom'

      END IF

  999 CONTINUE

C     Internal procedure for MULTICUT:

      CONTAINS

C        ---------------------------------------------------------------
C
         SUBROUTINE CUT (CASE, ZSLICE, IDIM, JDIM, X, Y, Z)
C
C        CUT makes one cut of all surface patches at the specified "Z"
C        (which may be X at the higher level).  Results are written in
C        2-column form suited to typical flow solvers.
C
C        ---------------------------------------------------------------

C        Arguments:

         CHARACTER
     >      CASE * 1  ! 'X'[cut] or 'Z'[cut]
         INTEGER
     >      IDIM, JDIM
         REAL
     >      ZSLICE, X(IDIM,JDIM,*), Y(IDIM,JDIM,*), Z(IDIM,JDIM,*)

C        Local variables:

         INTEGER
     >      IER, L, L1, L2, LINE, NLINES, NNODE, NSEG, NSEGN
         REAL
     >      FN

C        Execution:

         NSEG  = 0
         NNODE = 0

         DO N = 1, NBLK

            IF (.NOT. BLANK(N)) THEN

C              Add the 2-pt. line segments from this block for this cut:

               IL = NI(N)
               JL = NJ(N)
               NSEGN = NSEG

C              CUTSURF expects (x,y,z,f), so pass x for f:

               CALL CUTSURF (X(1,1,N), Y(1,1,N), Z(1,1,N), X(1,1,N),
     >                       1, IL, 1, JL, IDIM, JDIM,
     >                       SEG, NNODE, 2*MAXSEG, NDGCUT,
     >                       NSEG, MAXSEG, P, S,
     >                       EPSXYZ, EPSXYZ, IER)

               IF (IER .NE. 0) THEN
                  WRITE (LUNCRT, 1001) 'CUTSURF', IER, N, CASE,
     >                                 ZSLICE, NSEG, 'Proceeding.'
               END IF

CCCCC          write (luncrt, '(a, 2i6)')
CCCCC>            ' nseg, nsegn: ', nseg, nsegn

               IF (NSEG .GT. NSEGN) THEN

C                 Convert the 2-pt. line segments to contiguous curves:

                  CALL CUTORDER (NAD, NADDIM, XCRV, YCRV, ZCRV, FCRV,
     >                           MAXSEG, SEG, NNODE, 2*MAXSEG, NDGCUT,
     >                           NSEG, MAXSEG, EPSXYZ, EPSXYZ, IER)

                  IF (IER .NE. 0) THEN
                     WRITE (LUNCRT, 1001) 'CUTORDER', IER, N, CASE,
     >                                ZSLICE, NSEG, 'Aborting this cut.'
                     EXIT  ! Skip remaining blocks
                  END IF

CCCCC             write (luncrt, '(a, i10)') ' NAD(1): ', NAD(1)
               END IF

            END IF

         END DO  ! Next block


         NLINES = NAD(1)

         IF (NLINES .GT. 0) THEN

            NAD(1) = 1  ! So the following works for line segment 1

            DO LINE = 1, NLINES
               L1 = NAD(LINE)
               L2 = NAD(LINE+1) - 1
               FN = REAL (L2 - L1 + 1)

               IF (CASE .EQ. 'X') THEN
                  WRITE (LUNCUT, 1002) '  FNFP          XF  FSEC'
                  WRITE (LUNCUT, 1003) FN, ZSLICE, ONE
                  WRITE (LUNCUT, 1002) '           Z           Y'
               ELSE  ! CASE = 'Z'
                  WRITE (LUNCUT, 1002) '  FNWP          ZW  FSEC'
                  WRITE (LUNCUT, 1003) FN, ZSLICE, ONE
                  WRITE (LUNCUT, 1002) '           X           Y'
               END IF

               WRITE (LUNCUT, '(2F12.6)') (XCRV(L), YCRV(L), L = L1, L2)

            END DO  ! Next line segment

         END IF

  999    CONTINUE

C        Formats:

 1001    FORMAT (/, ' *** ', A, ' error:', I2, '   Block:', I3, 3X, A1,
     >           ' cut:', F10.3, '   NSEG:', I5, /, 1X, A)
 1002    FORMAT (A)
 1003    FORMAT (F6.0, F12.6, F6.0)

         END SUBROUTINE CUT

      END PROGRAM MULTICUT
