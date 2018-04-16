!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE TRIANGULATE_PATCHES (NPATCH, NIJ, XYZ, IP, ICOMPONENT,
     >                                REF_LEN, NPTS, NFACES, IV, ICOMP)
!
!     One-liner:  Triangulate structured surface patches.
!
!     TRIANGULATE_PATCHES converts a surface paneling (one or more structured
!     surface patches) to unstructured form.  The shorter diagonal of each
!     quadrilateral cell is used to define each pair of triangles.  Account
!     is taken of collapsed cell edges: all output triangles have nonzero
!     area.  However, no check is made for non-convex cells.
!
!     The patch points are expected to be in triplet form (possibly via
!     appropriate reading of an AEROSURF/PLOT3D form).  Thus there is no
!     need to copy (x,y,z)s.
!
!
!     10/12/00  DAS  Initial adaptation of Mark Rimlinger's Struct2Unstruct
!                    routine from UV_MAP.
!
!     Author:  David Saunders, ELORET/NASA Ames Research Center, Mtn. View, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (IN) ::
     >   NPATCH                ! Number of patches

      INTEGER, INTENT (IN) ::
     >   NIJ(2,NPATCH)         ! Patch dimensions

      REAL, INTENT (IN) ::
     >   XYZ(3,*)              ! Surface point coordinates (packed)

      INTEGER, INTENT (IN) ::
     >   IP(NPATCH),           ! Indices in XYZ(1:3,*) of the start of
     >   ICOMPONENT(NPATCH)    ! each patch, and geometry component #s

      REAL, INTENT (IN) ::
     >   REF_LEN               ! Reference length: cell edges less than
                               ! REF_LEN * 1.E-15 are considered to be
                               ! collapsed
      INTEGER, INTENT (OUT) ::
     >   NPTS                  ! Total mumber of points in XYZ(3,*):
                               ! all must be included if the triangulation
                               ! is written to disk in FAST format
      INTEGER, INTENT (OUT) ::
     >   NFACES                ! Number of triangles (tet. faces) found
                               ! (not counting collapsed triangles)
      INTEGER, INTENT (OUT) ::
     >   IV(3,*)               ! Pointers to the three vertices of
                               ! each triangle in XYZ(1:3,*)
      INTEGER, INTENT (OUT) ::
     >   ICOMP(*)              ! Flag for each triangle, set to the
                               ! appropriate ICOMPONENT element

!     Local constants:

      REAL, PARAMETER ::
     >   ZERO = 0.

!     Local variables:

      INTEGER
     >   I, I0, I1, I2, I3, I4, IC, J, J0, L, N, NI, NJ, NTRI

      REAL
     >   DIAG1, DIAG2, EDGE1, EDGE2, EDGE3, EDGE4, TOL

!     Execution:

      TOL  = 1.E-30 * (REF_LEN ** 2)
      N    = NPATCH
      NPTS = IP(N) + NIJ(1,N) * NIJ(2,N) - 1 ! Possibly IP(NPATCH+1)
      NTRI = 0

      DO N = 1, NPATCH

         I0 = IP(N) - 1 ! Offset in XYZ(1:3,*) for (i,j) elements of patch N
         IC = ICOMPONENT(N)
         NI = NIJ(1,N)
         NJ = NIJ(2,N)

         DO J = 1, NJ - 1

            J0 = I0 + (J - 1) * NI

            DO I = 1, NI - 1 ! For each quad. cell ...

               I1 = J0 + I  ! (i,j)
               I2 = I1 + 1  ! (i+1,j)
               I3 = I2 + NI ! (i+1,j+1)
               I4 = I3 - 1  ! (i,j+1)

               DIAG1 = ZERO
               DIAG2 = ZERO
               EDGE1 = ZERO
               EDGE2 = ZERO
               EDGE3 = ZERO
               EDGE4 = ZERO

               DO L = 1, 3
                  DIAG1 = DIAG1 + (XYZ(L,I1) - XYZ(L,I3)) ** 2
                  DIAG2 = DIAG2 + (XYZ(L,I2) - XYZ(L,I4)) ** 2
                  EDGE1 = EDGE1 + (XYZ(L,I1) - XYZ(L,I2)) ** 2
                  EDGE2 = EDGE2 + (XYZ(L,I2) - XYZ(L,I3)) ** 2
                  EDGE3 = EDGE3 + (XYZ(L,I3) - XYZ(L,I4)) ** 2
                  EDGE4 = EDGE4 + (XYZ(L,I4) - XYZ(L,I1)) ** 2
               END DO

               IF (DIAG1 <= DIAG2) THEN

                  IF (EDGE1 > TOL .AND. EDGE2 > TOL) THEN
                     NTRI = NTRI + 1
                     IV(1,NTRI)  = I1
                     IV(2,NTRI)  = I2
                     IV(3,NTRI)  = I3
                     ICOMP(NTRI) = IC
                  END IF

                  IF (EDGE3 > TOL .AND. EDGE4 > TOL) THEN
                     NTRI = NTRI + 1
                     IV(1,NTRI)  = I3
                     IV(2,NTRI)  = I4
                     IV(3,NTRI)  = I1
                     ICOMP(NTRI) = IC
                  END IF

               ELSE ! DIAG1 > DIAG2

                  IF (EDGE2 > TOL .AND. EDGE3 > TOL) THEN
                     NTRI = NTRI + 1
                     IV(1,NTRI)  = I2
                     IV(2,NTRI)  = I3
                     IV(3,NTRI)  = I4
                     ICOMP(NTRI) = IC
                  END IF

                  IF (EDGE4 > TOL .AND. EDGE1 > TOL) THEN
                     NTRI = NTRI + 1
                     IV(1,NTRI)  = I4
                     IV(2,NTRI)  = I1
                     IV(3,NTRI)  = I2
                     ICOMP(NTRI) = IC
                  END IF

               END IF

            END DO ! Next I

         END DO ! Next J

      END DO ! Next patch

      NFACES = NTRI

      END SUBROUTINE TRIANGULATE_PATCHES
