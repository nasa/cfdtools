!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine nearest_cell_point (itype, nvert, xvert, xtarget, xnear, coefs, &
                                  dsq)

!  One-liner:  Locate the nearest point in or on a grid cell of indicated type
!
!  Overview:
!
!     This is a generalization of the earlier utilities NEAREST_EDGE_POINT,
!  NEAREST_TRI_POINT, and NEAREST_TET_POINT as needed for rapid searching of
!  unstructured grids with mixed cell types, prompted by the US3D flow solver.
!  Certain US3D specifics here could be adjusted for some other nomenclature
!  if necessary.
!
!     The rapid part is done well by ADT techniques (Alternating Digital Trees).
!  An ADT search package is typically not affected by whether a target point is
!  within some mesh cell or not.  Apart from avoiding expensive testing of the
!  bulk of the cells via tree structures, the strategy for each cell that IS
!  tested is to return interpolation coefficients corresponding to the nearest
!  point not outside that cell.  Then, even if no cell is found to contain a
!  given target point, the search returns the best possible cell along with
!  interpolation coefficients to go with it that do not imply extrapolation.
!  This is typically what is really wanted in the presence of unintended
!  geometric discrepancies, or where best mesh boundary values are appropriate.
!
!     If a target point is not outside the current cell, then the distance to
!  the interpolated point is zero and the cell has been found.  More often, the
!  target point is outside a cell being tested.  For geometric applications at
!  least, perpendicular distances to cell boundaries are then clearly the
!  correct measure of closeness.  Sometimes the foot of the perpendicular to a
!  cell boundary lies within the finite edge or face.  More often, it is not and
!  it needs to be adjusted to provide the boundary point closest to the target.
!
!     Each combines an initial calculation of linear interpolation coefficients
!  with the adjustment that is most commonly needed to produce an interpolated
!  point not outside the given cell.  Note that for quadrilateral surface cells
!  hexahedral structured volume grid cells, the standard coefficients are non-
!  linear.  The author chose to apply the same basic linear method to all cells
!  out of curiosity, knowing that the nonlinear methods could always be reverted
!  to for quads and hexes if necessary.  Thus, at least initially, the computed
!  coefficients for any cell type are such that
!
!    Sum (i = 1:nvert) coefs(i)*xvert(i) = xtarget (or nearest cell point to it)
!  & Sum (i = 1:nvert) coefs(i)          = 1
!
!     Note that all cases of the linear method involve 3 equations for x, y, z,
!  and a fourth for the sum of the coefficients.  These can always be reduced to
!  just 3 equations for n - 1 coefficients, from which the nth coefficient is
!  deduced since the sum is 1 (where n = # vertices).  Thus, the system may be
!  overdetermined (triangle case) or square (tet and linear quad cases) or
!  underdetermined (all others).  Underdetermined least squares solutions for
!  the first step of the calculation provide the solution with shortest length
!  (smallest 2-norm), but then the second step of the algorithm may adjust
!  the first-step result if the target is outside the cell being processed.
!
!     Also, in the interest of efficiency, degenerate cells are not handled at
!  this low level, although a collapsed vertex may not cause a problem.
!
!  History:
!
!     08/27/04  D.A.Saunders  Initial implementation of NEAREST_TET_POINT.
!     06/06/13    "     "     Adapted NEAREST_TET_POINT as NEAREST_CELL_POINT.
!
!  Author:
!
!     David Saunders, ERC Inc./NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in) :: &            ! Fluent nomenclature for cell types:
      itype                             !   type               # vertices
                                        !   1 = triangle            3
                                        !   2 = tetrahedron         4
                                        !   3 = quadrilateral       4
                                        !   4 = hexahedron          8
                                        !   5 = pyramid             5
                                        !   6 = prism               6
   integer, intent (in) :: &
      nvert                             ! Number of vertices for this type

   real, intent (in) :: &
      xvert(3,nvert)                    ! (x,y,z)s for each vertex

   real, intent (in) :: &
      xtarget(3)                        ! Target point (x,y,z) to find the
                                        ! best cell for
   real, intent (out), dimension (3) :: &
      xnear                             ! The point nearest to xtarget that is
                                        ! not outside the cell
   real, intent (out) :: &
      coefs(nvert)                      ! Interpolation coefficients such that
                                        !   c1 x1 + c2 x2 + ... + cn xn = xnear;
                                        !   c1    + c2    + ... + cn    = 1
   real, intent (out) :: &
      dsq                               ! The corresponding shortest distance
                                        ! between xtarget and xnear, squared
!  Execution:

   select case (itype)

     case (1)  ! Triangle

       call nearest_tri_point (xvert(:,1), xvert(:,2), xvert(:,3), xtarget, &
                               xnear, coefs(1), coefs(2), coefs(3), dsq)

     case (2)  ! Tetrahedron

       call nearest_tet_point (xvert(:,1), xvert(:,2), xvert(:,3), xvert(:,4),&
                               xtarget, xnear, coefs(1), coefs(2), coefs(3),  &
                               coefs(4), dsq)

     case (3)  ! Quadrilateral

       call nearest_quad_point (xvert, xtarget, xnear, coefs, dsq)

     case (4)  ! Hexahedron

       call nearest_hex_point (xvert, xtarget, xnear, coefs, dsq)

     case (5)  ! Pyramid

       call nearest_pyramid_point (xvert, xtarget, xnear, coefs, dsq)

     case (6)  ! Prism

       call nearest_prism_point (xvert, xtarget, xnear, coefs, dsq)

   end select

   end subroutine nearest_cell_point
