!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module packages a data structure for grid cells supported by the
!  unstructured data format documented in the FieldView Reference Manual
!  along with a subroutine for reading such a dataset (single zone assumed).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module element_type_module  ! Derived data type for unstructured grid cells

   type element_type
      character (len=3) :: cell_type       ! Element/cell type id, e.g. 'tet'
      integer           :: nodes_per_cell  ! # nodes in element
      integer           :: ncell           ! # cells/elements of this type
      integer, pointer  :: c2n(:,:)        ! c2n(1:nodes_per_cell,i) contain
   end type element_type                   ! the node indices for cell i

   end module element_type_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module FieldView_module

   public :: read_fieldview

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_fieldview (lundat, lunerr, nvars, nnodes, x, y, z,       &
                                 flowvar, varname, nelem, element, ios)

!     Read a 32-bit unformatted file in a FieldView format written by FUN3D for
!     shadowgraph purposes.  The file is expected to contain (x,y,z) coordinates
!     for all grid points (only one unstructured grid in FieldView terminology),
!     along with density and its x/y/z partial derivatives, although the routine
!     will read any number of flow variables.
!
!     Up to 4 cell types are handled: tetrahedral, hexahedral, prism, & pyramid,
!     and they are expected to follow the node numbering convention that appears
!     in the FieldView Reference Manual.
!
!     The cell types are grouped in the sense that the indices for all tet cells
!     are contiguous, with unused integer header words preceding each quadruple.
!     Other cell types are grouped similarly.  The groups may be in any order.
!     Boundary face information follows the x,y,z data but is just skipped here,
!     as is any special boundary-only data following the flow data.
!
!     The flow variables in the file are assumed to be in PLOT3D-type order as
!     needed for possible rotation of the gradient vectors of variable 1 that
!     are (initially) expected to be variables 2, 3, 4 (although this reading
!     routine could handle any number of flow variables).
!
!     The file being read is expected to be open upon entry, and is not closed
!     here.  All storage is allocated here.
!
!     History:
!
!     07/01/10  J.R.Carlson, LaRC  Original reading program adapted from the
!                                  writing.
!     07/06/10  D.A.Saunders, ARC  Trapping read errors can avoid backspacing.
!     07/09/10      "     "        Subroutine version.
!     07/28/10      "     "        Moved the reading routine to a module to
!                                  avoid the explicit interface otherwise reqd.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use element_type_module  ! Derived data type for unstructured grid cells:
                               ! %cell_type,      e.g., 'tet'
                               ! %nodes_per_cell  e.g., m = 4
                               ! %ncell           e.g., n = 232382
                               ! %c2n(1:m,1:n)    indices in x/y/z/flowvar(:)
      implicit none

!     Arguments:

      integer, intent (in)                   :: lundat  ! For file being read
      integer, intent (in)                   :: lunerr  ! For error messages
      integer, intent (out)                  :: nvars   ! # flow variables found
      integer, intent (out)                  :: nnodes  ! # grid pts/nodes found
      real (kind=4), pointer, dimension(:)   :: x, y, z ! Grid point coordinates
      real (kind=4), pointer, dimension(:,:) :: flowvar ! Grid point flow values
      character (len=80),  pointer           :: varname(:) ! Flow variable names
      integer, intent (out)                  :: nelem   ! # element types found
      type (element_type), pointer           :: element(:) ! Grouped cell types
      integer, intent (out)                  :: ios     ! 0 => no error detected

!     Local constants:

      integer, parameter :: maxelem = 4           ! Tet, hex, prism &/or pyramid

!     Local variables:

!     Bit pattern used to detect byte ordering:
!     integer, parameter :: fv_magic = 66051      ! Not clear how we use it

!     Data section id codes:
!     integer, parameter :: fv_nodes      = 1001  ! These aren't used either
!     integer, parameter :: fv_faces      = 1002
!     integer, parameter :: fv_elements   = 1003
!     integer, parameter :: fv_variables  = 1004
!     integer, parameter :: fv_bndry_vars = 1006

!     Element type id codes:
!     integer, parameter :: fv_tet_elem_id   = 1  ! Unused iheader values
!     integer, parameter :: fv_hex_elem_id   = 2 
!     integer, parameter :: fv_prism_elem_id = 3 
!     integer, parameter :: fv_pyra_elem_id  = 4 

      character (len=3)  :: cell_type
      character (len=80) :: bc_name, txt

      integer            :: fv_magic
      integer            :: fv_nodes
      integer            :: fv_faces
      integer            :: fv_elements
      integer            :: fv_variables
      integer            :: fv_bndry_vars
      integer            :: iheader
      integer            :: ver1
      integer            :: ver2
      integer            :: number_of_grids
      integer            :: number_of_boundaries
      integer            :: ibound
      integer            :: ibound_data_flag, irighthand
      integer            :: nv
      integer            :: i_am_zero
      integer            :: nbnds
      integer            :: nbface
      integer            :: ncell, ncell1, ncell2, ncell3, ncell4
      integer            :: nodes_per_cell
      integer            :: i, j
      integer            :: ielem

      real (kind=4)      :: simulation_time  ! These could be output arguments
      real (kind=4)      :: xmach
      real (kind=4)      :: alpha
      real (kind=4)      :: re

!     Execution:

!     Read magic number to verify byte ordering (but we don't actually use it):

      read (lundat, iostat=ios) fv_magic
      if (ios /= 0) then
         write (lunerr, '(a, i6)') 'Error reading fv_magic. ios:', ios
         go to 99  ! Single stop philosophy
      end if
      write (lunerr, '(a, i10)') 'fv_magic:', fv_magic

!     Read header:

      read (lundat, iostat=ios) txt
      if (ios /= 0) then
         write (lunerr, '(a, i6)') 'Error reading header. ios:', ios
         go to 99
      end if
      write (lunerr, '(2a)') 'header: ', trim (txt)

!     Read FieldView version:

      read (lundat, iostat=ios) ver1, ver2
      if (ios /= 0) then
         write (lunerr, '(a, i6)') 'Error reading FieldView version. ios:', ios
         go to 99
      end if
      write (lunerr, '(a,i1,a,i1)') 'FieldView version #: ',ver1, '.',ver2

!     Time, free stream Mach number, angle of attack, Reynolds #:

      read (lundat, iostat=ios) simulation_time, xmach, alpha, re
      if (ios /= 0) then
         write (lunerr, '(a, i6)') 'Error reading t/M/Alpha/Re. ios:', ios
         go to 99
      end if
      write (lunerr, '(a, 1p, 4e16.7)') &
         't/M/Alpha/Re: ', simulation_time, xmach, alpha, re

!     Number of grids (assumed to be 1):

      read (lundat, iostat=ios) number_of_grids
      if (ios /= 0) then
         write (lunerr, '(a, i6)') 'Error reading # grids. ios:', ios
         go to 99
      end if
      write (lunerr, '(a, i4)') '# grids:', number_of_grids

!     Boundary table:

      read (lundat) number_of_boundaries
      if (ios /= 0) then
         write (lunerr, '(a, i6)') 'Error reading # boundaries. ios:', ios
         go to 99
      end if
      write (lunerr, '(a, i4)') '# boundaries:', number_of_boundaries

      write (lunerr, '(a)') ' Flag Hand    Name'
      do ibound = 1, number_of_boundaries
         read (lundat, iostat=ios) ibound_data_flag, irighthand, bc_name
         if (ios /= 0) then
            write (lunerr, '(a, i6)') 'Error reading boundary flags. ios:', ios
            go to 99
         end if
         write (lunerr, '(2i5, 4x, a)') &
            ibound_data_flag, irighthand, trim (bc_name)
      end do

!     Number of volume flow variables:

      read (lundat, iostat=ios) nvars
      if (ios /= 0) then
         write (lunerr, '(a, i6)') 'Error reading # flow variables. ios:', ios
         go to 99
      end if
      write (lunerr, '(a, i4)') '# flow variables:', nvars

      allocate (varname(nvars), stat=ios)
      if (ios /= 0) then
         write (lunerr, '(a, i6)') 'Error allocating flow var. names. ios:', ios
         go to 99
      end if

!     Flow variable names:

      write (lunerr, '(a)') 'Flow variable names:'
      do nv = 1, nvars
         read (lundat, iostat=ios) varname(nv)
         if (ios /= 0) then
            write (lunerr, '(a, i6)') 'Error reading variables names. ios:', ios
            go to 99
         end if
         write (lunerr, '(i2, 2x, a)') nv, trim (varname(nv))
      end do

!     Number of special boundary face flow variables (none expected):

      read (lundat, iostat=ios) i_am_zero
      if (ios /= 0) then
         write (lunerr, '(a, i6)') &
            'Error reading # special boundary face variables. ios:', ios
         go to 99
      end if

!     x, y, z coordinates for each node:

      read (lundat, iostat=ios) fv_nodes, nnodes
      if (ios /= 0) then
         write (lunerr, '(a, i6)') 'Error reading # nodes. ios:', ios
         go to 99
      end if
      write (lunerr, '(a, 2i10)') 'fv_nodes, nnodes:', fv_nodes, nnodes

      allocate (x(nnodes), y(nnodes), z(nnodes), stat=ios)
      if (ios /= 0) then
         write (lunerr, '(a, i6)') 'Error allocating node (x,y,z)s. ios:', ios
         go to 99
      end if

      read (lundat, iostat=ios) x, y, z
      if (ios /= 0) then
         write (lunerr, '(a, i6)') 'Error reading node (x,y,z)s. ios:', ios
         go to 99
      end if

!     Boundary face-to-node connectivity (skipped here):

      write (lunerr, '(a)') 'fv_faces  boundary  # boundary faces'

      do ibound = 1, number_of_boundaries

         read (lundat, iostat=ios) fv_faces, nbnds, nbface
         if (ios /= 0) then
            write (lunerr, '(a, i6)') &
               'Error reading face-to-node connectivity. ios:', ios
            go to 99
         end if
         write (lunerr, '(i8, 2i10)') fv_faces, nbnds, nbface

!!!      allocate (f2n(4,nbface))

!!!      read (lundat, iostat=ios) f2n
         read (lundat, iostat=ios)
         if (ios /= 0) then
            write (lunerr, '(a, i6)') &
!!!            'Error reading face-to-node connectivity. ios:', ios
               'Error skipping face-to-node connectivity. ios:', ios
            go to 99
         end if
!!!      deallocate (f2n)

      end do

!     Read cell-to-node connectivity for flow volume cells (elements).
!     Assume no more than 4 element types (but not all are necessarily present).
!     Assume a read error means no more element types, and that the
!     "fv_variables" code has been encountered instead.

      allocate (element(maxelem), stat=ios)  ! No good way to size it exactly

      if (ios /= 0) then
         write (lunerr, '(a, i5)') 'Error allocating element(:); ios:', ios
         go to 99
      end if

      nelem = 0

      do i = 1, maxelem  ! At most this many cell types, but some may be absent

         read (lundat, iostat=ios) fv_elements, ncell1, ncell2, ncell3, ncell4

         if (ios /= 0) then
            write (lunerr, '(a, i5, a)') &
               'Error reading cell-type header; fv_elements:', fv_elements, &
               '; proceeding to flow variables.'
            exit
         else
            nelem = nelem + 1
            if (ncell1 > 0) then
               ncell = ncell1;  cell_type = 'tet';  nodes_per_cell = 4
            else if (ncell2 > 0) then
               ncell = ncell2;  cell_type = 'hex';  nodes_per_cell = 8
            else if (ncell3 > 0) then
               ncell = ncell3;  cell_type = 'prz';  nodes_per_cell = 6
            else if (ncell4 > 0) then
               ncell = ncell4;  cell_type = 'pyr';  nodes_per_cell = 5
            end if

            element(nelem)%ncell          = ncell
            element(nelem)%cell_type      = cell_type
            element(nelem)%nodes_per_cell = nodes_per_cell

            write (lunerr, '(i1, 3a, i10)') &
               nelem, ': # ', cell_type, ' cells:', ncell

            allocate (element(nelem)%c2n(nodes_per_cell,ncell), stat=ios)

            if (ios /= 0) then
               write (lunerr, '(a, /, a, 3i10)') &
                  'Error allocating flow connectivity arrays.', &
                  'ios, # nodes per cell, # cells:', ios, nodes_per_cell, ncell
               go to 99
            end if

            read (lundat, iostat=ios) &
               (iheader, element(nelem)%c2n(:,j), j = 1, ncell)

            if (ios /= 0) then
               write (lunerr, '(a, /, a, 3i10)') &
                  'Error reading flow connectivity arrays.', &
                  'ios, # nodes per cell, # cells:', ios, nodes_per_cell, ncell
               go to 99
            end if

         end if

      end do  ! Next element type

!     Show the node numbers for the first cell of each element type found:

      do ielem = 1, nelem
         write (lunerr, '(3a, /, 8i9)') 'Nodes for the first ', &
            element(ielem)%cell_type, ' cell:', element(ielem)%c2n(:,1)
      end do

!     Volume data proper:

!!    read  (lundat, iostat=ios) fv_variables  ! This was error-trapped above

!!    write (lunerr, '(a, i6)') 'fv_variables ', fv_variables

      allocate (flowvar(nnodes,nvars), stat=ios)  ! PLOT3D function order

      if (ios /= 0) then
         write (lunerr, '(a, /, a, 3i10)') &
            'Error allocating flow variable arrays.', &
            'ios, # variables, # nodes:', ios, nvars, nnodes
         go to 99
      end if

      read (lundat, iostat=ios) flowvar

      if (ios /= 0) then
         write (lunerr, '(a, /, a, 3i10)') &
            'Error reading flow variables.', &
            'ios, # variables, # nodes:', ios, nvars, nnodes
         go to 99 
      end if

!     Boundary data (none expected):

      read  (lundat, iostat=ios) fv_bndry_vars
      write (lunerr, '(a, i6)') 'fv_bndry_vars:', fv_bndry_vars

99    continue

      end subroutine read_fieldview

   end module FieldView_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
