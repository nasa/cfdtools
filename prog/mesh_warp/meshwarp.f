c*******************************************************************************
c
c  PROGRAM: MESHWARP.v1
c
c  PURPOSE:
c    MESHWARP, a member of the GMORPH suite of utilities, is used to morph an
c    existing surface or volume mesh to reflect design shape changes in an
c    underlying geometry definition.  The warping procedure used by this utility
c    is based on the methodology developed in the SYN107-MB design code.
c    MESHWARP may be used to update unstructured or structured, surface or
c    volume meshes based on unstructured or structured geometry definitions.
c
c***
c  INPUTS:
c     Six sets of data are required as inputs to the MESH_WARP utility:
c     the original geometry representation used to generate the surface mesh,
c     the modified geometry representation, the surface or volume mesh points,
c     the mesh connectivity file, a relational mapping between the surface mesh
c     points and the geometry definition produced by the UV_MAP utility, and
c     the meshwarp.inp control file.  Except for the control file, which must
c     be named meshwarp.input, all file names and formats are specified in
c     meshwarp.inp.  All unformatted files must contain 64-bit words.
c
c  GEOMETRY DATABASES:
c     A discretized representation of the geometry must be supplied in
c     either structured or unstructured formats.  MESHWARP supports structured
c     definitions in either PLOT3D or AEROSURF formats in 3-D, multiple grid,
c     formatted or unformatted, IBLANK or no IBLANK styles.
c     Unstructured geometry databases must be in FAST format in 3-D, single
c     grid, formatted or unformatted style.
c
c     The original and modified geometry definitions must be of the same
c     format and consistent with the geometry used by UV_MAP to create the
c     relational mapping.
c
c  VOLUME/SURFACE MESH POINTS:
c     The volume or surface mesh points to be warped may also be supplied in
c     structured or unstructured formats.  Unstructured surface/volume meshes
c     must be in FAST format with 3-D, single grid, formatted or unformatted
c     attributes. The IFLAG associated with each triangle is used to group
c     the unstructured points into groups analogous to a block in a structured
c     mesh.  Structured meshes may be in PLOT3D format in 3-D, multiple grid,
c     formatted or unformatted, IBLANK or no IBLANK styles.
c
c     Currently, an unstructured volume warping method has not been implemented.
c     Using MESHWARP on an unstructured volume mesh will result in only the
c     surface points defined in the (u,v) map being translated to follow the
c     geometry perturbation.  This feature may be implemented in the future.
c
c     The volume meshes must be of the same format and structure as that used
c     by UV_MAP to create the relational mapping.
c
c  CONNECTIVITY FILE:
c     When using structured meshes, the connectivity file has the same format
c     and interpretation as that used by SYN107-MB.
c
c     Unstructured meshes require one boundary condition line for each patch
c     grouping specified through the IFLAG mechanism described in the section
c     on volume mesh inputs.  The boundaries of the unstructured mesh patches
c     are not ensured to remain coincident during warping.
c
c     Structured faces or unstructured patches for which a relational map has
c     been constructed should be activated in the HYPERFACE section of the
c     connectivity file by block number and face.  Structured surface meshes
c     should use face 5 for activation while volume meshes behave as a standard
c     SYN107-MB face numbering.  Unstructured patches should specify the block
c     field with the patch number and the face field as one.
c
c     Design faces may be specified in the connectivity file to force the motion
c     of one block face (the "design" face) relative to a face on a geometry
c     surface (the "hyperface").  Modification of the design face specification
c     does not require UV_MAP to be rerun.
c
c  RELATIONAL MAPPING:
c     A relational mapping produced by UV_MAP must be supplied in either
c     absolute, (p,q), or parameterized, (u,v), formats.  See the UV_MAP
c     documentation for more details.  The original geometry, volume/surface
c     mesh, and connectivity files used with MESHWARP should be the same as
c     those used with UV_MAP.  The (u,v) map has an implied order of block
c     faces on the surface - that of the "activation" list in the connectivity
c     file.
c
c  meshwarp.inp:
c     The control file description is elaborated in a separate section below.
c
c***
c  OUTPUT:
c     MESHWARP produces a single output at the completion of execution:
c     a perturbed volume/surface mesh in the same format as the unperturbed
c     volume/surface mesh.
c
c***
c  CONTROL FILE:
c     The control file, named meshwarp.inp, contains information on the I/O
c     file names/formats.
c
c     The file formats are described below with the appropriate character
c     string specifier.  The input strings are case insensitive.
c
c     TYPE           STRING           DESCRIPTION
c
c     FORM/UNFORM    form/unform      Formatted or unformatted data
c     STRCT/UNSTRCT  strct/unstrct    Structured or unstructured
c                                     geometry database or mesh
c     NOBLANK/BLANK  noblank/blank    IBLANKing present in file
c     AERO/P3D       aero/p3d         A structured geometry database
c                                     may be in AEROSURF or PLOT3D format
c     ABSOLUTE/PARAM  absolute/param  A relational mapping may be constructed
c                                     based on triangle indices or patch
c                                     arc-length parameterization
c
c     A sample input file is reproduced below.  The leading 'c' should be
c     deleted from an actual meshwarp.inp file.
c
c--Cut here--
c -----------------------------------------------------
c                  INPUT FILES
c -----------------------------------------------------
c ORIGINAL_GEO_DATA   FORM/UNFORM STRCT/UNSTRCT NOBLANK/BLANK  AERO/P3D
c 1_geometry.xyz         form        strct         noblank     aero
c
c MODIFIED_GEO_DATA  !The file format of the modified geometry is taken
c 2_geometry.xyz     !to be the same as the original geometry file.
c
c ORIGINAL_VOLUME_MESH FORM/UNFORM STRCT/UNSTRCT NOBLANK/BLANK
c 1_volume.xyz          unform       strct         blank
c
c CONNECTIVITY_FILE
c a.conn
c
c UV_MAP              FORM/UNFORM ABSOLUTE/PARAM
c uv.map                 form        absolute
c
c -----------------------------------------------------
c                  OUTPUT FILES
c -----------------------------------------------------
c MODIFIED_VOLUME_MESH FORM/UNFORM STRCT/UNSTRCT NOBLANK/BLANK
c 2_volume.xyz          unform       strct         blank
c
c--End Cut Here-
c****
c
c  PACKAGES, MODULES, and SOURCE FILES:
c     Input/Output library:
c       cfd_io_pacakge.f
c       cfd_io_package.mod
c     Dynamic binary tree management library:
c       tree_package.f
c       tree_package.mod
c     Triangle sorting library used with tree_package:
c       triangle_sort.f
c       triangle_sort.mod
c     Numutil library routines
c       tools.f (gathered for portability)
c     MPI overloaded routines for non-MPI environments:
c       mpi_dummy.f
c     Flush for VMS:
c       flush_vms.f
c
c***
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c
c     07/17/03 DAS  The MAPXYZ step made the undocumented assumption that the
c                   original and perturbed geometry patches had the same
c                   dimensions.  Even so, the perturbed patch should still be
c                   searched as well as the original for each (u,v).  Now,
c                   the patch dimensions do not have to match for this stand-
c                   alone tool (though they may well have to in an optimization
c                   application).
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c*******************************************************************************

c*******************************************************************************

c     Module Definitions:

      module Parameters
        integer, parameter ::
     .    MBCRD = 3, MBFACE = 6, MBEDGE = 12, MBCRN = 8,
     .    MLEDGE = 4, MLCRN = 4, MAXECNCT = 7, MAXCRNCNCT = 96,
     .    StdOut = 6, CntrlUnit = 1
      end module Parameters

      module TypeDefs

c--------

        type Geometry

c         Structured Representation
          real, pointer ::
     .      x(:,:)
          integer, pointer ::
     .      iBlank(:),
     .      ipx(:),
     .      ijk(:,:)
          integer ::
     .      nComp,
     .      nPnl,
     .      nPtch

c         Parameterized Representation
          real, pointer ::
     .      u_ptch(:), v_ptch(:)

c         Aerosurf Subpatch Information
          integer, pointer  ::
     .      pnl2Comp(:),
     .      pnl2Ptch(:,:),
     .      pnlNPtch(:),
     .      pnlNSubI(:), pnlNSubJ(:),
     .      pnlISub(:,:), pnlJSub(:,:),
     .      ptch2Pnl(:)

c         Unstructured Representation
          integer, pointer ::
     .      tri(:,:),
     .      iTriUnsrt(:),
     .      iTri2Ptch(:),
     .      iTri2PtchIJ(:,:),
     .      tet(:,:),
     .      ipTriPnl(:)

          integer ::
     .      nPt,
     .      nTri,
     .      nTet

        end type Geometry

c--------

        type AbstractMeshType

c         Structured Representation
          real, pointer ::
     .     x(:,:)
          integer, pointer ::
     .     iBlank(:),
     .     ipx(:),
     .     ijk(:,:)
          integer ::
     .     nBlock

c         Unstructured Representation
          integer, pointer ::
     .      tri(:,:),
     .      iTri(:),
     .      tet(:,:)
          integer ::
     .      nPt,
     .      nTri,
     .      nTet

        end type AbstractMeshType

c--------

        type BCData
          integer ::
     .      nBC
          integer, pointer ::
     .      typ(:,:),
     .      algn(:,:,:),
     .      is(:,:,:),             !Home block transf. matrix
     .      it(:,:,:,:),           !Neighbor block transf. matrix
     .      ir(:,:,:,:)            !Relative transf. matrix
        end type BCData

c--------

        type DesignFace

          integer, pointer ::
c           Source/target block/face specifications:
     .      srcBlock(:), srcFace(:),
     .      trgBlock(:), trgFace(:),

c           Source/target pointers to mesh x array:
     .      dfsrc2xp(:),
     .      dftrg2xp(:),
     .      dftrg2nbxp(:),
     .      ip(:),
     .      ijk(:,:),

c           Averaging corner/design face corner info:
     .      dfc2dfxp(:,:),
     .      dfc2accnct(:,:),
     .      dfc2ac(:,:),
     .      accnct(:),
     .      msc2ac(:),
     .      ac2igb(:,:),
     .      c_type(:,:),

c           Averaging edges/design face edge info:
     .      dfe2dfxp(:),
     .      dfeip(:,:),
     .      dfenPt(:,:),
     .      dfe2ae(:,:),
     .      dfe2aecnct(:,:),
     .      dfe2aexp(:),
     .      aecnct(:),
     .      aeip(:),
     .      mse2ae(:),
     .      ae2igb(:,:),
     .      e_type(:,:)

          integer ::
     .       acct, aect,
     .       nDF, nDFEpts

        end type DesignFace

c--------

        type HyperFace
          integer, pointer ::
c           Hyperface block/face specification:
     .      ihyp(:,:),
     .      hf2xp(:),
     .      ip(:)

          integer ::
     .      nHF
        end type HyperFace

c--------

        type MasterData

c         Master edge data:
          integer, pointer ::
     .      msep(:),
     .      msecnct(:),
     .      mse2igep(:,:),
     .      mse2igb(:,:),
     .      mse2ige(:,:),
     .      mse2ige_dir(:,:),
     .      mse2msc(:,:),
     .      ige2mse(:,:),
     .      ige2mse_dir(:,:),
     .      mse2ige_actHF(:),
     .      mse2ige_actDF(:),
     .      ige2mse_actHF(:,:),
     .      ige2mse_actDF(:,:)

          integer ::
     .      msect,
     .      mseptct

c         Master corner data
          integer, pointer ::
     .      msccnct(:),
     .      msc2igb(:,:),
     .      msc2ige(:,:),
     .      msc2igc(:,:),
     .      msc2igcp(:,:),
     .      igc2msc(:,:),
     .      msc2igc_actHF(:),
     .      msc2igc_actDF(:)

          integer ::
     .      mscct

c         Active block, face, edge lists:
          integer, pointer ::
     .      igb_act(:),
     .      igf_act(:,:),
     .      ige_act(:,:)

        end type MasterData

c--------

        type UnitBlockNumbering
          integer ::
     .      sng(4,4),    !local edge index ranges
     .      f2ei(4,6),   !block face composed of 4 block edges
     .      f2ci(4,6),   !block face composed of 4 block corners
     .      e2fi(2,12),  !block edge member of 2 block faces
     .      e2fj(2,12),  !block edge is an edge in the local face
     .      e2ci(2,12),  !block edge is composed of 2 block corners
     .      c2ei(3,8),   !block corner forms an end of 3 block edges
     .      c2ej(3,8)    !block corner is a beginning or end
                         !  of the three block edges.
        end type UnitBlockNumbering

c--------

        type UVMap
          real, pointer ::
     .      uv(:,:)

          integer, pointer ::
     .      iuv(:),
     .      ijk(:,:),
     .      ip(:)

          integer ::
     .      nHyp
        end type UVMap

c--------

      end module TypeDefs

c***********************************************************************
c
      program MESHWARP
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c
c***********************************************************************

c  Modules:

      use Parameters
      use TypeDefs
      use cfd_io_package

      implicit none

c  Local variables:

      type (AbstractMeshType) ::
     .  mesh0, meshN  !Initial and perturbed volume meshes.

      type (Geometry) ::
     .  geom0, geomN  !Initial and perturbed geometries.

      type (UVMap) ::
     .  uv            !(u,v) map relating geometry to mesh.

      type (BCData) ::
     .  bc            !Boundary condition data for mesh.

      type (DesignFace) ::
     .  df            !Design face specification.

      type (HyperFace) ::
     .  hf            !Hyperface specification.

      type (MasterData) ::
     .  mstr          !Master edge/corner lists.

      type (ioMeshFileType) ::
     .  iFileMesh0,
     .  iFileMeshN,
     .  iFileGeom0,
     .  iFileGeomN,
     .  iFileConn,
     .  iFileUV

      character (len=12) ::
     .  meshType,
     .  geomType,
     .  connType,
     .  uvType

      real ::
     .  mpi_wtime

      real ::
     .  cpu1, cpu2

c  Execution:

      cpu1 = mpi_wtime()

c     Read input control:
      call readInputControl (iFileMesh0, iFileMeshN, meshType,
     .                       iFileGeom0, iFileGeomN, geomType,
     .                       iFileConn, connType,
     .                       iFileUV, uvType)

c     Read first geometry database:
      call prepareGeometryData (geom0, iFileGeom0, geomType, uvType)

c     Read modified geometry database:
      call prepareGeometryData (geomN, iFileGeomN, geomType, uvType)

c     Read the initial volume mesh:
      call readMeshData (mesh0, iFileMesh0, meshType)

c     Copy the initial mesh into the soon to be perturbed mesh:
      call duplicateMesh (mesh0, meshN, meshType)

c     Read connectivity to obtain boundary conditions, hyperfaces,
c     and design faces:
      call readConnFile (bc, df, hf, mesh0, iFileConn, connType)

c     Read uv-map:
      call readUVFile (uv, hf%nHF, iFileUV, uvType)

c     Create master edge/corner lists:
      if (meshType == 'strct') then
        call createMasterInfo (mstr, hf, df, mesh0, bc)
      else
        call createMasterInfoUnstr (hf, mesh0, uv)
      end if

c     Perturb all hyperfaces:
      call mapXYZ (meshN, mesh0, geomN, geom0, hf, uv, uvType)

c     Perturb blocks based on hyperface perturbation:
      if (meshType == 'strct')
     .  call warpmb (meshN, mesh0, mstr, df)

c     Write perturbed mesh to file:
      call writeMeshData (meshN, iFileMeshN, meshType)

      cpu2 = mpi_wtime ()
      write (StdOut,'(1x,a,F10.2)') 'Execution time: ', cpu2 - cpu1

      end program MESHWARP

c+-----------------------------------------------------------------------
c
      subroutine readInputControl (iFileMesh0, iFileMeshN, meshType,
     .                             iFileGeom0, iFileGeomN, geomType,
     .                             iFileConn, connType,
     .                             iFileUV, uvType )
c
c  Description:  Read the I/O file names and formats from the
c  meshwarp.input control file.
c
c  meshType is set to one of unstrct/strct
c  geomType is set to one of unstrct/aero/p3d
c  connType is set to syn107.  No other connectivity files are supported.
c  uvType is set to one of absolute/param
c
c  The fields of the ioMeshFileType as reproduced below from the
c  cfd_io_package are set to the appropriate values.
c
c  ndim is set to 3, status is set to old for input files and unknown for
c  output files.
c
c  triples is set to .T.
c
c  type, public :: ioMeshFileType
c     character(len=80) :: name    !File name
c     character(len=11) :: form    !File format, 'unformatted'|'formatted'
c     character(len=7)  :: status  !File status, 'old'|'new'|'unknown'
c     logical           :: mg      !Multiple grid format, .T.|.F.
c     logical           :: triples !Shape of returned array:
c                                  !x(ndim,npts) -> .T.
c                                  !x(npts,ndim) -> .F.
c     integer           :: ndim    !Spatial dimensions, N -> N-D
c     logical           :: iblank  !iblanking .T.
c  end type ioMeshFileType
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c
c------------------------------------------------------------------------

c  Modules:

      use Parameters
      use cfd_io_package

      implicit none


c  Arguments:

      type (ioMeshFileType) ::
     .  iFileMesh0, iFileMeshN,
     .  iFileGeom0, iFileGeomN,
     .  iFileConn, iFileUV

      character (*) ::
     .  meshType, geomType,
     .  connType, uvType

c  Local variables:

      logical :: exist

      integer :: i

      character (len=12) ::
     .  input_file = 'meshwarp.inp'
      character (len=16) ::
     .  routine = 'readInputControl'
      character (len=7) ::
     .  structure, blanking, aerop3d

c  Execution:

      write (StdOut,*) 'Reading control file...'
      call flush (StdOut)

      inquire (file=input_file, exist=exist)
      if (.not. exist) goto 100

      open (unit=CntrlUnit, file=input_file, form='formatted')

c     Brush past INPUT comment lines:
      do i = 1, 3
        read (CntrlUnit,*)
      end do

c     Read original geometry data file name/format
      read (CntrlUnit,*)
      read (CntrlUnit,*) iFileGeom0%name, iFileGeom0%form,
     .                   structure, blanking, aerop3d
      read (CntrlUnit,*) !Clear blank line after original geometry

      call locase (iFileGeom0%form)
      call locase (structure)
      call locase (blanking)
      call locase (aerop3d)

      select case (iFileGeom0%form)
        case ('form')
          iFileGeom0%form = 'formatted'
        case ('unform')
          iFileGeom0%form = 'unformatted'
        case default
          goto 101
      end select

      select case (structure)
        case ('strct')
          select case (aerop3d)
            case ('aero')
              geomType = 'aero'
            case ('p3d')
              geomType = 'p3d'
            case default
              goto 103
          end select
        case ('unstrct')
          geomType = 'unstrct'
        case default
          goto 102
      end select

      select case (blanking)
        case ('blank')
          iFileGeom0%iblank = .TRUE.
        case ('noblank')
          iFileGeom0%iblank = .FALSE.
        case default
          goto 104
      end select

      iFileGeom0%mg      = .TRUE.
      iFileGeom0%status  = 'old'
      iFileGeom0%triples = .TRUE.
      iFileGeom0%ndim    = MBCRD

c     Read modified geometry data file name/format
c     The modified geometry geomType is taken to be the same as the original.
      read (CntrlUnit,*)
      read (CntrlUnit,*) iFileGeomN%name
      read (CntrlUnit,*) !Clear blank line after modified geometry

      iFileGeomN%form   = iFileGeom0%form
      iFileGeomN%iblank = iFileGeom0%iblank
      iFileGeomN%mg      = .TRUE.
      iFileGeomN%status  = 'old'
      iFileGeomN%triples = .TRUE.
      iFileGeomN%ndim    = MBCRD

c     Read original mesh data file name/format
      read (CntrlUnit,*)
      read (CntrlUnit,*) iFileMesh0%name, iFileMesh0%form,
     .                   structure, blanking
      read (CntrlUnit,*) !Clear blank line after original mesh

      call locase (iFileMesh0%form)
      call locase (structure)
      call locase (blanking)

      select case (iFileMesh0%form)
        case ('form')
          iFileMesh0%form = 'formatted'
        case ('unform')
          iFileMesh0%form = 'unformatted'
        case default
          goto 101
      end select

      select case (structure)
        case ('strct')
          meshType = 'strct'
        case ('unstrct')
          meshType = 'unstrct'
        case default
          goto 102
      end select

      select case (blanking)
        case ('blank')
          iFileMesh0%iblank = .TRUE.
        case ('noblank')
          iFileMesh0%iblank = .FALSE.
        case default
          goto 104
      end select

      iFileMesh0%mg      = .TRUE.
      iFileMesh0%status  = 'old'
      iFileMesh0%triples = .TRUE.
      iFileMesh0%ndim    = MBCRD

c     Read connectivity file name
      read (CntrlUnit,*)
      read (CntrlUnit,*) iFileConn%name
      read (CntrlUnit,*) !Clear blank line after connectivity file

      iFileConn%form = 'formatted'

      connType = 'syn107'
      iFileConn%status  = 'old'

c     Read UV_MAP file name
      read (CntrlUnit,*)
      read (CntrlUnit,*) iFileUV%name, iFileUV%form, uvType
      read (CntrlUnit,*) !Clear blank line after UV_MAP file

      call locase (iFileUV%form)
      call locase (uvType)

      select case (iFileUV%form)
        case ('form')
          iFileUV%form = 'formatted'
        case ('unform')
          iFileUV%form = 'unformatted'
        case default
          goto 101
      end select

      select case (uvType)
        case ('absolute')
        case ('param')
        case default
          goto 105
      end select

      iFileUV%mg      = .true.
      iFileUV%triples = .true.
      iFileUV%status  = 'old'
      iFileUV%ndim    = 2
      iFileUV%iblank  = .true.

c     Brush past OUTPUT comment lines:
      do i = 1, 3
        read (CntrlUnit,*)
      end do

c     Read modified volume mesh name/format
c     Force all options to be the same as the input mesh file
c     except for the file format.
      read (CntrlUnit,*)
      read (CntrlUnit,*) iFileMeshN%name, iFileMeshN%form

      call locase (iFileMeshN%form)

      select case (iFileMeshN%form)
        case ('form')
          iFileMeshN%form = 'formatted'
        case ('unform')
          iFileMeshN%form = 'unformatted'
        case default
          goto 101
      end select

      iFileMeshN%iblank  = iFileMesh0%iblank
      iFileMeshN%mg      = iFileMesh0%mg
      iFileMeshN%status  = 'unknown'
      iFileMeshN%triples = iFileMesh0%triples
      iFileMeshN%ndim    = MBCRD

      return

c  Error handling:

 100  call meshwarpGracefulExit (routine,
     .  'meshwarp.inp control file does not exist.')

 101  call meshwarpGracefulExit (routine,
     .  'File format must be one of: form/unform.')

 102  call meshwarpGracefulExit (routine,
     .  'File structure must be one of: strct/unstrct.')

 103  call meshwarpGracefulExit (routine,
     .  'Geometry file package must be one of: aero/p3d.')

 104  call meshwarpGracefulExit (routine,
     .  'Blanking format must be one of: blank/noblank.')

 105  call meshwarpGracefulExit (routine,
     .  'UVType must be one of: absolute/param.')

      end subroutine readInputControl

c---------------------------------------------------------------------------
c
c     Geometry Preparation Routines:  Note that the routines here
c     differ from the UV_MAP preparation in that the component->panel
c     arrays are not initialized.  MESHWARP does not care about components.
c
c+--------------------------------------------------------------------------
c
      subroutine prepareGeometryData (geom, iFileGeom, geomType, uvType)
c
c  Description:  Read geometry data from file in the format specified by
c  geomType.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c---------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs
      use cfd_io_package

      implicit none

c  Arguments:

      type (Geometry), intent(inout) ::
     .  geom

      type (ioMeshFileType), intent(in) ::
     .  iFileGeom

      character (*), intent (in) ::
     .  geomType, uvType

c  Local variables:

      integer ::
     .  ncall = -1,
     .  i, j

      integer, allocatable ::
     .  iTri2PtchTmp(:),
     .  iTri2PtchIJTmp(:,:)

      character (len=19) ::
     .  routine = 'prepareGeometryData'

      logical, save ::
     .  first = .true.

c  Execution:

      if (first) then
        write (StdOut,*) 'Preparing original geometry...'
        first = .false.
      else
        write (StdOut,*) 'Preparing modified geometry...'
      end if

c     Read Geometry Data from file.  If other flavors of
c     database arise in future, add them as more cases here.

      select case (geomType)
        case ('aero')
          call prepareAerosurfGeometry (geom, iFileGeom, uvType)
        case ('p3d')
          call prepareP3DGeometry (geom, iFileGeom, uvType)
        case ('unstrct')
          call prepareUnstructGeometry (geom, iFileGeom)
        case default
          goto 100
      end select

      return

c  Error handling:

 100  call meshwarpGracefulExit (routine,
     .      'Unknown type of geometry file.')

      end subroutine prepareGeometryData

c+---------------------------------------------------------------------
c
      subroutine prepareAerosurfGeometry (geom, iFileGeom, uvType)
c
c  Description:  Read Aerosurf Geometry, construct
c  component=>panel=>patch relationship and triangulate.
c  Assumes that the panels have been pre-ordered into component blocks.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c----------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs
      use cfd_io_package

      implicit none

c  Arguments:

      type (Geometry) ::
     .  geom

      type (ioMeshFileType) ::
     .  iFileGeom

      character (len=*) ::
     .  uvType

c  Local variables:

      integer ::
     .  i, j, l, m, i_sub, j_sub,
     .  ncall = -1,
     .  i_Ptch, i_Pnl, id_Pnl,
     .  nPts, lPts,
     .  l_cnct, l_Pnl, l_Ptch,
     .  id, jd, ip_s, ip_e

      integer, allocatable ::
     .  ijkPtch(:,:), ipPtch(:)

      real, allocatable ::
     .  xPtch(:,:), xtmp(:), ytmp(:), ztmp(:)

      type (ioMeshFileType), pointer :: iFileTmp

      character (len=23) :: routine = 'prepareAerosurfGeometry'

c  Execution:

      allocate (iFileTmp)
      iFileTmp = iFileGeom

c     Read data from file.

      call io_r_aerosurf (geom%x, geom%nPnl, geom%ijk, geom%ipx,
     .                    geom%pnl2Comp, geom%pnlNSubI, geom%pnlNSubJ,
     .                    geom%pnlISub, geom%pnlJSub,
     .                    ncall, iFileTmp)

c     Construct panel=>patch relationship.

c     Count the total number of patches:
      geom%nPtch = 0
      allocate (geom%pnlNPtch(geom%nPnl))
      geom%pnlNPtch = 0

      do i = 1, geom%nPnl
        j = geom%pnlNSubI(i)*geom%pnlNSubJ(i)
        geom%pnlNPtch(i) = j
        geom%nPtch = geom%nPtch + j
      end do

      allocate (geom%ptch2Pnl(geom%nPtch))
      geom%ptch2Pnl = 0
      allocate (geom%pnl2Ptch(geom%nPtch,geom%nPnl))
      geom%pnl2Ptch = 0

c     The x, ijk geometry arrays were read in based on panels.  Now
c     they will be redimensioned/shifted based on patches:

      allocate (ijkPtch(3,geom%nPtch), ipPtch(geom%nPtch+1))

c     Load the ijk/ip of each patch and count the total number of patch points:
      i_Ptch = 0
      nPts = 0
      do i_Pnl = 1, geom%nPnl
        l_cnct = 0
        do j = 1, geom%pnlNSubJ(i_Pnl)
        do i = 1, geom%pnlNSubI(i_Pnl)
          i_Ptch = i_Ptch + 1
          l_cnct = l_cnct + 1
          ipPtch(i_Ptch) = nPts + 1
          ijkPtch(1,i_Ptch) =
     .      geom%pnlISub(i+1,i_Pnl)-geom%pnlISub(i,i_Pnl)+1
          ijkPtch(2,i_Ptch) =
     .      geom%pnlJSub(j+1,i_Pnl)-geom%pnlJsub(j,i_Pnl)+1
          ijkPtch(3,i_Ptch) = 1
          lPts = ijkPtch(1,i_Ptch)*ijkPtch(2,i_Ptch)
          nPts = nPts + lPts
          geom%ptch2Pnl(i_Ptch) = i_Pnl
          geom%pnl2Ptch(l_cnct,i_Pnl) = i_Ptch
        end do
        end do
      end do

      ipPtch(geom%nPtch+1) = nPts+1

c     Allocate x array for patches and load from geometry x array:

      allocate (xPtch(3,nPts))

      i_Ptch = 0
      do i_Pnl = 1, geom%nPnl
        id_Pnl = geom%ijk(1,i_Pnl)
        l_Pnl  = geom%ipx(i_Pnl)
        do j_sub = 1, geom%pnlNSubJ(i_Pnl)
        do i_sub = 1, geom%pnlNSubI(i_Pnl)
          i_Ptch = i_Ptch + 1
          l_Ptch = ipPtch(i_Ptch)
          do j = geom%pnlJSub(j_sub,i_Pnl), geom%pnlJSub(j_sub+1,i_Pnl)
          do i = geom%pnlISub(i_sub,i_Pnl), geom%pnlISub(i_sub+1,i_Pnl)
            m = l_Pnl + (i-1) + id_Pnl*(j-1)
            xPtch(1:3,l_Ptch) = geom%x(1:3,m)
            l_Ptch = l_Ptch + 1
          end do
          end do
        end do
        end do
      end do

c     Copy restructured patches back into geometry arrays:

      deallocate (geom%x, geom%ipx, geom%ijk)

      allocate (geom%x(3,nPts), geom%ipx(geom%nPtch+1),
     .          geom%ijk(3,geom%nPtch), geom%iBlank(nPts))

      geom%x = xPtch
      geom%ipx = ipPtch
      geom%ijk = ijkPtch
      geom%iBlank = 1

      deallocate (xPtch, ipPtch, ijkPtch)

c     Convert to an unstructured mesh:

      if (uvType == 'absolute') then

        call struct2Unstruct (geom)

      else  !parameterized (u,v) map

c       Create parameterization of patches:

        allocate (geom%u_ptch(nPts), geom%v_ptch(nPts))

        do i_Ptch = 1, geom%nPtch

          id = geom%ijk(1,i_Ptch)
          jd = geom%ijk(2,i_Ptch)
          ip_s = geom%ipx(i_Ptch)
          ip_e = geom%ipx(i_Ptch+1)-1

          allocate (xtmp(id*jd), ytmp(id*jd), ztmp(id*jd))

          xtmp = geom%x(1,ip_s:ip_e)
          ytmp = geom%x(2,ip_s:ip_e)
          ztmp = geom%x(3,ip_s:ip_e)

          geom%u_ptch(ip_s) = 0. ! -999. suppresses normalization

          call param2d (id, jd, 1, id, 1, jd,
     .                  xtmp, ytmp, ztmp,
     .                  geom%u_ptch(ip_s), geom%v_ptch(ip_s))

          deallocate (xtmp, ytmp, ztmp)
        end do
      end if

      end subroutine prepareAerosurfGeometry

c+---------------------------------------------------------------------
c
      subroutine prepareP3DGeometry (geom, iFileGeom, uvType)
c
c  Description:  Read Plot3D geometry in blanked or unblanked format,
c  then create an unstructured representation.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c----------------------------------------------------------------------

c  Modules:

      use TypeDefs
      use cfd_io_package

      implicit none

c  Arguments:

      type (Geometry) ::
     .  geom

      type (ioMeshFileType) ::
     .  iFileGeom

      character (len=*) ::
     .  uvType

c  Local variables:

      integer ::
     .  i, j, l, m, i_sub, j_sub,
     .  ncall = -1,
     .  i_Ptch,
     .  nPts,
     .  id, jd, ip_s, ip_e

      real, allocatable ::
     .  xtmp(:), ytmp(:), ztmp(:)

      type (ioMeshFileType), pointer ::
     .  iFileTmp

      character (len=18) :: routine = 'prepareP3DGeometry'

c  Execution:

      allocate (iFileTmp)
      iFileTmp = iFileGeom

c     Read geometry data from file.

      if (iFileGeom%iblank) then
        call io_r_str_x_p3d_ib (geom%x, geom%nPnl, geom%ijk,
     .                          geom%ipx, geom%iBlank, ncall, iFileTmp)
      else
        call io_r_str_x_p3d (geom%x, geom%nPnl, geom%ijk,
     .                       geom%ipx, ncall, iFileTmp )
        allocate (geom%iBlank(geom%ipx(geom%nPnl+1)-1))
        geom%iBlank = 1
      end if

c     Construct panel=>patch relationship:

      allocate (geom%pnlNPtch(geom%nPnl))
      geom%pnlNPtch = 1

      allocate (geom%pnl2Ptch(1,geom%nPnl))

      geom%nPtch = geom%nPnl
      allocate (geom%ptch2Pnl(geom%nPtch))

      do i = 1, geom%nPnl
        geom%pnl2Ptch(1,i) = i
        geom%ptch2Pnl(i) = i
      end do

c     The x, ijk geometry arrays were read in based on panels. Panels
c     are equivalent to patches for P3D files so no reordering is necessary.
c     Convert to an unstructured mesh:

      if (uvType == 'absolute') then

        call struct2Unstruct (geom)

      else !parameterized (u,v) map

c       Create parameterization of patches:

        nPts = geom%ipx(geom%nPnl+1)-1
        allocate (geom%u_ptch(nPts), geom%v_ptch(nPts))

        do i_Ptch = 1, geom%nPtch
          id = geom%ijk(1,i_Ptch)
          jd = geom%ijk(2,i_Ptch)
          ip_s = geom%ipx(i_Ptch)
          ip_e = geom%ipx(i_Ptch+1)-1

          allocate (xtmp(id*jd), ytmp(id*jd), ztmp(id*jd))
          xtmp = geom%x(1,ip_s:ip_e)
          ytmp = geom%x(2,ip_s:ip_e)
          ztmp = geom%x(3,ip_s:ip_e)

          geom%u_ptch(ip_s) = 0. ! -999. suppresses normalization

          call param2d (id, jd, 1, id, 1, jd,
     .                  xtmp, ytmp, ztmp,
     .                  geom%u_ptch(ip_s), geom%v_ptch(ip_s))

          deallocate (xtmp, ytmp, ztmp)
        end do
      end if

      end subroutine prepareP3DGeometry

c+---------------------------------------------------------------------
c
      subroutine prepareUnstructGeometry (geom, iFileGeom)
c
c  Description:  Read unstructured geometry in FAST format.
c  I am assuming that the surface file is single grid. The integer flag
c  for each triangle should specify the component number to which the
c  triangle belongs. The component numbers should run sequentially from
c  1 to nGeomComponent Also, the triangles in the file should be pre-sorted
c  by component.
c
c  Sort the array of triangles into component numbers
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c----------------------------------------------------------------------

c  Modules:

      use TypeDefs
      use cfd_io_package

      implicit none

c  Arguments:

      type (Geometry) ::
     .  geom

      type (ioMeshFileType) ::
     .  iFileGeom

c  Local variables:

      integer ::
     .  i, i_Ptch, l_Ptch, iPrevPtch,
     .  ncall = -1

      type (ioMeshFileType), pointer ::
     .  iFileTmp

      character (len=23) ::
     .  routine = 'prepareUnstructGeometry'

c  Execution:

      allocate (iFileTmp)
      iFileTmp = iFileGeom

c     Read data from file.

      call io_r_unstr_x_fast (geom%x, geom%nPt,
     .                        geom%tri, geom%iTri2Ptch, geom%nTri,
     .                        geom%tet, geom%nTet, ncall, iFileTmp)

c     Determine the total number of components, and create ipTriPnl component
c     pointers into tri array.  Assume the patches are ordered in the file.

      geom%nComp = geom%iTri2Ptch(geom%nTri)
      geom%nPnl  = geom%nComp
      geom%nPtch = geom%nComp

      allocate (geom%ipTriPnl(geom%nPnl+1))

      l_Ptch = 0
      iPrevPtch = -1
      do i = 1, geom%nTri
        i_Ptch = geom%iTri2Ptch(i)
        if (i_Ptch /= iPrevPtch) then
          l_Ptch = l_Ptch + 1
          if (l_Ptch > geom%nPtch) goto 100
          if (l_Ptch /= i_Ptch) goto 101
          geom%ipTriPnl(l_Ptch) = i
          iPrevPtch = geom%iTri2Ptch(i)
        end if
      end do
      geom%ipTriPnl(geom%nPtch+1) = geom%nTri+1

      allocate (geom%iTriUnsrt(geom%nTri))

      do i = 1, geom%nTri
        geom%iTriUnsrt(i) = i
      end do

c     Create component, panel, patch correspondence.  Even though
c     the component/panel/patch are the same entities in an unstructured
c     mesh, filling the arrays allows Aerosurf, P3D, and unstructured
c     databases to be handled exactly the same way.

      allocate (geom%pnlNPtch(geom%nPnl),
     .          geom%pnl2Ptch(1,geom%nPnl),
     .          geom%ptch2Pnl(geom%nPtch))

      do i = 1, geom%nPtch
        geom%pnlNPtch(i)   = 1
        geom%pnl2Ptch(1,i) = i
        geom%ptch2Pnl(i)   = i
      end do

      return

c  Error handling:

 100  call meshwarpGracefulExit (routine,
     .  'Number of components found to be more than iTri2Ptch')

 101  call meshwarpGracefulExit (routine,
     .  'Components in file are not in counting order.')

      end subroutine prepareUnstructGeometry

c+--------------------------------------------------------------------------
c
      subroutine struct2Unstruct (geom)

c
c     Description:  Convert structured surface meshes into unstructured
c     surface meshes.  Return the mesh number in the iTriVerts flag,
c     and return zero tetrahedra.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c---------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs

      implicit none

c  Arguments:

      type (Geometry) ::
     .  geom

c  Local constants:

      integer, parameter :: THREE = 3, FOUR = 4

c  Local variables:

      integer ::
     .  lx, id,jd, nv,
     .  i,j,n, idb,
     .  i1,i2,i3,i4,
     .  i_Tri, i_Ptch, i_Pnl,
     .  l_Blank, l_Ptch

      real ::
     .  diag_1, diag_2,
     .  edge_1, edge_2,
     .  edge_3, edge_4,
     .  eps

      logical, save ::
     .  first = .true.

      type (Geometry), save ::
     .  geom_first

      character (len=15) ::
     .  routine = 'struct2Unstruct'

c  Execution:

c     Copy original triangulation to modified geometry:
      if (.not. first) then
        if (geom%nPtch /= geom_first%nPtch) goto 100
        do i = 1,geom%nPtch
          do j = 1, 2
            if (geom%ijk(j,i) /= geom_first%ijk(j,i)) goto 100
          end do
        end do
        geom%tri => geom_first%tri
        geom%iTriUnsrt => geom_first%iTriUnsrt
        geom%iTri2Ptch => geom_first%iTri2Ptch
        geom%iTri2PtchIJ => geom_first%iTri2PtchIJ
        geom%tet => geom_first%tet
        geom%ipTriPnl => geom_first%ipTriPnl
        geom%nPt = geom_first%nPt
        geom%nTri = geom_first%nTri
        geom%nTet = geom_first%nTet
        return
      end if

      eps = epsilon (eps)

c     Calculate total number of triangles in mesh and set up ip_tri indexing:

      allocate (geom%ipTriPnl(geom%nPnl+1))

      geom%nTri = 0
      do i = 1, geom%nPtch
       j = 2*(geom%ijk(1,i)-1)*(geom%ijk(2,i)-1)
       geom%nTri = geom%nTri + j
      end do

      allocate (geom%tri(THREE,geom%nTri),
     .          geom%iTriUnsrt(geom%nTri))
      geom%iTriUnsrt = -999

      allocate (geom%iTri2Ptch(geom%nTri),
     .          geom%iTri2PtchIJ(2,geom%nTri),
     .          geom%tet(FOUR,0))

      geom%nPt = geom%ipx(geom%nPtch+1)-1
      geom%nTet = 0

c     Triangulate

      i_Tri = 0
      geom%ipTriPnl(1) = 1
      do i_Pnl = 1, geom%nPnl
        do l_Ptch = 1, geom%pnlNPtch(i_Pnl)
          i_Ptch = geom%pnl2Ptch(l_Ptch,i_Pnl)
          id = geom%ijk(1,i_Ptch)
          jd = geom%ijk(2,i_Ptch)
          do i = 1, id-1
          do j = 1, jd-1

c           Establish corner indices of cell
            lx = geom%ipx(i_Ptch)-1
            i1 = lx + (i  )+id*(j-1)  !Lower left corner
            i2 = lx + (i+1)+id*(j-1)
            i3 = lx + (i+1)+id*(j  )
            i4 = lx + (i  )+id*(j  )

c           Establish triangle leg vectors.  Simple triangulation for now
c           Future improvement requires the handling of nonconvex cells with
c           the check_convex subroutine.

            edge_1 = 0.
            edge_2 = 0.
            edge_3 = 0.
            edge_4 = 0.
            diag_1 = 0.
            diag_2 = 0.
            do nv = 1,3
              edge_1 = edge_1 + (geom%x(nv,i2)-geom%x(nv,i1))**2
              edge_2 = edge_2 + (geom%x(nv,i3)-geom%x(nv,i2))**2
              edge_3 = edge_3 + (geom%x(nv,i4)-geom%x(nv,i3))**2
              edge_4 = edge_4 + (geom%x(nv,i1)-geom%x(nv,i4))**2
              diag_1 = diag_1 + (geom%x(nv,i3)-geom%x(nv,i1))**2
              diag_2 = diag_2 + (geom%x(nv,i4)-geom%x(nv,i2))**2
            end do

            if (diag_1 <= diag_2) then
              if (edge_1 > eps .and. edge_2 > eps) then
                i_Tri = i_Tri + 1
                geom%tri(1,i_Tri) = i1
                geom%tri(2,i_Tri) = i2
                geom%tri(3,i_Tri) = i3
                geom%iTriUnsrt(i_Tri) = i_Tri
            l_Blank = geom%iBlank(i1)*geom%iBlank(i2)*geom%iBlank(i3)
                if (l_Blank == 0) then
                  geom%iTri2Ptch(i_Tri) = 0
                  geom%iTri2PtchIJ(1:2,i_Tri) = 0
                else
                  geom%iTri2Ptch(i_Tri) = i_Ptch
                  geom%iTri2PtchIJ(1,i_Tri) = i
                  geom%iTri2PtchIJ(2,i_Tri) = j
                end if
              end if
              if (edge_3 > eps .and. edge_4 > eps) then
                i_Tri = i_Tri + 1
                geom%tri(1,i_Tri) = i3
                geom%tri(2,i_Tri) = i4
                geom%tri(3,i_Tri) = i1
                geom%iTriUnsrt(i_Tri) = i_Tri
            l_Blank = geom%iBlank(i1)*geom%iBlank(i4)*geom%iBlank(i3)
                if (l_Blank == 0) then
                  geom%iTri2Ptch(i_Tri) = 0
                  geom%iTri2PtchIJ(1:2,i_Tri) = 0
                else
                  geom%iTri2Ptch(i_Tri) = i_Ptch
                  geom%iTri2PtchIJ(1,i_Tri) = i
                  geom%iTri2PtchIJ(2,i_Tri) = j
                end if
              end if
            else
              if (edge_2 > eps .and. edge_3 > eps) then
                i_Tri = i_Tri + 1
                geom%tri(1,i_Tri) = i2
                geom%tri(2,i_Tri) = i3
                geom%tri(3,i_Tri) = i4
                geom%iTriUnsrt(i_Tri) = i_Tri
            l_Blank = geom%iBlank(i2)*geom%iBlank(i3)*geom%iBlank(i4)
                if (l_Blank == 0) then
                  geom%iTri2Ptch(i_Tri) = 0
                  geom%iTri2PtchIJ(1:2,i_Tri) = 0
                else
                  geom%iTri2Ptch(i_Tri) = i_Ptch
                  geom%iTri2PtchIJ(1,i_Tri) = i
                  geom%iTri2PtchIJ(2,i_Tri) = j
                end if
              end if
              if (edge_4 > eps .and. edge_1 > eps) then
                i_Tri = i_Tri + 1
                geom%tri(1,i_Tri) = i4
                geom%tri(2,i_Tri) = i1
                geom%tri(3,i_Tri) = i2
                geom%iTriUnsrt(i_Tri) = i_Tri
             l_Blank = geom%iBlank(i1)*geom%iBlank(i2)*geom%iBlank(i4)
                if (l_Blank == 0) then
                  geom%iTri2Ptch(i_Tri) = 0
                  geom%iTri2PtchIJ(1:2,i_Tri) = 0
                else
                  geom%iTri2Ptch(i_Tri) = i_Ptch
                  geom%iTri2PtchIJ(1,i_Tri) = i
                  geom%iTri2PtchIJ(2,i_Tri) = j
                end if
              end if
            end if
          end do
          end do
        end do
        geom%ipTriPnl(i_Pnl+1) = i_Tri+1
      end do
      geom%nTri = i_Tri

      if (first) then
        geom_first%ijk => geom%ijk
        geom_first%tri => geom%tri
        geom_first%iTriUnsrt => geom%iTriUnsrt
        geom_first%iTri2Ptch => geom%iTri2Ptch
        geom_first%iTri2PtchIJ => geom%iTri2PtchIJ
        geom_first%tet => geom%tet
        geom_first%ipTriPnl => geom%ipTriPnl
        geom_first%nPt = geom%nPt
        geom_first%nTri = geom%nTri
        geom_first%nTet = geom%nTet
        geom_first%nPtch = geom%nPtch
        first = .false.
      end if

      return

c  Error handling:

 100  call meshwarpGracefulExit (routine,
     .  'Absolute mapping requires fixed point count in geo-patches.')

      end subroutine struct2Unstruct

c----------------------------------------------------------------------------
c
c     Mesh, connectivity, uv input/output routines:
c
c+---------------------------------------------------------------------------
c
      subroutine readMeshData (mesh, iFile, meshType)
c
c  Description:  Read the volume/surface mesh from a file.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c----------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs
      use cfd_io_package

      implicit none

c  Arguments:

      type (AbstractMeshType) ::
     .  mesh

      type (ioMeshFileType) ::
     .  iFile

      character (*) ::
     .  meshType

c  Local variables:

      integer :: ncall = -1

      type (ioMeshFileType), pointer ::
     .  iFileTmp

      character (len=12) ::
     .  routine = 'readMeshData'

c  Execution:

      write (StdOut,*) 'Reading mesh...'
      call flush (StdOut)

      allocate (iFileTmp)
      iFileTmp = iFile

c     Read mesh from file in correct format:

      select case (meshType)
        case ('strct')
          if (iFile%iblank) then
            call io_r_str_x_p3d_ib (mesh%x, mesh%nBlock, mesh%ijk,
     .                              mesh%ipx, mesh%iBlank, ncall,
     .                              iFileTmp)
          else
            call io_r_str_x_p3d (mesh%x, mesh%nBlock, mesh%ijk,
     .                           mesh%ipx, ncall, iFileTmp)
            if (associated (mesh%iBlank)) deallocate (mesh%iBlank)
            allocate (mesh%iBlank(mesh%ipx(mesh%nBlock+1)-1))
            mesh%iBlank = 1
          end if
          mesh%nPt = mesh%ipx(mesh%nBlock+1)-1
        case ('unstrct')
          call io_r_unstr_x_fast (mesh%x, mesh%nPt,
     .                            mesh%tri, mesh%iTri, mesh%nTri,
     .                            mesh%tet, mesh%nTet, ncall, iFileTmp)
            if (associated (mesh%iBlank)) deallocate (mesh%iBlank)
            allocate (mesh%iBlank(mesh%nPt))
            mesh%iBlank = 1
            mesh%nBlock = 1
        case default
          call meshwarpGracefulExit (routine,
     .      'Unknown type of mesh file.')
      end select

 999  return

      end subroutine readMeshData

c+-----------------------------------------------------------------------
c
      subroutine readConnFile(bc, df, hf, mesh, iFile, connType)
c
c  Description:  Read bc connectivity from file.  Only syn107-mb
c  connectivity files are supported at this time.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs
      use cfd_io_package

      implicit none

c  Arguments:

      type (BCData) :: bc
      type (DesignFace) :: df
      type (HyperFace)  :: hf
      type (AbstractMeshType) :: mesh
      type (ioMeshFileType) :: iFile
      character (*) :: connType

c  Local variables:

      real, pointer ::
     .  rhand(:)

      integer ::
     .  ncall = -1

      type (ioMeshFileType), pointer ::
     .  iFileTmp

      character (len=12) ::
     .  routine = 'readConnFile'

c  Execution:

      write (StdOut,*) 'Reading connectivity file...'
      call flush (StdOut)

      allocate (iFileTmp)
      iFileTmp = iFile

      select case (connType)
        case ('syn107')
          call io_r_bc_flo107 (bc%typ, bc%algn, rhand,
     .                         bc%nBC, ncall, iFileTmp)

          call io_r_hyp_syn107 (hf%ihyp, hf%nHF, iFileTmp)

          call io_r_dsnf_syn107 (df%nDF, df%srcBlock, df%srcFace,
     .                           df%trgBlock, df%trgFace, iFileTmp)

        case default
          call meshwarpGracefulExit (routine,
     .      'Unknown type of connectivity file.')
      end select

      if (bc%nBC /= mesh%nBlock)
     .  call meshwarpGracefulExit (routine,
     .    'Number of bcs does not match number of grid blocks.')

      end subroutine readConnFile

c+-----------------------------------------------------------------------
c
      subroutine readUVFile (uv, nHF, iFile, uvType)
c
c     Description:  Read the uv data from a file. Need to implement
c     something for the AEROSURF subpatching method.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs
      use cfd_io_package

      implicit none

c  Arguments:

      integer ::
     .  nHF

      type(UVMap) ::
     .  uv

      type(ioMeshFileType) ::
     .  iFile

      character (*) ::
     .  uvType

c  Local variables:

      integer ::
     .  ncall = -1

      type (ioMeshFileType), pointer ::
     .  iFileTmp

      character (len=9) ::
     .  routine = 'readUVMap'

c  Execution:

      write (StdOut,*) 'Reading (u,v) map file...'
      call flush (StdOut)

      allocate (iFileTmp)
      iFileTmp = iFile

      call io_r_str_x_p3d_ib (uv%uv, uv%nHyp, uv%ijk, uv%ip,
     .                        uv%iuv, ncall, iFileTmp)

      if (uv%nHyp /= nHF) goto 100

      return

c  Error handling:

 100  call meshwarpGracefulExit (routine,
     .  'Discrepancy in no. of hyperfaces in conn and uv file.')

      end subroutine readUVFile

c+-----------------------------------------------------------------------
c
      subroutine writeMeshData (mesh, iFile, meshType)
c
c  Description:  Write a surface/volume mesh to a file.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs
      use cfd_io_package

      implicit none

c  Arguments:

      type(AbstractMeshType) ::
     .  mesh

      type(ioMeshFileType) ::
     .  iFile

      character (*) ::
     .  meshType

c  Local variables:

      integer :: ncall = -1

      type (ioMeshFileType), pointer :: iFileTmp

      character (len=12) ::
     .  routine = 'writeMeshData'

c  Execution:

      write (StdOut,'(a)') ' Writing warped mesh to file...'
      call flush (StdOut)

      allocate (iFileTmp)
      iFileTmp = iFile

c     Read mesh from file in correct format:

      select case (meshType)
        case ('strct')
          if (iFile%iblank) then
            call io_w_str_x_p3d_ib (mesh%x, mesh%nBlock, mesh%ijk,
     .                              mesh%ipx, mesh%iBlank, ncall,
     .                              iFileTmp)
          else
            call io_w_str_x_p3d (mesh%x, mesh%nBlock, mesh%ijk,
     .                           mesh%ipx, ncall, iFileTmp)
          end if
        case ('unstrct')
          call io_w_unstr_x_fast (mesh%x, mesh%nPt,
     .                            mesh%tri, mesh%iTri, mesh%nTri,
     .                            mesh%tet, mesh%nTet, ncall, iFileTmp)
        case default
          call meshwarpGracefulExit (routine, 'Unknown mesh file type.')
      end select

 999  return

      end subroutine writeMeshData

c+-----------------------------------------------------------------------------
c
      subroutine duplicateMesh (meshSrc, meshTrg, meshType)
c
c  Description:  Copy the x array of a source mesh into a target mesh.
c  Since the x array is the only data to change in the warping procedure,
c  link the data arrays directly to the source.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c------------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs

      implicit none

c  Arguments:

      type (AbstractMeshType) ::
     .  meshSrc, meshTrg

      character (*) :: meshType

c  Execution:

      if (meshType == 'strct') then
        allocate (meshTrg%x(3,meshSrc%nPt))
        meshTrg%x      =  meshSrc%x
        meshTrg%iBlank => meshSrc%iBlank
        meshTrg%ipx    => meshSrc%ipx
        meshTrg%ijk    => meshSrc%ijk
        meshTrg%nBlock =  meshSrc%nBlock
        meshTrg%nPt    =  meshSrc%nPt
      else
        allocate (meshTrg%x(3,meshSrc%nPt))
        meshTrg%x      =  meshSrc%x
        meshTrg%iBlank => meshSrc%iBlank
        meshTrg%tri    => meshSrc%tri
        meshTrg%iTri   => meshSrc%iTri
        meshTrg%tet    => meshSrc%tet
        meshTrg%nPt    =  meshSrc%nPt
        meshTrg%nTri   =  meshSrc%nTri
        meshTrg%nTet   =  meshSrc%nTet
      end if

      end subroutine duplicateMesh

c----------------------------------------------------------------------------
c
c     Master Edge/Corner Creation Routines
c
c
c+-----------------------------------------------------------------------
c
      subroutine activateMasterInfo (lc, nBlock, maxDim,
     .                               bcType, ir, df, hf,
     .                               msect, msecnct, mscct, msccnct,
     .                               ige2mse, igc2msc, mse2igb, mse2ige,
     .                               msc2igb, msc2ige, mse2ige_actHF,
     .                               mse2ige_actDF, msc2igc_actHF,
     .                               msc2igc_actDF, ige2mse_actHF,
     .                               ige2mse_actDF,
     .                               igb_act, igf_act, ige_act)
c
c  Description:  Flag the active edges and faces in the master list
c  corresponding to the hyper and design faces.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs

      implicit none

c  Arguments:

      integer, intent (in) ::
     .  nBlock,
     .  maxDim,
     .  msect,
     .  mscct

      integer, intent (in) ::
     .  bcType(MBFACE,nBlock),
     .  ir(MBCRD,MBCRD,MBFACE,nBlock),
     .  ige2mse(MBEDGE,nBlock),
     .  igc2msc(MBCRN,nBlock),
     .  msecnct(MBEDGE*nBlock),
     .  mse2igb(MAXECNCT,MBEDGE*nBlock),
     .  mse2ige(MAXECNCT,MBEDGE*nBlock),
     .  msccnct(MBCRN*nBlock),
     .  msc2igb(MAXCRNCNCT,MBEDGE*nBlock),
     .  msc2ige(MAXCRNCNCT,MBEDGE*nBlock)

      type (UnitBlockNumbering) ::
     .  lc

      type (DesignFace), intent (inout) ::
     .  df

      type (HyperFace),  intent (in) ::
     .  hf

      integer :: !intent (out)
     .  mse2ige_actHF(msect),
     .  mse2ige_actDF(msect),
     .  msc2igc_actHF(mscct),
     .  msc2igc_actDF(mscct),
     .  ige2mse_actHF(MBEDGE,nBlock),
     .  ige2mse_actDF(MBEDGE,nBlock),
     .  igb_act(nBlock),
     .  igf_act(MBFACE,nBlock),
     .  ige_act(MBEDGE,nBlock)

c  Local variables:

      integer ::
     .  i, j, k, l, m, n,
     .  i_b, i_f, i_e, i_c, j_e,
     .  mse, msc, i_eAdj,
     .  nb_b, nb_f, ltot

      integer, allocatable ::
     .  nbfacid(:,:)

c  Execution:

      write (StdOut,'(4x,a)') 'Activating master edges/corners...'
      call flush (StdOut)

c     Initialize activation lists:
      mse2ige_actHF = 0
      mse2ige_actDF = 0
      msc2igc_actHF = 0
      msc2igc_actDF = 0
      ige2mse_actHF = 0
      ige2mse_actDF = 0
      igb_act = 0
      igf_act = 0
      ige_act = 0

      df%c_type = 0
      df%e_type = 0

c     Initialize local variables:

      allocate (nbfacid(MBFACE,nBlock))
      nbfacid = 0

c     Establish the face id's for the neighboring mating faces on all blocks.
c     (I'm surprised this info hasn't been used before.  Check mstedg, mstcrn
c     for possible re-use.)

      do i_b = 1, nBlock
        do i_f = 1, MBFACE
          if (bcType(i_f,i_b) > 0) then
            if (abs(ir(1,1,i_f,i_b)) > 0) then
              if (ir(1,1,i_f,i_b) > 0) then
                nbfacid(i_f,i_b) = 2
              else
                nbfacid(i_f,i_b) = 1
              end if
            else if (abs(ir(1,2,i_f,i_b)) > 0) then
              if (ir(1,2,i_f,i_b) > 0) then
                nbfacid(i_f,i_b) = 4
              else
                nbfacid(i_f,i_b) = 3
              end if
            else if (abs(ir(1,3,i_f,i_b)) > 0) then
              if (ir(1,3,i_f,i_b) > 0) then
                nbfacid(i_f,i_b) = 6
              else
                nbfacid(i_f,i_b) = 5
              end if
            end if
          end if
        end do
      end do

c     For each hyperface, flag active master and global edges, master
c     corners, and set the type of active faces to type = 2.

      do i = 1, hf%nHF
        i_b = hf%ihyp(1,i)  !HypF block number
        i_f = hf%ihyp(2,i)  !HypF face number
        igf_act(i_f,i_b) = 2

        do j = 1, MLEDGE

c         Activate the four edges of the face.

          i_e = lc%f2ei(j,i_f)
          ige2mse_actHF(i_e,i_b) = 1
          mse = ige2mse(i_e,i_b)
          mse2ige_actHF(mse) = 1

c         Activate the two corners of this edge.

          i_c = lc%e2ci(1,i_e)
          msc = igc2msc(i_c,i_b)
          msc2igc_actHF(msc) = 1

          i_c = lc%e2ci(2,i_e)
          msc = igc2msc(i_c,i_b)
          msc2igc_actHF(msc) = 1

        end do
      end do

c     For each design face and their mating
c     neighbors, flag edge and face activity type=2.
c
      do i = 1, df%nDF
        i_b = df%trgBlock(i)
        i_f = df%trgFace(i)
        igf_act(i_f,i_b) = 2

c       Turn on the neighboring face:
        if (bcType(i_f,i_b) > 0) then
          nb_b = bcType(i_f,i_b)
          nb_f = nbfacid(i_f,i_b)
          igf_act(nb_f,nb_b) = 2
        end if

c       Modify the corner and face activation:
        do l = 1, MLCRN !ok
          i_c = lc%f2ci(l,i_f)    !Corner mods
          msc = igc2msc(i_c,i_b)
          msc2igc_actDF(msc) = 1
          if (msc2igc_actHF(msc) > 0)
     .      df%c_type(l,i) = 1

          i_e = lc%f2ei(l,i_f)    !Edge mods
          mse = ige2mse(i_e,i_b)
          ige2mse_actDF(i_e,i_b) = 1
          mse2ige_actDF(mse) = 1
          if (mse2ige_actHF(mse) > 0)
     .      df%e_type(l,i) = 1
        end do

      end do

c     Set all active design and hyper face edges
c     to type=2 by looping through the
c     master edges and detecting activity.

      do mse = 1, msect
        if (mse2ige_actHF(mse) == 1 .or.
     .      mse2ige_actDF(mse) == 1     ) then
          do i = 1, msecnct(mse)
            i_b = mse2igb(i,mse)
            i_e = mse2ige(i,mse)
            ige_act(i_e,i_b) = 2
          end do
        end if
      end do

c     Set all edges connected to an active master
c     corner or active master design corner to
c     be at least implicit edges; type=1.

      do msc = 1, mscct
        if (msc2igc_actHF(msc) == 1 .or.
     .      msc2igc_actDF(msc) == 1     ) then
          do i = 1, msccnct(msc)
            i_b = msc2igb(i,msc)
            i_e = msc2ige(i,msc)
            ige_act(i_e,i_b) = max (ige_act(i_e,i_b), 1)
          end do
        end if
      end do

c     For each block....

      do i_b = 1, nBlock

c       Loop through edges on the block and for those that are active, ensure
c       that their locally connected edges are set to at least type=1.

        do i_e = 1, MBEDGE
          if (ige_act(i_e,i_b) == 2) then
            do l = 1, 2
              i_f = lc%e2fi(l,i_e)
              j_e = lc%e2fj(l,i_e)
              if (j_e <= 2) then
                i_eAdj = lc%f2ei(3,i_f)
                ige_act(i_eAdj,i_b) = max (1, ige_act(i_eAdj,i_b))
                i_eAdj = lc%f2ei(4,i_f)
                ige_act(i_eAdj,i_b) = max (1, ige_act(i_eAdj,i_b))
              else
                i_eAdj = lc%f2ei(1,i_f)
                ige_act(i_eAdj,i_b) = max (1, ige_act(i_eAdj,i_b))
                i_eAdj = lc%f2ei(2,i_f)
                ige_act(i_eAdj,i_b) = max (1, ige_act(i_eAdj,i_b))
              end if
            end do
          end if
        end do

c       Loop through all faces on all blocks and test
c       whether any edges are to be moved on non-active faces.
c       For those faces that contain such edges, set the face type = 1.
c       If the block contains any active faces, set the block activity
c       to at least 1.  (I'm not sure at the moment how igb_act is supposed
c       to obtain values greater than 1!)
c
        do i_f = 1, MBFACE
          if (igf_act(i_f,i_b) == 0) then
            ltot = 0
            do l = 1,4
              i_e = lc%f2ei(l,i_f)
              ltot = ltot + ige_act(i_e,i_b)
            end do

            if (ltot > 0) igf_act(i_f,i_b) = 1

          end if
          if (igf_act(i_f,i_b) > 0) then
            igb_act(i_b) = max (1, igb_act(i_b))  !This seems to be the only
          end if                                  !place to set iblkptg so why
        end do                                    !is the max necessary?

      end do !Next block

      end subroutine activateMasterInfo

c+-----------------------------------------------------------------------
c
      subroutine createDesignFacePtr (df, lc,
     .                                nBlock, ijkx, ipx,
     .                                bcType, bcAlign, is, ir,
     .                                mscct, msect,
     .                                igc2msc, ige2mse, ige2mse_dir)
c
c     Description:  Establish pointers from design faces into the x
c     array of the mesh.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs

      implicit none

c  Arguments:

      type (DesignFace), intent (out) ::
     .  df

      type (UnitBlockNumbering), intent (in) ::
     .  lc

      integer, intent (in) ::
     .  nBlock,
     .  ijkx(3,nBlock),
     .  ipx(nBlock)

      integer, intent (in) ::
     .  bcType(MBFACE,nBlock),
     .  bcAlign(MBCRD,MBFACE,nBlock),
     .  is(MBCRD,MBCRD,MBFACE),
     .  ir(MBCRD,MBCRD,MBFACE,nBlock)

      integer, intent (in) ::
     .  mscct,
     .  msect,
     .  igc2msc(MBCRN,nBlock),
     .  ige2mse(MBEDGE,nBlock),
     .  ige2mse_dir(MBEDGE,nBlock)

c  Local variables:

      integer ::
     .  i, i1, i2, i_c, i_df, ip_s,
     .  p, ps, pe, pd,
     .  q, qs, qe, qd,
     .  r, rs, re,
     .  df_b, df_f, df_c,
     .  iPtr, lDFpt, nDFpt,
     .  a_c, ac_cnct,
     .  l_c, msc, ipDFE, i_dfe,
     .  ip_dfx,
     .  npAe, l, l_e, df_e,
     .  mse, mse_dir, a_e,
     .  i_s, i_e

      integer ::
     .  pqr(3)

c  Execution:

      write (StdOut,'(4x,a)') 'Assigning design face pointers...'
      call flush (StdOut)

      allocate (df%ijk(MBCRD,df%nDF))
      df%ijk = 0
      allocate (df%ip(df%nDF+1))
      df%ip = 0

c     Count the number of design face points, note that this assumes that design
c     face targets and sources are from the same block and opposite faces.
c
      df%ip(1) = 1
      do i_df = 1, df%nDF
        df_b = df%srcBlock(i_df)
        df_f = df%srcFace(i_df)
        pqr = matmul(is(1:MBCRD,1:MBCRD,df_f),ijkx(1:MBCRD,df_b))
        lDFpt = pqr(2)*pqr(3)
        df%ijk(1:MBCRD,i_df) = pqr(1:MBCRD)
        df%ip(i_df+1) = df%ip(i_df) + lDFpt
      end do

      nDFpt = df%ip(df%nDF+1)-1
      allocate (df%dfsrc2xp(nDFpt))
      df%dfsrc2xp = 0

      allocate (df%dftrg2xp(nDFpt))
      df%dftrg2xp = 0

      allocate (df%dftrg2nbxp(nDFpt))
      df%dftrg2nbxp = 0

c     Loop over source design faces and store the pointers to x array

      do i_df = 1, df%nDF
        df_b = df%srcBlock(i_df)
        df_f = df%srcFace(i_df)

c       Source design face-2-mesh face pointers

        call facxv (is, df_f, df%ip(i_df), ipx(df_b), df%ijk(1,i_df),
     .              ijkx(1,df_b), df%dfsrc2xp, nDFpt)

c       Load target design face2block pointers and
c       design face2neighbor block pointers.

        df_b = df%trgBlock(i_df)
        df_f = df%trgFace(i_df)

        call facxd (bcType, bcAlign, is, ir, nBlock, ipx, ijkx,
     .              df_f, df_b, df%ip(i_df), df%ijk(1,i_df), nDFpt,
     .              df%dftrg2xp, df%dftrg2nbxp)
      end do

c     For each corner of the design face, determine its location in
c     the design-face-x array:

      allocate (df%dfc2dfxp(MLCRN,df%nDF))
      df%dfc2dfxp = 0

      do i_df = 1, df%nDF
        df_b = df%trgBlock(i_df)
        df_f = df%trgFace(i_df)

        ip_s = df%ip(i_df)

        pd = df%ijk(2,i_df)
        qd = df%ijk(3,i_df)

        i_c = 1
        do q = 1, qd, qd-1
        do p = 1, pd, pd-1
          df%dfc2dfxp(i_c,i_df) = ip_s + (p-1) + pd*(q-1)
          i_c = i_c + 1
        end do
        end do
      end do

c     Establish a unique list of averaging corners and all design face
c     corners which are attached to the averaging corners.

      allocate (df%msc2ac(mscct))
      df%msc2ac = 0
      allocate (df%ac2igb(MAXCRNCNCT,MLCRN*df%nDF))  !ok Overdimensioned
      df%ac2igb = 0
      allocate (df%accnct(MLCRN*df%nDF)) !ok
      df%accnct = 0
      allocate (df%dfc2accnct(MLCRN,df%nDF)) !ok
      df%dfc2accnct = 0
      allocate (df%dfc2ac(MLCRN,df%nDF)) !ok
      df%dfc2ac = 0

      df%acct = 0

c     For each design face:
      do i_df = 1, df%nDF
        df_b = df%trgBlock(i_df)
        df_f = df%trgFace(i_df)

c       For each design face corner:
        do l_c = 1, MLCRN !ok

          df_c = lc%f2ci(l_c,df_f)  !unit block corner of this local corner
          msc  = igc2msc(df_c,df_b) !master corner this block corner
          a_c  = df%msc2ac(msc)     ! is attached to.

          if (a_c == 0) then !This master corner has not been associated with
                             ! an averaging corner.
            df%acct = df%acct + 1
            df%accnct(df%acct) = 1
            df%ac2igb(1,df%acct) = df_b
            df%msc2ac(msc) = df%acct

            df%dfc2accnct(l_c,i_df) = 1
            df%dfc2ac(l_c,i_df) = df%acct

          else  !An averaging corner has been created for this master corner

c           Check to see if this block corner has already been inserted
c           into this averaging corner list.
            ac_cnct = df%accnct(a_c)

            do i = 1, ac_cnct
              if (df_b == df%ac2igb(i,a_c)) then
                df%dfc2accnct(l_c,i_df) = i
                df%dfc2ac(l_c,i_df) = a_c
                goto 10
              end if
            end do

            ac_cnct = ac_cnct + 1
            df%accnct(a_c) = ac_cnct
            df%ac2igb(ac_cnct,a_c) = df_b

            df%dfc2accnct(l_c,i_df) = ac_cnct
            df%dfc2ac(l_c,i_df) = a_c

 10       end if

        end do
      end do

c     Establish averaging edge to design face edge connections,
c     count the total number of design face edge points:

      df%nDFEpts = 0
      do i_df = 1, df%nDf
        df%nDFEpts = df%nDFEpts +
     .               2*df%ijk(2,i_df) + 2*df%ijk(3,i_df)
      end do

      allocate (df%dfe2dfxp(df%nDFEpts))
      df%dfe2dfxp = 0

      allocate (df%dfeip(MLEDGE,df%nDF))
      df%dfeip = 0

      allocate (df%dfe2ae(MLEDGE,df%nDF))
      df%dfe2ae = 0

      allocate (df%dfe2aecnct(MLEDGE,df%nDF))
      df%dfe2aecnct = 0

      allocate (df%dfe2aexp(df%nDFEpts))
      df%dfe2aexp = 0

      allocate (df%aecnct(MLEDGE*df%nDF))
      df%aecnct = 0

      allocate (df%aeip(MLEDGE*df%nDF)) !Overdimensioned
      df%aeip = 0

      allocate (df%mse2ae(msect))
      df%mse2ae = 0

      allocate (df%ae2igb(MAXECNCT,MLEDGE*df%nDF)) !Overdimensioned
      df%ae2igb = 0

      allocate (df%dfenPt(MLEDGE,df%nDF))
      df%dfenPt = 0

      df%aect = 0

c     Construct pointers from the design face edges to their location
c     in the design face x array.

      ipDFE = 1
      do i_df = 1, df%nDF
        pd = df%ijk(2,i_df)
        qd = df%ijk(3,i_df)
        ip_dfx = df%ip(i_df)
        i_dfe = 1

        do q = 1,qd,qd-1
          df%dfeip(i_dfe,i_df) = ipDFE
          df%dfenPt(i_dfe,i_df) = pd
          do p = 1,pd,1
            df%dfe2dfxp(ipDFE) = ip_dfx + (p-1) + pd*(q-1)
            ipDFE = ipDFE + 1
          end do
          i_dfe = i_dfe + 1
        end do

        do p = 1,pd,pd-1
          df%dfeip(i_dfe,i_df) = ipDFE
          df%dfenPt(i_dfe,i_df) = qd
          do q = 1,qd,1
            df%dfe2dfxp(ipDFE) = ip_dfx + (p-1) + pd*(q-1)
            ipDFE = ipDFE + 1
          end do
          i_dfe = i_dfe + 1
        end do
      end do

      npAe = 1
      do i_df = 1, df%nDF
        df_b = df%trgBlock(i_df)
        df_f = df%trgFace(i_df)

        do l_e = 1, MLEDGE
          df_e = lc%f2ei(l_e,df_f)
          mse = ige2mse(df_e,df_b)
          mse_dir = ige2mse_dir(df_e,df_b)
          a_e = df%mse2ae(mse)

          if (a_e == 0) then !The master edge has not been associated with
                             !an averaging edge.
            df%aect = df%aect + 1
            df%aecnct(df%aect) = 1
            df%aeip(df%aect) = npAe
            df%ae2igb(1,df%aect) = df_b
            npAe = npAe + df%dfenPt(l_e,i_df)
            df%aeip(df%aect+1) = npAe

            df%mse2ae(mse) = df%aect
            df%dfe2ae(l_e,i_df) = df%aect
            df%dfe2aecnct(l_e,i_df) = 1

            i1 = df%aeip(df%aect)
            i2 = df%aeip(df%aect+1) - 1

            i_s = i1*(mse_dir+1)/2 + i2*(1-mse_dir)/2
            i_e = i1*(1-mse_dir)/2 + i2*(mse_dir+1)/2

            l = df%dfeip(l_e,i_df)
            do i = i_s, i_e, mse_dir
              df%dfe2aexp(l) = i
              l = l+1
            end do

          else

            do i = 1,df%aecnct(a_e)
              if (df_b == df%ae2igb(i,a_e)) then
                df%dfe2ae(l_e,i_df) = a_e
                df%dfe2aecnct(l_e,i_df) = i
                goto 20
              end if
            end do

            df%aecnct(a_e) = df%aecnct(a_e) + 1
            df%ae2igb(df%aecnct(a_e),a_e) = df_b
            df%dfe2ae(l_e,i_df) = a_e
            df%dfe2aecnct(l_e,i_df) = df%aecnct(a_e)

 20         i1 = df%aeip(a_e)
            i2 = df%aeip(a_e+1) - 1

            i_s = i1*(mse_dir+1)/2 + i2*(1-mse_dir)/2
            i_e = i1*(1-mse_dir)/2 + i2*(mse_dir+1)/2

            l = df%dfeip(l_e,i_df)
            do i = i_s, i_e, mse_dir
              df%dfe2aexp(l) = i
              l = l+1
            end do

          end if

        end do
      end do

      end subroutine createDesignFacePtr

c+-----------------------------------------------------------------------
c
      subroutine createHyperFacePtr (hf, is, nBlock, ijk, ipx)
c
c     Description:  Establish pointers from hyperaces into the x array.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs

      implicit none

c  Arguments:

      type (HyperFace) ::
     .  hf

      integer, intent (in) ::
     .  nBlock

      integer, intent (in) ::
     .  ijk(MBCRD,nBlock),
     .  ipx(nBlock+1),
     .  is(MBCRD,MBCRD,MBFACE)

c  Local variables:

      integer ::
     .  p, ps, pe,
     .  q, qs, qe,
     .  r, rs, re,
     .  hf_b, hf_f,
     .  iPtr, nHFpt,
     .  i_hf, ipx_s

      integer ::
     .  pqr(MBCRD,hf%nHF)

c  Execution:

      write (StdOut,'(4x,a)') 'Assigning hyperface pointers...'
      call flush (StdOut)

c     Count the number of hyper face points
      nHFpt = 0
      do i_hf = 1, hf%nHF
        hf_b = hf%ihyp(1,i_hf)
        hf_f = hf%ihyp(2,i_hf)
        pqr(1:MBCRD,i_hf) =
     .    matmul(is(1:MBCRD,1:MBCRD,hf_f),ijk(1:MBCRD,hf_b))
        nHFpt = nHFpt + pqr(2,i_hf)*pqr(3,i_hf)
      end do

      allocate (hf%hf2xp(nHFpt))
      hf%hf2xp = 0

      allocate (hf%ip(hf%nHF+1))
      hf%ip = 0

c     Loop over hyperfaces and store the pointers to x array

      iPtr = 1
      do i_hf = 1, hf%nHF
        hf_b = hf%ihyp(1,i_hf)
        hf_f = hf%ihyp(2,i_hf)

c       Set up indexing of hyperface

        ps = mod(hf_f-1,2)*(pqr(1,i_hf)-1) + 1
        qs = 1
        qe = pqr(2,i_hf)
        rs = 1
        re = pqr(3,i_hf)

        ipx_s = ipx(hf_b)

c       Load the pointers into hf2xp

        hf%ip(i_hf) = iPtr
        do r = rs, re
        do q = qs, qe
          hf%hf2xp(iPtr) = (ipx_s) +   (   ps*is(1,1,hf_f)
     .                                   + q *is(2,1,hf_f)
     .                                   + r *is(3,1,hf_f)-1)
     .       + ijk(1,hf_b)*            (   ps*is(1,2,hf_f)
     .                                   + q *is(2,2,hf_f)
     .                                   + r *is(3,2,hf_f)-1)
     .       + ijk(1,hf_b)*ijk(2,hf_b)*(   ps*is(1,3,hf_f)
     .                                   + q *is(2,3,hf_f)
     .                                   + r *is(3,3,hf_f)-1)
          iPtr = iPtr + 1
        end do
        end do

      end do

      hf%ip(hf%nHF+1) = iPtr

      end subroutine createHyperFacePtr

c+-----------------------------------------------------------------------
c
      subroutine createMasterInfoUnstr (hf, mesh, uv)
c
c  Description:  Create hyperface point list for an unstructured mesh.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs

      implicit none

c  Arguments:

      type (HyperFace), intent (inout) ::
     .  hf

      type (AbstractMeshType), intent (inout) ::
     .  mesh

      type (UVMap), intent (in) ::
     .  uv

c  Local variables:

      integer ::
     .  i, j, l, i_HF,
     .  nHFTri, nHFPt, nHF,
     .  hfPrev

      character (len=21) ::
     .  routine = 'createMasterInfoUnstr'

c  Execution:

c     Count the number of triangles and points in the hyperfaces:

      mesh%iBlank = 0

      nHFTri = 0
      nHF = 0
      do i = 1, mesh%nTri
        i_HF = mesh%iTri(i)
        nHF = max (nHF, i_HF)
        if (i_HF /= 0) then
          nHFTri = nHFTri + 1
          do j = 1, 3
            mesh%iBlank(mesh%tri(j,i)) = -i_HF  !Mark all x in surface tri's
          end do
        end if
      end do

      nHFPt = 0
      do i = 1, mesh%nPt
        if (mesh%iBlank(i) /= 0) nHFPT = nHFPt + 1
      end do

      allocate (hf%ip(nHF+1), hf%hf2xp(nHFPt))

      hfPrev = 0
      i_HF = 0
      l = 0
      do i = 1, mesh%nPt
        if (mesh%iBlank(i) /= 0) then
          l = l+1
          if (hfPrev /= mesh%iBlank(i)) then
            hfPrev = mesh%iBlank(i)
            i_HF = i_HF + 1
            hf%ip(i_HF)  = l
          end if
          hf%hf2xp(l) = i
        end if
      end do

      hf%ip(nHF+1) = l+1

c     Check for consistency with the (u,v) map:
      if (uv%nHyp /= hf%nHF) goto 100
      if (hf%ip(hf%nHF+1) /= uv%ip(uv%nHyp+1)) goto 101

      return

c  Error handling:

 100  call meshwarpGracefulExit (routine,
     .  'Number of hyperfaces does not match (u,v) map')
 101  call meshwarpGracefulExit (routine,
     .  'Number of hyperface points does not match (u,v) map')

      end subroutine createMasterInfoUnstr

c+-----------------------------------------------------------------------
c
      subroutine createMasterInfo (mstr, hf, df, mesh, bc)
c
c  Description:  Create master face/edge/corner information.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs

      implicit none

c  Arguments:

      type (AbstractMeshType), intent (in) ::
     .  mesh

      type (BCData), intent (inout) ::
     .  bc

      type (MasterData), intent (out) ::
     .  mstr

      type (HyperFace), intent (out) ::
     .  hf

      type (DesignFace), intent (out) ::
     .  df

c  Local variables:

      integer ::
     .  nBlock,
     .  maxDim

      type (UnitBlockNumbering) ::
     .  lc

c  Execution:

      write (StdOut,'(a)') ' Creating topological interfaces...'
      call flush (StdOut)

      nBlock = mesh%nBlock
      maxDim = maxval (mesh%ijk)

c     Establish the unit-block face/edge/corner numbering convention:
      call setLocalNumConvention (lc)

c     Allocate transformation matrices:

      allocate (bc%is(MBCRD,MBCRD,MBFACE),
     .          bc%it(MBCRD,MBCRD,MBFACE,nBlock),
     .          bc%ir(MBCRD,MBCRD,MBFACE,nBlock))

c     Create transformation matrices is, it, ir:

      call createTransformMat (bc%is, bc%it, bc%ir,
     .                         bc%typ, bc%algn, nBlock)

c     Allocate master edge arrays:

      allocate (mstr%msep(MBEDGE*nBlock),
     .          mstr%msecnct(MBEDGE*nBlock),
     .          mstr%mse2igep(MAXECNCT,MBEDGE*nBlock*maxDim),
     .          mstr%mse2igb(MAXECNCT,MBEDGE*nBlock),
     .          mstr%mse2ige(MAXECNCT,MBEDGE*nBlock),
     .          mstr%mse2ige_dir(MAXECNCT,MBEDGE*nBlock),
     .          mstr%ige2mse(MBEDGE,nBlock),
     .          mstr%ige2mse_dir(MBEDGE,nBlock))

c     Establish master edge arrays and pointers:

      call createMasterEdges (lc, nBlock, mesh%ijk,
     .                        mesh%ipx, maxDim,
     .                        bc%typ, bc%algn, bc%is, bc%ir,
     .                        mstr%msep, mstr%msecnct, mstr%mse2igep,
     .                        mstr%mse2igb, mstr%mse2ige,
     .                        mstr%mse2ige_dir, mstr%ige2mse,
     .                        mstr%ige2mse_dir, mstr%msect,
     .                        mstr%mseptct)

c     Allocate master corner arrays:

      allocate (mstr%mse2msc(2,mstr%msect),
     .          mstr%msccnct(MBCRN*nBlock),
     .          mstr%msc2igb(MAXCRNCNCT,MBCRN*nBlock),
     .          mstr%msc2ige(MAXCRNCNCT,MBCRN*nBlock),
     .          mstr%msc2igc(MAXCRNCNCT,MBCRN*nBlock),
     .          mstr%msc2igcp(MAXCRNCNCT,MBCRN*nBlock),
     .          mstr%igc2msc(MBCRN,nBlock))

c     Establish master corner arrays and pointers:

      call createMasterCorners (lc, nBlock, maxDim,
     .                          mstr%msect, mstr%msecnct, mstr%msep,
     .                          mstr%mse2igb, mstr%mse2ige,
     .                          mstr%mse2igep, mstr%mse2ige_dir,
     .                          mstr%mse2msc, mstr%ige2mse,
     .                          mstr%ige2mse_dir, mstr%mscct,
     .                          mstr%msccnct, mstr%msc2igb,
     .                          mstr%msc2ige, mstr%msc2igc,
     .                          mstr%msc2igcp, mstr%igc2msc)

c     Allocate activation lists:

      allocate (mstr%mse2ige_actHF(mstr%msect),
     .          mstr%mse2ige_actDF(mstr%msect),
     .          mstr%msc2igc_actHF(mstr%mscct),
     .          mstr%msc2igc_actDF(mstr%mscct),
     .          mstr%ige2mse_actHF(MBEDGE,nBlock),
     .          mstr%ige2mse_actDF(MBEDGE,nBlock),
     .          mstr%igb_act(nBlock),
     .          mstr%igf_act(MBFACE,nBlock),
     .          mstr%ige_act(MBEDGE,nBlock),
     .          df%c_type(MLCRN,df%nDF),
     .          df%e_type(MLEDGE,df%nDF))

c     Activate hyper and design faces:

      call activateMasterInfo (lc, nBlock, maxDim,
     .                         bc%typ, bc%ir, df, hf,
     .                         mstr%msect, mstr%msecnct,
     .                         mstr%mscct, mstr%msccnct,
     .                         mstr%ige2mse, mstr%igc2msc,
     .                         mstr%mse2igb, mstr%mse2ige,
     .                         mstr%msc2igb, mstr%msc2ige,
     .                         mstr%mse2ige_actHF, mstr%mse2ige_actDF,
     .                         mstr%msc2igc_actHF, mstr%msc2igc_actDF,
     .                         mstr%ige2mse_actHF, mstr%ige2mse_actDF,
     .                         mstr%igb_act, mstr%igf_act, mstr%ige_act)

c     Create hyperface pointers:

      call createHyperFacePtr (hf, bc%is, nBlock, mesh%ijk, mesh%ipx)

c     Create design face pointers:

      call createDesignFacePtr (df, lc, nBlock,
     .                          mesh%ijk, mesh%ipx,
     .                          bc%typ, bc%algn, bc%is, bc%ir,
     .                          mstr%mscct, mstr%msect, mstr%igc2msc,
     .                          mstr%ige2mse, mstr%ige2mse_dir)

      end subroutine createMasterInfo

c+-----------------------------------------------------------------------
c
      subroutine createTransformMat (is, it, ir,
     .                               bcType, bcAlign, nBlock)
c
c    Description:  Create the home, neighbor, and relative
c    transformation matrices which describe the face-face connectivity
c    and relative coordinate orientations.
c
c    Assumes six faces and three coordinate directions.
c
c    Home Transformation Matrix:  is(3,3,6)
c     Establishes the mapping of the global coordinates (i,j,k) to the
c     local face coordinates (l,m).
c     First Row:  nonzero element indicates constant computational
c       coordinate (cc) of the face.  E.g.  (0 1 0) -> J-Constant face
c     Second Row: non-zero element indicates the l local coordinate mapping.
c       E.g. (1 0 0) indicates that on a local face with local indices (l,m),
c       the I-coord is associated with l.
c     Third Row:  Same as second row but for the m local coordinate.
c
c    Neighbor Transformation Matrix:  it(3,3,6,nBlock)
c      From the boundary condition alignment data, maps which coordinate
c      (and direction by sign) of the neighbor block corresponds to the
c      home block.  Each of the three rows represents a home coordinate
c      direction: it(1,,,)->I, it(2,,,)->J, it(3,,,)->K.  For each row,
c      the non-zero column indicates the mapped neighbor coordinate direction.
c      E.g.  [0 1 0]   would be the matrix resulting from a bc: 2 -3 1
c            [0 0-1]   which is interpreted as: i_home -> j_nghbr
c            [1 0 0]                            j_home -> k_nghbr in opposing
c                                                           direction
c                                               k_home -> i_nghbr
c
c      Each face usually maps to a separate neighbor, hence six transformation
c      matrices are needed for each block.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs

      implicit none

c  Arguments:

      integer, intent(in) ::
     .  nBlock

      integer ::
     .  is(MBCRD,MBCRD,MBFACE),
     .  it(MBCRD,MBCRD,MBFACE,nBlock),
     .  ir(MBCRD,MBCRD,MBFACE,nBlock)

      integer, intent(in) ::
     .  bcType(MBFACE,nBlock),
     .  bcAlign(MBCRD,MBFACE,nBlock)

c  Local variables:

      integer ::
     .  Stat, k, l, m, n,
     .  mt(MBFACE)

      character (len=18) ::
     .  routine = 'createTransformMat'

c  Execution:

      write (StdOut,'(4x,a)') 'Assigning transformation matrices...'
      call flush (StdOut)

c     Create home transformation matrix:

      is = 0

      do k = 1, MBFACE
        if (k <= 2) then
          is(1,1,k) = 1
          is(2,2,k) = 1
          is(3,3,k) = 1
        else if (k<=4) then   !Note that j-faces are left handed
          is(1,2,k) = 1       !in local space. i.e. (l,m) => (i,k)
          is(2,1,k) = 1       !instead of (l,m) => (k,i)
          is(3,3,k) = 1
        else
          is(1,3,k) = 1
          is(2,1,k) = 1
          is(3,2,k) = 1
        end if
      end do

c     Create neighbor transformation matrix for flow-through faces:

      it = 0

      do n = 1, nBlock
        do k = 1, MBFACE
          if (bcType(k,n) > 0) then !Block-Block connection
            do m = 1,3
              l = abs(bcAlign(m,k,n))
              if (l<1 .or. l>3) goto 100
              it(m,l,k,n) = sign (1,bcAlign(m,k,n))
            end do
          else  !Specified boundary condition
            do m = 1,3
              it(m,m,k,n) = 1
            end do
          end if
        end do  !Faces
      end do  !Blocks

c     Create relative transformation matrix.

      do k = 1,5,2
        mt(k) = 1
        mt(k+1) = -1
      end do

      ir = 0

      do n = 1, nBlock
        do k = 1, MBFACE
          ir(1:MBCRD,1:MBCRD,k,n) =
     .      matmul (is(1:MBCRD,1:MBCRD,k),
     .              it(1:MBCRD,1:MBCRD,k,n))
          ir(1,1:MBCRD,k,n) = mt(k)*ir(1,1:MBCRD,k,n)
        end do
      end do

      return

c  Error handling:

 100  write (StdOut,*) 'Error in block: ', n
      call meshwarpGracefulExit (routine, 'Boundary conditions error')

      end subroutine createTransformMat

c+-----------------------------------------------------------------------
c
      subroutine createMasterCorners (lc, nBlock, maxDim, msect,
     .                                msecnct, msep, mse2igb, mse2ige,
     .                                mse2igep, mse2ige_dir, mse2msc,
     .                                ige2mse, ige2mse_dir,
     .                                mscct, msccnct, msc2igb, msc2ige,
     .                                msc2igc, msc2igcp, igc2msc)
c
c  Description:  Establish the master corner indexing based on the
c  master edge lists. The master edge and corner arrays are described
c  in the type MasterData module.  The descriptions of arrays being
c  initialized by mstcrn are repeated below with dimension information
c  included.
c
c  integer, pointer ::
c    msccnct(nBCrn*nBlock)
c      Number of block corners connected to master corner.
c    mse2msc(2,msect)
c      Each master edge contains two master corners.
c    igc2msc(nBCrn,nBlock)
c      Pointers from a block corner to the master corner.
c    msc2ige(MAXCRNCNCT,nBCrn*nBlock)
c      Each block corner connected to a master corner belongs to a
c      specific edge.
c    msc2igc(MAXCRNCNCT,nBCrn*nBlock)
c      Each master corner is connected to a series of block corners.
c    msc2igcp(MAXCRNCNCT,nBCrn*nBlock)
c      Pointers into x array for each connected block corner of
c      each master corner.
c    msc2igb(MAXCRNCNCT,nBCrn*nBlock)
c      Each block corner connected to a master corner
c      belongs to a specific block.
c
c  integer ::
c    mscct  !Total number of master corners.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs

      implicit none

c  Arguments:

      type (UnitBlockNumbering), intent (in) ::
     .  lc

      integer, intent (in) ::
     .  nBlock, maxDim

      integer, intent (in) ::
     .  msect,
     .  msecnct(MBEDGE*nBlock),
     .  mse2igep(MAXECNCT,maxDim),
     .  mse2igb(MAXECNCT,MBEDGE*nBlock),
     .  mse2ige(MAXECNCT,MBEDGE*nBlock),
     .  mse2ige_dir(MAXECNCT,MBEDGE*nBlock),
     .  ige2mse(MBEDGE,nBlock),
     .  ige2mse_dir(MBEDGE,nBlock),
     .  msep(MBEDGE*nBlock)

      integer :: !intent(out)
     .  mse2msc(2,msect),
     .  msccnct(MBCRN*nBlock),
     .  msc2igb(MAXCRNCNCT,MBCRN*nBlock),
     .  msc2ige(MAXCRNCNCT,MBCRN*nBlock),
     .  msc2igc(MAXCRNCNCT,MBCRN*nBlock),
     .  msc2igcp(MAXCRNCNCT,MBCRN*nBlock),
     .  igc2msc(MBCRN,nBlock)

      integer :: !intent(out)
     .  mscct

c  Local variables:

      integer ::
     .  n, j, l, msc, j_b, j_c, j_cnct,
     .  i_e1, i_end1, mse1, mse1_end, mse1_cnct, msc1,
     .  i_e2, i_end2, mse2, mse2_end, mse2_cnct, msc2,
     .  i_e3, i_end3, mse3, mse3_end, mse3_cnct, msc3

c  Execution:

      write (StdOut,'(4x,a)') 'Assigning master corners...'
      call flush (StdOut)

c     Initialize master corner arrays:
      mse2msc  = 0
      msccnct  = 0
      msc2igb  = 0
      msc2ige  = 0
      msc2igc  = 0
      msc2igcp = 0
      igc2msc  = 0
      mscct    = 0

      do n = 1, nBlock
      do l = 1, MBCRN

c       Determine the three master edges and master edge endpoints to which
c       this block corner is connected.

        i_e1       =  lc%c2ei(1,l)
        i_end1     =  lc%c2ej(1,l)
        mse1       =  ige2mse(i_e1,n)
        mse1_end   = -ige2mse_dir(i_e1,n)*i_end1
        mse1_cnct  =  msecnct(mse1)
        msc1       =  mse2msc( (mse1_end+3)/2, mse1)

        i_e2       =  lc%c2ei(2,l)
        i_end2     =  lc%c2ej(2,l)
        mse2       =  ige2mse(i_e2,n)
        mse2_end   = -ige2mse_dir(i_e2,n)*i_end2
        mse2_cnct  =  msecnct(mse2)
        msc2       =  mse2msc( (mse2_end+3)/2, mse2)

        i_e3       =  lc%c2ei(3,l)
        i_end3     =  lc%c2ej(3,l)
        mse3       =  ige2mse(i_e3,n)
        mse3_end   = -ige2mse_dir(i_e3,n)*i_end3
        mse3_cnct  =  msecnct(mse3)
        msc3       =  mse2msc( (mse3_end+3)/2, mse3)

        msc = 0
        if (msc1 > 0) then
          igc2msc(l,n) = msc1
          msc          = msc1
        else if (msc2 > 0) then
          igc2msc(l,n) = msc2
          msc          = msc2
        else if (msc3 > 0) then ! Mark omitted this test
          igc2msc(l,n) = msc3
          msc          = msc3
        end if

c       Test if the corner has previously been linked to the master
c       corner list.

        if (msc > 0) then

c         Check to see if the other corners are attached to a master
c         corner and whether they are attached to the same master corner.

          if (msc1 > 0 .and. msc2 > 0 .and. msc1 /= msc2) then

c           Move the pointers from the second master corner to
c           the first master corner

            do j = 1, msccnct(msc2)
              j_b = msc2igb(j,msc2)
              j_c = msc2igc(j,msc2)
              igc2msc(j_c,j_b) = msc  !reassign to new mscrn

              j_cnct = msccnct(msc) + 1
              msccnct(msc) = j_cnct

              msc2ige(j_cnct,msc)  = msc2ige(j,msc2)
              msc2igb(j_cnct,msc)  = msc2igb(j,msc2)
              msc2igc(j_cnct,msc)  = msc2igc(j,msc2)
              msc2igcp(j_cnct,msc) = msc2igcp(j,msc2)
            end do

            mse2msc((mse2_end+3)/2, mse2) = msc1
            msccnct(msc2) = 0
          end if

c         Move corner 3 to corner 1 if necessary

          if (msc1 > 0 .and. msc3 > 0 .and. msc1 /= msc3) then

c           Move the pointers from the third master corner to
c           the first master corner

            do j = 1, msccnct(msc3)
              j_b = msc2igb(j,msc3)
              j_c = msc2igc(j,msc3)
              igc2msc(j_c,j_b) = msc  !reassign to new mscrn
              j_cnct = msccnct(msc1) + 1
              msccnct(msc1) = j_cnct

              msc2ige(j_cnct,msc1)  = msc2ige(j,msc3)
              msc2igb(j_cnct,msc1)  = msc2igb(j,msc3)
              msc2igc(j_cnct,msc1)  = msc2igc(j,msc3)
              msc2igcp(j_cnct,msc1) = msc2igcp(j,msc3)
            end do

            mse2msc((mse3_end+3)/2, mse3) = msc1
            msccnct(msc3) = 0
          end if

c         Move corner 3 to corner 2 if necessary

          if (msc2 > 0 .and. msc3 > 0 .and. msc2 /= msc3 .and.
     .        msc1 == 0) then

c           Move the pointers from the third master corner to
c           the second master corner

            do j = 1, msccnct(msc3)
              j_b = msc2igb(j,msc3)
              j_c = msc2igc(j,msc3)
              igc2msc(j_c,j_b) = msc  !reassign to new mscrn
              j_cnct = msccnct(msc2) + 1
              msccnct(msc2) = j_cnct

              msc2ige(j_cnct,msc2)  = msc2ige(j,msc3)
              msc2igb(j_cnct,msc2)  = msc2igb(j,msc3)
              msc2igc(j_cnct,msc2)  = msc2igc(j,msc3)
              msc2igcp(j_cnct,msc2) = msc2igcp(j,msc3)
            end do

            mse2msc((mse3_end+3)/2, mse3) = msc2
            msccnct(msc3) = 0
          end if

c         Connect any corners which have not been previously
c         connected to a master corner.  Should this be done first?

          if (msc1 == 0) then
            call loadMasterCornerPtrs (lc, nBlock, maxDim, msect, msc,
     .                                 mse1, mse1_end, mse1_cnct, msep,
     .                                 mse2igb, mse2ige, mse2igep,
     .                                 mse2ige_dir, mse2msc,
     .                                 msc2igb, msc2ige, msc2igc,
     .                                 msc2igcp, msccnct)
          end if

          if (msc2 == 0) then
            call loadMasterCornerPtrs (lc, nBlock, maxDim, msect, msc,
     .                                 mse2, mse2_end, mse2_cnct, msep,
     .                                 mse2igb, mse2ige, mse2igep,
     .                                 mse2ige_dir, mse2msc,
     .                                 msc2igb, msc2ige, msc2igc,
     .                                 msc2igcp, msccnct)
          end if

          if (msc3 == 0) then
            call loadMasterCornerPtrs (lc, nBlock, maxDim, msect, msc,
     .                                 mse3, mse3_end, mse3_cnct, msep,
     .                                 mse2igb, mse2ige, mse2igep,
     .                                 mse2ige_dir, mse2msc,
     .                                 msc2igb, msc2ige, msc2igc,
     .                                 msc2igcp, msccnct)
          end if

        else  !No connections for this corner so create new master corner

          mscct = mscct + 1
          igc2msc(l,n) = mscct

c         Link all three master edges to this master corner
          call loadMasterCornerPtrs (lc, nBlock, maxDim, msect, mscct,
     .                               mse1, mse1_end, mse1_cnct, msep,
     .                               mse2igb, mse2ige, mse2igep,
     .                               mse2ige_dir, mse2msc, msc2igb,
     .                               msc2ige, msc2igc,
     .                               msc2igcp, msccnct)

          call loadMasterCornerPtrs (lc, nBlock, maxDim, msect, mscct,
     .                               mse2, mse2_end, mse2_cnct, msep,
     .                               mse2igb, mse2ige, mse2igep,
     .                               mse2ige_dir, mse2msc, msc2igb,
     .                               msc2ige, msc2igc,
     .                               msc2igcp, msccnct)

          call loadMasterCornerPtrs (lc, nBlock, maxDim, msect, mscct,
     .                               mse3, mse3_end, mse3_cnct, msep,
     .                               mse2igb, mse2ige, mse2igep,
     .                               mse2ige_dir, mse2msc, msc2igb,
     .                               msc2ige, msc2igc,
     .                               msc2igcp, msccnct)
        end if

       end do !MBCRN
       end do !nBlock

       end subroutine createMasterCorners

c+-----------------------------------------------------------------------
c
      subroutine createMasterEdges (lc, nBlock, ijk, ip, maxDim,
     .                              bcType, bcAlign, is, ir,
     .                              msep, msecnct, mse2igep,
     .                              mse2igb, mse2ige, mse2ige_dir,
     .                              ige2mse, ige2mse_dir,
     .                              msect, mseptct)
c
c
c  Description:  Create the master edge arrays based on the block
c  boundary connectivity information.  The type MasterData module contains
c  the full list of the variables to be set in this routine and is repeated
c  below with the array dimensions added:
c
c  integer ::
c    msep(nBEdge*nBlock)
c      Beginning index of a master edge in xedge and mse2igep arrays.
c    msecnct(nBEdge*nBlock)
c      Number of edges connected to a master edge.
c    mse2igep(MAXECNCT,nBEdge*nBlock*maxval (ijk))
c      Pointers to the master edge's connected block edges in the x array.
c    mse2igb(MAXECNCT,nBEdge*nBlock)
c      Each master edge's connected block edges belong to a specific
c      block number.
c    mse2ige(MAXECNCT,nBEdge*nBlock)
c      The local edge number of each block edge connected to a master edge.
c    mse2ige_dir(MAXECNCT,nBEdge*nBlock)
c      Each master edge may run in the same (+) or opposite (-) direction
c      than that of the connected block edge.
c    ige2mse(nBEdge,nBlock)
c      Specifies the master edge to which a particular block edge is
c      connected.
c    ige2mse_dir(nBEdge,nBlock)
c      Each block edge may run in the same (+) or opposite (-) direction
c      compared with the master edge.
c  integer ::
c    msect    Total number of master edges.
c    mseptct  Total number of master edge points.
c
c  Warnings:
c    Currently, the maximum number of edges connected to a master edge
c    has been hardcoded at 7. This should be made dynamic, but entails
c    a multiple pass strategy.
c
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs

      implicit none

c  Arguments:

c     Mesh indices:
      integer, intent (in) ::
     .  nBlock, maxDim,
     .  ijk(MBCRD,nBlock),
     .  ip(nBlock)

c     Transformation matrices:
      integer, intent (inout) ::
     .  bcType(MBFACE,nBlock)

      integer, intent (in) ::
     .  bcAlign(MBCRD,MBFACE,nBlock),
     .  is(MBCRD,MBCRD,MBFACE),
     .  ir(MBCRD,MBCRD,MBFACE,nBlock)

c     Unit block numbering convention:
      type (UnitBlockNumbering), intent(in) ::
     .  lc

c     Master Edge arrays:
      integer :: !intent(out)
     .  msep(MBEDGE*nBlock),
     .  msecnct(MBEDGE*nBlock),
     .  mse2igep(MAXECNCT,MBEDGE*nBlock*maxDim),
     .  mse2igb(MAXECNCT,MBEDGE*nBlock),
     .  mse2ige(MAXECNCT,MBEDGE*nBlock),
     .  mse2ige_dir(MAXECNCT,MBEDGE*nBlock),
     .  ige2mse(MBEDGE,nBlock),
     .  ige2mse_dir(MBEDGE,nBlock)

      integer ::
     .  msect,
     .  mseptct

c  Local variables:

      integer ::
     .  i, j, k, l, m, n, p,
     .  itmp, loc(1),
     .  trace,       !Trace of a given relative transform matrix
     .  itr1, itr2,  !Local face l,m indices correspond to a neighbor i,j,k
     .  mse1,        !Each block edge may be connected to two faces which
     .  mse2,        ! could lead to two master edge connections.
     .  mse3,        !The edge could also be degenerate in which case it
     .  mse,         ! has a third opportunity to be attached to a master edge.
     .  mse1_cnct,   !Number of connections for the first and second master
     .  mse2_cnct,   ! edges.
     .  mse3_cnct,
     .  mse_cnct,
     .  mse_p,
     .  mse_p0,
     .  mse2_p_s,
     .  mse2_p_e,
     .  mse3_p_s,
     .  mse3_p_e,
     .  nb_b1,       !Edge belongs to two faces which may be attached to
     .  nb_b2,       ! two neighboring blocks, or to its own block in
     .  nb_b3,       ! degenerate cases.
     .  nb_e1,       !Edge is connected to a neighbor edge.
     .  nb_e2,       ! or through degeneracy to its own block edge.
     .  nb_e3,
     .  lc_f1,       !Edge is in two local faces
     .  lc_f2,
     .  mseptct0,    !Temporary holder for master edge point count.
     .  i_b, i_f, i_e, iflip, !Target block, face, and edge
     .  j_b, j_f, j_e, jflip, !Source block, face, and edge
     .  kflip

      integer, allocatable ::
     .  ieflag(:,:), !Each edge of a block may be degenerate.
     .  irsum1(:,:), !Row sums of relative transform matrix.
     .  irsum2(:,:),
     .  nbf2e(:,:,:) !Each of the four local edges of each face of a given
                     !block is mated to a neighbor edge.

      character (len=6) ::
     .  routine = 'mstedg'

c  Execution:

      write (StdOut,'(4x,a)') 'Assigning master edges...'
      call flush (StdOut)

c     Allocate and initialize the local arrays:

      allocate (ieflag(MBEDGE,nBlock))
      ieflag = 0
      allocate (irsum1(MBFACE,nBlock))
      irsum1 = 0
      allocate (irsum2(MBFACE,nBlock))
      irsum2 = 0
      allocate (nbf2e(MLEDGE,MBFACE,nBlock))
      nbf2e = 0

c     Initialize master edge arrays:

      msep        = 0
      msecnct     = 0
      mse2igep    = 0
      mse2igb     = 0
      mse2ige     = 0
      mse2ige_dir = 0
      ige2mse     = 0
      ige2mse_dir = 0
      msect       = 0
      mseptct     = 0

c     For a given edge of a given face of a given block, determine
c     the connected neighbors and convert to the neighbors' local
c     edge numbering and orientation using the ir transformation matrix.

      do n = 1, nBlock

        do k = 1, MBFACE

c       If boundary condition indicates a block-block connection,
c       identify the mating edge numbers of the neighbor blocks.

          if (bcType(k,n) > 0) then

c           Find the non-zero element in the row indicating coordinate direction
            loc = maxval (maxloc (abs (ir(1,1:MBCRD,k,n))))
            l   = loc (1) !Easier to reference scalar

            p = 2*l- (1 - sign (1, ir(1,l,k,n)))/2 !Determine if max or min face
            nbf2e(1:4,k,n) = lc%f2ei(1:4,p)

            select case (l)   !Find trace of condensed 2x2 neighbor matrix
              case (1)        !in order to identify swapped local coords.
                itr1 = 2
                itr2 = 3
              case (2)      !But for this left handed case, could
                itr1 = 1    !get rid of this case statement.
                itr2 = 3
              case (3)
                itr1 = 1
                itr2 = 2
              case default
                goto 100
            end select

c           Determine whether the local coordinates of the home
c           and connected face are swapped (l,m)_home => (m,l)_neigh
c           Also test whether the edges are flipped. l => -l?
c           A swapped condition is revealed by the trace of the
c           condensed 2x2 neighbor matrix. A flipped condition is
c           revealed by a negative row sum.

            trace = abs (ir(2,itr1,k,n)*ir(3,itr2,k,n))
            irsum1(k,n) = ir(2,itr1,k,n) + ir(2,itr2,k,n) !Save these for
            irsum2(k,n) = ir(3,itr1,k,n) + ir(3,itr2,k,n) !future use.

            if (trace == 0) then
              itmp         = nbf2e(1,k,n)
              nbf2e(1,k,n) = nbf2e(3,k,n)
              nbf2e(3,k,n) = itmp

              itmp         = nbf2e(2,k,n)
              nbf2e(2,k,n) = nbf2e(4,k,n)
              nbf2e(4,k,n) = itmp
            end if

            if (irsum1(k,n) < 0) then
              itmp         = nbf2e(3,k,n)
              nbf2e(3,k,n) = nbf2e(4,k,n)
              nbf2e(4,k,n) = itmp
            end if

            if ( irsum2(k,n) < 0 ) then
              itmp          = nbf2e(1,k,n)
              nbf2e(1,k,n) = nbf2e(2,k,n)
              nbf2e(2,k,n) = itmp
            end if

c         The current local face may have collapsed edges, which cause
c         two or more of the edges of this face to be coincident.
c         Examples of this type of degeneration are faces that are
c         either line or point singularities.  Both cases must be
c         such that they may be treated properly when the master edge
c         list is created.  The input assumptions are as follows:
c
c             ITYPE = -11    -----        edge 1 and 2 are coincident
c             ITYPE = -12    -----        edge 3 and 4 are coincident
c             ITYPE = -13    -----        edge 1 and 2, and
c                                         edge 3 and 4 are coincident
c
c         Modify the boundary condition to be extrapolated.
c         Possibly unsuitable for viscous case?

          else if (bcType(k,n) == -11) then
            ieflag(lc%f2ei(1,k),n) = lc%f2ei(2,k)
            ieflag(lc%f2ei(2,k),n) = lc%f2ei(1,k)
            bcType(k,n) = -1
          else if (bcType(k,n) == -12) then
            ieflag(lc%f2ei(3,k),n) = lc%f2ei(4,k)
            ieflag(lc%f2ei(4,k),n) = lc%f2ei(3,k)
            bcType(k,n) = -1
          else if (bcType(k,n) == -13) then
            ieflag(lc%f2ei(1,k),n) = lc%f2ei(2,k)
            ieflag(lc%f2ei(2,k),n) = lc%f2ei(1,k)
            ieflag(lc%f2ei(3,k),n) = lc%f2ei(4,k)
            ieflag(lc%f2ei(4,k),n) = lc%f2ei(3,k)
            bcType(k,n) = -1
          end if  !bcType

        end do ! nFace
      end do !nBlock

c  *********************************************************************
c  End part 1, begin part 2.  Having determined neighboring edge
c  numbers, begin constructing master edge list.
c  *********************************************************************

      do n = 1, nBlock
        do l = 1, MBEDGE

          mse1 = 0
          mse2 = 0
          mse3 = 0

c         Each edge is a member of two faces, both of which may
c         be connected to a neighboring block.

          lc_f1 = lc%e2fi(1,l)
          nb_b1 = bcType(lc_f1,n)

          if (nb_b1 > 0) then
            nb_e1 = nbf2e(lc%e2fj(1,l),lc_f1,n)
            mse1  = ige2mse(nb_e1,nb_b1)
          end if

          lc_f2 = lc%e2fi(2,l)
          nb_b2 = bcType(lc_f2,n)

          if (nb_b2 > 0) then
            nb_e2 = nbf2e(lc%e2fj(2,l),lc_f2,n)
            mse2  = ige2mse(nb_e2,nb_b2)
          end if

          nb_e3 = ieflag(l,n)

          if (nb_e3 > 0) then !The edge is on a collapsed face
            nb_b3 = n
            mse3  = ige2mse(nb_e3,nb_b3)
          end if

c         Check to see whether this edge in one of its neighboring block
c         incarnations has already been attached to a master edge.  If
c         so, connect the current incarnation to the existing master edge.

          if (mse1 > 0 .or. mse2 > 0 .or. mse3 > 0) then

            if (mse1 > 0) then
              i_b = nb_b1
              i_f = lc_f1
              i_e = nb_e1
              mse = mse1
            else if (mse2 > 0) then
              i_b = nb_b2
              i_f = lc_f2
              i_e = nb_e2
              mse = mse2
            else
              i_b = nb_b3
              i_e = nb_e3
              mse = mse3
            end if

            ige2mse(l,n) = mse                !Connect edge to master list
            msecnct(mse) = msecnct(mse) + 1   !# edges attached to master edge

c           Determine the direction of the current edge relative to
c           the adjacent block edge.  If one looks at the unit cube
c           numbering convention, it turns out that on faces 1-2,
c           edges of numbers 1-8 are l-coord while 9-12 are m coords.
c           Similarly for faces 3-6, edges between 1-4 are l coords,
c           while 5-12 are m coords.  irsum1 represents the relative
c           direction of the l-coord, while irsum2 represents the m-coord.

            if (mse1 > 0 .or. mse2 > 0) then
              if (i_f <= 2) then
                if (l <= 8) then
                  iflip = irsum1(i_f,n)
                else
                  iflip = irsum2(i_f,n)
                end if
              else
                if (l <= 4) then
                  iflip = irsum1(i_f,n)
                else
                  iflip = irsum2(i_f,n)
                end if
              end if
            else
              iflip = 1
            end if

c           Set the direction in the master to global flip list
c           and the global to master flip list. Note that the
c           total allowable connections to a master edge is 7.

            ige2mse_dir(l,n) = iflip*ige2mse_dir(i_e,i_b)
            mse2ige_dir(msecnct(mse),mse) = ige2mse_dir(l,n)

c           Set the master to global edge and block pointers.

            mse2ige(msecnct(mse),mse) = l
            mse2igb(msecnct(mse),mse) = n

            mse_p  = msep(mse)
            mse_p0 = mse_p

            call loadMasterEdgePtrs (lc%e2fi(1,l), mse_p, mse2igep,
     .                               msecnct(mse), ige2mse_dir(l,n),
     .                               is, lc%sng(1,lc%e2fj(1,l)),
     .                               ijk(1,n), ip(n))

c           If the edge is connected to two blocks, check to ensure that both
c           connections refer to the same master edge.  If a redundant master
c           edge is discovered, consolidate it into a single master edge.

            if (mse1 > 0 .and. mse2 > 0 .and. mse1 /= mse2 ) then

c             Shuffle mse2 onto mse1.
              mse1_cnct = msecnct(mse1)
              mse2_cnct = msecnct(mse2)
              msecnct(mse1) = mse1_cnct + mse2_cnct

              if (msecnct(mse1) > MAXECNCT) then
                write (StdOut, '(1x, a, 2i6)')
     .            'Block, edge: ', n, l,
     .            'nb_b1, nb_b2:', nb_b1, nb_b2,
     .            'mse1, mse2:  ', mse1, mse2,
     .            'mse1_cnct, mse2_cnct:', mse1_cnct, mse2_cnct,
     .            'mse, msecnct(mse1):  ', mse, msecnct(mse1),
     .            'MAXECNCT exceeded.  Check the connectivity file.'
                go to 100
              end if

              mse2_p_s = msep(mse2) + 1
              mse2_p_e = msep(mse2+1)

c             Determine the orientation between the two master edges.

              if (lc_f2 <= 2) then
                if (l <= 8) then
                  jflip = irsum1(lc_f2,n)
                else
                  jflip = irsum2(lc_f2,n)
                end if
              else
                if (l <= 4) then
                  jflip = irsum1(lc_f2,n)
                else
                  jflip = irsum2(lc_f2,n)
                end if
              end if

              kflip = jflip*iflip

c             Shuffle the second master edge data onto the first.

              do m = 1, mse2_cnct
                j_b = mse2igb(m,mse2)
                j_e = mse2ige(m,mse2)
                ige2mse(j_e,j_b) = mse1

                ige2mse_dir(j_e,j_b) =    kflip*ige2mse_dir(j_e,j_b)
                mse2ige_dir(mse1_cnct+m,mse1) = ige2mse_dir(j_e,j_b)

                mse2igb(mse1_cnct+m,mse1) = j_b
                mse2ige(mse1_cnct+m,mse1) = j_e

                if (jflip == 1) then
                  j = mse2_p_s
                  do i = mse_p0+1, mse_p
                    mse2igep(mse1_cnct+m,i) = mse2igep(m,j)
                    j = j+1
                  end do
                else
                  j = mse2_p_e
                  do i = mse_p0+1, mse_p
                    mse2igep(mse1_cnct+m,i) = mse2igep(m,j)
                    j = j-1
                  end do
                end if
              end do  !mse2_cnct

              msecnct(mse2) = 0  !Nullify the second master edge connection
            end if  !mse1 /= mse2 /= 0

c           Check for the possibility that the edge is connected to a
c           neighboring block and also collapsed.  If the two edges are
c           connected to multiple master edges, combine the master edges.

            if ((mse1 > 0 .or. mse2 > 0) .and.
     .          (mse3 > 0              ) .and.
     .          (mse3 /= mse           )) then

c             Shuffle mse3 onto mse.
              mse_cnct  = msecnct(mse)
              mse3_cnct = msecnct(mse3)
              msecnct(mse) = mse_cnct + mse3_cnct

              mse3_p_s = msep(mse3) + 1
              mse3_p_e = msep(mse3+1)

              kflip = iflip

              do m = 1, mse3_cnct
                j_b = mse2igb(m,mse3)
                j_e = mse2ige(m,mse3)
                ige2mse(j_e,j_b) = mse

                ige2mse_dir(j_e,j_b) = kflip*ige2mse_dir(j_e,j_b)
                mse2ige_dir(mse_cnct+m,mse) = ige2mse_dir(j_e,j_b)

                mse2igb(mse_cnct+m,mse) = j_b
                mse2ige(mse_cnct+m,mse) = j_e

                j = mse3_p_s
                do i = mse_p0+1, mse_p
                 mse2igep(mse_cnct+m,i) = mse2igep(m,j)
                 j = j+1
                end do
              end do !mse3_cnct

              msecnct(mse3) = 0
            end if

          else

c           This edge does not appear in the master edge list, so add it now.

            msect = msect + 1
            ige2mse(l,n) = msect
            ige2mse_dir(l,n) = 1

            mseptct0 = mseptct

            msecnct(msect) = 1
            mse2ige(1,msect) = l
            mse2igb(1,msect) = n
            mse2ige_dir(1,msect) = 1

            call loadMasterEdgePtrs (lc%e2fi(1,l), mseptct, mse2igep,
     .                               msecnct(msect), ige2mse_dir(l,n),
     .                               is, lc%sng(1,lc%e2fj(1,l)),
     .                               ijk(1,n), ip(n))

            msep(msect+1) = mseptct

          end if   !msedi1 > 0 .or. mseid2 > 0 .or. mseid3 > 0

        end do !nBEdge

      end do !nBlock

      return

c  Error handling::

 100  call meshwarpGracefulExit (routine, 'Error in array (ir?).')

      end subroutine createMasterEdges

c+-----------------------------------------------------------------------
c
      subroutine facxd (bcType, bcAlign, is, ir,
     .                  nBlock, ipx, ijkx,
     .                  iDFace, iDBlock, ipdf_s, pqr_df, nDFpt,
     .                  dftrg2xp, dftrg2nbxp)
c
c     Description: Load the target design face x-pointer and the attached
c     face x-pointer.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c------------------------------------------------------------------------

c  Modules:

      use Parameters

      implicit none

c  Arguments:

      integer, intent (in) ::
     .  nBlock, iDFace, iDBlock,
     .  ipDF_s, nDFpt

      integer, intent (in) ::
     .  bcType(MBFACE,nBlock),
     .  bcAlign(MBCRD,MBFACE,nBlock),
     .  is(MBCRD,MBCRD,MBFACE),
     .  ir(MBCRD,MBCRD,MBFACE,nBlock),
     .  ipx(nBlock),
     .  ijkx(MBCRD,nBlock),
     .  pqr_df(MBCRD)

      integer, intent (out) ::
     .  dftrg2xp(nDFpt),
     .  dftrg2nbxp(nDFpt)

c   Local variables

      integer ::
     .  ps,
     .  q, qs, qe, qf,
     .  r, rs, re, rf,
     .  ip, ipx_s,
     .  pqr_nb(3),
     .  nb_Block

      integer :: one(3) = 1, ms(3)

      integer :: is_t(3,3)

c  Execution:

c     Set up indexing of block

      ps = mod(iDFace-1,2)*(pqr_df(1)-1) + 1

      qs = 1
      qe = pqr_df(2)
      rs = 1
      re = pqr_df(3)

c     In-lined the fload subroutine.  Load the pointers into iftem:

      ipx_s = ipx(iDBlock)
      ip = ipdf_s

      do r = rs, re
      do q = qs, qe
        dftrg2xp(ip) = (ipx_s) +         (  ps*is(1,1,iDFace)
     .                                    + q *is(2,1,iDFace)
     .                                    + r *is(3,1,iDFace)-1)
     . + ijkx(1,iDBlock)*                (  ps*is(1,2,iDFace)
     .                                    + q *is(2,2,iDFace)
     .                                    + r *is(3,2,iDFace)-1)
     . + ijkx(1,iDBlock)*ijkx(2,iDBlock)*(  ps*is(1,3,iDFace)
     .                                    + q *is(2,3,iDFace)
     .                                    + r *is(3,3,iDFace)-1)
        ip = ip + 1
      end do
      end do

c     Set up neighbor block pointers:
c     Find orientation:

      nb_Block = bcType(iDFace,iDBlock)
      ipx_s = ipx(nb_Block)

      pqr_nb = matmul(abs(ir(1:MBCRD,1:MBCRD,iDFace,iDBlock)),
     .                ijkx(1:3,nb_Block))

      ms = matmul(ir(1:MBCRD,1:MBCRD,iDFace,iDBlock),one)

      ps = (( ms(1) +1)/2)*(pqr_nb(1) -1) + 1

      qs = ((-ms(2) +1)/2)*(pqr_nb(2) -1) + 1
      qe = (( ms(2) +1)/2)*(pqr_nb(2) -1) + 1
      qf = ms(2)

      rs = ((-ms(3) +1)/2)*(pqr_nb(3) -1) + 1
      re = (( ms(3) +1)/2)*(pqr_nb(3) -1) + 1
      rf = ms(3)

      ipx_s = ipx(nb_Block)
      ip = ipdf_s

      is_t = abs(ir(1:MBCRD,1:MBCRD,iDFace,iDBlock))
      do r = rs, re, rf
      do q = qs, qe, qf
        dftrg2nbxp(ip) = (ipx_s) + (   ps*is_t(1,1)
     .                               + q *is_t(2,1)
     .                               + r *is_t(3,1)-1)
     .   + ijkx(1,nb_Block)*       (   ps*is_t(1,2)
     .                               + q *is_t(2,2)
     .                               + r *is_t(3,2)-1)
     .    + ijkx(1,nb_Block)*ijkx(2,nb_Block)*(   ps*is_t(1,3)
     .                               + q *is_t(2,3)
     .                               + r *is_t(3,3)-1)
        ip = ip + 1
      end do
      end do

      end subroutine facxd

c+-----------------------------------------------------------------------
c
      subroutine facxv(is, iFace, ipdf_s, ipx_s, pqr, ijk, iftem, nPt)
c
c     Description: Set xp indices of a design face and the connected
c     neighbor face.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c------------------------------------------------------------------------

c  Modules:

      use Parameters

      implicit none

c  Arguments:

      integer, intent (in) ::
     .  iFace,
     .  ipdf_s,
     .  ipx_s,
     .  nPt

      integer, intent (in) ::
     .  is(MBCRD,MBCRD,MBFACE)

      integer, intent (in) ::
     .  pqr(3),
     .  ijk(3)

      integer, intent (out) ::
     .  iftem(nPt)

c   Local variables:

      integer ::
     .  ps,
     .  q, qs, qe,
     .  r, rs, re,
     .  ip

c  Execution:

c     Set up indexing of block

      ps = mod(iFace-1,2)*(pqr(1)-1) + 1

      qs = 1
      qe = pqr(2)
      rs = 1
      re = pqr(3)

c     In-lined the fload subroutine.  Load the pointers into iftem:

      ip = ipdf_s
      do r = rs, re
      do q = qs, qe
        iftem(ip) = (ipx_s) +  (   ps*is(1,1,iFace)
     .                           + q *is(2,1,iFace)
     .                           + r *is(3,1,iFace)-1)
     .         + ijk(1)*       (   ps*is(1,2,iFace)
     .                           + q *is(2,2,iFace)
     .                           + r *is(3,2,iFace)-1)
     .         + ijk(1)*ijk(2)*(   ps*is(1,3,iFace)
     .                           + q *is(2,3,iFace)
     .                           + r *is(3,3,iFace)-1)
        ip = ip + 1
      end do
      end do

      end subroutine facxv

c+--------------------------------------------------------------------------
c
      subroutine loadMasterCornerPtrs (lc, nBlock, maxDim, msect,
     .                                 msc, mse, mse_end,
     .                                 mse_cnct, msep,
     .                                 mse2igb, mse2ige, mse2igep,
     .                                 mse2ige_dir, mse2msc,
     .                                 msc2igb, msc2ige, msc2igc,
     .                                 msc2igcp, msccnct)
c
c     Load msccnct:  master-corner connection counter
c          msc2ige:  master-corner-to-global-edge array
c          msc2igb:  master-corner-to-global-block array
c          msc2igcp: master-corner-to-global-xyz-index array
c          msc2igc:  master-corner-to-global-corner array
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c
c---------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs

      implicit none

c  Arguments:

      integer, intent (in) ::
     .  nBlock, maxDim, msect

      type (UnitBlockNumbering) ::
     .  lc

      integer, intent(in) ::
     .  msc, mse_cnct, mse, mse_end

      integer, intent(in) ::
     .  mse2ige(MAXECNCT,MBEDGE*nBlock),
     .  mse2igb(MAXECNCT,MBEDGE*nBlock),
     .  mse2igep(MAXECNCT,MBEDGE*nBlock*maxDim),
     .  mse2ige_dir(MAXECNCT,MBEDGE*nBlock),
     .  msep(MBEDGE*nBlock)

      integer, intent(inout) ::
     .  mse2msc(2,msect),
     .  msccnct(MBCRN*nBlock),
     .  msc2igb(MAXCRNCNCT,MBCRN*nBlock),
     .  msc2ige(MAXCRNCNCT,MBCRN*nBlock),
     .  msc2igc(MAXCRNCNCT,MBCRN*nBlock),
     .  msc2igcp(MAXCRNCNCT,MBCRN*nBlock)

c  Local variables:

      integer ::
     .  j, j_s, j_e, j_p,
     .  ige, ige_end, igc,
     .  j_cnct

      character (len=20) ::
     .  routine = 'loadMasterCornerPtrs'

c  Execution:

      do j = 1, mse_cnct
        msccnct(msc) = msccnct(msc) + 1
        j_cnct = msccnct(msc)

        if (j_cnct > MAXCRNCNCT) goto 100

        msc2ige(j_cnct,msc) = mse2ige(j,mse)
        msc2igb(j_cnct,msc) = mse2igb(j,mse)

        j_s = msep(mse) + 1
        j_e = msep(mse+1)
        j_p =   j_e*(mse_end+1)/2  !Grid point either at beginning
     .        + j_s*(1-mse_end)/2  !or end of master edge x pointer list.

        msc2igcp(j_cnct,msc) = mse2igep(j,j_p)

        ige     = mse2ige(j,mse)
        ige_end = mse2ige_dir(j,mse)*mse_end
        igc     = lc%e2ci( (ige_end+3)/2, ige)
        msc2igc(j_cnct,msc) = igc
      end do

c     Connect the master edge to the master corner

      mse2msc((mse_end+3)/2, mse) = msc

      return

c  Error handling:

 100  call meshwarpGracefulExit (routine, 'MSCCNCT(MSC) > MAXCRNCNCT')

      end subroutine loadMasterCornerPtrs

c+-----------------------------------------------------------------------
c
      subroutine loadMasterEdgePtrs (iFace, iPtr, mse2igep, iCnct,
     .                               iEdgFlp, is, isng, ijk, ip)
c
c  Description:   Calculate the pointers for all boundary edges for
c  block to block exchanges.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c------------------------------------------------------------------------

c  Modules:

      use Parameters

      implicit none

c  Arguments:

      integer ::
     .  iFace, iPtr, iCnct, iEdgFlp

      integer ::
     .  mse2igep(MAXECNCT,*),
     .  isng(MLEDGE),
     .  is(MBCRD,MBCRD,MBFACE),
     .  ijk(MBCRD),
     .  ip

c  Local variables:

      integer ::
     .  ps,
     .  q, qs, qe,
     .  r, rs, re,
     .  tmp

      integer ::
     .  pqr(MBCRD)

c  Execution:

c     Find orientation of coordinates on home face

      pqr = matmul(is(1:MBCRD,1:MBCRD,iFace),ijk)

c     Set up indexing of block

      ps = mod(iFace-1,2)*(pqr(1)-1) + 1

      qs = ((1+isng(1))  + pqr(2)*(1-isng(1)))/2
      qe = ((1+isng(2))  + pqr(2)*(1-isng(2)))/2
      rs = ((1+isng(3))  + pqr(3)*(1-isng(3)))/2
      re = ((1+isng(4))  + pqr(3)*(1-isng(4)))/2

c     Flip direction if necessary.

      if (iEdgFlp < 0) then
       tmp = qs
       qs  = qe
       qe  = tmp

       tmp = rs
       rs  = re
       re  = tmp
      end if

c     In-lined the eload subroutine.  Load the pointers into mstoigp.

      do r = rs, re, iEdgFlp
      do q = qs, qe, iEdgFlp
        iPtr = iPtr + 1
        mse2igep(iCnct,iPtr) = (ip) +  (   ps*is(1,1,iFace)
     .                                   + q *is(2,1,iFace)
     .                                   + r *is(3,1,iFace) - 1)
     .                 + ijk(1)*       (   ps*is(1,2,iFace)
     .                                   + q *is(2,2,iFace)
     .                                   + r *is(3,2,iFace) - 1)
     .                 + ijk(1)*ijk(2)*(   ps*is(1,3,iFace)
     .                                   + q *is(2,3,iFace)
     .                                   + r *is(3,3,iFace) - 1)
      end do
      end do

      end subroutine loadMasterEdgePtrs

c+-----------------------------------------------------------------------
c
      subroutine setLocalNumConvention(lc)
c
c  Description:  Establish the local face/edge/corner numbering convention.
c
c  Face Numbering:  IMin = 1, IMax = 2
c                   JMin = 3, JMax = 4
c                   KMin = 5, KMax = 6
c
c  Edge Numbering:  This is difficult to draw on the screen so I
c                   am going to furnish 2D views of each face and
c                   let the user draw it in 3D.
c
c  Face: 1
c           5------7
c           |  7   |
c           |9   11|
c      k,m  |  5   |
c      ^    1------3
c      |
c      -->j,l
c
c  Face: 2
c           6------8
c           |  8   |
c           |10  12|
c      k,m  |  6   |
c      ^    2------4
c      |
c      -->j,l
c
c  Face: 3
c           5------6
c           |  3   |
c           |9   10|
c      k,m  |  1   |
c      ^    1------2
c      |
c      -->i,l
c
c  Face: 4
c           7------8
c           |  4   |
c           |11  12|
c      k,m  |  2   |
c      ^    3------4
c      |
c      -->i,l
c
c  Face: 5
c           3------4
c           |  2   |
c           |5    6|
c      j,m  |  1   |
c      ^    1------2
c      |
c      -->i,l
c
c  Face: 6
c           7------8
c           |  4   |
c           |7    8|
c      j,m  |  3   |
c      ^    5------6
c      |
c      -->i,l
c
c  Original implementation J. Reuther
c  Adapted to MESHWARP by M. Rimlinger
c
c------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs

      implicit none

c  Arguments:

      type (UnitBlockNumbering) :: lc

c  Execution:

      write (StdOut,'(4x,a)') 'Assigning unit block convention...'
      call flush (StdOut)

c     For a local edge of a local face, indicate the indexing endpoints
c     in (l,m).
      lc%sng(1,1)    =  1  !Generic Edge 1 has a range of l=1,lmax
      lc%sng(2,1)    = -1  !                              m=1,1
      lc%sng(3,1)    =  1
      lc%sng(4,1)    =  1
      lc%sng(1,2)    =  1  !Edge 2 l=1,lmax
      lc%sng(2,2)    = -1  !       m = mmax,mmax
      lc%sng(3,2)    = -1
      lc%sng(4,2)    = -1
      lc%sng(1,3)    =  1
      lc%sng(2,3)    =  1
      lc%sng(3,3)    =  1
      lc%sng(4,3)    = -1
      lc%sng(1,4)    = -1
      lc%sng(2,4)    = -1
      lc%sng(3,4)    =  1
      lc%sng(4,4)    = -1

c     Each face is composed of four block edges
      lc%f2ei(1,1)   = 5
      lc%f2ei(2,1)   = 7
      lc%f2ei(3,1)   = 9
      lc%f2ei(4,1)   = 11
      lc%f2ei(1,2)   = 6
      lc%f2ei(2,2)   = 8
      lc%f2ei(3,2)   = 10
      lc%f2ei(4,2)   = 12
      lc%f2ei(1,3)   = 1
      lc%f2ei(2,3)   = 3
      lc%f2ei(3,3)   = 9
      lc%f2ei(4,3)   = 10
      lc%f2ei(1,4)   = 2
      lc%f2ei(2,4)   = 4
      lc%f2ei(3,4)   = 11
      lc%f2ei(4,4)   = 12
      lc%f2ei(1,5)   = 1
      lc%f2ei(2,5)   = 2
      lc%f2ei(3,5)   = 5
      lc%f2ei(4,5)   = 6
      lc%f2ei(1,6)   = 3
      lc%f2ei(2,6)   = 4
      lc%f2ei(3,6)   = 7
      lc%f2ei(4,6)   = 8

c     Each face is composed of four block corners
      lc%f2ci(1,1)   = 1
      lc%f2ci(2,1)   = 3
      lc%f2ci(3,1)   = 5
      lc%f2ci(4,1)   = 7
      lc%f2ci(1,2)   = 2
      lc%f2ci(2,2)   = 4
      lc%f2ci(3,2)   = 6
      lc%f2ci(4,2)   = 8
      lc%f2ci(1,3)   = 1
      lc%f2ci(2,3)   = 2
      lc%f2ci(3,3)   = 5
      lc%f2ci(4,3)   = 6
      lc%f2ci(1,4)   = 3
      lc%f2ci(2,4)   = 4
      lc%f2ci(3,4)   = 7
      lc%f2ci(4,4)   = 8
      lc%f2ci(1,5)   = 1
      lc%f2ci(2,5)   = 2
      lc%f2ci(3,5)   = 3
      lc%f2ci(4,5)   = 4
      lc%f2ci(1,6)   = 5
      lc%f2ci(2,6)   = 6
      lc%f2ci(3,6)   = 7
      lc%f2ci(4,6)   = 8

c     Each edge is a member of two block faces
      lc%e2fi(1,1)  = 5
      lc%e2fi(2,1)  = 3
      lc%e2fi(1,2)  = 5
      lc%e2fi(2,2)  = 4
      lc%e2fi(1,3)  = 6
      lc%e2fi(2,3)  = 3
      lc%e2fi(1,4)  = 6
      lc%e2fi(2,4)  = 4
      lc%e2fi(1,5)  = 5
      lc%e2fi(2,5)  = 1
      lc%e2fi(1,6)  = 5
      lc%e2fi(2,6)  = 2
      lc%e2fi(1,7)  = 6
      lc%e2fi(2,7)  = 1
      lc%e2fi(1,8)  = 6
      lc%e2fi(2,8)  = 2
      lc%e2fi(1,9)  = 3
      lc%e2fi(2,9)  = 1
      lc%e2fi(1,10) = 3
      lc%e2fi(2,10) = 2
      lc%e2fi(1,11) = 4
      lc%e2fi(2,11) = 1
      lc%e2fi(1,12) = 4
      lc%e2fi(2,12) = 2

c     Location of each block edge in the local face
      lc%e2fj(1,1)  = 1
      lc%e2fj(2,1)  = 1
      lc%e2fj(1,2)  = 2
      lc%e2fj(2,2)  = 1
      lc%e2fj(1,3)  = 1
      lc%e2fj(2,3)  = 2
      lc%e2fj(1,4)  = 2
      lc%e2fj(2,4)  = 2
      lc%e2fj(1,5)  = 3
      lc%e2fj(2,5)  = 1
      lc%e2fj(1,6)  = 4
      lc%e2fj(2,6)  = 1
      lc%e2fj(1,7)  = 3
      lc%e2fj(2,7)  = 2
      lc%e2fj(1,8)  = 4
      lc%e2fj(2,8)  = 2
      lc%e2fj(1,9)  = 3
      lc%e2fj(2,9)  = 3
      lc%e2fj(1,10) = 4
      lc%e2fj(2,10) = 3
      lc%e2fj(1,11) = 3
      lc%e2fj(2,11) = 4
      lc%e2fj(1,12) = 4
      lc%e2fj(2,12) = 4

c     Each edge is composed of two block corners
      lc%e2ci(1,1) = 1
      lc%e2ci(2,1) = 2
      lc%e2ci(1,2) = 3
      lc%e2ci(2,2) = 4
      lc%e2ci(1,3) = 5
      lc%e2ci(2,3) = 6
      lc%e2ci(1,4) = 7
      lc%e2ci(2,4) = 8
      lc%e2ci(1,5) = 1
      lc%e2ci(2,5) = 3
      lc%e2ci(1,6) = 2
      lc%e2ci(2,6) = 4
      lc%e2ci(1,7) = 5
      lc%e2ci(2,7) = 7
      lc%e2ci(1,8) = 6
      lc%e2ci(2,8) = 8
      lc%e2ci(1,9) = 1
      lc%e2ci(2,9) = 5
      lc%e2ci(1,10)= 2
      lc%e2ci(2,10)= 6
      lc%e2ci(1,11)= 3
      lc%e2ci(2,11)= 7
      lc%e2ci(1,12)= 4
      lc%e2ci(2,12)= 8

c     Each corner forms an end of three edges
      lc%c2ei(1,1) = 1
      lc%c2ei(2,1) = 5
      lc%c2ei(3,1) = 9
      lc%c2ei(1,2) = 1
      lc%c2ei(2,2) = 6
      lc%c2ei(3,2) = 10
      lc%c2ei(1,3) = 2
      lc%c2ei(2,3) = 5
      lc%c2ei(3,3) = 11
      lc%c2ei(1,4) = 2
      lc%c2ei(2,4) = 6
      lc%c2ei(3,4) = 12
      lc%c2ei(1,5) = 3
      lc%c2ei(2,5) = 7
      lc%c2ei(3,5) = 9
      lc%c2ei(1,6) = 3
      lc%c2ei(2,6) = 8
      lc%c2ei(3,6) = 10
      lc%c2ei(1,7) = 4
      lc%c2ei(2,7) = 7
      lc%c2ei(3,7) = 11
      lc%c2ei(1,8) = 4
      lc%c2ei(2,8) = 8
      lc%c2ei(3,8) = 12

c     Each corner is either the beginning or end of the
c     three edges
      lc%c2ej(1,1) = 1
      lc%c2ej(2,1) = 1
      lc%c2ej(3,1) = 1
      lc%c2ej(1,2) =-1
      lc%c2ej(2,2) = 1
      lc%c2ej(3,2) = 1
      lc%c2ej(1,3) = 1
      lc%c2ej(2,3) =-1
      lc%c2ej(3,3) = 1
      lc%c2ej(1,4) =-1
      lc%c2ej(2,4) =-1
      lc%c2ej(3,4) = 1
      lc%c2ej(1,5) = 1
      lc%c2ej(2,5) = 1
      lc%c2ej(3,5) =-1
      lc%c2ej(1,6) =-1
      lc%c2ej(2,6) = 1
      lc%c2ej(3,6) =-1
      lc%c2ej(1,7) = 1
      lc%c2ej(2,7) =-1
      lc%c2ej(3,7) =-1
      lc%c2ej(1,8) =-1
      lc%c2ej(2,8) =-1
      lc%c2ej(3,8) =-1

      end subroutine setLocalNumConvention

c----------------------------------------------------------------------------
c
c     Mesh Perturbation Routines
c
c+-----------------------------------------------------------------------
c
      subroutine dfaces (df, meshN, mesh0)
c
c     Description:  Update the design faces to follow the hyperface motion.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c
c------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs

      implicit none

c  Arguments:

      type (DesignFace) ::
     .  df

      type (AbstractMeshType) ::
     .  mesh0,
     .  meshN

c  Local variables:

      integer ::
     .  i, j, l, m, l_c, l_e, nv,
     .  ac, accnct, i_dfcxp,
     .  ae, aecnct, ip_ae, ip_de,
     .  id1, id2, id3,
     .  i_df, ipdf_s, ipdf_e, nPt,
     .  isrc, itrg, npts,
     .  ip, ip_s, ip_e

      real ::
     .  eps,
     .  total(3),
     .  dxc(3,4)

      real, allocatable ::
     .  xtemp(:), ytemp(:), ztemp(:),
     .  uDF0(:), vDF0(:),
     .  dx_dfsrc(:,:),
     .  dx_dftrg(:,:),
     .  ac_xden(:,:,:),
     .  ac_xnum(:,:,:),
     .  ac_xdes(:,:,:),
     .  ae_xden(:,:,:),
     .  ae_xnum(:,:,:),
     .  ae_xdes(:,:,:),
     .  dxej(:,:,:),
     .  dxek(:,:,:)

c  Execution:

      write (StdOut,'(a)')
     .  ' Updating design faces from modified hyperfaces...'
      call flush (StdOut)

      eps = epsilon(eps)

c     Allocate the temporary arrays for the source and target delta x:

      npts = df%ip(df%nDF+1)-1
      allocate(uDF0(npts))
      uDF0 = 0
      allocate(vDF0(npts))
      vDF0 = 0
      allocate(dx_dfsrc(3,npts))
      dx_dfsrc = 0
      allocate(dx_dftrg(3,npts))

c     Allocate the averaging corner/edge arrays:

      allocate (ac_xden(3,MAXCRNCNCT,df%acct))
      ac_xden = 0
      allocate (ac_xnum(3,MAXCRNCNCT,df%acct))
      ac_xnum = 0
      allocate (ac_xdes(3,MAXCRNCNCT,df%acct))
      ac_xdes = 0

      allocate (ae_xden(3,MAXECNCT,df%nDFEpts))
      ae_xden = 0
      allocate (ae_xnum(3,MAXECNCT,df%nDFEpts))
      ae_xnum = 0
      allocate (ae_xdes(3,MAXECNCT,df%nDFEpts))
      ae_xdes = 0

c     Load the averaging corner/edge arrays with contributions from all the
c     design faces:

      do i_df = 1, df%nDF

c       Indexing for this design face:
        id2 = df%ijk(2,i_df)
        id3 = df%ijk(3,i_df)
        ipdf_s = df%ip(i_df)
        ipdf_e = df%ip(i_df+1)-1
        npts = ipdf_e - ipdf_s + 1

c       Create a parameterization of the original design face:

        allocate (xtemp(npts), ytemp(npts), ztemp(npts))

        do l = ipdf_s, ipdf_e
          isrc = df%dfsrc2xp(l)
          m  = l-ipdf_s+1
          xtemp(m) = mesh0%x(1,isrc)
          ytemp(m) = mesh0%x(2,isrc)
          ztemp(m) = mesh0%x(3,isrc)
        end do

        uDF0(ipdf_s) = 0.  !-999. suppresses normalization

        call param2d (id2, id3, 1, id2, 1, id3,
     .                xtemp, ytemp, ztemp,
     .                uDF0(ipdf_s), vDF0(ipdf_s) )

        deallocate (xtemp, ytemp, ztemp)

c       Load the initial perturbations into the source/target arrays:

        do l = ipdf_s, ipdf_e
          isrc = df%dfsrc2xp(l)
          itrg = df%dftrg2xp(l)
          dx_dfsrc(1:MBCRD,l) =  meshN%x(1:MBCRD,isrc)
     .                         - mesh0%x(1:MBCRD,isrc)
          dx_dftrg(1:MBCRD,l) =  meshN%x(1:MBCRD,itrg)
     .                         - mesh0%x(1:MBCRD,itrg)
        end do

c       Load the averaging corner arrays from the perturbed faces:

        do l_c = 1, MLCRN !ok
          ac = df%dfc2ac(l_c,i_df)
          accnct = df%dfc2accnct(l_c,i_df)
          i_dfcxp = df%dfc2dfxp(l_c,i_df)

          if (df%c_type(l_c,i_df) == 1) then !Attached to hyperface corner.
            ac_xdes(1:MBCRD,accnct,ac) = dx_dftrg(1:MBCRD,i_dfcxp)
          else
            do nv = 1,3
             ac_xden(nv,accnct,ac) =
     .           ac_xden(nv,accnct,ac) + abs(dx_dfsrc(nv,i_dfcxp))
             ac_xnum(nv,accnct,ac) =
     .           ac_xnum(nv,accnct,ac) +
     .           dx_dfsrc(nv,i_dfcxp)*abs(dx_dfsrc(nv,i_dfcxp))
            end do
          end if
        end do

c       Load the averaging edge arrays from the perturbed faces:

        do l_e = 1, MLEDGE
          ae = df%dfe2ae(l_e,i_df)
          ip_ae = df%aeip(ae)
          aecnct = df%dfe2aecnct(l_e,i_df)
          ip_de = df%dfeip(l_e,i_df)

          if (df%e_type(l_e,i_df) == 1) then !Attached to hyperface edge.
            do l = 1, df%dfenPt(l_e,i_df)
              i = df%dfe2dfxp(ip_de)
              j = df%dfe2aexp(ip_ae)
              do nv = 1,3
                ae_xdes(nv,aecnct,j) = dx_dftrg(nv,i)
              end do
              ip_de = ip_de + 1
              ip_ae = ip_ae + 1
            end do
          else
            do l = 1, df%dfenPt(l_e,i_df)
              i = df%dfe2dfxp(ip_de)
              j = df%dfe2aexp(ip_ae)
              do nv = 1,3
                ae_xden(nv,aecnct,j) = ae_xden(nv,aecnct,j) +
     .                                abs(dx_dfsrc(nv,i))
                ae_xnum(nv,aecnct,j) = ae_xnum(nv,aecnct,j) +
     .                           dx_dfsrc(nv,i)*abs(dx_dfsrc(nv,i))
              end do
              ip_de = ip_de + 1
              ip_ae = ip_ae + 1
            end do
          end if
        end do
      end do  !Finished loading the averaging corner and edge arrays

c     Perform special averaging for all corner and edge values
c     to get a single contribution from each block.

c     Average corners:
      do i_df = 1, df%nDf
        do l_c = 1, MLCRN  !ok
          if (df%c_type(l_c,i_df) /= 1) then !average
            ac = df%dfc2ac(l_c,i_df)
            accnct = df%dfc2accnct(l_c,i_df)
            do nv = 1,3
              ac_xdes(nv,accnct,ac) = ac_xnum(nv,accnct,ac) /
     .                                max (ac_xden(nv,accnct,ac), eps)
            end do
          end if
        end do
      end do

      do ac = 1, df%acct
        total = 0
        do accnct = 1, df%accnct(ac)
          do nv = 1,3
            total(nv) = total(nv) + ac_xdes(nv,accnct,ac)
          end do
        end do
        do nv = 1,3
          total(nv) = total(nv)/real(df%accnct(ac))
          do accnct = 1, df%accnct(ac)
            ac_xdes(nv,accnct,ac) = total(nv)
          end do
        end do
      end do

c     Average edges:
      do i_df = 1, df%nDf
        do l_e = 1, MLEDGE
          if (df%e_type(l_e,i_df) /= 1) then
            ae = df%dfe2ae(l_e,i_df)
            ip_ae = df%aeip(ae)
            aecnct = df%dfe2aecnct(l_e,i_df)
c            ip_de = df%dfeip(l_e,i_df)

            do l = 1, df%dfenPt(l_e,i_df)
c              i = df%dfe2dfxp(ip_de)
              j = df%dfe2aexp(ip_ae)
              do nv = 1,3
                ae_xdes(nv,aecnct,j) = ae_xnum(nv,aecnct,j) /
     .                                 max (ae_xden(nv,aecnct,j), eps)
              end do
c              ip_de = ip_de + 1
              ip_ae = ip_ae + 1
            end do
          end if
        end do
      end do

      do ae = 1, df%aect
        do ip_ae = df%aeip(ae), df%aeip(ae+1)-1
          total = 0
          do aecnct = 1, df%aecnct(ae)
            do nv = 1,3
              total(nv) = total(nv) + ae_xdes(nv,aecnct,ip_ae)
            end do
          end do

          do nv = 1,3
            total(nv) = total(nv)/real(df%aecnct(ae))
            do aecnct = 1, df%aecnct(ae)
              ae_xdes(nv,aecnct,ip_ae) = total(nv)
            end do
          end do
        end do
      end do

c     Perform design face perturbation based on corner and edge
c     values:

      do i_df = 1, df%nDF

c       Pick up corners and edges:
        do l_c = 1, MLCRN
          ac = df%dfc2ac(l_c,i_df)
          accnct = df%dfc2accnct(l_c,i_df)
          do nv = 1,3
            dxc(nv,l_c) = ac_xdes(nv,accnct,ac)
          end do
        end do

        nPt = df%dfenPt(1,i_df)
        allocate (dxej(3,2,nPt))

        do l_e = 1,2
          ae = df%dfe2ae(l_e,i_df)
          ip_ae = df%aeip(ae)
          aecnct = df%dfe2aecnct(l_e,i_df)

          do l = 1, df%dfenPt(l_e,i_df)
            j = df%dfe2aexp(ip_ae)
            do nv = 1,3
              dxej(nv,l_e,l) = ae_xdes(nv,aecnct,j)
            end do
            ip_ae = ip_ae + 1
          end do
        end do

        nPt = df%dfenPt(3,i_df)
        allocate (dxek(3,2,nPt))

        do l_e = 3,4
          ae = df%dfe2ae(l_e,i_df)
          ip_ae = df%aeip(ae)
          aecnct = df%dfe2aecnct(l_e,i_df)

          do l = 1, df%dfenPt(l_e,i_df)
            j = df%dfe2aexp(ip_ae)
            do nv = 1,3
              dxek(nv,l_e-2,l) = ae_xdes(nv,aecnct,j)
            end do
            ip_ae = ip_ae + 1
          end do
        end do

c       Correct the design face:
        ipdf_s = df%ip(i_df)
        call warpface (dxc, dxej, dxek, df%ijk(2,i_df), df%ijk(3,i_df),
     .                 dx_dfsrc(1,ipdf_s), uDF0(ipdf_s), vDF0(ipdf_s))

c       Add the perturbations to the target faces:
        ipdf_e = df%ip(i_df+1)-1
        do l = ipdf_s, ipdf_e
          itrg = df%dftrg2xp(l)
          do nv = 1,3
            meshN%x(nv,itrg) = mesh0%x(nv,itrg) + dx_dfsrc(nv,l)
          end do
        end do
        deallocate (dxek, dxej)

      end do

c     Update block faces connected to the design faces:
      do i_df = 1, df%nDF
        ip_s = df%ip(i_df)
        ip_e = df%ip(i_df+1)-1

        do ip = ip_s, ip_e
          i = df%dftrg2nbxp(ip)
          j = df%dftrg2xp(ip)
          meshN%x(1:3,i) = meshN%x(1:3,j)
        end do
      end do

      end subroutine dfaces

c+-----------------------------------------------------------------------
c
      subroutine mapXYZ (meshN, mesh0, geomN, geom0, hf, uv, uvType)
c
c     Description:  Perturb hyperfaces based on the (u,v) map
c     and the geometries.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c     07/17/03 DAS  2D searches must be done for both the original and the
c                   perturbed patch (u,v)s.  Don't assume the same patch
c                   sizes either.
c-------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs

      implicit none

c  Arguments:

      type (AbstractMeshType), intent (in) ::
     .  mesh0

      type (AbstractMeshType), intent (out) ::
     .  meshN

      type (Geometry), intent (in) ::
     .  geom0,
     .  geomN

      type (HyperFace), intent (in) ::
     .  hf

      type (UVMap), intent (in) ::
     .  uv

      character (*) ::
     .  uvType

c  Local constants:

      real, parameter ::
     .  ONE = 1.

c  Local variables:

      real ::
     .  a(3,2),
     .  xpt_g0(3), xpt_gN(3)

      real ::
     .  eps, p, q, pm1, qm1

      integer ::
     .  ix, i_tri, l, nv,
     .  ptch, prev_ptch,
     .  im_ptch, jm_ptch, iptch_s,
     .  ig, jg, ig0, jg0, igN, jgN,
     .  imjm, ipjm, imjp, ipjp,
     .  ier, i_hf

      character (len=6) ::
     .  routine = 'mapXYZ'

c  Execution:

      write (StdOut,'(a)')
     .  ' Updating hyperfaces from modified geometry...'
      call flush (StdOut)

      if (uvType == 'absolute') then
        do i_hf = 1, hf%nHF                      !For each hyperface
          do l = hf%ip(i_hf), hf%ip(i_hf+1)-1    !loop over each point
            ix = hf%hf2xp(l)    !Point to be updated
            i_tri = uv%iuv(l)   !Geometry triangle to interpolate from
            p     = uv%uv(1,l)
            q     = uv%uv(2,l)

c           Interpolate (x,y,z) from original database triangle
            xpt_g0(1:3) =  (1.-p-q)*geom0%x(1:3,geom0%tri(2,i_tri))
     .                    +(   p  )*geom0%x(1:3,geom0%tri(1,i_tri))
     .                    +(     q)*geom0%x(1:3,geom0%tri(3,i_tri))

c           Interpolate (x,y,z) from updated database triangle
            xpt_gN(1:3) =  (1.-p-q)*geomN%x(1:3,geomN%tri(2,i_tri))
     .                    +(   p  )*geomN%x(1:3,geomN%tri(1,i_tri))
     .                    +(     q)*geomN%x(1:3,geomN%tri(3,i_tri))

c           Construct updated (x,y,z) from original mesh and
c           delta of database positions.
            meshN%x(1:3,ix) = mesh0%x(1:3,ix) + xpt_gN(1:3)-xpt_g0(1:3)
          end do
        end do

      else  !uvType is parametric

ccccc   eps = 10. * epsilon (eps)
        eps = 100.* epsilon (eps) ! Match UV_MAP; should scale by data range.

        prev_ptch = 0

        do i_hf = 1, hf%nHF
          do l = hf%ip(i_hf), hf%ip(i_hf+1)-1

            ix = hf%hf2xp(l)            !Point to be updated
            ptch = uv%iuv(l)            !Patch to interpolate from

c           Interpolate from original geometry first:

            iptch_s = geom0%ipx(ptch)   !Beginning index of patch
            im_ptch = geom0%ijk(1,ptch) !Dimensions of patch
            jm_ptch = geom0%ijk(2,ptch)

            if (ptch /= prev_ptch) then !Don't update prev_ptch yet
*****         ig = max (1, im_ptch/2)
*****         jg = max (1, jm_ptch/2)
              ig = nint (real (im_ptch - 2) * uv%uv(1,l) + ONE)
              jg = nint (real (jm_ptch - 2) * uv%uv(2,l) + ONE)
            else
              ig = ig0
              jg = jg0
            end if

c           Locate enclosing cell in original geometry:

            call ripple2d (im_ptch, jm_ptch, 1, im_ptch, 1, jm_ptch,
     .                     geom0%u_ptch(iptch_s), geom0%v_ptch(iptch_s),
     .                     uv%uv(1,l), uv%uv(2,l), ig, jg, eps, p,q,ier)

            if (ier /= 0) then
              write (StdOut,*)
     .          'MAPXYZ: Error finding geom0 cell containing',
     .          ' u = ', uv%uv(1,l), ' v = ', uv%uv(2,l),
     .          ' Patch #: ', ptch, ' ier: ', ier,
     .          ' Proceeding with cell ', ig, jg
ccccc           goto 100
            end if

            ig0  = ig
            jg0  = jg
            ipjp = iptch_s + ig + im_ptch*jg
            imjp = ipjp - 1
            ipjm = ipjp - im_ptch
            imjm = ipjm - 1
            pm1  = ONE - p
            qm1  = ONE - q

c           Bilinear interpolation of (x,y,z) at this (u,v) in original geom.:

            do nv = 1, 3
              xpt_g0(nv) =
     >          qm1*(pm1*geom0%x(nv, imjm) + p*geom0%x(nv, ipjm))
     >        + q  *(pm1*geom0%x(nv, imjp) + p*geom0%x(nv, ipjp))
            end do

c           Interpolate from perturbed geometry (may be different-sized patch):

            iptch_s = geomN%ipx(ptch)   !Beginning index of patch
            im_ptch = geomN%ijk(1,ptch) !Dimensions of patch
            jm_ptch = geomN%ijk(2,ptch)

            if (ptch /= prev_ptch) then
              prev_ptch = ptch
              ig = nint (real (im_ptch - 2) * uv%uv(1,l) + ONE)
              jg = nint (real (jm_ptch - 2) * uv%uv(2,l) + ONE)
            else
              ig = igN
              jg = jgN
            end if

c           Locate enclosing cell in perturbed geometry:

            call ripple2d (im_ptch, jm_ptch, 1, im_ptch, 1, jm_ptch,
     .                     geomN%u_ptch(iptch_s), geomN%v_ptch(iptch_s),
     .                     uv%uv(1,l), uv%uv(2,l), ig, jg, eps, p,q,ier)

            if (ier /= 0) then
              write (StdOut,*)
     .          'MAPXYZ: Error finding geomN cell containing',
     .          ' u = ', uv%uv(1,l), ' v = ', uv%uv(2,l),
     .          ' Patch #: ', ptch, ' ier: ', ier,
     .          ' Proceeding with cell ', ig, jg
ccccc           goto 100
            end if

            igN  = ig
            jgN  = jg
            ipjp = iptch_s + ig + im_ptch*jg
            imjp = ipjp - 1
            ipjm = ipjp - im_ptch
            imjm = ipjm - 1
            pm1  = ONE - p
            qm1  = ONE - q

c           Bilinear interpolation of (x,y,z) at this (u,v) in perturbed
c           geometry, and update of this surface mesh point:

            do nv = 1, 3
              xpt_gN(nv) =
     >          qm1*(pm1*geomN%x(nv, imjm) + p*geomN%x(nv, ipjm))
     >        + q  *(pm1*geomN%x(nv, imjp) + p*geomN%x(nv, ipjp))

              meshN%x(nv,ix) = mesh0%x(nv,ix)
     .                         + (xpt_gN(nv) - xpt_g0(nv))
            end do

          end do !Points in the hyperface
        end do   !nHF
      end if

      return

c  Error handling:

cc100 call meshwarpGracefulExit (routine,
cccc .  'Error in finding cell containing u,v')

      end subroutine mapxyz

c+-----------------------------------------------------------------------
c
      subroutine update_mse (meshN, mesh0, maxDim,
     .                       msect, msecnct, msep,
     .                       mse2igb, mse2ige, mse2igep,
     .                       mse2ige_act, ige2mse_act,
     .                       mse2msc, msc2igcp, msccnct)
c
c     Description:  Update all the master edge point locations due to
c     hyperface or design motion, then cascade this motion to all attached
c     block edges.
c
c     This routine could be modified if one created the following:
c       hfe2msexp, hfe2igexp,
c       dfe2msexp, dfe2igexp
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs

      implicit none

c  Arguments:

      type (AbstractMeshType), intent (inout) ::
     .  meshN

      type (AbstractMeshType), intent (in) ::
     .  mesh0

      integer ::
     .  msect, maxDim

      integer ::
     .  msecnct(MBEDGE*meshN%nBlock),
     .  msep(MBEDGE*meshN%nBlock),
     .  mse2igb(MAXECNCT,MBEDGE*meshN%nBlock),
     .  mse2ige(MAXECNCT,MBEDGE*meshN%nBlock),
     .  mse2igep(MAXECNCT,MBEDGE*meshN%nBlock*maxDim),
     .  mse2ige_act(msect),
     .  ige2mse_act(MBEDGE,meshN%nBlock),
     .  mse2msc(2,msect),
     .  msc2igcp(MAXCRNCNCT,MBCRN*meshN%nBlock),
     .  msccnct(MBCRN*meshN%nBlock)

c  Local variables:

      real, allocatable ::
     .  x_mse(:,:)

      integer ::
     .  mse, msc, hf_b, hf_e, i_mse,
     .  l_cnct, mse_s, mse_e, ixp

c  Execution:

      allocate (x_mse(3,msep(msect+1)-1))
      x_mse = 0

c     For each master edge that has a flag indicating that it is
c     attached to at least one active hyperface edge: determine
c     the active connection and fill the master edge x array with
c     the modified mesh x values for the active hyperface edge.

c     Update the master edges:

      do mse = 1,msect

        if (mse2ige_act(mse) == 1) then  !Master edge is connected
                                         !to at least one active hyperface
          do l_cnct = 1, msecnct(mse)
            hf_b = mse2igb(l_cnct,mse)
            hf_e = mse2ige(l_cnct,mse)
            if (ige2mse_act(hf_e,hf_b) == 1) exit
          end do

          mse_s = msep(mse) + 1
          mse_e = msep(mse+1)

          do i_mse = mse_s, mse_e
            ixp = mse2igep(l_cnct,i_mse)
            x_mse(1:3,i_mse) =   meshN%x(1:3,ixp)
     .                         - mesh0%x(1:3,ixp)
          end do

c         Having updated the master edge x array due to hyperface motion,
c         update all block edges attached to the master edge.

          do l_cnct = 1, msecnct(mse)

            do i_mse = mse_s, mse_e
              ixp = mse2igep(l_cnct,i_mse)
              meshN%x(1:3,ixp) =   x_mse(1:3,i_mse)
     .                           + mesh0%x(1:3,ixp)
            end do
          end do

c         Update all corners attached to the two master corners:

          msc = mse2msc(1,mse)

          do l_cnct = 1, msccnct(msc)
            ixp = msc2igcp(l_cnct,msc)
            meshN%x(1:3,ixp) =  x_mse(1:3,mse_s)
     .                        + mesh0%x(1:3,ixp)
          end do

          msc = mse2msc(2,mse)

          do l_cnct = 1, msccnct(msc)
            ixp = msc2igcp(l_cnct,msc)
            meshN%x(1:3,ixp) =   x_mse(1:3,mse_e)
     .                         + mesh0%x(1:3,ixp)
          end do

        end if

      end do

      deallocate (x_mse)

      end subroutine update_mse

c+-----------------------------------------------------------------------
c
      subroutine warpmb (meshN, mesh0, mstr, df)
c
c  Description:  Update the master edge/corners based on the mapXYZ
c  perturbations of the hyperfaces.  Update the positions of the design
c  faces based on the hyperface modifications, then redistribute the
c  results through the master edge/corner lists.  Finally, update the
c  block volumes which have faces or edges in the perturbed master
c  edge/corner lists.
c
c  AUTHOR:
c     Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  HISTORY:
c     10/1999  MJR  Version 1: MESHWARP created as a stand-alone version
c                   of SYN107-MB's mesh perturbation method.
c------------------------------------------------------------------------

c  Modules:

      use Parameters
      use TypeDefs

      implicit none

c  Arguments:

      type (AbstractMeshType) :: meshN, mesh0

      type (MasterData) :: mstr

      type (DesignFace) :: df

c  Local variables:

      integer ::
     .  ig, ip, maxDim

c  Execution:

c     Reflect the changes in the hyperface edges to the master edges
c     and then to the block edges connected to the master edges.

      maxDim = maxval (meshN%ijk)

      call update_mse (meshN, mesh0, maxDim,
     .                 mstr%msect, mstr%msecnct, mstr%msep,
     .                 mstr%mse2igb, mstr%mse2ige, mstr%mse2igep,
     .                 mstr%mse2ige_actHF, mstr%ige2mse_actHF,
     .                 mstr%mse2msc, mstr%msc2igcp, mstr%msccnct)

c     Modify the design faces to follow the hyperface motion.
      call dfaces (df, meshN, mesh0)

c     Distribute design face changes through master edge mechanism:
      call update_mse (meshN, mesh0, maxDim,
     .                 mstr%msect, mstr%msecnct, mstr%msep,
     .                 mstr%mse2igb, mstr%mse2ige, mstr%mse2igep,
     .                 mstr%mse2ige_actDF, mstr%ige2mse_actDF,
     .                 mstr%mse2msc, mstr%msc2igcp, mstr%msccnct)

c     Perturb blocks

      write (StdOut,'(a)') ' Warping blocks...'
      call flush (StdOut)

      do ig = 1, meshN%nBlock
        if (meshN%ijk(3,ig) == 1) cycle ! Surface grid only
        if (mstr%igb_act(ig) > 0) then
          ip = meshN%ipx(ig)
          call warpblk_gm (mstr%igf_act(1,ig), mstr%ige_act(1,ig),
     .                     meshN%ijk(1,ig), meshN%ijk(2,ig),
     .                     meshN%ijk(3,ig), mesh0%x(1,ip),
     .                     meshN%x(1,ip))
        end if
      end do

      end subroutine warpmb

c+---------------------------------------------------------------------------
c
      subroutine meshwarpGracefulExit (routine, message)
c
c  Description: Gracefully terminates the program execution by
c  writing an error message, calling mpi_finalize, flushing the standard
c  output buffer, and calling stop.  If not executing in an MPI environment,
c  the MPI_DUMMY library should be linked upon compilation.
c
c  Variables:
c     routine  : Routine requesting the exit.
c     message  : Reason for the termination.
c
c  Constants:
c     StdOut : Package constant = 6.
c
c  Author: Mark Rimlinger, Raytheon/NASA Ames Research Center
c
c  6/21/99  General ioPackage overhaul
c
c----------------------------------------------------------------------------

c  Modules:

      use Parameters

      implicit none

c  Arguments:

      character(*), intent(in) :: routine, message

c  Procedures:

      external mpi_finalize

c  Execution:

c     Write message, flush the buffer, finalize MPI, and quit.

      write (StdOut,'(/,(1x,a))') routine, message,
     .  'MESHWARP TERMINATING DUE TO ERROR'
      call flush (StdOut)
      call mpi_finalize (StdOut)
      stop

      end subroutine meshwarpGracefulExit
