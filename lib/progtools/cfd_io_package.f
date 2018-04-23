c*******************************************************************************
c
c                               CFD_IO_PACKAGE
c
c  CFD_IO_PACKAGE provides a library of routines to perform common file I/O
c  tasks encountered in the development of applications utilizing CFD datasets
c  in PLOT3D, FAST, OVERFLOW, and SYN107-MB formats.  Written in Fortran 90,
c  the package encapsulates all memory allocation through the use of dynamic
c  pointers, thus freeing the user from the onus of memory management.
c
c  Compilation:
c
c  'f90 cfd_io_package -c' produces a module, CFD_IO_PACKAGE.mod,
c  and an object file, cfd_io_package.o both of which must be present for
c  successful compilation of programs using the CFD_IO_PACKAGE routines.
c  In non-MPI environments, the MPI_DUMMY package must be compiled.  For VMS
c  systems, the FLUSH_VMS source dummy must also be compiled.
c
c  Modules:
c
c  The subroutines and derived data types have been encapsulated in a single
c  module named CFD_IO_PACKAGE.  Any user defined procedures calling routines
c  defined in this package must have the 'use cfd_io_package' statement.
c
c  Derived Data Types:
c
c  The ioMeshFileType has been defined below to describe the file and data
c  attributes that cannot be determined from the contents of the file itself.
c  Contents of the derived data type include the file name, file format, file
c  status, multiple vs. single grid, an option for the returned array to be
c  shaped as n-tuples (x(ndim,npts)) or contiguous dimension (x(npts,ndim)),
c  and the number of spatial dimensions (1-D,2-D,...). A generic routine that
c  interactively queries the user for these attributes can be activated by two
c  methods: leaving the file variable out of the calling argument list (optional
c  argument), or passing an unassociated ioMeshFileType variable.  The second
c  method should be used in cases where the file attributes will be used by
c  the calling program or by subsequent calls to I/O package utilities.
c
c  Notes on F90 Platform Compatibility:
c
c  SGI and VMS Fortran 90 differ in the handling of read (err=, iostat=).
c  On the SGI, an erroneous eof/eor does not cause branching to the err label,
c  but iostat will be loaded with a useful negative number.  A formatted read
c  from an unformatted file, however, is treated as a branching offense.  If the
c  err specifier is missing, the code discontinues execution in an ugly fashion.
c  VMS on the other hand, treats both an unexpected eof/eor and a format clash
c  as an error causing the read/write statement to branch.  Thus read/writes
c  have both err and iostat specifiers followed by an explicit test on the
c  returned iostat value.
c
c  History:
c
c     07/13/99  M. Rimlinger  Initial coding of package completed.
c     11/18/99  D. Saunders   VMS compiler bug involving pointers in implied
c                             I/O loops required work-arounds throughout.
c
c  Author: Mark J. Rimlinger, Raytheon/NASA Ames Research Center, Moffett Field.
c
c*******************************************************************************

      module cfd_io_package

        private   !All module variables and routines default to private

c       Public routines:

        public ::
     .    IO_R_AEROSURF,      !Read AEROSURF panels, subpatch information
     .    IO_R_BC_FLO107,     !Read FLO107 boundary conditions
     .    IO_R_FNC_P3D,       !Read PLOT3D function file
     .    IO_R_HYP_SYN107,    !Read SYN107 hyperfaces from connectivity file
     .    IO_R_DSNF_SYN107,   !Read SYN107 design faces from connectivity file
     .    IO_R_Q_P3D,         !Read PLOT3D solution file
     .    IO_R_STR_X_P3D,     !Read PLOT3D structured mesh without IBLANK
     .    IO_R_STR_X_P3D_IB,  !Read PLOT3D structured mesh with IBLANK
     .    IO_R_UNSTR_X_FAST,  !Read FAST unstructured mesh
     .    IO_R_WAKE_FLO107,   !Read FLO107 wake specification
     .    IO_W_FNC_P3D,       !Write PLOT3D function file
     .    IO_W_Q_P3D,         !Write PLOT3D solution file
     .    IO_W_STR_X_P3D,     !Write PLOT3D structured mesh without IBLANK
     .    IO_W_STR_X_P3D_IB,  !Write PLOT3D structured mesh with IBLANK
     .    IO_W_UNSTR_X_FAST   !Write FAST unstructured mesh

c       Module constants:

        character, parameter ::
     .    FORMATTED * 9 = 'formatted', UNFORMATTED * 11 = 'unformatted',
     .    UNKNOWN * 7   = 'unknown',   NEW * 3 = 'new', OLD * 3 = 'old'

        integer, parameter ::
     .    ioStdOut = 6, ioStdIn = 5

c       Module type definitions:

        type, public :: ioFileType
          character :: name * 80
          character :: form * 11
          character :: status * 7
        end type ioFileType

        type, public :: ioMeshFileType
          character :: name * 80  !File name
          character :: form * 11  !File format, 'unformatted'|'formatted'
          character :: status * 7 !File status, 'old'|'new'|'unknown'
          logical   :: mg         !Multiple grid format; .F. -> single block
          logical   :: triples    !Shape of array in memory
                                  !x(ndim,npts) -> .T.
                                  !x(npts,ndim) -> .F.
          integer   :: ndim       !# spatial dimensions, 2 | 3
          logical   :: iblank     !iblanking -> .T.
        end type ioMeshFileType

c       Internal procedures:

        contains

c*******************************************************************************
c                             Public subroutines
c*******************************************************************************

c+------------------------------------------------------------------------------
c
      subroutine IO_R_AEROSURF (x, ng, ijk, ip, pnl2Comp, nSubI, nSubJ,
     .                          iSub, jSub, ncall, usr_file)
c
c  Description: IO_R_AEROSURF reads an AEROSURF panel file with subpatch info.
c
c  The usr_file variable is an optional argument.  Absence, or the passing of an
c  unassociated variable in the actual argument, will cause a predefined set of
c  user queries to activate.
c
c  The ncall variable controls to what extent the file will be read.
c  Ncall = -2 Read only the subpatch info at the end of the file.
c  Ncall = -1 results in the entire file, header, panels, and subpatches,
c             being returned.
c  Ncall = 0  requests that only the header be read.
c  Ncall > 0  requires that the header has been read in a previous call and
c        < ng requests that only the next block's x data be returned.
c  Ncall > ng Read the subpatch information after the last data block
c             has been read sequentially.
c
c  Arguments: See declaration for the complete type specification.
c
c     Intent(inout)
c       usr_file       : File description variable, optional except for
c                        ncall > 0 which requires the file to be previously
c                        opened with the ncall = 0 option
c
c     If ncall = -1    : Read whole file
c      Intent(out)
c       ng             : Number of grids
c       ijk(ndim,ng)   : Computational dimensions of each grid
c       ip(ng+1)       : Index specifying beginning of mesh i in x array;
c                        ip(ng+1) is the total number of points + 1
c       x              : Packed mesh data for ng grids
c                      : x(ndim,npts) => prg_file%triples = .T.
c                      : x(npts,ndim) => prg_file%triples = .F.
c       pnl2Comp(ng)   : Index relating panel to component
c       nSubI(ng)      : Number of subpatches in the I-direction for each panel
c       nSubJ(ng)      : Number of subpatches in the J-direction for each panel
c       iSub(MAXSUBI+1,ng) : Starting I-index of each subpatch
c       jSub(MAXSUBJ+1,ng) : Starting J-index of each subpatch
c
c     If ncall = 0     : Read header of file
c      Intent(out)
c       ng             : Number of grids
c       ijk(ndim,ng)   : Computational dimensions of each grid
c       ip             : Returned as an unassociated variable
c       x              : Returned as an unassociated variable
c
c     If ncall > 0     : Read next data block
c      Intent(in)
c       ng             : Number of grids
c       ijk(ndim,ng)   : Computational dimensions of each grid
c
c      Intent(out)
c       ip(2)          : Degenerate case in which ip(1) = 1 and
c                        ip(2) = (npts in ncall_th block) + 1
c       x              : Packed mesh data for next grid
c                      : x(nq,npts) => prg_file%triples = .T.
c                      : x(npts,nq) => prg_file%triples = .F.
c
c     If ncall = ng+1  : Read subpatch info from end of file.  Assume file is
c                      : initially open & pointing to start of subpatch info.
c      Intent(in)
c       ng             : Number of grids
c      Intent(out)
c       pnl2Comp(ng)   : Index relating panel to component
c       nSubI(ng)      : Number of subpatches in the I-direction for each panel
c       nSubJ(ng)      : Number of subpatches in the J-direction for each panel
c       iSub(MAXSUBI+1,ng) : Starting I-index of each subpatch
c       jSub(MAXSUBJ+1,ng) : Starting J-index of each subpatch
c
c  Side Effects:
c     Since the data pointers are disassociated and allocated during execution,
c     actual argument memory references prior to the call will be lost.
c
c     The file will be closed at the time of return when ncall = -1, -2 or ng+1.
c     However, the file will remain open for ng = [0,ng].
c
c  Warnings:
c     Use of ncall >= 0, requires the presence of the usr_file argument.
c     This does not imply that the argument must be associated at the
c     time of the call when ncall = 0.
c
c     The value of ncall > 0 must correspond to the block number at the
c     current file pointer location.  The routine does not rewind/skip
c     ncall-1 times, but assumes that the current location is correct.
c
c  8/17/99  Initial coding of AEROSURF sub-patch info.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      real, pointer    :: x(:,:)
      integer, pointer :: ijk(:,:), ip(:)
      integer          :: ng, ncall
      integer, pointer :: pnl2Comp(:), nSubI(:), nSubJ(:),
     .                    iSub(:,:), jSub(:,:)
      type (ioMeshFileType), pointer, optional :: usr_file

c  Local constants:

      character, parameter :: routine * 13 = 'IO_R_AEROSURF'

c  Local variables:

      integer ::
     .  i, ioStat, ioUnit, ig, ig_s, ig_e, ig_tot, j, ndim, ni, nj,
     .  npts, maxsubi, maxsubj
      logical ::
     .  ascii, exist, tuples
      type (ioMeshFileType), pointer ::
     .  prg_file

c  Execution:

c     Test if the programmer is using the optional file argument,
c     the unassociated file argument, or preinitialized options.

      call io_Qry_Str_X_File (usr_file, prg_file, ncall, routine)

c     Check the file attributes for errors, then open the file:

      call io_Prep_Open_X_File (prg_file, ioUnit, ncall, routine)

      if (ncall <= 0) then
        inquire (file=prg_file%name, exist=exist)
        if (.not. exist) goto 103
        open (unit=ioUnit, file=prg_file%name,
     .        form=prg_file%form, status=prg_file%status)
      end if

c     If requested, read file header.

      ndim  = prg_file%ndim
      ascii = prg_file%form == FORMATTED

      if (ncall <= 0) then                 !Begin reading header
        if (prg_file%mg) then
          if (ascii) then
            read (ioUnit, *, err=100, iostat=ioStat) ng  !Number of grids
          else
            read (ioUnit, err=100, iostat=ioStat) ng
          end if
          if (ioStat /= 0) goto 100
        else
          ng = 1
        end if

        if (associated (ijk)) deallocate (ijk)  !Allocate ijk and
        allocate (ijk(ndim,ng), stat=ioStat)    !read dimensions
        if (ioStat /= 0) goto 200

        if (ascii) then
          read (ioUnit, *, err=101, iostat=ioStat) ijk
        else
          read (ioUnit, err=101, iostat=ioStat) ijk
        end if
        if (ioStat /= 0) goto 101

        if (ncall == 0) goto 999              !Finished reading header
      end if

      if (ncall > ng) goto 777

c     Set loop counters to the number of blocks being read.

      if (ncall < 0) then
        ig_s = 1
        ig_e = ng
        ig_tot = ng
      else
        ig_s = ncall
        ig_e = ncall
        ig_tot = 1
      end if

c     Allocate ip array.

      if (associated (ip)) deallocate (ip)
      allocate (ip(ig_tot+1), stat=ioStat)
      if (ioStat /= 0) goto 201

c     Set the ip index for the number of blocks being read.

      ip(1) = 1
      do ig = ig_s, ig_e
        npts = 1
        do j = 1, ndim
          npts = npts*ijk(j,ig)
        end do
        ip(ig-ig_s+2) = ip(ig-ig_s+1) + npts
      end do

c     Allocate x array.

      npts   = ip(ig_tot+1) - 1
      tuples = prg_file%triples

      if (tuples) then
        ni = ndim
        nj = npts
      else
        ni = npts
        nj = ndim
      end if

      if (associated (x)) deallocate (x)
      allocate (x(ni,nj), stat=ioStat)
      if (ioStat /= 0) goto 202

c     Read block data in correct format.

      do ig = 1, ig_tot
        call IO_R_BLOCK (ig_tot, ig, ip, ni, nj, 0, x, ip,
     .                   tuples, ascii, ioUnit, ioStat)
        if (ioStat /= 0) goto 102
      end do

 777  if ((ncall < 0) .or. (ncall > ng)) then
        if (associated (pnl2Comp)) deallocate (pnl2Comp)
        if (associated (nSubI))    deallocate (nSubI)
        if (associated (nSubJ))    deallocate (nSubJ)
        if (associated (iSub))     deallocate (iSub)
        if (associated (jSub))     deallocate (jSub)

        if (ascii) then
          allocate (pnl2Comp(ng), stat=ioStat)
          if (ioStat /= 0) goto 203
          read (ioUnit, *, err=104, iostat=ioStat)  !Skip descriptor
          read (ioUnit, *, err=104, iostat=ioStat) pnl2Comp
          if (ioStat /= 0) goto 104

          read (ioUnit, *, err=105, iostat=ioStat)
          read (ioUnit, *, err=105, iostat=ioStat) maxsubi, maxsubj
          if (ioStat /= 0) goto 105

          allocate (nSubI(ng), stat=ioStat)
          if (ioStat /= 0) goto 204
          read (ioUnit, *, err=106, iostat=ioStat)
          read (ioUnit, *, err=106, iostat=ioStat) nSubI
          if (ioStat /= 0) goto 106

          allocate (nSubJ(ng), stat=ioStat)
          if (ioStat /= 0) goto 205
          read (ioUnit, *, err=107, iostat=ioStat)
          read (ioUnit, *, err=107, iostat=ioStat) nSubJ
          if (ioStat /= 0) goto 107

          allocate (iSub(maxsubi+1,ng), stat=ioStat)
          if (ioStat /= 0) goto 206
          read (ioUnit, *, err=108, iostat=ioStat)
          read (ioUnit, *, err=108, iostat=ioStat) iSub
          if (ioStat /= 0) goto 108

          allocate (jSub(maxsubj+1,ng), stat=ioStat)
          if (ioStat /= 0) goto 207
          read (ioUnit, *, err=109, iostat=ioStat)
          read (ioUnit, *, err=109, iostat=ioStat) jSub
          if (ioStat /= 0) goto 109

        else

          allocate (pnl2Comp(ng), stat=ioStat)
          if (ioStat /= 0) goto 203
          read (ioUnit, err=104, iostat=ioStat) pnl2Comp
          if (ioStat /= 0) goto 104

          read (ioUnit, err=105, iostat=ioStat) maxsubi, maxsubj
          if (ioStat /= 0) goto 105

          allocate (nSubI(ng), stat=ioStat)
          if (ioStat /= 0) goto 204
          read (ioUnit, err=106, iostat=ioStat) nSubI
          if (ioStat /= 0) goto 106

          allocate (nSubJ(ng), stat=ioStat)
          if (ioStat /= 0) goto 205
          read (ioUnit, err=107, iostat=ioStat) nSubJ
          if (ioStat /= 0) goto 107

          allocate (iSub(maxsubi+1,ng), stat=ioStat)
          if (ioStat /= 0) goto 206
          read (ioUnit, err=108, iostat=ioStat) iSub
          if (ioStat /= 0) goto 108

          allocate (jSub(maxsubj+1,ng), stat=ioStat)
          if (ioStat /= 0) goto 207
          read (ioUnit, err=109, iostat=ioStat) jSub
          if (ioStat /= 0) goto 109
        end if

        close (ioUnit)
        goto 999
      end if

 999  return

c  Error handling:

 100  call ioExit (routine, 'Read ng returned error')
 101  call ioExit (routine, 'Read ijk returned error')
 102  call ioExit (routine, 'Read x returned error')
 103  call ioExit (routine, 'File does not exist')
 104  call ioExit (routine, 'Read component number returned error')
 105  call ioExit (routine, 'Read maxsubI, maxsubJ returned error')
 106  call ioExit (routine, 'Read nSubI returned error')
 107  call ioExit (routine, 'Read nSubJ returned error')
 108  call ioExit (routine, 'Read iSub returned error')
 109  call ioExit (routine, 'Read jSub returned error')

 200  call ioExit (routine, 'Error in allocation of ijk')
 201  call ioExit (routine, 'Error in allocation of ip')
 202  call ioExit (routine, 'Error in allocation of x')
 203  call ioExit (routine, 'Error in allocation of pnl2Comp')
 204  call ioExit (routine, 'Error in allocation of nSubI')
 205  call ioExit (routine, 'Error in allocation of nSubJ')
 206  call ioExit (routine, 'Error in allocation of iSub')
 207  call ioExit (routine, 'Error in allocation of jSub')

      end subroutine IO_R_AEROSURF

c+------------------------------------------------------------------------------
c
      subroutine IO_R_BC_FLO107 (ibc_type, ibc_orient, rhand,
     .                           ng, ncall, usr_file)
c
c  Description: IO_R_BC_FLO107 reads the block boundary conditions from a
c  FLO107-MB or SYN107-MB connectivity file.
c
c  The usr_file variable is an optional argument.  Absence, or the passing of an
c  unassociated variable in the actual argument, will cause a predefined set of
c  user queries to activate.
c
c  The ncall variable controls to what extent the file will be read.
c  Ncall = -1 results in the entire file, header and all BCs, being returned.
c  Ncall = 0  requests that only the header be read.
c  Ncall > 0  requires that the header has been read in a previous call and
c              requests that only the next block's bc line be returned.
c
c  Note that the header in a connectivity file is only the number of grids, ng.
c
c  Arguments: See declaration for the complete type specification.
c
c     Intent(inout)
c       usr_file        : File description variable, optional unless ncall > 0
c
c     If ncall = -1     : Read whole file
c      Intent(out)
c       ng              : Number of grids
c       ibc_type (6,ng) : Boundary condition for the six faces of each block
c       ibc_orient(3,6,ng) : Boundary coordinate equivalence and orientation
c                          : for the six faces of each block
c       rhand(ng)       : The right or left handedness of each block
c
c     If ncall = 0      : Read header of file
c      Intent(out)
c       ng              : Number of grids
c       ibc_type (:,:)  : Returned as an unassociated variable
c
c       ibc_orient(:,:) : Returned as an unassociated variable
c       rhand(:)        : Returned as an unassociated variable
c
c     If ncall > 0      : Read next data block
c      Intent(in)
c       ng              : Number of grids
c       usr_file        : File descriptor must be associated and file open
c      Intent(out)
c       ibc_type (6,1)  : Boundary condition for the six faces of this block
c       ibc_orient(3,6,1) : Boundary coordinate equivalence and orientation
c                         : for the six faces of this block
c       rhand(ng)       : The right or left handedness of the block
c
c  Side Effects:
c     Since the data pointers are disassociated and allocated during
c     execution, memory references prior to the subroutine call will
c     be lost.
c
c     The file will be closed at the time of return when ncall = -1 or
c     ng.  However, the file will remain open for ng = [0,ng-1].
c
c  Warnings:
c     Use of ncall >= 0, requires the presence of the usr_file argument.
c     This does not imply that the argument must be associated at the
c     time of the call when ncall = 0.
c
c     The value of ncall > 0 must correspond to the block number at the
c     current file pointer location.  The routine does not rewind/skip
c     ncall-1 times, but assumes that the current location is correct.
c
c  6/21/99  Updated ioPackage to allocate required memory for user using
c           the pointer attribute for the actual and dummy arguments.
c           Note: Cannot use intent and pointer in same attribute line!
c
c  7/1/99   Rewrote after adding ncall options.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      integer, pointer ::
     .  ibc_type(:,:),        !Boundary condition applied to each face
     .  ibc_orient(:,:,:)     !Interblock coordinate orientation data

      real, pointer ::
     .  rhand(:)              !Left or right (-1 or 1) handedness of grid

      integer ::
     .  ng,                   !Number of grids
     .  ncall                 !Extent of read: -1, entire file
                              !                 0, header
                              !                >0, next block

      type (ioMeshFileType), pointer, optional ::
     .  usr_file              !File description

c  Local constants:

      integer,   parameter :: nface = 6, ndim = 3
      character, parameter :: routine * 14 = 'IO_R_BC_FLO107'

c  Local variables:

      integer ::
     .  ihold, ig, ig_tot, ioStat, ioUnit, k, l, n
      logical ::
     .  exist
      type (ioMeshFileType), pointer ::
     .  prg_file

c  Execution:

c     Test if the programmer is using the optional file argument,
c     the unassociated file argument, or preinitialized options.

      call io_Qry_File (usr_file, prg_file, ncall, routine)

c     Check the file attributes for errors, then open the file:

      call io_Prep_Open_X_File (prg_file, ioUnit, ncall, routine)

      if (ncall <= 0) then
        inquire (file=prg_file%name, exist=exist)
        if (.not. exist) goto 103
        open (unit=ioUnit, file=prg_file%name,
     .        form=prg_file%form, status=prg_file%status)
      end if

c     If requested, read file header.

      if (ncall <= 0) then                       !Begin reading header
        if (prg_file%form == FORMATTED) then
          read (ioUnit,*, err=100, iostat=ioStat) ng
        else
          read (ioUnit,   err=100, iostat=ioStat) ng
        end if
        if (ioStat /= 0) goto 100
        if (ncall  == 0) goto 999                !Finished reading header
      end if

c     Set the loop counters to the number of blocks being read.

      if (ncall < 0) then
        ig_tot = ng
      else
        ig_tot = 1
      end if

c     Ensure that work space is deallocated prior to being allocated.

      if (associated (ibc_type))   deallocate (ibc_type)
      if (associated (ibc_orient)) deallocate (ibc_orient)
      if (associated (rhand))      deallocate (rhand)

c     Allocate the required workspace.

      allocate (ibc_type(nface,ig_tot), stat=ioStat)
      if (ioStat /= 0) goto 200
      allocate (ibc_orient(ndim,nface,ig_tot), stat=ioStat)
      if (ioStat /= 0) goto 201
      allocate (rhand(ig_tot), stat=ioStat)
      if (ioStat /= 0) goto 202

      if (prg_file%form == FORMATTED) then
        do ig = 1, ig_tot
          read (ioUnit, 500, err=102, iostat=ioStat)
     .      n,(ibc_type(k,ig),(ibc_orient(l,k,ig),l=1,ndim),k=1,nface),
     .      rhand(ig), ihold
          if (ioStat /= 0) goto 102
        end do
      else
        do ig = 1,ig_tot
          read (ioUnit, err=102, iostat=ioStat)
     .      n,(ibc_type(k,ig),(ibc_orient(l,k,ig),l=1,ndim),k=1,nface),
     .      rhand(ig), ihold
          if (ioStat /= 0) goto 102
        end do
      end if

      if (ncall < 0 .or. ncall == ng) close (ioUnit)

 999  return

c  Error handling:

 100  call ioExit (routine, 'Read ng returned error')
 102  call ioExit (routine, 'Read ibc_type et. al. returned error')
 103  call ioExit (routine, 'File does not exist')

 200  call ioExit (routine, 'Error in allocation of ibc_type')
 201  call ioExit (routine, 'Error in allocation of ibc_orient')
 202  call ioExit (routine, 'Error in allocation of rhand')

 500  format(I5,I3,3I2,5(1X,I3,3I2),2X,F2.0,3x,I2)

      end subroutine IO_R_BC_FLO107

c+------------------------------------------------------------------------------
c
      subroutine IO_R_FNC_P3D (q, ng, nfnc, ijk, ip, ncall, usr_file)
c
c  Description: IO_R_FNC_P3D reads a function file in PLOT3D format, but is
c  general enough to read files containing an arbitrary number of scalars.
c
c  The usr_file variable is an optional argument.  Absence, or the passing of an
c  unassociated variable in the actual argument, will cause an interactive user
c  query session to activate.
c
c  The ncall variable controls to what extent the file will be read.
c  Ncall = -1 results in the entire file, header and all BCs, being returned.
c  Ncall = 0  requests that only the header be read.
c  Ncall > 0  requires that the header has been read in a previous call and
c              requests that only the next block's function data be returned.
c
c  Note that the header of a PLOT3D function file consists of ng, nfnc, & ijk.
c
c  Arguments: See declaration for the complete type specification.
c
c     Intent(inout)
c       usr_file       : File description variable, optional except for
c                        ncall > 0 which requires the file to be previously
c                        opened with the ncall = 0 option
c
c     If ncall = -1    : Read whole file
c      Intent(out)
c       ng             : Number of grids
c       nfnc           : Number of scalars contained in file
c       ijk(ndim,ng)   : Computational dimensions of each grid
c       ip(ng+1)       : Index specifying beginning of mesh i in q array;
c                        ip(ng+1) is the total number of points + 1
c       q              : Packed scalar function data, for ng grids
c                      : q(nfnc,npts) => usr_file%triples = .T.
c                      : q(npts,nfnc) => usr_file%triples = .F.
c
c     If ncall = 0     : Read header of file
c      Intent(out)
c       ng             : Number of grids
c       nfnc           : Number of scalars contained in file
c       ijk(ndim,ng)   : Computational dimensions of each grid
c       ip             : Returned as an unassociated variable
c       q              : Returned as an unassociated variable
c
c     If ncall > 0     : Read next data block
c      Intent(in)
c       ng             : Number of grids
c       nfnc           : Number of scalars contained in file
c       ijk(ndim,ng)   : Computational dimensions of each grid
c
c      Intent(out)
c       ip(2)          : Degenerate case in which ip(1) = 1 and
c                        ip(2) = number of points in ncall_th block + 1
c       q              : Packed scalar function data for next grid
c                      : q(nfnc,npts) => usr_file%triples = .T.
c                      : q(npts,nfnc) => usr_file%triples = .F.
c
c  Side Effects:
c     Since the data pointers are disassociated and allocated during execution,
c     actual argument memory references prior to the call will be lost.
c
c     The file will be closed at the time of return when ncall = -1 or ng.
c     However, the file will remain open for ng = [0,ng-1].
c
c  Warnings:
c     Use of ncall >= 0, requires the presence of the usr_file argument.
c     This does not imply that the argument must be associated at the
c     time of the call when ncall = 0.
c
c     The value of ncall > 0 must correspond to the block number at the
c     current file pointer location.  The routine does not rewind/skip
c     ncall-1 times, but assumes that the current location is correct.
c
c  6/21/99  Updated ioPackage to allocate required memory for user using
c           the pointer attribute for the actual and dummy arguments.
c           Note: Cannot use intent and pointer in same attribute line!
c
c  7/1/99   Rewrote after adding ncall options.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      real, pointer    :: q(:,:)
      integer, pointer :: ijk(:,:), ip(:)
      integer          :: ng, nfnc, ncall
      type (ioMeshFileType), pointer, optional :: usr_file

c  Local constants:

      character, parameter :: routine * 12 = 'IO_R_FNC_P3D'

c  Local variables:

      integer ::
     .  i, ig, ig_s, ig_e, ig_tot, ioStat, ioUnit, j, ni, nj, ndim, npts
      logical ::
     .  ascii, exist, tuples
      type (ioMeshFileType), pointer ::
     .  prg_file

c  Execution:

c     Test if the programmer is using the optional file argument,
c     the unassociated file argument, or preinitialized options.

      call io_Qry_Str_X_File (usr_file, prg_file, ncall, routine)

c     Check the file attributes for errors, then open the file.

      call io_Prep_Open_X_File (prg_file, ioUnit, ncall, routine)

      if (ncall <= 0) then
        inquire (file=prg_file%name, exist=exist)
        if (.not. exist) goto 103
        open (unit=ioUnit, file=prg_file%name,
     .        form=prg_file%form, status=prg_file%status)
      end if

c     If requested, read the file header.

      ndim  = prg_file%ndim
      ascii = prg_file%form == FORMATTED

      if (ncall <= 0) then                       !Begin reading header
        if (prg_file%mg) then
          if (ascii) then
            read (ioUnit,*, err=100, iostat=ioStat) ng  !Read number of grids
          else
            read (ioUnit,   err=100, iostat=ioStat) ng
          end if
        else
          ng = 1
        end if
        if (ioStat /= 0) goto 100

        if (associated (ijk)) deallocate (ijk)
        allocate (ijk(ndim,ng), stat=ioStat)     !Read idim, jdim, kdim
        if (ioStat /= 0) goto 200

        if (ascii) then
          read (ioUnit, *, err=101, iostat=ioStat) ijk
        else
          read (ioUnit,    err=101, iostat=ioStat) ijk
        end if
        if (ioStat /= 0) goto 101

        if (ncall  == 0) goto 999                !Finished reading header
      end if

c     Set the loop counter to the number of blocks being read.

      if (ncall < 0) then
        ig_s = 1
        ig_e = ng
        ig_tot = ng
      else
        ig_s = ncall
        ig_e = ncall
        ig_tot = 1
      end if

c     Allocate ip array.

      if (associated (ip)) deallocate (ip)
      allocate (ip(ig_tot+1), stat=ioStat) !Size of the number of grids+1
      if (ioStat /= 0) goto 201

c     Set the ip index for the number of blocks being read.

      ip(1) = 1
      do ig = ig_s, ig_e
        npts = 1
        do j = 1, ndim
          npts = npts*ijk(j,ig)
        end do
        ip(ig-ig_s+2) = ip(ig-ig_s+1) + npts
      end do

c     Allocate q array.

      npts   = ip(ig_tot+1) - 1
      tuples = prg_file%triples

      if (tuples) then
        ni = nfnc
        nj = npts
      else
        ni = npts
        nj = nfnc
      end if

      if (associated (q)) deallocate (q)
      allocate (q(ni,nj), stat=ioStat)
      if (ioStat /= 0) goto 202

c     Read block data in the correct format.

      do ig = 1, ig_tot
        call IO_R_BLOCK (ig_tot, ig, ip, ni, nj, 0, q, ip,
     .                   tuples, ascii, ioUnit, ioStat)
        if (ioStat /= 0) goto 102
      end do

      if (ncall < 0 .or. ncall == ng) close (ioUnit)

 999  return

c  Error handling:

 100  call ioExit (routine, 'Read ng returned error')
 101  call ioExit (routine, 'Read ijk returned error')
 102  call ioExit (routine, 'Read q returned error')
 103  call ioExit (routine, 'File does not exist')

 200  call ioExit (routine, 'Error in allocation of ijk')
 201  call ioExit (routine, 'Error in allocation of ip')
 202  call ioExit (routine, 'Error in allocation of q')

      end subroutine IO_R_FNC_P3D

c+------------------------------------------------------------------------------
c
      subroutine IO_R_HYP_SYN107 (ihyp, nhyp, usr_file)
c
c  Note: This routine does not handle multiple hyperface groups.
c
c  Description: IO_R_HYP_SYN107 reads the hyperface specification from a
c  SYN107-MB connectivity file.  The boundary condition and wake data preceding
c  the hyperface information is skipped.  See IO_R_BC_FLO107 & IO_R_WAKE_FLO107
c  for extracting this information.
c
c  The usr_file variable is an optional argument.  Absence, or the passing of an
c  unassociated variable in the actual argument, will cause a predefined set of
c  user queries to activate.
c
c  Arguments: See declaration for the complete type specification.
c
c     Intent(inout)
c       usr_file       : File description variable (optional)
c
c      Intent(out)
c       nhyp           : Number of hyperfaces
c       ihyp(2,nhyp)   : Hyperface data where
c                      : ihyp(1,i) is the block number of ith hyperface
c                      : ihyp(2,i) is the face number of ith hyperface
c
c  Side Effects:
c     Since the data pointers are disassociated and allocated during execution,
c     memory references prior to the subroutine call will be lost.
c
c  7/7/99   Initial coding.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      integer, pointer :: ihyp(:,:)
      integer          :: nhyp
      type (ioMeshFileType), pointer, optional :: usr_file

c  Local constants:

      character, parameter :: routine * 15 = 'IO_R_HYP_SYN107'

c  Local variables:

      integer ::
     .  ig, ioStat, ioUnit, ng, ncall = -1
      logical ::
     .  exist
      type (ioMeshFileType), pointer ::
     .  prg_file

c  Execution:

c     Test if the programmer is using the optional file argument,
c     the unassociated file argument, or preinitialized options.

      call io_Qry_File (usr_file, prg_file, ncall, routine)

c     Check the file attributes for errors, then open the file:

      call io_Prep_Open_X_File (prg_file, ioUnit, ncall, routine)

      if (ncall <= 0) then
        inquire (file=prg_file%name, exist=exist)
        if (.not. exist) goto 100
        open (unit=ioUnit, file=prg_file%name,
     .        form=prg_file%form, status=prg_file%status)
      end if

c     Skip boundary condition and wake data to reach the hyperface section.

      if (prg_file%form == FORMATTED) then
        read (ioUnit, *, err=101, iostat=ioStat) ng    !BCs
        if (ioStat /= 0) goto 101
        do ig = 1, ng
          read (ioUnit, *, err=102, iostat=ioStat)
          if (ioStat /= 0) goto 102
        end do
        read (ioUnit, *, err=103, iostat=ioStat)       !Wakes
        read (ioUnit, *, err=104, iostat=ioStat) ng
        if (ioStat /= 0) goto 104
        do ig = 1,ng
          read (ioUnit, *, err=105, iostat=ioStat)
          if (ioStat /= 0) goto 105
        end do
        read (ioUnit, *, err=106, iostat=ioStat)       !# Hyperface regions
        if (ioStat /= 0) goto 106
        read (ioUnit, *, err=107, iostat=ioStat) ng
        if (ioStat /= 0) goto 107
        read (ioUnit, *, err=107, iostat=ioStat) nhyp
        if (ioStat /= 0) goto 107
      else                                             !Unformatted file
        read (ioUnit, err=101, iostat=ioStat) ng       !BCs
        if (ioStat /= 0) goto 101
        do ig = 1, ng
          read (ioUnit, err=102, iostat=ioStat)
          if (ioStat /= 0) goto 102
        end do
        read (ioUnit, err=103, iostat=ioStat)          !Wakes
        read (ioUnit, err=104, iostat=ioStat) ng
        do ig = 1,ng
          read (ioUnit, err=105, iostat=ioStat)
          if (ioStat /= 0) goto 105
        end do
        read (ioUnit, err=106, iostat=ioStat)          !# Hyperface regions
        if (ioStat /= 0) goto 106
        read (ioUnit, err=107, iostat=ioStat) ng
        if (ioStat /= 0) goto 107
        read (ioUnit, err=107, iostat=ioStat) nhyp
        if (ioStat /= 0) goto 107
      end if

c     Allocate workspace.

      if (associated (ihyp)) deallocate (ihyp)
      allocate (ihyp(2,nhyp), stat=ioStat)
      if (ioStat /= 0) goto 200

c     Read hyperface data

      if (prg_file%form == FORMATTED) then
        do ig = 1, nhyp
          read (ioUnit, *, err=108, iostat=ioStat) ihyp(1:2,ig)
          if (ioStat /= 0) goto 108
        end do
      else
        do ig = 1, nhyp
          read (ioUnit, err=108, iostat=ioStat) ihyp
          if (ioStat /= 0) goto 108
        end do
      end if

      close (ioUnit)

      return

c  Error handling:

 100  call ioExit (routine, 'File does not exist')
 101  call ioExit (routine, 'Error reading ng')
 102  call ioExit (routine, 'Error skipping bc data')
 103  call ioExit (routine, 'Error skipping wake title')
 104  call ioExit (routine, 'Error skipping number of wakes')
 105  call ioExit (routine, 'Error skipping a wake line')
 106  call ioExit (routine, 'Error skipping a hyper title line')
 107  call ioExit (routine, 'Error reading number of hyperfaces')
 108  call ioExit (routine, 'Error reading a hyperface line')

 200  call ioExit (routine, 'Error in allocation of ihyp')

      end subroutine IO_R_HYP_SYN107

c+------------------------------------------------------------------------------
c
      subroutine IO_R_DSNF_SYN107 (ndf, srcBlock, srcFace, targBlock,
     .                             targFace, usr_file)
c
c  Description: Read the design face info of a SYN107-MB connectivity file.
c
c  7/7/99   Initial coding.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      integer, pointer :: srcBlock(:), srcFace(:),
     .                    targBlock(:), targFace(:)
      integer          :: ndf
      type (ioMeshFileType), pointer, optional :: usr_file

c  Local constants:

      character, parameter :: routine * 16 = 'IO_R_DSNF_SYN107'

c  Local variables:

      integer ::
     .  ioStat, ioUnit, ig, l, ng, nhyp, ncall = -1
      logical ::
     .  exist
      type (ioMeshFileType), pointer ::
     .  prg_file

c  Execution:

c     Test if the programmer is using the optional file argument,
c     the unassociated file argument, or preinitialized options.

      call io_Qry_File (usr_file, prg_file, ncall, routine)

c     Check the file attributes for errors, then open the file:

      call io_Prep_Open_X_File (prg_file, ioUnit, ncall, routine)

      if (ncall <= 0) then
        inquire (file=prg_file%name, exist=exist)
        if (.not. exist) goto 100
        open (unit=ioUnit, file=prg_file%name,
     .        form=prg_file%form, status=prg_file%status)
      end if

c     Skip boundary condition and wake data to reach the hyperface section.

      if (prg_file%form == FORMATTED) then
        read (ioUnit, *, err=101, iostat=ioStat) ng    !BCs
        if (ioStat /= 0) goto 101
        do ig = 1, ng
          read (ioUnit, *, err=102, iostat=ioStat)
          if (ioStat /= 0) goto 102
        end do
        read (ioUnit, *, err=103, iostat=ioStat)       !Wakes
        read (ioUnit, *, err=104, iostat=ioStat) ng
        if (ioStat /= 0) goto 104
        do ig = 1,ng
          read (ioUnit, *, err=105, iostat=ioStat)
          if (ioStat /= 0) goto 105
        end do
        read (ioUnit, *, err=106, iostat=ioStat)       !Number of Hyperfaces
        if (ioStat /= 0) goto 106
        read (ioUnit, *, err=107, iostat=ioStat) nhyp
        if (ioStat /= 0) goto 107
        read (ioUnit, *, err=107, iostat=ioStat) ng
        if (ioStat /= 0) goto 107
        do ig = 1,ng
          read (ioUnit, *, err=108, iostat=ioStat)
        end do
        read (ioUnit, *, err=109, iostat=ioStat)
        if (ioStat /= 0) goto 109
        read (ioUnit, *, err=110, iostat=ioStat) ndf
        if (ioStat /= 0) goto 110

      else                                             !Unformatted File
        read (ioUnit, err=101, iostat=ioStat) ng       !BCs
        if (ioStat /= 0) goto 101
        do ig = 1, ng
          read (ioUnit, err=102, iostat=ioStat)
          if (ioStat /= 0) goto 102
        end do
        read (ioUnit, err=103, iostat=ioStat)          !Wakes
        read (ioUnit, err=104, iostat=ioStat) ng
        do ig = 1,ng
          read (ioUnit, err=105, iostat=ioStat)
          if (ioStat /= 0) goto 105
        end do

        read (ioUnit, err=106, iostat=ioStat)          !Number of hyperfaces
        if (ioStat /= 0) goto 106
        read (ioUnit, err=107, iostat=ioStat) nhyp
        if (ioStat /= 0) goto 107
        read (ioUnit, err=107, iostat=ioStat) ng
        if (ioStat /= 0) goto 107
        do ig = 1,ng
          read (ioUnit,err=108,iostat=ioStat)
        end do

        read (ioUnit, err=109, iostat=ioStat)
        if (ioStat /= 0) goto 109
        read (ioUnit, err=110, iostat=ioStat) ndf
        if (ioStat /= 0) goto 110

      end if

c     Allocate workspace:

      if (associated (srcBlock))  deallocate (srcBlock)
      if (associated (srcFace))   deallocate (srcFace)
      if (associated (targBlock)) deallocate (targBlock)
      if (associated (targFace))  deallocate (targFace)

      allocate (srcBlock(ndf), stat=ioStat)
      if (ioStat /= 0) goto 200
      allocate (targBlock(ndf), stat=ioStat)
      if (ioStat /= 0) goto 201
      allocate (srcFace(ndf), stat=ioStat)
      if (ioStat /= 0) goto 202
      allocate (targFace(ndf), stat=ioStat)
      if (ioStat /= 0) goto 203

c     Read hyperface data

      if (prg_file%form == FORMATTED) then
        do ig = 1, ndf
          read (ioUnit, *, err=111, iostat=ioStat)
     .      srcBlock(ig), srcFace(ig), targBlock(ig), targFace(ig)
          if (ioStat /= 0) goto 111
        end do
      else                                    !Unformatted File
        do ig = 1, ndf
          read (ioUnit, err=111, iostat=ioStat)
     .      srcBlock(ig), srcFace(ig), targBlock(ig), targFace(ig)
          if (ioStat /= 0) goto 111
        end do
      end if

      close (ioUnit)

      return

c  Error handling:

 100  call ioExit (routine, 'File does not exist')
 101  call ioExit (routine, 'Error reading ng')
 102  call ioExit (routine, 'Error skipping bc data')
 103  call ioExit (routine, 'Error skipping wake title')
 104  call ioExit (routine, 'Error skipping number of wakes')
 105  call ioExit (routine, 'Error skipping a wake line')
 106  call ioExit (routine, 'Error skipping a hyper title line')
 107  call ioExit (routine, 'Error reading nhyper')
 108  call ioExit (routine, 'Error skipping a hyperface')
 109  call ioExit (routine, 'Error reading a design title line')
 110  call ioExit (routine, 'Error reading number of design faces')
 111  call ioExit (routine, 'Error reading a design face')

 200  call ioExit (routine, 'Error in allocation of srcBlock')
 201  call ioExit (routine, 'Error in allocation of targBlock')
 202  call ioExit (routine, 'Error in allocation of srcFace')
 203  call ioExit (routine, 'Error in allocation of targFace')

      end subroutine IO_R_DSNF_SYN107

c+------------------------------------------------------------------------------
c
      subroutine IO_R_Q_P3D (q, ng, nq, ijk, ip,
     .                       fsmach, alpha, re, time, ncall, usr_file)
c
c  Description: IO_R_Q_P3D reads a solution file in PLOT3D format.
c  PLOT3D format stipulates that the solution file contain only pgr_file%ndim+2
c  variables, i.e., no species or turbulent quantities are expected.  A routine
c  for reading an OVERFLOW solution file may be implemented as IO_R_Q_OVERFLOW.
c
c  The usr_file variable is an optional argument.  Absence, or the passing
c  of an unassociated variable in the actual argument, will cause a predefined
c  set of user queries to activate.
c
c  The ncall variable controls to what extent the file will be read.
c  Ncall = -1 results in the entire file, header and all BCs, being returned.
c  Ncall = 0  requests that only the header be read.
c  Ncall > 0  requires that the header has been read in a previous call and
c              requests that only the next block's q-data be returned.
c
c  Note that the header of a PLOT3D solution file consists of ng, nq, and ijk.
c
c  Arguments: See declaration for the complete type specification.
c
c     Intent(inout)
c       usr_file       : File description variable, optional except for
c                        ncall > 0 which requires the file to be previously
c                        opened with the ncall = 0 option
c
c     If ncall = -1    : Read whole file
c      Intent(out)
c       ng             : Number of grids
c       nq             : Assumed to be ndim+2
c       ijk(ndim,ng)   : Computational dimensions of each grid
c       ip(ng+1)       : Index specifying beginning of mesh i in q array;
c                        ip(ng+1) is the total number of points + 1
c       q              : Packed scalar function data, for ng grids
c                      : q(nq,npts) => prg_file%triples = .T.
c                      : q(npts,nq) => prg_file%triples = .F.
c       fsmach, alpha  : Scalar quantities assumed to be the same
c       re, time       : in each q data set
c
c     If ncall = 0     : Read header of file
c      Intent(out)
c       ng             : Number of grids
c       nq             : Assumed to be ndim+2
c       ijk(ndim,ng)   : Computational dimensions of each grid
c       ip             : Returned as an unassociated variable
c      Intent(none)
c       q              : Returned as an unassociated variable
c       fsmach, alpha  : Returned as unitialized data
c       re, time       :
c
c     If ncall > 0     : Read next data block
c      Intent(in)
c       ng             : Number of grids
c       nq             : Assumed to be ndim+2
c       ijk(ndim,ng)   : Computational dimensions of each grid
c
c      Intent(out)
c       ip(2)          : Degenerate case in which ip(1) = 1 and
c                        ip(2) = (number of points in ncall_th block) + 1
c       q              : Packed scalar q data, for next grid
c                      : q(nq,npts) => prg_file%triples = .T.
c                      : q(npts,nq) => prg_file%triples = .F.
c       fsmach, alpha  : Scalar quantities assumed to be the same
c       re, time       : in each q data set
c
c  Side Effects:
c     Since the data pointers are disassociated and allocated during execution,
c     actual argument memory references prior to the call will be lost.
c
c     The file will be closed at the time of return when ncall = -1 or ng.
c     However, the file will remain open for ng = [0,ng-1].
c
c  Warnings:
c     Use of ncall >= 0, requires the presence of the usr_file argument.
c     This does not imply that the argument must be associated at the
c     time of the call when ncall = 0.
c
c     The value of ncall > 0 must correspond to the block number at the
c     current file pointer location.  The routine does not rewind/skip
c     ncall-1 times, but assumes that the current location is correct.
c
c  7/8/99   Rewrote after adding ncall options.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      real, pointer    :: q(:,:)
      integer, pointer :: ijk(:,:), ip(:)
      integer          :: ng, nq, ncall
      real             :: fsmach, alpha, re, time
      type (ioMeshFileType), pointer, optional :: usr_file

c  Local constants:

      character, parameter :: routine * 10 = 'IO_R_Q_P3D'

c  Local variables:

      integer ::
     .  i, ig, ig_s, ig_e, ig_tot, ioStat, ioUnit, j, ni, nj, ndim, npts
      logical ::
     .  ascii, exist, tuples
      type (ioMeshFileType), pointer ::
     .  prg_file

c  Execution:

c     Test if the programmer is using the optional file argument,
c     the unassociated file argument, or preinitialized options.

      call io_Qry_File (usr_file, prg_file, ncall, routine)

c     Check the file attributes for errors, then open the file:

      call io_Prep_Open_X_File (prg_file, ioUnit, ncall, routine)

      if (ncall <= 0) then
        inquire (file=prg_file%name, exist=exist)
        if (.not. exist) goto 103
        open (unit=ioUnit, file=prg_file%name,
     .        form=prg_file%form, status=prg_file%status)
      end if

c     Set the spatial dimension and length of vector in q-file.
c     This routine assumes a PLOT3D q-file without turbulence and species data.

      ndim  = prg_file%ndim
      ascii = prg_file%form == FORMATTED
      nq    = ndim + 2

      if (ncall <= 0) then          !Begin reading header
        if (prg_file%mg) then
          if (ascii) then
            read (ioUnit, *, err=100, iostat=ioStat) ng  !Number of grids
          else
            read (ioUnit,    err=100, iostat=ioStat) ng
          end if
          if (ioStat /= 0) goto 100
        else
          ng = 1
        end if

        if (associated (ijk)) deallocate (ijk)  !Allocate and read
        allocate (ijk(ndim,ng), stat=ioStat)    !idim, jdim, kdim
        if (ioStat /= 0) goto 200

        if (ascii) then
          read (ioUnit, *, err=101, iostat=ioStat) ijk
        else
          read (ioUnit,    err=101, iostat=ioStat) ijk
        end if
        if (ioStat /= 0) goto 101

        if (ncall == 0) goto 999               !Finished reading header
      end if

c     Set the loop counter to the number of blocks being read.

      if (ncall < 0) then
        ig_s = 1
        ig_e = ng
        ig_tot = ng
      else
        ig_s = ncall
        ig_e = ncall
        ig_tot = 1
      end if

c     Allocate ip array.

      if (associated (ip)) deallocate (ip)
      allocate (ip(ig_tot+1), stat=ioStat)
      if (ioStat /= 0) goto 201

c     Set the ip index for the number of blocks being read.
      ip(1) = 1
      do ig = ig_s, ig_e
        npts = 1
        do j = 1,ndim
          npts = npts*ijk(j,ig)
        end do
        ip(ig-ig_s+2) = ip(ig-ig_s+1) + npts
      end do

c     Allocate q array.

      npts   = ip(ig_tot+1) - 1
      tuples = prg_file%triples

      if (tuples) then
        ni = nq
        nj = npts
      else
        ni = npts
        nj = nq
      end if

      if (associated (q)) deallocate (q)
      allocate (q(ni,nj), stat=ioStat)
      if (ioStat /= 0) goto 202

c     Read block data in correct format.

      do ig = 1, ig_tot

        if (ascii) then
          read (ioUnit, *, err=104, iostat=ioStat)
     .      fsmach, alpha, re, time
        else
          read (ioUnit,    err=104, iostat=ioStat)
     .      fsmach, alpha, re, time
        end if
        if (ioStat /= 0) goto 104

        call IO_R_BLOCK (ig_tot, ig, ip, ni, nj, 0, q, ip,
     .                   tuples, ascii, ioUnit, ioStat)
        if (ioStat /= 0) goto 102

      end do

      if (ncall < 0 .or. ncall == ng) close (ioUnit)

 999  return

c  Error handling:

 100  call ioExit (routine, 'Read ng returned error')
 101  call ioExit (routine, 'Read ijk returned error')
 102  call ioExit (routine, 'Read q returned error')
 103  call ioExit (routine, 'File does not exist')
 104  call ioExit (routine, 'Read fsmach returned error')

 200  call ioExit (routine, 'Error in allocation of ijk')
 201  call ioExit (routine, 'Error in allocation of ip')
 202  call ioExit (routine, 'Error in allocation of q')

      end subroutine IO_R_Q_P3D

c+------------------------------------------------------------------------------
c
      subroutine IO_R_STR_X_P3D (x, ng, ijk, ip, ncall, usr_file)
c
c  Description: IO_R_STR_X_P3D reads an unblanked, structured mesh file in
c  PLOT3D format.  Refer to IO_R_STR_X_P3D_IB for the iblanked mesh format.
c
c  The usr_file variable is an optional argument.  Absence, or the passing of an
c  unassociated variable in the actual argument, will cause a predefined set of
c  user queries to activate.
c
c  The ncall variable controls to what extent the file will be read.
c  Ncall = -1 results in the entire file, header and all data, being returned.
c  Ncall = 0  requests that only the header be read.
c  Ncall > 0  requires that the header has been read in a previous call and
c             requests that only the next block's x data be returned.
c
c  Note that the header of a PLOT3D solution file consists of ng and ijk.
c
c  Arguments: See declaration for the complete type specification.
c
c     Intent(inout)
c       usr_file       : File description variable, optional except for
c                        ncall > 0 which requires the file to be previously
c                        opened with the ncall = 0 option.
c
c     If ncall = -1    : Read whole file
c      Intent(out)
c       ng             : Number of grids.
c       ijk(ndim,ng)   : Computational dimensions of each grid.
c       ip(ng+1)       : Index specifying beginning of mesh i in x array;
c                        ip(ng+1) is the total number of points + 1.
c       x              : Packed mesh data for ng grids.
c                      : x(ndim,npts) => usr_file%triples = .T.
c                      : x(npts,ndim) => usr_file%triples = .F.
c
c     If ncall = 0     : Read header of file.
c      Intent(out)
c       ng             : Number of grids.
c       ijk(ndim,ng)   : Computational dimensions of each grid.
c       ip             : Returned as an unassociated variable.
c       x              : Returned as an unassociated variable.
c
c     If ncall > 0     : Read next data block.
c      Intent(in)
c       ng             : Number of grids.
c       ijk(ndim,ng)   : Computational dimensions of each grid.
c
c      Intent(out)
c       ip(2)          : Degenerate case in which ip(1) = 1 and
c                        ip(2) = (npts in ncall_th block) + 1.
c       x              : Packed mesh data for next grid.
c                      : x(nq,npts) => usr_file%triples = .T.
c                      : x(npts,nq) => usr_file%triples = .F.
c
c  Side Effects:
c     Since the data pointers are disassociated and allocated during execution,
c     actual argument memory references prior to the call will be lost.
c
c     The file will be closed at the time of return when ncall = -1 or ng.
c     However, the file will remain open for ng = [0,ng-1].
c
c  Warnings:
c     Use of ncall >= 0, requires the presence of the usr_file argument.
c     This does not imply that the argument must be associated at the
c     time of the call when ncall = 0.
c
c     The value of ncall > 0 must correspond to the block number at the
c     current file pointer location.  The routine does not rewind/skip
c     ncall-1 times, but assumes that the current location is correct.
c
c  06/21/99   Initial coding with auto-memory handling.
c  07/01/99   Rewrote to add ncall option.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      real, pointer    :: x(:,:)
      integer, pointer :: ijk(:,:), ip(:)
      integer          :: ng, ncall
      type (ioMeshFileType), pointer, optional :: usr_file

c  Local constants:

      character, parameter :: routine * 14 = 'IO_R_STR_X_P3D'

c  Local variables:

      integer ::
     .  i, ioStat, ioUnit, ig, ig_s, ig_e, ig_tot, j, ndim, ni, nj, npts
      logical ::
     .  ascii, exist, tuples
      type (ioMeshFileType), pointer ::
     .  prg_file

c  Execution:

c     Test if the programmer is using the optional file argument,
c     the unassociated file argument, or pre-initialized options.

      call io_Qry_Str_X_File (usr_file, prg_file, ncall, routine)

c     Check the file attributes for errors, then open the file:

      call io_Prep_Open_X_File (prg_file, ioUnit, ncall, routine)

      if (ncall <= 0) then
        inquire (file=prg_file%name, exist=exist)
        if (.not. exist) goto 103
        open (unit=ioUnit, file=prg_file%name,
     .    form=prg_file%form, status=prg_file%status)
      end if

c     If requested, read file header.

      ndim  = prg_file%ndim
      ascii = prg_file%form == FORMATTED

      if (ncall <= 0) then
        if (prg_file%mg) then                 !Begin reading header
          if (ascii) then
            read (ioUnit, *, err=100, iostat=ioStat) ng  !Number of grids
          else
            read (ioUnit, err=100, iostat=ioStat) ng
          end if
          if (ioStat /= 0) goto 100
        else
          ng = 1
        end if

        if (associated (ijk)) deallocate (ijk)  !Allocate ijk and
        allocate (ijk(ndim,ng), stat=ioStat)    !read dimensions
        if (ioStat /= 0) goto 200

        if (ascii) then
          read (ioUnit, *, err=101, iostat=ioStat) ijk
        else
          read (ioUnit, err=101, iostat=ioStat) ijk
        end if
        if (ioStat /= 0) goto 101

        if (ncall == 0) goto 999              !Finished reading header
      end if

c     Set loop counters to the number of blocks being read.

      if (ncall < 0) then
        ig_s = 1
        ig_e = ng
        ig_tot = ng
      else
        ig_s = ncall
        ig_e = ncall
        ig_tot = 1
      end if

c     Allocate ip array.

      if (associated (ip)) deallocate (ip)
      allocate (ip(ig_tot+1), stat=ioStat)
      if (ioStat /= 0) goto 201

c     Set the ip index for the number of blocks being read.

      ip(1) = 1
      do ig = ig_s, ig_e
        npts = 1
        do j = 1, ndim
          npts = npts*ijk(j,ig)
        end do
        ip(ig-ig_s+2) = ip(ig-ig_s+1) + npts
      end do

c     Allocate x array.

      npts   = ip(ig_tot+1) - 1
      tuples = prg_file%triples

      if (tuples) then
        ni = ndim
        nj = npts
      else
        ni = npts
        nj = ndim
      end if

      if (associated (x)) deallocate (x)
      allocate (x(ni,nj), stat=ioStat)
      if (ioStat /= 0) goto 202

c     Read the data as packed blocks.

      do ig = 1, ig_tot
        call IO_R_BLOCK (ig_tot, ig, ip, ni, nj, 0, x, ip,
     .                   tuples, ascii, ioUnit, ioStat)
        if (ioStat /= 0) goto 102
      end do

      if (ncall < 0 .or. ncall == ng) close (ioUnit)

 999  return

c  Error handling:

 100  call ioExit  (routine, 'Read ng returned error')
 101  call ioExit  (routine, 'Read ijk returned error')
 102  call ioExit  (routine, 'Read x returned error')
 103  call ioExit  (routine, 'File does not exist')

 200  call ioExit  (routine, 'Error in allocation of ijk')
 201  call ioExit  (routine, 'Error in allocation of ip')
 202  call ioExit  (routine, 'Error in allocation of x')

      end subroutine IO_R_STR_X_P3D

c+------------------------------------------------------------------------------
c
      subroutine IO_R_STR_X_P3D_IB (x, ng, ijk, ip, iblank, ncall,
     .                              usr_file)
c
c  Description: IO_R_STR_X_P3D_IB reads a blanked, structured mesh file in
c  PLOT3D format.  Refer to IO_R_STR_X_P3D for the blanked mesh file format.
c
c  The usr_file variable is an optional argument.  Absence, or the passing of an
c  unassociated variable in the actual argument, will cause a predefined set of
c  user queries to activate.
c
c  The ncall variable controls to what extent the file will be read.
c  Ncall = -1 results in the entire file, header and all data, being returned.
c  Ncall = 0  requests that only the header be read.
c  Ncall > 0  requires that the header has been read in a previous call and
c             requests that only the next block's x & iblank data be returned.
c
c  Note that the header of a PLOT3D solution file consists of ng and ijk.
c
c  Arguments: See declaration for the complete type specification.
c
c     Intent(inout)
c       usr_file       : File description variable, optional except for
c                        ncall > 0 which requires the file to be previously
c                        opened with the ncall = 0 option
c
c     If ncall = -1    : Read whole file
c      Intent(out)
c       ng             : Number of grids
c       ijk(ndim,ng)   : Computational dimensions of each grid
c       ip(ng+1)       : Index specifying beginning of mesh i in x array;
c                        ip(ng+1) is the total number of points + 1
c       x              : Packed mesh data for ng grids
c                      : x(ndim,npts) => usr_file%triples = .T.
c                      : x(npts,ndim) => usr_file%triples = .F.
c       iblank(npts)   : Packed blanking data for each grid
c
c     If ncall = 0     : Read header of file
c      Intent(out)
c       ng             : Number of grids
c       ijk(ndim,ng)   : Computational dimensions of each grid
c       ip             : Returned as an unassociated variable
c       x              : Returned as an unassociated variable
c       iblank         : Returned as an unassociated variable
c
c     If ncall > 0     : Read next data block
c      Intent(in)
c       ng             : Number of grids
c       ijk(ndim,ng)   : Computational dimensions of each grid
c
c      Intent(out)
c       ip(2)          : Degenerate case in which ip(1) = 1 and
c                        ip(2) = (npts in ncall_th block) + 1
c       x              : Packed mesh data for next grid
c                      : x(nq,npts) => usr_file%triples = .T.
c                      : x(npts,nq) => usr_file%triples = .F.
c       iblank(npts)   : Packed blanking data for next grid
c
c  Side Effects:
c     Since the data pointers are disassociated and allocated during execution,
c     actual argument memory references prior to the call will be lost.
c
c     The file will be closed at the time of return when ncall = -1 or ng.
c     However, the file will remain open for ng = [0,ng-1].
c
c  Warnings:
c     Use of ncall >= 0, requires the presence of the usr_file argument.
c     This does not imply that the argument must be associated at the
c     time of the call when ncall = 0.
c
c     The value of ncall > 0 must correspond to the block number at the
c     current file pointer location.  The routine does not rewind/skip
c     ncall-1 times, but assumes that the current location is correct.
c
c  06/21/99   Initial coding with auto-memory handling.
c  07/01/99   Rewrote to add ncall option.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      real, pointer    :: x(:,:)
      integer, pointer :: ijk(:,:), ip(:), iblank(:)
      integer          :: ng, ncall
      type (ioMeshFileType), pointer, optional :: usr_file

c  Local constants:

      character, parameter :: routine * 17 = 'IO_R_STR_X_P3D_IB'

c  Local variables:

      integer ::
     .  i, ig, ig_s, ig_e, ig_tot, ioStat, ioUnit, j, ndim, ni, nj, npts
      logical ::
     .  ascii, exist, tuples
      type (ioMeshFileType), pointer ::
     .  prg_file

c  Execution:

c     Test if the programmer is using the optional file argument,
c     the unassociated file argument, or preinitialized options.

      call io_Qry_Str_X_File (usr_file, prg_file, ncall, routine)

c     Check the file attributes for errors, then open the file:

      call io_Prep_Open_X_File (prg_file, ioUnit, ncall, routine)

      if (ncall <= 0) then
        inquire (file=prg_file%name, exist=exist)
        if (.not. exist) goto 103
        open (unit=ioUnit, file=prg_file%name,
     .        form=prg_file%form, status=prg_file%status)
      end if

c     If requested, read file header.

      ndim  = prg_file%ndim
      ascii = prg_file%form == FORMATTED

      if (ncall <= 0) then
        if (prg_file%mg) then                 !Begin reading header
          if (ascii) then
            read (ioUnit, *, err=100, iostat=ioStat) ng  !Number of grids
          else
            read (ioUnit,    err=100, iostat=ioStat) ng
          end if
          if (ioStat /= 0) goto 100
        else
          ng = 1
        end if

        if (associated (ijk)) deallocate (ijk)  !Allocate ijk and
        allocate (ijk(ndim,ng), stat=ioStat)    !read dimensions
        if (ioStat /= 0) goto 200

        if (ascii) then
          read (ioUnit, *, err=101, iostat=ioStat) ijk
        else
          read (ioUnit,    err=101, iostat=ioStat) ijk
        end if
        if (ioStat /= 0) goto 101

        if (ncall == 0) goto 999              !Finished reading header
      end if

c     Set loop counters to the number of blocks being read.

      if (ncall < 0) then
        ig_s = 1
        ig_e = ng
        ig_tot = ng
      else
        ig_s = ncall
        ig_e = ncall
        ig_tot = 1
      end if

c     Allocate ip array.

      if (associated (ip)) deallocate (ip)
      allocate (ip(ig_tot+1), stat=ioStat)
      if (ioStat /= 0) goto 201

c     Set the ip index for the number of blocks being read.
        ip(1) = 1
        do ig = ig_s, ig_e
          npts = 1
          do j = 1, ndim
            npts = npts*ijk(j,ig)
          end do
          ip(ig-ig_s+2) = ip(ig-ig_s+1) + npts
        end do

c     Allocate x and iblank arrays.

      npts   = ip(ig_tot+1) - 1
      tuples = prg_file%triples

      if (tuples) then
        ni = ndim
        nj = npts
      else
        ni = npts
        nj = ndim
      end if

      if (associated (x)) deallocate (x)
      allocate (x(ni,nj), stat=ioStat)
      if (ioStat /= 0) goto 202

      if (associated (iblank)) deallocate (iblank)
      allocate (iblank(npts), stat=ioStat)
      if (ioStat /= 0) goto 203

c     Read block data in correct format.

      do ig = 1, ig_tot
        call IO_R_BLOCK (ig_tot, ig, ip, ni, nj, npts, x, iblank,
     .                   tuples, ascii, ioUnit, ioStat)
        if (ioStat /= 0) goto 102
      end do

      if (ncall < 0 .or. ncall == ng) close (ioUnit)

 999  return

c  Error handling:

 100  call ioExit (routine, 'Read ng returned error')
 101  call ioExit (routine, 'Read ijk returned error')
 102  call ioExit (routine, 'Read x returned error')
 103  call ioExit (routine, 'File does not exist')

 200  call ioExit (routine, 'Error in allocation of ijk')
 201  call ioExit (routine, 'Error in allocation of ip')
 202  call ioExit (routine, 'Error in allocation of x')
 203  call ioExit (routine, 'Error in allocation of iblank')

      end subroutine IO_R_STR_X_P3D_IB

c+------------------------------------------------------------------------------
c
      subroutine IO_R_UNSTR_X_FAST (x, nPoints, triVerts, iTriVerts,
     .                              nTris, tetVerts, nTets, ncall,
     .                              usr_file)
c
c  Description: IO_R_UNSTR_X_FAST reads a single-grid unstructured mesh in
c  FAST format.
c
c  The usr_file variable is an optional argument.  Absence, or the passing of an
c  unassociated variable in the actual argument, will cause a predefined set of
c  user queries to activate.
c
c  The ncall variable controls to what extent the file will be read.
c  Ncall = -1 results in the entire file, header and all BCs, being returned.
c  Ncall = 0  requests that only the header be read.
c  Ncall > 0  requires that the header has been read in a previous call and
c             requests that the data block be read.
c
c  Note that the header of a FAST unstructured mesh file consists of
c  nPoints, nTriangles, nTets.
c
c  Arguments: See declaration for the complete type specification.  Note that
c  an unstructured mesh is required to be 3-D.
c
c     Intent(inout)
c       usr_file       : File description variable, optional except for
c                        ncall > 0 which requires the file to be previously
c                        opened with the ncall = 0 option
c
c     If ncall = -1    : Read whole file
c      Intent(out)
c       nPoints        : Number of (x,y,z) points in mesh
c       nTri           : Number of triangles in mesh
c       nTet           : Number of tetrahedra in mesh
c
c       x              : Packed (x,y,z) data for ng grids
c                      : x(ndim,npts) => usr_file%triples = .T.
c                      : x(npts,ndim) => usr_file%triples = .F.
c       triVerts       : List of triangles in mesh
c                      : triVerts(3,nTri)
c       iTriVerts      : Integer flag associated with each triangle
c                      : iTriVerts(nTri)
c       tetVerts       : List of tetrahedra in mesh
c                      : tetVerts(4,nTet)
c
c     If ncall = 0     : Read header of file
c      Intent(out)
c       nPoints        : Number of (x,y,z) points in mesh
c       nTri           : Number of triangles in mesh
c       nTet           : Number of tetrahedra in mesh
c
c       x              : Returned as an unassociated variable
c       triVerts       : Returned as an unassociated variable
c       iTriVerts      : Returned as an unassociated variable
c       tetVerts       : Returned as an unassociated variable
c
c     If ncall > 0     : Read next data block
c      Intent(in)
c       nPoints        : Number of (x,y,z) points in mesh
c       nTri           : Number of triangles in mesh
c       nTet           : Number of tetrahedra in mesh
c
c      Intent(out)
c       x              : Packed (x,y,z) data for ng grids
c                      : x(ndim,npts) => usr_file%triples = .T.
c                      : x(npts,ndim) => usr_file%triples = .F.
c       triVerts       : List of triangles in mesh
c                      : triVerts(3,nTri)
c       iTriVerts      : Integer flag associated with each triangle
c                      : iTriVerts(nTri)
c       tetVerts       : List of tetrahedra in mesh
c                      : tetVerts(4,nTet)
c
c  Side Effects:
c     Since the data pointers are disassociated and allocated during execution,
c     actual argument memory references prior to the call will be lost.
c
c     The file will be closed at the time of return when ncall = -1 or
c     ncall > 0, but will remain open when ncall = 0.
c
c  Warnings:
c     Use of ncall >= 0, requires the presence of the usr_file argument.
c     This does not imply that the argument must be associated at the
c     time of the call when ncall = 0.
c
c     The routine assumes a single unstructured mesh.
c
c  08/08/99  MJR  Initial coding with auto-memory handling.
c  11/19/99  DAS  Avoid implied loops in the I/O.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      real, pointer    :: x(:,:)
      integer, pointer :: triVerts(:,:), tetVerts(:,:), iTriVerts(:)
      integer          :: nPoints, nTris, nTets
      integer          :: ncall
      type (ioMeshFileType), pointer, optional :: usr_file

c  Local constants:

      integer,   parameter :: THREE = 3, FOUR = 4
      character, parameter :: routine * 17 = 'IO_R_UNSTR_X_FAST'

c  Local variables:

      integer ::
     .  i, ioStat, ioUnit, ip(2), j, nblank, ndim, ni, nj
      integer, allocatable ::
     .  iblank(:)
      logical ::
     .  ascii, exist, tuples
      type (ioMeshFileType), pointer ::
     .  prg_file

c  Execution:

c     Test if the programmer is using the optional file argument,
c     the unassociated file argument, or preinitialized options.

      call io_Qry_Unstr_X_File (usr_file, prg_file, ncall, routine)

c     Check the file attributes for errors, then open the file:

      call io_Prep_Open_X_File (prg_file, ioUnit, ncall, routine)

      if (ncall <= 0) then
        inquire (file=prg_file%name, exist=exist)
        if (.not. exist) goto 102
        open (unit=ioUnit, file=prg_file%name,
     .        form=prg_file%form, status=prg_file%status)
      end if

c     If requested, read file header.

      ndim  = prg_file%ndim
      ascii = prg_file%form == FORMATTED

      if (ncall <= 0) then                    !Begin reading header
        if (ascii) then
          read (ioUnit, *, err=100, iostat=ioStat)
     .      nPoints, nTris, nTets
        else
          read (ioUnit,    err=100, iostat=ioStat)
     .      nPoints, nTris, nTets
        end if
        if (ioStat /= 0) goto 100

        if (ncall == 0) goto 999              !Finished reading header
      end if

c     Allocate x.

      tuples = prg_file%triples

      if (tuples) then
        ni = ndim
        nj = nPoints
      else
        ni = nPoints
        nj = ndim
      end if

      if (associated (x)) deallocate (x)
      allocate (x(ni,nj), stat=ioStat)
      if (ioStat /= 0) goto 200

c     Allocate vertex info storage in a form suited to IO_R_BLOCK.

      nblank = FOUR * (nTris + nTets)
      allocate (iblank(nblank))

c     Read grid and pointers in correct format.

      ip(1) = 1
      ip(2) = nPoints + 1

      call IO_R_BLOCK (1, 1, ip, ni, nj, nblank, x, iblank,
     .                 tuples, ascii, ioUnit, ioStat)
      if (ioStat /= 0) goto 101

c     Transfer iblank(*) to triVerts, iTriVerts, tetVerts.

      if (associated (triVerts)) deallocate (triVerts)
      allocate (triVerts(THREE,nTris), stat=ioStat)
      if (ioStat /= 0) goto 201

      i = 1
      do j = 1, nTris
         triVerts(1:THREE,j) = iblank(i:i+2)
         i = i + THREE
      end do

      if (associated (itriVerts)) deallocate (itriVerts)
      allocate (itriVerts(nTris), stat=ioStat)
      if (ioStat /= 0) goto 202

      j = i + nTris
      itriVerts = iblank(i:j-1)

      if (associated (tetVerts)) deallocate (tetVerts)
      allocate (tetVerts(FOUR,nTets), stat=ioStat)
      if (ioStat /= 0) goto 203

      i = j
      do j = 1, nTets
         tetVerts(1:FOUR,j) = iblank(i:i+THREE)
         i = i + FOUR
      end do

      deallocate (iblank)

      close (ioUnit)

 999  return

c  Error handling:

 100  call ioExit (routine, 'Read nPoints returned error')
 101  call ioExit (routine, 'Read x returned error')
 102  call ioExit (routine, 'File does not exist')

 200  call ioExit (routine, 'Error in allocation of x')
 201  call ioExit (routine, 'Error in allocation of triVerts')
 202  call ioExit (routine, 'Error in allocation of iTriVerts')
 203  call ioExit (routine, 'Error in allocation of tetVerts')

      end subroutine IO_R_UNSTR_X_FAST

c+------------------------------------------------------------------------------
c
      subroutine IO_R_WAKE_FLO107 (iw_bl, iw_face, nw, usr_file)
c
c  Description: IO_R_WAKE_FLO107 reads the wake boundary conditions from a
c  FLO107-MB connectivity file.  The boundary condition data preceding the
c  wake information is skipped.  See IO_R_BC_FLO107 to read this information.
c
c  The usr_file variable is an optional argument.  Absence, or the passing of an
c  unassociated variable in the actual argument, will cause a predefined set of
c  user queries to activate.
c
c  The ncall feature is not implemented for this routine.
c
c  Arguments: See declaration for the complete type specification.
c
c     Intent(inout)
c       usr_file       : File description variable
c
c     Intent(out)
c       nw             : Number of wakes specified in file
c       iw_bl(nw)      : List of blocks in which a wake BC exists
c       iw_face(6,nw)  : Wake specification for each face of the nw blocks
c
c  Side Effects:
c     Since the data pointers are disassociated and allocated during execution,
c     actual argument memory references prior to the call will be lost.
c
c  7/1/99   Recoded with auto-memory feature.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      integer, pointer :: iw_bl(:), iw_face(:,:)
      integer          :: nw
      type (ioMeshFileType), pointer, optional :: usr_file

c  Local constants:

      integer,   parameter :: nface = 6, ncall = -1
      character, parameter :: routine * 16 = 'IO_R_WAKE_FLO107'

c  Local variables:

      integer ::
     .  ig, ioStat, ioUnit, l, ng
      logical ::
     .  exist
      type (ioMeshFileType), pointer ::
     .  prg_file

c  Execution:

c     Test if the programmer is using the optional file argument,
c     the unassociated file argument, or preinitialized options.

      call io_Qry_File (usr_file, prg_file, ncall, routine)

c     Check the file attributes for errors, then open the file:

      call io_Prep_Open_X_File (prg_file, ioUnit, ncall, routine)

      inquire (file=prg_file%name, exist=exist)
      if (.not. exist) goto 103
      open (unit=ioUnit, file=prg_file%name,
     .      form=prg_file%form, status=prg_file%status)

c     Skip boundary condition data to reach the wake section.

      if (prg_file%form == FORMATTED) then           !Formatted File
        read (ioUnit, *, err=101, iostat=ioStat) ng    !BCs
        if (ioStat /= 0) goto 101
        do ig = 1, ng
          read (ioUnit, *, err=102, iostat=ioStat)
          if (ioStat /= 0) goto 102
        end do
        read (ioUnit, *, err=103, iostat=ioStat)       !Wake title
        read (ioUnit, *, err=104, iostat=ioStat) nw    !Number of wake BCs
        if (ioStat /= 0) goto 104
      else                                           !Unformatted file
        read (ioUnit, err=101, iostat=ioStat) ng       !BCs
        if (ioStat /= 0) goto 101
        do ig = 1, ng
          read (ioUnit, err=102, iostat=ioStat)
          if (ioStat /= 0) goto 102
        end do
        read (ioUnit, err=103, iostat=ioStat)          !Wake title
        read (ioUnit, err=104, iostat=ioStat) nw       !Number of wake BCs
        if (ioStat /= 0) goto 104
      end if

c     Allocate wake block and face arrays.

      if (associated (iw_bl)) deallocate (iw_bl)
      allocate (iw_bl(nw), stat=ioStat)
      if (ioStat /= 0) goto 200

      if (associated (iw_face)) deallocate (iw_face)
      allocate (iw_face(nface,nw), stat=ioStat)
      if (ioStat /= 0) goto 201

c     Read block and face data based on format.

      if (prg_file%form == FORMATTED) then
        do ig = 1,nw
          read (ioUnit, *, err=105, iostat=ioStat)
     .      iw_bl(ig), (iw_face(l,ig), l = 1, nface)
          if (ioStat /= 0) goto 105
        end do
      else
        do ig = 1,nw
          read (ioUnit, err=105, iostat=ioStat)
     .      iw_bl(ig), (iw_face(l,ig), l = 1, nface)
          if (ioStat /= 0) goto 105
        end do
      end if

      close (ioUnit)

      return

c  Error handling:

 101  call ioExit (routine, 'Error reading ng')
 102  call ioExit (routine, 'Error skipping bc data')
 103  call ioExit (routine, 'Error skipping wake title')
 104  call ioExit (routine, 'Error reading number of wakes')
 105  call ioExit (routine, 'Error reading a wake line')

 200  call ioExit (routine, 'Error in allocation of iw_bl')
 201  call ioExit (routine, 'Error in allocation of iw_face')

      end subroutine IO_R_WAKE_FLO107

c+------------------------------------------------------------------------------
c
      subroutine IO_W_FNC_P3D (q, ng, nfnc, ijk, ip, ncall, usr_file)
c
c  Description: IO_W_FNC_P3D writes a function file in PLOT3D format, but is
c  general enough to write files containing an arbitrary number of scalars.
c
c  The usr_file variable is an optional argument.  Absence, or the passing of an
c  unassociated variable in the actual argument, will cause an interactive user
c  query session to activate.
c
c  The ncall variable controls to what extent the file will be written.
c  Ncall = -1 results in the entire file, header and all BCs, being written.
c  Ncall = 0  requests that only the header be written.
c  Ncall > 0  requires that the header has been written in a previous call and
c             requests that only the next block's function data be written.
c
c  Note that the header of a PLOT3D function file consists of ng, nfnc, & ijk.
c
c  Arguments: See declaration for the complete type specification.
c
c     Intent(inout)
c       usr_file       : File description variable, optional except for
c                        ncall > 0 which requires the file to be previously
c                        opened with the ncall = 0 option
c
c     If ncall = -1    : Write entire file
c      Intent(in)
c       ng             : Number of grids
c       nfnc           : Number of scalars contained in file
c       ijk(ndim,ng)   : Computational dimensions of each grid
c       ip(ng+1)       : Index specifying beginning of mesh i in q array;
c                        ip(ng+1) is the total number of points + 1
c       q              : Packed scalar function data, for ng grids
c                      : q(nfnc,npts) => usr_file%triples = .T.
c                      : q(npts,nfnc) => usr_file%triples = .F.
c
c     If ncall = 0     : Write file header
c      Intent(in)
c       ng             : Number of grids
c       nfnc           : Number of scalars contained in file
c       ijk(ndim,ng)   : Computational dimensions of each grid
c      Intent(none)
c       ip             : Variable is not accessed
c       q              : Variable is not accessed
c
c     If ncall > 0     : Write next data block
c      Intent(in)
c       ng             : Number of grids
c       nfnc           : Number of scalars contained in file
c       ijk(ndim,ng)   : Computational dimensions of each grid
c       ip(2)          : Degenerate case in which ip(1) = 1 and
c                        ip(2) = number of points in ncall_th block + 1
c       q              : Packed scalar function data for next grid
c                      : q(nfnc,npts) => usr_file%triples = .T.
c                      : q(npts,nfnc) => usr_file%triples = .F.
c
c  Side Effects:
c     The file will be closed at the time of return when ncall = -1 or ng.
c     However, the file will remain open for ng = [0,ng-1].
c
c  Warnings:
c     Use of ncall >= 0, requires the presence of the usr_file argument.
c     This does not imply that the argument must be associated at the
c     time of the call when ncall = 0.
c
c     The value of ncall > 0 must correspond to the block number at the
c     current file pointer location.  The routine does not rewind/skip
c     ncall-1 times, but assumes that the current location is correct.
c
c  7/1/99   Rewrote after adding ncall options.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      real, pointer    :: q(:,:)
      integer, pointer :: ijk(:,:), ip(:)
      integer          :: ng, nfnc, ncall
      type (ioMeshFileType), pointer, optional :: usr_file

c  Local constants:

      character, parameter :: routine * 12 = 'IO_W_FNC_P3D'

c  Local variables:

      integer ::
     .  i, ig, ig_tot, ioStat, ioUnit, j, ndim, ni, nj, npts
      logical ::
     .  ascii, tuples
      type (ioMeshFileType), pointer ::
     .  prg_file

c  Execution:

c     Test if the programmer is using the optional file argument,
c     the unassociated file argument, or preinitialized options.

      call io_Qry_Str_X_File (usr_file, prg_file, ncall, routine)

c     Check the file attributes for errors, then open the file.

      call io_Prep_Open_X_File (prg_file, ioUnit, ncall, routine)

      ndim  = prg_file%ndim
      ascii = prg_file%form == FORMATTED

c     If requested, write the header information.

      if (ncall <= 0) then

        open (unit=ioUnit, file=prg_file%name,        !Agreement of status and
     .    form=prg_file%form, status=prg_file%status) !existence is checked in
                                                      !io_Prep_Open_X_File

        if (prg_file%mg) then                 !Begin writing header
          if (ascii) then
            write (ioUnit, *, err=100, iostat=ioStat) ng !Number of grids
          else
            write (ioUnit,    err=100, iostat=ioStat) ng
          end if
          if (ioStat /= 0) goto 100
        else
          ng = 1  !Ensure it
        end if

        if (.not. associated (ijk)) goto 200  !Check for proper association of
        if (size (ijk,1) /= ndim)   goto 201  !ijk data
        if (size (ijk,2) /= ng)     goto 202

        if (ascii) then  !Write ijk data
          write (ioUnit, '(3I5)', err=101, iostat=ioStat)
     .      ((ijk(j,ig), j = 1, ndim), nfnc, ig = 1, ng)
        else
          write (ioUnit, err=101, iostat=ioStat)
     .      ((ijk(j,ig), j = 1, ndim), nfnc, ig = 1, ng)
        end if

        if (ioStat /= 0) goto 101
        if (ncall  == 0) goto 999             !Finished writing header
      end if

c     Set the loop counter to the number of blocks being written.

      if (ncall < 0) then
        ig_tot = ng
      else
        ig_tot = 1
      end if

c     Test ip and q arrays for association and proper dimensions.

      if (.not. associated (q))  goto 203
      if (.not. associated (ip)) goto 204

      if (size (ip,1) /= ig_tot+1) goto 205

      npts = ip(ig_tot+1) - 1
      tuples = prg_file%triples

      if (tuples) then
        ni = ndim
        nj = npts
        i = 1
        j = 2
      else
        ni = npts
        nj = ndim
        i = 2
        j = 1
      end if

      if (size (q,i) /= nfnc) goto 206
      if (size (q,j) /= npts) goto 207

c     Write block data in correct format.

      do ig = 1, ig_tot
        call IO_W_BLOCK (ig_tot, ig, ip, ni, nj, 0, q, ip,
     .                   tuples, ascii, ioUnit, ioStat)
        if (ioStat /= 0) goto 102
      end do

      if (ncall < 0 .or. ncall == ng) close (ioUnit)

 999  return

c  Error handling:

 100  call ioExit (routine, 'Write ng returned error')
 101  call ioExit (routine, 'Write ijk returned error')
 102  call ioExit (routine, 'Write q returned error')

 200  call ioExit (routine, 'Working array not associated: ijk')
 201  call ioExit (routine, 'Size of ijk inconsistent with ndim')
 202  call ioExit (routine, 'Size of ijk inconsistent with ng')
 203  call ioExit (routine, 'Working array q not associated')
 204  call ioExit (routine, 'Working array ip not associated')
 205  call ioExit (routine, 'Size of ip inconsistent with # blocks')
 206  call ioExit (routine, 'Size of q inconsistent with nfnc')
 207  call ioExit (routine, 'Size of q inconsistent with ip(last)-1')

      end subroutine IO_W_FNC_P3D

c+------------------------------------------------------------------------------
c
      subroutine IO_W_Q_P3D (q, ng, nq, ijk, ip,
     .                       fsmach, alpha, re, time, ncall, usr_file)
c
c  Description: IO_W_Q_P3D writes a solution file in PLOT3D format.
c  PLOT3D format stipulates that the solution file contain only pgr_file%ndim+2
c  variables, i.e., no species or turbulent quantities are expected.  A routine
c  for writing an OVERFLOW solution file may be implemented as IO_W_Q_OVERFLOW.
c
c  The usr_file variable is an optional argument.  Absence, or the passing of an
c  unassociated variable in the actual argument, will cause a predefined set of
c  user queries to activate.
c
c  The ncall variable controls to what extent the file will be written.
c  Ncall = -1 results in the entire file, header and all BCs, being written.
c  Ncall = 0  requests that only the header be written.
c  Ncall > 0  requires that the header has been written in a previous call and
c             requests that only the next block's q-data be written.
c
c  Note that the header of a PLOT3D solution file consists of ng, nq, & ijk.
c
c  Arguments: See declaration for the complete type specification.
c
c     Intent(inout)
c       usr_file       : File description variable, optional except for
c                        ncall > 0 which requires the file to be previously
c                        opened with the ncall = 0 option
c
c     If ncall = -1    : Write whole file
c      Intent(in)
c       ng             : Number of grids
c       nq             : Assumed to be ndim+2
c       ijk(ndim,ng)   : Computational dimensions of each grid
c       ip(ng+1)       : Index specifying beginning of mesh i in q array;
c                        ip(ng+1) is the total number of points + 1
c       q              : Packed scalar function data, for ng grids
c                      : q(nq,npts) => usr_file%triples = .T.
c                      : q(npts,nq) => usr_file%triples = .F.
c       fsmach, alpha  : Scalar quantities assumed to be the same
c       re, time       : for each set of q data
c
c     If ncall = 0     : Write header of file
c      Intent(in)
c       ng             : Number of grids
c       nq             : Assumed to be ndim+2
c       ijk(ndim,ng)   : Computational dimensions of each grid
c
c      Intent(none)
c       ip             : Variable is not accessed
c       q              : Variable is not accessed
c       fsmach, alpha  : Variable is not accessed
c       re, time       :
c
c     If ncall > 0     : Write next data block
c      Intent(in)
c       ng             : Number of grids
c       nq             : Assumed to be ndim+2
c       ijk(ndim,ng)   : Computational dimensions of each grid
c       ip(2)          : Degenerate case in which ip(1) = 1 and
c                        ip(2) = (number of points in ncall_th block) + 1
c       q              : Packed scalar q data, for next grid
c                      : q(nq,npts) => usr_file%triples = .T.
c                      : q(npts,nq) => usr_file%triples = .F.
c       fsmach, alpha  : Scalar quantities assumed to be the same
c       re, time       : for each set of q data
c
c  Side Effects:
c     The file will be closed at the time of return when ncall = -1 or ng.
c     However, the file will remain open for ng = [0,ng-1].
c
c  Warnings:
c     Use of ncall >= 0, requires the presence of the usr_file argument.
c     This does not imply that the argument must be associated at the
c     time of the call when ncall = 0.
c
c     The value of ncall > 0 must correspond to the block number at the
c     current file pointer location.  The routine does not rewind/skip
c     ncall-1 times, but assumes that the current location is correct.
c
c  7/8/99   Rewrote after adding ncall options.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      real, pointer    :: q(:,:)
      integer, pointer :: ijk(:,:), ip(:)
      integer          :: ng, nq, ncall
      real             :: fsmach, alpha, re, time
      type (ioMeshFileType), pointer, optional :: usr_file

c  Local constants:

      character, parameter :: routine * 10 = 'IO_W_Q_P3D'

c  Local variables:

      integer ::
     .  i, ig, ig_tot, ioStat, ioUnit, j, ndim, ni, nj, npts
      logical ::
     .  ascii, tuples
      type (ioMeshFileType), pointer ::
     .  prg_file

c  Execution:

c     Test if the programmer is using the optional file argument,
c     the unassociated file argument, or preinitialized options.

      call io_Qry_File (usr_file, prg_file, ncall, routine)

c     Check the file attributes for errors, then open the file:

      call io_Prep_Open_X_File (prg_file, ioUnit, ncall, routine)

      ndim  = prg_file%ndim
      ascii = prg_file%form == FORMATTED

c     If requested, open file and write header.

      if (ncall <= 0) then

        open (unit=ioUnit, file=prg_file%name,        !Agreement of status and
     .    form=prg_file%form, status=prg_file%status) !existence is checked in
                                                      !io_Prep_Open_X

        if (prg_file%mg) then                       !Begin writing header
          if (ascii) then
            write (ioUnit, *, err=100, iostat=ioStat) ng !Number of grids
          else
            write (ioUnit,    err=100, iostat=ioStat) ng
          end if
          if (ioStat /= 0) goto 100
        else
          ng = 1  !Ensure it
        end if

        if (.not. associated (ijk)) goto 200      !Check association and
        if (size (ijk,1) /= ndim)   goto 201      !dimension of ijk before
        if (size (ijk,2) /= ng)     goto 202      !writing to file

        if (ascii) then
          write (ioUnit, '(3I5)', err=101, iostat=ioStat)
     .      ((ijk(j,ig), j = 1, ndim), ig = 1, ng)
        else
          write (ioUnit,          err=101, iostat=ioStat)
     .      ((ijk(j,ig), j = 1, ndim), ig = 1, ng)
        end if
        if (ioStat /= 0) goto 101

        if (ncall == 0) goto 999                  !Finished writing header
      end if

c     Set the loop counter to the number of blocks to be written.

      if (ncall < 0) then  !Set the grid block loop counter
        ig_tot = ng
      else
        ig_tot = 1
      end if

c     Check for proper association and dimension of q and ip arrays.

      if (.not. associated (q))  goto 203
      if (.not. associated (ip)) goto 204

      if (size (ip,1) /= ig_tot+1) goto 205

      npts = ip(ig_tot+1) - 1
      tuples = prg_file%triples

      if (tuples) then
        ni = ndim
        nj = npts
        i = 1
        j = 2
      else
        ni = npts
        nj = ndim
        i = 2
        j = 1
      end if

      if (size (q,i) /= nq)   goto 206
      if (size (q,j) /= npts) goto 207

c     Write block data in the correct format.

      do ig = 1, ig_tot

        if (ascii) then
          write (ioUnit, '(4F16.8)', err=103, iostat=ioStat)
     .      fsmach, alpha, re, time
        else
          write (ioUnit,             err=103, iostat=ioStat)
     .      fsmach, alpha, re, time
        end if
        if (ioStat /= 0) goto 103

        call IO_W_BLOCK (ig_tot, ig, ip, ni, nj, 0, q, ip,
     .                   tuples, ascii, ioUnit, ioStat)
      end do

      if (ncall < 0 .or. ncall == ng) close (ioUnit)

 999  return

c  Error handling:

 100  call ioExit (routine, 'Write ng returned error')
 101  call ioExit (routine, 'Write ijk returned error')
 103  call ioExit (routine, 'Write fsmach returned error')

 200  call ioExit (routine, 'IJK array not associated')
 201  call ioExit (routine, 'Size of ijk inconsistent with ndim')
 202  call ioExit (routine, 'Size of ijk inconsistent with ng')
 203  call ioExit (routine, 'Q array not associated')
 204  call ioExit (routine, 'IParray ip not associated')
 205  call ioExit (routine, 'Size of ip inconsistent with # blocks')
 206  call ioExit (routine, 'Size of q inconsistent with nq')
 207  call ioExit (routine, 'Size of q inconsistent with ip(last)-1')

      end subroutine IO_W_Q_P3D

c+------------------------------------------------------------------------------
c
      subroutine IO_W_STR_X_P3D (x, ng, ijk, ip, ncall, usr_file)
c
c  Description: IO_W_STR_X_P3D writes an unblanked, structured mesh file in
c  PLOT3D format.  Refer to IO_W_STR_X_P3D_IB for the iblanked mesh format.
c
c  The usr_file variable is an optional argument.  Absence, or the passing of an
c  unassociated variable in the actual argument, will cause a predefined set of
c  user queries to activate.
c
c  The ncall variable controls to what extent the file will be written.
c  Ncall = -1 results in the entire file, header and all BCs, being written.
c  Ncall = 0  requests that only the header be written.
c  Ncall > 0  requires that the header has been written in a previous call and
c             requests that only the next block's x data be written.
c
c  Note that the header of a PLOT3D solution file consists of ng and ijk.
c
c  Arguments: See declaration for the complete type specification.
c
c     Intent(inout)
c       usr_file       : File description variable, optional except for
c                        ncall > 0 which requires the file to be previously
c                        opened with the ncall = 0 option.
c
c     If ncall = -1    : Write whole file
c      Intent(in)
c       ng             : Number of grids.
c       ijk(ndim,ng)   : Computational dimensions of each grid.
c       ip(ng+1)       : Index specifying beginning of mesh i in x array;
c                        ip(ng+1) is the total number of points + 1.
c       x              : Packed mesh data for ng grids.
c                      : x(ndim,npts) => usr_file%triples = .T.
c                      : x(npts,ndim) => usr_file%triples = .F.
c
c     If ncall = 0     : Write header of file.
c      Intent(in)
c       ng             : Number of grids.
c       ijk(ndim,ng)   : Computational dimensions of each grid.
c      Intent(none)
c       ip             : Variable is not accessed.
c       x              : Variable is not accessed.
c
c     If ncall > 0     : Write next data block.
c      Intent(in)
c       ng             : Number of grids.
c       ijk(ndim,ng)   : Computational dimensions of each grid.
c       ip(2)          : Degenerate case in which ip(1) = 1 and
c                        ip(2) = (npts in ncall_th block) + 1.
c       x              : Packed mesh data for next grid.
c                      : x(nq,npts) => usr_file%triples = .T.
c                      : x(npts,nq) => usr_file%triples = .F.
c
c  Side Effects:
c     The file will be closed at the time of return when ncall = -1 or ng.
c     However, the file will remain open for ng = [0,ng-1].
c
c  Warnings:
c     Use of ncall >= 0, requires the presence of the usr_file argument.
c     This does not imply that the argument must be associated at the
c     time of the call when ncall = 0.
c
c     The value of ncall > 0 must correspond to the block number
c     at the current file pointer location.  The routine does not
c     rewind/skip ncall-1 times, but assumes that the current location
c     is correct.
c
c  06/21/99  MJR  Initial coding with auto-memory handling.
c  07/01/99  MJR  Rewrote after adding ncall options.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      real, pointer    :: x(:,:)
      integer, pointer :: ijk(:,:), ip(:)
      integer          :: ng, ncall
      type (ioMeshFileType), pointer, optional :: usr_file

c  Local constants:

      character :: routine * 14 = 'IO_W_STR_X_P3D'

c  Local variables:

      integer ::
     .  i, ig, ig_tot, ioStat, ioUnit, j, ndim, ni, nj, npts
      logical ::
     .  ascii, tuples
      type (ioMeshFileType), pointer ::
     .  prg_file

c  Execution:

c     Test if the programmer is using the optional file argument,
c     the unassociated file argument, or pre-initialized options.

      call io_Qry_Str_X_File (usr_file, prg_file, ncall, routine)

c     Check the file attributes for errors, then open the file:

      call io_Prep_Open_X_File (prg_file, ioUnit, ncall, routine)

      ndim  = prg_file%ndim
      ascii = prg_file%form == FORMATTED

c     If requested, open file and write header.

      if (ncall <= 0) then

        open (unit=ioUnit, file=prg_file%name,        !Agreement of status and
     .    form=prg_file%form, status=prg_file%status) !existence is checked
                                                      !in io_Prep_Open_X

        if (prg_file%mg) then                         !Begin writing header
          if (ascii) then
            write (ioUnit, *, err=100, iostat=ioStat) ng
          else
            write (ioUnit,    err=100, iostat=ioStat) ng
          end if
          if (ioStat /= 0) goto 100
        else
          ng = 1  !Ensure that ng=1
        end if

        if (.not. associated (ijk)) goto 200          !Check that ijk is
        if (size (ijk,1) /= ndim)   goto 201          !associated and properly
        if (size (ijk,2) /= ng)     goto 202          !dimensioned

        if (ascii) then
          write (ioUnit, '(3I5)', err=101, iostat=ioStat) ijk
        else
          write (ioUnit,          err=101, iostat=ioStat) ijk
        end if
        if (ioStat /= 0) goto 101

        if (ncall  == 0) goto 999                     !Finished writing header
      end if

c     Set the loop counter to the number of blocks being written.

      if (ncall < 0) then
        ig_tot = ng
      else
        ig_tot = 1
      end if

c     Check that x and ip arrays are associated and properly dimensioned.
c     Also set triples aliases for implicit loops.

      if (.not. associated (x))  goto 203
      if (.not. associated (ip)) goto 204

      if (size (ip,1) /= ig_tot+1) goto 205

      npts = ip(ig_tot+1) - 1
      tuples = prg_file%triples

      if (tuples) then
        ni = ndim
        nj = npts
        i = 1
        j = 2
      else
        ni = npts
        nj = ndim
        i = 2
        j = 1
      end if

      if (size (x,i) /= ndim) goto 206
      if (size (x,j) /= npts) goto 207

c     Write block data in correct format.

      do ig = 1, ig_tot
        call IO_W_BLOCK (ig_tot, ig, ip, ni, nj, 0, x, ip,
     .                   tuples, ascii, ioUnit, ioStat)
        if (ioStat /= 0) goto 102
      end do

      if (ncall < 0 .or. ncall == ng) close (ioUnit)

 999  return

c  Error handling:

 100  call ioExit  (routine, 'Write ng returned error')
 101  call ioExit  (routine, 'Write ijk returned error')
 102  call ioExit  (routine, 'Write x returned error')

 200  call ioExit  (routine, 'IJK array not associated')
 201  call ioExit  (routine, 'Size of ijk inconsistent w/ ndim')
 202  call ioExit  (routine, 'Size of ijk inconsistent with ng')
 203  call ioExit  (routine, 'X array not associated')
 204  call ioExit  (routine, 'IP array not associated')
 205  call ioExit  (routine, 'Size of ip inconsistent w/ # blks')
 206  call ioExit  (routine, 'Size of x inconsistent with ndim')
 207  call ioExit  (routine, 'Size of x inconsistent with ip(last)-1')

      end subroutine IO_W_STR_X_P3D

c+------------------------------------------------------------------------------
c
      subroutine IO_W_STR_X_P3D_IB (x, ng, ijk, ip, iblank, ncall,
     .                              usr_file)
c
c  Description: IO_W_STR_X_P3D_IB writes a blanked, structured mesh file in
c  PLOT3D format.  Refer to IO_W_STR_X_P3D for the unblanked mesh format.
c
c  The usr_file variable is an optional argument.  Absence, or the passing of an
c  unassociated variable in the actual argument, will cause a predefined set of
c  user queries to activate.
c
c  The ncall variable controls to what extent the file will be written.
c  Ncall = -1 results in the entire file, header and all BCs, being written.
c  Ncall = 0  requests that only the header be written.
c  Ncall > 0  requires that the header has been written in a previous call and
c             requests that only the next block's x data be written.
c
c  Note that the header of a PLOT3D solution file consists of ng and ijk.
c
c  Arguments: See declaration for the complete type specification.
c
c     Intent(inout)
c       usr_file       : File description variable, optional except for
c                        ncall > 0 which requires the file to be previously
c                        opened with the ncall = 0 option
c
c     If ncall = -1    : Write whole file
c      Intent(in)
c       ng             : Number of grids
c       ijk(ndim,ng)   : Computational dimensions of each grid
c       ip(ng+1)       : Index specifying beginning of mesh i in x array;
c                        ip(ng+1) is the total number of points + 1
c       x              : Packed mesh data for ng grids
c                      : x(ndim,npts) => usr_file%triples = .T.
c                      : x(npts,ndim) => usr_file%triples = .F.
c       iblank(npts)   : Packed iblank data for ng grids
c
c     If ncall = 0     : Write header of file
c      Intent(in)
c       ng             : Number of grids
c       ijk(ndim,ng)   : Computational dimensions of each grid
c      Intent(none)
c       ip             : Variable is not accessed
c       x              : Variable is not accessed
c       iblank         : Variable is not accessed
c
c     If ncall > 0     : Write next data block
c      Intent(in)
c       ng             : Number of grids
c       ijk(ndim,ng)   : Computational dimensions of each grid
c       ip(2)          : Degenerate case in which ip(1) = 1 and
c                        ip(2) = (npts in ncall_th block) + 1
c       x              : Packed mesh data for next grid;
c                      : x(nq,npts) => usr_file%triples = .T.
c                      : x(npts,nq) => usr_file%triples = .F.
c       iblank(npts)   : Packed iblank data
c
c  Side Effects:
c     The file will be closed at the time of return when ncall = -1 or ng.
c     However, the file will remain open for ng = [0,ng-1].
c
c  Warnings:
c     Use of ncall >= 0, requires the presence of the usr_file argument.
c     This does not imply that the argument must be associated at the
c     time of the call when ncall = 0.
c
c     The value of ncall > 0 must correspond to the block number at the
c     current file pointer location.  The routine does not rewind/skip
c     ncall-1 times, but assumes that the current location is correct.
c
c  06/21/99  MJR  Initial coding with auto-memory handling.
c  07/01/99  MJR  Rewrote after adding ncall options.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      real, pointer    :: x(:,:)
      integer, pointer :: ijk(:,:), ip(:), iblank(:)
      integer          :: ng, ncall
      type (ioMeshFileType), pointer, optional :: usr_file

c  Local constants:

      character, parameter :: routine * 14 = 'IO_W_STR_X_P3D'

c  Local variables:

      integer ::
     .  i, ig, ig_tot, ioStat, ioUnit, j, ndim, ni, nj, npts
      logical ::
     .  ascii, tuples
      type (ioMeshFileType), pointer ::
     .  prg_file

c  Execution:

c     Test if the programmer is using the optional file argument,
c     the unassociated file argument, or preinitialized options.

      call io_Qry_Str_X_File (usr_file, prg_file, ncall, routine)

c     Check the file attributes for errors, then open the file:

      call io_Prep_Open_X_File (prg_file, ioUnit, ncall, routine)

      ndim  = prg_file%ndim
      ascii = prg_file%form == FORMATTED

c     If requested, open file and write header.

      if (ncall <= 0) then

        open (unit=ioUnit, file=prg_file%name,        !Agreement of status and
     .    form=prg_file%form, status=prg_file%status) !existence is checked in
                                                      !io_Prep_Open_X

        if (prg_file%mg) then                         !Begin writing header
          if (ascii) then
            write (ioUnit, *, err=100, iostat=ioStat) ng
          else
            write (ioUnit,    err=100, iostat=ioStat) ng
          end if
          if (ioStat /= 0) goto 100
        else
            ng = 1  !Ensure it
        end if

        if (.not. associated (ijk)) goto 200          !Check that ijk is
        if (size (ijk,1) /= ndim)   goto 201          !associated and properly
        if (size (ijk,2) /= ng)     goto 202          !dimensioned

        if (ascii) then
          write (ioUnit, '(3I5)', err=101, iostat=ioStat)
     .      ((ijk(j,ig), j = 1, ndim), ig = 1, ng)
        else
          write (ioUnit, err=101, iostat=ioStat)
     .      ((ijk(j,ig), j = 1, ndim), ig = 1, ng)
        end if
        if (ioStat /= 0) goto 101

        if (ncall  == 0) goto 999                     !Finished writing header
      end if

c     Set the loop counter to the number of blocks being written.

      if (ncall < 0) then
        ig_tot = ng
      else
        ig_tot = 1
      end if

c     Check that x and ip arrays are associated and properly dimensioned.

      if (.not. associated (x))      goto 203
      if (.not. associated (ip))     goto 204
      if (.not. associated (iblank)) goto 208

      if (size (ip,1) /= ig_tot+1) goto 205

      npts = ip(ig_tot+1) - 1
      tuples = prg_file%triples

      if (tuples) then
        ni = ndim
        nj = npts
        i = 1
        j = 2
      else
        ni = npts
        nj = ndim
        i = 2
        j = 1
      end if

      if (size (x,i)      /= ndim) goto 206
      if (size (x,j)      /= npts) goto 207
      if (size (iblank,1) /= npts) goto 209

c     Write block data in correct format.

      do ig = 1, ig_tot
        call IO_W_BLOCK (ig_tot, ig, ip, ni, nj, npts, x, iblank,
     .                   tuples, ascii, ioUnit, ioStat)
        if (ioStat /= 0) goto 102
      end do

      if (ncall < 0 .or. ncall == ng) close (ioUnit)

 999  return

c  Error handling:

 100  call ioExit (routine, 'Write ng returned error')
 101  call ioExit (routine, 'Write ijk returned error')
 102  call ioExit (routine, 'Write x returned error')

 200  call ioExit (routine, 'IJK array not associated')
 201  call ioExit (routine, 'Size of ijk inconsistent with ndim')
 202  call ioExit (routine, 'Size of ijk inconsistent with ng')
 203  call ioExit (routine, 'X array not associated')
 204  call ioExit (routine, 'IP array not associated')
 205  call ioExit (routine, 'Size of ip inconsistent with # blocks')
 206  call ioExit (routine, 'Size of x inconsistent with ndim')
 207  call ioExit (routine, 'Size of x inconsistent with ip(last)-1')
 208  call ioExit (routine, 'IBLANK array not associated')
 209  call ioExit (routine, 'iblank size inconsistent with ip(last)-1')

      end subroutine IO_W_STR_X_P3D_IB

c+------------------------------------------------------------------------------
c
      subroutine IO_W_UNSTR_X_FAST (x, nPoints, triVerts, iTriVerts,
     .                              nTris, tetVerts, nTets, ncall,
     .                              usr_file)
c
c  Description: IO_W_UNSTR_X_FAST writes an unstructured mesh file in FAST
c  format.  This routine is suitable only for SINGLE-GRID cases.
c
c  The usr_file variable is an optional argument.  Absence, or the passing of an
c  unassociated variable in the actual argument, will cause a predefined set of
c  user queries to activate.
c
c  The ncall variable controls to what extent the file will be written.
c  Ncall = -1 results in the entire file, header and data, being written.
c  Ncall = 0  requests that only the header be written.
c  Ncall > 0  requires that the header has been written in a previous call and
c             requests that only the x, triVerts, and tetVerts data be written.
c
c  Note that the header of a FAST unformatted file consists of nPoints,
c  nTris, nTets.
c
c  Arguments: See declaration for the complete type specification.
c
c     Intent(inout)
c       usr_file       : File description variable, optional except for
c                        ncall > 0 which requires the file to be previously
c                        opened with the ncall = 0 option
c
c     If ncall = -1    : Write whole file
c      Intent(in)
c       nPoints        : Number of (x,y,z) points in mesh
c       nTri           : Number of triangles in mesh
c       nTet           : Number of tetrahedra in mesh
c
c       x              : Packed (x,y,z) data for ng grids;
c                      : x(ndim,npts) => usr_file%triples = .T.
c                      : x(npts,ndim) => usr_file%triples = .F.
c       triVerts       : List of triangles in mesh
c                      : triVerts(3,nTri)
c       iTriVerts      : Integer flag associated with each triangle
c                      : iTriVerts(nTri)
c       tetVerts       : List of tetrahedra in mesh
c                      : tetVerts(4,nTet)
c
c     If ncall = 0     : Write header of file
c      Intent(in)
c       nPoints        : Number of (x,y,z) points in mesh
c       nTri           : Number of triangles in mesh
c       nTet           : Number of tetrahedra in mesh
c      Intent(none)
c       x              : Returned as an unassociated variable
c       triVerts       : Returned as an unassociated variable
c       iTriVerts      : Returned as an unassociated variable
c       tetVerts       : Returned as an unassociated variable
c
c     If ncall > 0     : Write next data block
c      Intent(in)
c       nPoints        : Number of (x,y,z) points in mesh
c       nTri           : Number of triangles in mesh
c       nTet           : Number of tetrahedra in mesh
c
c       x              : Packed (x,y,z) data for ng grids;
c                      : x(ndim,npts) => usr_file%triples = .T.
c                      : x(npts,ndim) => usr_file%triples = .F.
c       triVerts       : List of triangles in mesh
c                      : triVerts(3,nTri)
c       iTriVerts      : Integer flag associated with each triangle
c                      : iTriVerts(nTri)
c       tetVerts       : List of tetrahedra in mesh
c                      : tetVerts(4,nTet)
c
c  Side Effects:
c     The file will be closed at the time of return when ncall = -1 or
c     ncall > 0.  However, the file will remain open for ncall = 0.
c
c  Warnings:
c     Use of ncall >= 0, requires the presence of the usr_file argument.
c     This does not imply that the argument must be associated at the
c     time of the call when ncall = 0.
c
c     This subroutine is designed for single-grid data only.
c
c  8/8/99   Initial coding with auto-memory handling.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      real, pointer    :: x(:,:)
      integer, pointer :: triVerts(:,:), tetVerts(:,:), iTriVerts(:)
      integer          :: nPoints, nTris, nTets
      integer          :: ncall
      type (ioMeshFileType), pointer, optional :: usr_file

c  Local constants:

      integer,   parameter :: THREE = 3, FOUR = 4
      character, parameter :: routine * 17 = 'IO_W_UNSTR_X_FAST'

c  Local variables:

      integer ::
     .  i, ioStat, ioUnit, ip(2), j, nblank, ndim, ni, nj
      integer, allocatable ::
     .  iblank(:)
      logical ::
     .  ascii, tuples
      type (ioMeshFileType), pointer ::
     .  prg_file

c  Execution:

c     Test if the programmer is using the optional file argument,
c     the unassociated file argument, or preinitialized options.

      call io_Qry_Unstr_X_File (usr_file, prg_file, ncall, routine)

c     Check the file attributes for errors, then open the file:

      call io_Prep_Open_X_File (prg_file, ioUnit, ncall, routine)

      ndim  = prg_file%ndim
      ascii = prg_file%form == FORMATTED

c     If requested, open file and write header.

      if (ncall <= 0) then

        open (unit=ioUnit, file=prg_file%name,        !Agreement of status and
     .    form=prg_file%form, status=prg_file%status) !existence is checked by
                                                      !IO_PREP...
        if (ascii) then                               !Write header line
          write (ioUnit, '(3I10)', err=100, iostat=ioStat)
     .      nPoints, nTris, nTets
        else
          write (ioUnit,           err=100, iostat=ioStat)
     .      nPoints, nTris, nTets
        end if
        if (ioStat /= 0) goto 100

        if (ncall  == 0) goto 999                     !Finished writing header
      end if

c     Check that x, triVerts, iTriVerts, and tetVerts are associated and
c     properly dimensioned.

      if (.not. associated (x))         goto 200
      if (.not. associated (triVerts))  goto 201
      if (.not. associated (iTriVerts)) goto 202
      if (.not. associated (tetVerts))  goto 203

      tuples = prg_file%triples

      if (tuples) then
        ni = ndim
        nj = nPoints
        i = 1
        j = 2
      else
        ni = nPoints
        nj = ndim
        i = 2
        j = 1
      end if

      if (size (x,i) /= ndim)    goto 204
      if (size (x,j) /= nPoints) goto 205

      if (size (triVerts,1) /= THREE) goto 206
      if (size (triVerts,2) /= nTris) goto 207

      if (size (iTriVerts)  /= nTris) goto 208

      if (size (tetVerts,1) /= FOUR)  goto 209
      if (size (tetVerts,2) /= nTets) goto 210

c     Kludge to use the block writing utility:

      nblank = FOUR * (nTris + nTets)
      allocate (iblank(nblank))

      i = 1
      do j = 1, nTris
         iblank(i:i+2) = triVerts(1:THREE,j)
         i = i + THREE
      end do

      j = i + nTris
      iblank(i:j-1) = itriVerts

      i = j
      do j = 1, nTets
         iblank(i:i+THREE) = tetVerts(1:FOUR,j)
         i = i + FOUR
      end do

c     Write block data in correct format.

      ip(1) = 1
      ip(2) = nPoints + 1

      call IO_W_BLOCK (1, 1, ip, ni, nj, nblank, x, iblank,
     .                 tuples, ascii, ioUnit, ioStat)

      deallocate (iblank)

      if (ioStat /= 0) goto 101

      close (ioUnit)

 999  return

c  Error handling:

 100  call ioExit (routine, 'Write nPoints returned error')
 101  call ioExit (routine, 'Write x returned error')

 200  call ioExit (routine, 'x array not associated')
 201  call ioExit (routine, 'triVerts array not associated')
 202  call ioExit (routine, 'iTriVerts array not associated')
 203  call ioExit (routine, 'tetVerts array not associated')
 204  call ioExit (routine, 'Size of x disagrees with ndim')
 205  call ioExit (routine, 'Size of x disagrees with nPoints')
 206  call ioExit (routine, 'Size of triVerts disagrees with THREE')
 207  call ioExit (routine, 'Size of triVerts disagrees with nTris')
 208  call ioExit (routine, 'Size of iTriVerts disagrees with nTris')
 209  call ioExit (routine, 'Size of tetVerts disagrees with FOUR')
 210  call ioExit (routine, 'Size of tetVerts disagrees with nTets')

      end subroutine IO_W_UNSTR_X_FAST

c*******************************************************************************
c                             Private subroutines
c*******************************************************************************

c+------------------------------------------------------------------------------
c
      subroutine IO_R_BLOCK (nblocks, iblock, ip, ni, nj, nblank, x,
     .                       iblank, tuples, ascii, ioUnit, ioStat)
c
c  Description: IO_R_BLOCK reads the indicated block of [un]blanked data in
c  PLOT3D-type order on disk (NOT tuples), formatted or not.  Results are
c  returned in either PLOT3D or tuples order.
c
c  11/06/99  DAS  Pushing the I/O down a level eliminates implied loops.
c  11/22/99   "   Had to push it down ANOTHER level for multiblock x(*,*) cases.
c  11/24/99   "   Had to work harder to include the unstructured case (1 block).
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      integer, intent (in)  :: nblocks,       ! Total # blocks
     .                         iblock,        ! Current block #
     .                         ip(nblocks+1), ! Pointers to start of each block
     .                         ni, nj,        ! Dimensions of x(*,*)
     .                         nblank         ! Length of iblank(*) (may be 0)
      real,    intent (out) :: x(ni,nj)       ! Packed block(s) and the ...
      integer, intent (out) :: iblank(nblank) ! ... blanking (or vertex info)
      logical, intent (in)  :: tuples         ! T = reorder the disk data
      logical, intent (in)  :: ascii          ! Formatted/unformatted flag
      integer, intent (in)  :: ioUnit         ! Unit number to read from
      integer, intent (out) :: ioStat         ! 0 = no read error;
                                              ! 1 = read error
c  Local variables:

      integer               :: i1, ibuf, nintegers, nitems, npoints
      integer, target       :: i, j
      integer, pointer      :: l, m
      real, allocatable     :: xbuf (:,:)

c  Execution:

      i1 = ip(iblock)               ! Start of block in relevant index of x
      npoints = ip(iblock+1) - i1   ! Length of block
      nintegers = nblank            ! If its only one block

      if (nblocks == 1) then        ! Covers the only unstructured case handled

        if (.not. tuples) then      ! Most efficient case

          call IO_R_P3D (npoints, nj, nintegers, x, iblank) ! Internal procedure

          go to 999
        end if

      else ! Multiblock: assume it is structured; no good way for unstructured

        if (nintegers /= 0) nintegers = npoints ! Since nblank is for all blocks

      end if

c     Preserve I/O efficiency by reading into a buffer.

      if (tuples) then
        nitems = ni
        l => j
        m => i
      else
        nitems = nj
        l => i
        m => j
      end if

      allocate (xbuf(npoints,nitems))

      call IO_R_P3D (npoints, nitems, nintegers, xbuf, iblank(i1))

      if (ioStat /= 0) go to 999

      do j = 1, nitems
        ibuf = 0
        do i = ip(iblock), ip(iblock+1) - 1
          ibuf = ibuf + 1
          x(l,m) = xbuf(ibuf,j)
        end do
      end do

      deallocate (xbuf)

 999  return

      contains ! Internal procedure for IO_R_BLOCK

!       -----------------------------------------------------------
        subroutine IO_R_P3D (npoints, nitems, nintegers, x, iblank)
!       -----------------------------------------------------------

!       IO_R_P3D reads one block of [un]blanked data in PLOT3D order on disk
!       (NOT tuples), formatted or not.  Results are returned in PLOT3D order.

!       Arguments:

        integer, intent (in)  :: npoints, nitems,   ! # pts. & # items per pt. &
     .                           nintegers          ! # iblanks for this block
        real,    intent (out) :: x(npoints,nitems)  ! Packed block ...
        integer, intent (out) :: iblank(nintegers)  ! ... and its blanking

!       Execution:

        ioStat = 1

        if (ascii) then
          if (nintegers == 0) then
            read (ioUnit, *, err=999, iostat=ioStat) x
          else
            read (ioUnit, *, err=999, iostat=ioStat) x, iblank
          end if
        else
          if (nintegers == 0) then
            read (ioUnit,    err=999, iostat=ioStat) x
          else
            read (ioUnit,    err=999, iostat=ioStat) x, iblank
          end if
        end if

        ioStat = 0

 999    return

        end subroutine IO_R_P3D

      end subroutine IO_R_BLOCK

c+------------------------------------------------------------------------------
c
      subroutine IO_W_BLOCK (nblocks, iblock, ip, ni, nj, nblank, x,
     .                       iblank, tuples, ascii, ioUnit, ioStat)
c
c  Description: IO_W_BLOCK reads the indicated block of [un]blanked data to
c  disk in PLOT3D-type order (NOT tuples), formatted or not.  Incoming data
c  may be in either PLOT3D or tuples order.
c
c  11/06/99  DAS  Pushing the I/O down a level eliminates implied loops.
c  11/22/99   "   Had to push it down ANOTHER level for multiblock x(*,*) cases.
c  11/24/99   "   Had to work harder to include the unstructured case (1 block).
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      integer, intent (in)  :: nblocks,       ! Total # blocks
     .                         iblock,        ! Current block #
     .                         ip(nblocks+1), ! Pointers to start of each block
     .                         ni, nj,        ! Dimensions of x(*,*)
     .                         nblank         ! Length of iblank(*) (may be 0)
      real,    intent (in)  :: x(ni,nj)       ! Packed block(s) and the ...
      integer, intent (in)  :: iblank(nblank) ! ... blanking (or vertex info)
      logical, intent (in)  :: tuples         ! T = reorder the disk data
      logical, intent (in)  :: ascii          ! Formatted/unformatted flag
      integer, intent (in)  :: ioUnit         ! Unit number to read from
      integer, intent (out) :: ioStat         ! 0 = no read error;
                                              ! 1 = read error
c  Local variables:

      integer               :: i1, ibuf, nintegers, nitems, npoints
      integer, target       :: i, j
      integer, pointer      :: l, m
      real, allocatable     :: xbuf (:,:)

c  Execution:

      i1 = ip(iblock)               ! Start of block in relevant index of x
      npoints = ip(iblock+1) - i1   ! Length of block
      nintegers = nblank            ! If its only one block

      if (nblocks == 1) then        ! Covers the only unstructured case handled

        if (.not. tuples) then      ! Most efficient case

          call IO_W_P3D (npoints, nj, nintegers, x, iblank) ! Internal procedure

          go to 999
        end if

      else ! Multiblock: assume it is structured; no good way for unstructured

        if (nintegers /= 0) nintegers = npoints ! Since nblank is for all blocks

      end if

c     Preserve I/O efficiency by reading into a buffer.

      if (tuples) then
        nitems = ni
        l => j
        m => i
      else
        nitems = nj
        l => i
        m => j
      end if

      allocate (xbuf(npoints,nitems))

      do j = 1, nitems
        ibuf = 0
        do i = ip(iblock), ip(iblock+1) - 1
          ibuf = ibuf + 1
          xbuf(ibuf,j) = x(l,m)
        end do
      end do

      call IO_W_P3D (npoints, nitems, nintegers, xbuf, iblank(i1))

      deallocate (xbuf)

 999  return

      contains ! Internal procedure for IO_R_BLOCK

!       -----------------------------------------------------------
        subroutine IO_W_P3D (npoints, nitems, nintegers, x, iblank)
!       -----------------------------------------------------------

!       IO_W_P3D writes one block of [un]blanked data in PLOT3D order.

!       Arguments:

        integer, intent (in)  :: npoints, nitems,   ! # pts. & # items per pt. &
     .                           nintegers          ! # iblanks for this block
        real,    intent (in)  :: x(npoints,nitems)  ! Packed block ...
        integer, intent (in)  :: iblank(nintegers)  ! ... and its blanking

!       Execution:

        ioStat = 1

        if (ascii) then
          if (nblank == 0) then
            write (ioUnit, '(1p, 6e19.11)', err=999, iostat=ioStat) x
          else
            write (ioUnit, '(1p, 6e19.11)', err=999, iostat=ioStat) x
            write (ioUnit, '(10i8)',   err=999, iostat=ioStat) iblank
          end if
        else
          if (nblank == 0) then
            write (ioUnit, err=999, iostat=ioStat) x
          else
            write (ioUnit, err=999, iostat=ioStat) x, iblank
          end if
        end if

        ioStat = 0

 999    return

        end subroutine IO_W_P3D

      end subroutine IO_W_BLOCK

c+------------------------------------------------------------------------------
c
      subroutine ioFindCleanUnitNumber (ioUnit)
c
c  Description:  In order to prevent the CFD_IO_PACKAGE routines from clashing
c  with the programmer's unit numbers, this routine identifies a unit number
c  that is not currently in use by employing the INQUIRE intrinsic.
c
c  Warning: The routine is hardwired to look for unit numbers in the range
c  [105,299].  Note that valid Cray unit numbers are: [1,4], [7,99], [105,299].
c
c  Fortran 90 allows two successive "open (unit=ioUnit)" statements, where
c  ioUnit is the same in both calls, without causing an error.  F90 will
c  simply close the original file and execute the new open statement.  Thus,
c  if the I/O routines are used in the ncall >= 0 modes, and the programmer
c  unwittingly opens another file with the same unit number between reads/
c  writes to the file, the I/O routine will report that the file has been
c  prematurely closed and exit gracefully.  An expedient but by no means
c  bulletproof work-around is to avoid opening heavily used unit numbers
c  within the package.  Hence only unit numbers within [105,299] are tested.
c
c  Arguments:
c    Intent(out)
c     ioUnit : a unit number not currently in use
c
c  06/21/99  General cleaning of ioPackage
c  07/09/99  Restrict to [105,299]
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      integer, intent(out) :: ioUnit

c  Local constants:

      integer, parameter :: is = 105, ie = 299

c  Local variables:

      integer :: i
      logical :: opened

c  Execution:

      do i = is, ie
       inquire (unit=i, opened=opened)
       if (.not. opened) then
         ioUnit = i
         exit
       end if
      end do

      end subroutine ioFindCleanUnitNumber

c+------------------------------------------------------------------------------
c
      subroutine ioExit (routine, message)
c
c  Description: Gracefully terminate program execution by writing an error
c  message, calling mpi_finalize, flushing the standard output buffer, and
c  stopping.  If not executing in an MPI environment, the MPI_DUMMY library
c  should be linked.
c
c  Arguments:
c     routine  : Routine requesting the exit.
c     message  : Reason for the termination.
c
c  Constants:
c     ioStdOut : Package constant = 6.
c
c  6/21/99  General ioPackage overhaul.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      character(*), intent(in) :: routine, message

c  Procedures:

      external mpi_finalize

c  Execution:

c     Write message, flush the buffer, finalize MPI, and leave.

      write (ioStdOut, '(/,(1x,a))') routine, message,
     .  'CFD_IO_PACKAGE TERMINATING EXECUTION DUE TO ERROR'
      call flush (ioStdOut)
      call mpi_finalize (ioStdOut)
      stop

      end subroutine ioExit

c+------------------------------------------------------------------------------
c
      subroutine io_Prep_Open_X_File (prg_file, ioUnit, ncall, caller)
c
c  Description: In preparation for opening an ioMeshFileType file, this
c  routine performs error checking on the fields of prg_file and checks
c  for possible conflicts with ncall and opened file states.  The routine
c  also returns either a safe unit number under which to open the file,
c  or the unit number the file is currently open under.
c
c  The error checks include:
c     file format string:  FORMATTED || UNFORMATTED
c     file status string:  OLD || NEW || UNKNOWN
c     agreement between file status and existence
c     agreement between ioUnit and previously opened file
c
c  Possible conflicts between ncall and the file state include:
c     For ncall <= 0:  A previously opened file will be closed on return.
c     For ncall >  0:  Conflicts between the prg_file%form attribute and
c                      the specifier used in a previous open statment
c                      will be corrected by modifying prg_file%form.
c  Side effects:
c     The character string attributes of the prg_file variable will be
c     changed to lower case.
c
c  Arguments:
c     prg_file  : File upon which to perform error checks.
c     caller    : Subroutine making the io_Prep_Open_X_File request.
c     ncall     : -1 Read header and all blocks.
c               :  0 Read header.
c               : >0 Read block at current file pointer location,
c                    which implies that the file should be open.
c     ioUnit    : safe unit to use in open statement
c
c  06/21/99  General ioPackage overhaul.
c  06/28/99  Added ncall feature, cleaned.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      type (ioMeshFileType), intent(inout) :: prg_file
      character, intent(in)                :: caller * (*)
      integer, intent(out)                 :: ioUnit
      integer, intent(in)                  :: ncall

c  Local constants:

      character, parameter :: yes * 3 = 'yes'

c  Local variables:

      logical   :: exist, opened
      character :: inq_unf * 3, inq_form * 3

c  Execution:

c     Convert strings to lower case:

      call locase (prg_file%form)
      call locase (prg_file%status)

      inquire (file=prg_file%name, opened=opened, number=ioUnit,
     .         formatted=inq_form, unformatted=inq_unf, exist=exist)

      call locase (inq_form)
      call locase (inq_unf)

c     For cases in which the file should not be open:

      if (ncall <= 0) then

c       Check that the file format is intelligible.

        select case (prg_file%form)
        case (FORMATTED)
        case (UNFORMATTED)
        case default
          goto 100
        end select

c       Check that file status is intelligible and agrees with file existence.

        select case (prg_file%status)
        case (OLD)
          if (.not. exist) goto 101
        case (NEW)
          if (exist) goto 102
        case (UNKNOWN)
        case default
          goto 103
        end select

        if (opened) then      !File is currently open, but agrees
          close (unit=ioUnit) !with expectation on status, so close
        else
          call ioFindCleanUnitNumber (ioUnit)  !Get a new unit number
        end if

      else  !For cases in which the file should be open (ncall > 0):

        if (opened) then             !Force file format to agree with
          if (inq_form == yes) then  !the results of the inquire statement
            prg_file%form = FORMATTED
          else if (inq_unf == yes) then
            prg_file%form = UNFORMATTED
          end if
        else
          goto 104  !Error, the file should be open already.
        end if

      end if

      return

c  Error handling:

 100  call ioExit (caller, 'File format is unintelligible')
 101  call ioExit (caller, 'File status is old but file does not exist')
 102  call ioExit (caller, 'File status is new but file exists')
 103  call ioExit (caller, 'File status is unintelligible')
 104  call ioExit (caller, 'File is closed but ncall > 0')

      end subroutine io_Prep_Open_X_File

c+------------------------------------------------------------------------------
c
      subroutine io_Qry_File (usr_file, prg_file, ncall, caller)
c
c  Description: io_Qry_File provides an optional file attribute querying
c  mechanism for the name, form, & status fields of an ioMeshFileType variable.
c
c  The interactive session is activated if usr_file is either not present or
c  not associated.  For cases in which usr_file is present and unassociated,
c  usr_file memory will be allocated, the user queried, and two distinct copies
c  of the results will be returned, usr_file and prg_file.  If usr_file is not
c  present, no usr_file memory will be allocated and a single copy of the query
c  results will be returned in prg_file.
c
c  If usr_file is present and associated, the interactive session is skipped
c  and a copy of usr_file is returned in prg_file.  Specifying ncall > 0
c  implies that the file has been previously opened, thus usr_file must be
c  present and associated.
c
c  F90 note: A "non present" optional argument may legally be passed to
c  another routine and will retain the "non present" status in the called
c  routine.  E.g., "call io_Qry_File (usr_file, prg_file, ncall, caller)" will
c  still treat usr_file as not present within io_Qry_File if it was not present
c  in the calling routine.
c
c  7/1/99   General cleaning of ioPackage
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      type (ioMeshFileType), pointer, optional ::
     .  usr_file                         !User's file variable
      type (ioMeshFileType), pointer ::
     .  prg_file                         !Program's version of usr_file
      integer ::
     .  ncall                            !Ncall = -1, Read entire file
                                         !Ncall =  0, Read header
                                         !Ncall >  0, Read next block
      character ::
     .  caller * (*)                     !Calling routine; character 4 is
                                         !checked for 'R', 'W', or <other>

c  Local variables:

      integer   :: i_response

c  Execution:

c     Determine whether usr_file name has been pre-set:

      allocate (prg_file)

      if (present (usr_file)) then
        if (associated (usr_file)) then
          prg_file = usr_file     !Present and associated: copy and return
          goto 999
        else
          if (ncall > 0) goto 301 !Present but unassociated: query, copy,
          allocate (usr_file)     !and return
        end if
      else
        if (ncall > 0) goto 300   !Not present: query and return
      end if

c     Since usr_file is either not present or not associated,
c     activate the generic user queries.

      select case (caller(4:4))   !Tell the user what's going on
      case ('R')
        write (ioStdOut, '(/,a)') ' ---File Query For Read---'
      case ('W')
        write (ioStdOut, '(/,a)') ' ---File Query For Write---'
      case default
        write (ioStdOut, '(/,a)') ' ---File Query---'
      end select

      write (ioStdOut, '(/,a)', advance='no')  !Prompt for name
     .  ' Enter mesh filename: '
      read (ioStdIn, '(a)') prg_file%name

 1    write (ioStdOut, '(/, (1x,a))')          !Prompt for format
     . '--------------------------------',
     . ' 1: formatted    2: unformatted',
     . '--------------------------------'
      write (ioStdOut, '(a)', advance='no') ' Enter a file format: '
      read (ioStdIn,*) i_response

      select case (i_response)
      case (1)
        prg_file%form = FORMATTED
      case (2)
        prg_file%form = UNFORMATTED
      case default
        write (ioStdOut, '(a)') ' Error: Select a number (1-2).'
        goto 1
      end select

      prg_file%status = UNKNOWN      !Calling routine will check existence

      if (present (usr_file)) usr_file = prg_file !Return filename if needed

 999  return

c  Error handling:

 300  call ioExit (caller, 'NCALL > 0, but file is optional')
 301  call ioExit (caller, 'NCALL > 0, but file is not associated')

      end subroutine io_Qry_File

c+------------------------------------------------------------------------------
c
      subroutine io_Qry_Str_X_File (usr_file, prg_file, ncall, caller)
c
c  Description: io_Qry_File_X_File provides an optional file attribute querying
c  mechanism for the name, form, status, mg, triples, and ndim fields of a
c  structured ioMeshFileType variable.
c
c  The interactive session is activated if usr_file is either not present or
c  not associated.  For cases in which usr_file is present and unassociated,
c  usr_file memory will be allocated, the user queried, and two distinct copies
c  of the results will be returned, usr_file and prg_file.  If usr_file is not
c  present, no usr_file memory will be allocated and a single copy of the query
c  results will be returned in prg_file.
c
c  If usr_file is present and associated, the interactive session is skipped
c  and a copy of usr_file is returned in prg_file.  Specifying ncall > 0
c  implies that the file has been previously opened, thus usr_file must be
c  present and associated.
c
c  The file%triples attribute refers to the shape of the array in memory,
c  not the format which has been used to store the data in the file.
c
c  F90 note:  A "non present" optional argument may legally be passed to
c  another routine and will retain the "non present" status in the called
c  routine.  E.g., "call io_Qry_File (usr_file, prg_file, ncall, caller)" will
c  still treat usr_file as not present within io_Qry_File if it was not present
c  in the calling routine.
c
c  7/1/99   General Cleaning of ioPackage.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      type (ioMeshFileType), pointer, optional ::
     .  usr_file                         !User's version of variable
      type (ioMeshFileType), pointer ::
     .  prg_file                         !Program's version of usr_file
      integer ::
     .  ncall                            !Ncall = -1, Read entire file
                                         !Ncall =  0, Read header
                                         !Ncall >  0, Read next block
      character ::
     .  caller * (*)                     !Calling routine; character 4 is
                                         !checked for 'R', 'W', or <other>

c  Local variables:

      integer   :: i_response

c  Execution:

c     Determine whether usr_file name is pre-set:

      allocate (prg_file)

      if (present (usr_file)) then
        if (associated (usr_file)) then
          prg_file = usr_file      !Present and associated: copy and return
          goto 999
        else
          if (ncall > 0) goto 301  !Present but unassociated: query, copy,
          allocate (usr_file)      !and return
        end if
      else
        if (ncall > 0) goto 300    !Not present: query and return
      end if

c     Since usr_file is either not present or not associated,
c     activate the generic user queries.

      select case (caller(4:4))    !Tell the user what is going on
      case ('R')
        write (ioStdOut, '(/,a)') ' ---File Query For Read---'
      case ('W')
        write (ioStdOut, '(/,a)') ' ---File Query For Write---'
      case default
        write (ioStdOut, '(/,a)') ' ---File Query---'
      end select

      write (ioStdOut, '(/,a)', advance='no') !Prompt for filename
     .  ' Enter mesh filename: '
      read (ioStdIn, '(a)') prg_file%name

 1    write (ioStdOut, '(/,(1x,a))')          !Prompt for format
     . '--------------------------------',
     . ' 1: mg/form/3d      2: mg/unf/3d',
     . ' 3: sg/form/3d      4: sg/unf/3d',
     . '                                ',
     . ' 5: mg/form/2d      6: mg/unf/2d',
     . ' 7: sg/form/2d      8: sg/unf/2d',
     . '--------------------------------'
      write (ioStdOut, '(a)', advance='no') ' Enter a file format: '
      read (ioStdIn, *) i_response

      select case (i_response)
      case (1)
        prg_file%form = FORMATTED
        prg_file%mg   = .true.
        prg_file%ndim = 3
      case (2)
        prg_file%form = UNFORMATTED
        prg_file%mg   = .true.
        prg_file%ndim = 3
      case (3)
        prg_file%form = FORMATTED
        prg_file%mg   = .false.
        prg_file%ndim = 3
      case (4)
        prg_file%form = UNFORMATTED
        prg_file%mg   = .false.
        prg_file%ndim = 3
      case (5)
        prg_file%form = FORMATTED
        prg_file%mg   = .true.
        prg_file%ndim = 2
      case (6)
        prg_file%form = UNFORMATTED
        prg_file%mg   = .true.
        prg_file%ndim = 2
      case (7)
        prg_file%form = FORMATTED
        prg_file%mg   = .false.
        prg_file%ndim = 2
      case (8)
        prg_file%form = UNFORMATTED
        prg_file%mg   = .false.
        prg_file%ndim = 2
      case default
        write (ioStdOut, '(a)') ' Error: Enter a number (1-8).'
        goto 1
      end select

 2    write (ioStdOut, '(/, (1x,a))')       !Prompt for shape of returned array
     . '--------------------------------',
     . ' 1: x(ndim,npts) 2: x(npts,ndim)',
     . '--------------------------------'
      write (ioStdOut, '(a)', advance='no')
     .  ' Specify the memory storage order: '
      read (ioStdIn, *) i_response

      select case (i_response)
      case (1)
        prg_file%triples = .true.
      case (2)
        prg_file%triples = .false.
      case default
        write (ioStdOut, '(a)') ' Error: Enter a number (1-2).'
        goto 2
      end select

      prg_file%status = UNKNOWN        !Let the calling routine check existence

      if (present (usr_file)) usr_file = prg_file  !Return filename if needed

 999  return

c  Error handling:

 300  call ioExit (caller, 'NCALL > 0, but file is optional.')
 301  call ioExit (caller, 'NCALL > 0, but file is not associated.')

      end subroutine io_Qry_Str_X_File

c+------------------------------------------------------------------------------
c
      subroutine io_Qry_Unstr_X_File (usr_file, prg_file, ncall, caller)
c
c  Description: io_Qry_Unstr_X_File provides an optional file attribute
c  querying mechanism for the name, form, status, mg, triples, and ndim
c  fields of an unstructured ioMeshFileType variable.
c
c  The interactive session is activated if usr_file is either not present
c  or not associated.  For cases in which usr_file is present and unassociated,
c  usr_file memory will be allocated, the user queried, and two distinct copies
c  of the results will be returned, usr_file and prg_file.  If usr_file is not
c  present, no usr_file memory will be allocated and a single copy of the query
c  results will be returned in prg_file.
c
c  If usr_file is present and associated, the interactive session is skipped
c  and a copy of usr_file is returned in prg_file.  Specifying ncall > 0
c  implies that the file has been previously opened, thus usr_file must be
c  present and associated.
c
c  The file%triples attribute refers to the shape of the array in memory,
c  not the format which has been used to store the data in the file.
c
c  F90 note: A "non present" optional argument may legally be passed to
c  another routine and will retain the "non present" status in the called
c  routine.  E.g., "call io_Qry_File (usr_file, prg_file, ncall, caller)" will
c  still treat usr_file as not present within io_Qry_File if it was not present
c  in the calling routine.
c
c  7/1/99   General Cleaning of ioPackage.
c
c  Author:  Mark Rimlinger, Raytheon/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c  Arguments:

      type (ioMeshFileType), pointer, optional ::
     .  usr_file                         !User's version of variable
      type (ioMeshFileType), pointer ::
     .  prg_file                         !Program's version of usr_file
      integer ::
     .  ncall                            !Ncall = -1, Read entire file
                                         !Ncall =  0, Read header
                                         !Ncall >  0, Read next block
      character ::
     .  caller * (*)                     !Calling routine; character 4 is
                                         !checked for 'R', 'W', or <other>

c  Local variables:

      integer   :: i_response

c  Execution:

c     Determine whether usr_file name is pre-set:

      allocate (prg_file)

      if (present (usr_file)) then
        if (associated (usr_file)) then
          prg_file = usr_file      !Present and associated: copy and return
          goto 999
        else
          if (ncall > 0) goto 301  !Present but unassociated: query, copy,
          allocate (usr_file)      !and return
        end if
      else
        if (ncall > 0) goto 300    !Not present: query and return
      end if

c     Since usr_file is either not present or not associated,
c     activate the generic user queries.

      select case (caller(4:4))    !Notify the user what is going on
      case ('R')
        write (ioStdOut, '(/,a)') ' ---File Query For Read---'
      case ('W')
        write (ioStdOut, '(/,a)') ' ---File Query For Write---'
      case default
        write (ioStdOut, '(/,a)') ' ---File Query---'
      end select

      write (ioStdOut, '(/,a)', advance='no')   !Prompt for filename
     .  ' Enter mesh filename: '
      read (ioStdIn, '(a)') prg_file%name

 1    write (ioStdOut, '(/,(1x,a))')            !Prompt for format
     . '--------------------------------',
     . ' 1: formatted     2: unformatted',
     . '--------------------------------'
      write (ioStdOut, '(a)', advance='no') ' Enter a file format: '
      read (ioStdIn, *) i_response

      select case (i_response)
      case (1)
        prg_file%form = FORMATTED
        prg_file%mg   = .false.
        prg_file%ndim = 3
      case (2)
        prg_file%form = UNFORMATTED
        prg_file%mg   = .false.
        prg_file%ndim = 3
      case default
        write (ioStdOut, '(a)') 'Error: Enter a number (1-2).'
        goto 1
      end select

 2    write (ioStdOut, '(/, (1x,a))')       !Prompt for shape of returned array
     . '--------------------------------',
     . ' 1: x(ndim,npts) 2: x(npts,ndim)',
     . '--------------------------------'
      write (ioStdOut, '(a)', advance='no')
     .  ' Specify the memory storage order: '
      read (ioStdIn, *) i_response

      select case (i_response)
      case (1)
        prg_file%triples = .true.
      case (2)
        prg_file%triples = .false.
      case default
        write (ioStdOut, '(a)') ' Error: Enter a number (1-2).'
        goto 2
      end select

      prg_file%status = UNKNOWN        !Let the calling routine check existence

      if (present (usr_file)) usr_file = prg_file  !Return filename if needed

 999  return

c  Error handling:

 300  call ioExit (caller, 'NCALL > 0, but file is optional.')
 301  call ioExit (caller, 'NCALL > 0, but file is not associated.')

      end subroutine io_Qry_Unstr_X_File

C+------------------------------------------------------------------------------
C
      SUBROUTINE LOCASE (STRING)
C
C PURPOSE:  LOCASE changes all upper case letters in the given
C           character string to lower case.
C
C METHOD:   Each character in STRING is treated in turn.  The intrinsic
C           function INDEX effectively allows a table lookup, with
C           the local strings LOW and UPP acting as two tables.
C           This method avoids the use of CHAR and ICHAR, which appear
C           to behave differently on ASCII and EBCDIC machines.
C
C ARGUMENTS:
C    ARG    TYPE I/O/S DESCRIPTION
C  STRING    C   I/O   Character string possibly containing some
C                      uppercase letters on input;
C                      strictly lowercase letters on output with no
C                      change to any non-alphabetic characters.
C
C HISTORY:
C
C  09/10/1985  M. Saunders   Rewrite of version by Hooper/Kennelly, 1983.
C  July 1999   M. Rimlinger  Adaptation of existing UPCASE.
C
C AUTHOR (UPCASE): Michael Saunders, Systems Optimization Lab., Stanford U.
C
C-------------------------------------------------------------------------------

C   Arguments:

      CHARACTER, INTENT (INOUT) :: STRING * (*)

C   Local constants:

      CHARACTER, PARAMETER :: LOW*26 = 'abcdefghijklmnopqrstuvwxyz',
     .                        UPP*26 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
C   Local variables:

      INTEGER    I, J
      CHARACTER  C*1

C   Execution:

      DO J = 1, LEN (STRING)
         C = STRING(J:J)
         IF (C >= 'A' .AND. C <= 'Z') THEN
            I = INDEX (UPP, C)
            IF (I > 0) STRING(J:J) = LOW(I:I)
         END IF
      END DO

      END SUBROUTINE LOCASE

      end module cfd_io_package
