add_executable(uslos uslos.f90)
target_link_libraries(uslos
  PRIVATE
    cfdtools_global
    adt
    eigenlib
    geomlib
    gridlib
    interp3d
    intrinsics
    linsys
    numodules
    optlib
    progtools
    triangulation_io
    ugridlib
    xyzq_io
)

#$F90 $FFLAGS -o uslos \
#   uslos.f90 \
#   grid_block_structure.o \
#   xyzq_io.o \
#   string_justify.o \
#   triangulation_io.o \
#   adt_utilities.o \
#   $ROOT/interp3d/intsec6.f90 \
#   $ROOT/interp3d/intsec8.f90 \
#   $ROOT/ugridlib/spherical_triangulation.f90 \
#   $ROOT/geomlib/f90geomlib.a \
#   $ROOT/gridlib/f90gridlib.a \
#   $ROOT/eigenlib/f90eigenlib.a \
#   $ROOT/optlib/f90optlib.a \
#   $ROOT/ugridlib/f90ugridlib.a \
#   $ROOT/linsys/f90linsys.a \
#   $ROOT/interp3d/f90interp3d.a \
#   $ROOT/searchlib/f90searchlib.a \
#   $ROOT/numodules/f90numodules.a \
#   $ROOT/progtools/libprmodules.a

