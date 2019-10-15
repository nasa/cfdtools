add_executable(slos slos.f90)
target_link_libraries(slos
  PRIVATE
    global
    adt
    eigenlib
    geomlib
    gridlib
    interp3d
    intrinsics
    linsys
    numodules
    progtools
    optlib
    searchlib
    triangulation_io
    ugridlib
    xyzq_io
)
install(TARGETS slos DESTINATION ${CMAKE_INSTALL_BINDIR})



#$F90 $FFLAGS -c $ROOT/xyzq_io/grid_block_structure.f90
#$F90 $FFLAGS -c $ROOT/xyzq_io/xyzq_io.f90
#$F90 $FFLAGS -c $ROOT/progtools/string_justify.f90
#$F90 $FFLAGS -c $ROOT/triangulation_io/triangulation_io.f90
#$F90 $FFLAGS -c $ROOT/adt/adt_utilities.f90
#
## intsec6.o needs the ADT module, so it's not in f90interp3d.a
#
#$F90 $FFLAGS -c $ROOT/interp3d/intsec6.f90
#$F90 $FFLAGS -o slos \
#   slos.f90 \
#   grid_block_structure.o \
#   xyzq_io.o \
#   string_justify.o \
#   triangulation_io.o \
#   adt_utilities.o \
#   intsec6.o \
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

