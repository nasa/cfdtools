add_executable(tri_diffs tri_diffs.f90)
target_link_libraries(tri_diffs
  PRIVATE
    cfdtools_global
    adt
    eigenlib
    geomlib
    linsys
    progtools
    triangulation_io
    ugridlib
    xyzq_io
)


# $F90 $FFLAGS -o tri_diffs \
#    $ROOT/progtools/string_justify.f90 \
#    $ROOT/triangulation_io/triangulation_io.f90 \
#    $ROOT/xyzq_io/grid_block_structure.f90 \
#    $ROOT/adt/adt_utilities.f90 \
#    tri_diffs.f90 \
#    $ROOT/geomlib/f90geomlib.a \
#    $ROOT/eigenlib/f90eigenlib.a \
#    $ROOT/ugridlib/f90ugridlib.a \
#    $ROOT/linsys/f90linsys.a \
#    $ROOT/progtools/libprmodules.a

