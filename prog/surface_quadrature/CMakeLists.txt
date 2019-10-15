add_executable(surface_quadrature surface_quadrature.f90)
target_link_libraries(surface_quadrature
  PRIVATE
    global
    geomlib
    integratelib
    progtools
    tecplot_io
)
install(TARGETS surface_quadrature DESTINATION ${CMAKE_INSTALL_BINDIR})

# F90="ifort"
# #FFLAGS="-r8 -C -traceback -O0 -extend-source 132 -warn general"
# FFLAGS="-r8 -O2 -traceback -extend-source 132 -warn general"
# ROOT="/share/apps/cfdtools"
#
# $F90 $FFLAGS -o surface_quadrature \
#    $ROOT/tecplot_io/tecplot_io.f90 \
#    $ROOT/tecplot_io/tec360_dummy.f90 \
#    surface_quadrature.f90 \
#    $ROOT/integrate/f90integrate.a \
#    $ROOT/geomlib/f90geomlib.a \
#    $ROOT/progtools/libprmodules.a

