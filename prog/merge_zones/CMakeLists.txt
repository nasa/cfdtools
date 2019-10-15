add_executable(merge_zones merge_zones.f90)
target_link_libraries(merge_zones
  PRIVATE
    global
    progtools
    tecplot_io
)
install(TARGETS merge_zones DESTINATION ${CMAKE_INSTALL_BINDIR})

#$F90 $FFLAGS -c $ROOT/tecplot_io/tecplot_io.f90
#$F90 $FFLAGS -c $ROOT/tecplot_io/tec360_dummy.f90
#$F90 $FFLAGS -o merge_zones \
#   merge_zones.f90 \
#   tecplot_io.o \
#   tec360_dummy.o \
#   $ROOT/progtools/libprmodules.a

