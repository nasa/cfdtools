add_executable(extract_blocks_2d extract_blocks_2d.f90)
target_link_libraries(extract_blocks_2d
  PRIVATE
    global
    progtools
    xyzq_io
)
install(TARGETS extract_blocks_2d DESTINATION ${CMAKE_INSTALL_BINDIR})

#$F90 $FFLAGS -c $ROOT/xyzq_io/grid_block_structure.f90
#$F90 $FFLAGS -c $ROOT/xyzq_io/xyq_io.f90
#$F90 $FFLAGS -o extract_blocks_2d \
#   grid_block_structure.o \
#   xyq_io.o \
#   extract_blocks_2d.f90 \
#   $ROOT/progtools/libprmodules.a

