add_executable(iblank_move iblank_move.f90)
target_link_libraries(iblank_move
  PRIVATE
    global
    gridlib
    xyzq_io
)
install(TARGETS iblank_move DESTINATION ${CMAKE_INSTALL_BINDIR})

# $F90 $FFLAGS -c $ROOT/xyzq_io/grid_block_structure.f90
# $F90 $FFLAGS -c $ROOT/xyzq_io/xyzq_io.f90
# $F90 $FFLAGS -o iblank_move \
#    grid_block_structure.o \
#    xyzq_io.o \
#    iblank_move.f90 \
#    $ROOT/gridlib/f90gridlib.a

