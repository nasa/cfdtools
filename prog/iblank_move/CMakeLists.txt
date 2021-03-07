cfdtools_add_executable(iblank_move iblank_move.f90)
target_link_libraries(cfdtools_iblank_move
    cfdtools_gridlib
    cfdtools_xyzq_io
)
