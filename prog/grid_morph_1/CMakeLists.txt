add_executable(grid_morph_1
  morph_utilities.f
  grid_morph_1.f90
  read_surface_patch.f90
)
target_link_libraries(grid_morph_1 cfdtools_global)
install(TARGETS grid_morph_1 DESTINATION ${CMAKE_INSTALL_BINDIR})
