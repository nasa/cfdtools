add_executable(insert_z
  insert_z.f90
)
target_link_libraries(insert_z global grid_block_structure)
install(TARGETS insert_z DESTINATION ${CMAKE_INSTALL_BINDIR})
