add_executable(adjust_flow
  adjust_flow.f90
)
target_link_libraries(adjust_flow global)
target_link_libraries(adjust_flow
  grid_block_structure
  gridlib
  progtools
  geomlib
  xyzq_io
)
install(FILES README DESTINATION doc RENAME adjust_flow)
install(TARGETS adjust_flow DESTINATION ${CMAKE_INSTALL_BINDIR})
