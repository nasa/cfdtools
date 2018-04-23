add_executable(shadowgraph
  FieldView_Data.f90
  flow_gradients.f90
  grid_block_structure.f90
  process_hex_cell.f90
  process_prism.f90
  process_pyramid.f90
  process_tet_cell.f90
  read_fieldview.f90
  shadowgraph.f90
  test_process_prism.f90
  test_process_pyramid.f90
  test_process_tet_cell.f90
  triangulate_quad_face.f90
  VTK_Data.f90
)
target_link_libraries(shadowgraph global)
target_link_libraries(shadowgraph
  xyzq_io
  interp3d
  numodules
  progtools
  grid_block_utilities
  linsys
  SHADOWGRAPH
  gridlib
  geomlib
)
install(FILES README DESTINATION doc RENAME shadowgraph)
install(TARGETS shadowgraph DESTINATION ${CMAKE_INSTALL_BINDIR})
