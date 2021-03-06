add_executable(refine_grid
  densify_grid_block.f90
  refine_grid.f90
)
target_link_libraries(refine_grid cfdtools_global)
target_link_libraries(refine_grid
  interp3d
  interplib
  numodules
  linsys
  optlib
  gridlib
  xyzq_io
  searchlib
  geomlib
)
install(FILES README DESTINATION doc RENAME refine_grid)
install(TARGETS refine_grid DESTINATION ${CMAKE_INSTALL_BINDIR})

if(CFDTOOLS_TEST)

  add_test(
    NAME refine_grid_cube3d
    COMMAND sh -c "refine_grid <cube3d.inp"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test
  )
  set_tests_properties(refine_grid_cube3d PROPERTIES
    ENVIRONMENT "PATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{PATH}"
    PASS_REGULAR_EXPRESSION  "Processing block 1:     9   5  11  -->    9    9   11"
  )

  add_test(
    NAME refine_grid_cube2d
    COMMAND sh -c "refine_grid <cube2d.inp"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test
  )
  set_tests_properties(refine_grid_cube2d PROPERTIES
    ENVIRONMENT "PATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{PATH}"
    PASS_REGULAR_EXPRESSION  "The file appears to be two dimensional"
  )

endif()
