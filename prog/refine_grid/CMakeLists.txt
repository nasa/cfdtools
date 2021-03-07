cfdtools_add_executable(refine_grid
  densify_grid_block.f90
  refine_grid.f90
)
target_link_libraries(cfdtools_refine_grid
  cfdtools_interp3d
  cfdtools_interplib
  cfdtools_numodules
  cfdtools_linsys
  cfdtools_optlib
  cfdtools_gridlib
  cfdtools_xyzq_io
  cfdtools_searchlib
  cfdtools_geomlib
)

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
