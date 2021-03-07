cfdtools_add_executable(radial_interp
  radial_interp.f90
)
target_link_libraries(cfdtools_radial_interp
  cfdtools_xyzq_io
  cfdtools_interp3d
  cfdtools_numodules
  cfdtools_adt
  cfdtools_linsys
  cfdtools_optlib
  cfdtools_gridlib
  cfdtools_geomlib
  cfdtools_searchlib
)

if(CFDTOOLS_BUILD_TESTING)

  add_test(
    NAME radial_interp_grid
    COMMAND test.sh ${CMAKE_CURRENT_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/test_grid
  )

  add_test(
    NAME radial_interp_grid_redist
    COMMAND test.sh ${CMAKE_CURRENT_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/test_grid_redist
  )

  add_test(
    NAME radial_interp_soln
    COMMAND test.sh ${CMAKE_CURRENT_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/test_soln
  )

  add_test(
    NAME radial_interp_soln_redist
    COMMAND test.sh ${CMAKE_CURRENT_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/test_soln_redist
  )

endif()
