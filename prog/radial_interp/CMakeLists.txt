add_executable(radial_interp
  radial_interp.f90
)
target_link_libraries(radial_interp cfdtools_global)
target_link_libraries(radial_interp
  xyzq_io
  interp3d
  numodules
  adt
  linsys
  optlib
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME radial_interp)
install(TARGETS radial_interp DESTINATION ${CMAKE_INSTALL_BINDIR})

if(CFDTOOLS_TEST)

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
