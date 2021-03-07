add_executable(nozzle_throat_conditions
  nozzle_throat_conditions.f90
)
target_link_libraries(nozzle_throat_conditions
  cfdtools_global
  cfdtools_lewis
  cfdtools_progtools
  cfdtools_integratelib
  cfdtools_interplib
  cfdtools_searchlib
  cfdtools_numodules
  cfdtools_optlib
  cfdtools_linsys
)
install(FILES README DESTINATION doc RENAME nozzle_throat_conditions)
install(TARGETS nozzle_throat_conditions DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(
    NAME    nozzle_throat_conditions
    COMMAND sh -c "nozzle_throat_conditions < example.inp"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test
)
set_tests_properties(nozzle_throat_conditions PROPERTIES
    ENVIRONMENT "PATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{PATH}"
    PASS_REGULAR_EXPRESSION "bulk enthalpy:        1.82"
)
