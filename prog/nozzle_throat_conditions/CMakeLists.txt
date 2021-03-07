cfdtools_add_executable(nozzle_throat_conditions
  nozzle_throat_conditions.f90
)
target_link_libraries(cfdtools_nozzle_throat_conditions
  cfdtools_lewis
  cfdtools_progtools
  cfdtools_integratelib
  cfdtools_interplib
  cfdtools_searchlib
  cfdtools_numodules
  cfdtools_optlib
  cfdtools_linsys
)

add_test(
    NAME    nozzle_throat_conditions
    COMMAND sh -c "nozzle_throat_conditions < example.inp"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test
)
set_tests_properties(nozzle_throat_conditions PROPERTIES
    ENVIRONMENT "PATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{PATH}"
    PASS_REGULAR_EXPRESSION "bulk enthalpy:        1.82"
)
