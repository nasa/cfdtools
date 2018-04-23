add_library(lewis
    equilibrium_composition.f90
    equilibrium_gas.c
    lewis.c
)
# Can't link to global; adds bad options to C files. Manually propagate
# global compile flags to f90 files.
#target_link_libraries(lewis PRIVATE global)
# TODO: This doesn't work. Need to figur out why
#get_target_property(GLOBAL_OPTS global INTERFACE_COMPILE_OPTIONS)
#set_source_files_properties(equilibrium_composition.f90
#  PROPERTIES
#    COMPILE_FLAGS ${GLOBAL_OPTS}
#)
target_include_directories(lewis PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
