add_library(cfdtools_clewis equilibrium_gas.c lewis.c)
install(TARGETS cfdtools_clewis
    EXPORT cfdtools-config
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

cfdtools_add_library(lewis equilibrium_composition.f90)
target_link_libraries(cfdtools_lewis PRIVATE cfdtools_clewis)

