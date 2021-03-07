cfdtools_add_library(lewis equilibrium_composition.f90)
add_library(cfdtools_clewis equilibrium_gas.c lewis.c)
install(TARGETS cfdtools_clewis EXPORT cfdtools-config)
target_link_libraries(cfdtools_lewis PRIVATE cfdtools_clewis)

