cfdtools_add_library(lewis equilibrium_composition.f90)
cfdtools_add_library(clewis equilibrium_gas.c lewis.c)
target_link_libraries(cfdtools_lewis PRIVATE cfdtools_clewis)

