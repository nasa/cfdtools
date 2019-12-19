add_library(lewis equilibrium_composition.f90)
add_library(clewis equilibrium_gas.c lewis.c)
target_link_libraries(lewis PRIVATE global clewis)

