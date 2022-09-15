cfdtools_add_executable(usinterp usinterp.f90)
target_link_libraries(cfdtools_optinterp PRIVATE
    cfdtools_linsys
    cfdtools_lapacksubset
    cfdtools_progtools
    cfdtools_numodules
)
