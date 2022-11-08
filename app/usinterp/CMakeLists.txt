cfdtools_add_executable(usinterp usinterp.f90)
target_link_libraries(cfdtools_usinterp
    cfdtools_linsys
    cfdtools_lapacksubset
    cfdtools_optinterp
    cfdtools_progtools
    cfdtools_numodules
    cfdtools_table_io
)
