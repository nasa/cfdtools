cfdtools_add_executable(table_arithmetic table_arithmetic.f90)
target_link_libraries(cfdtools_table_arithmetic
    cfdtools_numodules
    cfdtools_progtools
    cfdtools_searchlib
    cfdtools_table_io
)
