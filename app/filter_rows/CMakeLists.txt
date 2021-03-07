cfdtools_add_executable(filter_rows filter_rows.f90)
target_link_libraries(cfdtools_filter_rows
    cfdtools_progtools
    cfdtools_table_io
)
