cfdtools_add_executable(free_stream_conditions free_stream_conditions.f90)
target_link_libraries(cfdtools_free_stream_conditions
    cfdtools_progtools
    cfdtools_table_io
)
