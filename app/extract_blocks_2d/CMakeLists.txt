cfdtools_add_executable(extract_blocks_2d extract_blocks_2d.f90)
target_link_libraries(cfdtools_extract_blocks_2d
    cfdtools_progtools
    cfdtools_xyzq_io
)

