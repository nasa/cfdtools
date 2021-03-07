cfdtools_add_executable(prepare_fiat_data
  prepare_fiat_data.f90
)
target_link_libraries(cfdtools_prepare_fiat_data
  cfdtools_table_io
  cfdtools_progtools
)
