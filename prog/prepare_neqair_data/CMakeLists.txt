cfdtools_add_executable(prepare_neqair_data
  prepare_neqair_data.f90
)
target_link_libraries(cfdtools_prepare_neqair_data
  cfdtools_gridlib
  cfdtools_progtools
)
