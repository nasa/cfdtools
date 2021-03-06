add_library(tecplot_io
  tec360_dummy.f90
  tecplot_io.f90
  # tecplot_tri_io.f90
)
target_link_libraries(tecplot_io PRIVATE cfdtools_global progtools)
target_include_directories(tecplot_io PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc RENAME tecplot_io)
