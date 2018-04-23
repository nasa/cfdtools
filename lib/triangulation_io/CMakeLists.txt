add_library(triangulation_io
  triangulation_io.f90
)
target_link_libraries(triangulation_io PRIVATE global eigenlib progtools)
target_include_directories(triangulation_io PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc)
