add_library(special_functions
  error_function.f
  gaussian.f90
  normal_distribution.f90
)
target_link_libraries(special_functions global)
target_include_directories(special_functions PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc)
