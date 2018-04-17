add_library(optinterp
  optimal_interpolation.f90
)
target_link_libraries(optinterp PRIVATE global)
target_include_directories(optinterp PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc)
