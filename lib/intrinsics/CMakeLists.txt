add_library(intrinsics
  trig_functions_in_degrees.f90
)
target_link_libraries(intrinsics PRIVATE global)
target_include_directories(intrinsics PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION intrinsics)
