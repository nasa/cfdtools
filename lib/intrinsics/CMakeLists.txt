add_library(intrinsics
  trig_functions_in_degrees.f90
)
target_link_libraries(intrinsics global)
install(FILES README DESTINATION doc)
