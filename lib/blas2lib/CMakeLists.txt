add_library(blas2lib
  blas2.f
)
target_link_libraries(blas2lib PRIVATE global)
target_include_directories(blas2lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc)
