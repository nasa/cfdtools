add_library(blaslib
  blas1a.f
  blas1b.f
)
target_link_libraries(blaslib PRIVATE cfdtools_global)
target_include_directories(blaslib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc RENAME blaslib)
