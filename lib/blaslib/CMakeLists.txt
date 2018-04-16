add_library(blaslib
  blas1a.f
  blas1b.f
)
target_link_libraries(blaslib global)
install(FILES README DESTINATION doc)
