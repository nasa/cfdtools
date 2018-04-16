add_library(blas2lib
  blas2.f
)
target_link_libraries(blas2lib global)
install(FILES README DESTINATION doc)
