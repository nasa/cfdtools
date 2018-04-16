add_library(eigenlib
  pythag.f
  rs.f
  tql1.f
  tql2.f
  tred1.f
  tred2.f
)
target_link_libraries(eigenlib global)
install(FILES README DESTINATION doc)
