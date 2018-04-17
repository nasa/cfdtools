add_library(obsolete
  curv1.f
  htdis.f
  splint.f
  xformx.f
)
target_link_libraries(obsolete global)
target_include_directories(obsolete PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc)
