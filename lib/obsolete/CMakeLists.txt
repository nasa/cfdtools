add_library(obsolete
  clusgeo.f.no-longer-in-use
  curv1.f
  htdis.f
  splint.f
  xformx.f
)
target_link_libraries(obsolete global)
install(FILES README DESTINATION doc)
