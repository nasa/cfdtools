add_library(searchlib
  asrch3d.f
  bsrchr.f
  hsortcc.f
  hsortii.f
  hsortiind.f
  hsortri.f
  hsortrr.f
  insertxy.f90
  interval.f
  interval2.f
  interval_i.f
  interval_i2.f
  merger.f
  merger2.f
  ripple2d.f
  ripple2p.f
  ripple3d.f
  ripple_surf.f
  ssrchr.f
  unique_ri.f90
  xsort.f
)
target_link_libraries(searchlib PRIVATE global interp2d interp3d)
target_include_directories(searchlib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc RENAME searchlib)
