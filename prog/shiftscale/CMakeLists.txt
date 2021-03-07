cfdtools_add_executable(shiftscale
  shiftscale.f
)
target_link_libraries(cfdtools_shiftscale
  cfdtools_numodules
  cfdtools_progtools
)
