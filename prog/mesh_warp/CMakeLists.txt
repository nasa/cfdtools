add_executable(mesh_warp
  meshwarp.f
  tools_meshwarp.f
)
target_link_libraries(mesh_warp cfdtools_global progtools)
install(TARGETS mesh_warp DESTINATION ${CMAKE_INSTALL_BINDIR})
