cmake_minimum_required (VERSION 3.5)

# Specify targets
add_library(surfseg STATIC graph.cpp MaxFlow.cpp surface_segment.cpp manifold_mesh.cpp util.cpp subdivided_icosahedron.cpp tet_mesh_4d.cpp surface_segment_4d.cpp)

target_compile_definitions (surfseg PRIVATE NO_OPENGL)

set_property(TARGET surfseg PROPERTY POSITION_INDEPENDENT_CODE ON)

target_link_libraries(surfseg GEL)