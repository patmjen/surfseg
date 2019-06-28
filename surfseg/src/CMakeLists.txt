cmake_minimum_required (VERSION 3.5)

# Specify targets
add_library(surfseg STATIC graph.cpp MaxFlow.cpp surface_segment.cpp manifold_mesh.cpp util.cpp subdivided_icosahedron.cpp tet_mesh_4d.cpp surface_segment_4d.cpp)

target_compile_definitions (SurfSegment PRIVATE NO_OPENGL)
target_compile_definitions (SmoothMesh PRIVATE NO_OPENGL)
target_compile_definitions (surfseg PRIVATE NO_OPENGL)

set_property(TARGET surfseg PROPERTY POSITION_INDEPENDENT_CODE ON)

# Setup linking for all targets
foreach (EXE ObjView SubIcoView DscView DscSmooth SurfSegment SmoothMesh surfseg)
	target_link_libraries(${EXE} GEL)
	if (LINUX)
		target_link_libraries (${EXE} GL GLU glut)
	elseif (WIN32)
		target_link_libraries (${EXE} glew32 glut32)
	else ()
		message (FATAL_ERROR "Apple not supported -- you can maybe use the linux target")
	endif ()
endforeach ()
