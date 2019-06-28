find_package(Matlab)

set (MESH_LIB_ROOT_DIR ${PROJECT_SOURCE_DIR}/../MESH)

if (Matlab_FOUND)
	include_directories (${Matlab_INCLUDE_DIRS})
	include_directories (${PROJECT_SOURCE_DIR}/surfseg/include ${PROJECT_SOURCE_DIR}/surfseg/src)
	include_directories (${GEL_LIB_ROOT_DIR}/src ${MESH_LIB_ROOT_DIR} ${MESH_LIB_ROOT_DIR}/lib3d/include)
	if (WIN32)
		link_directories (${MESH_LIB_ROOT_DIR}/build/x64-${CMAKE_BUILD_TYPE})
	else ()
		link_directories (${MESH_LIB_ROOT_DIR}/build)
	endif()

	set (CMAKE_POSITION_INDEPENDENT_CODE ON)

	# 3D Surface cuts
	matlab_add_mex (NAME mex_surfcut SRC mex_surfcut.cpp matlab_util.cpp LINK_TO surfseg)
	matlab_add_mex (NAME mex_ksurfcut SRC mex_ksurfcut.cpp matlab_util.cpp LINK_TO surfseg)
	matlab_add_mex (NAME mex_surfcut_planesep SRC mex_surfcut_planesep.cpp matlab_util.cpp LINK_TO surfseg)
	matlab_add_mex (NAME mex_surfcut_planesep_dual SRC mex_surfcut_planesep_dual.cpp matlab_util.cpp LINK_TO surfseg)
	matlab_add_mex (NAME mex_surfcut_planesep_qpbo SRC mex_surfcut_planesep_qpbo.cpp matlab_util.cpp LINK_TO surfseg)

	# 4D surface cuts
	matlab_add_mex (NAME mex_surfcut_4d SRC mex_surfcut_4d.cpp matlab_util.cpp LINK_TO surfseg)

	# Mesh utils
	matlab_add_mex (NAME mex_gauss_curvature SRC mex_gauss_curvature.cpp matlab_util.cpp LINK_TO surfseg)
	matlab_add_mex (NAME mex_star_intersect SRC mex_star_intersect.cpp matlab_util.cpp LINK_TO surfseg)
	matlab_add_mex (NAME mex_hausdorff SRC mex_hausdorff.cpp matlab_util.cpp LINK_TO surfseg mesh)
	matlab_add_mex (NAME mex_subdiv_icosahedron SRC mex_subdiv_icosahedron.cpp matlab_util.cpp LINK_TO surfseg mesh)

	if (WIN32)
		# For some reason we have to specify this manually on Windows
		foreach (EXE mex_surfcut mex_ksurfcut mex_surfcut_planesep mex_surfcut_planesep_dual mex_surfcut_planesep_qpbo mex_surfcut_4d mex_gauss_curvature mex_star_intersect mex_hausdorff mex_subdiv_icosahedron)
			target_link_libraries (${EXE} ${Matlab_ROOT_DIR}/extern/lib/win64/microsoft/libmat.lib ${Matlab_ROOT_DIR}/extern/lib/win64/microsoft/libmx.lib ${Matlab_ROOT_DIR}/extern/lib/win64/microsoft/libut.lib)
			set_target_properties (${EXE} PROPERTIES LINK_FLAGS "/export:mexFunction /debug")
		endforeach ()
	endif ()
endif ()