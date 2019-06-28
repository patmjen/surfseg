cmake_minimum_required (VERSION 3.5)

project ("surfseg")

set (GEL_LIB_ROOT_DIR ${PROJECT_SOURCE_DIR}/../GEL)

set (CMAKE_CXX_STANDARD 14)

if (UNIX AND NOT APPLE)
	set (LINUX true)
endif ()

if (WIN32)
	# For Windows we manually specify the location of the header- and lib directories
    set (EXT_LIB_ROOT_DIR ${PROJECT_SOURCE_DIR}/../../../lib)
	if (CMAKE_SIZEOF_VOID_P EQUAL 8)
		# 64 bits
		link_directories (${EXT_LIB_ROOT_DIR}/freeglut/lib/x64 ${EXT_LIB_ROOT_DIR}/glew-2.1.0/lib/Release/x64 ${GEL_LIB_ROOT_DIR}/GEL_WIN/x64/${CMAKE_BUILD_TYPE})
	elseif (CMAKE_SIZEOF_VOID_P EQUAL 4)
		# 32 bits
		link_directories (${EXT_LIB_ROOT_DIR}/freeglut/lib ${EXT_LIB_ROOT_DIR}/glew-2.1.0/lib/Release/Win32 ${GEL_LIB_ROOT_DIR}/GEL_WIN/${CMAKE_BUILD_TYPE})
	endif ()
    include_directories (${EXT_LIB_ROOT_DIR}/freeglut/include ${EXT_LIB_ROOT_DIR}/glew-2.1.0/include)
endif ()
include_directories (${PROJECT_SOURCE_DIR}/surfseg/include ${PROJECT_SOURCE_DIR}/surfseg/src ${GEL_LIB_ROOT_DIR}/src)

# Include sub-projects.
add_subdirectory ("surfseg")
