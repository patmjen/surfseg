cmake_minimum_required (VERSION 3.5)

project ("surfseg")

if (DEFINED ENV{GEL_LIB_ROOT_DIR})
	set (GEL_LIB_ROOT_DIR "$ENV{GEL_LIB_ROOT_DIR}")
else ()
	set (GEL_LIB_ROOT_DIR ${PROJECT_SOURCE_DIR}/../GEL)
endif ()

set (CMAKE_CXX_STANDARD 14)

if (UNIX AND NOT APPLE)
	set (LINUX true)
endif ()

if (WIN32)
	# For Windows we manually specify the location of the header- and lib directories
	# For other OSs we assume the compiler can find them on its own
	if (CMAKE_SIZEOF_VOID_P EQUAL 8)
		# 64 bits
		link_directories (${GEL_LIB_ROOT_DIR}/GEL_WIN/x64/${CMAKE_BUILD_TYPE})
	elseif (CMAKE_SIZEOF_VOID_P EQUAL 4)
		# 32 bits
		link_directories (${GEL_LIB_ROOT_DIR}/GEL_WIN/${CMAKE_BUILD_TYPE})
	endif ()
endif ()
include_directories (${PROJECT_SOURCE_DIR}/surfseg/include ${PROJECT_SOURCE_DIR}/surfseg/src ${GEL_LIB_ROOT_DIR}/src)

# Include sub-projects.
add_subdirectory ("surfseg")
