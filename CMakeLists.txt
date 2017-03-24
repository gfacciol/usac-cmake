cmake_minimum_required(VERSION 2.8)
project(USAC)

# Are we using gcc?
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	# GCC on MacOs needs this option to use the clang assembler
	if (APPLE)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wa,-q")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-q")
	endif ()
	# Optimize to the current CPU and enable warnings
	#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -std=c99 -O3")
	#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -g")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif ()

# Default build is in Release mode
IF(NOT CMAKE_BUILD_TYPE AND NOT MSVC)
	SET(CMAKE_BUILD_TYPE "Release")
ENDIF(NOT CMAKE_BUILD_TYPE AND NOT MSVC)
SET(CMAKE_BUILD_TYPE "Debug")


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

option(USAC_BUILD_TEST "Build the USAC test program." OFF)

add_definitions(-DLIBCONFIGXX_STATIC)
find_package(Config++)
find_package(Lapack)

add_subdirectory(src)
#if(USAC_BUILD_TEST)
#	add_subdirectory(RunSingleTest)
#endif(USAC_BUILD_TEST)