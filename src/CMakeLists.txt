cmake_minimum_required(VERSION 3.6)
set(MODULE_NAME USAC)
project(${MODULE_NAME})

## Source files -----------------------------------------------------------------------------------

set (SOURCE_FILES
		estimators/AggHomogEstimator.h
		estimators/FundmatrixEstimator.h
		estimators/HomogEstimator.h
		estimators/OriginalAggHomogEstimator.h
		estimators/USAC.h
		config/ConfigFileReader.cpp
		config/ConfigParams.cpp
		config/ConfigParamsHomog.cpp
		utils/HomographyFunctions.cpp
		utils/FundmatrixFunctions.cpp
		utils/MathFunctions.cpp
		utils/Timer.h)
set (SOURCE_IIO
		iio/iio.c
		iio/iio.h)
set (SOURCE_SIFT
		sift/lib_sift.c
		sift/lib_sift.h
		sift/lib_discrete.c
		sift/lib_discrete.h
		sift/lib_matching.c
		sift/lib_matching.h
		sift/lib_keypoint.c
		sift/lib_keypoint.h
		sift/lib_scalespace.c
		sift/lib_scalespace.h
		sift/lib_util.c
		sift/lib_util.h
		sift/lib_sift_anatomy.c
		sift/lib_sift_anatomy.h
		sift/lib_description.c
		sift/lib_description.h)


#file(GLOB_RECURSE sources     "*.cpp")
#file(GLOB_RECURSE headers_hpp "*.hpp")
#file(GLOB_RECURSE headers_h   "*.h")
#list(APPEND headers ${headers_h} ${headers_hpp})


## Project setup ----------------------------------------------------------------------------------

include_directories(
	${CONFIG++_INCLUDE_DIR}
)

link_directories(
	${LAPACK_LIBRARY_DIRS}
)

add_executable(${MODULE_NAME} ${SOURCE_FILES} ${SOURCE_IIO} ${SOURCE_SIFT} RansaacDemo.cpp)

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})

find_package(TIFF REQUIRED)
include_directories(${TIFF_INCLUDE_DIR})

find_package(JPEG REQUIRED)
include_directories(${JPEG_INCLUDE_DIR})

target_link_libraries(${MODULE_NAME}
		${PNG_LIBRARY}
		${TIFF_LIBRARY}
		${JPEG_LIBRARY}
		${CONFIG++_LIBRARY}
		${LAPACK_LIBRARIES}
)

if(MSVC)
	add_definitions(/D_CRT_SECURE_NO_WARNINGS /MP)
endif(MSVC)

