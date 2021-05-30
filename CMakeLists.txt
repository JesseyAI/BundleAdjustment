cmake_minimum_required(VERSION 2.8)
project(TestBA)

set(DEFAULT_BUILD_TYPE "Debug")

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall")

# Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

# Sophus
find_package(Sophus REQUIRED)
include_directories( ${Sophus_INCLUDE_DIRS} )

# OpenCV
find_package(OpenCV 3 REQUIRED)
include_directories( ${OpenCV_INCLUDE_DIRS})

# add libarary
add_library(sfm_shared SHARED src/Sfm.cpp)
add_library(camera_shared SHARED src/Camera.cpp)
add_library(mappoint_shared SHARED src/MapPoint.cpp)
add_library(costfunction_shared SHARED src/CostFunction.cpp)
add_library(bundleadjustment_shared SHARED src/BundleAdjustment.cpp)

target_link_libraries(sfm_shared
                      camera_shared
                      mappoint_shared
                      costfunction_shared
                      bundleadjustment_shared
                      ${OpenCV_LIBS})

add_executable(testBA testBA.cpp)

target_link_libraries(testBA sfm_shared)
target_link_libraries(testBA ${OpenCV_LIBS})
target_link_libraries(testBA camera_shared)
target_link_libraries(testBA mappoint_shared)
target_link_libraries(testBA costfunction_shared)
target_link_libraries(testBA bundleadjustment_shared)