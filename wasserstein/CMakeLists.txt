project                     (wasserstein)
cmake_minimum_required      (VERSION 2.8.9)

# Default to Release
if                          (NOT CMAKE_BUILD_TYPE)
    set                     (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif                       (NOT CMAKE_BUILD_TYPE)

# Boost
find_package                (Boost)

include_directories         (${CMAKE_CURRENT_SOURCE_DIR}/src
                             ${CMAKE_CURRENT_SOURCE_DIR}/src/include
                             SYSTEM ${Boost_INCLUDE_DIR})
add_definitions             (-std=c++11)

add_executable              (wasserstein src/wasserstein.cpp src/main.cpp src/basic_defs.cpp)
