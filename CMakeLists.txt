cmake_minimum_required(VERSION 3.30 FATAL_ERROR)

project(CHEMKIN LANGUAGES CXX)

add_library(chemkin SHARED)
add_library(lab109::chemkin ALIAS chemkin)

target_include_directories(
    chemkin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                   $<INSTALL_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_compile_features(chemkin PUBLIC cxx_std_23)

add_subdirectory(include)
add_subdirectory(source)
