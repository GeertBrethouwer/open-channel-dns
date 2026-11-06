cmake_minimum_required(VERSION 3.22.1)
project(ocdns CXX)

# The solver itself lives in native/ocdns.hpp at the repository root (shared with the desktop test).
add_library(ocdns SHARED ocdns_jni.cpp)
target_include_directories(ocdns PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../../../../native)
target_compile_features(ocdns PRIVATE cxx_std_17)
# Fast code in every build type; no -ffast-math, so the blow-up check (isfinite) stays reliable.
target_compile_options(ocdns PRIVATE -O3)
target_link_libraries(ocdns PRIVATE log)
