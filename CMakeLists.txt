find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_library(homography_lib INTERFACE)

target_include_directories(homography_lib INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
