cmake_minimum_required(VERSION 3.8)
project(homology_relation)
set(CMAKE_CXX_STANDARD 17)

include_directories("include")
add_subdirectory("src")

#add_executable(homology src/main.cpp)
add_executable(homology src/perf_main.cpp)
target_link_libraries(homology lib_homology_construction)

add_executable(org_homology src/new_main.cpp)
target_link_libraries(org_homology lib_homology_construction)

add_executable(forget_homology src/forget_main.cpp)
target_link_libraries(forget_homology lib_homology_construction)

add_executable(meta_homology src/meta_main.cpp)
target_link_libraries(meta_homology lib_homology_construction)

#####################
# TEST
#####################
enable_testing()
add_executable(unit_tests src/test/main.cpp src/test/test_relation_set.cpp)
target_link_libraries(unit_tests lib_homology_construction)
add_test(unit_tests unit_tests)
