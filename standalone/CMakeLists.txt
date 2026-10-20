cmake_minimum_required(VERSION 3.20)

project(without_sgems_hosc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(COMMON_SOURCES
    src/common/config.cpp
    src/common/dense_grid.cpp
    src/common/grid_io.cpp
    src/common/text.cpp)

set(HOSC2008_SOURCES
    vendor/hosc2008_core/AVL.cpp
    vendor/hosc2008_core/CumulantCal.cpp
    vendor/hosc2008_core/FunctionCall.cpp
    vendor/hosc2008_core/Initialization.cpp
    vendor/hosc2008_core/SearchNode.cpp
    vendor/hosc2008_core/Tools.cpp)

set(HOSCWINDOW_SOURCES
    vendor/hoscwindow_core/AVL.cpp
    vendor/hoscwindow_core/CumulantCal.cpp
    vendor/hoscwindow_core/FunctionCall.cpp
    vendor/hoscwindow_core/Initialization.cpp
    vendor/hoscwindow_core/SearchNode.cpp
    vendor/hoscwindow_core/Tools.cpp)

add_library(weighted_moment_core STATIC vendor/weighted_moment_core/WeightedMomentCore.cpp)
target_include_directories(weighted_moment_core PUBLIC vendor/weighted_moment_core)

add_library(hosc2008_core STATIC ${HOSC2008_SOURCES})
target_include_directories(hosc2008_core PUBLIC vendor/compat vendor/hosc2008_core)

add_library(hoscwindow_core STATIC ${HOSCWINDOW_SOURCES})
target_include_directories(hoscwindow_core PUBLIC vendor/compat vendor/hoscwindow_core)

add_executable(template_cli src/template/main.cpp ${COMMON_SOURCES})
target_include_directories(template_cli PRIVATE src/common vendor/compat vendor/hosc2008_core)
target_link_libraries(template_cli PRIVATE hosc2008_core weighted_moment_core)

add_executable(movewindow_cli
    src/window/main.cpp
    ${COMMON_SOURCES})
target_include_directories(movewindow_cli PRIVATE src/common vendor/compat vendor/hoscwindow_core)
target_link_libraries(movewindow_cli PRIVATE hoscwindow_core weighted_moment_core)

enable_testing()
add_executable(weighted_core_checks tests/compare_weighted_core.cpp)
target_link_libraries(weighted_core_checks PRIVATE weighted_moment_core)
add_test(NAME weighted_core_checks COMMAND weighted_core_checks)

add_executable(window_avl_core_checks tests/compare_window_avl_core.cpp)
target_link_libraries(window_avl_core_checks PRIVATE hoscwindow_core)
target_include_directories(window_avl_core_checks PRIVATE vendor/compat vendor/hoscwindow_core)
add_test(NAME window_avl_core_checks COMMAND window_avl_core_checks)

set(MOVING_WINDOW_FIXTURE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures")
set(MOVING_WINDOW_SMOKE_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/test_output/movewindow")
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/movewindow_cli_smoke.ini.in
    ${CMAKE_CURRENT_BINARY_DIR}/tests/movewindow_cli_smoke.ini
    @ONLY)
add_test(NAME movewindow_cli_smoke
    COMMAND movewindow_cli ${CMAKE_CURRENT_BINARY_DIR}/tests/movewindow_cli_smoke.ini)

set(TEMPLATE_EXHAUSTIVE_SMOKE_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/test_output/template_exhaustive")
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/template_exhaustive_smoke.ini.in
    ${CMAKE_CURRENT_BINARY_DIR}/tests/template_exhaustive_smoke.ini
    @ONLY)
add_test(NAME template_cli_exhaustive_integer_smoke
    COMMAND template_cli ${CMAKE_CURRENT_BINARY_DIR}/tests/template_exhaustive_smoke.ini)

# The formal directional workflow relies on signed integer offsets rather than
# legacy angular search.  Keep this as a process-level test so it exercises the
# production template_cli parser, exhaustive sampler, and AVL evaluator.
add_executable(template_cli_geometry_checks tests/template_cli_geometry_checks.cpp)
add_dependencies(template_cli_geometry_checks template_cli)
add_test(NAME template_cli_exhaustive_integer_geometry
    COMMAND template_cli_geometry_checks $<TARGET_FILE:template_cli>)
set_tests_properties(template_cli_exhaustive_integer_geometry PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
