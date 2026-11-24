cmake_minimum_required(VERSION 3.20)

project(without_sgems_hosc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(COMMON_SOURCES
    src/common/config.cpp
    src/common/dense_grid.cpp
    src/common/grid_io.cpp
    src/common/point_set.cpp
    src/common/cell_declustering.cpp
    src/common/spatial_declustering.cpp
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

add_executable(template_cli src/template/main.cpp src/template/point_set_template.cpp ${COMMON_SOURCES})
target_include_directories(template_cli PRIVATE src/common vendor/compat vendor/hosc2008_core)
target_link_libraries(template_cli PRIVATE hosc2008_core weighted_moment_core)

add_executable(movewindow_cli
    src/window/main.cpp
    ${COMMON_SOURCES})
target_include_directories(movewindow_cli PRIVATE src/common vendor/compat vendor/hoscwindow_core)
target_link_libraries(movewindow_cli PRIVATE hoscwindow_core weighted_moment_core)

enable_testing()
add_executable(order_cap_checks tests/order_cap_checks.cpp)
target_link_libraries(order_cap_checks PRIVATE hosc2008_core weighted_moment_core)
add_test(NAME order_cap_checks COMMAND order_cap_checks)
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

add_executable(point_set_io_checks tests/point_set_io_checks.cpp ${COMMON_SOURCES})
target_include_directories(point_set_io_checks PRIVATE src/common)
add_test(NAME point_set_io_checks COMMAND point_set_io_checks)

add_executable(template_pointset_checks tests/template_pointset_checks.cpp)
add_dependencies(template_pointset_checks template_cli)
add_test(NAME template_pointset_checks COMMAND template_pointset_checks $<TARGET_FILE:template_cli>)
set_tests_properties(template_pointset_checks PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(window_pointset_checks tests/window_pointset_checks.cpp)
add_dependencies(window_pointset_checks movewindow_cli)
add_test(NAME window_pointset_checks COMMAND window_pointset_checks $<TARGET_FILE:movewindow_cli>)
set_tests_properties(window_pointset_checks PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(cell_declustering_checks tests/cell_declustering_checks.cpp ${COMMON_SOURCES})
target_include_directories(cell_declustering_checks PRIVATE src/common)
add_test(NAME cell_declustering_checks COMMAND cell_declustering_checks)

add_executable(template_declustering_checks tests/template_declustering_checks.cpp)
add_dependencies(template_declustering_checks template_cli)
add_test(NAME template_declustering_checks COMMAND template_declustering_checks $<TARGET_FILE:template_cli>)
set_tests_properties(template_declustering_checks PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(cartesian_declustering_checks tests/cartesian_declustering_checks.cpp)
add_dependencies(cartesian_declustering_checks template_cli)
add_test(NAME cartesian_declustering_checks COMMAND cartesian_declustering_checks $<TARGET_FILE:template_cli>)
set_tests_properties(cartesian_declustering_checks PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(spatial_declustering_checks tests/spatial_declustering_checks.cpp
    src/common/spatial_declustering.cpp)
target_include_directories(spatial_declustering_checks PRIVATE src/common)
add_test(NAME spatial_declustering_checks COMMAND spatial_declustering_checks)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    add_test(NAME external_weight_checks COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/external_weight_checks.py
        --template $<TARGET_FILE:template_cli> --window $<TARGET_FILE:movewindow_cli>
        --out ${CMAKE_CURRENT_BINARY_DIR}/test_output/external_weights)
    add_test(NAME point_window_grid_checks COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/point_window_grid_checks.py
        --cli $<TARGET_FILE:movewindow_cli>
        --out ${CMAKE_CURRENT_BINARY_DIR}/test_output/point_window_grid)
    add_test(NAME window_declustering_checks COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/window_declustering_checks.py
        --cli $<TARGET_FILE:movewindow_cli>
        --out ${CMAKE_CURRENT_BINARY_DIR}/test_output/window_declustering)
endif()
