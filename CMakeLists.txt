cmake_minimum_required(VERSION 3.20)
project(ascii_protein VERSION 0.9.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ASCII_PROTEIN_BUILD_TESTS "Build the test suite" ON)

add_library(ascii_protein_core
    src/structure.cpp
    src/pdb_parser.cpp
    src/mmcif_parser.cpp
    src/structure_loader.cpp
    src/framebuffer.cpp
    src/renderer.cpp
    src/cli.cpp
    src/terminal.cpp
    src/rcsb.cpp
    src/command.cpp
    src/selection.cpp
    src/secondary_assignment.cpp
    src/trajectory.cpp
    src/xyz_trajectory_source.cpp
)
target_include_directories(ascii_protein_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(ascii_protein_core PUBLIC cxx_std_20)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(ascii_protein_core PRIVATE -Wall -Wextra -Wpedantic)
elseif(MSVC)
    target_compile_options(ascii_protein_core PRIVATE /W4 /permissive-)
endif()

add_executable(ascii-protein src/main.cpp)
target_link_libraries(ascii-protein PRIVATE ascii_protein_core)
target_compile_definitions(ascii-protein PRIVATE ASCII_PROTEIN_VERSION="${PROJECT_VERSION}")

include(GNUInstallDirs)
install(TARGETS ascii-protein RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(ASCII_PROTEIN_BUILD_TESTS)
    enable_testing()
    add_executable(ascii_protein_tests tests/tests.cpp)
    target_link_libraries(ascii_protein_tests PRIVATE ascii_protein_core)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        target_compile_options(ascii_protein_tests PRIVATE -Wall -Wextra -Wpedantic)
    elseif(MSVC)
        target_compile_options(ascii_protein_tests PRIVATE /W4 /permissive-)
    endif()
    add_test(NAME ascii_protein_tests COMMAND ascii_protein_tests)
endif()