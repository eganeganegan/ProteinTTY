if(NOT DEFINED BUILD_DIR OR NOT DEFINED SOURCE_DIR OR NOT DEFINED OUTPUT_DIR OR
   NOT DEFINED VERSION OR NOT DEFINED PLATFORM)
    message(FATAL_ERROR "BUILD_DIR, SOURCE_DIR, OUTPUT_DIR, VERSION, and PLATFORM are required")
endif()

set(executable "${BUILD_DIR}/ascii-protein")
if(WIN32)
    set(executable "${BUILD_DIR}/Release/ascii-protein.exe")
    if(NOT EXISTS "${executable}")
        set(executable "${BUILD_DIR}/ascii-protein.exe")
    endif()
endif()
if(NOT EXISTS "${executable}")
    message(FATAL_ERROR "release executable not found: ${executable}")
endif()

string(REGEX REPLACE "^v" "" package_version "${VERSION}")
set(package_name "ascii-protein-${package_version}-${PLATFORM}")
set(staging_root "${OUTPUT_DIR}/staging")
set(package_dir "${staging_root}/${package_name}")
file(REMOVE_RECURSE "${staging_root}")
file(MAKE_DIRECTORY "${package_dir}")
file(COPY "${executable}" "${SOURCE_DIR}/README.md" "${SOURCE_DIR}/LICENSE" DESTINATION "${package_dir}")

if(WIN32)
    set(archive "${OUTPUT_DIR}/${package_name}.zip")
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E tar cf "${archive}" --format=zip "${package_name}"
        WORKING_DIRECTORY "${staging_root}"
        RESULT_VARIABLE archive_result
    )
else()
    set(archive "${OUTPUT_DIR}/${package_name}.tar.gz")
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E tar czf "${archive}" --format=gnutar "${package_name}"
        WORKING_DIRECTORY "${staging_root}"
        RESULT_VARIABLE archive_result
    )
endif()
if(NOT archive_result EQUAL 0)
    message(FATAL_ERROR "failed to create release archive")
endif()

file(SHA256 "${archive}" checksum)
get_filename_component(archive_name "${archive}" NAME)
file(WRITE "${archive}.sha256" "${checksum}  ${archive_name}\n")
file(REMOVE_RECURSE "${staging_root}")
message(STATUS "Created ${archive}")
