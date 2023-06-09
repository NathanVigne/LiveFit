function(setup_vcpck_before_project)
    if(DEFINED ENV{VCPCK_ROOT})
        set(vcpkg_toolchain_path "$ENV{VCPCK_ROOT}/scripts/buildsystems/vcpkg.cmake")
        set(CMAKE_TOOLCHAIN_FILE "${vcpkg_toolchain_path}" CACHE STRING "" FORCE)
    endif()

    if(DEFINED ENV{VCPKG_DEFAULT_TRIPLET} AND NOT VCPK_TARGET_TRIPLET)
        set(VCPK_TARGET_TRIPLET "ENV{VCPKG_DEFAULT_TRIPLET}" CACHE STRING "")
    endif()
endfunction()


