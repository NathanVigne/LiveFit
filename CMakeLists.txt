set(ENV{VCPCK_ROOT} C:/dev/sdks/vcpkg)
set(ENV{VCPKG_DEFAULT_TRIPLET} x64-mingw-static)
set(ENV{GSL_ROOT_DIR} C:/dev/sdks/vcpkg/installed/x64-mingw-static/)

cmake_minimum_required(VERSION 3.14)
include(cmake/app_utils.cmake)
setup_vcpck_before_project()

project(testGSL LANGUAGES C CXX)

set(CMAKE_AUTOMOC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GSL REQUIRED)
find_package(Qt6 REQUIRED COMPONENTS Core Widgets Charts Gui)

add_executable(testGSL
  main.cpp
  mainwidget.h mainwidget.cpp
  worker.h worker.cpp
)

target_link_libraries(testGSL PRIVATE GSL::gsl GSL::gslcblas Qt6::Gui Qt6::Widgets Qt6::Core Qt6::Charts)


install(TARGETS testGSL
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
