cmake_minimum_required(VERSION 3.16)

project(FDenoise VERSION 0.0.1 LANGUAGES CXX)

add_compile_definitions(${PROJECT_NAME} PRIVATE
    TARGET_NAME="${PROJECT_NAME}"
    TARGET_VER="${PROJECT_VERSION}"
)

file(GLOB_RECURSE ${PROJECT_NAME}_LIBS
    "${CMAKE_CURRENT_SOURCE_DIR}/libs/*.dll"
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/.install")

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt6 COMPONENTS Core Concurrent REQUIRED)

set(LIBS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/libs")

set(FFTW_INCLUDE_DIR "${LIBS_ROOT}/fftw3")
set(FFTW_LIBRARY "${LIBS_ROOT}/fftw3/libfftw3-3.lib")
set(FFTW_DLL "${LIBS_ROOT}/fftw3/libfftw3-3.dll")

set(SNDFILE_INCLUDE_DIR "${LIBS_ROOT}/libsnd/include")
set(SNDFILE_LIBRARY "${LIBS_ROOT}/libsnd/lib/sndfile.lib")
set(SNDFILE_DLL "${LIBS_ROOT}/libsnd/bin/sndfile.dll")

add_executable(${PROJECT_NAME} main.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE
    ${FFTW_INCLUDE_DIR}
    ${SNDFILE_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME} PRIVATE
    Qt6::Core
    Qt6::Concurrent
    ${FFTW_LIBRARY}
    ${SNDFILE_LIBRARY}
)

if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /wd4996 /wd4244)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -w -fpermissive)
endif()

install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION .
)

if(WIN32)
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/libs/libsnd/bin/sndfile.dll"
            DESTINATION .
    )
    install(FILES ${${PROJECT_NAME}_LIBS}
        DESTINATION .
    )
endif()

qt_generate_deploy_script(
  OUTPUT_SCRIPT deploy_script
  TARGET ${PROJECT_NAME}
  CONTENT "
    qt_deploy_runtime_dependencies(
      EXECUTABLE $<TARGET_FILE:${PROJECT_NAME}>
      BIN_DIR .
      LIB_DIR .
      PLUGINS_DIR .
      NO_TRANSLATIONS
    )
  "
)
install(SCRIPT ${deploy_script})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} --install . --config $<CONFIG>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
