cmake_minimum_required(VERSION 3.1)
project(Morpher)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)

find_package(LIBIGL REQUIRED QUIET)

#dqs
include_directories(${PROJECT_SOURCE_DIR}/dqs)
file(GLOB dqs ${PROJECT_SOURCE_DIR}/dqs/*)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Add your project files
file(GLOB SRCFILES *.cpp *.cc)
add_executable(${PROJECT_NAME} ${SRCFILES} ${dqs})
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw ${Geant4_LIBRARIES})
