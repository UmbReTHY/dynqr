include(ExternalProject)
find_package(Git REQUIRED)

ExternalProject_Add(
    catch_project
    PREFIX ${CMAKE_BINARY_DIR}/catch
    GIT_REPOSITORY https://github.com/philsquared/Catch.git
    TIMEOUT 10
    UPDATE_COMMAND ${GIT_EXECUTABLE} pull
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
   )

# Expose required variable (CATCH_INCLUDE_DIR) to parent scope
ExternalProject_Get_Property(catch_project source_dir)
add_library(catch INTERFACE)
target_include_directories(catch INTERFACE ${source_dir}/single_include)