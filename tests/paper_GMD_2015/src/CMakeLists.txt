find_package(libmpdata++)

if (NOT libmpdataxx_FOUND)
  message(FATAL_ERROR "libmpdata++ not found")
endif()

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  # make the Release flags the default
  set(CMAKE_CXX_FLAGS ${libmpdataxx_CXX_FLAGS_RELEASE})
  # ignore CMake defaults
  set(CMAKE_CXX_FLAGS_RELEASE "")
else()
  set(CMAKE_CXX_FLAGS_DEBUG ${libmpdataxx_CXX_FLAGS_DEBUG})
endif()

# std::future
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("#include <future>\n int main() { std::future<void> f; }\n" STD_FUTURE_WORKS)
if(${STD_FUTURE_WORKS})
  add_definitions(-DSTD_FUTURE_WORKS)
endif()

add_executable(icicle icicle.cpp)

target_link_libraries(icicle ${libmpdataxx_LIBRARIES})
target_include_directories(icicle PUBLIC ${libmpdataxx_INCLUDE_DIRS})

find_package(Boost COMPONENTS thread iostreams system timer program_options filesystem REQUIRED)
target_link_libraries(icicle ${Boost_LIBRARIES})

# search for a local version of libcloudphxx
find_library(libcloudphxx_LIBRARIES cloudphxx_lgrngn ${CMAKE_SOURCE_DIR}/../../build/src NO_DEFAULT_PATH)
# if no local version found, look for an installed one
find_library(libcloudphxx_LIBRARIES cloudphxx_lgrngn)

if (NOT libcloudphxx_LIBRARIES)
  message(FATAL_ERROR "libclouphxx++ not found")
endif()

target_include_directories(icicle PUBLIC ${CMAKE_SOURCE_DIR}/../../include)
target_link_libraries(icicle ${libcloudphxx_LIBRARIES})

install(TARGETS icicle DESTINATION bin)
