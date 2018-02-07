# std::future
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("#include <future>\n int main() { std::future<void> f; }\n" STD_FUTURE_WORKS)
if(${STD_FUTURE_WORKS})
  add_definitions(-DSTD_FUTURE_WORKS)
endif()

add_executable(icicle icicle.cpp)

target_link_libraries(icicle ${libmpdataxx_LIBRARIES})
target_include_directories(icicle PUBLIC ${libmpdataxx_INCLUDE_DIRS})

target_link_libraries(icicle ${libcloudphxx_LIBRARIES})

find_package(Boost COMPONENTS program_options REQUIRED)
target_link_libraries(icicle ${Boost_LIBRARIES})

install(TARGETS icicle DESTINATION bin)
