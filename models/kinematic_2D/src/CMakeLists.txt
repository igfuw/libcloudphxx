# std::future
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("#include <future>\n int main() { std::future<void> f; }\n" STD_FUTURE_WORKS)
if(${STD_FUTURE_WORKS})
  add_definitions(-DSTD_FUTURE_WORKS)
endif()

add_executable(icicle icicle.cpp)

# libmpdata requires std_14, but sets it with a CXX_FLAG
# which then gets overriden by libcloudphxx std_11 compile_feature
# hence we require std_14 feature here
# TODO: add an exported target to libmpdata++ that requires std_14 through compile_feaures
#       and not through a CXX_FLAG
target_compile_features(icicle PRIVATE cxx_std_14)

target_link_libraries(icicle PRIVATE ${libmpdataxx_LIBRARIES})
target_include_directories(icicle PRIVATE ${libmpdataxx_INCLUDE_DIRS})
target_link_libraries(icicle PRIVATE clphxx::cloudphxx_lgrngn)
find_package(Boost COMPONENTS program_options REQUIRED)
target_link_libraries(icicle PRIVATE ${Boost_LIBRARIES})

install(TARGETS icicle DESTINATION bin)
