add_subdirectory(prtcls)

install(
  FILES 
    const_cp.hpp
    hydrostatic.hpp
    macros.hpp
    moist_air.hpp
    theta_dry.hpp
    theta_std.hpp
    units.hpp
    zip.hpp
  DESTINATION 
    include/libcloudph++/common
)
