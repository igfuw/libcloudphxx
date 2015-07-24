// Python binding for libcloudph++
//
// author[s]: Sylwester Arabas, Dorota Jarecka
// licensing: GPU GPL v3
// copyright: University of Warsaw

// all libcloudph++'s includes which potentially include <cassert> or <assert.h>
// needed here as assert.h redefines assert() every time it is included
#include <libcloudph++/blk_1m/extincl.hpp>
#include <libcloudph++/blk_2m/extincl.hpp>
#include <libcloudph++/lgrngn/extincl.hpp> // TODO: still missing some includes as asserts give aborts :(

// turning asserts into exceptions
#undef assert
#define assert_str(s) #s
// TODO: for some reason it's not working everywhere 
//#if defined(__GNUC_PREREQ) && __GNUC_PREREQ(2,6)
//#  define assert_fun __PRETTY_FUNCTION__
//#else
#  define assert_fun __func__
//#endif
#define assert(cond) { \
  if (!(cond)) std::cerr <<  \
    __FILE__ ":" + std::to_string(__LINE__) + ": " + assert_fun + ": " \
    "Assertion `" assert_str(cond) "' failed." \
  ; \
}
