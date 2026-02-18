#
# Check that the compiler can compile the flags specified
# for the compiler.
# This should be a first check for the compiler
# to ensure that the compiler will handle all cases.
#
# This module only exposes the SIESTA_HAS_IPO variable
# to signal whether it can be used or not.

include(CheckCCompilerFlag)
include(CheckCSourceRuns)

include(CheckCXXCompilerFlag)
include(CheckCXXSourceRuns)

include(CheckFortranCompilerFlag)
include(CheckFortranSourceRuns)

include(CheckIPOSupported)

get_property(project_languages GLOBAL PROPERTY ENABLED_LANGUAGES)

# Filter out the C, CXX and Fortran languages
list(FILTER project_languages INCLUDE REGEX "^(CXX|C|Fortran)$")


# When we are doing cross-compilation, we have to specify the exact run-time
# flags, and we have to bypass the `try_run` runs by manually setting
# the return values.
# Manually specify flags for the try-runs, otherwise it'll fail
# This should be called like this:
#  siesta_try_run_default(run_test 0)
#  try_run(run_test ...)
# Which will then setup the run_test to a cache variable
# when cross compiling.
macro(siesta_try_run_default var)
  if( ${ARGC} GREATER 0 )
    set(${var} ${ARGV0} CACHE INTERNAL "Default value for the variable ${var}")
  else()
    set(${var} 0 CACHE INTERNAL "Default value for the variable ${var}")
  endif()
endmacro()


set(compiler_ok TRUE)
# Start checks
message(CHECK_START "Checking default compiler usability and flags.")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

# Now check details
set(SIESTA_HAS_IPO TRUE)
foreach(lang IN LISTS project_languages)
  message(CHECK_START "checking ${lang} flags")
  list(APPEND CMAKE_MESSAGE_INDENT "  ")

  set(compiler_${lang}_ok TRUE)

  # just check if the current flags are *good*
  cmake_language(CALL "check_${lang}_compiler_flag" "" flags_${lang}_default)
  if(NOT flags_${lang}_default)
    message(STATUS "Default flags for ${lang} failed, please check your default flags")
    set(compiler_${lang}_ok FALSE)
  endif()

  # check if a simple program can be compiled
  string(TOLOWER ${lang} lang_lower)
  siesta_try_run_default(runs_${lang}_default TRUE)

  if(lang_lower MATCHES "fortran")
    check_fortran_source_runs([==[
program t
print *, "hello world"
end
]==] runs_${lang}_default)
  elseif(lang_lower MATCHES "c")
    check_c_source_runs([==[
#include <stdio.h>
int main() {
   printf("Hello, World!");
   return 0;
}]==] runs_${lang}_default)
  elseif(lang_lower MATCHES "cxx")
    check_cxx_source_runs([==[
#include <iostream>
int main() {
    std::cout << "Hello World!";
    return 0;
}
]==] runs_${lang}_default)

  endif()

  if(NOT runs_${lang}_default)
    message(STATUS "Default program for ${lang} failed, please check your compiler")
    set(compiler_${lang}_ok FALSE)
  endif()

  # Store in the SIESTA_<lang>_HAS_IPO variable whether it can be enabled
  check_ipo_supported(RESULT SIESTA_${lang}_HAS_IPO LANGUAGES ${lang})
  if(NOT SIESTA_${lang}_HAS_IPO)
    # This is not fatal
    message(STATUS "IPO for ${lang} is not available")
    set(SIESTA_HAS_IPO FALSE)
  endif()

  if(compiler_${lang}_ok)
    message(CHECK_PASS "passed")
  else()
    message(CHECK_FAIL "failed")
    set(compiler_ok FALSE)
  endif()
  list(POP_BACK CMAKE_MESSAGE_INDENT)

endforeach()

if(SIESTA_HAS_IPO)
  message(STATUS "Can use IPO for all enabled languages")
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)
if(compiler_ok)
  message(CHECK_PASS "passed")
else()
  message(CHECK_FAIL "failed")
  message(FATAL_ERROR [==[
The compilers for the specified languages are setup incorrectly!

Please correct the compiler and flags specified!
]==])
endif()
