#ifndef CONFIGURE_HPP
#define CONFIGURE_HPP

#define ENABLE_DEBUG true

#include <variant_configure.hpp>
#include <stdint.h>

// Exit Codes
enum ExitCode : uint8_t {
  SUCCESS = 0,
  VERIFICATION_FAIL = 255,
  COMMAND_LINE_ERROR = 254,
  ASSERTION_FAILURE = 134 // TODO Change to non-literal
};

// Program Options struct
typedef struct struct_ProgramOptions {
  int T;
  int nx;
  int ny;
  int nz;
  double epsilon;
  unsigned int seed;
  bool verify;
} ProgramOptions;

ProgramOptions parseArguments( char** argv, int argc );

#endif
