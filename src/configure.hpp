#ifndef CONFIGURE_HPP
#define CONFIGURE_HPP

#define ENABLE_DEBUG true

#define LONG_OPTIONS_TERMINAL {0,0,0,0}

#include <variant_configure.hpp>
#include <stdint.h>

// Exit Codes
enum ExitCode : uint8_t {
  SUCCESS = 0,
  VERIFICATION_FAIL = 255,
  COMMAND_LINE_ERROR = 254,
  ASSERTION_FAILURE = 134 // TODO Change to non-literal
};

static const char* bar_flag = (char*) "--" ;

typedef struct struct_ArgsParameters {
  int argc;
  char** argv;
} ArgsParameters;

// Program Options struct
typedef struct struct_ProgramOptions {
  int T;
  int nx;
  int ny;
  int nz;
  double epsilon;
  unsigned int seed;
  bool verify;
  VariantOptions variant_options;
} ProgramOptions;

static const ProgramOptions ProgramOptions_Defaults {
  .T = 100,
  .nx = 100,
  .ny = 100,
  .nz = 100,
  .epsilon = 0.00001,
  .seed = 123456789,
  .verify = false,
  .variant_options = VariantOptions_Default
};

ArgsParameters createProgramArgs( int argc, char** argv );

ArgsParameters createVariantArgs( int argc, char** argv );

ProgramOptions parseProgramOptions( int argc, char** argv );

VariantOptions parseVariantOptions( int argc, char** argv );

#endif
