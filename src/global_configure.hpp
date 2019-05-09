#ifndef CONFIGURE_HPP
#define CONFIGURE_HPP

#include <stdint.h>
#include <getopt.h>
#include <stdio.h>

#define ENABLE_DEBUG true

#define LONG_OPTIONS_TERMINAL {0,0,0,0}
#define IS_LONG_OPTIONS_TERMINAL(x) ( ((void*)(x).name == nullptr) && ((int)(x).has_arg == 0) && ((void*)(x).flag == nullptr) && ((int)(x).val == 0) )

// Exit Codes
typedef uint8_t ExitType;
enum GeneralExitCode : ExitType {
  SUCCESS            = 0,
  VERIFICATION_FAIL  = 255,
  COMMAND_LINE_ERROR = 254,
  ASSERTION_FAILURE  = 134 // TODO Change to non-literal (defined by system)
};

static const char* bar_flag = (char*) "--" ;

typedef struct struct_ArgsParameters {
  int argc;
  char** argv;
} ArgsParameters;

char* createShortOptionsFromLongOptions( const struct option* long_options );

ArgsParameters createProgramArgs( int argc, char** argv );

ArgsParameters createVariantArgs( int argc, char** argv );

void printGeneralInformationMessage( FILE* stream );
void printVariantInformationMessage( FILE* stream );
void printGeneralOptionsMessage( FILE* stream );
void printVariantOptionsMessage( FILE* stream );
void printGeneralExitCodeMessage( FILE* stream );
void printVariantExitCodeMessage( FILE* stream );
void printHelpMessage( FILE* stream );

#include <variant_configure.hpp>

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

static const struct option general_long_options[] = {
  {"size",      required_argument, nullptr, 's'},
  {"seed",      required_argument, nullptr, 'S'},
  {"x_size",    required_argument, nullptr, 'x'},
  {"y_size",    required_argument, nullptr, 'y'},
  {"z_size",    required_argument, nullptr, 'z'},
  {"timesteps", required_argument, nullptr, 'T'},
  {"epsilon",   required_argument, nullptr, 'e'},
  {"verify",    no_argument,       nullptr, 'V'}, // TODO Should this use the flagging process of getopt_long?
  {"help",      no_argument,       nullptr, 'h'},
  LONG_OPTIONS_TERMINAL
};

ProgramOptions parseProgramOptions( int argc, char** argv );

VariantOptions parseVariantOptions( int argc, char** argv );



#endif
