#include <configure.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <ctype.h>
#include <stdbool.h>

ArgsParameters createProgramArgs( int argc, char** argv ){
  // find the bar ("--") flag
  int program_argc = 0;
  for( /*program_argc = 0*/;
       program_argc < argc && strcmp(argv[program_argc], bar_flag) != 0;
       program_argc += 1
     ){ /* do nothing */ }
     // Note: This does not include the bar flag (Trust me -Ian)

   ArgsParameters args {
     .argc = program_argc,
     .argv = argv
   };

  return args;
}

ArgsParameters createVariantArgs( int argc, char** argv ){

  // find the bar ("--") flag
  int program_argc = 0;
  for( /*program_argc = 0*/;
       program_argc < argc && strcmp(argv[program_argc], bar_flag) != 0;
       program_argc += 1
     ){ /* do nothing */ }

  // conditional values if there is a bar flag
  bool more_args = argc > program_argc;
  ArgsParameters args = {
    .argc = (more_args)? argc - (program_argc): 0,
    .argv = (more_args)? &argv[program_argc] : nullptr
  };

  return args;
}

char* createShortOptionsFromLongOptions( const struct option* long_options ){

  // Count how many characters are required
  int count = 0;
  for( int i = 0; ! IS_LONG_OPTIONS_TERMINAL( long_options[i] ); ++i ){
    if( long_options[i].flag != nullptr ) continue;
    count += 1;
    if( long_options[i].has_arg == required_argument ) count += 1;
  }

  char* short_options = (char*) calloc(sizeof(char), count+1);

  char* char_pos = short_options;
  for( int i = 0; ! IS_LONG_OPTIONS_TERMINAL( long_options[i] ); ++i ){
    if( long_options[i].flag != nullptr ) continue;
    *char_pos = (char) long_options[i].val;
    ++char_pos;
    if( long_options[i].has_arg == required_argument ){
      *char_pos = ':';
      ++char_pos;
    }
  }
  *char_pos = '\0';

  return short_options;
}

void printGeneralInformationMessage( FILE* stream ){
  /* Nothing */
}

void printGeneralOptionsMessage( FILE* stream ){
  fprintf( stream,
    "General Program Command Line Options:\n"
    "  --verify\n"
    "   -V\n"
    "    Perform verification after experiment to check that experimental output is correct.\n"
    "    Program returns VERIFICATION_FAIL if verification fails.\n"
    "  --x_size <N : uint>\n"
    "   -x      <N : uint>\n"
    "    Set size of domain in x direction.\n"
    "  --y_size <N : uint>\n"
    "   -y      <N : uint>\n"
    "    Set size of domain in y direction.\n"
    "  --z_size <N : int>\n"
    "   -z      <N : uint>\n"
    "    Set size of domain in z direction.\n"
    "  --size <N : uint>\n"
    "   -s    <N : uint>\n"
    "    Set size of domain in x, y, z directions to the same value.\n"
    "    --size 10 is equivalent to --x_size 10 --y_size 10 --z_size 10.\n"
    "  --epsilon <epsilon : float>\n"
    "   -e       <epsilon : float>\n"
    "    Set acceptable error bound used for verification.\n"
    "    An absolute difference between experimental and verification output \n"
    "    grids less than epsilon is considered 'passing'\n"
    "  --seed <N : uint>\n"
    "   -S    <N : uint>\n"
    "    Set the seed used to generate grid-seeds.\n"
    "  --help\n"
    "   -h\n"
    "    Print this message.\n"
  );
}

void printGeneralExitCodeMessage( FILE* stream ){
  fprintf( stream,
    "General Program Exit Codes:\n"
    "  SUCCESS:            %u\n"
    "  VERIFICATION_FAIL:  %u\n"
    "  COMMAND_LINE_ERROR: %u\n"
    "  ASSERTION_FAILURE:  %u\n",
    SUCCESS,
    VERIFICATION_FAIL,
    COMMAND_LINE_ERROR,
    ASSERTION_FAILURE
  );
}

void printHelpMessage( FILE* stream ){
  printGeneralInformationMessage( stream );
  printVariantInformationMessage( stream );
  printGeneralOptionsMessage( stream );
  printVariantOptionsMessage( stream );
  printGeneralExitCodeMessage( stream );
  printVariantExitCodeMessage( stream );
}

// Program arguments parser
ProgramOptions parseProgramOptions( int argc, char** argv ){
  ProgramOptions opts = ProgramOptions_Defaults;

  ArgsParameters args = createProgramArgs( argc, argv );

  char* short_options = createShortOptionsFromLongOptions( general_long_options );

  char c;
  optind = 1;
  opterr = 0;
  int option_index = 0;

  while( (c = getopt_long(args.argc, args.argv, short_options, general_long_options, &option_index)) != -1 ){
    switch (c){
      case 's':
      {
        int size = atoi( optarg );
        opts.nx = size;
        opts.ny = size;
        opts.nz = size;
        break;
      }

      case 'S':
        opts.seed = strtoul( optarg, NULL, 10 );
        break;

      case 'x':
        opts.nx = atoi( optarg );
        break;

      case 'y':
        opts.ny = atoi( optarg );
        break;

      case 'z':
        opts.nz = atoi( optarg );
        break;

      case 'T':
        opts.T = atoi( optarg );
        break;

      case 'e':
        opts.epsilon = strtod( optarg, NULL );
        break;

      case 'V':
        opts.verify = true;
        break;

      case 'h':
        printHelpMessage( stdout );
        exit(SUCCESS);

      default:
        fprintf(stderr, "Error while parsing general command line arguments! (case %c, 0x%x)\n", c, c);
        printHelpMessage( stderr );
        exit(COMMAND_LINE_ERROR);
    }
  }

  free( short_options );

  // Set the variant's options.
  // Note: Pass full arguments to parseVariantOptions
  opts.variant_options = parseVariantOptions( argc, argv );

  return opts;
}
