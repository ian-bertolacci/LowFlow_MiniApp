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

void printVariantInformationMessage( FILE* stream ){
  fprintf( stream, "This is the serial implementation of the LowFlow mini-app\n" );
}

void printVariantOptionsMessage( FILE* stream ){
  fprintf( stream,
    "Variant Specific Command Line Options:\n"
    "  --help\n"
    "   -h\n"
    "    Print this message.\n"
  );
}

void printVariantExitCodeMessage( FILE* stream ){
  /* Nothing */
}

VariantOptions parseVariantOptions( int argc, char** argv ){
  VariantOptions opts = VariantOptions_Default;
  ArgsParameters args = createVariantArgs( argc, argv );

  // if there are no flags do not parse.
  if( args.argc > 0 ){
    // create short options string
    char* short_options = createShortOptionsFromLongOptions( variant_long_options );

    // Parse arguments
    char c;
    int option_index = 0;
    while( (c = getopt_long(args.argc, args.argv, short_options, variant_long_options, &option_index)) != -1 ){

      switch (c){
        case 'h':
          printHelpMessage( stdout );
          exit(SUCCESS);

        default:
          fprintf(stderr, "Error while parsing command line arguments for variant\n");
          printHelpMessage( stdout );
          exit(COMMAND_LINE_ERROR);
      }
    }
    free( short_options );
  }

  return opts;
}
