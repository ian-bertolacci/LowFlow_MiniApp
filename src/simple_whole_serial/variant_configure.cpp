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

VariantOptions parseVariantOptions( int argc, char** argv ){
  VariantOptions opts = VariantOptions_Default;

  ArgsParameters args = createVariantArgs( argc, argv );

  char* options = (char*) "h";

  // if there are no flags do not parse.
  if( args.argc > 0 ){
    // Parse arguments
    char c;
    optind = 1; // reset getopt process
    opterr = 0;
    while( (c = getopt(args.argc, args.argv, options)) != -1){
      switch (c){
        case 'h':
          printf( "The simple whole serial variant has no arguments.\n" );
          exit(SUCCESS);

        case '?':
        {
          char* position = strchr( options, optopt );
          // if option not in options string
          if( position == NULL ){
            // if printable
            if( isprint(optopt) ) fprintf(stderr, "Unknown option -%c'.\n", optopt);
            // if not printable
            else fprintf(stderr, "Unknown option character 0x%x'.\n", optopt);
            exit(COMMAND_LINE_ERROR);
          }
          // if option in string...
          // if next character is the options modifier (':')
          else if( position[1] == ':' ){
            fprintf(stderr, "Option -%c requires an argument.\n", optopt);
            exit(COMMAND_LINE_ERROR);
          }
          // Then theres been some kind of parsing error,
        }

        default:
        {
          fprintf(stderr, "Error while parsing command line arguments, case (%x: '%c')\n", c, c);
          exit(COMMAND_LINE_ERROR);
        }

      }
    }
  }

  return opts;
}
