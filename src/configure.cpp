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

static const ProgramOptions ProgramOptions_Defaults {
  .T = 100,
  .nx = 100,
  .ny = 100,
  .nz = 100,
  .epsilon = 0.00001,
  .seed = 123456789,
  .verify = false
};

// Program arguments parser
ProgramOptions parseArguments( char** argv, int argc ){
  ProgramOptions opts = ProgramOptions_Defaults;

  char* options = (char*) "x:y:z:s:S:e:Vh-";

  char c;
  opterr = 0;
  bool continue_loop = true;

  while( continue_loop && (c = getopt(argc, argv, options)) != -1){
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

      case '-':
        continue_loop = false;
        break;

      case 'h':
        printf(
          "Command Line Options:\n"
          "  -V : perform verification after experiment to check that experimental output is correct.\n"
          "       Program returns 255 if verification fails.\n"
          "  -x <N : int> : Set size of domain in x direction.\n"
          "  -y <N : int> : Set size of domain in y direction.\n"
          "  -z <N : int> : Set size of domain in z direction.\n"
          "  -s <N : int> : Set size of domain in x, y, z directions to the same value.\n"
          "  -s 10 is equivalent to -x 10 -y 10 -z 10.\n"
          "  -e <epsilon : float> : Set acceptable error bound used for verification. An absolute difference between\n"
          "                         experimental and verification output grids less than epsilon is considered 'passing'\n"
          "  -S <N : uint> : Set the seed used to generate grid-seeds.\n"
          "  -h : Print this message.\n"
          "\nExit Codes:\n"
          "  SUCCESS: %u\n"
          "  VERIFICATION_FAIL: %u\n"
          "  COMMAND_LINE_ERROR: %u\n"
          "  ASSERTION_FAILURE: %u\n",
          SUCCESS,
          VERIFICATION_FAIL,
          COMMAND_LINE_ERROR,
          ASSERTION_FAILURE
        );
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

  return opts;
}
