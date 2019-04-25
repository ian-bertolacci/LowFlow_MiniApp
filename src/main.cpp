#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>

#include <omp.h>

#include <configure.hpp>
#include <types.hpp>
#include <macros.hpp>
#include <science.hpp>
#include <util.hpp>

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

typedef struct struct_ExperimentalResults {
  double elapsed;
  Variant_Domain* variant_domain;
  Variant_Grid* fp;
  Variant_Grid* vx;
  Variant_Grid* vy;
  Variant_Grid* vz;
} ExperimentalResults;

// Program arguments parser
ProgramOptions parseArguments( char** argv, int argc ){
  ProgramOptions opts = {
    .T = 100,
    .nx = 100,
    .ny = 100,
    .nz = 100,
    .epsilon = 0.00001,
    .seed = 123456789,
    .verify = false
  };

  char* options = (char*) "x:y:z:s:S:e:Vh";

  int c;
  opterr = 0;
  while( (c = getopt(argc, argv, options)) != -1){
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

bool verify( Basic_Domain* domain, ExperimentalResults results, double epsilon, unsigned int seed ){
  // Perform allocation in bulk, name later
  const int input_grid_count = 16;
  const int output_grid_count = 4;
  Basic_Grid* output_grids[output_grid_count];
  Basic_Grid* input_grids[input_grid_count];
  unsigned int seeds[input_grid_count];

  // Create seeds for input grids
  // use seeds in different loop avoids seed mixing.
  srand( seed );
  for( int i = 0; i < input_grid_count; i += 1 ){
    seeds[i] = rand();
  }
  // Allocate and populate input grids
  for( int i = 0; i < input_grid_count; i += 1 ){
    input_grids[i] = Basic_Grid_alloc( domain );
    Basic_Grid_populate_seeded( domain, input_grids[i], seeds[i] );
  }

  // Zero allocate output grids
  for( int i = 0; i < output_grid_count; i += 1 ){
    output_grids[i] = Basic_Grid_alloc( domain );
    Basic_Grid_populate_zero( domain, output_grids[i] );
  }

  // Give grids names
  Basic_Grid* fp = output_grids[0];
  Basic_Grid* vx = output_grids[1];
  Basic_Grid* vy = output_grids[2];
  Basic_Grid* vz = output_grids[3];

  Basic_Grid* dp         = input_grids[ 0];
  Basic_Grid* et         = input_grids[ 1];
  Basic_Grid* odp        = input_grids[ 2];
  Basic_Grid* opp        = input_grids[ 3];
  Basic_Grid* osp        = input_grids[ 4];
  Basic_Grid* permxp     = input_grids[ 5];
  Basic_Grid* permyp     = input_grids[ 6];
  Basic_Grid* permzp     = input_grids[ 7];
  Basic_Grid* pop        = input_grids[ 8];
  Basic_Grid* pp         = input_grids[ 9];
  Basic_Grid* rpp        = input_grids[10];
  Basic_Grid* sp         = input_grids[11];
  Basic_Grid* ss         = input_grids[12];
  Basic_Grid* z_mult_dat = input_grids[13];
  Basic_Grid* x_ssl_dat  = input_grids[14];
  Basic_Grid* y_ssl_dat  = input_grids[15];

  // Do baseline scientific kernel
  {
    // NlFunctionEval:261 analogue
    {
      int x, y, z;
      Basic_Domain_loop_interior(domain, x,y,z,
        {
          Basic_Grid_access(fp, x,y,z) =
           (  Basic_Grid_access(sp,  x,y,z)
            * Basic_Grid_access(dp,  x,y,z)
            - Basic_Grid_access(osp, x,y,z)
            * Basic_Grid_access(odp, x,y,z)
           )
           * Basic_Grid_access(pop, x,y,z)
           * Basic_Grid_access(z_mult_dat, x,y,z);
        }
      );
    }

    // NlFunctionEval:338 analogue
    {
      int x, y, z;
      Basic_Domain_loop_interior(domain, x,y,z,
        {
          Basic_Grid_access(fp, x,y,z) +=
              Basic_Grid_access(ss, x,y,z)
            * Basic_Grid_access(z_mult_dat, x,y,z)
            * (   Basic_Grid_access(pp, x,y,z)
                * Basic_Grid_access(sp, x,y,z)
                * Basic_Grid_access(dp, x,y,z)
                - Basic_Grid_access(opp, x,y,z)
                * Basic_Grid_access(osp, x,y,z)
                * Basic_Grid_access(odp, x,y,z)
              );
        }
      );
    }

    // NlFunctionEval:416 analogue
    {
      int x, y, z;

      Basic_Domain_loop_interior(domain, x,y,z,
        {
          Basic_Grid_access(fp, x,y,z) -=
              Basic_Grid_access(z_mult_dat, x,y,z)
            * (   Basic_Grid_access(sp, x,y,z)
                * Basic_Grid_access(et, x,y,z)
              );
        }
      );
    }

    // NlFunctionEval:551 analogue
    {
      int x, y, z;

      Basic_Domain_loop_interior(domain, x,y,z,
        {
          double x_dir_g   = ArithmeticMean( Basic_Grid_access( x_ssl_dat, x, y, 0), Basic_Grid_access( x_ssl_dat, x+1,  y, 0 ) );
          double x_dir_g_c = ArithmeticMean( Basic_Grid_access( x_ssl_dat, x, y, 0), Basic_Grid_access( x_ssl_dat, x+1,  y, 0 ) );
          double y_dir_g   = ArithmeticMean( Basic_Grid_access( y_ssl_dat, x, y, 0), Basic_Grid_access( y_ssl_dat, x,  y+1, 0 ) );
          double y_dir_g_c = ArithmeticMean( Basic_Grid_access( y_ssl_dat, x, y, 0), Basic_Grid_access( y_ssl_dat, x,  y+1, 0 ) );

          double diff_right = Basic_Grid_access(pp, x,y,z) - Basic_Grid_access(pp, x+1, y,   z);
          double diff_front = Basic_Grid_access(pp, x,y,z) - Basic_Grid_access(pp, x,   y+1, z);

          double updir_right = diff_right * x_dir_g_c - x_dir_g;
          double updir_front = diff_front * y_dir_g_c - y_dir_g;

          double sep = ArithmeticMean( Basic_Grid_access(z_mult_dat, x,y,z),  Basic_Grid_access(z_mult_dat, x, y, z+1) );

          double lower_cond =
            Basic_Grid_access(pp, x,y,z) / sep
          - (   Basic_Grid_access(z_mult_dat, x,y,z)
              / ( Basic_Grid_access(z_mult_dat, x,y,z)
                + Basic_Grid_access(z_mult_dat, x, y, z+1)
                )
            )
          * Basic_Grid_access(dp, x,y,z);

          double upper_cond =
            Basic_Grid_access(pp, x, y, z+1) / sep
          + (   Basic_Grid_access(z_mult_dat, x, y, z+1)
              / ( Basic_Grid_access(z_mult_dat, x,y,z)
                + Basic_Grid_access(z_mult_dat, x, y, z+1)
                )
            )
          * Basic_Grid_access(dp, x,y,z+1);

          double diff_upper = lower_cond - upper_cond;

          /*
          NOTE!
          Originally the harmonic mean was
          PMean( Basic_Grid_access(pp, x,y,z),
                 Basic_Grid_access(pp, x+1, y, z),
                 permxp[ip],
                 permxp[ip + 1]
          )
          However! PMean(a,b,c,d) in parflow is simply HarmonicMean(c, d)
          so the first two terms have been removed
          */

          double u_right =
            (
               Basic_Grid_access(z_mult_dat, x,y,z)
             * HarmonicMean( Basic_Grid_access(permxp, x,   y, z),
                             Basic_Grid_access(permxp, x+1, y, z) )
             * diff_right * x_dir_g_c
             * UpstreamMean(
                 updir_right,
                 0.0,
                 Basic_Grid_access(rpp, x,   y, z) * Basic_Grid_access(dp, x,   y, z),
                 Basic_Grid_access(rpp, x+1, y, z) * Basic_Grid_access(dp, x+1, y, z)
               )
           ) + (
              Basic_Grid_access(z_mult_dat, x,y,z)
             * HarmonicMean( Basic_Grid_access(permxp, x,   y, z),
                             Basic_Grid_access(permxp, x+1, y, z) )
             * (-x_dir_g)
             * UpstreamMean(
                 updir_right,
                 0.0,
                 Basic_Grid_access(rpp, x,   y, z) * Basic_Grid_access(dp, x,   y, z),
                 Basic_Grid_access(rpp, x+1, y, z) * Basic_Grid_access(dp, x+1, y, z)
               )
           );

           double u_front =
             (
                Basic_Grid_access(z_mult_dat, x,y,z)
              * HarmonicMean( Basic_Grid_access(permyp, x,   y, z),
                              Basic_Grid_access(permyp, x+1, y, z) )
              * diff_front * x_dir_g_c
              * UpstreamMean(
                  updir_front,
                  0.0,
                  Basic_Grid_access(rpp, x,   y, z) * Basic_Grid_access(dp, x,   y, z),
                  Basic_Grid_access(rpp, x+1, y, z) * Basic_Grid_access(dp, x+1, y, z)
                )
            ) + (
               Basic_Grid_access(z_mult_dat, x,y,z)
              * HarmonicMean( Basic_Grid_access(permyp, x,   y, z),
                              Basic_Grid_access(permyp, x+1, y, z) )
              * (-x_dir_g)
              * UpstreamMean(
                  updir_front,
                  0.0,
                  Basic_Grid_access(rpp, x,   y, z) * Basic_Grid_access(dp, x,   y, z),
                  Basic_Grid_access(rpp, x+1, y, z) * Basic_Grid_access(dp, x+1, y, z)
                )
            );

          double u_upper =
                      HarmonicMeanDZ(
                          Basic_Grid_access(permzp,     x, y, z  ),
                          Basic_Grid_access(permzp,     x, y, z+1),
                          Basic_Grid_access(z_mult_dat, x, y, z  ),
                          Basic_Grid_access(z_mult_dat, x, y, z+1)
                      )
                    * diff_upper
                    * UpstreamMean(
                        lower_cond,
                        upper_cond,
                        Basic_Grid_access(rpp, x, y, z  ) * Basic_Grid_access(dp, x, y, z  ),
                        Basic_Grid_access(rpp, x, y, z+1) * Basic_Grid_access(dp, x, y, z+1)
                      );

          Basic_Grid_access(vx, x,y,z) = u_right;
          Basic_Grid_access(vy, x,y,z) = u_front;
          Basic_Grid_access(vz, x,y,z) = u_upper;

          Basic_Grid_access(fp, x  , y  , z  ) += u_right * u_front * u_upper;
          Basic_Grid_access(fp, x+1, y  , z  ) += u_right;
          Basic_Grid_access(fp, x  , y+1, z  ) -= u_front;
          Basic_Grid_access(fp, x  , y  , z+1) -= u_upper;
        }
      );
    }

    // NlFunctionEval:747 analogue placholder
    // TODO Loop NlFunctionEval:747
    /*
    {

    }
    */
    // old jacobi placeholder
    // Basic_Domain_loop_interior(domain, x,y,z,
    //   {
    //     int idx_000  = Basic_Domain_idx(domain, x,   y,   z  );
    //     int idx_p100 = Basic_Domain_idx(domain, x+1, y,   z  );
    //     int idx_p010 = Basic_Domain_idx(domain, x,   y+1, z  );
    //     int idx_p001 = Basic_Domain_idx(domain, x,   y,   z+1);
    //     int idx_n100 = Basic_Domain_idx(domain, x-1, y,   z  );
    //     int idx_n010 = Basic_Domain_idx(domain, x,   y-1, z  );
    //     int idx_n001 = Basic_Domain_idx(domain, x,   y,   z-1);
    //
    //     output_data[idx_000] = (1.0/7.0)*(
    //       source_data[idx_000] +
    //       source_data[idx_p100] + source_data[idx_p010] + source_data[idx_p001] +
    //       source_data[idx_n100] + source_data[idx_n010] + source_data[idx_n001]
    //     );
    //   }
    // );
  }

  char* grid_names[output_grid_count] = { (char*)"fp", (char*)"vx", (char*)"vy", (char*)"vz" };
  Variant_Grid*  experimental_grids[output_grid_count] = { results.fp, results.vx, results.vy, results.vz };

  // Compare each experimental grid
  bool passed = true;
  for( int i = 0; i < output_grid_count; i += 1 ){
    printf( "Verifying experimental %s: ", grid_names[i] );
    bool this_passed = true;
    int x,y,z;
    Basic_Domain_loop_whole(domain, x,y,z,
      {
        double compare_value = Variant_Grid_access(experimental_grids[i], x,y,z);
        double baseline_value = Basic_Grid_access(output_grids[i], x,y,z);
        double delta = fabs( compare_value - baseline_value );
        if( delta >= epsilon ){
          if( this_passed ) printf("Failed!\n(x, y, z): | Comparison - Baseline | = absolute delta >= epsilon\n");
          printf("(%d, %d, %d): |%15.7f - %15.7f| = %-10.7f >= %-10.5f \n",
            x, y, z,
            compare_value,
            baseline_value,
            delta,
            epsilon
          );
          this_passed = false;
        }
      }
    );
    passed &= this_passed;
    if( this_passed ) printf("Passed!\n");
  }

  if( passed ){
    printf("All Passed!\n");
  }

  // Deallocate verification input grids
  for( int i = 0; i < input_grid_count; i += 1 ){
    Basic_Grid_dealloc( input_grids[i] );
  }

  // Deallocate verification output grids
  for( int i = 0; i < output_grid_count; i += 1 ){
    Basic_Grid_dealloc( output_grids[i] );
  }

  return passed;
}

ExperimentalResults experiment( Basic_Domain* domain, unsigned int seed ){
  // Perform allocation in bulk, name later
  const int input_grid_count = 16;
  const int output_grid_count = 4;
  Variant_Grid* output_grids[output_grid_count];
  Variant_Grid* input_grids[input_grid_count];
  unsigned int seeds[input_grid_count];

  Variant_Domain* variant_domain = Variant_Domain_alloc( Basic_Domain_nx(domain), Basic_Domain_ny(domain),Basic_Domain_nz(domain) );

  // Create seeds for input grids
  // use seeds in different loop avoids seed mixing.
  srand( seed );
  for( int i = 0; i < input_grid_count; i += 1 ){
    seeds[i] = rand();
  }
  // Allocate and populate input grids
  for( int i = 0; i < input_grid_count; i += 1 ){
    input_grids[i] = Variant_Grid_alloc( domain );
    Variant_Grid_populate_seeded( variant_domain, input_grids[i], seeds[i] );
  }

  // Zero allocate output grids
  for( int i = 0; i < output_grid_count; i += 1 ){
    output_grids[i] = Variant_Grid_alloc( variant_domain );
    Variant_Grid_populate_zero( variant_domain, output_grids[i] );
  }

  // Give grids names
  // Output Grids
  Variant_Grid* fp = output_grids[0];
  Variant_Grid* vx = output_grids[1];
  Variant_Grid* vy = output_grids[2];
  Variant_Grid* vz = output_grids[3];

  // Input Grids
  Variant_Grid* dp         = input_grids[ 0];
  Variant_Grid* et         = input_grids[ 1];
  Variant_Grid* odp        = input_grids[ 2];
  Variant_Grid* opp        = input_grids[ 3];
  Variant_Grid* osp        = input_grids[ 4];
  Variant_Grid* permxp     = input_grids[ 5];
  Variant_Grid* permyp     = input_grids[ 6];
  Variant_Grid* permzp     = input_grids[ 7];
  Variant_Grid* pop        = input_grids[ 8];
  Variant_Grid* pp         = input_grids[ 9];
  Variant_Grid* rpp        = input_grids[10];
  Variant_Grid* sp         = input_grids[11];
  Variant_Grid* ss         = input_grids[12];
  Variant_Grid* z_mult_dat = input_grids[13];
  Variant_Grid* x_ssl_dat  = input_grids[14];
  Variant_Grid* y_ssl_dat  = input_grids[15];

  // Time and Perform science
  double start = omp_get_wtime();
  science( domain, fp, vx, vy, vz, dp, et, odp, opp, osp, permxp, permyp, permzp, pop, pp, rpp, sp, ss, z_mult_dat, x_ssl_dat, y_ssl_dat );
  double end = omp_get_wtime();

  // Deallocate input grids
  for( int i = 0; i < input_grid_count; i += 1 ){
    Variant_Grid_dealloc( input_grids[i] );
  }

  // Return results
  ExperimentalResults results = {
    .elapsed = end - start,
    .variant_domain = variant_domain,
    .fp = fp,
    .vx = vx,
    .vy = vy,
    .vz = vz
  };

  return results;
}


int main( int argc, char** argv ){

  ProgramOptions opts = parseArguments( argv, argc );

  printf( "(nx, ny, nz): (%d, %d, %d)\n", opts.nx, opts.ny, opts.nz );

  Basic_Domain* base_domain = Basic_Domain_alloc( opts.nx, opts.ny, opts.nz );

  ExperimentalResults results = experiment( base_domain, opts.seed );

  printf( "Elapsed: %fs\n", results.elapsed );

  bool passed = true;
  if( opts.verify ){
    passed = verify( base_domain, results, opts.epsilon, opts.seed );
  }

  Basic_Domain_dealloc( base_domain );
  Variant_Grid_dealloc( results.fp );
  Variant_Grid_dealloc( results.vx );
  Variant_Grid_dealloc( results.vy );
  Variant_Grid_dealloc( results.vz );
  Variant_Domain_dealloc( results.variant_domain );

  return (!passed)?VERIFICATION_FAIL:SUCCESS;
}
