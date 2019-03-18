#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>

#include <omp.h>

#include <types.hpp>
#include <macros.hpp>
#include <util.hpp>

#define ENABLE_DEBUG true

// Program Options struct
typedef struct struct_ProgramOptions {
  int T;
  int nx;
  int ny;
  int nz;
  double epsilon;
  bool verify;
} ProgramOptions;

// Program arguments parser
ProgramOptions parseArguments( char** argv, int argc ){
  ProgramOptions opts;
  opts.nx = 100;
  opts.ny = 100;
  opts.nz = 100;
  opts.T = 100;
  opts.verify = false;
  opts.epsilon = 0.00001;

  int c;
  opterr = 0;

  while( (c = getopt (argc, argv, "x:y:z:s:T:e:v")) != -1){
    switch (c){
      case 's':
      {
        int size = atoi( optarg );
        opts.nx = size;
        opts.ny = size;
        opts.nz = size;
        break;
      }

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

      case 'v':
        opts.verify = true;
        break;

      case '?':
        printf( "?\n" );
        if( optopt == 'x' || optopt == 'y' || optopt == 'z' || optopt == 'T' )
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);

        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);

        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);

        abort();
      default:
        printf( "default\n" );
        abort ();
    }
  }

  return opts;
}


bool verify( Basic_Domain* domain, Variant_Grid* compare, int T, bool verbose, double epsilon ){
  if( verbose ) printf("Verifying: ");

  // allocate ping-pong grids
  Basic_Grid* source = Basic_Grid_alloc( domain );
  Basic_Grid* target = Basic_Grid_alloc( domain );
  Basic_Grid* result = source;
  int x, y, z;

  Basic_Grid_populate( domain, source );

  double start = omp_get_wtime();

  for( int t = 0; t < T; t += 1 ){
    double* source_data = Variant_Grid_data(source);
    double* target_data = Variant_Grid_data(target);
    Basic_Domain_loop_interior(domain, x,y,z,
      {
          int idx_000  = Basic_Domain_idx(domain, x,   y,   z  );
          int idx_p100 = Basic_Domain_idx(domain, x+1, y,   z  );
          int idx_p010 = Basic_Domain_idx(domain, x,   y+1, z  );
          int idx_p001 = Basic_Domain_idx(domain, x,   y,   z+1);
          int idx_n100 = Basic_Domain_idx(domain, x-1, y,   z  );
          int idx_n010 = Basic_Domain_idx(domain, x,   y-1, z  );
          int idx_n001 = Basic_Domain_idx(domain, x,   y,   z-1);

          target_data[idx_000] = (1.0/7.0)*(
            source_data[idx_000] +
            source_data[idx_p100] + source_data[idx_p010] + source_data[idx_p001] +
            source_data[idx_n100] + source_data[idx_n010] + source_data[idx_n001]
          );
      }
    );
    result = target;
    target = source;
    source = result;
  }

  bool passed = true;

  double* compare_data = Variant_Grid_data( compare );
  double* result_data = Basic_Grid_data( result );

  Basic_Domain_loop_whole(domain, x,y,z,
    {
      double compare_value = compare_data[Variant_Domain_idx(Variant_Grid_domain(compare), x,y,z)];
      double base_value = result_data[Basic_Domain_idx(domain, x,y,z)];
      double delta = abs( compare_value - base_value );
      if( delta >= epsilon ){
        if( verbose ){
          if( passed ) printf("Failed!\n");
          printf("(%d, %d, %d): %f - %f = %f >= %f \n",
            x, y, z,
            compare_value,
            base_value,
            delta,
            epsilon
          );
        }
        passed = false;
      }
      if( !passed && !verbose ) goto verify_loop_quick_exit;
    }
  );

  verify_loop_quick_exit:
  if( passed && verbose ){
    printf("Passed!\n");
  }

  return passed;
}

Variant_Grid* jacobi( Variant_Domain* domain, Variant_Grid* source, Variant_Grid* target, int T ){
  if( ENABLE_DEBUG ){
    assert( domain != nullptr && "Error during jacobi: domain null" );
    assert( source != nullptr && "Error during jacobi: source grid null" );
    assert( target != nullptr && "Error during jacobi: target grid null" );
    assert( source != target && "Error during jacobi: source and target point to same grid" );
    assert( Variant_Domain_equal(domain, Variant_Grid_domain(source)) && "Error during jacobi: Domain and Source's domain are different" );
    assert( Variant_Domain_equal(domain, Variant_Grid_domain(target)) && "Error during jacobi: Domain and Target's domain are different" );
    assert( Variant_Grid_data(source) != nullptr && "Error during jacobi: Source grid has Null data array" );
    assert( Variant_Grid_data(target) != nullptr && "Error during jacobi: Target grid has Null data array" );
  }

  int x, y, z, t;
  Variant_Grid* result;

  for( t = 0; t < T; t += 1 ){
    double* source_data = Variant_Grid_data(source);
    double* target_data = Variant_Grid_data(target);

    Variant_Domain_fast_loop_interior(domain, x,y,z,
      {
        int idx_000t = Variant_Grid_idx(target, x,   y,   z  );
        int idx_000  = Variant_Grid_idx(source, x,   y,   z  );
        int idx_p100 = Variant_Grid_idx(source, x+1, y,   z  );
        int idx_p010 = Variant_Grid_idx(source, x,   y+1, z  );
        int idx_p001 = Variant_Grid_idx(source, x,   y,   z+1);
        int idx_n100 = Variant_Grid_idx(source, x-1, y,   z  );
        int idx_n010 = Variant_Grid_idx(source, x,   y-1, z  );
        int idx_n001 = Variant_Grid_idx(source, x,   y,   z-1);

        target_data[idx_000t] = (1.0/7.0)*(
          source_data[idx_000] +
          source_data[idx_p100] + source_data[idx_p010] + source_data[idx_p001] +
          source_data[idx_n100] + source_data[idx_n010] + source_data[idx_n001]
        );
      }
    );
    // swap and set result grid;
    result = target;
    target = source;
    source = result;
  }

  return result;
}

int main( int argc, char** argv ){

  ProgramOptions opts = parseArguments( argv, argc );

  printf( "(nx, ny, nz): (%d, %d, %d)\nT: %d\n", opts.nx, opts.ny, opts.nz, opts.T );

  Basic_Domain* base_domain = Basic_Domain_alloc( opts.nx, opts.ny, opts.nz );
  Variant_Domain* variant_domain = Variant_Domain_alloc( opts.nx, opts.ny, opts.nz );
  Variant_Grid* source_grid = Variant_Grid_alloc( variant_domain );
  Variant_Grid* target_grid = Variant_Grid_alloc( variant_domain );
  Variant_Grid_populate( variant_domain, source_grid );

  double start = omp_get_wtime();

  // Note: result_grid is either source or target grid, do not free.
  Variant_Grid* result_grid = jacobi( variant_domain, source_grid, target_grid, opts.T );

  double end = omp_get_wtime();
  double elapsed = end - start;

  printf( "Elapsed: %fs\n", elapsed );

  if( opts.verify ){
    bool passed = verify( base_domain, result_grid, opts.T, true, opts.epsilon );
  }

  Variant_Grid_dealloc( source_grid );
  Variant_Grid_dealloc( target_grid );
  Variant_Domain_dealloc( variant_domain );
  Basic_Domain_dealloc( base_domain );

  return 0;
}
