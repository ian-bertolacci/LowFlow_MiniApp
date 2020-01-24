# Build
Uses CMake (developed using 3.13.3).
Assume that `$LOWFLOW_ROOT` is the path to this directory.
1. Create build directory somewhere (we'll call it $BUILD_ROOT)
2. `cd $BUILD_ROOT`
3. `cmake $LOWFLOW_ROOT`
4. `make`
5. `make install`
  - This installs the executable into `$LOWFLOW_ROOT/bin` (TODO: make install prefix)

## Configuration
There are multiple build types:
+ Release
  - Enables lots of optimization flags.
+ Debug
  - Disables optimization flags.
  - Sets gdb debug flag

They are specified at cmake configure time with `-DCMAKE_BUILD_TYPE=<build-type>`

It's worth noting that attempting to change the compiler for any target is (as of writing) frustratingly impossible.
Setting the compiler needs to be done during cmake configuration with `-DCMAKE_CXX_COMPILER=<compiler>` where compiler is a valid command, or is a path (not sure if it needs to be absolute) to the compiler executable.

Below are flags used to enable/disable certain families of builds
- `-D WITH_OPENCL=true`
  + Enables OpenCL variants.
- `-D WITH_KOKKOS=true`
  + Enables Kokkos variants (namely simple_whole_kokkos).
- `-D WITH_CUDA=true`
  + Enables CUDA variants.

## Prerequisite libraries
### OpenCL
OpenCL is required for variants using it.
All OpenCL variants are disabled by default.
They can be Enabled during configuration with `-D WITH_OPENCL=true`.

Note that there are different versions of the API, particularly 1.x vs 2.x.
The preferred API is 2.0, however because NVIDIA refuses to upgrade their implementation
to fully support this OpenCL 2.0 standard (as ratified in late 2013), OpenCL 1.2
must be also supported.

### Kokkos
[Kokkos](https://github.com/kokkos/kokkos) is required for variants using it.
All Kokkos variants are disabled by default.
They can be enabled during configuration with `-D WITH_KOKKOS=true`.
Also enabling CUDA (with `-D WITH_CUDA=true`) will also enable the Kokkos CUDA variant.

# Running
Each program is a different, self-contained variant, and will be run in roughly the same way.

Command Line Options:
+ `-V` : perform verification after experiment to check that experimental output is correct.
  Program returns `255` if verification fails.
+ `-x <N : int>` : Set size of domain in x direction.
+ `-y <N : int>` : Set size of domain in y direction.
+ `-z <N : int>` : Set size of domain in z direction.
+ `-s <N : int>` : Set size of domain in x, y, z directions to the same value.
  `-s 10` is equivalent to `-x 10 -y 10 -z 10`.
+ `-e <epsilon : float>` : Set acceptable error bound used for verification.
  An absolute difference between experimental and verification output grids less than epsilon is considered 'passing'
+ `-S <N : uint>` : Set the seed used to generate grid-seeds.
+ `-h` : Prints this message (including variant's message) and exit codes.
+ `--` : End of common command line arguments, beginning of variant specific command line arguments.

Variants can also have their own command line options, which are also printed in the help message.

## Exit Codes
+ SUCCESS (0) : Successful run.
+ VERIFICATION_FAIL (255) : Verification failed (only if `-V` flag enabled).
+ COMMAND_LINE_ERROR (254) : Invalid command line options.
+ ASSERTION_FAILURE (134) : Assertion caused program abort.

Variant can also have their own exit codes, which are also printed in the help message.


# Creating a new Variant

In an (possibly misguided) attempt to simplify the development of each variant,
all variants share a common source (files directrly under `src/`), and implement
common interfaces with their variant specific source.
+ Type system
+ Macro system
+ Configuration system
+ Metric system
+ Setup/Teardown system
+ Science system

## Type System
Each variant must define the following types:
+ `Variant_Domain`
  - This represents a 3D rectilinear index space starting at (0,0,0) and ending
    at (nx,ny,nz), as provided by the constructor.
  - This can be implemented in any reasonable manner. In many cases, it simply
    mirrors the Basic_Domain struct.
  - (At the time of writing) a Variant_Domain instance is not manipulated or
    queried, and is usually only used in the Variant_Domain_loop macros.
+ `Variant_Grid`
  - This represents a 3D rectilinear array whose index space is defined by a domain
    provided to the constructor.
  - It may be represented in any way.
  - When facing 'global' code (particularly code in `main.cpp`) all indices in
    the domain must be accessible at any time and must have their correct value.

Each variant must also define the following functions:
+ `Variant_Domain* Variant_Domain_alloc( int nx, int ny, int nz )`
  - Allocate a Variant_Domain instance of size nx, ny, nx
+ `void Variant_Domain_dealloc( Variant_Domain* domain )`
  - Free Variant_Domain instance.
+ `Variant_Grid* Variant_Grid_alloc( Variant_Domain* domain )`
  - Allocate Variant_Grid instance given a Variant_Domain instance.
+ `void Variant_Grid_dealloc( Variant_Grid* grid )`
  - Free Variant_Grid instance.
  + `void Variant_Grid_populate_zero( Variant_Domain* domain, Variant_Grid* grid )`
  - Populate a Variant_Grid with 0.0 at all locations.
+ `void Variant_Grid_populate_seeded( Variant_Domain* domain, Variant_Grid* grid, unsigned int seed )`
  - Populate a grid (in lexicographic order) with values from standard `random` random number generator seeded with value `seed`.
+ `void Variant_Grid_populate_int_increment( Variant_Domain* domain, Variant_Grid* grid )`
  - Populate a Variant_Grid with values starting at (0,0,0) with 1.0, such each index is one greater than the previous index (in lexicographic order).

## Macro System
Each variant must define the following macros:
+ `Variant_Grid_access(grid,x,y,z)`
  - 'Returns' the element at position (x,y,z) in the grid instance.
+ `Variant_Domain_loop_whole(domain, iter_x, iter_y, iter_z, body)`
  - Iterate over *all* indices in the domain, and invoking body at every iteration.
  - The iterator values are stored into the variables passed into `iter_x`, `iter_y`, `iter_z`.
    These variables must be declared before invoking the macro.
+ `Variant_Domain_loop_interior(domain, iter_x, iter_y, iter_z, body)`
  - Iterate over all *interior* indices in the domain (indices which are not at the faces of the domain), and invoking body at every iteration.
  - The iterator values are stored into the variables passed into `iter_x`, `iter_y`, `iter_z`.
    These variables must be declared before invoking the macro

### Global Macro System
Similarly, the global macro system has the macros:
+ `Basic_Grid_access(grid,x,y,z)`
  - 'Returns' the element at position (x,y,z) in the grid instance.
+ `Basic_Domain_loop_whole(domain, iter_x, iter_y, iter_z, body)`
  - Iterate over *all* indices in the domain, and invoking body at every iteration.
  - The iterator values are stored into the variables passed into `iter_x`, `iter_y`, `iter_z`.
    These variables must be declared before invoking the macro.
+ `Basic_Domain_loop_interior(domain, iter_x, iter_y, iter_z, body)`
  - Iterate over all *interior* indices in the domain (indices which are not at the faces of the domain), and invoking body at every iteration.
  - The iterator values are stored into the variables passed into `iter_x`, `iter_y`, `iter_z`.
    These variables must be declared before invoking the macro

As well as the below macros (which, while mandatory, are often also implemented on the variant side):
+ `Basic_Domain_nx(domain)`
  - 'Returns' domain's size in the X direction.
+ `Basic_Domain_ny(domain)`
  - 'Returns' evaluates to a domain's size in the Y direction.
+ `Basic_Domain_nz(domain)`
  - 'Returns' evaluates to a domain's size in the Z direction.
+ `Basic_Domain_idx(domain, x,y,z)`
  - 'Returns' a flat index from the 3D tuple (x,y,z).
+ `Basic_Grid_data(grid)`
  - 'Returns' a Basic_Grid's backing array.
+ `Basic_Grid_domain(grid)`
  - 'Returns' a Basic_Grid's backing domain.
+ `Basic_Grid_idx(grid, x,y,z)`
  - Same as Basic_Domain_idx using the Grid's backing domain.
+ `Basic_Grid_access_index(grid, idx)`
  - Similar to Basic_Grid_access, but using a flat index (such as one produced by Basic_Domain_idx).
+ `Basic_Domain_equal(domain_a, domain_b)`
  - Evalutates to true of the two domains are equal.
+ `Basic_Domain_fast_loop_whole(domain, iter_x, iter_y, iter_z, body)`
  - (As time of writing) Mostly unused, is the same as Basic_Domain_loop_whole.
+ `Basic_Domain_fast_loop_interior(domain, iter_x, iter_y, iter_z, body)`
  - (As time of writing) Mostly unused, is the same as Basic_Domain_loop_interior.

(Yes, many of these could also be functions, but we'd like the expressions they
represent to be embedded directly, so that the compiler can very easily perform
expression optimizations such as common-subexpression-elimination.)

## Configuration System

There are several components to defining command line arguments for a variant.
While it is not necessary to provide additional command line arguments, the basic skeleton does need to be implemented.

Several types must be defined:
+ `VariantOptions`
  - This stores the program options resulting from parsing the command line arguments.
  - Usually a struct, typedef-ed to simply be the typename `VariantOptions`.
    * Can be empty.
+ (optional) `enum VariantExitCode : ExitType`
  - Defines exit codes that are specific to this variant (this is rare).
  - Actual values should be explicitly defined.
    * See `GeneralExitCode` to avoid overlap with previously defined exit code values.
  - `ExitType` is a `uint8_t`.

Several global static variables need to be defined.
+ `static const VariantOptions VariantOptions_Default`
  - Stores the default values for all options of the variant.
+ static const struct option variant_long_options[]
  - This is an array of `stuct option` objects.
    * the option struct has the following compoenents:
      + `const char* name` : the full name of the option; will become the string of the command line flag.
      + `int has_arg` : indicates if the flag takes an argument.
        - Legal values for `has_arg` (these are defined by the getopt library):
          * `no_argument` : flag takes no arguments.
          * `required_argument` : flag takes one, mandatory argument.
          * `optional_argument` : flag takes one, optional argument.
      + `int* flag` : (Just make use `nullptr`)
      + `int val` : (Just use a char literal indicating the 'short' version of the flag (e.g. `'h'` for `--help`))
    * `struct option` comes from getopt. Please see the [GNU getopt-long documentation](https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Options.html) for more details.
  - It must be terminated with the `LONG_OPTIONS_TERMINAL` macro.

This is all done in place with a [compound literal](https://en.cppreference.com/w/c/language/compound_literal).

Several functions must also be defined:
+ `VariantOptions parseVariantOptions( int argc, char** argv )`
  - Parses (ideally using `getopt_long`) the variant's command line options (and *only* the variant's) and returns an `VariantOptions` instance.
  - If no arguments are present, the returned `VariantOptions` should be equivalent to `VariantOptions_Default`
+ `void printVariantInformationMessage( FILE* stream )`
  - Prints a quick, 1-2 sentence description of the variant.
+ `void printVariantOptionsMessage( FILE* stream )`
  - Prints a help message listing the variant's command line options (if any) and their purpose.
+ `void printVariantExitCodeMessage( FILE* stream )`
  - Prints a message listing exit codes introduced by this variant (if any) and their integer values.

If the variant has no flags actual implementation for the message function can simply be to do absolutely nothing.

Below is a template/example for the configuration code.

The below should go in the header (variant_configure.hpp), and defines types and defaults.
```cpp
#ifndef VARIANT_CONFIGURE_HPP
#define VARIANT_CONFIGURE_HPP

#include <global_configure.hpp>

typedef struct struct_VariantOptions {
  bool verbose;
} VariantOptions;

static const VariantOptions VariantOptions_Default {
  .verbose = false
};

static const struct option variant_long_options[] = {
  {"help",           no_argument, nullptr, 'h'},
  {"verbose",        no_argument, nullptr, 'v'},
  LONG_OPTIONS_TERMINAL
};

enum VariantExitCode : ExitType {
  TURBO_ENCABULATOR_DESYNCRONIZATION_EVENT = 123
};

#endif
```

The below should go in a cpp files (such as variant_configure.cpp), and defines the functions.
```cpp
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
  fprintf( stream, "This is the manual OpenCL with temprary storage implementation of the LowFlow mini-app\n" );
}

void printVariantOptionsMessage( FILE* stream ){
  fprintf( stream,
    "Variant Specific Command Line Options:\n"
    "  --help\n"
    "   -h\n"
    "    Print this message.\n"
    "  --verbose\n"
    "   -v\n"
    "    Enable verbose printing during experiment run.\n"
  );
}

void printVariantExitCodeMessage( FILE* stream ){
  fprintf( stream,
    "Variant Exit Codes:\n"
    "-  TURBO_ENCABULATOR_DESYNCRONIZATION_EVENT:  %u\n"
    "  + Turbo-Encabulator has encountered an issue where it can no longer automatically synchronize cardinal grammeters"
    TURBO_ENCABULATOR_DESYNCRONIZATION_EVENT
  );
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
    optind = 1; // Reset getopt_long
    opterr = 0;
    while( (c = getopt_long(args.argc, args.argv, short_options, variant_long_options, &option_index)) != -1 ){
      switch (c){
        case 'h':
          printHelpMessage( stdout );
          exit(SUCCESS);

        case 'v':
          opts.verbose =  true;
          break;

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
```

## Metrics system
Variants each have their own metrics that they collect during the program's run (for example, the runtime of a specific loop).

It is required to define the type: `Variant_Metrics`.
This is usually a struct, typedef-ed to simply be the typename `Variant_Metrics`.
All the fields of this type should be a different metric being measured.

It is also required to define the function: `void printVariantMetricInformation( FILE* stream, Variant_Metrics* metrics )`
This should write the name and values of each field in the `Variant_Metrics` instance to the file stream.
This is invoked automatically at the end of the program.

Several timing macros are provided to simplify timing of code sections (defined in metrics.hpp):
+ `START_TIMER( timer_variable )`: Start the timer `timer_variable`.
+ `STOP_TIMER( timer_variable )`: Stop the timer `timer_variable`; `timer_variable` contains the elapsed time in seconds.
+ `TIMEIT( timer_variable, body )`: Time the execution of the body; `timer_variable` contains the elapsed time in seconds.
+ `TIMEIT_ACCUMULATE( timer_variable, body )`: Time the execution of the body; `timer_variable` contains the elapsed time in seconds plus its original value.
`timer_variable` must be a `double` type lvalue.

When `ENABLE_METRICS` is *not* defined, only the body is executed (if the macro is provided one); the timing code is not, and so the timer_variables do not contain the elapsed time.
When `ENABLE_METRICS` *is* defined, *both* the body and the timing code is executed.

(NOTE: `TIMEIT` and `TIMEIT_ACCUMULATE` are implemented in a way that allows un-parenthesized commas, such as with multiple variable declaration and cuda kernel invocations).

Below is a template/example for the metrics code.

The below should go in the header (variant_metrics.hpp), and defines types and defaults.
```cpp
#ifndef VARIANT_METRICS_HPP
#define VARIANT_METRICS_HPP

#include <configure.hpp>
#include <omp.h>

typedef struct struct_Variant_Metrics {
  double elapsed_216;
  double elapsed_338;
  double elapsed_416;
  double elapsed_551;
  double elapsed_551_reduce;
  double elapsed_setup;
  double elapsed_compile;
  double elapsed_copy_host_to_device;
  double elapsed_copy_device_to_host;
  double elapsed_exec_setup;
  double elapsed_exec;
  double elapsed_teardown;
} Variant_Metrics;

#endif
```

The below should go in a cpp file (such as variant_metrics.cpp), and defines the functions.
```cpp
#include <metrics.hpp>

void printVariantMetricInformation( FILE* stream, Variant_Metrics* metrics ){
    fprintf( stream,
      "Elapsed 261: %f\n"
      "Elapsed 338: %f\n"
      "Elapsed 416: %f\n"
      "Elapsed 551: %f\n"
      "Elapsed 551 reduce: %f\n"
      "Elapsed setup: %f\n"
      "Elapsed compile: %f\n"
      "Elapsed host->device: %f\n"
      "Elapsed device->host: %f\n"
      "Elapsed setup execution: %f\n"
      "Elapsed execution: %f\n"
      "Elapsed teardown: %f\n",

      metrics->elapsed_216,
      metrics->elapsed_338,
      metrics->elapsed_416,
      metrics->elapsed_551,
      metrics->elapsed_551_reduce,
      metrics->elapsed_setup,
      metrics->elapsed_compile,
      metrics->elapsed_copy_host_to_device,
      metrics->elapsed_copy_device_to_host,
      metrics->elapsed_exec_setup,
      metrics->elapsed_exec,
      metrics->elapsed_teardown
    );
}
```

## Setup/Teardown system
Some variants may require different setup and tear-down processes that must be called at the beginning and end of the program.

This is done by defining the (mandatory) functions:
+ `void programSetup( ProgramOptions program_options )`
  - Invoked at the beginning of the program after parsing command line arguments.
+ `void programTeardown( ProgramOptions program_options )`
  - Invoked at the end of the program.

They are provided the program options from parsing the command line arguments via `ProgramOptions`

`ProgramOptions` is a struct with the following fields:
+ `ArgsParameters all_argument_parameters`
  - Contains the `main`'s `argc` and `argv` values
    * Access via fields `int argc` and `char** argv`
+ `int nx`
  - Size of all grids in the X direction.
+ `int ny`
  - Size of all grids in the Y direction.
+ `int nz`
  - Size of all grids in the Z direction.
+ `double epsilon`
  - The maximum allowed difference between experiemental and baseline output grids
+ `unsigned int seed`
  - Seed value for creating random grid values
+ `bool verify`
  - True if `--verify` flag used.
+ `VariantOptions variant_options`
  - The variant options struct, containing its own fields (defined by you).

## Science system

The core of the variant implementation is the `science` function.
This is adapted from the [Parflow project](https://github.com/parflow/parflow), specifically, [NlFunctionEval](https://github.com/parflow/parflow/blob/master/pfsimulator/parflow_lib/nl_function_eval.c#L109).

(Comment: I, Ian Bertolacci, have no idea what science or mathematics are happening, so I will not describe it).

The following grids are inputs:
+ dp
+ et
+ odp
+ opp
+ osp
+ permxp
+ permyp
+ permzp
+ pop
+ pp
+ rpp
+ sp
+ ss
+ z_mult_dat
+ x_ssl_dat
+ y_ssl_dat

The following grids are outputs:
+ fp
+ vx
+ vy
+ vz

Below is a template/example implementation.
```cpp
void science(
  Variant_Domain* domain,
  Variant_Grid* fp,
  Variant_Grid* vx,
  Variant_Grid* vy,
  Variant_Grid* vz,
  Variant_Grid* dp,
  Variant_Grid* et,
  Variant_Grid* odp,
  Variant_Grid* opp,
  Variant_Grid* osp,
  Variant_Grid* permxp,
  Variant_Grid* permyp,
  Variant_Grid* permzp,
  Variant_Grid* pop,
  Variant_Grid* pp,
  Variant_Grid* rpp,
  Variant_Grid* sp,
  Variant_Grid* ss,
  Variant_Grid* z_mult_dat,
  Variant_Grid* x_ssl_dat,
  Variant_Grid* y_ssl_dat,
  VariantOptions options,
  Variant_Metrics* metrics
){
  if( ENABLE_DEBUG ){
    assert( domain != nullptr     && "Error during science: domain is null" );
    assert( fp != nullptr         && "Error during science: output grid fp is null" );
    assert( vx != nullptr         && "Error during science: output grid vx is null" );
    assert( vy != nullptr         && "Error during science: output grid vy is null" );
    assert( vz != nullptr         && "Error during science: output grid vz is null" );
    assert( dp != nullptr         && "Error during science: input grid dp is null" );
    assert( et != nullptr         && "Error during science: input grid et is null" );
    assert( odp != nullptr        && "Error during science: input grid odp is null" );
    assert( opp != nullptr        && "Error during science: input grid opp is null" );
    assert( osp != nullptr        && "Error during science: input grid osp is null" );
    assert( permxp != nullptr     && "Error during science: input grid permxp is null" );
    assert( permyp != nullptr     && "Error during science: input grid permyp is null" );
    assert( permzp != nullptr     && "Error during science: input grid permzp is null" );
    assert( pop != nullptr        && "Error during science: input grid pop is null" );
    assert( pp != nullptr         && "Error during science: input grid pp is null" );
    assert( rpp != nullptr        && "Error during science: input grid rpp is null" );
    assert( sp != nullptr         && "Error during science: input grid sp is null" );
    assert( ss != nullptr         && "Error during science: input grid ss is null" );
    assert( z_mult_dat != nullptr && "Error during science: input grid z_mult_dat is null" );
    assert( x_ssl_dat != nullptr  && "Error during science: input grid x_ssl_dat is null" );
    assert( y_ssl_dat != nullptr  && "Error during science: input grid y_ssl_dat is null" );
  }

  // Do baseline scientific kernel
  // NlFunctionEval:261 analogue
  TIMEIT( metrics->elapsed_216,
    {
      int x;
      int y;
      int z;
      Variant_Domain_fast_loop_interior(domain, x,y,z,
        {
          Variant_Grid_access(fp, x,y,z) =
           (  Variant_Grid_access(sp,  x,y,z)
            * Variant_Grid_access(dp,  x,y,z)
            - Variant_Grid_access(osp, x,y,z)
            * Variant_Grid_access(odp, x,y,z)
           )
           * Variant_Grid_access(pop, x,y,z)
           * Variant_Grid_access(z_mult_dat, x,y,z);
        }
      );
    }
  )

  // NlFunctionEval:338 analogue
  TIMEIT( metrics->elapsed_338,
    {
      int x;
      int y;
      int z;
      Variant_Domain_fast_loop_interior(domain, x,y,z,
        {
          Variant_Grid_access(fp, x,y,z) +=
              Variant_Grid_access(ss, x,y,z)
            * Variant_Grid_access(z_mult_dat, x,y,z)
            * (   Variant_Grid_access(pp, x,y,z)
                * Variant_Grid_access(sp, x,y,z)
                * Variant_Grid_access(dp, x,y,z)
                - Variant_Grid_access(opp, x,y,z)
                * Variant_Grid_access(osp, x,y,z)
                * Variant_Grid_access(odp, x,y,z)
              );
        }
      );
    }
  )
  // NlFunctionEval:416 analogue
  TIMEIT( metrics->elapsed_416,
    {
      int x;
      int y;
      int z;
      Variant_Domain_fast_loop_interior(domain, x,y,z,
        {
          Variant_Grid_access(fp, x,y,z) -=
              Variant_Grid_access(z_mult_dat, x,y,z)
            * (   Variant_Grid_access(sp, x,y,z)
                * Variant_Grid_access(et, x,y,z)
              );
        }
      );
    }
  )

  // NlFunctionEval:551 analogue
  TIMEIT( metrics->elapsed_551,
    {
      int x;
      int y;
      int z;
      Variant_Domain_fast_loop_interior(domain, x,y,z,
        {

          double x_dir_g   = ArithmeticMean( Variant_Grid_access( x_ssl_dat, x, y, 0), Variant_Grid_access( x_ssl_dat, x+1,  y, 0 ) );
          double x_dir_g_c = ArithmeticMean( Variant_Grid_access( x_ssl_dat, x, y, 0), Variant_Grid_access( x_ssl_dat, x+1,  y, 0 ) );
          double y_dir_g   = ArithmeticMean( Variant_Grid_access( y_ssl_dat, x, y, 0), Variant_Grid_access( y_ssl_dat, x,  y+1, 0 ) );
          double y_dir_g_c = ArithmeticMean( Variant_Grid_access( y_ssl_dat, x, y, 0), Variant_Grid_access( y_ssl_dat, x,  y+1, 0 ) );

          double diff_right = Variant_Grid_access(pp, x,y,z) - Variant_Grid_access(pp, x+1, y,   z);
          double diff_front = Variant_Grid_access(pp, x,y,z) - Variant_Grid_access(pp, x,   y+1, z);

          double updir_right = diff_right * x_dir_g_c - x_dir_g;
          double updir_front = diff_front * y_dir_g_c - y_dir_g;

          double sep = ArithmeticMean( Variant_Grid_access(z_mult_dat, x,y,z),  Variant_Grid_access(z_mult_dat, x, y, z+1) );

          double lower_cond =
            Variant_Grid_access(pp, x,y,z) / sep
          - (   Variant_Grid_access(z_mult_dat, x,y,z)
              / ( Variant_Grid_access(z_mult_dat, x,y,z)
                + Variant_Grid_access(z_mult_dat, x, y, z+1)
                )
            )
          * Variant_Grid_access(dp, x,y,z);

          double upper_cond =
            Variant_Grid_access(pp, x, y, z+1) / sep
          + (   Variant_Grid_access(z_mult_dat, x, y, z+1)
              / ( Variant_Grid_access(z_mult_dat, x,y,z)
                + Variant_Grid_access(z_mult_dat, x, y, z+1)
                )
            )
          * Variant_Grid_access(dp, x,y,z+1);

          double diff_upper = lower_cond - upper_cond;


          // NOTE!
          // Originally the harmonic mean was
          // PMean( Variant_Grid_access(pp, x,y,z),
          //        Variant_Grid_access(pp, x+1, y, z),
          //        permxp[ip],
          //        permxp[ip + 1]
          // )
          // However! PMean(a,b,c,d) in parflow is simply HarmonicMean(c, d)
          // so the first two terms have been removed

          double u_right =
            (
               Variant_Grid_access(z_mult_dat, x,y,z)
             * HarmonicMean( Variant_Grid_access(permxp, x,   y, z),
                             Variant_Grid_access(permxp, x+1, y, z) )
             * diff_right * x_dir_g_c
             * UpstreamMean(
                 updir_right,
                 0.0,
                 Variant_Grid_access(rpp, x,   y, z) * Variant_Grid_access(dp, x,   y, z),
                 Variant_Grid_access(rpp, x+1, y, z) * Variant_Grid_access(dp, x+1, y, z)
               )
           ) + (
              Variant_Grid_access(z_mult_dat, x,y,z)
             * HarmonicMean( Variant_Grid_access(permxp, x,   y, z),
                             Variant_Grid_access(permxp, x+1, y, z) )
             * (-x_dir_g)
             * UpstreamMean(
                 updir_right,
                 0.0,
                 Variant_Grid_access(rpp, x,   y, z) * Variant_Grid_access(dp, x,   y, z),
                 Variant_Grid_access(rpp, x+1, y, z) * Variant_Grid_access(dp, x+1, y, z)
               )
           );

           double u_front =
             (
                Variant_Grid_access(z_mult_dat, x,y,z)
              * HarmonicMean( Variant_Grid_access(permyp, x,   y, z),
                              Variant_Grid_access(permyp, x+1, y, z) )
              * diff_front * x_dir_g_c
              * UpstreamMean(
                  updir_front,
                  0.0,
                  Variant_Grid_access(rpp, x,   y, z) * Variant_Grid_access(dp, x,   y, z),
                  Variant_Grid_access(rpp, x+1, y, z) * Variant_Grid_access(dp, x+1, y, z)
                )
            ) + (
               Variant_Grid_access(z_mult_dat, x,y,z)
              * HarmonicMean( Variant_Grid_access(permyp, x,   y, z),
                              Variant_Grid_access(permyp, x+1, y, z) )
              * (-x_dir_g)
              * UpstreamMean(
                  updir_front,
                  0.0,
                  Variant_Grid_access(rpp, x,   y, z) * Variant_Grid_access(dp, x,   y, z),
                  Variant_Grid_access(rpp, x+1, y, z) * Variant_Grid_access(dp, x+1, y, z)
                )
            );

          double u_upper =
                      HarmonicMeanDZ(
                          Variant_Grid_access(permzp,     x, y, z  ),
                          Variant_Grid_access(permzp,     x, y, z+1),
                          Variant_Grid_access(z_mult_dat, x, y, z  ),
                          Variant_Grid_access(z_mult_dat, x, y, z+1)
                      )
                    * diff_upper
                    * UpstreamMean(
                        lower_cond,
                        upper_cond,
                        Variant_Grid_access(rpp, x, y, z  ) * Variant_Grid_access(dp, x, y, z  ),
                        Variant_Grid_access(rpp, x, y, z+1) * Variant_Grid_access(dp, x, y, z+1)
                      );

          Variant_Grid_access(vx, x,y,z) = u_right;
          Variant_Grid_access(vy, x,y,z) = u_front;
          Variant_Grid_access(vz, x,y,z) = u_upper;

          Variant_Grid_access(fp, x  , y  , z  ) += u_right * u_front * u_upper;
          Variant_Grid_access(fp, x+1, y  , z  ) += u_right;
          Variant_Grid_access(fp, x  , y+1, z  ) -= u_front;
          Variant_Grid_access(fp, x  , y  , z+1) -= u_upper;
        }
      );
    }
  )
}
```

## Files
All variants share the code directly under `src/`.
Each variant has it's own directory under `src/` with the following files:
1. CMakeLists.txt
  + Defines compilation targets for this variant.
    This is further explained below.
2. Header files
  + Required files are
    - variant_configure.hpp
      * Should define varaint configuration system as described above.
    - variant_types.hpp
      * Should define variant type system as described above.
    - variant_macros.hpp
      * Should define the variant macro system as described above.
    - variant_metrics.hpp
      Should define the variant metrics system as described above.
  + These are included statically in the various global*.hpp headers
3. Source C++ files
  + Their names can be whatever, as they are specified to the build system in
    the variant's CMakeLists.txt file.

## Build System
There are two parts to adding to the CMake build process
First, a new CMakeLists.txt is created in the variants source directory.
Here is an example/template for a variant's CMakeLists.txt file:

```cmake
# Mandatory
cmake_minimum_required(VERSION 3.5)
include(${CMAKE_DIR}/utilities.txt)
include(${CMAKE_DIR}/build_variant.txt)

# Not necessary but useful
# This is this varaint's base target name (executable filename)
# Different complation configuraitons will have enabled macros appended to it.
set(VARIANT_NAME template_variant_name)
# This is the list of source files for this variant, relative to the variant's source directory
set(VARIANT_SOURCE_FILES variant_types.cpp variant_science.cpp variant_metrics.cpp variant_configure.cpp variant_setup_teardown.cpp)

# This creates a single target
CREATE_VARIANT_TARGET(
  # Name of the target
  NAME ${VARIANT}
  # The list of global source files (this has been defined by src/CMakeLists.txt)
  BASE_SRC_FILES ${BASE_SOURCE_FILES}
  # The list of the variant's source files
  VARIANT_SRC_FILES ${VARIANT_SOURCE_FILES}
  # Where the variant should be installed (Please do not change this.)
  INSTALL_DIR ${INSTALL_BIN_DIR}
)

# This creates another, single target.
# In this case, we are creating a target with the compile time configuration
# macro `ENABLE_METRICS` enabled.
CREATE_VARIANT_TARGET(
  # Base name of target. NOTE! It's the same as the above, but will automatically
  # have the configuration macro appended to it, creating a new name.
  NAME ${VARIANT}
  # The list of global source files (this has been defined by src/CMakeLists.txt)
  BASE_SRC_FILES ${BASE_SOURCE_FILES}
  # The list of the variant's source files
  VARIANT_SRC_FILES ${VARIANT_SOURCE_FILES}
  # Where the variant should be installed
  INSTALL_DIR ${INSTALL_BIN_DIR}
  # List of one or more macro tokens that will be enabled in compilation (using -DTOKEN_NAME)
  # for this target
  CONFIG_MACRO_DEFNS ENABLE_METRICS
)
```

CREATE_VARIANT_TARGET function can change various aspects of a target's compilation using
the following arguments (see cmake_parse_arguments for deeper explanation of how this works):

+ Boolean arguments (flags that take no additional argument parameters)
  - VARIANT_CXX_OVERWRITE_FLAGS
    * If enabled, overwrites the compilation flags with from VARIANT_CXX_FLAGS argument (see below)
  - DEBUG_CMAKE_CONFIGURE
    * if enabled, prints additional verbose debugging information during cmake configuration stage (not compilation)
+ Single arguments (take exactly one argument parameter)
  - NAME <name : string>
    * Sets base name of target. If configuration macro tokens are provided (using CONFIG_MACRO_DEFNS, see below), those tokens are appended to this name for this compilation target.
  - INSTALL_DIR <path : string>
    * Sets the installation location for this target.
+ List arguments (take one or more argument parameters)
  - BASE_SRC_FILES <paths : list<string>>
    * Takes list of paths to source files common to all variants (as defined globally in BASE_SOURCE_FILES variable).
  - VARIANT_SRC_FILES <paths : list<string>>
    * Takes list of paths to source files for this variant (ideally the same for all targets in this variant family).
  - VARIANT_INCLUDE_DIRS <paths : list<string>>
    * If necessary, takes list of paths to include directories necessary to build this target.
  - VARIANT_LINK_LIBS  <path : list<string>>
    * If necessary, takes list of library names necessary to build this target.
  - VARIANT_LINK_DIRS <path : list<string>>
    * If necessary, takes list of paths to directories for .
  - CONFIG_MACRO_DEFNS  <token : list<string>>
    * List of macros to be 'defined' (using -DTOKEN_NAME) to the compiler,
      to enabled code sections using the `#ifdef TOKEN ... #endif` style.
  - VARIANT_CXX_FLAGS <compiler-flag : list<string>>
    * List of compiler flags to be used for this specific target.

To set an argument, first specify the argument name (such as `CONFIG_MACRO_DEFNS`)
followed by the argument parameters, if any.
For arguments with single parameters, you shouldn't have to do anything.
For arguments with a list of 2 or more elements, please use a list instance (see https://cmake.org/cmake/help/v3.5/command/list.html for more details)

It's worth noting that attempting to change the compiler for any target is (as of writing) frustratingly impossible.
Setting the compiler needs to be done during cmake configuration with `-DCMAKE_CXX_COMPILER=<compiler>` where compiler is a valid command, or is a path (not sure if it needs to be absolute) to the compiler executable.
