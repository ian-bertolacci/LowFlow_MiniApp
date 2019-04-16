#include <science.hpp>
#include <macros.hpp>
#include <configure.hpp>
#include <util.hpp>

#include <assert.h>
#include <CL/opencl.h>

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
  Variant_Grid* y_ssl_dat
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

  const int num_kernels = 4;
  const int num_preambles = 1;

  int err;
  cl_device_id device_id;                          /* compute device id */
  cl_context context;                              /* compute context */
  cl_command_queue commands;                       /* compute command queue */
  cl_program program;                              /* compute program */

  char* kernel_sources[num_kernels+num_preambles]; /* kernel source */
  char* kernel_names[num_kernels];                 /* names of kernel functions */
  cl_kernel kernels[num_kernels];                  /* compute kernel */

  cl_mem input;                                    /* device memory used for the input array */
  cl_mem output;                                   /* device memory used for the output array */

  /* Connect to a compute device */
  int gpu = 1;
  err = clGetDeviceIDs(nullptr, gpu ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, 1, &device_id, nullptr);
  if( err != CL_SUCCESS ){
    printf("Error: Failed to create a device group! (Error code %d)\n", err);
  }

  /* Create a compute context */
  context = clCreateContext(0, 1, &device_id, nullptr, nullptr, &err);
  if( !context ){
    printf("Error: Failed to create a compute context! (Error code %d)\n", err);
  }

  /* Create a command commands */
  commands = clCreateCommandQueueWithProperties(context, device_id, nullptr, &err);
  if( !commands ){
    printf("Error: Failed to create a command commands! (Error code %d)\n", err);
  }


  kernel_sources[0] = (char*)
  "typedef struct struct_Basic_Domain { \
      int nx, ny, nz; \
    } Basic_Domain; \
    typedef struct struct_Basic_Grid { \
      __constant Basic_Domain* domain; \
      __global   double* data; \
    } Basic_Grid; \
    typedef Basic_Domain Variant_Domain; \
    typedef Basic_Grid Variant_Grid;";

  // Do baseline scientific kernel
  // NlFunctionEval:261 analogue
  kernel_names[0] = (char*)"function_NFE261";
  kernel_sources[0+num_preambles] = (char*) STRINGIZE(
    __kernel void function_NFE261(
      __constant Variant_Domain* domain,
      __global   Variant_Grid*   fp,
      __constant Variant_Grid*   dp,
      __constant Variant_Grid*   odp,
      __constant Variant_Grid*   osp,
      __constant Variant_Grid*   pop,
      __constant Variant_Grid*   sp,
      __constant Variant_Grid*   z_mult_dat
    ){
      int x;
      int y;
      int z;
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

  // NlFunctionEval:338 analogue
  kernel_names[1] = (char*)"function_NFE338";
  kernel_sources[1+num_preambles] = (char*) STRINGIZE(
    __kernel void function_NFE338(
      __constant Variant_Domain* domain,
      __constant Variant_Grid*   fp,
      __global   Variant_Grid*   dp,
      __global   Variant_Grid*   odp,
      __global   Variant_Grid*   opp,
      __global   Variant_Grid*   osp,
      __global   Variant_Grid*   pp,
      __global   Variant_Grid*   sp,
      __global   Variant_Grid*   ss,
      __global   Variant_Grid*   z_mult_dat
    ){
      int x;
      int y;
      int z;
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

  // NlFunctionEval:416 analogue
  kernel_names[2] = (char*)"function_NFE416";
  kernel_sources[2+num_preambles] = (char*) STRINGIZE(
    __kernel void function_NFE416(
      __constant Variant_Domain* domain,
      __global   Variant_Grid*   fp,
      __constant Variant_Grid*   et,
      __constant Variant_Grid*   sp,
      __constant Variant_Grid*   z_mult_dat
    ){
      int x;
      int y;
      int z;
      Variant_Grid_access(fp, x,y,z) -=
          Variant_Grid_access(z_mult_dat, x,y,z)
        * (   Variant_Grid_access(sp, x,y,z)
            * Variant_Grid_access(et, x,y,z)
          );
    }
  );

  // NlFunctionEval:551 analogue
  kernel_names[3] = (char*)"function_NFE551";
  kernel_sources[3+num_preambles] = (char*) STRINGIZE(
    __kernel void function_NFE551(
      __constant Variant_Domain* domain,
      __global   Variant_Grid*   fp,
      __global   Variant_Grid*   vx,
      __global   Variant_Grid*   vy,
      __global   Variant_Grid*   vz,
      __constant Variant_Grid*   dp,
      __constant Variant_Grid*   permxp,
      __constant Variant_Grid*   permyp,
      __constant Variant_Grid*   permzp,
      __constant Variant_Grid*   pp,
      __constant Variant_Grid*   rpp,
      __constant Variant_Grid*   z_mult_dat,
      __constant Variant_Grid*   x_ssl_dat,
      __constant Variant_Grid*   y_ssl_dat
    ){
      int x;
      int y;
      int z;
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

  for( int i = 0; i < num_preambles; ++i ){
    printf("==============================\nPreamble %i\n------------------------------\n%s\n", i+1, kernel_sources[i]);
  }

  for( int i = 0; i < num_kernels; i += 1 ){
    printf("==============================\n%s\n------------------------------\n%s\n", kernel_names[i], kernel_sources[num_preambles+i]);
  }

  /* Create the compute program from the source buffer */
  program = clCreateProgramWithSource(context, num_kernels+num_preambles, (const char**)kernel_sources, nullptr, &err);

  char* options = (char*)"-cl-std=CL2.0";

  err = clBuildProgram( program, 1, &device_id, options, nullptr, nullptr);

  if( !program || err != CL_SUCCESS ){
    printf("Error: Failed to create compute program! (Error code %d)\n", err);

    cl_build_status build_status;
    clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &build_status, NULL);

    switch( build_status ){
      case CL_BUILD_NONE: {
        printf("Build Status: CL_BUILD_NONE\n");
        break;
      }
      case CL_BUILD_ERROR: {
        printf("Build CL_BUILD_ERROR\n");
        break;
      }
      case CL_BUILD_SUCCESS: {
        printf("Build Status: CL_BUILD_SUCCESS\n");
        break;
      }
      case CL_BUILD_IN_PROGRESS: {
        printf("Build Status: CL_BUILD_IN_PROGRESS\n");
        break;
      }
    }


    const int infos = 2;
    unsigned int keys[infos] = { CL_PROGRAM_BUILD_OPTIONS, CL_PROGRAM_BUILD_LOG };
    char* names[infos] = { (char*)"CL_PROGRAM_BUILD_OPTIONS", (char*)"CL_PROGRAM_BUILD_LOG" };

    for(int i = 0; i < infos; ++i ){
      // Determine the size of the logCL_PROGRAM_BUILD_STATUS
      size_t log_size;
      clGetProgramBuildInfo(program, device_id, keys[i], 0, NULL, &log_size);

      // Allocate memory for the log
      char *log = (char *) malloc(log_size);

      // Get the log
      clGetProgramBuildInfo(program, device_id, keys[i], log_size, log, NULL);

      // Print the log
      printf("==============================\n%s\n------------------------------\n%s\n", names[i], log);
      free(log);
    }

    exit(-1);
  }

  for( int i = 0; i < num_kernels; ++i ){
    printf("%s\n", kernel_names[i]);
    kernels[i] = clCreateKernel(program, kernel_names[i], &err);
    if( !kernels[i] || err != CL_SUCCESS ){
      printf("Error: Failed to create compute kernel \"%s\"! (Error code %d)\n", kernel_names[i], err);
      exit(1);
    }
  }


  /* Need host->device copy */
  /* Need exec */
  /* Need device->host copy */

}
