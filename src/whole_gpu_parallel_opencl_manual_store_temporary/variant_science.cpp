#include <science.hpp>
#include <macros.hpp>
#include <configure.hpp>
#include <util.hpp>

#include <assert.h>
#include <CL/opencl.h>
#include <math.h>
#include <stdio.h>

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

  const int num_domains = 24;

  const int num_compute_kernels = 4 /*main kernels*/ + 1 /*reduce kernels*/;
  const int num_util_kernels = 0;
  const int num_preambles = 1 + num_domains;

  // Error code returned/out-var'ed by OpenCL functions
  int err;
  // compute device id
  cl_device_id device_id;
  // compute context
  cl_context context;
  // compute command queue
  cl_command_queue commands;
  // compute program
  cl_program program;

  // all source
  char* sources[num_preambles+num_compute_kernels+num_util_kernels];
  // preamble sources
  char** preamble_sources = &sources[0];
  // only compute kernel sources
  char** compute_kernel_sources = &sources[num_preambles];
  // only utility kernel sources
  char** util_kernel_sources = &sources[num_preambles+num_compute_kernels];
  // names of kernel functions
  char* kernel_names[num_compute_kernels];
  // compute kernels
  cl_kernel compute_kernels[num_compute_kernels];
  // utility function kernels
  cl_kernel util_kernels[num_compute_kernels];
  // array of the number of args for each kernel
  int kernel_num_args[num_compute_kernels];
  // array of arrays of args in order for each kernel : kernel_idx - > ( arg_idx -> cl_mem* )
  cl_mem** kernel_args[num_compute_kernels];

  // work_group_size domain size available to each kernel
  size_t work_group_size[num_compute_kernels];
  // Total work (domain size) done by each kernel
  size_t total_work[num_compute_kernels][3];
  // Number of work units per group for each kernel;
  size_t local[num_compute_kernels][3];

  /* Connect to a compute device */
  bool use_gpu = true;
  err = clGetDeviceIDs(nullptr, use_gpu ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, 1, &device_id, nullptr);
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

  Variant_Grid* u_right = Variant_Grid_alloc( domain );
  Variant_Grid* u_front = Variant_Grid_alloc( domain );
  Variant_Grid* u_upper = Variant_Grid_alloc( domain );

  Variant_Grid_populate_zero( Variant_Grid_domain(u_right), u_right );
  Variant_Grid_populate_zero( Variant_Grid_domain(u_right), u_front );
  Variant_Grid_populate_zero( Variant_Grid_domain(u_right), u_upper );

  // Device memories for each of the function args
  cl_mem device_fp         = clCreateBuffer( context, CL_MEM_READ_WRITE, Variant_Grid_sizeof(fp),         nullptr, &err );
  cl_mem device_vx         = clCreateBuffer( context, CL_MEM_WRITE_ONLY, Variant_Grid_sizeof(vx),         nullptr, &err );
  cl_mem device_vy         = clCreateBuffer( context, CL_MEM_WRITE_ONLY, Variant_Grid_sizeof(vy),         nullptr, &err );
  cl_mem device_vz         = clCreateBuffer( context, CL_MEM_WRITE_ONLY, Variant_Grid_sizeof(vz),         nullptr, &err );
  cl_mem device_dp         = clCreateBuffer( context, CL_MEM_READ_ONLY,  Variant_Grid_sizeof(dp),         nullptr, &err );
  cl_mem device_et         = clCreateBuffer( context, CL_MEM_READ_ONLY,  Variant_Grid_sizeof(et),         nullptr, &err );
  cl_mem device_odp        = clCreateBuffer( context, CL_MEM_READ_ONLY,  Variant_Grid_sizeof(odp),        nullptr, &err );
  cl_mem device_opp        = clCreateBuffer( context, CL_MEM_READ_ONLY,  Variant_Grid_sizeof(opp),        nullptr, &err );
  cl_mem device_osp        = clCreateBuffer( context, CL_MEM_READ_ONLY,  Variant_Grid_sizeof(osp),        nullptr, &err );
  cl_mem device_permxp     = clCreateBuffer( context, CL_MEM_READ_ONLY,  Variant_Grid_sizeof(permxp),     nullptr, &err );
  cl_mem device_permyp     = clCreateBuffer( context, CL_MEM_READ_ONLY,  Variant_Grid_sizeof(permyp),     nullptr, &err );
  cl_mem device_permzp     = clCreateBuffer( context, CL_MEM_READ_ONLY,  Variant_Grid_sizeof(permzp),     nullptr, &err );
  cl_mem device_pop        = clCreateBuffer( context, CL_MEM_READ_ONLY,  Variant_Grid_sizeof(pop),        nullptr, &err );
  cl_mem device_pp         = clCreateBuffer( context, CL_MEM_READ_ONLY,  Variant_Grid_sizeof(pp),         nullptr, &err );
  cl_mem device_rpp        = clCreateBuffer( context, CL_MEM_READ_ONLY,  Variant_Grid_sizeof(rpp),        nullptr, &err );
  cl_mem device_sp         = clCreateBuffer( context, CL_MEM_READ_ONLY,  Variant_Grid_sizeof(sp),         nullptr, &err );
  cl_mem device_ss         = clCreateBuffer( context, CL_MEM_READ_ONLY,  Variant_Grid_sizeof(ss),         nullptr, &err );
  cl_mem device_z_mult_dat = clCreateBuffer( context, CL_MEM_READ_ONLY,  Variant_Grid_sizeof(z_mult_dat), nullptr, &err );
  cl_mem device_x_ssl_dat  = clCreateBuffer( context, CL_MEM_READ_ONLY,  Variant_Grid_sizeof(x_ssl_dat),  nullptr, &err );
  cl_mem device_y_ssl_dat  = clCreateBuffer( context, CL_MEM_READ_ONLY,  Variant_Grid_sizeof(y_ssl_dat),  nullptr, &err );
  cl_mem device_u_right    = clCreateBuffer( context, CL_MEM_READ_WRITE, Variant_Grid_sizeof(u_right),    nullptr, &err );
  cl_mem device_u_front    = clCreateBuffer( context, CL_MEM_READ_WRITE, Variant_Grid_sizeof(u_front),    nullptr, &err );
  cl_mem device_u_upper    = clCreateBuffer( context, CL_MEM_READ_WRITE, Variant_Grid_sizeof(u_upper),    nullptr, &err );

  // Below is a collection of arrays that map between the host and device memory objects
  const int num_device_read_buffers = 20;
  const int num_device_write_buffers = 4;
  // Device side grid buffers marked as read
  cl_mem* device_read_buffers[num_device_read_buffers] = {
    &device_fp,
    &device_dp,
    &device_et,
    &device_odp,
    &device_opp,

    &device_osp,
    &device_permxp,
    &device_permyp,
    &device_permzp,
    &device_pop,

    &device_pp,
    &device_rpp,
    &device_sp,
    &device_ss,
    &device_z_mult_dat,

    &device_x_ssl_dat,
    &device_y_ssl_dat,
    &device_u_right,
    &device_u_front,
    &device_u_upper
  };

  // Associated host side grid buffers for device side grid buffers marked as read
  double* host_read_buffers[num_device_read_buffers] = {
    Variant_Grid_data(fp),
    Variant_Grid_data(dp),
    Variant_Grid_data(et),
    Variant_Grid_data(odp),
    Variant_Grid_data(opp),

    Variant_Grid_data(osp),
    Variant_Grid_data(permxp),
    Variant_Grid_data(permyp),
    Variant_Grid_data(permzp),
    Variant_Grid_data(pop),

    Variant_Grid_data(pp),
    Variant_Grid_data(rpp),
    Variant_Grid_data(sp),
    Variant_Grid_data(ss),
    Variant_Grid_data(z_mult_dat),

    Variant_Grid_data(x_ssl_dat),
    Variant_Grid_data(y_ssl_dat),
    Variant_Grid_data(u_right),
    Variant_Grid_data(u_front),
    Variant_Grid_data(u_upper)
  };

  size_t sizeof_read_buffers[num_device_read_buffers] = {
    Variant_Grid_sizeof(fp),
    Variant_Grid_sizeof(dp),
    Variant_Grid_sizeof(et),
    Variant_Grid_sizeof(odp),
    Variant_Grid_sizeof(opp),
    Variant_Grid_sizeof(osp),
    Variant_Grid_sizeof(permxp),
    Variant_Grid_sizeof(permyp),
    Variant_Grid_sizeof(permzp),
    Variant_Grid_sizeof(pop),
    Variant_Grid_sizeof(pp),
    Variant_Grid_sizeof(rpp),
    Variant_Grid_sizeof(sp),
    Variant_Grid_sizeof(ss),
    Variant_Grid_sizeof(z_mult_dat),
    Variant_Grid_sizeof(x_ssl_dat),
    Variant_Grid_sizeof(y_ssl_dat),
    Variant_Grid_sizeof(u_right),
    Variant_Grid_sizeof(u_front),
    Variant_Grid_sizeof(u_upper)
  };

  // Device side grid buffers marked as write
  cl_mem* device_write_buffers[num_device_write_buffers] = {
    &device_fp,
    &device_vx,
    &device_vy,
    &device_vz
  };

  // Associated host side grid buffers for device side grid buffers marked as write
  double* host_write_buffers[num_device_write_buffers] = {
    Variant_Grid_data(fp),
    Variant_Grid_data(vx),
    Variant_Grid_data(vy),
    Variant_Grid_data(vz)
  };

  size_t sizeof_write_buffers[num_device_read_buffers] = {
    Variant_Grid_sizeof(fp),
    Variant_Grid_sizeof(vx),
    Variant_Grid_sizeof(vy),
    Variant_Grid_sizeof(vz)
  };

  Variant_Domain* domains[num_domains] = {
    domain,
    Variant_Grid_domain(fp),
    Variant_Grid_domain(dp),
    Variant_Grid_domain(et),
    Variant_Grid_domain(odp),

    Variant_Grid_domain(opp),
    Variant_Grid_domain(osp),
    Variant_Grid_domain(permxp),
    Variant_Grid_domain(permyp),
    Variant_Grid_domain(permzp),

    Variant_Grid_domain(pop),
    Variant_Grid_domain(pp),
    Variant_Grid_domain(rpp),
    Variant_Grid_domain(sp),
    Variant_Grid_domain(ss),

    Variant_Grid_domain(z_mult_dat),
    Variant_Grid_domain(x_ssl_dat),
    Variant_Grid_domain(y_ssl_dat),
    Variant_Grid_domain(vx),
    Variant_Grid_domain(vy),

    Variant_Grid_domain(vz),
    Variant_Grid_domain(u_right),
    Variant_Grid_domain(u_front),
    Variant_Grid_domain(u_upper)
  };

  char* domain_names[num_domains] = {
    (char*) STRINGIZE(domain),
    (char*) STRINGIZE(fp),
    (char*) STRINGIZE(dp),
    (char*) STRINGIZE(et),
    (char*) STRINGIZE(odp),

    (char*) STRINGIZE(opp),
    (char*) STRINGIZE(osp),
    (char*) STRINGIZE(permxp),
    (char*) STRINGIZE(permyp),
    (char*) STRINGIZE(permzp),

    (char*) STRINGIZE(pop),
    (char*) STRINGIZE(pp),
    (char*) STRINGIZE(rpp),
    (char*) STRINGIZE(sp),
    (char*) STRINGIZE(ss),

    (char*) STRINGIZE(z_mult_dat),
    (char*) STRINGIZE(x_ssl_dat),
    (char*) STRINGIZE(y_ssl_dat),
    (char*) STRINGIZE(vx),
    (char*) STRINGIZE(vy),

    (char*) STRINGIZE(vz),
    (char*) STRINGIZE(u_right),
    (char*) STRINGIZE(u_front),
    (char*) STRINGIZE(u_upper)
  };

  preamble_sources[0] = (char*)
  "typedef struct struct_Basic_Domain { \
      int nx, ny, nz; \
    } Basic_Domain; \
    typedef Basic_Domain Variant_Domain;";

  {
    const int max_str_len = 512;
    for( int i = 0; i < num_domains; i += 1){

      preamble_sources[1+i] = (char*) calloc(max_str_len, sizeof(char));
      sprintf(
        preamble_sources[i+1],
        "const int %1$s_domain_nx=%2$d;"\
        "const int %1$s_domain_ny=%3$d;"\
        "const int %1$s_domain_nz=%4$d;",
        domain_names[i],
        Variant_Domain_nx( domains[i] ),
        Variant_Domain_ny( domains[i] ),
        Variant_Domain_nz( domains[i] )
      );
    }
  }
    // typedef struct struct_Basic_Grid { \
    //   __global Basic_Domain* domain; \
    //   __global double* data; \
    // } Basic_Grid; \
    // typedef Basic_Domain Variant_Domain; \
    // typedef Basic_Grid Variant_Grid;";

  int kernel_idx;
  // Do baseline scientific kernel
  // NlFunctionEval:261 analogue
  kernel_idx = 0;
  // Setup arguments
  kernel_num_args[kernel_idx] = 7;
  kernel_args[kernel_idx] = (cl_mem**) calloc( kernel_num_args[kernel_idx], sizeof(cl_mem*) );
  // Set each argument's input device memory object
  kernel_args[kernel_idx][ 0] = &device_fp;
  kernel_args[kernel_idx][ 1] = &device_dp;
  kernel_args[kernel_idx][ 2] = &device_odp;
  kernel_args[kernel_idx][ 3] = &device_osp;
  kernel_args[kernel_idx][ 4] = &device_pop;
  kernel_args[kernel_idx][ 5] = &device_sp;
  kernel_args[kernel_idx][ 6] = &device_z_mult_dat;
  // Set kernel's global and local work-group sizes
  total_work[kernel_idx][0] = Variant_Domain_nx(domain)-2;
  total_work[kernel_idx][1] = Variant_Domain_ny(domain)-2;
  total_work[kernel_idx][2] = Variant_Domain_nz(domain)-2;
  // Kernel name and source text
  kernel_names[kernel_idx] = (char*)"function_NFE261";
  compute_kernel_sources[kernel_idx] = (char*) STRINGIZE(
    __kernel void function_NFE261(
      __global   double* fp,
      __constant double* dp,
      __constant double* odp,
      __constant double* osp,
      __constant double* pop,
      __constant double* sp,
      __constant double* z_mult_dat
    ){
      int x = get_global_id(0)+1;
      int y = get_global_id(1)+1;
      int z = get_global_id(2)+1;
      // printf("function_NFE261 (%d %d %d) %f\n", x,y,z, Variant_Grid_fast_access(sp,x,y,z));
      Variant_Grid_fast_access(fp, x,y,z) =
       (  Variant_Grid_fast_access(sp,  x,y,z)
        * Variant_Grid_fast_access(dp,  x,y,z)
        - Variant_Grid_fast_access(osp, x,y,z)
        * Variant_Grid_fast_access(odp, x,y,z)
       )
       * Variant_Grid_fast_access(pop, x,y,z)
       * Variant_Grid_fast_access(z_mult_dat, x,y,z);
    }
  );

  // NlFunctionEval:338 analogue
  kernel_idx = 1;
  // Setup arguments
  kernel_num_args[kernel_idx] = 9;
  kernel_args[kernel_idx] = (cl_mem**) calloc( kernel_num_args[kernel_idx], sizeof(cl_mem*) );
  // Set each argument's input device memory object
  kernel_args[kernel_idx][ 0] = &device_fp;
  kernel_args[kernel_idx][ 1] = &device_dp;
  kernel_args[kernel_idx][ 2] = &device_odp;
  kernel_args[kernel_idx][ 3] = &device_opp;
  kernel_args[kernel_idx][ 4] = &device_osp;
  kernel_args[kernel_idx][ 5] = &device_pp;
  kernel_args[kernel_idx][ 6] = &device_sp;
  kernel_args[kernel_idx][ 7] = &device_ss;
  kernel_args[kernel_idx][ 8] = &device_z_mult_dat;
  // Set kernel's global and local work-group sizes
  total_work[kernel_idx][0] = Variant_Domain_nx(domain)-2;
  total_work[kernel_idx][1] = Variant_Domain_ny(domain)-2;
  total_work[kernel_idx][2] = Variant_Domain_nz(domain)-2;
  // Kernel name and source text
  kernel_names[kernel_idx] = (char*)"function_NFE338";
  compute_kernel_sources[kernel_idx] = (char*) STRINGIZE(
    __kernel void function_NFE338(
      __global   double* fp,
      __constant double* dp,
      __constant double* odp,
      __constant double* opp,
      __constant double* osp,
      __constant double* pp,
      __constant double* sp,
      __constant double* ss,
      __constant double* z_mult_dat
    ){
      int x = get_global_id(0)+1;
      int y = get_global_id(1)+1;
      int z = get_global_id(2)+1;
      // printf("function_NFE338 (%d %d %d)\n", x,y,z);
      Variant_Grid_fast_access(fp, x,y,z) +=
          Variant_Grid_fast_access(ss, x,y,z)
        * Variant_Grid_fast_access(z_mult_dat, x,y,z)
        * (   Variant_Grid_fast_access(pp, x,y,z)
            * Variant_Grid_fast_access(sp, x,y,z)
            * Variant_Grid_fast_access(dp, x,y,z)
            - Variant_Grid_fast_access(opp, x,y,z)
            * Variant_Grid_fast_access(osp, x,y,z)
            * Variant_Grid_fast_access(odp, x,y,z)
          );
    }
  );

  // NlFunctionEval:416 analogue
  // Setup arguments
  kernel_idx = 2;
  kernel_num_args[kernel_idx] = 4;
  kernel_args[kernel_idx] = (cl_mem**) calloc( kernel_num_args[kernel_idx], sizeof(cl_mem*) );
  // Set each argument's input device memory object
  kernel_args[kernel_idx][0] = &device_fp;
  kernel_args[kernel_idx][1] = &device_et;
  kernel_args[kernel_idx][2] = &device_sp;
  kernel_args[kernel_idx][3] = &device_z_mult_dat;
  // Set kernel's global and local work-group sizes
  total_work[kernel_idx][0] = Variant_Domain_nx(domain)-2;
  total_work[kernel_idx][1] = Variant_Domain_ny(domain)-2;
  total_work[kernel_idx][2] = Variant_Domain_nz(domain)-2;
  // Kernel name and source text
  kernel_names[kernel_idx] = (char*)"function_NFE416";
  compute_kernel_sources[kernel_idx] = (char*) STRINGIZE(
    __kernel void function_NFE416(
      __global   double* fp,
      __constant double* et,
      __constant double* sp,
      __constant double* z_mult_dat
    ){
      int x = get_global_id(0)+1;
      int y = get_global_id(1)+1;
      int z = get_global_id(2)+1;
      // printf("function_NFE416 (%d %d %d)\n", x,y,z);
      Variant_Grid_fast_access(fp, x,y,z) -=
          Variant_Grid_fast_access(z_mult_dat, x,y,z)
        * (   Variant_Grid_fast_access(sp, x,y,z)
            * Variant_Grid_fast_access(et, x,y,z)
          );
    }
  );

  // NlFunctionEval:551 analogue
  // Setup arguments
  kernel_idx = 3;
  kernel_num_args[kernel_idx] = 16;
  kernel_args[kernel_idx] = (cl_mem**) calloc( kernel_num_args[kernel_idx], sizeof(cl_mem*) );
  // Set each argument's input device memory object
  kernel_args[kernel_idx][ 0] = &device_u_right;
  kernel_args[kernel_idx][ 1] = &device_u_front;
  kernel_args[kernel_idx][ 2] = &device_u_upper;
  kernel_args[kernel_idx][ 3] = &device_fp;
  kernel_args[kernel_idx][ 4] = &device_vx;
  kernel_args[kernel_idx][ 5] = &device_vy;
  kernel_args[kernel_idx][ 6] = &device_vz;
  kernel_args[kernel_idx][ 7] = &device_dp;
  kernel_args[kernel_idx][ 8] = &device_permxp;
  kernel_args[kernel_idx][ 9] = &device_permyp;
  kernel_args[kernel_idx][10] = &device_permzp;
  kernel_args[kernel_idx][11] = &device_pp;
  kernel_args[kernel_idx][12] = &device_rpp;
  kernel_args[kernel_idx][13] = &device_z_mult_dat;
  kernel_args[kernel_idx][14] = &device_x_ssl_dat;
  kernel_args[kernel_idx][15] = &device_y_ssl_dat;
  // Set kernel's global and local work-group sizes
  total_work[kernel_idx][0] = Variant_Domain_nx(domain)-2;
  total_work[kernel_idx][1] = Variant_Domain_ny(domain)-2;
  total_work[kernel_idx][2] = Variant_Domain_nz(domain)-2;
  // Kernel name and source text
  kernel_names[kernel_idx] = (char*)"function_NFE551";
  compute_kernel_sources[kernel_idx] = (char*) STRINGIZE(
    __kernel void function_NFE551(
      __global   double* u_right,
      __global   double* u_front,
      __global   double* u_upper,
      __global   double* fp,
      __global   double* vx,
      __global   double* vy,
      __global   double* vz,
      __constant double* dp,
      __constant double* permxp,
      __constant double* permyp,
      __constant double* permzp,
      __constant double* pp,
      __constant double* rpp,
      __constant double* z_mult_dat,
      __constant double* x_ssl_dat,
      __constant double* y_ssl_dat
    ){
      int x = get_global_id(0)+1;
      int y = get_global_id(1)+1;
      int z = get_global_id(2)+1;

      // printf("global (%u, %u, %u) size (%u, %u, %u), local (%u, %u, %u) size (%u, %u, %u) := (%d, %d, %d)\n",
      //   get_global_id(0), get_global_id(1), get_global_id(2),
      //   get_global_size(0), get_global_size(1), get_global_size(2),
      //   get_local_id(0), get_local_id(1), get_local_id(2),
      //   get_local_size(0), get_local_size(1), get_local_size(2),
      //   x,y,z
      // );
      // printf("function_NFE551 (%d %d %d)\n", x,y,z);
      double x_dir_g   = ArithmeticMean( Variant_Grid_fast_access( x_ssl_dat, x, y, 0), Variant_Grid_fast_access( x_ssl_dat, x+1,  y, 0 ) );
      double x_dir_g_c = ArithmeticMean( Variant_Grid_fast_access( x_ssl_dat, x, y, 0), Variant_Grid_fast_access( x_ssl_dat, x+1,  y, 0 ) );
      double y_dir_g   = ArithmeticMean( Variant_Grid_fast_access( y_ssl_dat, x, y, 0), Variant_Grid_fast_access( y_ssl_dat, x,  y+1, 0 ) );
      double y_dir_g_c = ArithmeticMean( Variant_Grid_fast_access( y_ssl_dat, x, y, 0), Variant_Grid_fast_access( y_ssl_dat, x,  y+1, 0 ) );

      double diff_right = Variant_Grid_fast_access(pp, x,y,z) - Variant_Grid_fast_access(pp, x+1, y,   z);
      double diff_front = Variant_Grid_fast_access(pp, x,y,z) - Variant_Grid_fast_access(pp, x,   y+1, z);

      double updir_right = diff_right * x_dir_g_c - x_dir_g;
      double updir_front = diff_front * y_dir_g_c - y_dir_g;

      double sep = ArithmeticMean( Variant_Grid_fast_access(z_mult_dat, x,y,z),  Variant_Grid_fast_access(z_mult_dat, x, y, z+1) );

      double lower_cond =
        Variant_Grid_fast_access(pp, x,y,z) / sep
      - (   Variant_Grid_fast_access(z_mult_dat, x,y,z)
          / ( Variant_Grid_fast_access(z_mult_dat, x,y,z)
            + Variant_Grid_fast_access(z_mult_dat, x, y, z+1)
            )
        )
      * Variant_Grid_fast_access(dp, x,y,z);

      double upper_cond =
        Variant_Grid_fast_access(pp, x, y, z+1) / sep
      + (   Variant_Grid_fast_access(z_mult_dat, x, y, z+1)
          / ( Variant_Grid_fast_access(z_mult_dat, x,y,z)
            + Variant_Grid_fast_access(z_mult_dat, x, y, z+1)
            )
        )
      * Variant_Grid_fast_access(dp, x,y,z+1);

      double diff_upper = lower_cond - upper_cond;

      double u_right_local =
        (
           Variant_Grid_fast_access(z_mult_dat, x,y,z)
         * HarmonicMean( Variant_Grid_fast_access(permxp, x,   y, z),
                         Variant_Grid_fast_access(permxp, x+1, y, z) )
         * diff_right * x_dir_g_c
         * UpstreamMean(
             updir_right,
             0.0,
             Variant_Grid_fast_access(rpp, x,   y, z) * Variant_Grid_fast_access(dp, x,   y, z),
             Variant_Grid_fast_access(rpp, x+1, y, z) * Variant_Grid_fast_access(dp, x+1, y, z)
           )
       ) + (
          Variant_Grid_fast_access(z_mult_dat, x,y,z)
         * HarmonicMean( Variant_Grid_fast_access(permxp, x,   y, z),
                         Variant_Grid_fast_access(permxp, x+1, y, z) )
         * (-x_dir_g)
         * UpstreamMean(
             updir_right,
             0.0,
             Variant_Grid_fast_access(rpp, x,   y, z) * Variant_Grid_fast_access(dp, x,   y, z),
             Variant_Grid_fast_access(rpp, x+1, y, z) * Variant_Grid_fast_access(dp, x+1, y, z)
           )
       );

       double u_front_local =
         (
            Variant_Grid_fast_access(z_mult_dat, x,y,z)
          * HarmonicMean( Variant_Grid_fast_access(permyp, x,   y, z),
                          Variant_Grid_fast_access(permyp, x+1, y, z) )
          * diff_front * x_dir_g_c
          * UpstreamMean(
              updir_front,
              0.0,
              Variant_Grid_fast_access(rpp, x,   y, z) * Variant_Grid_fast_access(dp, x,   y, z),
              Variant_Grid_fast_access(rpp, x+1, y, z) * Variant_Grid_fast_access(dp, x+1, y, z)
            )
        ) + (
           Variant_Grid_fast_access(z_mult_dat, x,y,z)
          * HarmonicMean( Variant_Grid_fast_access(permyp, x,   y, z),
                          Variant_Grid_fast_access(permyp, x+1, y, z) )
          * (-x_dir_g)
          * UpstreamMean(
              updir_front,
              0.0,
              Variant_Grid_fast_access(rpp, x,   y, z) * Variant_Grid_fast_access(dp, x,   y, z),
              Variant_Grid_fast_access(rpp, x+1, y, z) * Variant_Grid_fast_access(dp, x+1, y, z)
            )
        );

      double u_upper_local =
                  HarmonicMeanDZ(
                      Variant_Grid_fast_access(permzp,     x, y, z  ),
                      Variant_Grid_fast_access(permzp,     x, y, z+1),
                      Variant_Grid_fast_access(z_mult_dat, x, y, z  ),
                      Variant_Grid_fast_access(z_mult_dat, x, y, z+1)
                  )
                * diff_upper
                * UpstreamMean(
                    lower_cond,
                    upper_cond,
                    Variant_Grid_fast_access(rpp, x, y, z  ) * Variant_Grid_fast_access(dp, x, y, z  ),
                    Variant_Grid_fast_access(rpp, x, y, z+1) * Variant_Grid_fast_access(dp, x, y, z+1)
                  );
      Variant_Grid_fast_access(fp,      x, y, z) += u_right_local * u_front_local * u_upper_local;

      Variant_Grid_fast_access(vx,      x, y, z)  = u_right_local;
      Variant_Grid_fast_access(vy,      x, y, z)  = u_front_local;
      Variant_Grid_fast_access(vz,      x, y, z)  = u_upper_local;

      Variant_Grid_fast_access(u_right, x, y, z)  = u_right_local;
      Variant_Grid_fast_access(u_front, x, y, z)  = u_front_local;
      Variant_Grid_fast_access(u_upper, x, y, z)  = u_upper_local;
    }
  );


  // NlFunctionEval:551 reduction portion
  // Setup arguments
  kernel_idx = 4;
  kernel_num_args[kernel_idx] = 4;
  kernel_args[kernel_idx] = (cl_mem**) calloc( kernel_num_args[kernel_idx], sizeof(cl_mem*) );
  // Set each argument's input device memory object
  kernel_args[kernel_idx][0] = &device_u_upper;
  kernel_args[kernel_idx][1] = &device_u_front;
  kernel_args[kernel_idx][2] = &device_u_right,
  kernel_args[kernel_idx][3] = &device_fp;
  // Set kernel's global and local work-group sizes
  total_work[kernel_idx][0] = Variant_Domain_nx(domain)-1;
  total_work[kernel_idx][1] = Variant_Domain_ny(domain)-1;
  total_work[kernel_idx][2] = Variant_Domain_nz(domain)-1;
  // Kernel name and source text
  kernel_names[kernel_idx] = (char*)"function_NFE551_reduce";
  compute_kernel_sources[kernel_idx] = (char*) STRINGIZE(
    __kernel void function_NFE551_reduce(
      __constant double* u_upper,
      __constant double* u_front,
      __constant double* u_right,
      __global   double* fp
    ){
      int x = get_global_id(0)+1;
      int y = get_global_id(1)+1;
      int z = get_global_id(2)+1;

      Variant_Grid_fast_access(fp, x, y, z) =
            Variant_Grid_fast_access(fp,      x,   y,   z  )
        +   Variant_Grid_fast_access(u_right, x-1, y  , z  )
        + (-Variant_Grid_fast_access(u_front, x  , y-1, z  ))
        + (-Variant_Grid_fast_access(u_upper, x  , y  , z-1));
    }
  );


  // for( int i = 0; i < num_preambles; ++i ){
  //   printf("==============================\nPreamble %i\n------------------------------\n%s\n", i+1, preamble_sources[i]);
  // }
  //
  // for( int i = 0; i < num_compute_kernels; i += 1 )
  // {
  //   printf("==============================\n%s\n------------------------------\n%s\n", kernel_names[i], compute_kernel_sources[i]);
  // }

  /* Create the compute program from the source buffer */
  program = clCreateProgramWithSource(context, num_compute_kernels+num_preambles, (const char**)sources, nullptr, &err);

  char* options = (char*)"-cl-std=CL2.0";

  err = clBuildProgram( program, 1, &device_id, options, nullptr, nullptr);

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

  if( !program || err != CL_SUCCESS || build_status != CL_BUILD_SUCCESS){
    printf("Error: Failed to create compute program! (Error code %d)\n", err);
    exit(-1);
  }

  for( int i = 0; i < num_compute_kernels; ++i ){
    printf("clCreateKernel %s\n", kernel_names[i]);
    compute_kernels[i] = clCreateKernel(program, kernel_names[i], &err);
    if( !compute_kernels[i] || err != CL_SUCCESS ){
      printf("Error: Failed to create compute kernel \"%s\"! (Error code %d)\n", kernel_names[i], err);
      exit(1);
    }
  }

  // Copy read buffers
  for( int i = 0; i < num_device_read_buffers; i += 1 ){
    printf("clEnqueueWriteBuffer %d\n", i);
    err = clEnqueueWriteBuffer( commands, *device_read_buffers[i], CL_FALSE, 0, sizeof_read_buffers[i], host_read_buffers[i], 0, nullptr, nullptr );
    if( err != CL_SUCCESS ){
      printf("Error: Failed enqueue write-buffer command (Error code %d)\n", err);
      exit(1);
    }
  }

  // Setup args for each kernel
  for( int k = 0; k < num_compute_kernels; k += 1 ){
    for( int i = 0; i < kernel_num_args[k]; i += 1 ){
      printf("clSetKernelArg %s, %d\n", kernel_names[k], i);
      err = clSetKernelArg( compute_kernels[k], i, sizeof(cl_mem), kernel_args[k][i] );
      if( err != CL_SUCCESS ){
        printf("Error: Failed set kernel %d args %d (Error code %d)\n", k, i, err);
        exit(1);
      }
    }
  }

  // Wait for memory copies (reminder: writes currently non-blocking)
  clFinish( commands );

  // Invoke kernels
  for( int k = 0; k < num_compute_kernels; k += 1 ){
    printf("clEnqueueNDRangeKernel %s\n", kernel_names[k]);
    err = clEnqueueNDRangeKernel( commands, compute_kernels[k], 3, nullptr, total_work[k], nullptr, 0, nullptr, nullptr );
    if( err != CL_SUCCESS ){
      printf("Error: Failed to execute kernel %s! (Error code %d)\n", kernel_names[k], err);
      exit(1);
    }
  }

  // Wait for all kernels to complete
  clFinish( commands );

  // Read from out grids
  for( int i = 0; i < num_device_write_buffers; i += 1 ){
    printf("clEnqueueReadBuffer %d\n", i);
    err = clEnqueueReadBuffer( commands, *device_write_buffers[i], CL_FALSE, 0, sizeof_write_buffers[i], host_write_buffers[i], 0, nullptr, nullptr );
    if( err != CL_SUCCESS ){
      printf("Error: Failed enqueue read-buffer command (Error code %d)\n", err);
      exit(1);
    }
  }

  // Wait for all memory copies (reminder: reads currently non-blocking)
  clFinish( commands );

  // Free structures
  // TODO FREE EVERYTHING INCLUDING OPENCL STRUCTURES
  for( int k = 0; k < num_compute_kernels; k += 1 ){
    free(kernel_args[k]);
  }
  for( int i = 0; i < num_domains; i += 1 ){
    free( preamble_sources[1+i] );
  }

}
