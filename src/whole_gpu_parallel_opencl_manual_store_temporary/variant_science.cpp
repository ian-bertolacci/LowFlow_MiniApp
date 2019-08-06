#include <science.hpp>
#include <macros.hpp>
#include <configure.hpp>
#include <util.hpp>

#include <assert.h>
#include <CL/cl.h>
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

  START_TIMER( metrics->elapsed_setup );

  const int num_domains = 23;

  const int num_compute_kernels = 4 /*main kernels*/ + 1 /*reduce kernels*/;
  const int num_util_kernels = 0;
  const int num_preambles = 1 ;

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
    exit(OPENCL_DEVICE_CREATION_FAILURE);
  }

  /* Create a compute context */
  context = clCreateContext(0, 1, &device_id, nullptr, nullptr, &err);
  if( !context ){
    printf("Error: Failed to create a compute context! (Error code %d)\n", err);
    exit(OPENCL_DEVICE_CREATION_FAILURE);
  }

  /* Create a command commands */
  // Because is trying to monopolize the GPU programming space,
  // NVIDIA does not support OpenCL version > 1.2
  #if defined(CL_VERSION_2_0) or defined(CL_VERSION_2_1) or defined(CL_VERSION_2_2)
    commands = clCreateCommandQueueWithProperties(context, device_id, nullptr, &err);
  #elif defined(CL_VERSION_1_0) or defined(CL_VERSION_1_1) or defined(CL_VERSION_1_2)
    commands = clCreateCommandQueue(context, device_id, 0, &err);
  #endif
  if( !commands ){
    printf("Error: Failed to create a command commands! (Error code %d)\n", err);
    exit(OPENCL_DEVICE_CREATION_FAILURE);
  }

  if( !commands ){
    printf("Error: Failed to create a command commands! (Error code %d)\n", err);
    exit(OPENCL_DEVICE_CREATION_FAILURE);
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

  // Domains
  cl_mem device_fp_domain         = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_dp_domain         = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_et_domain         = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_odp_domain        = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_opp_domain        = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_osp_domain        = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_permxp_domain     = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_permyp_domain     = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_permzp_domain     = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_pop_domain        = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_pp_domain         = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_rpp_domain        = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_sp_domain         = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_ss_domain         = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_z_mult_dat_domain = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_x_ssl_dat_domain  = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_y_ssl_dat_domain  = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_vx_domain         = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_vy_domain         = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_vz_domain         = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_u_right_domain    = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_u_front_domain    = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_u_upper_domain    = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );

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

  Variant_Domain* host_domains[num_domains] = {
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

  cl_mem* device_domains[num_domains] = {
    &device_fp_domain,
    &device_dp_domain,
    &device_et_domain,
    &device_odp_domain,
    &device_opp_domain,

    &device_osp_domain,
    &device_permxp_domain,
    &device_permyp_domain,
    &device_permzp_domain,
    &device_pop_domain,

    &device_pp_domain,
    &device_rpp_domain,
    &device_sp_domain,
    &device_ss_domain,
    &device_z_mult_dat_domain,

    &device_x_ssl_dat_domain,
    &device_y_ssl_dat_domain,
    &device_vx_domain,
    &device_vy_domain,
    &device_vz_domain,

    &device_u_right_domain,
    &device_u_front_domain,
    &device_u_upper_domain
  };

  preamble_sources[0] = (char*) STRINGIZE(
    typedef struct struct_Basic_Domain {
      int nx;
      int ny;
      int nz;
    } Basic_Domain;
    typedef Basic_Domain Variant_Domain;
  );

  int kernel_idx, arg_idx;
  // Do baseline scientific kernel
  // NlFunctionEval:261 analogue
  kernel_idx = 0;
  arg_idx = 0;
  // Setup arguments
  kernel_num_args[kernel_idx] = 7*2; /* 7 array-domain pairs */
  kernel_args[kernel_idx] = (cl_mem**) calloc( kernel_num_args[kernel_idx], sizeof(cl_mem*) );
  // Set each argument's input device memory object
  kernel_args[kernel_idx][arg_idx++] = &device_fp;
  kernel_args[kernel_idx][arg_idx++] = &device_fp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_dp;
  kernel_args[kernel_idx][arg_idx++] = &device_dp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_odp;
  kernel_args[kernel_idx][arg_idx++] = &device_odp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_osp;
  kernel_args[kernel_idx][arg_idx++] = &device_osp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_pop;
  kernel_args[kernel_idx][arg_idx++] = &device_pop_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_sp;
  kernel_args[kernel_idx][arg_idx++] = &device_sp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_z_mult_dat;
  kernel_args[kernel_idx][arg_idx++] = &device_z_mult_dat_domain;
  // Set kernel's global and local work-group sizes
  total_work[kernel_idx][0] = Variant_Domain_nx(domain)-2;
  total_work[kernel_idx][1] = Variant_Domain_ny(domain)-2;
  total_work[kernel_idx][2] = Variant_Domain_nz(domain)-2;
  // Kernel name and source text
  kernel_names[kernel_idx] = (char*)"function_NFE261";
  compute_kernel_sources[kernel_idx] = (char*) STRINGIZE(
    __kernel void function_NFE261(
      __global   double*         fp,
      __constant Variant_Domain* fp_domain,
      __constant double*         dp,
      __constant Variant_Domain* dp_domain,
      __constant double*         odp,
      __constant Variant_Domain* odp_domain,
      __constant double*         osp,
      __constant Variant_Domain* osp_domain,
      __constant double*         pop,
      __constant Variant_Domain* pop_domain,
      __constant double*         sp,
      __constant Variant_Domain* sp_domain,
      __constant double*         z_mult_dat,
      __constant Variant_Domain* z_mult_dat_domain
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
  arg_idx = 0;
  // Setup arguments
  kernel_num_args[kernel_idx] = 9*2; /* 9 array-domain pairs */
  kernel_args[kernel_idx] = (cl_mem**) calloc( kernel_num_args[kernel_idx], sizeof(cl_mem*) );
  // Set each argument's input device memory object
  kernel_args[kernel_idx][arg_idx++] = &device_fp;
  kernel_args[kernel_idx][arg_idx++] = &device_fp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_dp;
  kernel_args[kernel_idx][arg_idx++] = &device_dp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_odp;
  kernel_args[kernel_idx][arg_idx++] = &device_odp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_opp;
  kernel_args[kernel_idx][arg_idx++] = &device_opp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_osp;
  kernel_args[kernel_idx][arg_idx++] = &device_osp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_pp;
  kernel_args[kernel_idx][arg_idx++] = &device_pp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_sp;
  kernel_args[kernel_idx][arg_idx++] = &device_sp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_ss;
  kernel_args[kernel_idx][arg_idx++] = &device_ss_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_z_mult_dat;
  kernel_args[kernel_idx][arg_idx++] = &device_z_mult_dat_domain;
  // Set kernel's global and local work-group sizes
  total_work[kernel_idx][0] = Variant_Domain_nx(domain)-2;
  total_work[kernel_idx][1] = Variant_Domain_ny(domain)-2;
  total_work[kernel_idx][2] = Variant_Domain_nz(domain)-2;
  // Kernel name and source text
  kernel_names[kernel_idx] = (char*)"function_NFE338";
  compute_kernel_sources[kernel_idx] = (char*) STRINGIZE(
    __kernel void function_NFE338(
      __global   double*         fp,
      __constant Variant_Domain* fp_domain,
      __constant double*         dp,
      __constant Variant_Domain* dp_domain,
      __constant double*         odp,
      __constant Variant_Domain* odp_domain,
      __constant double*         opp,
      __constant Variant_Domain* opp_domain,
      __constant double*         osp,
      __constant Variant_Domain* osp_domain,
      __constant double*         pp,
      __constant Variant_Domain* pp_domain,
      __constant double*         sp,
      __constant Variant_Domain* sp_domain,
      __constant double*         ss,
      __constant Variant_Domain* ss_domain,
      __constant double*         z_mult_dat,
      __constant Variant_Domain* z_mult_dat_domain
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
  arg_idx = 0;
  kernel_num_args[kernel_idx] = 4*2 /* 4 array-domain object pairs */;
  kernel_args[kernel_idx] = (cl_mem**) calloc( kernel_num_args[kernel_idx], sizeof(cl_mem*) );
  // Set each argument's input device memory object
  kernel_args[kernel_idx][arg_idx++] = &device_fp;
  kernel_args[kernel_idx][arg_idx++] = &device_fp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_et;
  kernel_args[kernel_idx][arg_idx++] = &device_et_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_sp;
  kernel_args[kernel_idx][arg_idx++] = &device_sp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_z_mult_dat;
  kernel_args[kernel_idx][arg_idx++] = &device_z_mult_dat_domain;
  // Set kernel's global and local work-group sizes
  total_work[kernel_idx][0] = Variant_Domain_nx(domain)-2;
  total_work[kernel_idx][1] = Variant_Domain_ny(domain)-2;
  total_work[kernel_idx][2] = Variant_Domain_nz(domain)-2;
  // Kernel name and source text
  kernel_names[kernel_idx] = (char*)"function_NFE416";
  compute_kernel_sources[kernel_idx] = (char*) STRINGIZE(
    __kernel void function_NFE416(
      __global   double*         fp,
      __constant Variant_Domain* fp_domain,
      __constant double*         et,
      __constant Variant_Domain* et_domain,
      __constant double*         sp,
      __constant Variant_Domain* sp_domain,
      __constant double*         z_mult_dat,
      __constant Variant_Domain* z_mult_dat_domain
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
  arg_idx = 0;
  kernel_num_args[kernel_idx] = 16*2; /* 16 array-domain pairs */
  kernel_args[kernel_idx] = (cl_mem**) calloc( kernel_num_args[kernel_idx], sizeof(cl_mem*) );
  // Set each argument's input device memory object
  kernel_args[kernel_idx][arg_idx++] = &device_u_right;
  kernel_args[kernel_idx][arg_idx++] = &device_u_right_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_u_front;
  kernel_args[kernel_idx][arg_idx++] = &device_u_front_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_u_upper;
  kernel_args[kernel_idx][arg_idx++] = &device_u_upper_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_fp;
  kernel_args[kernel_idx][arg_idx++] = &device_fp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_vx;
  kernel_args[kernel_idx][arg_idx++] = &device_vx_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_vy;
  kernel_args[kernel_idx][arg_idx++] = &device_vy_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_vz;
  kernel_args[kernel_idx][arg_idx++] = &device_vz_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_dp;
  kernel_args[kernel_idx][arg_idx++] = &device_dp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_permxp;
  kernel_args[kernel_idx][arg_idx++] = &device_permxp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_permyp;
  kernel_args[kernel_idx][arg_idx++] = &device_permyp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_permzp;
  kernel_args[kernel_idx][arg_idx++] = &device_permzp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_pp;
  kernel_args[kernel_idx][arg_idx++] = &device_pp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_rpp;
  kernel_args[kernel_idx][arg_idx++] = &device_rpp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_z_mult_dat;
  kernel_args[kernel_idx][arg_idx++] = &device_z_mult_dat_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_x_ssl_dat;
  kernel_args[kernel_idx][arg_idx++] = &device_x_ssl_dat_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_y_ssl_dat;
  kernel_args[kernel_idx][arg_idx++] = &device_y_ssl_dat_domain;
  // Set kernel's global and local work-group sizes
  total_work[kernel_idx][0] = Variant_Domain_nx(domain)-2;
  total_work[kernel_idx][1] = Variant_Domain_ny(domain)-2;
  total_work[kernel_idx][2] = Variant_Domain_nz(domain)-2;
  // Kernel name and source text
  kernel_names[kernel_idx] = (char*)"function_NFE551";
  compute_kernel_sources[kernel_idx] = (char*) STRINGIZE(
    __kernel void function_NFE551(
      __global   double*         u_right,
      __constant Variant_Domain* u_right_domain,
      __global   double*         u_front,
      __constant Variant_Domain* u_front_domain,
      __global   double*         u_upper,
      __constant Variant_Domain* u_upper_domain,
      __global   double*         fp,
      __constant Variant_Domain* fp_domain,
      __global   double*         vx,
      __constant Variant_Domain* vx_domain,
      __global   double*         vy,
      __constant Variant_Domain* vy_domain,
      __global   double*         vz,
      __constant Variant_Domain* vz_domain,
      __constant double*         dp,
      __constant Variant_Domain* dp_domain,
      __constant double*         permxp,
      __constant Variant_Domain* permxp_domain,
      __constant double*         permyp,
      __constant Variant_Domain* permyp_domain,
      __constant double*         permzp,
      __constant Variant_Domain* permzp_domain,
      __constant double*         pp,
      __constant Variant_Domain* pp_domain,
      __constant double*         rpp,
      __constant Variant_Domain* rpp_domain,
      __constant double*         z_mult_dat,
      __constant Variant_Domain* z_mult_dat_domain,
      __constant double*         x_ssl_dat,
      __constant Variant_Domain* x_ssl_dat_domain,
      __constant double*         y_ssl_dat,
      __constant Variant_Domain* y_ssl_dat_domain
    ){
      int x = get_global_id(0)+1;
      int y = get_global_id(1)+1;
      int z = get_global_id(2)+1;

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
  arg_idx = 0;
  kernel_num_args[kernel_idx] = 4*2; /* 4 array-domain pairs */
  kernel_args[kernel_idx] = (cl_mem**) calloc( kernel_num_args[kernel_idx], sizeof(cl_mem*) );
  // Set each argument's input device memory object
  kernel_args[kernel_idx][arg_idx++] = &device_fp;
  kernel_args[kernel_idx][arg_idx++] = &device_fp_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_u_upper;
  kernel_args[kernel_idx][arg_idx++] = &device_u_upper_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_u_front;
  kernel_args[kernel_idx][arg_idx++] = &device_u_front_domain;
  kernel_args[kernel_idx][arg_idx++] = &device_u_right;
  kernel_args[kernel_idx][arg_idx++] = &device_u_right_domain;
  // Set kernel's global and local work-group sizes
  total_work[kernel_idx][0] = Variant_Domain_nx(domain)-1;
  total_work[kernel_idx][1] = Variant_Domain_ny(domain)-1;
  total_work[kernel_idx][2] = Variant_Domain_nz(domain)-1;
  // Kernel name and source text
  kernel_names[kernel_idx] = (char*)"function_NFE551_reduce";
  compute_kernel_sources[kernel_idx] = (char*) STRINGIZE(
    __kernel void function_NFE551_reduce(
      __global   double*         fp,
      __constant Variant_Domain* fp_domain,
      __constant double*         u_upper,
      __constant Variant_Domain* u_upper_domain,
      __constant double*         u_front,
      __constant Variant_Domain* u_front_domain,
      __constant double*         u_right,
      __constant Variant_Domain* u_right_domain
    ){
      int x = get_global_id(0)+1;
      int y = get_global_id(1)+1;
      int z = get_global_id(2)+1;

      double fp_val = Variant_Grid_fast_access(fp, x, y, z);
      // In transformation from scatter to gather, we need to not gather from
      // locations that did not scatter to us.
      // Basic rules:
      //   1. Gather from axis i if iterator i is > 1;
      //   2. Gather from *only* axis i if iterator i is == to Ni,
      //      where Ni is the size in the i axis
      double u_right_val = ( (1 < x && y < fp_domain->ny && z < fp_domain->nz ) ?  Variant_Grid_fast_access(u_right, x-1, y  , z  ) : 0.0 );
      double u_front_val = ( (1 < y && x < fp_domain->nx && z < fp_domain->nz ) ? -Variant_Grid_fast_access(u_front, x  , y-1, z  ) : 0.0 );
      double u_upper_val = ( (1 < z && x < fp_domain->nx && y < fp_domain->ny ) ? -Variant_Grid_fast_access(u_upper, x  , y  , z-1) : 0.0 );
      Variant_Grid_fast_access(fp, x, y, z) = fp_val + u_right_val + u_front_val + u_upper_val;
    }
  );


  if( options.verbose ){
    for( int i = 0; i < num_preambles; ++i ){
      printf("==============================\nPreamble %i\n------------------------------\n%s\n", i+1, preamble_sources[i]);
    }

    for( int i = 0; i < num_compute_kernels; i += 1 ){
      printf("==============================\n%s\n------------------------------\n%s\n", kernel_names[i], compute_kernel_sources[i]);
    }
  }

  STOP_TIMER( metrics->elapsed_setup );


  START_TIMER( metrics->elapsed_compile );
  /* Create the compute program from the source buffer */
  program = clCreateProgramWithSource(context, num_compute_kernels+num_preambles, (const char**)sources, nullptr, &err);
  if( err != SUCCESS ){
    printf("Error: Failed to create compute program! (Error code %d)\n", err);
  }

  char* compile_options = (char*)
  #if defined(INTEL_OPENCL)
    "-cl-std=CL2.0 -O3 -cl-mad-enable -cl-fast-relaxed-math";
  #elif defined(NVIDIA_OPENCL)
    "-cl-std=CL2.0 -cl-mad-enable -cl-fast-relaxed-math";
  #else
    "-cl-std=CL2.0";
  #endif

  err = clBuildProgram( program, 1, &device_id, compile_options, nullptr, nullptr);


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
      if( options.verbose ) printf("Build Status: CL_BUILD_SUCCESS\n");
      break;
    }
    case CL_BUILD_IN_PROGRESS: {
      printf("Build Status: CL_BUILD_IN_PROGRESS\n");
      break;
    }
  }

  if( options.verbose || !program || err != CL_SUCCESS || build_status != CL_BUILD_SUCCESS ){
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

    if( !program || err != CL_SUCCESS || build_status != CL_BUILD_SUCCESS ){
      printf("Error: Failed to build compute program! (Error code %d)\n", err);
      exit(OPENCL_COMPILATION_FAILURE);
    }
  }

  for( int i = 0; i < num_compute_kernels; ++i ){
    if( options.verbose ) printf("clCreateKernel %s\n", kernel_names[i]);
    compute_kernels[i] = clCreateKernel(program, kernel_names[i], &err);
    if( !compute_kernels[i] || err != CL_SUCCESS ){
      printf("Error: Failed to create compute kernel \"%s\"! (Error code %d)\n", kernel_names[i], err);
      exit(OPENCL_KERNEL_CREATION_FAILURE);
    }
  }

  STOP_TIMER( metrics->elapsed_compile );

  START_TIMER( metrics->elapsed_copy_host_to_device );

  // Copy read buffers
  for( int i = 0; i < num_device_read_buffers; i += 1 ){
    if( options.verbose ) printf("clEnqueueWriteBuffer read array %d\n", i);
    err = clEnqueueWriteBuffer( commands, *device_read_buffers[i], CL_FALSE, 0, sizeof_read_buffers[i], host_read_buffers[i], 0, nullptr, nullptr );
    if( err != CL_SUCCESS ){
      printf("Error: Failed enqueue write-buffer command (Error code %d)\n", err);
      exit(OPENCL_BUFFER_WRITE_FAILURE);
    }
  }

  // Copy Domains
  for( int i = 0; i < num_domains; i += 1 ){
    if( options.verbose ) printf("clEnqueueWriteBuffer domain %d (%p)\n", i, host_domains[i] );
    err = clEnqueueWriteBuffer( commands, *device_domains[i], CL_FALSE, 0, sizeof(Variant_Domain), host_domains[i], 0, nullptr, nullptr );
    if( err != CL_SUCCESS ){
      printf("Error: Failed enqueue write-buffer command (Error code %d)\n", err);
      exit(OPENCL_BUFFER_WRITE_FAILURE);
    }
  }

  START_TIMER( metrics->elapsed_exec_setup );

  // Setup args for each kernel
  for( int k = 0; k < num_compute_kernels; k += 1 ){
    for( int i = 0; i < kernel_num_args[k]; i += 1 ){
      if( options.verbose ) printf("clSetKernelArg %s, %d\n", kernel_names[k], i);
      err = clSetKernelArg( compute_kernels[k], i, sizeof(cl_mem), kernel_args[k][i] );
      if( err != CL_SUCCESS ){
        printf("Error: Failed set kernel %d args %d (Error code %d)\n", k, i, err);
        exit(OPENCL_KERNEL_SETUP_FAILURE);
      }
    }
  }


  // Wait for memory copies (reminder: writes currently non-blocking)
  clFinish( commands );
  STOP_TIMER( metrics->elapsed_copy_host_to_device );
  STOP_TIMER( metrics->elapsed_exec_setup );

  START_TIMER( metrics->elapsed_exec );
  #ifdef ENABLE_VARIANT_METRICS
  double* timers[num_compute_kernels] = { &metrics->elapsed_216, &metrics->elapsed_338, &metrics->elapsed_416, &metrics->elapsed_551, &metrics->elapsed_551_reduce };
  #endif
  // Invoke kernels
  for( int k = 0; k < num_compute_kernels; k += 1 ){
    if( options.verbose ) printf("clEnqueueNDRangeKernel %s\n", kernel_names[k]);
    START_TIMER( *timers[k] );
    err = clEnqueueNDRangeKernel( commands, compute_kernels[k], 3, nullptr, total_work[k], nullptr, 0, nullptr, nullptr );
    STOP_TIMER( *timers[k] );
    if( err != CL_SUCCESS ){
      printf("Error: Failed to execute kernel %s! (Error code %d)\n", kernel_names[k], err);
      exit(OPENCL_KERNEL_EXECUTION_FAILURE);
    }
  }

  // Wait for all kernels to complete
  clFinish( commands );
  STOP_TIMER( metrics->elapsed_exec );

  START_TIMER( metrics->elapsed_copy_device_to_host );
  // Read from out grids
  for( int i = 0; i < num_device_write_buffers; i += 1 ){
    if( options.verbose ) printf("clEnqueueReadBuffer %d\n", i);
    err = clEnqueueReadBuffer( commands, *device_write_buffers[i], CL_FALSE, 0, sizeof_write_buffers[i], host_write_buffers[i], 0, nullptr, nullptr );
    if( err != CL_SUCCESS ){
      printf("Error: Failed enqueue read-buffer command (Error code %d)\n", err);
      exit(OPENCL_BUFFER_READ_FAILURE);
    }
  }

  // Wait for all memory copies (reminder: reads currently non-blocking)
  clFinish( commands );
  STOP_TIMER( metrics->elapsed_copy_device_to_host );

  START_TIMER( metrics->elapsed_teardown );
  // Free structures
  // TODO FREE EVERYTHING INCLUDING OPENCL STRUCTURES
  for( int k = 0; k < num_compute_kernels; k += 1 ){
    free(kernel_args[k]);
  }
  // for( int i = 0; i < num_domains; i += 1 ){
  //   free( preamble_sources[1+i] );
  // }
  STOP_TIMER( metrics->elapsed_teardown );
}
