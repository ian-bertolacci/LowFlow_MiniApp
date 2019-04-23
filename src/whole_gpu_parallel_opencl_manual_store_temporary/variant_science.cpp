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

  const int num_compute_kernels = 4;
  const int num_util_kernels = 0;
  const int num_preambles = 1;

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
  // global domain size for kernal invocation
  size_t global[3];
  // local domain size for kernal invocation
  size_t local[3];

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


  // Note: domain is not a grid, and so it must be treated differently
  // TODO: Set this all ups so it doesnt have to be (void* and sizeof array)
  cl_mem device_domain            = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_fp_domain         = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_vx_domain         = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_vy_domain         = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
  cl_mem device_vz_domain         = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(Variant_Domain), nullptr, &err );
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

  // Below is a collection of arrays that map between the host and device memory objects
  const int num_device_read_buffers = 17;
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
    &device_y_ssl_dat
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
    Variant_Grid_data(y_ssl_dat)
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
    Variant_Grid_sizeof(y_ssl_dat)
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

  // Below is a collection of arrays that map between the host and device memory objects
  const int num_domains = 21;
  // Device side grid buffers marked as read
  cl_mem* device_domains[num_domains] = {
    &device_domain,
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
  };

  Variant_Domain* host_domains[num_domains] = {
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
    Variant_Grid_domain(vz)
  };

  preamble_sources[0] = (char*)
  "typedef struct struct_Basic_Domain { \
      int nx, ny, nz; \
    } Basic_Domain; \
    typedef Basic_Domain Variant_Domain;";
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
  kernel_num_args[kernel_idx] = 15;
  kernel_args[kernel_idx] = (cl_mem**) calloc( kernel_num_args[kernel_idx], sizeof(cl_mem*) );
  kernel_args[kernel_idx][ 0] = &device_domain;
  kernel_args[kernel_idx][ 1] = &device_fp;
  kernel_args[kernel_idx][ 2] = &device_fp_domain;
  kernel_args[kernel_idx][ 3] = &device_dp;
  kernel_args[kernel_idx][ 4] = &device_dp_domain;
  kernel_args[kernel_idx][ 5] = &device_odp;
  kernel_args[kernel_idx][ 6] = &device_odp_domain;
  kernel_args[kernel_idx][ 7] = &device_osp;
  kernel_args[kernel_idx][ 8] = &device_osp_domain;
  kernel_args[kernel_idx][ 9] = &device_pop;
  kernel_args[kernel_idx][10] = &device_pop_domain;
  kernel_args[kernel_idx][11] = &device_sp;
  kernel_args[kernel_idx][12] = &device_sp_domain;
  kernel_args[kernel_idx][13] = &device_z_mult_dat;
  kernel_args[kernel_idx][14] = &device_z_mult_dat_domain;
  kernel_names[kernel_idx] = (char*)"function_NFE261";
  compute_kernel_sources[kernel_idx] = (char*) STRINGIZE(
    __kernel void function_NFE261(
      __global Variant_Domain* domain,
      __global double*         fp,
      __global Variant_Domain* fp_domain,
      __global double*         dp,
      __global Variant_Domain* dp_domain,
      __global double*         odp,
      __global Variant_Domain* odp_domain,
      __global double*         osp,
      __global Variant_Domain* osp_domain,
      __global double*         pop,
      __global Variant_Domain* pop_domain,
      __global double*         sp,
      __global Variant_Domain* sp_domain,
      __global double*         z_mult_dat,
      __global Variant_Domain* z_mult_dat_domain
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
  kernel_num_args[kernel_idx] = 19;
  kernel_args[kernel_idx] = (cl_mem**) calloc( kernel_num_args[kernel_idx], sizeof(cl_mem*) );
  kernel_args[kernel_idx][ 0] = &device_domain;
  kernel_args[kernel_idx][ 1] = &device_fp;
  kernel_args[kernel_idx][ 2] = &device_fp_domain;
  kernel_args[kernel_idx][ 3] = &device_dp;
  kernel_args[kernel_idx][ 4] = &device_dp_domain;
  kernel_args[kernel_idx][ 5] = &device_odp;
  kernel_args[kernel_idx][ 6] = &device_odp_domain;
  kernel_args[kernel_idx][ 7] = &device_opp;
  kernel_args[kernel_idx][ 8] = &device_opp_domain;
  kernel_args[kernel_idx][ 9] = &device_osp;
  kernel_args[kernel_idx][10] = &device_osp_domain;
  kernel_args[kernel_idx][11] = &device_pp;
  kernel_args[kernel_idx][12] = &device_pp_domain;
  kernel_args[kernel_idx][13] = &device_sp;
  kernel_args[kernel_idx][14] = &device_sp_domain;
  kernel_args[kernel_idx][15] = &device_ss;
  kernel_args[kernel_idx][16] = &device_ss_domain;
  kernel_args[kernel_idx][17] = &device_z_mult_dat;
  kernel_args[kernel_idx][18] = &device_z_mult_dat_domain;
  kernel_names[kernel_idx] = (char*)"function_NFE338";
  compute_kernel_sources[kernel_idx] = (char*) STRINGIZE(
    __kernel void function_NFE338(
      __global Variant_Domain* domain,
      __global double*         fp,
      __global Variant_Domain* fp_domain,
      __global double*         dp,
      __global Variant_Domain* dp_domain,
      __global double*         odp,
      __global Variant_Domain* odp_domain,
      __global double*         opp,
      __global Variant_Domain* opp_domain,
      __global double*         osp,
      __global Variant_Domain* osp_domain,
      __global double*         pp,
      __global Variant_Domain* pp_domain,
      __global double*         sp,
      __global Variant_Domain* sp_domain,
      __global double*         ss,
      __global Variant_Domain* ss_domain,
      __global double*         z_mult_dat,
      __global Variant_Domain* z_mult_dat_domain
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
  kernel_idx = 2;
  kernel_num_args[kernel_idx] = 9;
  kernel_args[kernel_idx] = (cl_mem**) calloc( kernel_num_args[kernel_idx], sizeof(cl_mem*) );
  kernel_args[kernel_idx][0] = &device_domain;
  kernel_args[kernel_idx][1] = &device_fp;
  kernel_args[kernel_idx][2] = &device_fp_domain;
  kernel_args[kernel_idx][3] = &device_et;
  kernel_args[kernel_idx][4] = &device_et_domain;
  kernel_args[kernel_idx][5] = &device_sp;
  kernel_args[kernel_idx][6] = &device_sp_domain;
  kernel_args[kernel_idx][7] = &device_z_mult_dat;
  kernel_args[kernel_idx][8] = &device_z_mult_dat_domain;
  kernel_names[kernel_idx] = (char*)"function_NFE416";
  compute_kernel_sources[kernel_idx] = (char*) STRINGIZE(
    __kernel void function_NFE416(
      __global Variant_Domain* domain,
      __global double*         fp,
      __global Variant_Domain* fp_domain,
      __global double*         et,
      __global Variant_Domain* et_domain,
      __global double*         sp,
      __global Variant_Domain* sp_domain,
      __global double*         z_mult_dat,
      __global Variant_Domain* z_mult_dat_domain
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
  kernel_idx = 3;
  kernel_num_args[kernel_idx] = 27;
  kernel_args[kernel_idx] = (cl_mem**) calloc( kernel_num_args[kernel_idx], sizeof(cl_mem*) );
  kernel_args[kernel_idx][ 0] = &device_domain;
  kernel_args[kernel_idx][ 1] = &device_fp;
  kernel_args[kernel_idx][ 2] = &device_fp_domain;
  kernel_args[kernel_idx][ 3] = &device_vx;
  kernel_args[kernel_idx][ 4] = &device_vx_domain;
  kernel_args[kernel_idx][ 5] = &device_vy;
  kernel_args[kernel_idx][ 6] = &device_vy_domain;
  kernel_args[kernel_idx][ 7] = &device_vz;
  kernel_args[kernel_idx][ 8] = &device_vz_domain;
  kernel_args[kernel_idx][ 9] = &device_dp;
  kernel_args[kernel_idx][10] = &device_dp_domain;
  kernel_args[kernel_idx][11] = &device_permxp;
  kernel_args[kernel_idx][12] = &device_permxp_domain;
  kernel_args[kernel_idx][13] = &device_permyp;
  kernel_args[kernel_idx][14] = &device_permyp_domain;
  kernel_args[kernel_idx][15] = &device_permzp;
  kernel_args[kernel_idx][16] = &device_permzp_domain;
  kernel_args[kernel_idx][17] = &device_pp;
  kernel_args[kernel_idx][18] = &device_pp_domain;
  kernel_args[kernel_idx][19] = &device_rpp;
  kernel_args[kernel_idx][20] = &device_rpp_domain;
  kernel_args[kernel_idx][21] = &device_z_mult_dat;
  kernel_args[kernel_idx][22] = &device_z_mult_dat_domain;
  kernel_args[kernel_idx][23] = &device_x_ssl_dat;
  kernel_args[kernel_idx][24] = &device_x_ssl_dat_domain;
  kernel_args[kernel_idx][25] = &device_y_ssl_dat;
  kernel_args[kernel_idx][26] = &device_y_ssl_dat_domain;
  kernel_names[kernel_idx] = (char*)"function_NFE551";
  compute_kernel_sources[kernel_idx] = (char*) STRINGIZE(
    __kernel void function_NFE551(
      __global Variant_Domain* domain,
      __global double*         fp,
      __global Variant_Domain* fp_domain,
      __global double*         vx,
      __global Variant_Domain* vx_domain,
      __global double*         vy,
      __global Variant_Domain* vy_domain,
      __global double*         vz,
      __global Variant_Domain* vz_domain,
      __global double*         dp,
      __global Variant_Domain* dp_domain,
      __global double*         permxp,
      __global Variant_Domain* permxp_domain,
      __global double*         permyp,
      __global Variant_Domain* permyp_domain,
      __global double*         permzp,
      __global Variant_Domain* permzp_domain,
      __global double*         pp,
      __global Variant_Domain* pp_domain,
      __global double*         rpp,
      __global Variant_Domain* rpp_domain,
      __global double*         z_mult_dat,
      __global Variant_Domain* z_mult_dat_domain,
      __global double*         x_ssl_dat,
      __global Variant_Domain* x_ssl_dat_domain,
      __global double*         y_ssl_dat,
      __global Variant_Domain* y_ssl_dat_domain
    ){
      int x = get_global_id(0)+1;
      int y = get_global_id(1)+1;
      int z = get_global_id(2)+1;
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

      double u_right =
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

       double u_front =
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

      double u_upper =
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

      Variant_Grid_fast_access(vx, x,y,z) = u_right;
      Variant_Grid_fast_access(vy, x,y,z) = u_front;
      Variant_Grid_fast_access(vz, x,y,z) = u_upper;

      Variant_Grid_fast_access(fp, x  , y  , z  ) += u_right * u_front * u_upper;
      Variant_Grid_fast_access(fp, x+1, y  , z  ) += u_right;
      Variant_Grid_fast_access(fp, x  , y+1, z  ) -= u_front;
      Variant_Grid_fast_access(fp, x  , y  , z+1) -= u_upper;
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
  }

  for( int i = 0; i < num_compute_kernels; ++i ){
    compute_kernels[i] = clCreateKernel(program, kernel_names[i], &err);
    if( !compute_kernels[i] || err != CL_SUCCESS ){
      printf("Error: Failed to create compute kernel \"%s\"! (Error code %d)\n", kernel_names[i], err);
      exit(1);
    }
  }

  // Copy read buffers
  // TODO: Make non-blocking.
  for( int i = 0; i < num_device_read_buffers; i += 1 ){
    err = clEnqueueWriteBuffer( commands, *device_read_buffers[i], CL_TRUE, 0, sizeof_read_buffers[i], host_read_buffers[i], 0, nullptr, nullptr );
    if( err != CL_SUCCESS ){
      printf("Error: Failed enqueue write-buffer command (Error code %d)\n", err);
      exit(1);
    }

    // cl_mem_flags flag;
    // err = clGetMemObjectInfo( *device_read_buffers[i], CL_MEM_FLAGS, sizeof(cl_mem_flags), &flag, nullptr );
    // if( err != CL_SUCCESS ){
    //   printf("Error: clGetMemObjectInfo (Error code %d)\n", err);
    //   exit(1);
    // }
    // cl_mem_flags flags[6] =  { CL_MEM_READ_WRITE, CL_MEM_WRITE_ONLY, CL_MEM_READ_ONLY, CL_MEM_USE_HOST_PTR, CL_MEM_ALLOC_HOST_PTR, CL_MEM_COPY_HOST_PTR };
    // char* flag_names[6] =  { "CL_MEM_READ_WRITE", "CL_MEM_WRITE_ONLY", "CL_MEM_READ_ONLY", "CL_MEM_USE_HOST_PTR", "CL_MEM_ALLOC_HOST_PTR", "CL_MEM_COPY_HOST_PTR" };
    // for( int j = 0; j < 6; ++j ){
    //   if( flag & flags[j] ) printf("%s ", flag_names[j]);
    // }
    // printf("\n");
  }

  // Copy Domains
  for( int i = 0; i < num_domains; i += 1 ){
    err = clEnqueueWriteBuffer( commands, *device_domains[i], CL_TRUE, 0, sizeof(Variant_Domain), host_domains[i], 0, nullptr, nullptr );
    if( err != CL_SUCCESS ){
      printf("Error: Failed enqueue write-buffer command (Error code %d)\n", err);
      exit(1);
    }
  }

  // Setup args for each kernel
  for( int k = 0; k < num_compute_kernels; k += 1 ){
    for( int i = 0; i < kernel_num_args[k]; i += 1 ){
      err = clSetKernelArg( compute_kernels[k], i, sizeof(cl_mem), kernel_args[k][i] );
      if( err != CL_SUCCESS ){
        printf("Error: Failed set kernel %d args %d (Error code %d)\n", k, i, err);
        exit(1);
      }
    }
  }

  // Get the maximum work group size for executing the kernel on the device
  for( int k = 0; k < num_compute_kernels; k += 1 ){
    err = clGetKernelWorkGroupInfo( compute_kernels[k], device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(work_group_size[k]), &work_group_size[k], nullptr );
    if( err != CL_SUCCESS ){
        printf("Error: Failed to retrieve kernel %d work group info! %d\n", k, err);
        exit(1);
    }
  }


  global[0] = 1;//Variant_Domain_nx(domain)-2;
  global[1] = 1;//Variant_Domain_ny(domain)-2;
  global[2] = 1;//Variant_Domain_nz(domain)-2;
  for( int i = 0; i < 3; i += 1 ) Variant_Domain_nx(domain)-2;

  // TODO Is this necessary?
  clFinish( commands );

  for( int k = 0; k < num_compute_kernels; k += 1 ){
    err = clEnqueueNDRangeKernel( commands, compute_kernels[k], 3, nullptr, global, local, 0, nullptr, nullptr );
    clFinish( commands );
  }


  for( int i = 0; i < num_device_write_buffers; i += 1 ){
    err = clEnqueueReadBuffer( commands, *device_write_buffers[i], CL_TRUE, 0, sizeof_write_buffers[i], host_write_buffers[i], 0, nullptr, nullptr );
    if( err != CL_SUCCESS ){
      printf("Error: Failed enqueue write-buffer command (Error code %d)\n", err);
      exit(1);
    }
  }

  clFinish( commands );


}
