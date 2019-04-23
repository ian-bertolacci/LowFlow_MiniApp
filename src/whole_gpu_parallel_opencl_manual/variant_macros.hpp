#ifndef VARIANT_MACROS_HPP
#define VARIANT_MACROS_HPP

#include <global_macros.hpp>

// Basic Domain
#define Variant_Domain_nx(domain) Basic_Domain_nx(domain)
#define Variant_Domain_ny(domain) Basic_Domain_ny(domain)
#define Variant_Domain_nz(domain) Basic_Domain_nz(domain)
#define Variant_Domain_idx(domain,x,y,z) Basic_Domain_idx(domain,x,y,z)
#define Variant_Domain_count(domain) (Variant_Domain_nx(domain)*Variant_Domain_ny(domain)*Variant_Domain_nz(domain))

// Basic Grid
#define Variant_Grid_data(grid) Basic_Grid_data(grid)
#define Variant_Grid_domain(grid) Basic_Grid_domain(grid)
#define Variant_Grid_sizeof(grid) (sizeof(double)*Variant_Domain_count(Variant_Grid_domain(grid)))
// The below should be unused
// #define Basic_Grid_access_index(grid,idx) (Basic_Grid_data(grid)[idx])
// #define Basic_Grid_access(grid,x,y,z) (Basic_Grid_access_index(grid,Basic_Domain_idx(Basic_Grid_domain(grid),x,y,z)))
#define Variant_Grid_idx(grid, x,y,z) Basic_Grid_idx(grid, x,y,z)
#define Variant_Grid_access_index(grid, idx) Basic_Grid_access_index(grid, idx)
#define Variant_Grid_access(grid, x,y,z) Basic_Grid_access(grid, x,y,z)
#define Variant_Grid_fast_access(grid, x,y,z) (grid[Basic_Domain_idx(grid##_domain,x,y,z)])

#define Variant_Domain_equal(domain_a, domain_b) Basic_Domain_equal(domain_a, domain_b)

#define Variant_Domain_loop_whole(domain, iter_x, iter_y, iter_z, body) Basic_Domain_loop_whole(domain, iter_x, iter_y, iter_z, body)

#define Variant_Domain_fast_loop_whole(domain, iter_x, iter_y, iter_z, body) \
  { \
    char* source = STRINGIZE(\
    for( iter_x = 0; iter_x < Basic_Domain_nx(domain); iter_x += 1 ){ \
      for( iter_y = 0; iter_y < Basic_Domain_ny(domain); iter_y += 1 ){ \
        for( iter_z = 0; iter_z < Basic_Domain_nz(domain); iter_z += 1 ){ \
          body; \
        } \
      } \
    }) \
    printf( "Source %s\n", source ); \
  }

#define Variant_Domain_loop_interior(domain, iter_x, iter_y, iter_z, body) Basic_Domain_loop_interior(domain, iter_x, iter_y, iter_z, body)

#define Variant_Domain_fast_loop_interior(domain, iter_x, iter_y, iter_z, body) \
  { \
    int err; \
    cl_device_id device_id;             /* compute device id */ \
    cl_context context;                 /* compute context */ \
    cl_command_queue commands;          /* compute command queue */ \
    cl_program program;                 /* compute program */ \
    cl_kernel kernel;                   /* compute kernel */ \
    cl_mem input;                       /* device memory used for the input array */ \
    cl_mem output;                      /* device memory used for the output array */ \
    /* Connect to a compute device */ \
    int gpu = 1; \
    err = clGetDeviceIDs(nullptr, gpu ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, 1, &device_id, nullptr); \
    if (err != CL_SUCCESS) \
    { \
        printf("Error: Failed to create a device group! (Error code %d)\n", err); \
    } \
    /* Create a compute context */ \
    context = clCreateContext(0, 1, &device_id, nullptr, nullptr, &err); \
    if (!context) \
    { \
        printf("Error: Failed to create a compute context! (Error code %d)\n", err); \
    } \
    /* Create a command commands */ \
    commands = clCreateCommandQueueWithProperties(context, device_id, nullptr, &err); \
    if (!commands) \
    { \
        printf("Error: Failed to create a command commands! (Error code %d)\n", err); \
    } \
    const char* KernelSource = STRINGIZE( \
      __kernel void __kernel__function__( ){ \
        body; \
      } \
    ); \
    printf( "Kernel source:\n%s\n", KernelSource ); \
    /* Create the compute program from the source buffer */ \
    program = clCreateProgramWithSource(context, 1, & KernelSource, nullptr, &err); \
    if (!program) \
    { \
      printf("Error: Failed to create compute program! (Error code %d)\n", err); \
    } \
    kernel = clCreateKernel(program, "square", &err);\
    if (!kernel || err != CL_SUCCESS)\
    {\
        printf("Error: Failed to create compute kernel! (Error code %d)\n", err);\
        exit(1);\
    }\
    /* Need host->device copy */ \
    /* Need exec */ \
    /* Need device->host copy */ \
  }

#endif
