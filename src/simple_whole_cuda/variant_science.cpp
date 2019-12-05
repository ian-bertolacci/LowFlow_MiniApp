#include <science.hpp>
#include <macros.hpp>
#include <configure.hpp>
#include <util.hpp>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

//NlFunctionEval:216 CUDA kernel 
__global__ void NLFE216Kernel(int xmax, int ymax, int zmax, Variant_Grid *fp, Variant_Grid *sp, Variant_Grid *dp, Variant_Grid *osp, Variant_Grid *odp, Variant_Grid *pop, Variant_Grid *z_mult_dat) {
  //Get position in kernel grid
  //+1 is added to coordinates to account for offset in each direction
  int x = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int y = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int z = blockIdx.z * blockDim.z + threadIdx.z + 1;
  
  //Make sure that the thread is in bounds
  if (x == 0 || x > xmax || y == 0 || y > ymax || z == 0 || z > zmax) {
    return;
  }

  //Perform variant grid access
  Variant_Grid_access(fp,x,y,z) = (Variant_Grid_access(sp,x,y,z) * Variant_Grid_access(dp,x,y,z) - Variant_Grid_access(osp,x,y,z) * Variant_Grid_access(odp, x,y,z)) * Variant_Grid_access(pop, x,y,z) * Variant_Grid_access(z_mult_dat, x,y,z);
}

//NlFunctionEval:338 CUDA kernel 
__global__ void NLFE338Kernel(int xmax, int ymax, int zmax, Variant_Grid *fp, Variant_Grid *ss, Variant_Grid *z_mult_dat, Variant_Grid *pp, Variant_Grid *sp, Variant_Grid *dp, Variant_Grid *opp, Variant_Grid *osp, Variant_Grid *odp) {
  //Get position in kernel grid
  //+1 is added to coordinates to account for offset in each direction
  int x = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int y = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int z = blockIdx.z * blockDim.z + threadIdx.z + 1;

  //Make sure that the thread is in bounds
  if (x == 0 || x > xmax || y == 0 || y > ymax || z == 0 || z > zmax) {
    return;
  }

  //Perform variant grid access
  Variant_Grid_access(fp, x,y,z) += Variant_Grid_access(ss, x,y,z) * Variant_Grid_access(z_mult_dat, x,y,z) * (Variant_Grid_access(pp, x,y,z) * Variant_Grid_access(sp, x,y,z) * Variant_Grid_access(dp, x,y,z) - Variant_Grid_access(opp, x,y,z) * Variant_Grid_access(osp, x,y,z) * Variant_Grid_access(odp, x,y,z));
}

//NlFunctionEval:416 CUDA kernel 
__global__ void NLFE416Kernel(int xmax, int ymax, int zmax, Variant_Grid *fp, Variant_Grid *z_mult_dat, Variant_Grid *sp, Variant_Grid *et) {
  //Get position in kernel grid
  //+1 is added to coordinates to account for offset in each direction
  int x = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int y = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int z = blockIdx.z * blockDim.z + threadIdx.z + 1;

  //Make sure that the thread is in bounds
  if (x == 0 || x > xmax || y == 0 || y > ymax || z == 0 || z > zmax) {
    return;
  }

  //Perform variant grid access
  Variant_Grid_access(fp, x,y,z) -= Variant_Grid_access(z_mult_dat, x,y,z) * (Variant_Grid_access(sp, x,y,z) * Variant_Grid_access(et, x,y,z));
}

//NlFunctionEval:551 CUDA kernel 
__global__ void NLFE551Kernel(int xmax, int ymax, int zmax, Variant_Grid *x_ssl_dat, Variant_Grid *y_ssl_dat, Variant_Grid *pp, Variant_Grid *z_mult_dat, Variant_Grid *dp, Variant_Grid *permxp, Variant_Grid *rpp, Variant_Grid *permyp, Variant_Grid *permzp, Variant_Grid *vx, Variant_Grid *vy, Variant_Grid *vz, Variant_Grid *u_right, Variant_Grid *u_front, Variant_Grid *u_upper, Variant_Grid *fp) {
  //Get position in kernel grid
  //+1 is added to coordinates to account for offset in each direction
  int x = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int y = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int z = blockIdx.z * blockDim.z + threadIdx.z + 1;

  //Make sure that the thread is in bounds
  if (x == 0 || x > xmax || y == 0 || y > ymax || z == 0 || z > zmax) {
    return;
  }

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

  /*
  NOTE!
  Originally the harmonic mean was
  PMean( Variant_Grid_access(pp, x,y,z),
          Variant_Grid_access(pp, x+1, y, z),
          permxp[ip],
          permxp[ip + 1]
  )
  However! PMean(a,b,c,d) in parflow is simply HarmonicMean(c, d)
  so the first two terms have been removed
  */

  double u_right_val =
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

    double u_front_val =
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

  double u_upper_val =
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

  Variant_Grid_access(vx, x, y, z) = u_right_val;
  Variant_Grid_access(vy, x, y, z) = u_front_val;
  Variant_Grid_access(vz, x, y, z) = u_upper_val;

  /* Note:
  A reasonable person might ask:
  "Ian, Why are we not just using vx, vy, and vz for the value of
  u_right, u_front, u_upper in the reduction?", which is totally valid.
  The answer is that, in the original parflow code, vx, vy, and vz are
  not set to the u_*_value, but rather a constant expression involving
  u_*_value. Now, maybe it's worth it to delay that computation until
  later to reduce the temporary storage requirement at the cost of
  having to reload that data, but it's complicated.
  */

  Variant_Grid_access(u_right, x, y, z) = u_right_val;
  Variant_Grid_access(u_front, x, y, z) = u_front_val;
  Variant_Grid_access(u_upper, x, y, z) = u_upper_val;

  Variant_Grid_access(fp, x, y, z) += u_right_val * u_front_val * u_upper_val;
}

__global__ void NLFE551ReductionKernel(int xmax, int ymax, int zmax, Variant_Grid *u_right, Variant_Grid *u_front, Variant_Grid *u_upper, Variant_Grid *fp) {
  //Get position in kernel grid
  //+1 is added to coordinates to account for offset in each direction
  int x = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int y = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int z = blockIdx.z * blockDim.z + threadIdx.z + 1;

  //Make sure that the thread is in bounds
  if (x == 0 || x > xmax || y == 0 || y > ymax || z == 0 || z > zmax) {
    return;
  }

  double u_right_val =  Variant_Grid_access(u_right, x-1, y  , z  );
  double u_front_val = -Variant_Grid_access(u_front, x  , y-1, z  );
  double u_upper_val = -Variant_Grid_access(u_upper, x  , y  , z-1);
  Variant_Grid_access(fp, x, y, z) += u_right_val + u_front_val + u_upper_val;
}

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

  // Create domain for reduction portion. It is 1 larger than the normal domain.
  Variant_Domain* reduction_domain = Variant_Domain_alloc( Variant_Domain_nx(domain)+1, Variant_Domain_nz(domain)+1, Variant_Domain_ny(domain)+1 );

  // Create temporary grids for reduction portion.
  Variant_Grid* u_right = Variant_Grid_alloc( domain );
  Variant_Grid* u_front = Variant_Grid_alloc( domain );
  Variant_Grid* u_upper = Variant_Grid_alloc( domain );
  
/* ------------------------------ NLFE216 ------------------------------ */

  //Variables used for timing
  double NLFE216_CUDAAllocTime;
  double NLFE216_CUDAHTDTime;
  double NLFE216_CUDAKernelTime;

  //Determine dimensions of kernel grid
  //A static block size of 16x16x4 threads is used
  int blockx = 16;
  int blocky = 16;
  int blockz = 4;
  int gridx = (int) ceil(((float) (Basic_Domain_nx(domain) - 2)) / blockx);
  int gridy = (int) ceil(((float) (Basic_Domain_ny(domain) - 2)) / blocky);
  int gridz = (int) ceil(((float) (Basic_Domain_nz(domain) - 2)) / blockz);
  dim3 block = dim3(blockx, blocky, blockz);
  dim3 grid = dim3(gridx, gridy, gridz);

  //Device pointers for the grids
  Variant_Grid *fpCUDA, *spCUDA, *dpCUDA, *ospCUDA, *odpCUDA, *popCUDA, *z_mult_datCUDA;

  //Time the allocation time for NLFE216
  TIMEIT(NLFE216_CUDAAllocTime, {
    assert(cudaMalloc(&fpCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&spCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&dpCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&ospCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&odpCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&popCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&z_mult_datCUDA, sizeof(Variant_Grid)) == cudaSuccess);
  })

  //Time the Host to Device copy time
  TIMEIT(NLFE216_CUDAHTDTime, {
    assert(cudaMemcpy(fpCUDA, fp, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(spCUDA, sp, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(dpCUDA, dp, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(ospCUDA, osp, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(odpCUDA, odp, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(popCUDA, pop, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(z_mult_datCUDA, z_mult_dat, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
  })

  //Time the kernel execution time
  TIMEIT(NLFE216_CUDAKernelTime, {
    NLFE216Kernel<<<grid, block>>>(Basic_Domain_nx(domain) - 2, Basic_Domain_ny(domain) - 2, Basic_Domain_nz(domain) -2, fpCUDA, spCUDA, dpCUDA, ospCUDA, odpCUDA, popCUDA, z_mult_datCUDA);
    cudaDeviceSynchronize();
  })

/* ------------------------------ NLFE338 ------------------------------ */

  //Variables used for timing
  double NLFE338_CUDAAllocTime;
  double NLFE338_CUDAHTDTime;
  double NLFE338_CUDAKernelTime;

  //Device pointers for the grids not already on the device
  Variant_Grid *ssCUDA, *ppCUDA, *oppCUDA;

  //Time the allocation time for NLFE338
  TIMEIT(NLFE338_CUDAAllocTime, {
    assert(cudaMalloc(&ssCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&ppCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&oppCUDA, sizeof(Variant_Grid)) == cudaSuccess);
  })

  //Time the Host to Device copy time
  TIMEIT(NLFE338_CUDAHTDTime, {
    assert(cudaMemcpy(ssCUDA, ss, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(ppCUDA, pp, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(oppCUDA, opp, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
  })

  //Time the kernel execution time
  TIMEIT(NLFE338_CUDAKernelTime, {
    NLFE338Kernel<<<grid, block>>>(Basic_Domain_nx(domain) - 2, Basic_Domain_ny(domain) - 2, Basic_Domain_nz(domain) -2, fpCUDA, ssCUDA, z_mult_datCUDA, ppCuda, spCUDA, dpCUDA, oppCUDA, ospCUDA, odpCUDA);
    cudaDeviceSynchronize();
  })

/* ------------------------------ NLFE416 ------------------------------ */

  //Variables used for timing
  double NLFE416_CUDAAllocTime;
  double NLFE416_CUDAHTDTime;
  double NLFE416_CUDAKernelTime;

  //Device pointers for the grids not already on the device
  Variant_Grid *etCUDA;

  //Time the allocation time for NLFE416
  TIMEIT(NLFE416_CUDAAllocTime, {
    assert(cudaMalloc(&etCUDA, sizeof(Variant_Grid)) == cudaSuccess);
  })

  //Time the Host to Device copy time
  TIMEIT(NLFE416_CUDAHTDTime, {
    assert(cudaMemcpy(etCUDA, et, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
  })

  //Time the kernel execution time
  TIMEIT(NLFE416_CUDAKernelTime, {
    NLFE416Kernel<<<grid, block>>>(Basic_Domain_nx(domain) - 2, Basic_Domain_ny(domain) - 2, Basic_Domain_nz(domain) -2, fpCUDA, z_mult_datCUDA, spCUDA, etCUDA);
    cudaDeviceSynchronize();
  })

/* ------------------------------ NLFE551 ------------------------------ */

  //Variables used for timing
  double NLFE551_CUDAAllocTime;
  double NLFE551_CUDAHTDTime;
  double NLFE551_CUDAKernelTime;

  //Device pointers for the grids not already on the device
  Variant_Grid *x_ssl_datCUDA, *y_ssl_datCUDA, *permxpCUDA, *rppCUDA, *permypCUDA, *permzpCUDA, *vxCUDA, *vyCUDA, *vzCUDA, *u_rightCUDA, *u_frontCUDA, *u_upperCUDA;

  //Time the allocation time for NLFE551
  TIMEIT(NLFE551_CUDAAllocTime, {
    assert(cudaMalloc(&x_ssl_datCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&y_ssl_datCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&permxpCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&rppCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&permypCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&permzpCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&vxCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&vyCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&vzCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&u_rightCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&u_frontCUDA, sizeof(Variant_Grid)) == cudaSuccess);
    assert(cudaMalloc(&u_upperCUDA, sizeof(Variant_Grid)) == cudaSuccess);
  })

  //Time the Host to Device copy time
  TIMEIT(NLFE551_CUDAHTDTime, {
    assert(cudaMemcpy(x_ssl_datCUDA, x_ssl_dat, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(y_ssl_datCUDA, y_ssl_dat, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(permxpCUDA, permxp, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(rppCUDA, rpp, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(permypCUDA, permyp, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(permzpCUDA, permzp, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(vxCUDA, vx, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(vyCUDA, vy, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(vzCUDA, vz, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(u_rightCUDA, u_right, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(u_frontCUDA, u_front, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
    assert(cudaMemcpy(u_upperCUDA, u_upper, sizeof(Variant_Grid), cudaMemcpyHostToDevice) == cudaSuccess);
  })
  
  //Time the kernel execution time
  TIMEIT(NLFE551_CUDAKernelTime, {
    NLFE551Kernel<<<grid, block>>>(Basic_Domain_nx(domain) - 2, Basic_Domain_ny(domain) - 2, Basic_Domain_nz(domain) -2, x_ssl_datCUDA, y_ssl_datCUDA, ppCUDA, z_mult_datCUDA, dpCUDA, permxpCUDA, rppCUDA, permypCUDA, permzpCUDA, vxCUDA, vyCUDA, vzCUDA, u_rightCUDA, u_frontCUDA, u_upperCUDA, fpCUDA);
    cudaDeviceSynchronize();
  })

/* ------------------------------ NLFE551 Reduction ------------------------------ */

  //Variables used for timing
  double NLFE551Reduction_CUDAKernelTime;

  //Time the kernel execution time
  TIMEIT(NLFE551_CUDAKernelTime, {
    NLFE551ReductionKernel<<<grid, block>>>(Basic_Domain_nx(reduction_domain) - 2, Basic_Domain_ny(reduction_domain) - 2, Basic_Domain_nz(reduction_domain) -2, u_rightCUDA, u_frontCUDA, u_upperCUDA, fpCUDA);
    cudaDeviceSynchronize();
  })

/* ------------------------------ Cleanup ------------------------------ */

  //Variables used for timing
  double Cleanup_CUDADTHTime;
  double Cleanup_CUDAFreeTime;

  //Time the Device to Host transfer time for the final grids
  TIMEIT(Cleanup_CUDADTHTime, {
    assert(cudaMemcpy(vx, vxCUDA, sizeof(Variant_Grid), cudaMemcpyDeviceToHost) == cudaSuccess);
    assert(cudaMemcpy(vy, vyCUDA, sizeof(Variant_Grid), cudaMemcpyDeviceToHost) == cudaSuccess);
    assert(cudaMemcpy(vz, vzCUDA, sizeof(Variant_Grid), cudaMemcpyDeviceToHost) == cudaSuccess);
    assert(cudaMemcpy(fp, fpCUDA, sizeof(Variant_Grid), cudaMemcpyDeviceToHost) == cudaSuccess);
  })

  //Time the freeing of all device grids
  TIMEIT(Cleanup_CUDAFreeTime, {
    assert(cudaFree(fpCUDA) == cudaSuccess);
    assert(cudaFree(vxCUDA) == cudaSuccess);
    assert(cudaFree(vyCUDA) == cudaSuccess);
    assert(cudaFree(vzCUDA) == cudaSuccess);
    assert(cudaFree(dpCUDA) == cudaSuccess);
    assert(cudaFree(etCUDA) == cudaSuccess);
    assert(cudaFree(odpCUDA) == cudaSuccess);
    assert(cudaFree(oppCUDA) == cudaSuccess);
    assert(cudaFree(ospCUDA) == cudaSuccess);
    assert(cudaFree(permxpCUDA) == cudaSuccess);
    assert(cudaFree(permypCUDA) == cudaSuccess);
    assert(cudaFree(permzpCUDA) == cudaSuccess);
    assert(cudaFree(popCUDA) == cudaSuccess);
    assert(cudaFree(ppCUDA) == cudaSuccess);
    assert(cudaFree(rppCUDA) == cudaSuccess);
    assert(cudaFree(spCUDA) == cudaSuccess);
    assert(cudaFree(ssCUDA) == cudaSuccess);
    assert(cudaFree(z_mult_datCUDA) == cudaSuccess);
    assert(cudaFree(x_ssl_datCUDA) == cudaSuccess);
    assert(cudaFree(y_ssl_datCUDA) == cudaSuccess);
  })

  //Deallocate local grids and domains
  Variant_Domain_dealloc( reduction_domain );
  Variant_Grid_dealloc(u_right);
  Variant_Grid_dealloc(u_front);
  Variant_Grid_dealloc(u_upper);
}