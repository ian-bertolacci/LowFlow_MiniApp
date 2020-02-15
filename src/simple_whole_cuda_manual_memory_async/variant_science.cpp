#include <science.hpp>
#include <macros.hpp>
#include <configure.hpp>
#include <util.hpp>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

size_t gridSize, domainSize, dataSize;

inline void checkCUDAError(cudaError_t code, const char *file, int line) {
   if (code != cudaSuccess) {
      fprintf(stderr,"\nCUDA Error: %s\nFile: %s\nLine: %d\n", cudaGetErrorString(code), file, line);
      exit(code);
   }
}

#define check(expr) { \
  cudaError_t __e__ = expr; \
  if (ENABLE_DEBUG) { checkCUDAError(__e__, __FILE__, __LINE__); } \
}

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

void prepareDeviceGrid(Variant_Grid *host, Variant_Grid **device, cudaStream_t stream) {
  Variant_Grid temp;

  check(cudaMalloc(&temp.domain, domainSize));
  check(cudaMalloc(&temp.data, dataSize));

  check(cudaMemcpyAsync(temp.domain, host->domain, domainSize, cudaMemcpyHostToDevice, stream));
  check(cudaMemcpyAsync(temp.data, host->data, dataSize, cudaMemcpyHostToDevice, stream));

  check(cudaMalloc(device, gridSize));
  check(cudaMemcpyAsync(*device, &temp, gridSize, cudaMemcpyHostToDevice, stream));
}

void copyDeviceToHost(Variant_Grid *host, Variant_Grid **device) {
  Variant_Grid temp;

  check(cudaMemcpy(&temp, *device, gridSize, cudaMemcpyDeviceToHost));
  check(cudaMemcpy(host->data, temp.data, dataSize, cudaMemcpyDeviceToHost));

  check(cudaFree(temp.domain));
  check(cudaFree(temp.data));
  check(cudaFree(*device));
}

void freeDeviceGrid(Variant_Grid **device) {
  Variant_Grid temp;

  check(cudaMemcpy(&temp, *device, gridSize, cudaMemcpyDeviceToHost));

  check(cudaFree(temp.domain));
  check(cudaFree(temp.data));
  check(cudaFree(*device));
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
  Variant_Grid *u_right, *u_front, *u_upper;
  u_right = Variant_Grid_alloc( domain );
  u_front = Variant_Grid_alloc( domain );
  u_upper = Variant_Grid_alloc( domain );

  //Determine dimensions of kernel grid
  //A static block size of 16x16x4 threads is used
  int blockx = 16;
  int blocky = 16;
  int blockz = 4;
  int gridx = (int) ceil((float) Basic_Domain_nx(domain) / blockx);
  int gridy = (int) ceil((float) Basic_Domain_ny(domain) / blocky);
  int gridz = (int) ceil((float) Basic_Domain_nz(domain) / blockz);
  dim3 block = dim3(blockx, blocky, blockz);
  dim3 grid = dim3(gridx, gridy, gridz);

  //Get the size info of the grids
  gridSize = sizeof(Variant_Grid);
  domainSize = sizeof(Variant_Domain);
  dataSize = Basic_Domain_nx(domain) * Basic_Domain_ny(domain) * Basic_Domain_nz(domain) * sizeof(double);

  //Create the CUDA streams needed for async execution
  cudaStream_t s1, s2, s3, s4, s5;
  cudaStreamCreate(&s1);
  cudaStreamCreate(&s2);
  cudaStreamCreate(&s3);
  cudaStreamCreate(&s4);
  cudaStreamCreate(&s5);

  //Create grids on the device for this kernel call
  Variant_Grid *fpCUDA, *spCUDA, *dpCUDA, *ospCUDA, *odpCUDA, *popCUDA, *z_mult_datCUDA;
  Variant_Grid *ssCUDA, *ppCUDA, *oppCUDA;
  Variant_Grid *etCUDA;
  Variant_Grid *x_ssl_datCUDA, *y_ssl_datCUDA, *permxpCUDA, *rppCUDA, *permypCUDA, *permzpCUDA, *vxCUDA, *vyCUDA, *vzCUDA, *u_rightCUDA, *u_frontCUDA, *u_upperCUDA;

  //s1 grids
  prepareDeviceGrid(fp, &fpCUDA, s1);
  prepareDeviceGrid(sp, &spCUDA, s1);
  prepareDeviceGrid(dp, &dpCUDA, s1);
  prepareDeviceGrid(osp, &ospCUDA, s1);
  prepareDeviceGrid(odp, &odpCUDA, s1);
  prepareDeviceGrid(pop, &popCUDA, s1);
  prepareDeviceGrid(z_mult_dat, &z_mult_datCUDA, s1);

  //s2 grids
  prepareDeviceGrid(ss, &ssCUDA, s2);
  prepareDeviceGrid(pp, &ppCUDA, s2);
  prepareDeviceGrid(opp, &oppCUDA, s2);

  //s3 grids
  prepareDeviceGrid(et, &etCUDA, s3);

  //s4 grids
  prepareDeviceGrid(x_ssl_dat, &x_ssl_datCUDA, s4);
  prepareDeviceGrid(y_ssl_dat, &y_ssl_datCUDA, s4);
  prepareDeviceGrid(permxp, &permxpCUDA, s4);
  prepareDeviceGrid(rpp, &rppCUDA, s4);
  prepareDeviceGrid(permyp, &permypCUDA, s4);
  prepareDeviceGrid(permzp, &permzpCUDA, s4);
  prepareDeviceGrid(vx, &vxCUDA, s4);
  prepareDeviceGrid(vy, &vyCUDA, s4);
  prepareDeviceGrid(vz, &vzCUDA, s4);
  prepareDeviceGrid(u_right, &u_rightCUDA, s4);
  prepareDeviceGrid(u_front, &u_frontCUDA, s4);
  prepareDeviceGrid(u_upper, &u_upperCUDA, s4);

  //s1 kernel
  NLFE216Kernel<<<grid, block, 0, s1>>>(Basic_Domain_nx(domain) - 2, Basic_Domain_ny(domain) - 2, Basic_Domain_nz(domain) -2, fpCUDA, spCUDA, dpCUDA, ospCUDA, odpCUDA, popCUDA, z_mult_datCUDA);
  check(cudaStreamSynchronize(s1));

  //s2 kernel
  NLFE338Kernel<<<grid, block, 0, s2>>>(Basic_Domain_nx(domain) - 2, Basic_Domain_ny(domain) - 2, Basic_Domain_nz(domain) -2, fpCUDA, ssCUDA, z_mult_datCUDA, ppCUDA, spCUDA, dpCUDA, oppCUDA, ospCUDA, odpCUDA);
  check(cudaStreamSynchronize(s2));

  //s3 kernel
  NLFE416Kernel<<<grid, block, 0, s3>>>(Basic_Domain_nx(domain) - 2, Basic_Domain_ny(domain) - 2, Basic_Domain_nz(domain) -2, fpCUDA, z_mult_datCUDA, spCUDA, etCUDA);
  check(cudaStreamSynchronize(s3));

  //s4 kernel
  NLFE551Kernel<<<grid, block, 0, s4>>>(Basic_Domain_nx(domain) - 2, Basic_Domain_ny(domain) - 2, Basic_Domain_nz(domain) -2, x_ssl_datCUDA, y_ssl_datCUDA, ppCUDA, z_mult_datCUDA, dpCUDA, permxpCUDA, rppCUDA, permypCUDA, permzpCUDA, vxCUDA, vyCUDA, vzCUDA, u_rightCUDA, u_frontCUDA, u_upperCUDA, fpCUDA);
  check(cudaStreamSynchronize(s4));

  //s5 kernel
  NLFE551ReductionKernel<<<grid, block, 0, s5>>>(Basic_Domain_nx(reduction_domain) - 2, Basic_Domain_ny(reduction_domain) - 2, Basic_Domain_nz(reduction_domain) -2, u_rightCUDA, u_frontCUDA, u_upperCUDA, fpCUDA);
  check(cudaStreamSynchronize(s5));
  
  //Destroy streams
  check(cudaStreamDestroy(s1));
  check(cudaStreamDestroy(s2));
  check(cudaStreamDestroy(s3));
  check(cudaStreamDestroy(s4));
  check(cudaStreamDestroy(s5));

  //Copy over result grids
  copyDeviceToHost(fp, &fpCUDA);
  copyDeviceToHost(vx, &vxCUDA);
  copyDeviceToHost(vy, &vyCUDA);
  copyDeviceToHost(vz, &vzCUDA);

  //Free all other device grids and their members
  freeDeviceGrid(&dpCUDA);
  freeDeviceGrid(&etCUDA);
  freeDeviceGrid(&odpCUDA);
  freeDeviceGrid(&oppCUDA);
  freeDeviceGrid(&ospCUDA);
  freeDeviceGrid(&permxpCUDA);
  freeDeviceGrid(&permypCUDA);
  freeDeviceGrid(&permzpCUDA);
  freeDeviceGrid(&popCUDA);
  freeDeviceGrid(&ppCUDA);
  freeDeviceGrid(&rppCUDA);
  freeDeviceGrid(&spCUDA);
  freeDeviceGrid(&ssCUDA);
  freeDeviceGrid(&z_mult_datCUDA);
  freeDeviceGrid(&x_ssl_datCUDA);
  freeDeviceGrid(&y_ssl_datCUDA);
  freeDeviceGrid(&u_rightCUDA);
  freeDeviceGrid(&u_frontCUDA);
  freeDeviceGrid(&u_upperCUDA);

  //Cleanup temp domain and grids
  Variant_Domain_dealloc( reduction_domain );
  Variant_Grid_dealloc(u_right);
  Variant_Grid_dealloc(u_front);
  Variant_Grid_dealloc(u_upper);
}