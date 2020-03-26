#include <science.hpp>
#include <macros.hpp>
#include <configure.hpp>
#include <util.hpp>
#include <math.h>
#include <assert.h>
#include <stdio.h>

#define CHUNKS 10
#define STREAMS 5

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
__global__ void NLFE216Kernel(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, Variant_Grid *fp, Variant_Grid *sp, Variant_Grid *dp, Variant_Grid *osp, Variant_Grid *odp, Variant_Grid *pop, Variant_Grid *z_mult_dat) {
  //Get position in kernel grid
  //+1 is added to coordinates to account for offset in each direction
  int x = blockIdx.x * blockDim.x + threadIdx.x + xmin;
  int y = blockIdx.y * blockDim.y + threadIdx.y + ymin;
  int z = blockIdx.z * blockDim.z + threadIdx.z + zmin;
  
  //Make sure that the thread is in bounds
  if (x > xmax || y > ymax || z > zmax) {
    return;
  }

  //Perform variant grid access
  Variant_Grid_access(fp,x,y,z) = (Variant_Grid_access(sp,x,y,z) * Variant_Grid_access(dp,x,y,z) - Variant_Grid_access(osp,x,y,z) * Variant_Grid_access(odp, x,y,z)) * Variant_Grid_access(pop, x,y,z) * Variant_Grid_access(z_mult_dat, x,y,z);
}

//NlFunctionEval:338 CUDA kernel 
__global__ void NLFE338Kernel(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, Variant_Grid *fp, Variant_Grid *ss, Variant_Grid *z_mult_dat, Variant_Grid *pp, Variant_Grid *sp, Variant_Grid *dp, Variant_Grid *opp, Variant_Grid *osp, Variant_Grid *odp) {
  //Get position in kernel grid
  //+1 is added to coordinates to account for offset in each direction
  int x = blockIdx.x * blockDim.x + threadIdx.x + xmin;
  int y = blockIdx.y * blockDim.y + threadIdx.y + ymin;
  int z = blockIdx.z * blockDim.z + threadIdx.z + zmin;

  //Make sure that the thread is in bounds
  if (x > xmax || y > ymax || z > zmax) {
    return;
  }

  //Perform variant grid access
  Variant_Grid_access(fp, x,y,z) += Variant_Grid_access(ss, x,y,z) * Variant_Grid_access(z_mult_dat, x,y,z) * (Variant_Grid_access(pp, x,y,z) * Variant_Grid_access(sp, x,y,z) * Variant_Grid_access(dp, x,y,z) - Variant_Grid_access(opp, x,y,z) * Variant_Grid_access(osp, x,y,z) * Variant_Grid_access(odp, x,y,z));
}

//NlFunctionEval:416 CUDA kernel 
__global__ void NLFE416Kernel(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, Variant_Grid *fp, Variant_Grid *z_mult_dat, Variant_Grid *sp, Variant_Grid *et) {
  //Get position in kernel grid
  //+1 is added to coordinates to account for offset in each direction
  int x = blockIdx.x * blockDim.x + threadIdx.x + xmin;
  int y = blockIdx.y * blockDim.y + threadIdx.y + ymin;
  int z = blockIdx.z * blockDim.z + threadIdx.z + zmin;

  //Make sure that the thread is in bounds
  if (x > xmax || y > ymax || z > zmax) {
    return;
  }

  //Perform variant grid access
  Variant_Grid_access(fp, x,y,z) -= Variant_Grid_access(z_mult_dat, x,y,z) * (Variant_Grid_access(sp, x,y,z) * Variant_Grid_access(et, x,y,z));
}

//NlFunctionEval:551 CUDA kernel 
__global__ void NLFE551Kernel(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, Variant_Grid *x_ssl_dat, Variant_Grid *y_ssl_dat, Variant_Grid *pp, Variant_Grid *z_mult_dat, Variant_Grid *dp, Variant_Grid *permxp, Variant_Grid *rpp, Variant_Grid *permyp, Variant_Grid *permzp, Variant_Grid *vx, Variant_Grid *vy, Variant_Grid *vz, Variant_Grid *u_right, Variant_Grid *u_front, Variant_Grid *u_upper, Variant_Grid *fp) {
  //Get position in kernel grid
  //+1 is added to coordinates to account for offset in each direction
  int x = blockIdx.x * blockDim.x + threadIdx.x + xmin;
  int y = blockIdx.y * blockDim.y + threadIdx.y + ymin;
  int z = blockIdx.z * blockDim.z + threadIdx.z + zmin;

  //Make sure that the thread is in bounds
  if (x > xmax || y > ymax || z > zmax) {
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

__global__ void NLFE551ReductionKernel(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, Variant_Grid *u_right, Variant_Grid *u_front, Variant_Grid *u_upper, Variant_Grid *fp) {
  //Get position in kernel grid
  //+1 is added to coordinates to account for offset in each direction
  int x = blockIdx.x * blockDim.x + threadIdx.x + xmin;
  int y = blockIdx.y * blockDim.y + threadIdx.y + ymin;
  int z = blockIdx.z * blockDim.z + threadIdx.z + zmin;

  //Make sure that the thread is in bounds
  if (x > xmax || y > ymax || z > zmax) {
    return;
  }

  double u_right_val =  Variant_Grid_access(u_right, x-1, y  , z  );
  double u_front_val = -Variant_Grid_access(u_front, x  , y-1, z  );
  double u_upper_val = -Variant_Grid_access(u_upper, x  , y  , z-1);
  Variant_Grid_access(fp, x, y, z) += u_right_val + u_front_val + u_upper_val;
}

void allocateDeviceGrid(Variant_Grid *host, Variant_Grid **device) {
  Variant_Grid temp;

  check(cudaMalloc(&temp.domain, domainSize));
  check(cudaMalloc(&temp.data, dataSize));

  check(cudaMemcpy(temp.domain, host->domain, domainSize, cudaMemcpyHostToDevice));

  check(cudaMalloc(device, gridSize));
  check(cudaMemcpy(*device, &temp, gridSize, cudaMemcpyHostToDevice));
}

void freeDeviceGrid(Variant_Grid **device) {
  Variant_Grid temp;

  check(cudaMemcpy(&temp, *device, gridSize, cudaMemcpyDeviceToHost));

  check(cudaFree(temp.domain));
  check(cudaFree(temp.data));
  check(cudaFree(*device));
}

void copyChunk(Variant_Grid *host, Variant_Grid *device, int memoryOffset, int copySize, cudaStream_t stream) {
  Variant_Grid temp;

  check(cudaMemcpyAsync(&temp, device, gridSize, cudaMemcpyDeviceToHost, stream));

  check(cudaMemcpyAsync(temp.data + memoryOffset, host->data + memoryOffset, copySize, cudaMemcpyHostToDevice, stream));
}

void copyDeviceToHost(Variant_Grid *host, Variant_Grid *device) {
  Variant_Grid temp;

  check(cudaMemcpy(&temp, device, gridSize, cudaMemcpyDeviceToHost));
  check(cudaMemcpy(host->data, temp.data, dataSize, cudaMemcpyDeviceToHost));
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

  //Get the size info of the grids
  gridSize = sizeof(Variant_Grid);
  domainSize = sizeof(Variant_Domain);
  dataSize = Basic_Domain_nx(domain) * Basic_Domain_ny(domain) * Basic_Domain_nz(domain) * sizeof(double);

  /*
  * Basic idea, split thread grid into chunks, only copy over/execute for that chunk per iteration
  * Everything is allocated in the beginning, but memcpys happen asynchronously during execution
  * The reduction kernel can be done last
  */
  
  //These grids will be allocated on the device
  Variant_Grid *fpCUDA, *spCUDA, *dpCUDA, *ospCUDA, *odpCUDA, *popCUDA, *z_mult_datCUDA;
  Variant_Grid *ssCUDA, *ppCUDA, *oppCUDA;
  Variant_Grid *etCUDA;
  Variant_Grid *x_ssl_datCUDA, *y_ssl_datCUDA, *permxpCUDA, *rppCUDA, *permypCUDA, *permzpCUDA, *vxCUDA, *vyCUDA, *vzCUDA, *u_rightCUDA, *u_frontCUDA, *u_upperCUDA;

  //Allocate each of the grids listed above on the device
  allocateDeviceGrid(fp, &fpCUDA);
  allocateDeviceGrid(sp, &spCUDA);
  allocateDeviceGrid(dp, &dpCUDA);
  allocateDeviceGrid(osp, &ospCUDA);
  allocateDeviceGrid(odp, &odpCUDA);
  allocateDeviceGrid(pop, &popCUDA);
  allocateDeviceGrid(z_mult_dat, &z_mult_datCUDA);
  allocateDeviceGrid(ss, &ssCUDA);
  allocateDeviceGrid(pp, &ppCUDA);
  allocateDeviceGrid(opp, &oppCUDA);
  allocateDeviceGrid(et, &etCUDA);
  allocateDeviceGrid(x_ssl_dat, &x_ssl_datCUDA);
  allocateDeviceGrid(y_ssl_dat, &y_ssl_datCUDA);
  allocateDeviceGrid(permxp, &permxpCUDA);
  allocateDeviceGrid(rpp, &rppCUDA);
  allocateDeviceGrid(permyp, &permypCUDA);
  allocateDeviceGrid(permzp, &permzpCUDA);
  allocateDeviceGrid(vx, &vxCUDA);
  allocateDeviceGrid(vy, &vyCUDA);
  allocateDeviceGrid(vz, &vzCUDA);
  allocateDeviceGrid(u_right, &u_rightCUDA);
  allocateDeviceGrid(u_front, &u_frontCUDA);
  allocateDeviceGrid(u_upper, &u_upperCUDA);

  //Determine dimensions of kernel grid
  int blockx = 32;
  int blocky = 32;
  int blockz = 1;
  int gridx = (int) ceil((float) Basic_Domain_nx(domain) / blockx);
  int gridy = (int) ceil((float) Basic_Domain_ny(domain) / blocky);
  int gridz = (int) ceil((float) Basic_Domain_nz(domain) / CHUNKS);
  dim3 block = dim3(blockx, blocky, blockz);
  dim3 grid = dim3(gridx, gridy, gridz);

  //Create the streams
  cudaStream_t streams[STREAMS];  
  for (int i = 0; i < STREAMS; i++) {
    check(cudaStreamCreate(&streams[i]));
  }

  //Loop once for each chunk, running all kernels besides the reduction kernel
  int memoryOffset, copySize, zCopyAmount;
  int faceSize = Basic_Domain_nx(domain) * Basic_Domain_ny(domain);
  int faceMemorySize = faceSize * sizeof(double);
  int x1 = 1;
  int x2 = Basic_Domain_nx(domain) - 2;
  int y1 = 1;
  int y2 = Basic_Domain_ny(domain) - 2;
  int z1, z2;
  cudaStream_t stream;

  for (int i = 0; i < CHUNKS; i++) {
    //Get values for offsets and memcpy size
    //memoryOffset is how many doubles forward to jump in array
    //copySize is size in bytes to copy over, which could be a whole stride or whatever is left to compute
    //z1/z2 are the lower and upper bounds (inclusive) of the kernel in the z dimension, z2 - z1 should be <= gridz
    z1 = i * gridz + 1;
    z2 = min(z1 + gridz - 1, Basic_Domain_nz(domain) - 2);
    if (z1 > z2) {
      break;
    }
    if (i == 0) {
      zCopyAmount = gridz + 1;                                          //Also copy over z = 0
    }
    else if (i == CHUNKS - 1) {
      zCopyAmount = max(Basic_Domain_nz(domain) - z1, 0);          //Copy over whatever is left
    }
    else {
      zCopyAmount = gridz;                                              //We are somewhere in the middle, copy current section only
    }
    memoryOffset = faceSize * (i == 0 ? 0 : z1);
    copySize = faceMemorySize * zCopyAmount;
    stream = streams[i % STREAMS];

    printf("Chunk %d Info: {z-range: [%d, %d], offset: %d, size: %d}\n", i, z1, z2, memoryOffset, copySize);

    //Loop 1 memcpys and exec
    copyChunk(fp, fpCUDA, memoryOffset, copySize, stream);
    copyChunk(sp, spCUDA, memoryOffset, copySize, stream);
    copyChunk(dp, dpCUDA, memoryOffset, copySize, stream);
    copyChunk(osp, ospCUDA, memoryOffset, copySize, stream);
    copyChunk(odp, odpCUDA, memoryOffset, copySize, stream);
    copyChunk(pop, popCUDA, memoryOffset, copySize, stream);
    copyChunk(z_mult_dat, z_mult_datCUDA, memoryOffset, copySize, stream);
    NLFE216Kernel<<<grid, block, 0, stream>>>(x1, x2, y1, y2, z1, z2, fpCUDA, spCUDA, dpCUDA, ospCUDA, odpCUDA, popCUDA, z_mult_datCUDA);

    //Loop 2 memcpys and exec
    copyChunk(ss, ssCUDA, memoryOffset, copySize, stream);
    copyChunk(pp, ppCUDA, memoryOffset, copySize, stream);
    copyChunk(opp, oppCUDA, memoryOffset, copySize, stream);
    NLFE338Kernel<<<grid, block, 0, stream>>>(x1, x2, y1, y2, z1, z2, fpCUDA, ssCUDA, z_mult_datCUDA, ppCUDA, spCUDA, dpCUDA, oppCUDA, ospCUDA, odpCUDA);

    //Loop 3 memcpys and exec
    copyChunk(et, etCUDA, memoryOffset, copySize, stream);
    NLFE416Kernel<<<grid, block, 0, stream>>>(x1, x2, y1, y2, z1, z2, fpCUDA, z_mult_datCUDA, spCUDA, etCUDA);

    //Loop 4 memcpys and exec
    copyChunk(x_ssl_dat, x_ssl_datCUDA, memoryOffset, copySize, stream);
    copyChunk(y_ssl_dat, y_ssl_datCUDA, memoryOffset, copySize, stream);
    copyChunk(permxp, permxpCUDA, memoryOffset, copySize, stream);
    copyChunk(rpp, rppCUDA, memoryOffset, copySize, stream);
    copyChunk(permyp, permypCUDA, memoryOffset, copySize, stream);
    copyChunk(permzp, permzpCUDA, memoryOffset, copySize, stream);
    copyChunk(vx, vxCUDA, memoryOffset, copySize, stream);
    copyChunk(vy, vyCUDA, memoryOffset, copySize, stream);
    copyChunk(vz, vzCUDA, memoryOffset, copySize, stream);
    copyChunk(u_right, u_rightCUDA, memoryOffset, copySize, stream);
    copyChunk(u_front, u_frontCUDA, memoryOffset, copySize, stream);
    copyChunk(u_upper, u_upperCUDA, memoryOffset, copySize, stream);
    NLFE551Kernel<<<grid, block, 0, stream>>>(x1, x2, y1, y2, z1, z2, x_ssl_datCUDA, y_ssl_datCUDA, ppCUDA, z_mult_datCUDA, dpCUDA, permxpCUDA, rppCUDA, permypCUDA, permzpCUDA, vxCUDA, vyCUDA, vzCUDA, u_rightCUDA, u_frontCUDA, u_upperCUDA, fpCUDA);
  }

  //Wait for all of the above loop executions to finish before proceeding
  for (int i = 0; i < STREAMS; i++) {
    check(cudaStreamSynchronize(streams[i]));
  }

  //Run the reduction kernel, splitting similarly to above
  //All data has already been copied over, no need to copy chunks, but can split into chunks still
  gridx = (int) ceil((float) Basic_Domain_nx(reduction_domain) / blockx);
  gridy = (int) ceil((float) Basic_Domain_ny(reduction_domain) / blocky);
  gridz = (int) ceil((float) Basic_Domain_nz(reduction_domain) / CHUNKS);
  grid = dim3(gridx, gridy, gridz);
  x2 = Basic_Domain_nx(reduction_domain) - 2;
  y2 = Basic_Domain_ny(reduction_domain) - 2;
  for (int i = 0; i < CHUNKS; i++) {
    //Get the z bounds
    z1 = i * gridz + 1;
    z2 = min(z1 + gridz - 1, Basic_Domain_nz(reduction_domain) - 2);
    if (z1 > z2) {
      break;
    }
    printf("Chunk %d Info: {z-range: [%d, %d]}\n", i, z1, z2);
    //Execute the kernel on a different stream than the last
    NLFE551ReductionKernel<<<grid, block, 0, streams[i % STREAMS]>>>(x1, x2, y1, y2, z1, z2, u_rightCUDA, u_frontCUDA, u_upperCUDA, fpCUDA);
  }

  //Wait for all of the above loop executions to finish before proceeding
  for (int i = 0; i < STREAMS; i++) {
    check(cudaStreamSynchronize(streams[i]));
  }

  //Copy the output grids back to the host
  copyDeviceToHost(fp, fpCUDA);
  copyDeviceToHost(vx, vxCUDA);
  copyDeviceToHost(vy, vyCUDA);
  copyDeviceToHost(vz, vzCUDA);
  
  //Destroy the streams
  for (int i = 0; i < STREAMS; i++) {
    check(cudaStreamDestroy(streams[i]));
  }
  
  //Free the device grids
  freeDeviceGrid(&fpCUDA);
  freeDeviceGrid(&vxCUDA);
  freeDeviceGrid(&vyCUDA);
  freeDeviceGrid(&vzCUDA);
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