#include <science.hpp>
#include <macros.hpp>
#include <configure.hpp>
#include <util.hpp>
#include <math.h>
#include <assert.h>
#include <stdio.h>

size_t domainSize, dataSize;

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
__global__ void NLFE216Kernel(int xmax, int ymax, int zmax, Variant_Domain *fp_Domain, double *fp, Variant_Domain *sp_Domain, double *sp, Variant_Domain *dp_Domain, double *dp, Variant_Domain *osp_Domain, double *osp, Variant_Domain *odp_Domain, double *odp, Variant_Domain *pop_Domain, double *pop, Variant_Domain *z_mult_dat_Domain, double *z_mult_dat) {
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
  Variant_Grid_access_extracted(fp,x,y,z) = (Variant_Grid_access_extracted(sp,x,y,z) * Variant_Grid_access_extracted(dp,x,y,z) - Variant_Grid_access_extracted(osp,x,y,z) * Variant_Grid_access_extracted(odp, x,y,z)) * Variant_Grid_access_extracted(pop, x,y,z) * Variant_Grid_access_extracted(z_mult_dat, x,y,z);
}

//NlFunctionEval:338 CUDA kernel 
__global__ void NLFE338Kernel(int xmax, int ymax, int zmax, Variant_Domain *fp_Domain, double *fp, Variant_Domain *ss_Domain, double *ss, Variant_Domain *z_mult_dat_Domain, double *z_mult_dat, Variant_Domain *pp_Domain, double *pp, Variant_Domain *sp_Domain, double *sp, Variant_Domain *dp_Domain, double *dp, Variant_Domain *opp_Domain, double *opp, Variant_Domain *osp_Domain, double *osp, Variant_Domain *odp_Domain, double *odp) {
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
  Variant_Grid_access_extracted(fp, x,y,z) += Variant_Grid_access_extracted(ss, x,y,z) * Variant_Grid_access_extracted(z_mult_dat, x,y,z) * (Variant_Grid_access_extracted(pp, x,y,z) * Variant_Grid_access_extracted(sp, x,y,z) * Variant_Grid_access_extracted(dp, x,y,z) - Variant_Grid_access_extracted(opp, x,y,z) * Variant_Grid_access_extracted(osp, x,y,z) * Variant_Grid_access_extracted(odp, x,y,z));
}

//NlFunctionEval:416 CUDA kernel 
__global__ void NLFE416Kernel(int xmax, int ymax, int zmax, Variant_Domain *fp_Domain, double *fp, Variant_Domain *z_mult_dat_Domain, double *z_mult_dat, Variant_Domain *sp_Domain, double *sp, Variant_Domain *et_Domain, double *et) {
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
  Variant_Grid_access_extracted(fp, x,y,z) -= Variant_Grid_access_extracted(z_mult_dat, x,y,z) * (Variant_Grid_access_extracted(sp, x,y,z) * Variant_Grid_access_extracted(et, x,y,z));
}

//NlFunctionEval:551 CUDA kernel 
__global__ void NLFE551Kernel(int xmax, int ymax, int zmax, Variant_Domain *x_ssl_dat_Domain, double *x_ssl_dat, Variant_Domain *y_ssl_dat_Domain, double *y_ssl_dat, Variant_Domain *pp_Domain, double *pp, Variant_Domain *z_mult_dat_Domain, double *z_mult_dat, Variant_Domain *dp_Domain, double *dp, Variant_Domain *permxp_Domain, double *permxp, Variant_Domain *rpp_Domain, double *rpp, Variant_Domain *permyp_Domain, double *permyp, Variant_Domain *permzp_Domain, double *permzp, Variant_Domain *vx_Domain, double *vx, Variant_Domain *vy_Domain, double *vy, Variant_Domain *vz_Domain, double *vz, Variant_Domain *u_right_Domain, double *u_right, Variant_Domain *u_front_Domain, double *u_front, Variant_Domain *u_upper_Domain, double *u_upper, Variant_Domain *fp_Domain, double *fp) {
  //Get position in kernel grid
  //+1 is added to coordinates to account for offset in each direction
  int x = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int y = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int z = blockIdx.z * blockDim.z + threadIdx.z + 1;

  //Make sure that the thread is in bounds
  if (x == 0 || x > xmax || y == 0 || y > ymax || z == 0 || z > zmax) {
    return;
  }

  double x_dir_g   = ArithmeticMean( Variant_Grid_access_extracted( x_ssl_dat, x, y, 0), Variant_Grid_access_extracted( x_ssl_dat, x+1,  y, 0 ) );
  double x_dir_g_c = ArithmeticMean( Variant_Grid_access_extracted( x_ssl_dat, x, y, 0), Variant_Grid_access_extracted( x_ssl_dat, x+1,  y, 0 ) );
  double y_dir_g   = ArithmeticMean( Variant_Grid_access_extracted( y_ssl_dat, x, y, 0), Variant_Grid_access_extracted( y_ssl_dat, x,  y+1, 0 ) );
  double y_dir_g_c = ArithmeticMean( Variant_Grid_access_extracted( y_ssl_dat, x, y, 0), Variant_Grid_access_extracted( y_ssl_dat, x,  y+1, 0 ) );

  double diff_right = Variant_Grid_access_extracted(pp, x,y,z) - Variant_Grid_access_extracted(pp, x+1, y,   z);
  double diff_front = Variant_Grid_access_extracted(pp, x,y,z) - Variant_Grid_access_extracted(pp, x,   y+1, z);

  double updir_right = diff_right * x_dir_g_c - x_dir_g;
  double updir_front = diff_front * y_dir_g_c - y_dir_g;

  double sep = ArithmeticMean( Variant_Grid_access_extracted(z_mult_dat, x,y,z),  Variant_Grid_access_extracted(z_mult_dat, x, y, z+1) );

  double lower_cond =
    Variant_Grid_access_extracted(pp, x,y,z) / sep
  - (   Variant_Grid_access_extracted(z_mult_dat, x,y,z)
      / ( Variant_Grid_access_extracted(z_mult_dat, x,y,z)
        + Variant_Grid_access_extracted(z_mult_dat, x, y, z+1)
        )
    )
  * Variant_Grid_access_extracted(dp, x,y,z);

  double upper_cond =
    Variant_Grid_access_extracted(pp, x, y, z+1) / sep
  + (   Variant_Grid_access_extracted(z_mult_dat, x, y, z+1)
      / ( Variant_Grid_access_extracted(z_mult_dat, x,y,z)
        + Variant_Grid_access_extracted(z_mult_dat, x, y, z+1)
        )
    )
  * Variant_Grid_access_extracted(dp, x,y,z+1);

  double diff_upper = lower_cond - upper_cond;

  /*
  NOTE!
  Originally the harmonic mean was
  PMean( Variant_Grid_access_extracted(pp, x,y,z),
          Variant_Grid_access_extracted(pp, x+1, y, z),
          permxp[ip],
          permxp[ip + 1]
  )
  However! PMean(a,b,c,d) in parflow is simply HarmonicMean(c, d)
  so the first two terms have been removed
  */

  double u_right_val =
    (
        Variant_Grid_access_extracted(z_mult_dat, x,y,z)
      * HarmonicMean( Variant_Grid_access_extracted(permxp, x,   y, z),
                      Variant_Grid_access_extracted(permxp, x+1, y, z) )
      * diff_right * x_dir_g_c
      * UpstreamMean(
          updir_right,
          0.0,
          Variant_Grid_access_extracted(rpp, x,   y, z) * Variant_Grid_access_extracted(dp, x,   y, z),
          Variant_Grid_access_extracted(rpp, x+1, y, z) * Variant_Grid_access_extracted(dp, x+1, y, z)
        )
    ) + (
      Variant_Grid_access_extracted(z_mult_dat, x,y,z)
      * HarmonicMean( Variant_Grid_access_extracted(permxp, x,   y, z),
                      Variant_Grid_access_extracted(permxp, x+1, y, z) )
      * (-x_dir_g)
      * UpstreamMean(
          updir_right,
          0.0,
          Variant_Grid_access_extracted(rpp, x,   y, z) * Variant_Grid_access_extracted(dp, x,   y, z),
          Variant_Grid_access_extracted(rpp, x+1, y, z) * Variant_Grid_access_extracted(dp, x+1, y, z)
        )
    );

    double u_front_val =
      (
        Variant_Grid_access_extracted(z_mult_dat, x,y,z)
      * HarmonicMean( Variant_Grid_access_extracted(permyp, x,   y, z),
                      Variant_Grid_access_extracted(permyp, x+1, y, z) )
      * diff_front * x_dir_g_c
      * UpstreamMean(
          updir_front,
          0.0,
          Variant_Grid_access_extracted(rpp, x,   y, z) * Variant_Grid_access_extracted(dp, x,   y, z),
          Variant_Grid_access_extracted(rpp, x+1, y, z) * Variant_Grid_access_extracted(dp, x+1, y, z)
        )
    ) + (
        Variant_Grid_access_extracted(z_mult_dat, x,y,z)
      * HarmonicMean( Variant_Grid_access_extracted(permyp, x,   y, z),
                      Variant_Grid_access_extracted(permyp, x+1, y, z) )
      * (-x_dir_g)
      * UpstreamMean(
          updir_front,
          0.0,
          Variant_Grid_access_extracted(rpp, x,   y, z) * Variant_Grid_access_extracted(dp, x,   y, z),
          Variant_Grid_access_extracted(rpp, x+1, y, z) * Variant_Grid_access_extracted(dp, x+1, y, z)
        )
    );

  double u_upper_val =
              HarmonicMeanDZ(
                  Variant_Grid_access_extracted(permzp,     x, y, z  ),
                  Variant_Grid_access_extracted(permzp,     x, y, z+1),
                  Variant_Grid_access_extracted(z_mult_dat, x, y, z  ),
                  Variant_Grid_access_extracted(z_mult_dat, x, y, z+1)
              )
            * diff_upper
            * UpstreamMean(
                lower_cond,
                upper_cond,
                Variant_Grid_access_extracted(rpp, x, y, z  ) * Variant_Grid_access_extracted(dp, x, y, z  ),
                Variant_Grid_access_extracted(rpp, x, y, z+1) * Variant_Grid_access_extracted(dp, x, y, z+1)
              );

  Variant_Grid_access_extracted(vx, x, y, z) = u_right_val;
  Variant_Grid_access_extracted(vy, x, y, z) = u_front_val;
  Variant_Grid_access_extracted(vz, x, y, z) = u_upper_val;

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

  Variant_Grid_access_extracted(u_right, x, y, z) = u_right_val;
  Variant_Grid_access_extracted(u_front, x, y, z) = u_front_val;
  Variant_Grid_access_extracted(u_upper, x, y, z) = u_upper_val;

  Variant_Grid_access_extracted(fp, x, y, z) += u_right_val * u_front_val * u_upper_val;
}

__global__ void NLFE551ReductionKernel(int xmax, int ymax, int zmax, Variant_Domain *u_right_Domain, double *u_right, Variant_Domain *u_front_Domain, double *u_front, Variant_Domain *u_upper_Domain, double *u_upper, Variant_Domain *fp_Domain, double *fp) {
  //Get position in kernel grid
  //+1 is added to coordinates to account for offset in each direction
  int x = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int y = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int z = blockIdx.z * blockDim.z + threadIdx.z + 1;

  //Make sure that the thread is in bounds
  if (x == 0 || x > xmax || y == 0 || y > ymax || z == 0 || z > zmax) {
    return;
  }

  double u_right_val =  Variant_Grid_access_extracted(u_right, x-1, y  , z  );
  double u_front_val = -Variant_Grid_access_extracted(u_front, x  , y-1, z  );
  double u_upper_val = -Variant_Grid_access_extracted(u_upper, x  , y  , z-1);
  Variant_Grid_access_extracted(fp, x, y, z) += u_right_val + u_front_val + u_upper_val;
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
  TIMEIT(metrics->elapsed_temp_alloc, {
    u_right = Variant_Grid_alloc( domain );
    u_front = Variant_Grid_alloc( domain );
    u_upper = Variant_Grid_alloc( domain );
  });

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
  domainSize = sizeof(Variant_Domain);
  dataSize = Basic_Domain_nx(domain) * Basic_Domain_ny(domain) * Basic_Domain_nz(domain) * sizeof(double);

  //Create streams for async execution
  cudaStream_t s1, s2, s3, s4, s5;
  check(cudaStreamCreate(&s1));
  check(cudaStreamCreate(&s2));
  check(cudaStreamCreate(&s3));
  check(cudaStreamCreate(&s4));
  check(cudaStreamCreate(&s5));

  //s1 vars
  Variant_Domain *fp_CUDA_Domain, *sp_CUDA_Domain, *dp_CUDA_Domain, *osp_CUDA_Domain, *odp_CUDA_Domain, *pop_CUDA_Domain, *z_mult_dat_CUDA_Domain;
  double *fp_CUDA, *sp_CUDA, *dp_CUDA, *osp_CUDA, *odp_CUDA, *pop_CUDA, *z_mult_dat_CUDA;

  //s2 vars
  Variant_Domain *ss_CUDA_Domain, *pp_CUDA_Domain, *opp_CUDA_Domain;
  double *ss_CUDA, *pp_CUDA, *opp_CUDA;

  //s3 vars
  Variant_Domain *et_CUDA_Domain;
  double *et_CUDA;

  //s4 vars
  Variant_Domain *x_ssl_dat_CUDA_Domain, *y_ssl_dat_CUDA_Domain, *permxp_CUDA_Domain, *rpp_CUDA_Domain, *permyp_CUDA_Domain, *permzp_CUDA_Domain, *vx_CUDA_Domain, *vy_CUDA_Domain, *vz_CUDA_Domain, *u_right_CUDA_Domain, *u_front_CUDA_Domain, *u_upper_CUDA_Domain;
  double *x_ssl_dat_CUDA, *y_ssl_dat_CUDA, *permxp_CUDA, *rpp_CUDA, *permyp_CUDA, *permzp_CUDA, *vx_CUDA, *vy_CUDA, *vz_CUDA, *u_right_CUDA, *u_front_CUDA, *u_upper_CUDA;

  //Allocate s1 vars
  check(cudaMalloc(&fp_CUDA_Domain, domainSize));
  check(cudaMalloc(&sp_CUDA_Domain, domainSize));
  check(cudaMalloc(&dp_CUDA_Domain, domainSize));
  check(cudaMalloc(&osp_CUDA_Domain, domainSize));
  check(cudaMalloc(&odp_CUDA_Domain, domainSize));
  check(cudaMalloc(&pop_CUDA_Domain, domainSize));
  check(cudaMalloc(&z_mult_dat_CUDA_Domain, domainSize));
  check(cudaMalloc(&fp_CUDA, dataSize));
  check(cudaMalloc(&sp_CUDA, dataSize));
  check(cudaMalloc(&dp_CUDA, dataSize));
  check(cudaMalloc(&osp_CUDA, dataSize));
  check(cudaMalloc(&odp_CUDA, dataSize));
  check(cudaMalloc(&pop_CUDA, dataSize));
  check(cudaMalloc(&z_mult_dat_CUDA, dataSize));

  //Allocate s2 vars
  check(cudaMalloc(&ss_CUDA_Domain, domainSize));
  check(cudaMalloc(&pp_CUDA_Domain, domainSize));
  check(cudaMalloc(&opp_CUDA_Domain, domainSize));
  check(cudaMalloc(&ss_CUDA, dataSize));
  check(cudaMalloc(&pp_CUDA, dataSize));
  check(cudaMalloc(&opp_CUDA, dataSize));

  //Allocate s3 vars
  check(cudaMalloc(&et_CUDA_Domain, domainSize));
  check(cudaMalloc(&et_CUDA, dataSize));

  //Allocate s4 vars
  check(cudaMalloc(&x_ssl_dat_CUDA_Domain, domainSize));
  check(cudaMalloc(&y_ssl_dat_CUDA_Domain, domainSize));
  check(cudaMalloc(&permxp_CUDA_Domain, domainSize));
  check(cudaMalloc(&rpp_CUDA_Domain, domainSize));
  check(cudaMalloc(&permyp_CUDA_Domain, domainSize));
  check(cudaMalloc(&permzp_CUDA_Domain, domainSize));
  check(cudaMalloc(&vx_CUDA_Domain, domainSize));
  check(cudaMalloc(&vy_CUDA_Domain, domainSize));
  check(cudaMalloc(&vz_CUDA_Domain, domainSize));
  check(cudaMalloc(&u_right_CUDA_Domain, domainSize));
  check(cudaMalloc(&u_front_CUDA_Domain, domainSize));
  check(cudaMalloc(&u_upper_CUDA_Domain, domainSize));
  check(cudaMalloc(&x_ssl_dat_CUDA, dataSize));
  check(cudaMalloc(&y_ssl_dat_CUDA, dataSize));
  check(cudaMalloc(&permxp_CUDA, dataSize));
  check(cudaMalloc(&rpp_CUDA, dataSize));
  check(cudaMalloc(&permyp_CUDA, dataSize));
  check(cudaMalloc(&permzp_CUDA, dataSize));
  check(cudaMalloc(&vx_CUDA, dataSize));
  check(cudaMalloc(&vy_CUDA, dataSize));
  check(cudaMalloc(&vz_CUDA, dataSize));
  check(cudaMalloc(&u_right_CUDA, dataSize));
  check(cudaMalloc(&u_front_CUDA, dataSize));
  check(cudaMalloc(&u_upper_CUDA, dataSize));

  //Async copy s1 vars
  check(cudaMemcpyAsync(fp_CUDA_Domain, fp->domain, domainSize, cudaMemcpyHostToDevice, s1));
  check(cudaMemcpyAsync(sp_CUDA_Domain, sp->domain, domainSize, cudaMemcpyHostToDevice, s1));
  check(cudaMemcpyAsync(dp_CUDA_Domain, dp->domain, domainSize, cudaMemcpyHostToDevice, s1));
  check(cudaMemcpyAsync(osp_CUDA_Domain, osp->domain, domainSize, cudaMemcpyHostToDevice, s1));
  check(cudaMemcpyAsync(odp_CUDA_Domain, odp->domain, domainSize, cudaMemcpyHostToDevice, s1));
  check(cudaMemcpyAsync(pop_CUDA_Domain, pop->domain, domainSize, cudaMemcpyHostToDevice, s1));
  check(cudaMemcpyAsync(z_mult_dat_CUDA_Domain, z_mult_dat->domain, domainSize, cudaMemcpyHostToDevice, s1));
  check(cudaMemcpyAsync(fp_CUDA, fp->data, dataSize, cudaMemcpyHostToDevice, s1));
  check(cudaMemcpyAsync(sp_CUDA, sp->data, dataSize, cudaMemcpyHostToDevice, s1));
  check(cudaMemcpyAsync(dp_CUDA, dp->data, dataSize, cudaMemcpyHostToDevice, s1));
  check(cudaMemcpyAsync(osp_CUDA, osp->data, dataSize, cudaMemcpyHostToDevice, s1));
  check(cudaMemcpyAsync(odp_CUDA, odp->data, dataSize, cudaMemcpyHostToDevice, s1));
  check(cudaMemcpyAsync(pop_CUDA, pop->data, dataSize, cudaMemcpyHostToDevice, s1));
  check(cudaMemcpyAsync(z_mult_dat_CUDA, z_mult_dat->data, dataSize, cudaMemcpyHostToDevice, s1));

  //Async copy s2 vars
  check(cudaMemcpyAsync(ss_CUDA_Domain, ss->domain, domainSize, cudaMemcpyHostToDevice, s2));
  check(cudaMemcpyAsync(pp_CUDA_Domain, pp->domain, domainSize, cudaMemcpyHostToDevice, s2));
  check(cudaMemcpyAsync(opp_CUDA_Domain, opp->domain, domainSize, cudaMemcpyHostToDevice, s2));
  check(cudaMemcpyAsync(ss_CUDA, ss->data, dataSize, cudaMemcpyHostToDevice, s2));
  check(cudaMemcpyAsync(pp_CUDA, pp->data, dataSize, cudaMemcpyHostToDevice, s2));
  check(cudaMemcpyAsync(opp_CUDA, opp->data, dataSize, cudaMemcpyHostToDevice, s2));

  //Async copy s3 vars
  check(cudaMemcpyAsync(et_CUDA_Domain, et->domain, domainSize, cudaMemcpyHostToDevice, s3));
  check(cudaMemcpyAsync(et_CUDA, et->data, dataSize, cudaMemcpyHostToDevice, s3));

  //Async copy s4 vars
  check(cudaMemcpyAsync(x_ssl_dat_CUDA_Domain, x_ssl_dat->domain, domainSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(y_ssl_dat_CUDA_Domain, y_ssl_dat->domain, domainSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(permxp_CUDA_Domain, permxp->domain, domainSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(rpp_CUDA_Domain, rpp->domain, domainSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(permyp_CUDA_Domain, permyp->domain, domainSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(permzp_CUDA_Domain, permzp->domain, domainSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(vx_CUDA_Domain, vx->domain, domainSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(vy_CUDA_Domain, vy->domain, domainSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(vz_CUDA_Domain, vz->domain, domainSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(u_right_CUDA_Domain, u_right->domain, domainSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(u_front_CUDA_Domain, u_front->domain, domainSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(u_upper_CUDA_Domain, u_upper->domain, domainSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(x_ssl_dat_CUDA, x_ssl_dat->data, dataSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(y_ssl_dat_CUDA, y_ssl_dat->data, dataSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(permxp_CUDA, permxp->data, dataSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(rpp_CUDA, rpp->data, dataSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(permyp_CUDA, permyp->data, dataSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(permzp_CUDA, permzp->data, dataSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(vx_CUDA, vx->data, dataSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(vy_CUDA, vy->data, dataSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(vz_CUDA, vz->data, dataSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(u_right_CUDA, u_right->data, dataSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(u_front_CUDA, u_front->data, dataSize, cudaMemcpyHostToDevice, s4));
  check(cudaMemcpyAsync(u_upper_CUDA, u_upper->data, dataSize, cudaMemcpyHostToDevice, s4));

  NLFE216Kernel<<<grid, block, 0, s1>>>(Basic_Domain_nx(domain) - 2, Basic_Domain_ny(domain) - 2, Basic_Domain_nz(domain) -2, fp_CUDA_Domain, fp_CUDA, sp_CUDA_Domain, sp_CUDA, dp_CUDA_Domain, dp_CUDA, osp_CUDA_Domain, osp_CUDA, odp_CUDA_Domain, odp_CUDA, pop_CUDA_Domain, pop_CUDA, z_mult_dat_CUDA_Domain, z_mult_dat_CUDA);
  check(cudaStreamSynchronize(s1));

  NLFE338Kernel<<<grid, block, 0, s2>>>(Basic_Domain_nx(domain) - 2, Basic_Domain_ny(domain) - 2, Basic_Domain_nz(domain) -2, fp_CUDA_Domain, fp_CUDA, ss_CUDA_Domain, ss_CUDA, z_mult_dat_CUDA_Domain, z_mult_dat_CUDA, pp_CUDA_Domain, pp_CUDA, sp_CUDA_Domain, sp_CUDA, dp_CUDA_Domain, dp_CUDA, opp_CUDA_Domain, opp_CUDA, osp_CUDA_Domain, osp_CUDA, odp_CUDA_Domain, odp_CUDA);
  check(cudaStreamSynchronize(s2));

  NLFE416Kernel<<<grid, block, 0, s3>>>(Basic_Domain_nx(domain) - 2, Basic_Domain_ny(domain) - 2, Basic_Domain_nz(domain) -2, fp_CUDA_Domain, fp_CUDA, z_mult_dat_CUDA_Domain, z_mult_dat_CUDA, sp_CUDA_Domain, sp_CUDA, et_CUDA_Domain, et_CUDA);
  check(cudaStreamSynchronize(s3));

  NLFE551Kernel<<<grid, block, 0, s4>>>(Basic_Domain_nx(domain) - 2, Basic_Domain_ny(domain) - 2, Basic_Domain_nz(domain) -2, x_ssl_dat_CUDA_Domain, x_ssl_dat_CUDA, y_ssl_dat_CUDA_Domain, y_ssl_dat_CUDA, pp_CUDA_Domain, pp_CUDA, z_mult_dat_CUDA_Domain, z_mult_dat_CUDA, dp_CUDA_Domain, dp_CUDA, permxp_CUDA_Domain, permxp_CUDA, rpp_CUDA_Domain, rpp_CUDA, permyp_CUDA_Domain, permyp_CUDA, permzp_CUDA_Domain, permzp_CUDA, vx_CUDA_Domain, vx_CUDA, vy_CUDA_Domain, vy_CUDA, vz_CUDA_Domain, vz_CUDA, u_right_CUDA_Domain, u_right_CUDA, u_front_CUDA_Domain, u_front_CUDA, u_upper_CUDA_Domain, u_upper_CUDA, fp_CUDA_Domain, fp_CUDA);
  check(cudaStreamSynchronize(s4));

  NLFE551ReductionKernel<<<grid, block, 0, s5>>>(Basic_Domain_nx(reduction_domain) - 2, Basic_Domain_ny(reduction_domain) - 2, Basic_Domain_nz(reduction_domain) -2, u_right_CUDA_Domain, u_right_CUDA, u_front_CUDA_Domain, u_front_CUDA, u_upper_CUDA_Domain, u_upper_CUDA, fp_CUDA_Domain, fp_CUDA);
  check(cudaStreamSynchronize(s5));

/* ------------------------------ Cleanup ------------------------------ */

  //Destroy async streams
  check(cudaStreamDestroy(s1));
  check(cudaStreamDestroy(s2));
  check(cudaStreamDestroy(s3));
  check(cudaStreamDestroy(s4));
  check(cudaStreamDestroy(s5));

  //Copy result arrays back over
  check(cudaMemcpy(fp->data, fp_CUDA, dataSize, cudaMemcpyDeviceToHost));
  check(cudaMemcpy(vx->data, vx_CUDA, dataSize, cudaMemcpyDeviceToHost));
  check(cudaMemcpy(vy->data, vy_CUDA, dataSize, cudaMemcpyDeviceToHost));
  check(cudaMemcpy(vz->data, vz_CUDA, dataSize, cudaMemcpyDeviceToHost));

  //Free all allocated device vars
  check(cudaFree(fp_CUDA_Domain));
  check(cudaFree(vx_CUDA_Domain));
  check(cudaFree(vy_CUDA_Domain));
  check(cudaFree(vz_CUDA_Domain));
  check(cudaFree(dp_CUDA_Domain));
  check(cudaFree(et_CUDA_Domain));
  check(cudaFree(odp_CUDA_Domain));
  check(cudaFree(opp_CUDA_Domain));
  check(cudaFree(osp_CUDA_Domain));
  check(cudaFree(permxp_CUDA_Domain));
  check(cudaFree(permyp_CUDA_Domain));
  check(cudaFree(permzp_CUDA_Domain));
  check(cudaFree(pop_CUDA_Domain));
  check(cudaFree(pp_CUDA_Domain));
  check(cudaFree(rpp_CUDA_Domain));
  check(cudaFree(sp_CUDA_Domain));
  check(cudaFree(ss_CUDA_Domain));
  check(cudaFree(z_mult_dat_CUDA_Domain));
  check(cudaFree(x_ssl_dat_CUDA_Domain));
  check(cudaFree(y_ssl_dat_CUDA_Domain));
  check(cudaFree(u_right_CUDA_Domain));
  check(cudaFree(u_front_CUDA_Domain));
  check(cudaFree(u_upper_CUDA_Domain));
  check(cudaFree(fp_CUDA));
  check(cudaFree(vx_CUDA));
  check(cudaFree(vy_CUDA));
  check(cudaFree(vz_CUDA));
  check(cudaFree(dp_CUDA));
  check(cudaFree(et_CUDA));
  check(cudaFree(odp_CUDA));
  check(cudaFree(opp_CUDA));
  check(cudaFree(osp_CUDA));
  check(cudaFree(permxp_CUDA));
  check(cudaFree(permyp_CUDA));
  check(cudaFree(permzp_CUDA));
  check(cudaFree(pop_CUDA));
  check(cudaFree(pp_CUDA));
  check(cudaFree(rpp_CUDA));
  check(cudaFree(sp_CUDA));
  check(cudaFree(ss_CUDA));
  check(cudaFree(z_mult_dat_CUDA));
  check(cudaFree(x_ssl_dat_CUDA));
  check(cudaFree(y_ssl_dat_CUDA));
  check(cudaFree(u_right_CUDA));
  check(cudaFree(u_front_CUDA));
  check(cudaFree(u_upper_CUDA));
  
  //Cleanup temp domain and grids
  Variant_Domain_dealloc( reduction_domain );
  TIMEIT(metrics->elapsed_temp_dealloc, {
    Variant_Grid_dealloc(u_right);
    Variant_Grid_dealloc(u_front);
    Variant_Grid_dealloc(u_upper);
  });
}