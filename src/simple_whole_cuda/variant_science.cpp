#include <science.hpp>
#include <macros.hpp>
#include <configure.hpp>
#include <util.hpp>
#include <math.h>
#include <assert.h>
#include <stdio.h>

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
  
  Variant_Grid* u_right;
  Variant_Grid* u_front;
  Variant_Grid* u_upper;

  // Create temporary grids for reduction portion.
  TIMEIT(metrics->elapsed_temp_grid_alloc, {
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

/* ------------------------------ NLFE216 ------------------------------ */

  TIMEIT(metrics->elapsed_216, {
    NLFE216Kernel<<<grid, block>>>(Basic_Domain_nx(domain) - 2, Basic_Domain_ny(domain) - 2, Basic_Domain_nz(domain) -2, fp, sp, dp, osp, odp, pop, z_mult_dat);
    cudaDeviceSynchronize();
  });

/* ------------------------------ NLFE338 ------------------------------ */

  TIMEIT(metrics->elapsed_338, {
    NLFE338Kernel<<<grid, block>>>(Basic_Domain_nx(domain) - 2, Basic_Domain_ny(domain) - 2, Basic_Domain_nz(domain) -2, fp, ss, z_mult_dat, pp, sp, dp, opp, osp, odp);
    cudaDeviceSynchronize();
  });

/* ------------------------------ NLFE416 ------------------------------ */

  TIMEIT(metrics->elapsed_416, {
    NLFE416Kernel<<<grid, block>>>(Basic_Domain_nx(domain) - 2, Basic_Domain_ny(domain) - 2, Basic_Domain_nz(domain) -2, fp, z_mult_dat, sp, et);
    cudaDeviceSynchronize();
  });

/* ------------------------------ NLFE551 ------------------------------ */
  
  TIMEIT(metrics->elapsed_551, {
    NLFE551Kernel<<<grid, block>>>(Basic_Domain_nx(domain) - 2, Basic_Domain_ny(domain) - 2, Basic_Domain_nz(domain) -2, x_ssl_dat, y_ssl_dat, pp, z_mult_dat, dp, permxp, rpp, permyp, permzp, vx, vy, vz, u_right, u_front, u_upper, fp);
    cudaDeviceSynchronize();
  });
  
/* ------------------------------ NLFE551 Reduction ------------------------------ */

  TIMEIT(metrics->elapsed_551_reduce, {
    NLFE551ReductionKernel<<<grid, block>>>(Basic_Domain_nx(reduction_domain) - 2, Basic_Domain_ny(reduction_domain) - 2, Basic_Domain_nz(reduction_domain) -2, u_right, u_front, u_upper, fp);
    cudaDeviceSynchronize();
  });

/* ------------------------------ Cleanup ------------------------------ */

  Variant_Domain_dealloc( reduction_domain );

  TIMEIT(metrics->elapsed_temp_grid_dealloc, {
    Variant_Grid_dealloc(u_right);
    Variant_Grid_dealloc(u_front);
    Variant_Grid_dealloc(u_upper);
  });
}