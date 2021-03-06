#include <science.hpp>
#include <macros.hpp>
#include <configure.hpp>
#include <util.hpp>

#include <assert.h>
#include <stdio.h>

// NOTE: THIS IS REQUIRED FOR ROSE AND LCTOOL
void* nullptr = NULL;

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

  // Do baseline scientific kernel
  int x, y, z;
  int nx = Variant_Domain_nx( domain );
  int ny = Variant_Domain_ny( domain );
  int nz = Variant_Domain_nz( domain );

  #pragma omplc loopchain schedule( fuse() )
  {
    // NlFunctionEval:261 analogue
    #pragma omplc for domain(0:nx-1, 0:ny-1, 0:nz-1) \
      with (x,y,z) \
        read sp { (x,y,z) }, \
        read dp { (x,y,z) }, \
        read osp { (x,y,z) }, \
        read odp { (x,y,z) }, \
        read pop { (x,y,z) }, \
        read z_mult_dat { (x,y,z) }, \
        write fp { (x,y,z) }
    Variant_Domain_fast_loop_interior(domain, x,y,z,
      {
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
    #pragma omplc for domain(0:nx-1, 0:ny-1, 0:nz-1) \
      with (x,y,z) \
        read fp { (x,y,z) }, \
        read ss { (x,y,z) }, \
        read z_mult_dat { (x,y,z) }, \
        read pp { (x,y,z) }, \
        read sp { (x,y,z) }, \
        read dp { (x,y,z) }, \
        read opp { (x,y,z) }, \
        read osp { (x,y,z) }, \
        read odp { (x,y,z) }, \
        write fp { (x,y,z) }
    Variant_Domain_fast_loop_interior(domain, x,y,z,
      {
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
    #pragma omplc for domain(0:nx-1, 0:ny-1, 0:nz-1) \
      with (x,y,z) \
        read fp { (x,y,z) }, \
        read z_mult_dat { (x,y,z) }, \
        read sp { (x,y,z) }, \
        read et { (x,y,z) }, \
        write fp { (x,y,z) }
    Variant_Domain_fast_loop_interior(domain, x,y,z,
      {
        Variant_Grid_access(fp, x,y,z) -=
            Variant_Grid_access(z_mult_dat, x,y,z)
          * (   Variant_Grid_access(sp, x,y,z)
              * Variant_Grid_access(et, x,y,z)
            );
      }
    );

  // NlFunctionEval:551 analogue
    #pragma omplc for domain(0:nx-1, 0:ny-1, 0:nz-1) \
      with (x,y,z) \
        read fp { (x,y,z) }, \
        read x_ssl_dat { (x,y,0), (x+1, y,   0) }, \
        read y_ssl_dat { (x,y,0), (x,   y+1, 0) }, \
        read pp { (x,y,z), (x+1,y,z), (x,y+1,z), (x,y,z+1) }, \
        read z_mult_dat { (x,y,z), (x,y,z+1) }, \
        read dp { (x,y,z), (x,y,z+1), (x+1,y,z) }, \
        read permxp { (x,y,z), (x+1,y,z) }, \
        read permzp { (x,y,z), (x+1,y,z) }, \
        read rpp { (x,y,z), (x+1,y,z) }, \
        write vx { (x,y,z) }, \
        write vy { (x,y,z) }, \
        write vz { (x,y,z) }, \
        write fp { (x,y,z), (x+1,y,z), (x,y+1,z), (x,y,z+1) }
    Variant_Domain_fast_loop_interior(domain, x,y,z,
      {

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
  }
}
