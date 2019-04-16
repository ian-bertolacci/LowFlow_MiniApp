#include <science.hpp>
#include <macros.hpp>
#include <configure.hpp>
#include <util.hpp>

#include <assert.h>

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

  // Define all variables ahead like in ParFlow
  int x, y, z, idx, idx_alt, x_ssl_dat_xm, y_ssl_dat_ym, fp_xm, fp_ym, fp_zm;
  double x_dir_g, x_dir_g_c, y_dir_g, y_dir_g_c, diff_right, diff_front, updir_right, updir_front, sep, lower_cond, upper_cond, diff_upper, u_right, u_front, u_upper;

  double *fp_data, *vx_data, *vy_data, *vz_data, *dp_data, *et_data, *odp_data, *opp_data, *osp_data, *permxp_data, *permyp_data, *permzp_data, *pop_data, *pp_data, *rpp_data, *sp_data, *ss_data, *z_mult_dat_data, *x_ssl_dat_data, *y_ssl_dat_data;

  // Do baseline scientific kernel
  // NlFunctionEval:261 analogue
  {
    fp_data = Variant_Grid_data(fp);
    vx_data = Variant_Grid_data(vx);
    vy_data = Variant_Grid_data(vy);
    vz_data = Variant_Grid_data(vz);
    dp_data = Variant_Grid_data(dp);
    et_data = Variant_Grid_data(et);
    odp_data = Variant_Grid_data(odp);
    opp_data = Variant_Grid_data(opp);
    osp_data = Variant_Grid_data(osp);
    permxp_data = Variant_Grid_data(permxp);
    permyp_data = Variant_Grid_data(permyp);
    permzp_data = Variant_Grid_data(permzp);
    pop_data = Variant_Grid_data(pop);
    pp_data = Variant_Grid_data(pp);
    rpp_data = Variant_Grid_data(rpp);
    sp_data = Variant_Grid_data(sp);
    ss_data = Variant_Grid_data(ss);
    z_mult_dat_data = Variant_Grid_data(z_mult_dat);
    x_ssl_dat_data = Variant_Grid_data(x_ssl_dat);
    y_ssl_dat_data = Variant_Grid_data(y_ssl_dat);

    Variant_Domain_fast_loop_interior(domain, x,y,z,
      {
        idx = Variant_Grid_idx(fp, x,y,z);
        fp_data[idx] =
         (  sp_data[idx]
          * dp_data[idx]
          - osp_data[idx]
          * odp_data[idx]
         )
         * pop_data[idx]
         * z_mult_dat_data[idx];
      }
    );
  }

  // sNlFunctionEval:338 analogue
  {
    fp_data = Variant_Grid_data(fp);
    vx_data = Variant_Grid_data(vx);
    vy_data = Variant_Grid_data(vy);
    vz_data = Variant_Grid_data(vz);
    dp_data = Variant_Grid_data(dp);
    et_data = Variant_Grid_data(et);
    odp_data = Variant_Grid_data(odp);
    opp_data = Variant_Grid_data(opp);
    osp_data = Variant_Grid_data(osp);
    permxp_data = Variant_Grid_data(permxp);
    permyp_data = Variant_Grid_data(permyp);
    permzp_data = Variant_Grid_data(permzp);
    pop_data = Variant_Grid_data(pop);
    pp_data = Variant_Grid_data(pp);
    rpp_data = Variant_Grid_data(rpp);
    sp_data = Variant_Grid_data(sp);
    ss_data = Variant_Grid_data(ss);
    z_mult_dat_data = Variant_Grid_data(z_mult_dat);
    x_ssl_dat_data = Variant_Grid_data(x_ssl_dat);
    y_ssl_dat_data = Variant_Grid_data(y_ssl_dat);
    Variant_Domain_fast_loop_interior(domain, x,y,z,
      {
        idx = Variant_Grid_idx(fp, x,y,z);
        fp_data[idx] +=
            ss_data[idx]
          * z_mult_dat_data[idx]
          * (   pp_data[idx]
              * sp_data[idx]
              * dp_data[idx]
              - opp_data[idx]
              * osp_data[idx]
              * odp_data[idx]
            );
      }
    );
  }

  // NlFunctionEval:416 analogue
  {
    fp_data = Variant_Grid_data(fp);
    vx_data = Variant_Grid_data(vx);
    vy_data = Variant_Grid_data(vy);
    vz_data = Variant_Grid_data(vz);
    dp_data = Variant_Grid_data(dp);
    et_data = Variant_Grid_data(et);
    odp_data = Variant_Grid_data(odp);
    opp_data = Variant_Grid_data(opp);
    osp_data = Variant_Grid_data(osp);
    permxp_data = Variant_Grid_data(permxp);
    permyp_data = Variant_Grid_data(permyp);
    permzp_data = Variant_Grid_data(permzp);
    pop_data = Variant_Grid_data(pop);
    pp_data = Variant_Grid_data(pp);
    rpp_data = Variant_Grid_data(rpp);
    sp_data = Variant_Grid_data(sp);
    ss_data = Variant_Grid_data(ss);
    z_mult_dat_data = Variant_Grid_data(z_mult_dat);
    x_ssl_dat_data = Variant_Grid_data(x_ssl_dat);
    y_ssl_dat_data = Variant_Grid_data(y_ssl_dat);
    Variant_Domain_fast_loop_interior(domain, x,y,z,
      {
        idx = Variant_Grid_idx(fp, x,y,z);
        fp_data[idx] -=
            z_mult_dat_data[idx]
          * (   sp_data[idx]
              * et_data[idx]
            );
      }
    );
  }

  // NlFunctionEval:551 analogue
  {
    fp_data = Variant_Grid_data(fp);
    vx_data = Variant_Grid_data(vx);
    vy_data = Variant_Grid_data(vy);
    vz_data = Variant_Grid_data(vz);
    dp_data = Variant_Grid_data(dp);
    et_data = Variant_Grid_data(et);
    odp_data = Variant_Grid_data(odp);
    opp_data = Variant_Grid_data(opp);
    osp_data = Variant_Grid_data(osp);
    permxp_data = Variant_Grid_data(permxp);
    permyp_data = Variant_Grid_data(permyp);
    permzp_data = Variant_Grid_data(permzp);
    pop_data = Variant_Grid_data(pop);
    pp_data = Variant_Grid_data(pp);
    rpp_data = Variant_Grid_data(rpp);
    sp_data = Variant_Grid_data(sp);
    ss_data = Variant_Grid_data(ss);
    z_mult_dat_data = Variant_Grid_data(z_mult_dat);
    x_ssl_dat_data = Variant_Grid_data(x_ssl_dat);
    y_ssl_dat_data = Variant_Grid_data(y_ssl_dat);

    x_ssl_dat_xm = 1;
    y_ssl_dat_ym = Variant_Domain_nx( Variant_Grid_domain(y_ssl_dat) );
    fp_xm = 1;
    fp_ym = Variant_Domain_nx( Variant_Grid_domain(fp) );
    fp_zm = Variant_Domain_nx( Variant_Grid_domain(fp) ) * Variant_Domain_ny( Variant_Grid_domain(fp) );
    Variant_Domain_fast_loop_interior(domain, x,y,z,
      {
        idx  = Variant_Grid_idx(fp, x,y,z);
        idx_alt = Variant_Grid_idx(x_ssl_dat, x,y,0);
        x_dir_g   = ArithmeticMean( x_ssl_dat_data[idx_alt], x_ssl_dat_data[idx_alt+x_ssl_dat_xm] ); //Variant_Grid_access( x_ssl_dat, x+1,  y, 0 ) );
        x_dir_g_c = ArithmeticMean( x_ssl_dat_data[idx_alt], x_ssl_dat_data[idx_alt+x_ssl_dat_xm] ); //Variant_Grid_access( x_ssl_dat, x+1,  y, 0 ) );
        y_dir_g   = ArithmeticMean( y_ssl_dat_data[idx_alt], y_ssl_dat_data[idx_alt+y_ssl_dat_ym] ); //Variant_Grid_access( y_ssl_dat, x,  y+1, 0 ) );
        y_dir_g_c = ArithmeticMean( y_ssl_dat_data[idx_alt], y_ssl_dat_data[idx_alt+y_ssl_dat_ym] ); //Variant_Grid_access( y_ssl_dat, x,  y+1, 0 ) );

        diff_right = pp_data[idx] - pp_data[idx+fp_xm];
        diff_front = pp_data[idx] - pp_data[idx+fp_ym];

        updir_right = diff_right * x_dir_g_c - x_dir_g;
        updir_front = diff_front * y_dir_g_c - y_dir_g;

        sep = ArithmeticMean( z_mult_dat_data[idx],  z_mult_dat_data[idx+fp_zm] );

        lower_cond =
          pp_data[idx] / sep
        - (   z_mult_dat_data[idx]
            / ( z_mult_dat_data[idx]
              + z_mult_dat_data[idx+fp_zm]
              )
          )
        * dp_data[idx];

        upper_cond =
          pp_data[idx+fp_zm] / sep
        + (   z_mult_dat_data[idx+fp_zm]
            / ( z_mult_dat_data[idx]
              + z_mult_dat_data[idx+fp_zm]
              )
          )
        * dp_data[idx+fp_zm];

        diff_upper = lower_cond - upper_cond;

        /*
        NOTE!
        Originally the harmonic mean was
        PMean( pp_data[idx],
               Variant_Grid_access(pp, x+1, y, z),
               permxp[ip],
               permxp[ip + 1]
        )
        However! PMean(a,b,c,d) in parflow is simply HarmonicMean(c, d)
        so the first two terms have been removed
        */

        u_right =
          (
             z_mult_dat_data[idx]
           * HarmonicMean( permxp_data[idx],
                           permxp_data[idx+fp_xm] )
           * diff_right * x_dir_g_c
           * UpstreamMean(
               updir_right,
               0.0,
               rpp_data[idx] * dp_data[idx],
               rpp_data[idx+fp_xm] * dp_data[idx+fp_xm]
             )
         ) + (
            z_mult_dat_data[idx]
           * HarmonicMean( permxp_data[idx],
                           permxp_data[idx+fp_xm] )
           * (-x_dir_g)
           * UpstreamMean(
               updir_right,
               0.0,
               rpp_data[idx] * dp_data[idx],
               rpp_data[idx+fp_xm] * dp_data[idx+fp_xm]
             )
         );

         u_front =
           (
              z_mult_dat_data[idx]
            * HarmonicMean( permyp_data[idx],
                            permyp_data[idx+fp_xm] )
            * diff_front * x_dir_g_c
            * UpstreamMean(
                updir_front,
                0.0,
                rpp_data[idx] * dp_data[idx],
                rpp_data[idx+fp_xm] * dp_data[idx+fp_xm]
              )
          ) + (
             z_mult_dat_data[idx]
            * HarmonicMean( permyp_data[idx],
                            permyp_data[idx+fp_xm] )
            * (-x_dir_g)
            * UpstreamMean(
                updir_front,
                0.0,
                rpp_data[idx] * dp_data[idx],
                rpp_data[idx+fp_xm] * dp_data[idx+fp_xm]
              )
          );

        u_upper =
                    HarmonicMeanDZ(
                        permzp_data[idx],
                        permzp_data[idx+fp_zm],
                        z_mult_dat_data[idx],
                        z_mult_dat_data[idx+fp_zm]
                    )
                  * diff_upper
                  * UpstreamMean(
                      lower_cond,
                      upper_cond,
                      rpp_data[idx] * dp_data[idx],
                      rpp_data[idx+fp_zm] * dp_data[idx+fp_zm]
                    );

        vx_data[idx] = u_right;
        vy_data[idx] = u_front;
        vz_data[idx] = u_upper;

        fp_data[idx] += u_right * u_front * u_upper;
        fp_data[idx+fp_xm] += u_right;
        fp_data[idx+fp_ym] -= u_front;
        fp_data[idx+fp_zm] -= u_upper;
      }
    );
  }
}
