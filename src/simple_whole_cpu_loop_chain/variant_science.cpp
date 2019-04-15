#include <assert.h>
#include <configure.hpp>
#include <macros.hpp>
#include <science.hpp>
#include <util.hpp>

void science(Variant_Domain *domain, Variant_Grid *fp, Variant_Grid *vx,
             Variant_Grid *vy, Variant_Grid *vz, Variant_Grid *dp,
             Variant_Grid *et, Variant_Grid *odp, Variant_Grid *opp,
             Variant_Grid *osp, Variant_Grid *permxp, Variant_Grid *permyp,
             Variant_Grid *permzp, Variant_Grid *pop, Variant_Grid *pp,
             Variant_Grid *rpp, Variant_Grid *sp, Variant_Grid *ss,
             Variant_Grid *z_mult_dat, Variant_Grid *x_ssl_dat,
             Variant_Grid *y_ssl_dat) {
  if (true) {
    domain != nullptr && "Error during science: domain is null"
        ? (static_cast<void>(0))
        : __assert_fail(
              "domain != nullptr && \"Error during science: domain is null\"",
              "variant_science.in.cpp", 35, __PRETTY_FUNCTION__);
    fp != nullptr && "Error during science: output grid fp is null"
        ? (static_cast<void>(0))
        : __assert_fail("fp != nullptr && \"Error during science: output grid "
                        "fp is null\"",
                        "variant_science.in.cpp", 36, __PRETTY_FUNCTION__);
    vx != nullptr && "Error during science: output grid vx is null"
        ? (static_cast<void>(0))
        : __assert_fail("vx != nullptr && \"Error during science: output grid "
                        "vx is null\"",
                        "variant_science.in.cpp", 37, __PRETTY_FUNCTION__);
    vy != nullptr && "Error during science: output grid vy is null"
        ? (static_cast<void>(0))
        : __assert_fail("vy != nullptr && \"Error during science: output grid "
                        "vy is null\"",
                        "variant_science.in.cpp", 38, __PRETTY_FUNCTION__);
    vz != nullptr && "Error during science: output grid vz is null"
        ? (static_cast<void>(0))
        : __assert_fail("vz != nullptr && \"Error during science: output grid "
                        "vz is null\"",
                        "variant_science.in.cpp", 39, __PRETTY_FUNCTION__);
    dp != nullptr && "Error during science: input grid dp is null"
        ? (static_cast<void>(0))
        : __assert_fail("dp != nullptr && \"Error during science: input grid "
                        "dp is null\"",
                        "variant_science.in.cpp", 40, __PRETTY_FUNCTION__);
    et != nullptr && "Error during science: input grid et is null"
        ? (static_cast<void>(0))
        : __assert_fail("et != nullptr && \"Error during science: input grid "
                        "et is null\"",
                        "variant_science.in.cpp", 41, __PRETTY_FUNCTION__);
    odp != nullptr && "Error during science: input grid odp is null"
        ? (static_cast<void>(0))
        : __assert_fail("odp != nullptr && \"Error during science: input grid "
                        "odp is null\"",
                        "variant_science.in.cpp", 42, __PRETTY_FUNCTION__);
    opp != nullptr && "Error during science: input grid opp is null"
        ? (static_cast<void>(0))
        : __assert_fail("opp != nullptr && \"Error during science: input grid "
                        "opp is null\"",
                        "variant_science.in.cpp", 43, __PRETTY_FUNCTION__);
    osp != nullptr && "Error during science: input grid osp is null"
        ? (static_cast<void>(0))
        : __assert_fail("osp != nullptr && \"Error during science: input grid "
                        "osp is null\"",
                        "variant_science.in.cpp", 44, __PRETTY_FUNCTION__);
    permxp != nullptr && "Error during science: input grid permxp is null"
        ? (static_cast<void>(0))
        : __assert_fail("permxp != nullptr && \"Error during science: input "
                        "grid permxp is null\"",
                        "variant_science.in.cpp", 45, __PRETTY_FUNCTION__);
    permyp != nullptr && "Error during science: input grid permyp is null"
        ? (static_cast<void>(0))
        : __assert_fail("permyp != nullptr && \"Error during science: input "
                        "grid permyp is null\"",
                        "variant_science.in.cpp", 46, __PRETTY_FUNCTION__);
    permzp != nullptr && "Error during science: input grid permzp is null"
        ? (static_cast<void>(0))
        : __assert_fail("permzp != nullptr && \"Error during science: input "
                        "grid permzp is null\"",
                        "variant_science.in.cpp", 47, __PRETTY_FUNCTION__);
    pop != nullptr && "Error during science: input grid pop is null"
        ? (static_cast<void>(0))
        : __assert_fail("pop != nullptr && \"Error during science: input grid "
                        "pop is null\"",
                        "variant_science.in.cpp", 48, __PRETTY_FUNCTION__);
    pp != nullptr && "Error during science: input grid pp is null"
        ? (static_cast<void>(0))
        : __assert_fail("pp != nullptr && \"Error during science: input grid "
                        "pp is null\"",
                        "variant_science.in.cpp", 49, __PRETTY_FUNCTION__);
    rpp != nullptr && "Error during science: input grid rpp is null"
        ? (static_cast<void>(0))
        : __assert_fail("rpp != nullptr && \"Error during science: input grid "
                        "rpp is null\"",
                        "variant_science.in.cpp", 50, __PRETTY_FUNCTION__);
    sp != nullptr && "Error during science: input grid sp is null"
        ? (static_cast<void>(0))
        : __assert_fail("sp != nullptr && \"Error during science: input grid "
                        "sp is null\"",
                        "variant_science.in.cpp", 51, __PRETTY_FUNCTION__);
    ss != nullptr && "Error during science: input grid ss is null"
        ? (static_cast<void>(0))
        : __assert_fail("ss != nullptr && \"Error during science: input grid "
                        "ss is null\"",
                        "variant_science.in.cpp", 52, __PRETTY_FUNCTION__);
    z_mult_dat != nullptr &&
            "Error during science: input grid z_mult_dat is null"
        ? (static_cast<void>(0))
        : __assert_fail("z_mult_dat != nullptr && \"Error during science: "
                        "input grid z_mult_dat is null\"",
                        "variant_science.in.cpp", 53, __PRETTY_FUNCTION__);
    x_ssl_dat != nullptr && "Error during science: input grid x_ssl_dat is null"
        ? (static_cast<void>(0))
        : __assert_fail("x_ssl_dat != nullptr && \"Error during science: input "
                        "grid x_ssl_dat is null\"",
                        "variant_science.in.cpp", 54, __PRETTY_FUNCTION__);
    y_ssl_dat != nullptr && "Error during science: input grid y_ssl_dat is null"
        ? (static_cast<void>(0))
        : __assert_fail("y_ssl_dat != nullptr && \"Error during science: input "
                        "grid y_ssl_dat is null\"",
                        "variant_science.in.cpp", 55, __PRETTY_FUNCTION__);
  }
  // Do baseline scientific kernel
  int x;
  int y;
  int z;
  int nx = domain->nx;
  int ny = domain->ny;
  int nz = domain->nz;
  {
    if (nx >= 1 && ny >= 1 && nz >= 1) {
      {
        for (int omplc_gen_iter_2 = 0; omplc_gen_iter_2 < ny;
             omplc_gen_iter_2 = omplc_gen_iter_2 + 1) {
          for (int omplc_gen_iter_3 = 0; omplc_gen_iter_3 < nz;
               omplc_gen_iter_3 = omplc_gen_iter_3 + 1) {
            {{int z = omplc_gen_iter_3;
            int y = omplc_gen_iter_2;
            int x = 0;
            {
              fp->data[x + fp->domain->nx * y +
                       fp->domain->nx * fp->domain->ny * z] =
                  (sp->data[x + sp->domain->nx * y +
                            sp->domain->nx * sp->domain->ny * z] *
                       dp->data[x + dp->domain->nx * y +
                                dp->domain->nx * dp->domain->ny * z] -
                   osp->data[x + osp->domain->nx * y +
                             osp->domain->nx * osp->domain->ny * z] *
                       odp->data[x + odp->domain->nx * y +
                                 odp->domain->nx * odp->domain->ny * z]) *
                  pop->data[x + pop->domain->nx * y +
                            pop->domain->nx * pop->domain->ny * z] *
                  z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                   z_mult_dat->domain->nx *
                                       z_mult_dat->domain->ny * z];
            };
          }
          {
            int z = omplc_gen_iter_3;
            int y = omplc_gen_iter_2;
            int x = 0;
            {
              fp->data[x + fp->domain->nx * y +
                       fp->domain->nx * fp->domain->ny * z] +=
                  ss->data[x + ss->domain->nx * y +
                           ss->domain->nx * ss->domain->ny * z] *
                  z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                   z_mult_dat->domain->nx *
                                       z_mult_dat->domain->ny * z] *
                  (pp->data[x + pp->domain->nx * y +
                            pp->domain->nx * pp->domain->ny * z] *
                       sp->data[x + sp->domain->nx * y +
                                sp->domain->nx * sp->domain->ny * z] *
                       dp->data[x + dp->domain->nx * y +
                                dp->domain->nx * dp->domain->ny * z] -
                   opp->data[x + opp->domain->nx * y +
                             opp->domain->nx * opp->domain->ny * z] *
                       osp->data[x + osp->domain->nx * y +
                                 osp->domain->nx * osp->domain->ny * z] *
                       odp->data[x + odp->domain->nx * y +
                                 odp->domain->nx * odp->domain->ny * z]);
            };
          }
        }
        {
          int z = omplc_gen_iter_3;
          int y = omplc_gen_iter_2;
          int x = 0;
          {
            fp->data[x + fp->domain->nx * y +
                     fp->domain->nx * fp->domain->ny * z] -=
                z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                 z_mult_dat->domain->nx *
                                     z_mult_dat->domain->ny * z] *
                (sp->data[x + sp->domain->nx * y +
                          sp->domain->nx * sp->domain->ny * z] *
                 et->data[x + et->domain->nx * y +
                          et->domain->nx * et->domain->ny * z]);
          };
        }
      }
    }
  }
  {
    for (int omplc_gen_iter_1 = 1; omplc_gen_iter_1 < nx;
         omplc_gen_iter_1 = omplc_gen_iter_1 + 1) {
      {
        for (int omplc_gen_iter_3 = 0; omplc_gen_iter_3 < nz;
             omplc_gen_iter_3 = omplc_gen_iter_3 + 1) {
          {{int z = omplc_gen_iter_3;
          int y = 0;
          int x = omplc_gen_iter_1;
          {
            fp->data[x + fp->domain->nx * y +
                     fp->domain->nx * fp->domain->ny * z] =
                (sp->data[x + sp->domain->nx * y +
                          sp->domain->nx * sp->domain->ny * z] *
                     dp->data[x + dp->domain->nx * y +
                              dp->domain->nx * dp->domain->ny * z] -
                 osp->data[x + osp->domain->nx * y +
                           osp->domain->nx * osp->domain->ny * z] *
                     odp->data[x + odp->domain->nx * y +
                               odp->domain->nx * odp->domain->ny * z]) *
                pop->data[x + pop->domain->nx * y +
                          pop->domain->nx * pop->domain->ny * z] *
                z_mult_dat
                    ->data[x + z_mult_dat->domain->nx * y +
                           z_mult_dat->domain->nx * z_mult_dat->domain->ny * z];
          };
        }
        {
          int z = omplc_gen_iter_3;
          int y = 0;
          int x = omplc_gen_iter_1;
          {
            fp->data[x + fp->domain->nx * y +
                     fp->domain->nx * fp->domain->ny * z] +=
                ss->data[x + ss->domain->nx * y +
                         ss->domain->nx * ss->domain->ny * z] *
                z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                 z_mult_dat->domain->nx *
                                     z_mult_dat->domain->ny * z] *
                (pp->data[x + pp->domain->nx * y +
                          pp->domain->nx * pp->domain->ny * z] *
                     sp->data[x + sp->domain->nx * y +
                              sp->domain->nx * sp->domain->ny * z] *
                     dp->data[x + dp->domain->nx * y +
                              dp->domain->nx * dp->domain->ny * z] -
                 opp->data[x + opp->domain->nx * y +
                           opp->domain->nx * opp->domain->ny * z] *
                     osp->data[x + osp->domain->nx * y +
                               osp->domain->nx * osp->domain->ny * z] *
                     odp->data[x + odp->domain->nx * y +
                               odp->domain->nx * odp->domain->ny * z]);
          };
        }
      }
      {
        int z = omplc_gen_iter_3;
        int y = 0;
        int x = omplc_gen_iter_1;
        {
          fp->data[x + fp->domain->nx * y +
                   fp->domain->nx * fp->domain->ny * z] -=
              z_mult_dat
                  ->data[x + z_mult_dat->domain->nx * y +
                         z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] *
              (sp->data[x + sp->domain->nx * y +
                        sp->domain->nx * sp->domain->ny * z] *
               et->data[x + et->domain->nx * y +
                        et->domain->nx * et->domain->ny * z]);
        };
      }
    }
  }
  {
    for (int omplc_gen_iter_2 = 1; omplc_gen_iter_2 < ny;
         omplc_gen_iter_2 = omplc_gen_iter_2 + 1) {
      {{{int z = 0;
      int y = omplc_gen_iter_2;
      int x = omplc_gen_iter_1;
      {
        fp->data[x + fp->domain->nx * y + fp->domain->nx * fp->domain->ny * z] =
            (sp->data[x + sp->domain->nx * y +
                      sp->domain->nx * sp->domain->ny * z] *
                 dp->data[x + dp->domain->nx * y +
                          dp->domain->nx * dp->domain->ny * z] -
             osp->data[x + osp->domain->nx * y +
                       osp->domain->nx * osp->domain->ny * z] *
                 odp->data[x + odp->domain->nx * y +
                           odp->domain->nx * odp->domain->ny * z]) *
            pop->data[x + pop->domain->nx * y +
                      pop->domain->nx * pop->domain->ny * z] *
            z_mult_dat
                ->data[x + z_mult_dat->domain->nx * y +
                       z_mult_dat->domain->nx * z_mult_dat->domain->ny * z];
      };
    }
    {
      int z = 0;
      int y = omplc_gen_iter_2;
      int x = omplc_gen_iter_1;
      {
        fp->data[x + fp->domain->nx * y +
                 fp->domain->nx * fp->domain->ny * z] +=
            ss->data[x + ss->domain->nx * y +
                     ss->domain->nx * ss->domain->ny * z] *
            z_mult_dat
                ->data[x + z_mult_dat->domain->nx * y +
                       z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] *
            (pp->data[x + pp->domain->nx * y +
                      pp->domain->nx * pp->domain->ny * z] *
                 sp->data[x + sp->domain->nx * y +
                          sp->domain->nx * sp->domain->ny * z] *
                 dp->data[x + dp->domain->nx * y +
                          dp->domain->nx * dp->domain->ny * z] -
             opp->data[x + opp->domain->nx * y +
                       opp->domain->nx * opp->domain->ny * z] *
                 osp->data[x + osp->domain->nx * y +
                           osp->domain->nx * osp->domain->ny * z] *
                 odp->data[x + odp->domain->nx * y +
                           odp->domain->nx * odp->domain->ny * z]);
      };
    }
  }
  {
    int z = 0;
    int y = omplc_gen_iter_2;
    int x = omplc_gen_iter_1;
    {
      fp->data[x + fp->domain->nx * y + fp->domain->nx * fp->domain->ny * z] -=
          z_mult_dat
              ->data[x + z_mult_dat->domain->nx * y +
                     z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] *
          (sp->data[x + sp->domain->nx * y +
                    sp->domain->nx * sp->domain->ny * z] *
           et->data[x + et->domain->nx * y +
                    et->domain->nx * et->domain->ny * z]);
    };
  }
}
{
  for (int omplc_gen_iter_3 = 1; omplc_gen_iter_3 < nz;
       omplc_gen_iter_3 = omplc_gen_iter_3 + 1) {
    {{{int z = omplc_gen_iter_3;
    int y = omplc_gen_iter_2;
    int x = omplc_gen_iter_1;
    {
      fp->data[x + fp->domain->nx * y + fp->domain->nx * fp->domain->ny * z] =
          (sp->data[x + sp->domain->nx * y +
                    sp->domain->nx * sp->domain->ny * z] *
               dp->data[x + dp->domain->nx * y +
                        dp->domain->nx * dp->domain->ny * z] -
           osp->data[x + osp->domain->nx * y +
                     osp->domain->nx * osp->domain->ny * z] *
               odp->data[x + odp->domain->nx * y +
                         odp->domain->nx * odp->domain->ny * z]) *
          pop->data[x + pop->domain->nx * y +
                    pop->domain->nx * pop->domain->ny * z] *
          z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                           z_mult_dat->domain->nx * z_mult_dat->domain->ny * z];
    };
  }
  {
    int z = omplc_gen_iter_3;
    int y = omplc_gen_iter_2;
    int x = omplc_gen_iter_1;
    {
      fp->data[x + fp->domain->nx * y + fp->domain->nx * fp->domain->ny * z] +=
          ss->data[x + ss->domain->nx * y +
                   ss->domain->nx * ss->domain->ny * z] *
          z_mult_dat
              ->data[x + z_mult_dat->domain->nx * y +
                     z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] *
          (pp->data[x + pp->domain->nx * y +
                    pp->domain->nx * pp->domain->ny * z] *
               sp->data[x + sp->domain->nx * y +
                        sp->domain->nx * sp->domain->ny * z] *
               dp->data[x + dp->domain->nx * y +
                        dp->domain->nx * dp->domain->ny * z] -
           opp->data[x + opp->domain->nx * y +
                     opp->domain->nx * opp->domain->ny * z] *
               osp->data[x + osp->domain->nx * y +
                         osp->domain->nx * osp->domain->ny * z] *
               odp->data[x + odp->domain->nx * y +
                         odp->domain->nx * odp->domain->ny * z]);
    };
  }
}
{
  int z = omplc_gen_iter_3;
  int y = omplc_gen_iter_2;
  int x = omplc_gen_iter_1;
  {
    fp->data[x + fp->domain->nx * y + fp->domain->nx * fp->domain->ny * z] -=
        z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                         z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] *
        (sp->data[x + sp->domain->nx * y +
                  sp->domain->nx * sp->domain->ny * z] *
         et->data[x + et->domain->nx * y +
                  et->domain->nx * et->domain->ny * z]);
  };
}
}
{
  int z = omplc_gen_iter_3 - 1;
  int y = omplc_gen_iter_2 - 1;
  int x = omplc_gen_iter_1 - 1;
  {
    double x_dir_g =
        0.5 *
        (x_ssl_dat->data[x + x_ssl_dat->domain->nx * y +
                         x_ssl_dat->domain->nx * x_ssl_dat->domain->ny * 0] +
         x_ssl_dat->data[x + 1 + x_ssl_dat->domain->nx * y +
                         x_ssl_dat->domain->nx * x_ssl_dat->domain->ny * 0]);
    double x_dir_g_c =
        0.5 *
        (x_ssl_dat->data[x + x_ssl_dat->domain->nx * y +
                         x_ssl_dat->domain->nx * x_ssl_dat->domain->ny * 0] +
         x_ssl_dat->data[x + 1 + x_ssl_dat->domain->nx * y +
                         x_ssl_dat->domain->nx * x_ssl_dat->domain->ny * 0]);
    double y_dir_g =
        0.5 *
        (y_ssl_dat->data[x + y_ssl_dat->domain->nx * y +
                         y_ssl_dat->domain->nx * y_ssl_dat->domain->ny * 0] +
         y_ssl_dat->data[x + y_ssl_dat->domain->nx * (y + 1) +
                         y_ssl_dat->domain->nx * y_ssl_dat->domain->ny * 0]);
    double y_dir_g_c =
        0.5 *
        (y_ssl_dat->data[x + y_ssl_dat->domain->nx * y +
                         y_ssl_dat->domain->nx * y_ssl_dat->domain->ny * 0] +
         y_ssl_dat->data[x + y_ssl_dat->domain->nx * (y + 1) +
                         y_ssl_dat->domain->nx * y_ssl_dat->domain->ny * 0]);
    double diff_right =
        pp->data[x + pp->domain->nx * y + pp->domain->nx * pp->domain->ny * z] -
        pp->data[x + 1 + pp->domain->nx * y +
                 pp->domain->nx * pp->domain->ny * z];
    double diff_front =
        pp->data[x + pp->domain->nx * y + pp->domain->nx * pp->domain->ny * z] -
        pp->data[x + pp->domain->nx * (y + 1) +
                 pp->domain->nx * pp->domain->ny * z];
    double updir_right = diff_right * x_dir_g_c - x_dir_g;
    double updir_front = diff_front * y_dir_g_c - y_dir_g;
    double sep =
        0.5 *
        (z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                          z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] +
         z_mult_dat
             ->data[x + z_mult_dat->domain->nx * y +
                    z_mult_dat->domain->nx * z_mult_dat->domain->ny * (z + 1)]);
    double lower_cond =
        pp->data[x + pp->domain->nx * y + pp->domain->nx * pp->domain->ny * z] /
            sep -
        z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                         z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] /
            (z_mult_dat
                 ->data[x + z_mult_dat->domain->nx * y +
                        z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] +
             z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                              z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                                  (z + 1)]) *
            dp->data[x + dp->domain->nx * y +
                     dp->domain->nx * dp->domain->ny * z];
    double upper_cond =
        pp->data[x + pp->domain->nx * y +
                 pp->domain->nx * pp->domain->ny * (z + 1)] /
            sep +
        z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                         z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                             (z + 1)] /
            (z_mult_dat
                 ->data[x + z_mult_dat->domain->nx * y +
                        z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] +
             z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                              z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                                  (z + 1)]) *
            dp->data[x + dp->domain->nx * y +
                     dp->domain->nx * dp->domain->ny * (z + 1)];
    double diff_upper = lower_cond - upper_cond;
    double u_right =
        z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                         z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] *
            ((((bool)(permxp
                          ->data[x + permxp->domain->nx * y +
                                 permxp->domain->nx * permxp->domain->ny * z] +
                      permxp
                          ->data[x + 1 + permxp->domain->nx * y +
                                 permxp->domain->nx * permxp->domain->ny * z]))
                  ? 2.0 *
                        permxp->data[x + permxp->domain->nx * y +
                                     permxp->domain->nx * permxp->domain->ny *
                                         z] *
                        permxp->data[x + 1 + permxp->domain->nx * y +
                                     permxp->domain->nx * permxp->domain->ny *
                                         z] /
                        (permxp->data[x + permxp->domain->nx * y +
                                      permxp->domain->nx * permxp->domain->ny *
                                          z] +
                         permxp->data[x + 1 + permxp->domain->nx * y +
                                      permxp->domain->nx * permxp->domain->ny *
                                          z])
                  : ((double)0))) *
            diff_right * x_dir_g_c *
            ((updir_right - 0.0 >= ((double)0)
                  ? rpp->data[x + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * z] *
                        dp->data[x + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * z]
                  : rpp->data[x + 1 + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * z] *
                        dp->data[x + 1 + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * z])) +
        z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                         z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] *
            ((((bool)(permxp
                          ->data[x + permxp->domain->nx * y +
                                 permxp->domain->nx * permxp->domain->ny * z] +
                      permxp
                          ->data[x + 1 + permxp->domain->nx * y +
                                 permxp->domain->nx * permxp->domain->ny * z]))
                  ? 2.0 *
                        permxp->data[x + permxp->domain->nx * y +
                                     permxp->domain->nx * permxp->domain->ny *
                                         z] *
                        permxp->data[x + 1 + permxp->domain->nx * y +
                                     permxp->domain->nx * permxp->domain->ny *
                                         z] /
                        (permxp->data[x + permxp->domain->nx * y +
                                      permxp->domain->nx * permxp->domain->ny *
                                          z] +
                         permxp->data[x + 1 + permxp->domain->nx * y +
                                      permxp->domain->nx * permxp->domain->ny *
                                          z])
                  : ((double)0))) *
            -x_dir_g *
            ((updir_right - 0.0 >= ((double)0)
                  ? rpp->data[x + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * z] *
                        dp->data[x + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * z]
                  : rpp->data[x + 1 + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * z] *
                        dp->data[x + 1 + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * z]));
    double u_front =
        z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                         z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] *
            ((((bool)(permyp
                          ->data[x + permyp->domain->nx * y +
                                 permyp->domain->nx * permyp->domain->ny * z] +
                      permyp
                          ->data[x + 1 + permyp->domain->nx * y +
                                 permyp->domain->nx * permyp->domain->ny * z]))
                  ? 2.0 *
                        permyp->data[x + permyp->domain->nx * y +
                                     permyp->domain->nx * permyp->domain->ny *
                                         z] *
                        permyp->data[x + 1 + permyp->domain->nx * y +
                                     permyp->domain->nx * permyp->domain->ny *
                                         z] /
                        (permyp->data[x + permyp->domain->nx * y +
                                      permyp->domain->nx * permyp->domain->ny *
                                          z] +
                         permyp->data[x + 1 + permyp->domain->nx * y +
                                      permyp->domain->nx * permyp->domain->ny *
                                          z])
                  : ((double)0))) *
            diff_front * x_dir_g_c *
            ((updir_front - 0.0 >= ((double)0)
                  ? rpp->data[x + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * z] *
                        dp->data[x + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * z]
                  : rpp->data[x + 1 + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * z] *
                        dp->data[x + 1 + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * z])) +
        z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                         z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] *
            ((((bool)(permyp
                          ->data[x + permyp->domain->nx * y +
                                 permyp->domain->nx * permyp->domain->ny * z] +
                      permyp
                          ->data[x + 1 + permyp->domain->nx * y +
                                 permyp->domain->nx * permyp->domain->ny * z]))
                  ? 2.0 *
                        permyp->data[x + permyp->domain->nx * y +
                                     permyp->domain->nx * permyp->domain->ny *
                                         z] *
                        permyp->data[x + 1 + permyp->domain->nx * y +
                                     permyp->domain->nx * permyp->domain->ny *
                                         z] /
                        (permyp->data[x + permyp->domain->nx * y +
                                      permyp->domain->nx * permyp->domain->ny *
                                          z] +
                         permyp->data[x + 1 + permyp->domain->nx * y +
                                      permyp->domain->nx * permyp->domain->ny *
                                          z])
                  : ((double)0))) *
            -x_dir_g *
            ((updir_front - 0.0 >= ((double)0)
                  ? rpp->data[x + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * z] *
                        dp->data[x + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * z]
                  : rpp->data[x + 1 + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * z] *
                        dp->data[x + 1 + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * z]));
    double u_upper =
        ((((bool)(z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                   z_mult_dat->domain->nx *
                                       z_mult_dat->domain->ny * z] *
                      permzp->data[x + permzp->domain->nx * y +
                                   permzp->domain->nx * permzp->domain->ny *
                                       (z + 1)] +
                  permzp->data[x + permzp->domain->nx * y +
                               permzp->domain->nx * permzp->domain->ny * z] *
                      z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                       z_mult_dat->domain->nx *
                                           z_mult_dat->domain->ny * (z + 1)]))
              ? (z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                  z_mult_dat->domain->nx *
                                      z_mult_dat->domain->ny * z] +
                 z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                  z_mult_dat->domain->nx *
                                      z_mult_dat->domain->ny * (z + 1)]) *
                    permzp->data[x + permzp->domain->nx * y +
                                 permzp->domain->nx * permzp->domain->ny * z] *
                    permzp->data[x + permzp->domain->nx * y +
                                 permzp->domain->nx * permzp->domain->ny *
                                     (z + 1)] /
                    (permzp->data[x + permzp->domain->nx * y +
                                  permzp->domain->nx * permzp->domain->ny *
                                      (z + 1)] *
                         z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                          z_mult_dat->domain->nx *
                                              z_mult_dat->domain->ny * z] +
                     permzp->data[x + permzp->domain->nx * y +
                                  permzp->domain->nx * permzp->domain->ny * z] *
                         z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                          z_mult_dat->domain->nx *
                                              z_mult_dat->domain->ny * (z + 1)])
              : ((double)0))) *
        diff_upper *
        ((lower_cond - upper_cond >= ((double)0)
              ? rpp->data[x + rpp->domain->nx * y +
                          rpp->domain->nx * rpp->domain->ny * z] *
                    dp->data[x + dp->domain->nx * y +
                             dp->domain->nx * dp->domain->ny * z]
              : rpp->data[x + rpp->domain->nx * y +
                          rpp->domain->nx * rpp->domain->ny * (z + 1)] *
                    dp->data[x + dp->domain->nx * y +
                             dp->domain->nx * dp->domain->ny * (z + 1)]));
    vx->data[x + vx->domain->nx * y + vx->domain->nx * vx->domain->ny * z] =
        u_right;
    vy->data[x + vy->domain->nx * y + vy->domain->nx * vy->domain->ny * z] =
        u_front;
    vz->data[x + vz->domain->nx * y + vz->domain->nx * vz->domain->ny * z] =
        u_upper;
    fp->data[x + fp->domain->nx * y + fp->domain->nx * fp->domain->ny * z] +=
        u_right * u_front * u_upper;
    fp->data[x + 1 + fp->domain->nx * y +
             fp->domain->nx * fp->domain->ny * z] += u_right;
    fp->data[x + fp->domain->nx * (y + 1) +
             fp->domain->nx * fp->domain->ny * z] -= u_front;
    fp->data[x + fp->domain->nx * y +
             fp->domain->nx * fp->domain->ny * (z + 1)] -= u_upper;
  };
}
}
}
{
  int z = nz - 1;
  int y = omplc_gen_iter_2 - 1;
  int x = omplc_gen_iter_1 - 1;
  {
    double x_dir_g =
        0.5 *
        (x_ssl_dat->data[x + x_ssl_dat->domain->nx * y +
                         x_ssl_dat->domain->nx * x_ssl_dat->domain->ny * 0] +
         x_ssl_dat->data[x + 1 + x_ssl_dat->domain->nx * y +
                         x_ssl_dat->domain->nx * x_ssl_dat->domain->ny * 0]);
    double x_dir_g_c =
        0.5 *
        (x_ssl_dat->data[x + x_ssl_dat->domain->nx * y +
                         x_ssl_dat->domain->nx * x_ssl_dat->domain->ny * 0] +
         x_ssl_dat->data[x + 1 + x_ssl_dat->domain->nx * y +
                         x_ssl_dat->domain->nx * x_ssl_dat->domain->ny * 0]);
    double y_dir_g =
        0.5 *
        (y_ssl_dat->data[x + y_ssl_dat->domain->nx * y +
                         y_ssl_dat->domain->nx * y_ssl_dat->domain->ny * 0] +
         y_ssl_dat->data[x + y_ssl_dat->domain->nx * (y + 1) +
                         y_ssl_dat->domain->nx * y_ssl_dat->domain->ny * 0]);
    double y_dir_g_c =
        0.5 *
        (y_ssl_dat->data[x + y_ssl_dat->domain->nx * y +
                         y_ssl_dat->domain->nx * y_ssl_dat->domain->ny * 0] +
         y_ssl_dat->data[x + y_ssl_dat->domain->nx * (y + 1) +
                         y_ssl_dat->domain->nx * y_ssl_dat->domain->ny * 0]);
    double diff_right =
        pp->data[x + pp->domain->nx * y + pp->domain->nx * pp->domain->ny * z] -
        pp->data[x + 1 + pp->domain->nx * y +
                 pp->domain->nx * pp->domain->ny * z];
    double diff_front =
        pp->data[x + pp->domain->nx * y + pp->domain->nx * pp->domain->ny * z] -
        pp->data[x + pp->domain->nx * (y + 1) +
                 pp->domain->nx * pp->domain->ny * z];
    double updir_right = diff_right * x_dir_g_c - x_dir_g;
    double updir_front = diff_front * y_dir_g_c - y_dir_g;
    double sep =
        0.5 *
        (z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                          z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] +
         z_mult_dat
             ->data[x + z_mult_dat->domain->nx * y +
                    z_mult_dat->domain->nx * z_mult_dat->domain->ny * (z + 1)]);
    double lower_cond =
        pp->data[x + pp->domain->nx * y + pp->domain->nx * pp->domain->ny * z] /
            sep -
        z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                         z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] /
            (z_mult_dat
                 ->data[x + z_mult_dat->domain->nx * y +
                        z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] +
             z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                              z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                                  (z + 1)]) *
            dp->data[x + dp->domain->nx * y +
                     dp->domain->nx * dp->domain->ny * z];
    double upper_cond =
        pp->data[x + pp->domain->nx * y +
                 pp->domain->nx * pp->domain->ny * (z + 1)] /
            sep +
        z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                         z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                             (z + 1)] /
            (z_mult_dat
                 ->data[x + z_mult_dat->domain->nx * y +
                        z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] +
             z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                              z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                                  (z + 1)]) *
            dp->data[x + dp->domain->nx * y +
                     dp->domain->nx * dp->domain->ny * (z + 1)];
    double diff_upper = lower_cond - upper_cond;
    double u_right =
        z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                         z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] *
            ((((bool)(permxp
                          ->data[x + permxp->domain->nx * y +
                                 permxp->domain->nx * permxp->domain->ny * z] +
                      permxp
                          ->data[x + 1 + permxp->domain->nx * y +
                                 permxp->domain->nx * permxp->domain->ny * z]))
                  ? 2.0 *
                        permxp->data[x + permxp->domain->nx * y +
                                     permxp->domain->nx * permxp->domain->ny *
                                         z] *
                        permxp->data[x + 1 + permxp->domain->nx * y +
                                     permxp->domain->nx * permxp->domain->ny *
                                         z] /
                        (permxp->data[x + permxp->domain->nx * y +
                                      permxp->domain->nx * permxp->domain->ny *
                                          z] +
                         permxp->data[x + 1 + permxp->domain->nx * y +
                                      permxp->domain->nx * permxp->domain->ny *
                                          z])
                  : ((double)0))) *
            diff_right * x_dir_g_c *
            ((updir_right - 0.0 >= ((double)0)
                  ? rpp->data[x + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * z] *
                        dp->data[x + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * z]
                  : rpp->data[x + 1 + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * z] *
                        dp->data[x + 1 + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * z])) +
        z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                         z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] *
            ((((bool)(permxp
                          ->data[x + permxp->domain->nx * y +
                                 permxp->domain->nx * permxp->domain->ny * z] +
                      permxp
                          ->data[x + 1 + permxp->domain->nx * y +
                                 permxp->domain->nx * permxp->domain->ny * z]))
                  ? 2.0 *
                        permxp->data[x + permxp->domain->nx * y +
                                     permxp->domain->nx * permxp->domain->ny *
                                         z] *
                        permxp->data[x + 1 + permxp->domain->nx * y +
                                     permxp->domain->nx * permxp->domain->ny *
                                         z] /
                        (permxp->data[x + permxp->domain->nx * y +
                                      permxp->domain->nx * permxp->domain->ny *
                                          z] +
                         permxp->data[x + 1 + permxp->domain->nx * y +
                                      permxp->domain->nx * permxp->domain->ny *
                                          z])
                  : ((double)0))) *
            -x_dir_g *
            ((updir_right - 0.0 >= ((double)0)
                  ? rpp->data[x + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * z] *
                        dp->data[x + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * z]
                  : rpp->data[x + 1 + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * z] *
                        dp->data[x + 1 + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * z]));
    double u_front =
        z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                         z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] *
            ((((bool)(permyp
                          ->data[x + permyp->domain->nx * y +
                                 permyp->domain->nx * permyp->domain->ny * z] +
                      permyp
                          ->data[x + 1 + permyp->domain->nx * y +
                                 permyp->domain->nx * permyp->domain->ny * z]))
                  ? 2.0 *
                        permyp->data[x + permyp->domain->nx * y +
                                     permyp->domain->nx * permyp->domain->ny *
                                         z] *
                        permyp->data[x + 1 + permyp->domain->nx * y +
                                     permyp->domain->nx * permyp->domain->ny *
                                         z] /
                        (permyp->data[x + permyp->domain->nx * y +
                                      permyp->domain->nx * permyp->domain->ny *
                                          z] +
                         permyp->data[x + 1 + permyp->domain->nx * y +
                                      permyp->domain->nx * permyp->domain->ny *
                                          z])
                  : ((double)0))) *
            diff_front * x_dir_g_c *
            ((updir_front - 0.0 >= ((double)0)
                  ? rpp->data[x + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * z] *
                        dp->data[x + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * z]
                  : rpp->data[x + 1 + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * z] *
                        dp->data[x + 1 + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * z])) +
        z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                         z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] *
            ((((bool)(permyp
                          ->data[x + permyp->domain->nx * y +
                                 permyp->domain->nx * permyp->domain->ny * z] +
                      permyp
                          ->data[x + 1 + permyp->domain->nx * y +
                                 permyp->domain->nx * permyp->domain->ny * z]))
                  ? 2.0 *
                        permyp->data[x + permyp->domain->nx * y +
                                     permyp->domain->nx * permyp->domain->ny *
                                         z] *
                        permyp->data[x + 1 + permyp->domain->nx * y +
                                     permyp->domain->nx * permyp->domain->ny *
                                         z] /
                        (permyp->data[x + permyp->domain->nx * y +
                                      permyp->domain->nx * permyp->domain->ny *
                                          z] +
                         permyp->data[x + 1 + permyp->domain->nx * y +
                                      permyp->domain->nx * permyp->domain->ny *
                                          z])
                  : ((double)0))) *
            -x_dir_g *
            ((updir_front - 0.0 >= ((double)0)
                  ? rpp->data[x + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * z] *
                        dp->data[x + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * z]
                  : rpp->data[x + 1 + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * z] *
                        dp->data[x + 1 + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * z]));
    double u_upper =
        ((((bool)(z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                   z_mult_dat->domain->nx *
                                       z_mult_dat->domain->ny * z] *
                      permzp->data[x + permzp->domain->nx * y +
                                   permzp->domain->nx * permzp->domain->ny *
                                       (z + 1)] +
                  permzp->data[x + permzp->domain->nx * y +
                               permzp->domain->nx * permzp->domain->ny * z] *
                      z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                       z_mult_dat->domain->nx *
                                           z_mult_dat->domain->ny * (z + 1)]))
              ? (z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                  z_mult_dat->domain->nx *
                                      z_mult_dat->domain->ny * z] +
                 z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                  z_mult_dat->domain->nx *
                                      z_mult_dat->domain->ny * (z + 1)]) *
                    permzp->data[x + permzp->domain->nx * y +
                                 permzp->domain->nx * permzp->domain->ny * z] *
                    permzp->data[x + permzp->domain->nx * y +
                                 permzp->domain->nx * permzp->domain->ny *
                                     (z + 1)] /
                    (permzp->data[x + permzp->domain->nx * y +
                                  permzp->domain->nx * permzp->domain->ny *
                                      (z + 1)] *
                         z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                          z_mult_dat->domain->nx *
                                              z_mult_dat->domain->ny * z] +
                     permzp->data[x + permzp->domain->nx * y +
                                  permzp->domain->nx * permzp->domain->ny * z] *
                         z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                          z_mult_dat->domain->nx *
                                              z_mult_dat->domain->ny * (z + 1)])
              : ((double)0))) *
        diff_upper *
        ((lower_cond - upper_cond >= ((double)0)
              ? rpp->data[x + rpp->domain->nx * y +
                          rpp->domain->nx * rpp->domain->ny * z] *
                    dp->data[x + dp->domain->nx * y +
                             dp->domain->nx * dp->domain->ny * z]
              : rpp->data[x + rpp->domain->nx * y +
                          rpp->domain->nx * rpp->domain->ny * (z + 1)] *
                    dp->data[x + dp->domain->nx * y +
                             dp->domain->nx * dp->domain->ny * (z + 1)]));
    vx->data[x + vx->domain->nx * y + vx->domain->nx * vx->domain->ny * z] =
        u_right;
    vy->data[x + vy->domain->nx * y + vy->domain->nx * vy->domain->ny * z] =
        u_front;
    vz->data[x + vz->domain->nx * y + vz->domain->nx * vz->domain->ny * z] =
        u_upper;
    fp->data[x + fp->domain->nx * y + fp->domain->nx * fp->domain->ny * z] +=
        u_right * u_front * u_upper;
    fp->data[x + 1 + fp->domain->nx * y +
             fp->domain->nx * fp->domain->ny * z] += u_right;
    fp->data[x + fp->domain->nx * (y + 1) +
             fp->domain->nx * fp->domain->ny * z] -= u_front;
    fp->data[x + fp->domain->nx * y +
             fp->domain->nx * fp->domain->ny * (z + 1)] -= u_upper;
  };
}
}
}
{
  for (int omplc_gen_iter_3 = 1; omplc_gen_iter_3 <= nz;
       omplc_gen_iter_3 = omplc_gen_iter_3 + 1) {
    {
      int z = omplc_gen_iter_3 - 1;
      int y = ny - 1;
      int x = omplc_gen_iter_1 - 1;
      {
        double x_dir_g =
            0.5 *
            (x_ssl_dat
                 ->data[x + x_ssl_dat->domain->nx * y +
                        x_ssl_dat->domain->nx * x_ssl_dat->domain->ny * 0] +
             x_ssl_dat
                 ->data[x + 1 + x_ssl_dat->domain->nx * y +
                        x_ssl_dat->domain->nx * x_ssl_dat->domain->ny * 0]);
        double x_dir_g_c =
            0.5 *
            (x_ssl_dat
                 ->data[x + x_ssl_dat->domain->nx * y +
                        x_ssl_dat->domain->nx * x_ssl_dat->domain->ny * 0] +
             x_ssl_dat
                 ->data[x + 1 + x_ssl_dat->domain->nx * y +
                        x_ssl_dat->domain->nx * x_ssl_dat->domain->ny * 0]);
        double y_dir_g =
            0.5 *
            (y_ssl_dat
                 ->data[x + y_ssl_dat->domain->nx * y +
                        y_ssl_dat->domain->nx * y_ssl_dat->domain->ny * 0] +
             y_ssl_dat
                 ->data[x + y_ssl_dat->domain->nx * (y + 1) +
                        y_ssl_dat->domain->nx * y_ssl_dat->domain->ny * 0]);
        double y_dir_g_c =
            0.5 *
            (y_ssl_dat
                 ->data[x + y_ssl_dat->domain->nx * y +
                        y_ssl_dat->domain->nx * y_ssl_dat->domain->ny * 0] +
             y_ssl_dat
                 ->data[x + y_ssl_dat->domain->nx * (y + 1) +
                        y_ssl_dat->domain->nx * y_ssl_dat->domain->ny * 0]);
        double diff_right = pp->data[x + pp->domain->nx * y +
                                     pp->domain->nx * pp->domain->ny * z] -
                            pp->data[x + 1 + pp->domain->nx * y +
                                     pp->domain->nx * pp->domain->ny * z];
        double diff_front = pp->data[x + pp->domain->nx * y +
                                     pp->domain->nx * pp->domain->ny * z] -
                            pp->data[x + pp->domain->nx * (y + 1) +
                                     pp->domain->nx * pp->domain->ny * z];
        double updir_right = diff_right * x_dir_g_c - x_dir_g;
        double updir_front = diff_front * y_dir_g_c - y_dir_g;
        double sep =
            0.5 *
            (z_mult_dat
                 ->data[x + z_mult_dat->domain->nx * y +
                        z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] +
             z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                              z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                                  (z + 1)]);
        double lower_cond =
            pp->data[x + pp->domain->nx * y +
                     pp->domain->nx * pp->domain->ny * z] /
                sep -
            z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                             z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                                 z] /
                (z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                  z_mult_dat->domain->nx *
                                      z_mult_dat->domain->ny * z] +
                 z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                  z_mult_dat->domain->nx *
                                      z_mult_dat->domain->ny * (z + 1)]) *
                dp->data[x + dp->domain->nx * y +
                         dp->domain->nx * dp->domain->ny * z];
        double upper_cond =
            pp->data[x + pp->domain->nx * y +
                     pp->domain->nx * pp->domain->ny * (z + 1)] /
                sep +
            z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                             z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                                 (z + 1)] /
                (z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                  z_mult_dat->domain->nx *
                                      z_mult_dat->domain->ny * z] +
                 z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                  z_mult_dat->domain->nx *
                                      z_mult_dat->domain->ny * (z + 1)]) *
                dp->data[x + dp->domain->nx * y +
                         dp->domain->nx * dp->domain->ny * (z + 1)];
        double diff_upper = lower_cond - upper_cond;
        double u_right =
            z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                             z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                                 z] *
                ((((bool)(permxp->data[x + permxp->domain->nx * y +
                                       permxp->domain->nx * permxp->domain->ny *
                                           z] +
                          permxp->data[x + 1 + permxp->domain->nx * y +
                                       permxp->domain->nx * permxp->domain->ny *
                                           z]))
                      ? 2.0 *
                            permxp->data[x + permxp->domain->nx * y +
                                         permxp->domain->nx *
                                             permxp->domain->ny * z] *
                            permxp->data[x + 1 + permxp->domain->nx * y +
                                         permxp->domain->nx *
                                             permxp->domain->ny * z] /
                            (permxp->data[x + permxp->domain->nx * y +
                                          permxp->domain->nx *
                                              permxp->domain->ny * z] +
                             permxp->data[x + 1 + permxp->domain->nx * y +
                                          permxp->domain->nx *
                                              permxp->domain->ny * z])
                      : ((double)0))) *
                diff_right * x_dir_g_c *
                ((updir_right - 0.0 >= ((double)0)
                      ? rpp->data[x + rpp->domain->nx * y +
                                  rpp->domain->nx * rpp->domain->ny * z] *
                            dp->data[x + dp->domain->nx * y +
                                     dp->domain->nx * dp->domain->ny * z]
                      : rpp->data[x + 1 + rpp->domain->nx * y +
                                  rpp->domain->nx * rpp->domain->ny * z] *
                            dp->data[x + 1 + dp->domain->nx * y +
                                     dp->domain->nx * dp->domain->ny * z])) +
            z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                             z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                                 z] *
                ((((bool)(permxp->data[x + permxp->domain->nx * y +
                                       permxp->domain->nx * permxp->domain->ny *
                                           z] +
                          permxp->data[x + 1 + permxp->domain->nx * y +
                                       permxp->domain->nx * permxp->domain->ny *
                                           z]))
                      ? 2.0 *
                            permxp->data[x + permxp->domain->nx * y +
                                         permxp->domain->nx *
                                             permxp->domain->ny * z] *
                            permxp->data[x + 1 + permxp->domain->nx * y +
                                         permxp->domain->nx *
                                             permxp->domain->ny * z] /
                            (permxp->data[x + permxp->domain->nx * y +
                                          permxp->domain->nx *
                                              permxp->domain->ny * z] +
                             permxp->data[x + 1 + permxp->domain->nx * y +
                                          permxp->domain->nx *
                                              permxp->domain->ny * z])
                      : ((double)0))) *
                -x_dir_g *
                ((updir_right - 0.0 >= ((double)0)
                      ? rpp->data[x + rpp->domain->nx * y +
                                  rpp->domain->nx * rpp->domain->ny * z] *
                            dp->data[x + dp->domain->nx * y +
                                     dp->domain->nx * dp->domain->ny * z]
                      : rpp->data[x + 1 + rpp->domain->nx * y +
                                  rpp->domain->nx * rpp->domain->ny * z] *
                            dp->data[x + 1 + dp->domain->nx * y +
                                     dp->domain->nx * dp->domain->ny * z]));
        double u_front =
            z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                             z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                                 z] *
                ((((bool)(permyp->data[x + permyp->domain->nx * y +
                                       permyp->domain->nx * permyp->domain->ny *
                                           z] +
                          permyp->data[x + 1 + permyp->domain->nx * y +
                                       permyp->domain->nx * permyp->domain->ny *
                                           z]))
                      ? 2.0 *
                            permyp->data[x + permyp->domain->nx * y +
                                         permyp->domain->nx *
                                             permyp->domain->ny * z] *
                            permyp->data[x + 1 + permyp->domain->nx * y +
                                         permyp->domain->nx *
                                             permyp->domain->ny * z] /
                            (permyp->data[x + permyp->domain->nx * y +
                                          permyp->domain->nx *
                                              permyp->domain->ny * z] +
                             permyp->data[x + 1 + permyp->domain->nx * y +
                                          permyp->domain->nx *
                                              permyp->domain->ny * z])
                      : ((double)0))) *
                diff_front * x_dir_g_c *
                ((updir_front - 0.0 >= ((double)0)
                      ? rpp->data[x + rpp->domain->nx * y +
                                  rpp->domain->nx * rpp->domain->ny * z] *
                            dp->data[x + dp->domain->nx * y +
                                     dp->domain->nx * dp->domain->ny * z]
                      : rpp->data[x + 1 + rpp->domain->nx * y +
                                  rpp->domain->nx * rpp->domain->ny * z] *
                            dp->data[x + 1 + dp->domain->nx * y +
                                     dp->domain->nx * dp->domain->ny * z])) +
            z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                             z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                                 z] *
                ((((bool)(permyp->data[x + permyp->domain->nx * y +
                                       permyp->domain->nx * permyp->domain->ny *
                                           z] +
                          permyp->data[x + 1 + permyp->domain->nx * y +
                                       permyp->domain->nx * permyp->domain->ny *
                                           z]))
                      ? 2.0 *
                            permyp->data[x + permyp->domain->nx * y +
                                         permyp->domain->nx *
                                             permyp->domain->ny * z] *
                            permyp->data[x + 1 + permyp->domain->nx * y +
                                         permyp->domain->nx *
                                             permyp->domain->ny * z] /
                            (permyp->data[x + permyp->domain->nx * y +
                                          permyp->domain->nx *
                                              permyp->domain->ny * z] +
                             permyp->data[x + 1 + permyp->domain->nx * y +
                                          permyp->domain->nx *
                                              permyp->domain->ny * z])
                      : ((double)0))) *
                -x_dir_g *
                ((updir_front - 0.0 >= ((double)0)
                      ? rpp->data[x + rpp->domain->nx * y +
                                  rpp->domain->nx * rpp->domain->ny * z] *
                            dp->data[x + dp->domain->nx * y +
                                     dp->domain->nx * dp->domain->ny * z]
                      : rpp->data[x + 1 + rpp->domain->nx * y +
                                  rpp->domain->nx * rpp->domain->ny * z] *
                            dp->data[x + 1 + dp->domain->nx * y +
                                     dp->domain->nx * dp->domain->ny * z]));
        double u_upper =
            ((((bool)(z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                       z_mult_dat->domain->nx *
                                           z_mult_dat->domain->ny * z] *
                          permzp->data[x + permzp->domain->nx * y +
                                       permzp->domain->nx * permzp->domain->ny *
                                           (z + 1)] +
                      permzp->data[x + permzp->domain->nx * y +
                                   permzp->domain->nx * permzp->domain->ny *
                                       z] *
                          z_mult_dat
                              ->data[x + z_mult_dat->domain->nx * y +
                                     z_mult_dat->domain->nx *
                                         z_mult_dat->domain->ny * (z + 1)]))
                  ? (z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                      z_mult_dat->domain->nx *
                                          z_mult_dat->domain->ny * z] +
                     z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                      z_mult_dat->domain->nx *
                                          z_mult_dat->domain->ny * (z + 1)]) *
                        permzp->data[x + permzp->domain->nx * y +
                                     permzp->domain->nx * permzp->domain->ny *
                                         z] *
                        permzp->data[x + permzp->domain->nx * y +
                                     permzp->domain->nx * permzp->domain->ny *
                                         (z + 1)] /
                        (permzp->data[x + permzp->domain->nx * y +
                                      permzp->domain->nx * permzp->domain->ny *
                                          (z + 1)] *
                             z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                              z_mult_dat->domain->nx *
                                                  z_mult_dat->domain->ny * z] +
                         permzp->data[x + permzp->domain->nx * y +
                                      permzp->domain->nx * permzp->domain->ny *
                                          z] *
                             z_mult_dat
                                 ->data[x + z_mult_dat->domain->nx * y +
                                        z_mult_dat->domain->nx *
                                            z_mult_dat->domain->ny * (z + 1)])
                  : ((double)0))) *
            diff_upper *
            ((lower_cond - upper_cond >= ((double)0)
                  ? rpp->data[x + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * z] *
                        dp->data[x + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * z]
                  : rpp->data[x + rpp->domain->nx * y +
                              rpp->domain->nx * rpp->domain->ny * (z + 1)] *
                        dp->data[x + dp->domain->nx * y +
                                 dp->domain->nx * dp->domain->ny * (z + 1)]));
        vx->data[x + vx->domain->nx * y + vx->domain->nx * vx->domain->ny * z] =
            u_right;
        vy->data[x + vy->domain->nx * y + vy->domain->nx * vy->domain->ny * z] =
            u_front;
        vz->data[x + vz->domain->nx * y + vz->domain->nx * vz->domain->ny * z] =
            u_upper;
        fp->data[x + fp->domain->nx * y +
                 fp->domain->nx * fp->domain->ny * z] +=
            u_right * u_front * u_upper;
        fp->data[x + 1 + fp->domain->nx * y +
                 fp->domain->nx * fp->domain->ny * z] += u_right;
        fp->data[x + fp->domain->nx * (y + 1) +
                 fp->domain->nx * fp->domain->ny * z] -= u_front;
        fp->data[x + fp->domain->nx * y +
                 fp->domain->nx * fp->domain->ny * (z + 1)] -= u_upper;
      };
    }
  }
}
}
}
{
  for (int omplc_gen_iter_2 = 1; omplc_gen_iter_2 <= ny;
       omplc_gen_iter_2 = omplc_gen_iter_2 + 1) {
    for (int omplc_gen_iter_3 = 1; omplc_gen_iter_3 <= nz;
         omplc_gen_iter_3 = omplc_gen_iter_3 + 1) {
      {
        int z = omplc_gen_iter_3 - 1;
        int y = omplc_gen_iter_2 - 1;
        int x = nx - 1;
        {
          double x_dir_g =
              0.5 *
              (x_ssl_dat
                   ->data[x + x_ssl_dat->domain->nx * y +
                          x_ssl_dat->domain->nx * x_ssl_dat->domain->ny * 0] +
               x_ssl_dat
                   ->data[x + 1 + x_ssl_dat->domain->nx * y +
                          x_ssl_dat->domain->nx * x_ssl_dat->domain->ny * 0]);
          double x_dir_g_c =
              0.5 *
              (x_ssl_dat
                   ->data[x + x_ssl_dat->domain->nx * y +
                          x_ssl_dat->domain->nx * x_ssl_dat->domain->ny * 0] +
               x_ssl_dat
                   ->data[x + 1 + x_ssl_dat->domain->nx * y +
                          x_ssl_dat->domain->nx * x_ssl_dat->domain->ny * 0]);
          double y_dir_g =
              0.5 *
              (y_ssl_dat
                   ->data[x + y_ssl_dat->domain->nx * y +
                          y_ssl_dat->domain->nx * y_ssl_dat->domain->ny * 0] +
               y_ssl_dat
                   ->data[x + y_ssl_dat->domain->nx * (y + 1) +
                          y_ssl_dat->domain->nx * y_ssl_dat->domain->ny * 0]);
          double y_dir_g_c =
              0.5 *
              (y_ssl_dat
                   ->data[x + y_ssl_dat->domain->nx * y +
                          y_ssl_dat->domain->nx * y_ssl_dat->domain->ny * 0] +
               y_ssl_dat
                   ->data[x + y_ssl_dat->domain->nx * (y + 1) +
                          y_ssl_dat->domain->nx * y_ssl_dat->domain->ny * 0]);
          double diff_right = pp->data[x + pp->domain->nx * y +
                                       pp->domain->nx * pp->domain->ny * z] -
                              pp->data[x + 1 + pp->domain->nx * y +
                                       pp->domain->nx * pp->domain->ny * z];
          double diff_front = pp->data[x + pp->domain->nx * y +
                                       pp->domain->nx * pp->domain->ny * z] -
                              pp->data[x + pp->domain->nx * (y + 1) +
                                       pp->domain->nx * pp->domain->ny * z];
          double updir_right = diff_right * x_dir_g_c - x_dir_g;
          double updir_front = diff_front * y_dir_g_c - y_dir_g;
          double sep =
              0.5 *
              (z_mult_dat
                   ->data[x + z_mult_dat->domain->nx * y +
                          z_mult_dat->domain->nx * z_mult_dat->domain->ny * z] +
               z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                z_mult_dat->domain->nx *
                                    z_mult_dat->domain->ny * (z + 1)]);
          double lower_cond =
              pp->data[x + pp->domain->nx * y +
                       pp->domain->nx * pp->domain->ny * z] /
                  sep -
              z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                               z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                                   z] /
                  (z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                    z_mult_dat->domain->nx *
                                        z_mult_dat->domain->ny * z] +
                   z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                    z_mult_dat->domain->nx *
                                        z_mult_dat->domain->ny * (z + 1)]) *
                  dp->data[x + dp->domain->nx * y +
                           dp->domain->nx * dp->domain->ny * z];
          double upper_cond =
              pp->data[x + pp->domain->nx * y +
                       pp->domain->nx * pp->domain->ny * (z + 1)] /
                  sep +
              z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                               z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                                   (z + 1)] /
                  (z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                    z_mult_dat->domain->nx *
                                        z_mult_dat->domain->ny * z] +
                   z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                    z_mult_dat->domain->nx *
                                        z_mult_dat->domain->ny * (z + 1)]) *
                  dp->data[x + dp->domain->nx * y +
                           dp->domain->nx * dp->domain->ny * (z + 1)];
          double diff_upper = lower_cond - upper_cond;
          double u_right =
              z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                               z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                                   z] *
                  ((((bool)(permxp->data[x + permxp->domain->nx * y +
                                         permxp->domain->nx *
                                             permxp->domain->ny * z] +
                            permxp->data[x + 1 + permxp->domain->nx * y +
                                         permxp->domain->nx *
                                             permxp->domain->ny * z]))
                        ? 2.0 *
                              permxp->data[x + permxp->domain->nx * y +
                                           permxp->domain->nx *
                                               permxp->domain->ny * z] *
                              permxp->data[x + 1 + permxp->domain->nx * y +
                                           permxp->domain->nx *
                                               permxp->domain->ny * z] /
                              (permxp->data[x + permxp->domain->nx * y +
                                            permxp->domain->nx *
                                                permxp->domain->ny * z] +
                               permxp->data[x + 1 + permxp->domain->nx * y +
                                            permxp->domain->nx *
                                                permxp->domain->ny * z])
                        : ((double)0))) *
                  diff_right * x_dir_g_c *
                  ((updir_right - 0.0 >= ((double)0)
                        ? rpp->data[x + rpp->domain->nx * y +
                                    rpp->domain->nx * rpp->domain->ny * z] *
                              dp->data[x + dp->domain->nx * y +
                                       dp->domain->nx * dp->domain->ny * z]
                        : rpp->data[x + 1 + rpp->domain->nx * y +
                                    rpp->domain->nx * rpp->domain->ny * z] *
                              dp->data[x + 1 + dp->domain->nx * y +
                                       dp->domain->nx * dp->domain->ny * z])) +
              z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                               z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                                   z] *
                  ((((bool)(permxp->data[x + permxp->domain->nx * y +
                                         permxp->domain->nx *
                                             permxp->domain->ny * z] +
                            permxp->data[x + 1 + permxp->domain->nx * y +
                                         permxp->domain->nx *
                                             permxp->domain->ny * z]))
                        ? 2.0 *
                              permxp->data[x + permxp->domain->nx * y +
                                           permxp->domain->nx *
                                               permxp->domain->ny * z] *
                              permxp->data[x + 1 + permxp->domain->nx * y +
                                           permxp->domain->nx *
                                               permxp->domain->ny * z] /
                              (permxp->data[x + permxp->domain->nx * y +
                                            permxp->domain->nx *
                                                permxp->domain->ny * z] +
                               permxp->data[x + 1 + permxp->domain->nx * y +
                                            permxp->domain->nx *
                                                permxp->domain->ny * z])
                        : ((double)0))) *
                  -x_dir_g *
                  ((updir_right - 0.0 >= ((double)0)
                        ? rpp->data[x + rpp->domain->nx * y +
                                    rpp->domain->nx * rpp->domain->ny * z] *
                              dp->data[x + dp->domain->nx * y +
                                       dp->domain->nx * dp->domain->ny * z]
                        : rpp->data[x + 1 + rpp->domain->nx * y +
                                    rpp->domain->nx * rpp->domain->ny * z] *
                              dp->data[x + 1 + dp->domain->nx * y +
                                       dp->domain->nx * dp->domain->ny * z]));
          double u_front =
              z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                               z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                                   z] *
                  ((((bool)(permyp->data[x + permyp->domain->nx * y +
                                         permyp->domain->nx *
                                             permyp->domain->ny * z] +
                            permyp->data[x + 1 + permyp->domain->nx * y +
                                         permyp->domain->nx *
                                             permyp->domain->ny * z]))
                        ? 2.0 *
                              permyp->data[x + permyp->domain->nx * y +
                                           permyp->domain->nx *
                                               permyp->domain->ny * z] *
                              permyp->data[x + 1 + permyp->domain->nx * y +
                                           permyp->domain->nx *
                                               permyp->domain->ny * z] /
                              (permyp->data[x + permyp->domain->nx * y +
                                            permyp->domain->nx *
                                                permyp->domain->ny * z] +
                               permyp->data[x + 1 + permyp->domain->nx * y +
                                            permyp->domain->nx *
                                                permyp->domain->ny * z])
                        : ((double)0))) *
                  diff_front * x_dir_g_c *
                  ((updir_front - 0.0 >= ((double)0)
                        ? rpp->data[x + rpp->domain->nx * y +
                                    rpp->domain->nx * rpp->domain->ny * z] *
                              dp->data[x + dp->domain->nx * y +
                                       dp->domain->nx * dp->domain->ny * z]
                        : rpp->data[x + 1 + rpp->domain->nx * y +
                                    rpp->domain->nx * rpp->domain->ny * z] *
                              dp->data[x + 1 + dp->domain->nx * y +
                                       dp->domain->nx * dp->domain->ny * z])) +
              z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                               z_mult_dat->domain->nx * z_mult_dat->domain->ny *
                                   z] *
                  ((((bool)(permyp->data[x + permyp->domain->nx * y +
                                         permyp->domain->nx *
                                             permyp->domain->ny * z] +
                            permyp->data[x + 1 + permyp->domain->nx * y +
                                         permyp->domain->nx *
                                             permyp->domain->ny * z]))
                        ? 2.0 *
                              permyp->data[x + permyp->domain->nx * y +
                                           permyp->domain->nx *
                                               permyp->domain->ny * z] *
                              permyp->data[x + 1 + permyp->domain->nx * y +
                                           permyp->domain->nx *
                                               permyp->domain->ny * z] /
                              (permyp->data[x + permyp->domain->nx * y +
                                            permyp->domain->nx *
                                                permyp->domain->ny * z] +
                               permyp->data[x + 1 + permyp->domain->nx * y +
                                            permyp->domain->nx *
                                                permyp->domain->ny * z])
                        : ((double)0))) *
                  -x_dir_g *
                  ((updir_front - 0.0 >= ((double)0)
                        ? rpp->data[x + rpp->domain->nx * y +
                                    rpp->domain->nx * rpp->domain->ny * z] *
                              dp->data[x + dp->domain->nx * y +
                                       dp->domain->nx * dp->domain->ny * z]
                        : rpp->data[x + 1 + rpp->domain->nx * y +
                                    rpp->domain->nx * rpp->domain->ny * z] *
                              dp->data[x + 1 + dp->domain->nx * y +
                                       dp->domain->nx * dp->domain->ny * z]));
          double u_upper =
              ((((bool)(z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                         z_mult_dat->domain->nx *
                                             z_mult_dat->domain->ny * z] *
                            permzp->data[x + permzp->domain->nx * y +
                                         permzp->domain->nx *
                                             permzp->domain->ny * (z + 1)] +
                        permzp->data[x + permzp->domain->nx * y +
                                     permzp->domain->nx * permzp->domain->ny *
                                         z] *
                            z_mult_dat
                                ->data[x + z_mult_dat->domain->nx * y +
                                       z_mult_dat->domain->nx *
                                           z_mult_dat->domain->ny * (z + 1)]))
                    ? (z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                        z_mult_dat->domain->nx *
                                            z_mult_dat->domain->ny * z] +
                       z_mult_dat->data[x + z_mult_dat->domain->nx * y +
                                        z_mult_dat->domain->nx *
                                            z_mult_dat->domain->ny * (z + 1)]) *
                          permzp->data[x + permzp->domain->nx * y +
                                       permzp->domain->nx * permzp->domain->ny *
                                           z] *
                          permzp->data[x + permzp->domain->nx * y +
                                       permzp->domain->nx * permzp->domain->ny *
                                           (z + 1)] /
                          (permzp->data[x + permzp->domain->nx * y +
                                        permzp->domain->nx *
                                            permzp->domain->ny * (z + 1)] *
                               z_mult_dat
                                   ->data[x + z_mult_dat->domain->nx * y +
                                          z_mult_dat->domain->nx *
                                              z_mult_dat->domain->ny * z] +
                           permzp->data[x + permzp->domain->nx * y +
                                        permzp->domain->nx *
                                            permzp->domain->ny * z] *
                               z_mult_dat
                                   ->data[x + z_mult_dat->domain->nx * y +
                                          z_mult_dat->domain->nx *
                                              z_mult_dat->domain->ny * (z + 1)])
                    : ((double)0))) *
              diff_upper *
              ((lower_cond - upper_cond >= ((double)0)
                    ? rpp->data[x + rpp->domain->nx * y +
                                rpp->domain->nx * rpp->domain->ny * z] *
                          dp->data[x + dp->domain->nx * y +
                                   dp->domain->nx * dp->domain->ny * z]
                    : rpp->data[x + rpp->domain->nx * y +
                                rpp->domain->nx * rpp->domain->ny * (z + 1)] *
                          dp->data[x + dp->domain->nx * y +
                                   dp->domain->nx * dp->domain->ny * (z + 1)]));
          vx->data[x + vx->domain->nx * y +
                   vx->domain->nx * vx->domain->ny * z] = u_right;
          vy->data[x + vy->domain->nx * y +
                   vy->domain->nx * vy->domain->ny * z] = u_front;
          vz->data[x + vz->domain->nx * y +
                   vz->domain->nx * vz->domain->ny * z] = u_upper;
          fp->data[x + fp->domain->nx * y +
                   fp->domain->nx * fp->domain->ny * z] +=
              u_right * u_front * u_upper;
          fp->data[x + 1 + fp->domain->nx * y +
                   fp->domain->nx * fp->domain->ny * z] += u_right;
          fp->data[x + fp->domain->nx * (y + 1) +
                   fp->domain->nx * fp->domain->ny * z] -= u_front;
          fp->data[x + fp->domain->nx * y +
                   fp->domain->nx * fp->domain->ny * (z + 1)] -= u_upper;
        };
      }
    }
  }
}
}
}
}
// NlFunctionEval:261 analogue
