#ifndef VARIANT_MACROS_HPP
#define VARIANT_MACROS_HPP

#include <global_macros.hpp>

// Basic Domain
#define Variant_Domain_nx(domain) Basic_Domain_nx(domain)
#define Variant_Domain_ny(domain) Basic_Domain_ny(domain)
#define Variant_Domain_nz(domain) Basic_Domain_nz(domain)
#define Variant_Domain_idx(domain,x,y,z) Basic_Domain_idx(domain,x,y,z)

// Basic Grid
#define Variant_Grid_data(grid) Basic_Grid_data(grid)
#define Variant_Grid_domain(grid) Basic_Grid_domain(grid)
// The below should be unused
// #define Basic_Grid_access_index(grid,idx) (Basic_Grid_data(grid)[idx])
// #define Basic_Grid_access(grid,x,y,z) (Basic_Grid_access_index(grid,Basic_Domain_idx(Basic_Grid_domain(grid),x,y,z)))
#define Variant_Grid_idx(grid, x,y,z) Basic_Grid_idx(grid, x,y,z)
#define Variant_Grid_access_index(grid, idx) Basic_Grid_access_index(grid, idx)
#define Variant_Grid_access(grid, x,y,z) Basic_Grid_access(grid, x,y,z)

#define Variant_Domain_equal(domain_a, domain_b) Basic_Domain_equal(domain_a, domain_b)

#define Variant_Domain_loop_whole(domain, iter_x, iter_y, iter_z, body) Basic_Domain_loop_whole(domain, iter_x, iter_y, iter_z, body)

#define Variant_Domain_fast_loop_whole(domain, iter_x, iter_y, iter_z, body) \
  _Pragma( STRINGIZE(omp parallel for private(iter_x,iter_y,iter_z)) ) \
  for( iter_x = 0; iter_x < Variant_Domain_nx(domain); iter_x += 1 ){ \
    for( iter_y = 0; iter_y < Variant_Domain_ny(domain); iter_y += 1 ){ \
      for( iter_z = 0; iter_z < Variant_Domain_nz(domain); iter_z += 1 ){ \
        body; \
      } \
    } \
  }

#define Variant_Domain_loop_interior(domain, iter_x, iter_y, iter_z, body) Basic_Domain_loop_interior(domain, iter_x, iter_y, iter_z, body)

#define Variant_Domain_fast_loop_interior(domain, iter_x, iter_y, iter_z, body) \
  _Pragma( STRINGIZE(omp parallel for private(iter_x,iter_y,iter_z)) ) \
  for( iter_x = 1; iter_x < Variant_Domain_nx(domain)-1; iter_x += 1 ){ \
    for( iter_y = 1; iter_y < Variant_Domain_ny(domain)-1; iter_y += 1 ){ \
      for( iter_z = 1; iter_z < Variant_Domain_nz(domain)-1; iter_z += 1 ){ \
        body; \
      } \
    } \
  }

#endif
