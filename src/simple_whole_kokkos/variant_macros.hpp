#ifndef VARIANT_MACROS_HPP
#define VARIANT_MACROS_HPP

#include <global_macros.hpp>
#include <iostream>
#include <Kokkos_Core.hpp>
// #include <Kokkos_OpenMP.hpp>
// #include <Kokkos_Serial.hpp>

// Variant Domain
#define Variant_Domain_nx(domain) ((domain)->nx)
#define Variant_Domain_ny(domain) ((domain)->ny)
#define Variant_Domain_nz(domain) ((domain)->nz)
#define Variant_Domain_idx(domain, x,y,z) ((x)+(Variant_Domain_nx(domain)*(y))+(Variant_Domain_nx(domain)*Variant_Domain_ny(domain)*(z)))

// Variant Grid
#define Variant_Grid_data(grid) ((grid)->data)
#define Variant_Grid_domain(grid) ((grid)->domain)
#define Variant_Grid_idx(grid, x,y,z) (Variant_Domain_idx(Variant_Grid_domain(grid), x,y,z))
// #define Variant_Grid_access_index(grid,idx) (Variant_Grid_data(grid)[idx])
#define Variant_Grid_access_index(grid, idx) error "Do not use Variant_Grid_access_index"
#define Variant_Grid_access(grid,x,y,z) ((*Variant_Grid_data(grid))(x,y,z))

#define Variant_Domain_equal(domain_a, domain_b) \
  (Variant_Domain_nx(domain_a) == Variant_Domain_nx(domain_b) \
&& Variant_Domain_ny(domain_a) == Variant_Domain_ny(domain_b) \
&& Variant_Domain_nz(domain_a) == Variant_Domain_nz(domain_b))

#define Variant_Domain_loop_whole(domain, iter_x, iter_y, iter_z, body) \
  for( iter_x = 0; iter_x < Variant_Domain_nx(domain); iter_x += 1 ){ \
    for( iter_y = 0; iter_y < Variant_Domain_ny(domain); iter_y += 1 ){ \
      for( iter_z = 0; iter_z < Variant_Domain_nz(domain); iter_z += 1 ){ \
        body; \
      } \
    } \
  }

#define Variant_Domain_loop_interior(domain, iter_x, iter_y, iter_z, body) \
  for( iter_x = 1; iter_x < Variant_Domain_nx(domain)-1; iter_x += 1 ){ \
    for( iter_y = 1; iter_y < Variant_Domain_ny(domain)-1; iter_y += 1 ){ \
      for( iter_z = 1; iter_z < Variant_Domain_nz(domain)-1; iter_z += 1 ){ \
        body; \
      } \
    } \
  }

#define Variant_Domain_fast_loop_in_bounds(domain, \
  iter_x, lower_bound_x, upper_bound_x, \
  iter_y, lower_bound_y, upper_bound_y, \
  iter_z, lower_bound_z, upper_bound_z, \
  body) \
  { \
    /* std::cout << "Hello from Variant_Domain_fast_loop_in_bounds" << std::endl; */ \
    LowFlow_3D_Execution_Fast_Policy __local_Variant_Domain_fast_loop_in_bounds_policy( \
      {{lower_bound_x, lower_bound_y, lower_bound_z}}, \
      {{upper_bound_x, upper_bound_y, upper_bound_z}} \
    ); \
    Kokkos::parallel_for( "__local_Variant_Domain_fast_loop_in_bounds_parallel_for", __local_Variant_Domain_fast_loop_in_bounds_policy, \
      KOKKOS_LAMBDA (const int iter_x, const int iter_y, const int iter_z ){ \
        /* std::cout << "Iteration " << iter_x << ", " << iter_y << ", " << iter_z << std::endl; */ \
        body; \
      } \
    ); \
  }

#define Variant_Domain_fast_loop_whole(domain, iter_x, iter_y, iter_z, body) \
  Variant_Domain_fast_loop_in_bounds( domain, \
    iter_x, 0, Variant_Domain_nx(domain), \
    iter_y, 0, Variant_Domain_ny(domain), \
    iter_z, 0, Variant_Domain_nz(domain), \
    body \
  )



#define Variant_Domain_fast_loop_interior(domain, iter_x, iter_y, iter_z, body) \
  Variant_Domain_fast_loop_in_bounds( domain, \
    iter_x, 1, Variant_Domain_nx(domain)-1, \
    iter_y, 1, Variant_Domain_ny(domain)-1, \
    iter_z, 1, Variant_Domain_nz(domain)-1, \
    body \
  )
#endif
