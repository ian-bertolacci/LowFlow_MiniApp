#ifndef VARIANT_MACROS_HPP
#define VARIANT_MACROS_HPP

#include <global_macros.hpp>
#include <iostream>
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

// #include <Kokkos_OpenMP.hpp>
// #include <Kokkos_Serial.hpp>

// Variant Domain
#define Variant_Domain_nx(domain) ((domain)->nx)
#define Variant_Domain_ny(domain) ((domain)->ny)
#define Variant_Domain_nz(domain) ((domain)->nz)
#define Variant_Domain_idx(domain, x,y,z) ((x)+(Variant_Domain_nx(domain)*(y))+(Variant_Domain_nx(domain)*Variant_Domain_ny(domain)*(z)))
#define Variant_Grid_data(grid) ((grid)->data)

// Variant Grid
// #ifdef USE_CUDA
// #define Variant_Grid_unwrap_data(grid) ( is_on_cpu ? Variant_Grid_data(grid)->template view<view_3D_double_t::host_mirror_space>() : Variant_Grid_data(grid)->template view<view_3D_double_t::memory_space>() )
// #else
//   #define Variant_Grid_unwrap_data(grid) (*((grid)->data))
// #endif

#define Variant_Grid_domain(grid) ((grid)->domain)
#define Variant_Grid_idx(grid, x,y,z) (Variant_Domain_idx(Variant_Grid_domain(grid), x,y,z))(grid)->data
// #define Variant_Grid_access_index(grid,idx) (Variant_Grid_data(grid)[idx])
#define Variant_Grid_access_index(grid, idx) error "Do not use Variant_Grid_access_index"
#define Variant_Grid_access(grid,x,y,z) ( is_on_cpu ? Variant_Grid_data(grid)->template view<view_3D_double_t::host_mirror_space>()(x,y,z) : Variant_Grid_data(grid)->template view<view_3D_double_t::memory_space>()(x,y,z) )

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

#define unwrap_to_device(arg_name) arg_name->data->template view<view_3D_double_t::memory_space>();
#define unwrap_to_host(arg_name) arg_name->data->template view<view_3D_double_t::host_mirror_space>();
#define sync( arg_name, from, to ) \
  arg_name->data->modify<view_3D_double_t::from>(); \
  arg_name->data->sync<view_3D_double_t::to>();
#define sync_from_device_to_host( arg_name ) sync( arg_name, memory_space, host_mirror_space)

#define writes(...) DELAY(__VA_ARGS__)
#define reads(...) DELAY(__VA_ARGS__)

#define Variant_Domain_fast_loop_in_bounds( domain, \
  iter_x, lower_bound_x, upper_bound_x, \
  iter_y, lower_bound_y, upper_bound_y, \
  iter_z, lower_bound_z, upper_bound_z, \
  _writes, \
  _reads, \
  body \
) \
{ \
  FOR_EACH( unwrap_to_device, _reads ) \
  { \
    FOR_EACH( unwrap_to_device, _writes ) \
    LowFlow_3D_Execution_Fast_Policy __local_Variant_Domain_fast_loop_in_bounds_policy( \
      {{lower_bound_x, lower_bound_y, lower_bound_z}}, \
      {{upper_bound_x, upper_bound_y, upper_bound_z}} \
    ); \
    bool is_on_cpu = false; \
    Kokkos::parallel_for( "__local_Variant_Domain_fast_loop_in_bounds_parallel_for", __local_Variant_Domain_fast_loop_in_bounds_policy, \
      KOKKOS_LAMBDA (const int iter_x, const int iter_y, const int iter_z ){ \
        body; \
      } \
    ); \
    FOR_EACH( sync_from_device_to_host, _writes ) \
  } \
}


#define Variant_Domain_fast_loop_whole(domain, iter_x, iter_y, iter_z, _writes, _reads, body) \
  Variant_Domain_fast_loop_in_bounds( domain, \
    iter_x, 0, Variant_Domain_nx(domain), \
    iter_y, 0, Variant_Domain_ny(domain), \
    iter_z, 0, Variant_Domain_nz(domain), \
    DELAY(_writes), \
    DELAY(_reads), \
    body \
  )



#define Variant_Domain_fast_loop_interior(domain, iter_x, iter_y, iter_z, _writes, _reads, body) \
  Variant_Domain_fast_loop_in_bounds( domain, \
    iter_x, 1, Variant_Domain_nx(domain)-1, \
    iter_y, 1, Variant_Domain_ny(domain)-1, \
    iter_z, 1, Variant_Domain_nz(domain)-1, \
    DELAY(_writes), \
    DELAY(_reads), \
    body \
  )
#endif
