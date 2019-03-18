#ifndef GLOBAL_MACROS_HPP
#define GLOBAL_MACROS_HPP

// Basic Domain
#define Basic_Domain_nx(domain) (domain->nx)
#define Basic_Domain_ny(domain) (domain->ny)
#define Basic_Domain_nz(domain) (domain->nz)
#define Basic_Domain_idx(domain, x,y,z) (x+(Basic_Domain_nx(domain)*y)+(Basic_Domain_nx(domain)*Basic_Domain_ny(domain)*z))


// Basic Grid
#define Basic_Grid_data(grid) (grid->data)
#define Basic_Grid_domain(grid) (grid->domain)
// The below should be unused
// #define Basic_Grid_access_index(grid,idx) (Basic_Grid_data(grid)[idx])
// #define Basic_Grid_access(grid,x,y,z) (Basic_Grid_access_index(grid,Basic_Domain_idx(Basic_Grid_domain(grid),x,y,z)))
#define Basic_Grid_idx(grid, x,y,z) (Basic_Domain_idx(Basic_Grid_domain(grid), x,y,z))
#define Basic_Grid_access_index(grid, idx) error "Do not use Basic_Grid_access_index"
#define Basic_Grid_access(grid, x,y,z) error "Do not use Basic_Grid_access"

#define Basic_Domain_equal(domain_a, domain_b) \
  (Basic_Domain_nx(domain_a) == Basic_Domain_nx(domain_b) \
&& Basic_Domain_ny(domain_a) == Basic_Domain_ny(domain_b) \
&& Basic_Domain_nz(domain_a) == Basic_Domain_nz(domain_b))

#define Basic_Domain_loop_whole(domain, iter_x, iter_y, iter_z, body) \
  for( iter_x = 0; iter_x < Basic_Domain_nx(domain); iter_x += 1 ){ \
    for( iter_y = 0; iter_y < Basic_Domain_ny(domain); iter_y += 1 ){ \
      for( iter_z = 0; iter_z < Basic_Domain_nz(domain); iter_z += 1 ){ \
        body; \
      } \
    } \
  }

#define Basic_Domain_fast_loop_whole(domain, iter_x, iter_y, iter_z, body) \
  Basic_Domain_loop_whole(domain, iter_x, iter_y, iter_z, body)

#define Basic_Domain_loop_interior(domain, iter_x, iter_y, iter_z, body) \
  for( iter_x = 1; iter_x < Basic_Domain_nx(domain)-1; iter_x += 1 ){ \
    for( iter_y = 1; iter_y < Basic_Domain_ny(domain)-1; iter_y += 1 ){ \
      for( iter_z = 1; iter_z < Basic_Domain_nz(domain)-1; iter_z += 1 ){ \
        body; \
      } \
    } \
  }

#define Basic_Domain_fast_loop_interior(domain, iter_x, iter_y, iter_z, body) \
  Basic_Domain_loop_interior(domain, iter_x, iter_y, iter_z, body)

#endif
