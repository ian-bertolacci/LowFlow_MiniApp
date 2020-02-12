#ifndef VARIANT_MACROS_HPP
#define VARIANT_MACROS_HPP

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

#define Variant_Domain_loop_interior(domain, iter_x, iter_y, iter_z, body) Basic_Domain_loop_interior(domain, iter_x, iter_y, iter_z, body)

#define Variant_Domain_fast_loop_whole(domain, iter_x, iter_y, iter_z, body) Basic_Domain_loop_whole(domain, iter_x, iter_y, iter_z, body)

#define Variant_Domain_fast_loop_interior(domain, body) { \
	int blockx = 16; \
	int blocky = 16; \
	int blockz = 4; \
	int gridx = (int) ceil((float) Basic_Domain_nx(domain) / blockx); \
	int gridy = (int) ceil((float) Basic_Domain_ny(domain) / blocky); \
	int gridz = (int) ceil((float) Basic_Domain_nz(domain) / blockz); \
	kernelWrapper<<<dim3(gridx, gridy, gridz),  dim3(blockx, blocky, blockz)>>>([=] __device__ { \
		int x = blockIdx.x * blockDim.x + threadIdx.x + 1; \
		int y = blockIdx.y * blockDim.y + threadIdx.y + 1; \
		int z = blockIdx.z * blockDim.z + threadIdx.z + 1; \
		if (x == 0 || x > Basic_Domain_nx(domain) - 2 || y == 0 || y > Basic_Domain_ny(domain) - 2 || z == 0 || z > Basic_Domain_nz(domain) - 2) { return; } \
		body; \
	}); \
  	cudaDeviceSynchronize(); \
} 	

#endif
