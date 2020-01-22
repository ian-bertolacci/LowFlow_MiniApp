#include <variant_types.hpp>
#include <configure.hpp>
#include <assert.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

Variant_Domain* Variant_Domain_alloc( int nx, int ny, int nz ){
  Variant_Domain *result = nullptr;
  cudaError_t x = cudaMallocManaged(&result, sizeof(Variant_Domain));
  assert(x == cudaSuccess);

  //Fails here, when assigning nx/ny/nz. The above assert passes, so the memory should
  //be allocated correctly. Am I missing something?
  result->nx = nx;
  result->ny = ny;
  result->nz = nz;

  return result;
}

void Variant_Domain_dealloc( Variant_Domain* domain ){
  assert(cudaFree(domain) == cudaSuccess);
}

Variant_Grid* Variant_Grid_alloc( Variant_Domain* domain ){
  Variant_Grid *result = nullptr;
  cudaError_t x = cudaMallocManaged(&result, sizeof(Variant_Grid));
  assert(x == cudaSuccess);
  
  result->domain = domain;
  x = cudaMallocManaged(&result->data, (domain->nx * domain->ny * domain->nz) * sizeof(double));
  assert(x == cudaSuccess);

  return result;
}


void Variant_Grid_dealloc( Variant_Grid* grid ){
  assert(cudaFree(grid->data) == cudaSuccess);
  assert(cudaFree(grid) == cudaSuccess);
}


void Variant_Grid_populate_int_increment( Variant_Domain* domain, Variant_Grid* grid ){
  return Basic_Grid_populate_int_increment( (Basic_Domain*) domain, (Basic_Grid*) grid );
}

void Variant_Grid_populate_zero( Variant_Domain* domain, Variant_Grid* grid ){
  return Basic_Grid_populate_zero( (Basic_Domain*) domain, (Basic_Grid*) grid );
}

void Variant_Grid_populate_seeded( Variant_Domain* domain, Variant_Grid* grid, unsigned int seed ){
  return Basic_Grid_populate_seeded( (Basic_Domain*) domain, (Basic_Grid*) grid, seed );
}
