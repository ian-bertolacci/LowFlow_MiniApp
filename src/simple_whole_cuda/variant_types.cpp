#include <variant_types.hpp>
#include <configure.hpp>

Variant_Domain* Variant_Domain_alloc( int nx, int ny, int nz ){
  Variant_Domain *result = nullptr;
  if (cudaMallocManaged(&result, sizeof(Variant_Domain)) != cudaSuccess) {
    return nullptr;
  }
  
  result->nx = nx;
  result->ny = ny;
  result->nz = nz;

  return result;
}

void Variant_Domain_dealloc( Variant_Domain* domain ){
  cudaFree(domain);
}

Variant_Grid* Variant_Grid_alloc( Variant_Domain* domain ){
  Variant_Grid *result = nullptr;
  if (cudaMallocManaged(&result, sizeof(Variant_Grid)) != cudaSuccess) {
    return nullptr;
  };
  
  result->domain = domain;
  if (cudaMallocManaged(&result->data, (domain->nx * domain->ny * domain->nz) * sizeof(double)) != cudaSuccess) {
    return nullptr;
  }

  return result;
}


void Variant_Grid_dealloc( Variant_Grid* grid ){
  cudaFree(grid->data);
  cudaFree(grid);
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
