#include <variant_types.hpp>
#include <configure.hpp>

Variant_Domain* Variant_Domain_alloc( int nx, int ny, int nz ){
  return (Variant_Domain*) Basic_Domain_alloc( nx, ny, nz );
}


void Variant_Domain_dealloc( Variant_Domain* domain ){
  return Basic_Domain_dealloc( (Basic_Domain*) domain );
}



Variant_Grid* Variant_Grid_alloc( Variant_Domain* domain ){
  return (Variant_Grid*) Basic_Grid_alloc( (Basic_Domain*) domain );
}


void Variant_Grid_dealloc( Variant_Grid* grid ){
  return Basic_Grid_dealloc( (Basic_Grid*) grid );
}


void Variant_Grid_populate( Variant_Domain* domain, Variant_Grid* grid ){
  return Basic_Grid_populate( (Basic_Domain*) domain, (Basic_Grid*) grid );
}

void Variant_Grid_populate_seeded( Variant_Domain* domain, Variant_Grid* grid, unsigned int seed ){
  return Basic_Grid_populate_seeded( (Basic_Domain*) domain, (Basic_Grid*) grid, seed );
}
