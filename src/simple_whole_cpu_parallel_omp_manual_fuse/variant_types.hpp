#ifndef VARIANT_TYPES_HPP
#define VARIANT_TYPES_HPP

#include <global_types.hpp>

typedef Basic_Domain Variant_Domain;
typedef Basic_Grid Variant_Grid;

Variant_Domain* Variant_Domain_alloc( int nx, int ny, int nz );
void Variant_Domain_dealloc( Variant_Domain* domain );

Variant_Grid* Variant_Grid_alloc( Variant_Domain* domain );
void Variant_Grid_dealloc( Variant_Grid* grid );
void Variant_Grid_populate_int_increment( Variant_Domain* domain, Variant_Grid* grid );
void Variant_Grid_populate_zero( Variant_Domain* domain, Variant_Grid* grid );
void Variant_Grid_populate_seeded( Variant_Domain* domain, Variant_Grid* grid, unsigned int seed );

#endif
