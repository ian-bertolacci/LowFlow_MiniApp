#include <variant_types.hpp>
#include <macros.hpp>
#include <configure.hpp>
#include <cstring>
#include <cassert>

#include <global_types.hpp>
#include <global_macros.hpp>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>

Variant_Domain* Variant_Domain_alloc( int nx, int ny, int nz ){
  Variant_Domain* result = (Variant_Domain*) malloc( sizeof(Variant_Domain) );
  assert( result != nullptr && "Error during Variant_Domain_Alloc: Unable to allocate Variant_Domain struct" );

  Variant_Domain_nx(result) = nx;
  Variant_Domain_ny(result) = ny;
  Variant_Domain_nz(result) = nz;

  return result;
}

void Variant_Domain_dealloc( Variant_Domain* domain ){
  assert( domain != nullptr && "Error during Variant_Domain_dealloc: Variant_Domain pointer is nullptr" );
  free(domain);
}

Variant_Grid* Variant_Grid_alloc( Variant_Domain* domain ){
  Variant_Grid* result = (Variant_Grid*) malloc( sizeof(Variant_Grid) );
  assert( result != nullptr && "Error during Variant_Grid_Alloc: Unable to allocate Variant_Grid struct" );

  Variant_Grid_domain(result) = domain;

  result->data = new view_3D_double_t( "data", Variant_Domain_nx(domain), Variant_Domain_ny(domain), Variant_Domain_nz(domain) );

  assert( Variant_Grid_data(result) != nullptr && "Error during Variant_Grid_alloc: Unable to allocate Variant_Grid data" );

  return result;
}

void Variant_Grid_dealloc( Variant_Grid* grid ){
  assert( grid != nullptr && "Error during Variant_Grid_dealloc: Variant_Grid pointer is nullptr" );
  assert( Variant_Grid_data(grid) != nullptr && "Error during Variant_Grid_dealloc: Variant_Grid has Null data array");

  delete Variant_Grid_data(grid);

  free( grid );
}

void Variant_Grid_populate_int_increment( Variant_Domain* domain, Variant_Grid* grid ){
  assert( domain != nullptr && "Error during Variant_Grid_populate: Variant_Domain pointer is nullptr" );
  assert( grid != nullptr && "Error during Variant_Grid_populate: Variant_Grid pointer is nullptr" );
  double value = 1.0;
  int x, y, z;
  Variant_Domain_loop_whole(domain, x,y,z,
    {
      Variant_Grid_access( grid, x, y, z ) = (double) value;
      value += 1;
    }
  );
}

void Variant_Grid_populate_zero( Variant_Domain* domain, Variant_Grid* grid ){
  assert( domain != nullptr && "Error during Variant_Grid_populate: Variant_Domain pointer is nullptr" );
  assert( grid != nullptr && "Error during Variant_Grid_populate: Variant_Grid pointer is nullptr" );
  int x, y, z;
  Variant_Domain_loop_whole(domain, x,y,z,
    {
      Variant_Grid_access( grid, x, y, z) = 0.0;
    }
  );
}

void Variant_Grid_populate_seeded( Variant_Domain* domain, Variant_Grid* grid, unsigned int seed ){
  assert( domain != nullptr && "Error during Variant_Grid_populate_seeded: Variant_Domain pointer is nullptr" );
  assert( grid != nullptr && "Error during Variant_Grid_populate_seeded: Variant_Grid pointer is nullptr" );
  srand(seed);
  int x, y, z;
  Variant_Domain_loop_whole(domain, x,y,z,
    {
      // get signed value from unsigned source
      int rand_value = (int) rand();
      double sign = (double)((rand_value>=0)?1:-1);
      int abs_val = (int) abs(rand_value);
      int mod_factor = 100000;
      double mod = (double)(abs_val % mod_factor/2);

      double value = sign*(mod/(mod_factor));
      Variant_Grid_access(grid, x,y,z) = value;
    }
  );
}
