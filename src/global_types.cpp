#include <global_types.hpp>
#include <global_macros.hpp>
#include <stdlib.h>
#include <assert.h>

Basic_Domain* Basic_Domain_alloc( int nx, int ny, int nz ){
  Basic_Domain* result = (Basic_Domain*) malloc( sizeof(Basic_Domain) );
  assert( result != nullptr && "Error during Basic_Domain_Alloc: Unable to allocate Basic_Domain struct" );

  Basic_Domain_nx(result) = nx;
  Basic_Domain_ny(result) = ny;
  Basic_Domain_nz(result) = nz;

  return result;
}

void Basic_Domain_dealloc( Basic_Domain* domain ){
  assert( domain != nullptr && "Error during Basic_Domain_dealloc: Basic_Domain pointer is nullptr" );
  free(domain);
}

Basic_Grid* Basic_Grid_alloc( Basic_Domain* domain ){
  Basic_Grid* result = (Basic_Grid*) malloc( sizeof(Basic_Grid) );
  assert( result != nullptr && "Error during Basic_Grid_Alloc: Unable to allocate Basic_Grid struct" );

  Basic_Grid_domain(result) = domain;

  Basic_Grid_data(result) = (double*) malloc( (Basic_Domain_nx(domain)*Basic_Domain_ny(domain)*Basic_Domain_nz(domain)) * sizeof( double ) );
  assert( Basic_Grid_data(result) != nullptr && "Error during Basic_Grid_alloc: Unable to allocate Basic_Grid data" );

  return result;
}

void Basic_Grid_dealloc( Basic_Grid* grid ){
  assert( grid != nullptr && "Error during Basic_Grid_dealloc: Basic_Grid pointer is nullptr" );
  assert( Basic_Grid_data(grid) != nullptr && "Error during Basic_Grid_dealloc: Basic_Grid has Null data array");

  free( Basic_Grid_data(grid) );
  free( grid );
}

void Basic_Grid_populate( Basic_Domain* domain, Basic_Grid* grid ){
  assert( domain != nullptr && "Error during Basic_Grid_populate: Basic_Domain pointer is nullptr" );
  assert( grid != nullptr && "Error during Basic_Grid_populate: Basic_Grid pointer is nullptr" );
  double val = 1.0;
  int x, y, z;
  double* grid_data = Basic_Grid_data(grid);
  Basic_Domain_loop_whole(domain, x,y,z,
    {
      int idx = Basic_Domain_idx(domain, x,y,z);
      grid_data[idx] = val;
      val += 1.0;
    }
  );
}
