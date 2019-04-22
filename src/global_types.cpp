#include <global_types.hpp>
#include <global_macros.hpp>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>

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

void Basic_Grid_populate_int_increment( Basic_Domain* domain, Basic_Grid* grid ){
  assert( domain != nullptr && "Error during Basic_Grid_populate: Basic_Domain pointer is nullptr" );
  assert( grid != nullptr && "Error during Basic_Grid_populate: Basic_Grid pointer is nullptr" );
  double value = 1.0;
  int x, y, z;
  double* grid_data = Basic_Grid_data(grid);
  Basic_Domain_loop_whole(domain, x,y,z,
    {
      int idx = Basic_Domain_idx(domain, x,y,z);
      grid_data[idx] = (double) value;
      value += 1;
    }
  );
}

void Basic_Grid_populate_zero( Basic_Domain* domain, Basic_Grid* grid ){
  assert( domain != nullptr && "Error during Basic_Grid_populate: Basic_Domain pointer is nullptr" );
  assert( grid != nullptr && "Error during Basic_Grid_populate: Basic_Grid pointer is nullptr" );
  int x, y, z;
  double* grid_data = Basic_Grid_data(grid);
  Basic_Domain_loop_whole(domain, x,y,z,
    {
      int idx = Basic_Domain_idx(domain, x,y,z);
      grid_data[idx] = 0.0;
    }
  );
}

void Basic_Grid_populate_seeded( Basic_Domain* domain, Basic_Grid* grid, unsigned int seed ){
  assert( domain != nullptr && "Error during Basic_Grid_populate_seeded: Basic_Domain pointer is nullptr" );
  assert( grid != nullptr && "Error during Basic_Grid_populate_seeded: Basic_Grid pointer is nullptr" );
  srand(seed);
  int x, y, z;
  double* grid_data = Basic_Grid_data(grid);
  // printf("hi\n");
  Basic_Domain_loop_whole(domain, x,y,z,
    {
      int idx = Basic_Domain_idx(domain, x,y,z);
      // get signed value from unsigned source
      int rand_value = (int) rand();
      double sign = (double)((rand_value>=0)?1:-1);
      int abs_val = (int) abs(rand_value);
      int mod_factor = 100000;
      double mod = (double)(abs_val % mod_factor/2);

      double value = sign*(mod/(mod_factor));
      Basic_Grid_access(grid, x,y,z) = value;
      // printf( "(%d, %d, %d) %f\n", x,y,z, Basic_Grid_access(grid, x,y,z) );
    }
  );
}
