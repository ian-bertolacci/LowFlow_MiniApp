#ifndef GLOBAL_TYPES_HPP
#define GLOBAL_TYPES_HPP

typedef struct struct_Basic_Domain {
  // int x, y, z; // TODO Include positions? Would be used for multi-grid situations (like distributed solves)
  int nx, ny, nz;
} Basic_Domain;

typedef struct struct_Basic_Grid {
  Basic_Domain* domain;
  double* data;
} Basic_Grid;

Basic_Domain* Basic_Domain_alloc( int nx, int ny, int nz );
void Basic_Domain_dealloc( Basic_Domain* domain );

Basic_Grid* Basic_Grid_alloc( Basic_Domain* domain );
void Basic_Grid_dealloc( Basic_Grid* grid );
void Basic_Grid_populate( Basic_Domain* domain, Basic_Grid* grid );
void Basic_Grid_populate_seeded( Basic_Domain* domain, Basic_Grid* grid, unsigned int seed );

#endif
