#ifndef VARIANT_TYPES_HPP
#define VARIANT_TYPES_HPP

#include <global_types.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <configure.hpp>

static int is_on_cpu = true;

typedef struct struct_Variant_Domain {
  // int x, y, z; // TODO Include positions? Would be used for multi-grid situations (like distributed solves)
  int nx, ny, nz;
} Variant_Domain;

// Note: The Execution Space (kokkos_execution_space) symbol is defined in variant_configure.hpp
typedef typename Kokkos::Experimental::MDRangePolicy< kokkos_execution_space, Kokkos::Experimental::Rank<3, Kokkos::Experimental::Iterate::Left, Kokkos::Experimental::Iterate::Left> > LowFlow_3D_Execution_Fast_Policy;

#ifdef USE_CUDA
  typedef Kokkos::DualView<double***> view_3D_double_t;
#else
  typedef Kokkos::View<double***> view_3D_double_t;
#endif

typedef struct struct_Variant_Grid {
  Variant_Domain* domain;
  view_3D_double_t* data;
} Variant_Grid;

Variant_Domain* Variant_Domain_alloc( int nx, int ny, int nz );
void Variant_Domain_dealloc( Variant_Domain* domain );

Variant_Grid* Variant_Grid_alloc( Variant_Domain* domain );
void Variant_Grid_dealloc( Variant_Grid* grid );
void Variant_Grid_populate_int_increment( Variant_Domain* domain, Variant_Grid* grid );
void Variant_Grid_populate_zero( Variant_Domain* domain, Variant_Grid* grid );
void Variant_Grid_populate_seeded( Variant_Domain* domain, Variant_Grid* grid, unsigned int seed );


#endif
