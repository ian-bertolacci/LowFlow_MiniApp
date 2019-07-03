#include <setup_teardown.hpp>
#include <Kokkos_Core.hpp>

void programSetup( ProgramOptions program_options ){
  Kokkos::initialize( program_options.all_argument_parameters.argc, program_options.all_argument_parameters.argv );
}

void programTeardown( ProgramOptions program_options ){
  Kokkos::finalize();
}
