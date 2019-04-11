# Build
Uses CMake (developed using 3.13.3).
Assume that `$LOWFLOW_ROOT` is the path to this directory.
1. Create build directory somewhere (we'll call it $BUILD_ROOT)
2. `cd $BUILD_ROOT`
3. `cmake $LOWFLOW_ROOT`
4. `make`
5. `make install`
  - This installs the executables into `$LOWFLOW_ROOT/bin` (TODO: make install prefix)

# Running
Each program is a different, self-contained variant, and will be run in the same way.

Command Line Options:
+ `-V` : perform verification after experiment to check that experimental output is correct.
  Program returns `255` if verification fails.
+ `-x <N : int>` : Set size of domain in x direction.
+ `-y <N : int>` : Set size of domain in y direction.
+ `-z <N : int>` : Set size of domain in z direction.
+ `-s <N : int>` : Set size of domain in x, y, z directions to the same value.
  `-s 10` is equivalent to `-x 10 -y 10 -z 10`.
+ `-e <epsilon : float>` : Set acceptable error bound used for verification.
  An absolute difference between experimental and verification output grids less than epsilon is considered 'passing'
+ `-S <N : uint>` : Set the seed used to generate grid-seeds.
+ `-h` : Print this message, and exit codes.

## Exit Codes
+ SUCCESS (0) : Successful run.
+ VERIFICATION_FAIL (255) : Verification failed (only if `-V` flag enabled).
+ COMMAND_LINE_ERROR (254) : Invalid command line options.
+ ASSERTION_FAILURE (134) : Assertion caused program abort.
