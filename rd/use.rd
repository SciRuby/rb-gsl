=begin
= Using Ruby/GSL
== Installation
See ((<here|URL:index.html>)).

== Load the library
Put at the head of your scripts,
    
    require("gsl")

== Naming conventions, C Structs and Ruby Classes

Most of the GSL data types, functions or constants are named as (({gsl_xxx})) or (({GSL_XXX})). 
In Ruby/GSL, the prefix (({gsl_})) is replaced by the module identifier (({GSL::})), 
where (({GSL})) is the top level module of Ruby/GSL,
and the Ruby classes are defined for each of the GSL C structs under the (({GSL})) module.
According to the Ruby manner, the name of each class begins with a capital. For example,

  * Struct (({gsl_vector})) ---> Ruby class (({GSL::Vector}))
    * Function call as (({v = gsl_vector_alloc(5);})) ---> Class method (({v = GSL::Vector.alloc(5)}))
    * Function call as (({gsl_vector_set(v, i, 1.5);})) ---> Method (({v.set(i, 1.5)})) or (({v[i] = 1.5}))
  * Constant (({GSL_SUCCESS})) ---> Ruby module constant (({GSL::SUCCESS}))
  * Function (({gsl_sf_bessel_J0(x)})) --->
    * Submodule (({GSL::Sf}))
    * Module function (({GSL::Sf::bessel_J0(x)}))
      * (({GSL::Sf::Bessel::J0(x)})) is also OK, where (({J0(x)})) is a module function of the sub-sub-module (({GSL::Sf::Bessel})).
  * Function (({gsl_linalg_LU_decomp})) --->
    * Module (({GSL::Linalg::LU}))
      * Singleton method (({GSL::Linalg::LU_decomp}))
      * Submodule (({GSL::Linalg::LU}))
        * Singleton method (({GSL::Linalg::LU.decomp}))
    * Method (({GSL::Matrix#LU_decomp}))

== Examples
See the directories "examples/" and "tests/".

Some of the examples use the (({graph})) utility to show the results. The (({graph}))
utility is included in the ((<GNU plotutils|URL:http://www.gnu.org/software/plotutils/plotutils.html>)) package. Windows-cygwin binaries of (({GNU plotutils})) and 
related packages are available from 
((<here|URL:http://rustam.uwp.edu/support>)).

== Modules and Classes
The following is the list of Ruby/GSL modules and classes, <Name> (<Module or Class>)

* GSL (Module) 
  * Complex (Class) 
  * Poly (Class) 
    * Workspace (Class)
    * DividedDifferenceRepresentation (Class)
    * Taylor (Class)
  * Sf (Module) 
    * Result (Class)
  * Block (Class)
    * Int (Class)
    * Byte (Class)
    * Index < Permutation
  * Vector (Class) 
    * View < Vector
    * Complex (Class)
      * View < Vector::Complex
  * Matrix (Class) 
    * View < Matrix (Class)
    * Complex (Class)
      * View < Matrix::Complex
  * Permutation (Class) 
  * Combination (Class) 
  * Linalg (Module) 
    * LU (Module)
    * QR (Module)
    * QRPT (Module)
    * LQ (Module)
    * LQPT (Module)
    * SV (Module)
    * Cholesky (Module)
    * Symmtd (Module)
    * HH (Module)
  * Eigen (Module)
    * EigenValues < Vector
    * EigenVectors < Matrix
    * Symm (Module)
      * Workspace (Class)
    * Symmv (Module)
      * Workspace (Class)
    * Unsymm (Module)
      * Workspace (Class)
    * Unsymmv (Module)
      * Workspace (Class)
    * Herm (Module)
      * Workspace (Class)
    * Hermv (Module)
      * Workspace (Class)
  * FFT (Module) 
    * Complex (Module)
      * PackedArray (Class)
      * Wavetable (Class)
      * Workspace (Class)
    * Real (Module)
      * Wavetable (Class)
      * Workspace (Class)      
    * HalfComplex (Module)
      * Wavetable (Class)
    * Wavetable (Class)
    * WavetableFactor (Class)
    * Workspace (Class)
  * Wavelet (Class)
      * Wavelet2d < Wavelet
  * Function (Class)
  * Function_fdf (Class)
  * Integration (Module)
    * Workspace
    * QAWS_Table
    * QAWO_Table
  * Rng (Class) 
  * QRng (Class)
  * Ran (Module)
  * Stats (Module) 
  * Histogram (Class) 
    * Integral < Histogram
    * Pdf (Class)
  * Histogram2d (Class)
  * N-tuples 
    * SelectFn (Class)
    * ValueFn (Class)
  * Monte (Module)
    * Function (Class)
    * Plain (Class)
    * Miser (Class)
    * Vegas (Class)
  * Siman (Module)
    * Efunc (Class)
    * Step (Class)
    * Metric (Class)
    * Print (Class)
    * Params (Class)
    * Solver (Class)
  * Odeiv (Module)
    * Control (Class)
    * Evolve (Class)
    * System (Class)
    * Solver (Class)
  * Interp (Class)
    * Accel (Class)
  * Spline (Class)
  * Diff (Module) 
  * Deriv (Module)
  * Cheb (Class) 
  * Sum (Module) 
    * Levin_u (Class)
    * Levin_utrunc (Class)
  * Dht (Class) 
  * Root (Module)
    * Solver (Class)
    * FdfSolver (Class)
  * Min(Module) 
    * FMinimizer (Class)
  * MultiRoot (Module)
    * Function (Class)
    * FSolver (Class)
    * Function_fdf (Class)
    * FdfSolver (Class)
  * MultiMin (Module) 
    * Function (Class)
    * FMinimizer (Class)
    * Function_fdf (Class)
    * FdfMinimizer (Class)
  * Fit (Module) 
  * MultiFit (Module)
    * Workspace (Class)
    * Solver (Class)
    * Function_fdf (Class)
    * FdfSolver (Class)
  * CONST (Module) 
    * MKSA (Module)
    * CGSM (Module)
    * NUM (Module)

((<prev|URL:intro.html>))
((<next|URL:ehandling.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
