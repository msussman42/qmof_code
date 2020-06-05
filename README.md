# Quadratic Moment of Fluid Interface Construction
## Higher-order Quadrature on Implicitly Defined Domains
The files `saye_algorithm` and `saye_utils` contain objects and methods that perform integration along domains defined implicitly along one or more level set functions. Examples of the usage of these functions is provided in the files `integral2d_test.cpp` and `integral3d_test.cpp`:
```c++
#include <stdio.h>
#include "saye_utils.h"
#include "saye_algorithm.h"

int main()
{
    // Test integral in 2D by integrating various functions over the intersection
    //  of the unit circle and unit square
    
    // Set up initial level set function
    double coeff[6] = {1, 1, 0, 0, 0, -1};

    // Create Phi2D object
    Phi2D phi0( coeff );

    // Initialize lists of functions
    vector<Psi> psis{ Psi( &phi0 ) };
    vector<int> ss{ -1 };

    // Initialize box
    double xL[2] = {0, 0};
    double xU[2] = {1, 1};
    Box U( xL, xU, 2 );

    // Set S to be false because we are computing volume
    bool S = false;

    // Define number of quadrature nodes
    int q = 15;
    
    // Define base function objects
    Unit vof;   // f(x,y) = 1
    Cx cx;      // f(x,y) = x
    Mx2 mx2;    // f(x,y) = x*x

    double volume = I( psis, ss, U, &vof, false, q );
    printf( "int{ }:    %4.10f\n", volume );

    double moment_x = I( psis, ss, U, &cx, false, q );
    printf( "int{ x }:  %4.10f\n", moment_x/volume );
    
    double moment_xx = I( psis, ss, U, &mx2, false, q );
    printf( "int{ xx }: %4.10f\n", moment_xx/volume );
    
    return 0;
}
```
## Parabolic Moment-of-fluid Interface Construction
The file `grad_descent.cpp` contains a number of methods that calculate the optimal 2nd order level set function that best fits a given collection of moment data. This primarily through the functions `quadratic_MOF` and `parabolic_MOF`, the latter of which restricts the possible level sets to only those who generate a parabolic zero level set.
```c++
#include <vector>
#include <stdio.h>
#include "saye_utils.h"
#include "saye_algorithm.h"
#include "grad_descent.h"

int main()
{
    // Test the QMOF algorithm by approximating a circle on the unit square
    //  with a parabola

    int max_iter = 2500;
    char guess_type = 'l';
    double error_tols[4] = {1e-7, 1e-7, 1e-9, 1e-7};

    int NN_data[2] = {0, 0}; // For use with a NN, not needed when initial guess is linear

    double xL[2] = { 0, 0 };
    double xU[2] = { 1, 1 };
    Box U( xL, xU, 2 );

    double refs[6] = { 0.0625, 0.0625, 0.03978873, 0.21220659, 0.21220695, 0.19634954 };
    double init_coeffs[6];
    double final_coeffs[6];

    parabolic_MOF( refs, U, guess_type, error_tols, final_coeffs, max_iter, init_coeffs, NN_data );
    
    printf( "Final coefficients:\n" );
    printf( "%.10f, %.10f, %.10f, %.10f, %.10f, %.10f\n", final_coeffs[0], final_coeffs[1], 
                                                          final_coeffs[2], final_coeffs[3],
                                                          final_coeffs[4], final_coeffs[5] );

    return 0;
}
```
These files can be compiled and executed with the included `Makefile` with minimal additional library support. C++11 is needed for certain functions in the `gradient_descent` file, but these do not impact performance significantly. Otherwise, only the `math` library is needed. 

