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

