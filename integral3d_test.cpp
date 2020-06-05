#include <stdio.h>
#include "saye_utils.h"
#include "saye_algorithm.h"

int main()
{
    // Test integral in 3D by integrating over the intersection
    //  of the unit sphere and unit cube
    
    // Set up initial level set function
    double coeff[10] = {1, 1, 1, 0, 0, 0, 0, 0, 0, -1};

    // Create Phi2D object
    Phi3D phi0( coeff );

    // Initialize lists of functions
    vector<Psi> psis{ Psi( &phi0 ) };
    vector<int> ss{ -1 };

    // Initialize box
    double xL[3] = {0, 0, 0};
    double xU[3] = {1, 1, 1};
    Box U( xL, xU, 3 );

    // Set S to be false because we are computing volume
    bool S = false;

    // Define number of quadrature nodes
    int q = 15;
    
    // Define base function objects
    Mxyz2 f;   // f(x,y) = x*y*z^2

    double value = I( psis, ss, U, &f, false, q );
    printf( "int{ x*y*sin(z) }:    %4.10f\n", value );

    return 0;
    
}

