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
