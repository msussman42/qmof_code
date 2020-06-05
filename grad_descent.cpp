// Functions to perform 2D QMOF
#include <stdio.h>
#include <cmath>
#include <random>
#include "saye_utils.h"
#include "saye_algorithm.h"
#include "grad_descent.h"

using namespace std;

// Initialize several global instances of base functions used.
// This is done so that the objects can be added into an array, and looped over easily
Mx2 mx2;
My2 my2;
Mxy mxy;
Cx cx;
Cy cy;

F_params * bases[5];

// Find the cost of a linear level set function using only centroid moment data
// Data is provided in parabola format (alpha, h, k, theta)
// ref is length 5 to flush with quadratic construction
//  ref[5] = centroid x-coordinate
//  ref[4] = centroid y-coordinate
//  ref[3] = volume fraction
double cost_PLIC( double a[], double ref[], Box U )
{
    double the_cost = 0;

    // Convert parabola to coefficients
    double coeffs[6];
    poly_coefficients( coeffs, a );
    
    // Set up integral
    Phi2D phi0( coeffs );
    vector<Psi> psis{ Psi( &phi0 ) };
    vector<int> ss{-1};
    
    int q = 10;
    
    double temp;
    
    temp = ( I( psis, ss, U, &cx, false, q ) / ref[5] - ref[3] );
    the_cost += temp*temp;

    temp = ( I( psis, ss, U, &cy, false, q ) / ref[5] - ref[4] );
    the_cost += temp*temp;
    
    return the_cost;
}

// Evaluate the cost function of PQIC 
//  a (R^6) is coefficients for phi
//  a0 (R^6) is reference data (a0[5] is vof)
double cost_PQIC( double a[], double a0[], Box U )
{
    // Set up array of base functions
    init_bases();
    double the_cost = 0;

    Phi2D phi0( a );
    vector<Psi> psis{ Psi( &phi0 ) };
    vector<int> ss{ -1 };

    int q = 10;
    
    for( int i = 0; i < 5; i++ )
    {
        double temp = I( psis, ss, U, bases[i], false, q ) / a0[5];
        the_cost += ( temp - a0[i] )*( temp - a0[i] );
    }

    return the_cost;
}

// Evaluate the cost function of a parabola by casting it as a quadratic first
double cost_PPIC( double p[], double ref[], Box U )
{
    double a[6];
    poly_coefficients( a, p );

    return cost_PQIC( a, ref, U );
}

// Used in hill climbing on alpha (curvature) parameter
double alpha_forward_fdm( double cn[], double refs[], Box U, double h, double flood_tol )
{
    double original;
    double grad = 0;
    double an[6];

    original = cn[0];

    cn[0] = original + h;
    poly_coefficients( an, cn );
    flooding( an, refs[5], U, flood_tol);
    grad += cost_PQIC( an, refs, U );

    cn[0] = original;
    poly_coefficients( an, cn );
    flooding( an, refs[5], U, flood_tol );
    grad += -cost_PQIC( an, refs, U );
    
    grad = grad / h;
    cn[0] = original;

    return grad;
}

// Used in hill climbing on alpha (curvature) parameter
double alpha_backward_fdm( double cn[], double refs[], Box U, double h, double flood_tol )
{
    double original;
    double grad = 0;
    double an[6];

    original = cn[0];

    cn[0] = original;
    poly_coefficients( an, cn );
    flooding( an, refs[5], U, flood_tol);
    grad += cost_PQIC( an, refs, U );

    cn[0] = original + h;
    poly_coefficients( an, cn );
    flooding( an, refs[5], U, flood_tol );
    grad += -cost_PQIC( an, refs, U );
    
    grad = grad / h;
    cn[0] = original;

    return grad;
}

// Used in hill climbing on alpha (curvature) parameter
double alpha_centered_fdm( double cn[], double refs[], Box U, double h, double flood_tol )
{
    double original;
    double grad = 0;
    double an[6];

    original = cn[0];

    cn[0] = original + h;
    poly_coefficients( an, cn );
    flooding( an, refs[5], U, flood_tol);
    grad += cost_PQIC( an, refs, U );

    cn[0] = original - h;
    poly_coefficients( an, cn );
    flooding( an, refs[5], U, flood_tol );
    grad += -cost_PQIC( an, refs, U );
    
    grad = grad / 2 / h;
    cn[0] = original;

    return grad;
}




// Perform linear moment of fluid on U, using 'refs' as reference data.
//  'error_tols[]' is an array of 4 error tolerances 
//  (initial guess, finite difference precision, flooding algorithm, steepest descent)
void linear_MOF( double refs[], Box U, double error_tols[], double final_coeffs[], int max_iter )
{
    // Transform all data to the unit square
    double unit_refs[6];
    double unit_coeffs[6];
    double unit_init_coeffs[6];
    transform_refs_to_unit( refs, unit_refs, U );

    // IN [0,1]^2 MODE

    // Set initial guess as a straight line with correct volume
    double init_params[4] = {0, 1 - unit_refs[5], 0, 0};

    // If the box is full or empty, set level set manually
    if( unit_refs[5] <= 1e-5 )
    {
        poly_coefficients( unit_coeffs, init_params );
        transform_coeffs_to_nonunit( unit_coeffs, final_coeffs, U );
        final_coeffs[0] = final_coeffs[1] = final_coeffs[2] = 0;
        final_coeffs[3] = final_coeffs[4] = 0; 
        final_coeffs[5] = 1;
        return;
    }
    if( unit_refs[5] >= 1 - 1e-5 )
    {
        poly_coefficients( unit_coeffs, init_params );
        transform_coeffs_to_nonunit( unit_coeffs, final_coeffs, U );
        final_coeffs[0] = final_coeffs[1] = final_coeffs[2] = 0;
        final_coeffs[3] = final_coeffs[4] = 0; 
        final_coeffs[5] = -1;
        return;
    }
    
    // Calculate the optimal straight line with bisection method
    linear_MOF_bisection_unit( unit_refs, error_tols, init_params, max_iter );
    poly_coefficients( unit_coeffs, init_params ); // gets coefficients from parameters
    transform_coeffs_to_nonunit( unit_coeffs, final_coeffs, U ); // Transform coefficients back to original domain
}

// Calculate the optimal straight line MOF using a golden section search method
void linear_MOF_bisection_unit( double unit_refs[], double error_tols[], double final_params[], int max_iter )
{
    double error_tol = error_tols[0];
    double flood_tol = error_tols[2];

    double xL[2] = {0, 0};
    double xU[2] = {1, 1};
    Box UnitU( xL, xU, 2 );
    
    // Get dummy initial values, will be overwritten by flooding algorithm
    double init_params[4] = {0, 0.5, 0.5, 0};

    // Set up parameters and bounds for golden section search
    double R = 0.5 * (sqrt(5) - 1);
    double a = -2*atan(1);      // Overshoot the bounds to avoid unwanted symmetries
    double b = 10*atan(1);
    double d = R * (b - a);

    double x1 = a + d, x2 = b - d;
    
    // After adjusting the angle each time, do the flooding algorithm to maintain volume
    init_params[3] = x1;
    para_flooding( init_params, unit_refs[5], UnitU, flood_tol );
    
    double Cx1 = -cost_PPIC( init_params, unit_refs, UnitU );

    init_params[0] = 0;
    init_params[1] = 0.5;
    init_params[2] = 0.5;
    init_params[3] = x2;
    para_flooding( init_params, unit_refs[5], UnitU, flood_tol );
    double Cx2 = -cost_PPIC( init_params, unit_refs, UnitU );
       
    double interval, theta;
    // Perform the maximum number of golden section iterations (roughly 35 are needed)
    for( int i = 0; i < max_iter; i++ )
    {
        d = R * d;
        if( Cx1 > Cx2 )
        {
            a = x2;
            x2 = x1; Cx2 = Cx1;
            x1 = a + d;
            init_params[0] = 0;
            init_params[1] = 0.5;
            init_params[2] = 0.5;
            init_params[3] = x1;
            para_flooding( init_params, unit_refs[5], UnitU, flood_tol );
            Cx1 = -cost_PPIC( init_params, unit_refs, UnitU );
        }
        else
        {
            b = x1;
            x1 = x2; Cx1 = Cx2;
            x2 = b - d;
            init_params[0] = 0;
            init_params[1] = 0.5;
            init_params[2] = 0.5;
            init_params[3] = x2;
            para_flooding( init_params, unit_refs[5], UnitU, flood_tol );
            Cx2 = -cost_PPIC( init_params, unit_refs, UnitU );
        }
   
        theta = 0.5 * (a + b);
        interval = b - a;
        //printf("%d | interval length: %.6f | theta: %.6f\n", i, interval, theta );
        if( interval < error_tol )
            break;
    }
    
    final_params[0] = 0.0;
    final_params[1] = 0.5;
    final_params[2] = 0.5;
    final_params[3] = 0.5 * (a + b);
    para_flooding( final_params, unit_refs[5], UnitU, flood_tol );
}

// Main function to perform parabolic moment of fluid.
//  refs[]: 6-long array of reference data in order (xx, yy, xy, x, y, volume)
//  U: the original domain along which the MOF is performed
//  guess_type: 'l' for optimal line, 'h' for optimal line with hill climb algorithm, 'p' for parabolic NN guess (under construction)
//  error_tols[]: 4-array of various error tolerances needed
//  final_coeffs[]: array of the final coefficients of the optimal parabola
//  init_coeffs[]: array of the initial coefficients after using method defined by 'guess_type'
//  NN_data[]: 2-array denoting the structure of the neural network to do the initial guess (under construction)
// Reference data is in U, so we must transform it to [0,1]^d first, then transform back to U
void parabolic_MOF( double refs[], Box U, char guess_type, double error_tols[], double final_coeffs[], int max_iter, double init_coeffs[], int NN_data[] )
{
    // Define error tolerances for each
    double guess_tol = error_tols[0];
    double gradient_tol = error_tols[1];
    double flood_tol = error_tols[2];
    double steepest_tol = error_tols[3];

    // Transform coefficients to unit square
    double unit_refs[6];
    double unit_coeffs[6];
    double unit_init_coeffs[6];
    transform_refs_to_unit( refs, unit_refs, U );

    double xL[2] = {0, 0};
    double xU[2] = {1, 1};
    Box UnitU( xL, xU, 2 );

    // IN [0,1]^2 MODE

    double init_params[4] = {0, 1 - unit_refs[5], 0, 0};

    // Dont do any of this if the box is full or empty
    if( unit_refs[5] <= 1e-5 )
    {
        poly_coefficients( unit_coeffs, init_params );
        transform_coeffs_to_nonunit( unit_coeffs, final_coeffs, U );
        final_coeffs[0] = final_coeffs[1] = final_coeffs[2] = 0;
        final_coeffs[3] = final_coeffs[4] = 0; 
        final_coeffs[5] = 1;
        init_coeffs[5] = 1;
        return;
    }
    if( unit_refs[5] >= 1 - 1e-5 )
    {
        poly_coefficients( unit_coeffs, init_params );
        transform_coeffs_to_nonunit( unit_coeffs, final_coeffs, U );
        final_coeffs[0] = final_coeffs[1] = final_coeffs[2] = 0;
        final_coeffs[3] = final_coeffs[4] = 0; 
        final_coeffs[5] = -1;
        init_coeffs[5] = -1;
        return;
    }
    
    double coeff[6];

    double gamma = 1.0, gamma_num, gamma_den;

    vector<double> moment_errs;
    double this_moment_err, this_symdiff, this_curvature;

    // Store coefficients
    double an[4];
    double anp1[4];
    double anm1[4];

    // Store gradients
    double gradn[4];
    double gradnm1[4];

    if( guess_type == 'l' )
    {
        double init_coeffs[6];
        linear_MOF_bisection_unit( unit_refs, error_tols, init_params, max_iter ); // get parameters for best line
        para_flooding( init_params, unit_refs[5], UnitU, flood_tol );
        poly_coefficients( unit_init_coeffs, init_params ); // gets coefficients from parameters
        line_to_param( unit_init_coeffs, init_params ); // centers the vertex in the box
        init_params[3] = fmod( init_params[3], 8*atan(1) ); // shrinks angle
    }

    if( guess_type == 'h' )
    {
        int j;
        double init_coeffs[6];
        linear_MOF_bisection_unit( unit_refs, error_tols, init_params, max_iter ); // get parameters for best line
        para_flooding( init_params, unit_refs[5], UnitU, flood_tol );
        poly_coefficients( unit_init_coeffs, init_params ); // gets coefficients from parameters
        line_to_param( unit_init_coeffs, init_params ); // centers the vertex in the box
        init_params[3] = fmod( init_params[3], 8*atan(1) ); // shrinks angle

        // Set the original parameters so they can be changed
        double op[4] = {0, init_params[1], init_params[2], init_params[3]};
        double best_a, best_cost = 10000;
        
        // Do 3 hill climbs with these as the initial guesses
        double trials[3] = {-5, 0, 5};
        double as[100];
        double grads[100];
        double costs[100];
        
        // Loop over the initial guesses 
        for( int i = 0; i < 3; i++ )
        {
            init_params[0] = as[0] = trials[i];
            init_params[1] = op[1];
            init_params[2] = op[2];
            init_params[3] = op[3];
            para_flooding( init_params, unit_refs[5], UnitU, flood_tol );
            grads[0] = alpha_centered_fdm( init_params, unit_refs, UnitU, gradient_tol, flood_tol );
            costs[0] = cost_PPIC( init_params, unit_refs, UnitU );

            // Do first iteration manually
            init_params[0] = as[1] = as[0] - grads[0];
            init_params[1] = op[1];
            init_params[2] = op[2];
            init_params[3] = op[3];
            para_flooding( init_params, unit_refs[5], UnitU, flood_tol );
            grads[1] = alpha_centered_fdm( init_params, unit_refs, UnitU, gradient_tol, flood_tol );
            costs[1] = cost_PPIC( init_params, unit_refs, UnitU );

            // Fix the search length for the second iteration 
            gamma = 10;

            // Perform remaining 99 iterations of the hill climb
            for( j = 1; j < 99; j++ )
            {
                // Stop if we aren't moving anywhere
                if( grads[j] == grads[j-1] )
                    break;
                
                // Perform steepest descent iteration
                init_params[0] = as[j+1] = as[j] - gamma * grads[j];
                init_params[1] = op[1];
                init_params[2] = op[2];
                init_params[3] = op[3];
                para_flooding( init_params, unit_refs[5], UnitU, flood_tol );
        
                costs[j+1] = cost_PPIC( init_params, unit_refs, UnitU );
                grads[j+1] = alpha_centered_fdm( init_params, unit_refs, UnitU, 1e-6, flood_tol );
                
                // If we didn't descrease the cost, start over with a smaller step size
                if( costs[j+1] > costs[j] )
                {
                    j--;
                    gamma = gamma * 0.5;
                    continue;
                }

                gamma = dabs( (as[j+1] - as[j] ) / ( grads[j+1] - grads[j] ) );
                //printf("%.1f, %d: alpha: %f grad: %f cost: %f\n", trials[i], j+1, as[j+1], grads[j+1], costs[j+1] );
            }

            // Select the best alpha out of the three
            if( costs[j] < best_cost )
            {
                best_cost = costs[j];
                best_a = as[j];
            }
        }
        
        init_params[0] = best_a;
        init_params[1] = op[1];
        init_params[2] = op[2];
        init_params[3] = op[3];
        para_flooding( init_params, unit_refs[5], UnitU, flood_tol );
    }

    /* Calculate the inital guess using the Neural Network.
       Not very effective right now, commented out for repairs
    if( guess_type == 'p' )
    {
        int n_test_data = NN_data[0];
        int n_hidden_nodes = NN_data[1];

        // Save filenames so you dont need to recompute evrytime
        char input_fn[50];
        char output_fn[50]; 
        char centers_fn[50]; 
        char weights_fn[50];

        sprintf( input_fn, "rbf_data2/input%d.csv", n_test_data);
        sprintf( output_fn, "rbf_data2/output%d.csv", n_test_data);
        sprintf( centers_fn, "rbf_data2/centers%d_%d.csv", n_test_data, n_hidden_nodes);
        sprintf( weights_fn, "rbf_data2/weights%d_%d.csv", n_test_data, n_hidden_nodes);
        
        
        get_training_data_para( n_test_data, input_fn, output_fn ); 
        get_centers( n_hidden_nodes, input_fn, output_fn, centers_fn );
        get_weights( input_fn, output_fn, centers_fn, weights_fn ); 
        

        eval_network( unit_refs, init_params, centers_fn, weights_fn );

        para_flooding( init_params, unit_refs[5], UnitU, flood_tol );
    }
    */

    // Store the initial condition coefficients
    poly_coefficients( unit_init_coeffs, init_params );
    transform_coeffs_to_nonunit( unit_init_coeffs, init_coeffs, U );
    
    
    for( int i = 0; i < 4; i++ )
        an[i] = init_params[i];

    // Keep track of the error at each iteration
    this_moment_err = cost_PPIC( an, unit_refs, UnitU );
    moment_errs.push_back( this_moment_err );
    //printf( "0 %.10f | %f %f %f %f | %f\n", this_moment_err, an[0], an[1], an[2], an[3], gamma );

    double break_cond, gamma1, gamma2, old_err;
    int bs;
    for( int i = 1; i < max_iter; i++ )
    {
        old_err = moment_errs.back();

        // Calculate the gradient at the current iteration
        cost_fdm_PPIC( an, unit_refs, UnitU, gradn, gradient_tol, flood_tol );  

        bs = 0;
        gamma2 = 100.0; // fix the initial search length for each backstepping phase
        for( bs = 0; bs < 50; bs++ ) 
        {
            // Do the steepest descent
            for( int j = 0; j < 4; j++ )
                 anp1[j] = an[j] - gamma2*gradn[j];

            para_flooding( anp1, unit_refs[5], UnitU, flood_tol );
        
            this_moment_err = cost_PPIC( anp1, unit_refs, UnitU ); 
            
            //printf( "\t%d step: %.10f | %f %f %f %f\n", i, this_moment_err, anp1[0], anp1[1], anp1[2], anp1[3] );

            // If the cost doesnt decrease, try again with a smaller gamma
            if( this_moment_err > old_err )
                gamma2 = gamma2 * 0.5;
            else
                break;
        }
        
        moment_errs.push_back( this_moment_err );
        printf( "%d %.10f | %f %f %f %f | %f in %d\n", i, this_moment_err, an[0], an[1], an[2], an[3], gamma2, bs );
        //printf("\n");
        
        char energy_str[50];
        char filename_str[50];
        poly_coefficients( coeff, anp1 );

        // Refresh the old information
        for( int j = 0; j < 4; j++ )
        {
            anm1[j] = an[j];
            an[j] = anp1[j];
            gradnm1[j] = gradn[j];
        }
       
        // if the gradient is close to zero, stop iterations
        break_cond = 0;
        for( int j = 0; j < 4; j++ )
            break_cond += gradn[j]*gradn[j];
        if( sqrt( break_cond ) < steepest_tol )
            break;
    }
    
    poly_coefficients( coeff, anp1 );
    transform_coeffs_to_nonunit( coeff, final_coeffs, U );

    poly_coefficients( unit_init_coeffs, init_params );
    transform_coeffs_to_nonunit( unit_init_coeffs, init_coeffs, U );
    
    init_coeffs[0] = anp1[1]*( U[0][1] - U[0][0] ) + U[0][0];
    init_coeffs[1] = anp1[2]*( U[1][1] - U[1][0] ) + U[1][0];

}


// Calculate the optimal second order approximation, but optimize along the coefficients instead of the parameters of the parabola
void quadratic_MOF( double refs[], Box U, char guess_type, double error_tols[], double final_coeffs[], int max_iter, double init_coeffs[] )
{
    // Define error tolerances for each
    double guess_tol = error_tols[0];
    double gradient_tol = error_tols[1];
    double flood_tol = error_tols[2];
    double steepest_tol = error_tols[3];

    // Transform coefficients to unit square
    double unit_refs[6];
    double unit_coeffs[6];
    double unit_init_coeffs[6];
    transform_refs_to_unit( refs, unit_refs, U );

    double xL[2] = {0, 0};
    double xU[2] = {1, 1};
    Box UnitU( xL, xU, 2 );

    // IN [0,1]^2 MODE

    unit_init_coeffs[0] = unit_init_coeffs[1] = unit_init_coeffs[2] = 0;
    unit_init_coeffs[3] = -1;
    unit_init_coeffs[4] = 0;
    unit_init_coeffs[5] = 1 - unit_refs[5];

    // Dont do any of this if the box is full or empty
    if( unit_refs[5] <= 0 )
    {
        transform_coeffs_to_nonunit( unit_coeffs, final_coeffs, U );
        final_coeffs[5] = 1;
        init_coeffs[5] = 1;
        return;
    }
    if( unit_refs[5] >= 1 )
    {
        transform_coeffs_to_nonunit( unit_coeffs, final_coeffs, U );
        final_coeffs[5] = -1;
        init_coeffs[5] = -1;
        return;
    }
    
    double coeff[6];

    double gamma = 0.001, gamma_num, gamma_den;

    vector<double> moment_errs;
    double this_moment_err;

    // Store coefficients
    double an[6];
    double anp1[6];
    double anm1[6];

    // Store gradients
    double gradn[6];
    double gradnm1[6];

    if( guess_type == 'l' )
    {
        double init_params[4] = {0, 1 - unit_refs[5], 0, 0};
        linear_MOF_bisection_unit( unit_refs, error_tols, init_params, max_iter );
        para_flooding( init_params, unit_refs[5], UnitU, flood_tol );
        poly_coefficients( unit_init_coeffs, init_params );
        scaling( unit_init_coeffs );
    }

    /*
    if( guess_type == 'p' )
    {
        int n_test_data = 1000;
        int n_hidden_nodes = 1000;

        // Save filenames so you dont need to recompute evrytime
        char input_fn[50];
        char output_fn[50]; 
        char centers_fn[50]; 
        char weights_fn[50];

        sprintf( input_fn, "rbf_data2/quad_input%d.csv", n_test_data);
        sprintf( output_fn, "rbf_data2/quad_output%d.csv", n_test_data);
        sprintf( centers_fn, "rbf_data2/quad_centers%d_%d.csv", n_test_data, n_hidden_nodes);
        sprintf( weights_fn, "rbf_data2/quad_weights%d_%d.csv", n_test_data, n_hidden_nodes);
        
        
        get_training_data_quad( n_test_data, input_fn, output_fn ); 
        get_centers( n_hidden_nodes, input_fn, output_fn, centers_fn );
        get_weights( input_fn, output_fn, centers_fn, weights_fn ); 
        

        eval_network( unit_refs, unit_init_coeffs, centers_fn, weights_fn );

        flooding( unit_init_coeffs, unit_refs[5], UnitU, flood_tol );
    }
    */

    for( int i = 0; i < 4; i++ )
        an[i] = unit_init_coeffs[i];

    //printf( "0 %f %f | %f %f %f %f %f %f\n", this_moment_err, gamma, an[0], an[1], an[2], an[3], an[4], an[5] );
    this_moment_err = cost_PQIC( an, unit_refs, UnitU );
    moment_errs.push_back( this_moment_err );

    // Calculate first gradient
    cost_fdm_PQIC( an, unit_refs, UnitU, gradn, gradient_tol, flood_tol );  

    // Do the first steepest descent (no backstepping)
    for( int j = 0; j < 6; j++ )
        anp1[j] = an[j] - gamma*gradn[j];
    flooding( anp1, unit_refs[5], UnitU, flood_tol );
    
    //printf( "1 %f %f | %f %f %f %f %f %f\n", this_moment_err, gamma, anp1[0], anp1[1], anp1[2], anp1[3], anp1[4], anp1[5] );
    
    this_moment_err = cost_PQIC( anp1, unit_refs, UnitU );
    moment_errs.push_back( this_moment_err );
    
    // Refresh the old information
    for( int j = 0; j < 6; j++ )
    {
        anm1[j] = an[j];
        an[j] = anp1[j];
        gradnm1[j] = gradn[j];
    }
    
    double break_cond;
    for( int i = 2; i < max_iter; i++ )
    {

        // Calculate gradient
        cost_fdm_PQIC( an, unit_refs, UnitU, gradn, gradient_tol, flood_tol );  

        // Calculate gamma each iteration according to the Barzilai-Borwein method
        gamma_num = gamma_den = 0;
        for( int j = 0; j < 6; j++ )
        {
            gamma_num += ( an[j] - anm1[j] ) * ( gradn[j] - gradnm1[j] );
            gamma_den += ( gradn[j] - gradnm1[j] ) * ( gradn[j] - gradnm1[j] );
        }
        gamma = dabs( gamma_num ) / gamma_den;
        //printf( "%d %f %f | %f %f %f %f %f %f \n", i, this_moment_err, gamma, an[0], an[1], an[2], an[3], an[4], an[5]);

        // Do the steepest descent
        for( int j = 0; j < 6; j++ )
            anp1[j] = an[j] - gamma*gradn[j];
        flooding( anp1, unit_refs[5], UnitU, flood_tol );

        this_moment_err = cost_PQIC( anp1, unit_refs, UnitU ); 
        moment_errs.push_back( this_moment_err );

        // Refresh the old information
        for( int j = 0; j < 4; j++ )
        {
            anm1[j] = an[j];
            an[j] = anp1[j];
            gradnm1[j] = gradn[j];
        }
        
        break_cond = 0;
        for( int j = 0; j < 6; j++ )
            break_cond += gradn[j]*gradn[j];
        if( sqrt( break_cond ) < steepest_tol )
            break;
    }
    
    transform_coeffs_to_nonunit( anp1, final_coeffs, U );
    transform_coeffs_to_nonunit( unit_init_coeffs, init_coeffs, U );
}

// Take partial derivatives using 5 point stencil
void cost_fdm_PQIC( double an[], double refs[], Box U, double grad[], double h, double flood_tol )
{
    double original;

    for( int i = 0; i < 5; i++ )
    {
        grad[i] = 0;
        original = an[i];

        an[i] = original + 2*h; 
        flooding( an, refs[5], U, flood_tol );
        grad[i] += -cost_PQIC( an, refs, U );

        an[i] = original + h;
        flooding( an, refs[5], U, flood_tol );
        grad[i] += 8*cost_PQIC( an, refs, U );
        
        an[i] = original - h;
        flooding( an, refs[5], U, flood_tol );
        grad[i] += -8*cost_PQIC( an, refs, U );
        
        an[i] = original + 2*h;
        flooding( an, refs[5], U, flood_tol );
        grad[i] += cost_PQIC( an, refs, U );

        grad[i] = grad[i] / 12*h;
        an[i] = original;
    }

    grad[5] = 0;
}

// Calculate 5 point centered difference on each of the 4 parameters of a parabola
void cost_fdm_PPIC( double cn[], double refs[], Box U, double grad[], double h, double flood_tol )
{
    double original;
    double an[6];

    for( int i = 0; i < 4; i++ )
    {
        grad[i] = 0;
        original = cn[i];

        cn[i] = original + 2*h;
        poly_coefficients( an, cn );
        flooding( an, refs[5], U, flood_tol);
        grad[i] += -cost_PQIC( an, refs, U );

        cn[i] = original + h;
        poly_coefficients( an, cn );
        flooding( an, refs[5], U, flood_tol );
        grad[i] += 8*cost_PQIC( an, refs, U );
        
        cn[i] = original - h;
        poly_coefficients( an, cn );
        flooding( an, refs[5], U, flood_tol );
        grad[i] += -8*cost_PQIC( an, refs, U );
        
        cn[i] = original + 2*h;
        poly_coefficients( an, cn );
        flooding( an, refs[5], U, flood_tol );
        grad[i] += cost_PQIC( an, refs, U );

        grad[i] /= 12*h;
        cn[i] = original;
    }
}

// Wrapper function for flooding used on parabolas, adjusting the parameters afterwards
void para_flooding( double c[], double vof, Box U, double tol )
{
    double a[6];

    // Get terms in coefficient terms
    poly_coefficients( a, c );
    double new_const = a[5];
    
    // Do flooding on the coefficients
    flooding( a, vof, U, tol );
    double const_diff = a[5] - new_const;

    // Adjust parameters to match coefficient change
    c[1] += const_diff * cos( c[3] );
    c[2] += const_diff * sin( c[3] );
}

// Use the bisection method on the volume of the level set interior to the domain.
// Because the volume is naturally bounded by 0 and the volume of U, find coefficients
//  that give such a shape, and use them as the initial values
void flooding( double a[], double vof, Box U, double tol )
{
    vector<int> ss{ -1 };
    int q = 10;
    Unit unit;

    Phi2D phi0( a );
    vector<Psi> psis{ Psi( &phi0 ) };

    double low, high, end;
    double start = I( psis, ss, U, &unit, false, q );

    double test_vol;

    // get bounds for bisection method
    if( start < vof ) // Need to decrease a5
    {
        high = a[5];
        low = high - 1;
        while( true )
        {
            a[5] = low;
            phi0 = Phi2D( a );
            psis[0] = Psi( &phi0 );
            
            end = I( psis, ss, U, &unit, false, q );
            
            if( end >= vof )
                break;
            else
                low -= 1;
        }
    }
    else // need to increase a5
    {
        low = a[5];
        high = low + 1;
        while( true )
        {
            a[5] = high;
            phi0 = Phi2D( a );
            psis[0] = Psi( &phi0 );

            end = I( psis, ss, U, &unit, false, q );

            if( end <= vof )
                break;
            else
                high += 1;
        }
    }

    // Bisection method part
    test_vol = end;
    while( dabs( test_vol - vof ) > tol )
    {
        a[5] = 0.5*(low + high);
        
        phi0 = Phi2D( a );
        psis[0] = Psi( &phi0 );

        test_vol = I( psis, ss, U, &unit, false, q );
        
        if( test_vol > vof )
            low = a[5];
        else
            high = a[5];
    }
}

// Calculate the curvature of the level set function at a point (x, y)
double calc_curvature( double coeffs[], double x, double y )
{
    double f = 2 * coeffs[0] * x + coeffs[2] * y + coeffs[3];
    double fx = 2 * coeffs[0];
    double fy = coeffs[2];

    double g = 2 * coeffs[1] * y + coeffs[2] * x + coeffs[4];
    double gx = coeffs[2];
    double gy = 2 * coeffs[1];

    double h = sqrt( f * f + g * g );
    double hx = ( f * fx + g * gx ) / h;
    double hy = ( f * fy + g * gy ) / h;

    return ( h * fx - f * hx + h * gy - g * hy ) / h / h;
}

// Convert from parabola parameters to level set coefficients
// c is ( alpha (1/a), h, k, theta )
void poly_coefficients( double a[], double c[] ) 
{
    double sinth = sin(c[3]);
    double costh = cos(c[3]);
    double cosmsin = c[2]*costh - c[1]*sinth;
    
    a[0] =  c[0]*sinth*sinth/4.0;
    a[1] =  c[0]*costh*costh/4.0;
    a[2] = -c[0]*costh*sinth/2.0;
    a[3] =  c[0]/2.0*sinth*cosmsin - costh;
    a[4] = -c[0]/2.0*costh*cosmsin - sinth;
    a[5] =  c[1]*costh + c[2]*sinth + c[0]/4.0*cosmsin*cosmsin;
}

// Convert from level set coefficients known to represent a straight line
//  to parabola parameters. The vertex of the parabola is the center of the line
void line_to_param( double coeff[], double param[] )
{
    double a = coeff[3];
    double b = coeff[4];
    double c = coeff[5];

    double xs[2], ys[2]; 
    int i=0;
    
    if( -c/a > 0 && -c/a < 1 )
    {
        xs[i] = -c/a;
        ys[i] = 0;
        i++;
    }

    if( -c/b >= 0 && -c/b <= 1 )
    {
        xs[i] = 0;
        ys[i] = -c/b;
        i++;
    }

    if( -(c+a)/b >= 0 && -(c+a)/b <= 1 )
    {
        xs[i] = 1;
        ys[i] = -(a+c)/b;
        i++;
    }

    if( -(c+b)/a > 0 && -(c+b)/a < 1)
    {
        xs[i] = -(c+b)/a;
        ys[i] = 1;
        i++;
    }

    param[0] = 0.0;
    param[1] = 0.5 * ( xs[0] + xs[1] );
    param[2] = 0.5 * ( ys[0] + ys[1] );
    param[3] = atan(1)*4 + atan2( b, a );
}

// Scale each value in the array by the value of the maximum
void scaling( double a[] )
{
    double the_max = 0;
    for(int i = 0; i < 6; i++ )
        the_max = ( the_max > dabs(a[i]) ) ? the_max : dabs(a[i]);
    for(int i = 0; i < 6; i++ )
        a[i] = a[i] / the_max;
}

// Given a level set function, find each of its reference data
void get_refs( Box U, double coeffs[], double refs[] )
{
    Phi2D phi0( coeffs );
    vector<Psi> psis{ Psi( &phi0 ) };
    vector<int> ss{ -1 };

    bool S = false;
    int q = 15;

    Unit unit;
    refs[5] = I( psis, ss, U, &unit, S, q );

    if( refs[5] <= 0 || refs[5] >= U.volume() )
    {
        refs[0] = refs[1] = refs[2] = refs[3] = refs[4];
        return;
    }

    Mx2 mx2;
    refs[0] = I( psis, ss, U, &mx2, S, q ) / refs[5]; 
    
    My2 my2;
    refs[1] = I( psis, ss, U, &my2, S, q ) / refs[5]; 
    
    Mxy mxy;
    refs[2] = I( psis, ss, U, &mxy, S, q ) / refs[5]; 
    
    Cx cx;
    refs[3] = I( psis, ss, U, &cx, S, q ) / refs[5]; 
    
    Cy cy;
    refs[4] = I( psis, ss, U, &cy, S, q ) / refs[5]; 
}

// transform reference data from an arbitrary domain to the unit square
void transform_refs_to_unit( double refs[], double refs_t[], Box U )
{
    double a = U[0][0];
    double b = U[0][1];
    double c = U[1][0];
    double d = U[1][1];

    // vof
    refs_t[5] = refs[5] / (b-a) / (d-c);

    // Cx
    refs_t[3] = ( refs[3] - a ) / ( b - a );

    // Cy
    refs_t[4] = ( refs[4] - c ) / ( d - c );

    // Mxx
    refs_t[0] = ( refs[0] - 2*a*refs[3] + a*a ) / (b-a) / (b-a);
    
    // Myy
    refs_t[1] = ( refs[1] - 2*c*refs[4] + c*c ) / (d-c) / (d-c);

    // Mxy
    refs_t[2] = ( refs[2] - c*refs[3] - a*refs[4] + a*c ) / (b-a) / (d-c);
}

// transform reference data from the unit square to an arbitrary domain
void transform_coeffs_to_nonunit( double coeffs[], double coeffs_t[], Box U )
{
    double a = U[0][0];
    double b = U[0][1];
    double c = U[1][0];
    double d = U[1][1];

    coeffs_t[0] = coeffs[0] / (b-a) / (b-a);

    coeffs_t[1] = coeffs[1] / (d-c) / (d-c);

    coeffs_t[2] = coeffs[2] / (b-a) / (d-c);

    coeffs_t[3] = -2*coeffs[0]*a / (b-a)/(b-a) - coeffs[2]*c / (b-a)/(d-c) + coeffs[3] / (b-a);

    coeffs_t[4] = -2*coeffs[1]*c / (d-c)/(d-c) - coeffs[2]*a / (b-a)/(d-c) + coeffs[4] / (d-c);

    coeffs_t[5] = coeffs[0]*a*a/(b-a)/(b-a) + coeffs[1]*c*c/(d-c)/(d-c) + 
                  coeffs[2]*a*c/(b-a)/(d-c) - coeffs[3]*a/(b-a) - coeffs[4]*c/(d-c) + coeffs[5];

}

// Calculate the symmetric difference between two level set functions
double quad_symdiff( double a[], double b[], Box U, int q )
{
    Phi2D phia( a );
    Phi2D phib( b );

    vector<Psi> psis{ Psi( &phia ), Psi( &phib ) };
    vector<int> ssa{ -1, 1 };
    vector<int> ssb{ 1, -1 };

    bool S = false;

    Unit vof;

    return I( psis, ssa, U, &vof, S, q ) + I( psis, ssb, U, &vof, S, q );
}


void init_bases()
{
    bases[0] = &mx2;
    bases[1] = &my2;
    bases[2] = &mxy;
    bases[3] = &cx;
    bases[4] = &cy;
}

