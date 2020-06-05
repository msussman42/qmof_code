// saye_algorithm.cpp //
// Contains an implementation of Saye's algorithm
#include <stdio.h>
#include <set>
#include <vector>
#include <cmath>
#include "saye_utils.h"
#include "saye_algorithm.h"
using namespace std;

#define MAX_DEPTH 16

double qs[21], ws[21]; 
#include "quadrature.h"

// Evaluate thet maximum of the function |phi(x) - phi(x_c)| for x in U
//  Because the function parabolic, this maximum occurs on one of the vertices of U
double corner_max_phi( const Box& U, Psi psi, Point xc )
{
    int d = U.getD();
    double the_max = -1;
    double test_max;
    double vals[3];
    Point x;

    // Build point using 'vals', in 3D if necessary
    for( int i1 = 0; i1 < 2; i1++ )
    {
        vals[0] = U[0][i1];
        for( int i2 = 0; i2 < 2; i2++ )
        {
            vals[1] = U[1][i2];
            if( d == 3 )
            {
                for( int i3 = 0; i3 < 2; i3++ )
                {
                    vals[2] = U[2][i3];
                    x = Point( vals, d );
                    test_max = dabs( psi.eval(x) - psi.eval(xc) );
                    if( test_max > the_max )
                        the_max = test_max;
                }
            }
            else
            {
                x = Point( vals, d );
                test_max = dabs( psi.eval(x) - psi.eval(xc) );
                if( test_max > the_max )
                    the_max = test_max;
            }
        }
    }

    return the_max;
}

// Find maximum of function |psi_x1(x) - g| for x in U.
//  Because psi_x1 is a hyperplane, the maximum is on a vertex of U
double corner_max_grad( const Box& U, Psi psi, double g, int k )
{
    int d = U.getD();
    double the_max = -1;
    double test_max;
    double vals[3];
    Point x;

    for( int i1 = 0; i1 < 2; i1++ )
    {
        vals[0] = U[0][i1];
        for( int i2 = 0; i2 < 2; i2++ )
        {
            vals[1] = U[1][i2];
            if( d == 3 )
            {
                for( int i3 = 0; i3 < 2; i3++ )
                {
                    vals[2] = U[2][i3];
                    x = Point( vals, d );
                    test_max = dabs( psi.eval_k(x, k) - g );
                    if( test_max > the_max )
                        the_max = test_max;
                }
            }
            else
            {
                x = Point( vals, d );
                test_max = dabs( psi.eval_k(x, k) - g );
                if( test_max > the_max )
                    the_max = test_max;
            }
        }
    }

    return the_max;
}

// Compute the full, Gaussian tensor product of F on U in 1, 2, or 3D
double gauss_tensor_product( F_params * fp, const Box& U, int q )
{
    init_quadrature( q );
    double I = 0;
    int d = U.getD();

    if( d == 1 )
        for( int i1 = 0; i1 < q; i1++ )
        {
            Point x;
            // Build the quadrature node for each dimension
            x.saye( U[0][0] + (U[0][1] - U[0][0])*qs[i1], 0 );
            I += ws[i1]*fp->eval( x );
        }

    if( d == 2 )
        for( int i1 = 0; i1 < q; i1++ )
            for( int i2 = 0; i2 < q; i2++ )
            {
                Point x;
                // Build the quadrature node for each dimension
                x.saye( U[0][0] + (U[0][1] - U[0][0])*qs[i1], 0 );
                x.saye( U[1][0] + (U[1][1] - U[1][0])*qs[i2], 1 );
                I += ws[i1]*ws[i2]*fp->eval( x );
            }

    if( d == 3 )
        for( int i1 = 0; i1 < q; i1++ )
            for( int i2 = 0; i2 < q; i2++ )
                for( int i3 = 0; i3 < q; i3++ )
                {
                    Point x;
                    // Build the quadrature node for each dimension
                    x.saye( U[0][0] + (U[0][1] - U[0][0])*qs[i1], 0 );
                    x.saye( U[1][0] + (U[1][1] - U[1][0])*qs[i2], 1 );
                    x.saye( U[2][0] + (U[2][1] - U[2][0])*qs[i3], 2 );
                    I += ws[i1]*ws[i2]*ws[i3]*fp->eval( x );
                }

    return I;
}

/////////////////////////////////////
// Function to handle the recursive properties of the algorithm
//  When 'eval' is called on an F_params object, it runs through Saye's 
//  algorithm on the lower dimensional function. The possible base level
//  recursive cases are listed in the header file 'base_functions.h'

// Base class declaration to prevent automatic type conversion
F_params::F_params( )
{  }

// Copy constructor
F_params::F_params( const F_params & nfp )
{   
    fp = nfp.fp; 
    psis = nfp.psis;
    ss = nfp.ss;
    x1 = nfp.x1;
    x2 = nfp.x2;
    k = nfp.k;
    q = nfp.q;
    S = nfp.S;
}

// Basic constructor
F_params::F_params( vector<Psi>& n_psis, vector<int>& n_ss, 
                    F_params* n_fp, double n_x1, double n_x2, 
                    int n_k, int n_q, bool n_S )
{
    fp = n_fp; 
    psis = n_psis;
    ss = n_ss;
    x1 = n_x1;
    x2 = n_x2;
    k = n_k;
    q = n_q;
    S = n_S;
}

// Evaluate the function, which changes depending on if the surface integral
//  is being calculated
double F_params::eval( Point x )
{
    if( S == false )
        return F( x, fp, psis, ss, x1, x2, k, q );
    else
        return F_surf( x, fp, psis[0], x1, x2, k );
}

/////////////////////////////////////////

// First function described by Saye (2015), listed in the paper as Algorithm 1.
//  Evaluates the integrand of the integral evaluated in Algorithm 3.
double F( Point x, F_params * fp, vector<Psi>& psis, vector<int>& ss,
          double x1, double x2, int k, int q )
{
    int bad_sign = 0;

    int d = x.getD();
    int n = psis.size(); 
    init_quadrature( q );

    Point xc;
    double I = 0;
    double L;

    set<double> R{ x1, x2 };
    set<double> R_temp;

    // Get roots along [x1, x2] for each psi 
    for( int i = 0; i < n; i++ )
    {
        R_temp = psis[i].roots( x, x1, x2, k );
        // Only use of C++11, used for simplicity of syntax
        for( auto & elem : R_temp )
            R.insert( elem );
    }

    // Turn the set into indexable array
    vector<double> R_list( R.begin(), R.end() );

    // Loop over the partition of [x1, x2] defined by the roots
    for( int j = 0; j < R_list.size() - 1; j++ )
    {
        L = R_list[j+1] - R_list[j];
        xc = x;
        xc.saye( 0.5*( R_list[j+1] + R_list[j] ), k );
        
        // Check if partition is interior to the level-set function
        bad_sign = 0;
        for( int i = 0; i < n; i++ )
        {
            // If any of the psis are exterior, don't add to the integral
            if( ss[i]*psis[i].eval( xc ) < 0 ) 
            {
                bad_sign = 1;
                break;
            }
        }
        if( !bad_sign )
        {
            // Otherwise, evaluate for each quadrature node
            for( int i = 0; i < q; i++ )
            {
                xc = x;
                xc.saye( R_list[j] + L*qs[i], k );
                I += L * ws[i]*fp->eval( xc );
            }
        }
    }

    return I;
}

// Compute the value of Grad( phi(x) )
double gradient( Point x, Phi* phi )
{
    double grad = 0;

    for( int j = 0; j < x.getD(); j++ )
        grad += phi->eval_k( x, j )*phi->eval_k( x, j );

    return sqrt( grad );
}

// Algorithm 2: analog for algorithm 1 for use in surface integrals
double F_surf( Point x, F_params * fp, Psi phi,
          double x1, double x2, int k )
{    
    int d = x.getD();
    // Only one possible root to find
    set<double> R = phi.roots( x, x1, x2, k );
    double r;
    double grad = 0;

    if( R.size() == 1 )
    {
        // Select the first root
        r = *R.begin(); 
        x.saye( r, k );
        // Compute the gradient of phi at x
        for( int j = 0; j < x.getD(); j++ )
            grad += phi.eval_k( x, j )*phi.eval_k( x, j );
        
        // Compute the value of the integrand
        return fp->eval( x ) * sqrt( grad ) / 
                dabs( phi.eval_k( x, k ) );
    }
    else
        return 0;
}

// Algorithm 3: Compute the integral whose domain is interior to the set of level set 
//  functions 'psis', where the "interior" is defined by the sign functions in 'ss'.
//  'fp' represents the integrand of the function; when evaluated, either algorithm 1 or 2
//  is called depending on the value of 'S'. 
double I( vector<Psi> psis, vector<int> ss, const Box& U, 
          F_params * fp, bool S, int q, int depth )
{
    int d = U.getD();
    int n = psis.size();
    int k = 0;

    Point tmp, xc;
    
    vector<Psi> psis_t;
    vector<int> ss_t;

    vector<double> g;
    vector<double> delta_k;
    double jacobian;
    // assert s == false

    Box U1, U2, U3;

    F_params fpt;
    
    // Evaluate the base case scenario
    if( d == 1 )
        return F( tmp, fp, psis, ss, U[0][0], U[0][1], 0, q );

    // Calculate the center of U
    double xcs[3] = {0.5*(U[0][0] + U[0][1]),
                     0.5*(U[1][0] + U[1][1]),
                     0.5*(U[2][0] + U[2][1])};
    xc = Point( xcs, d );

    double delta;

    // Check if the domain U is either entirely bounded by the psis.
    //  This means that the integral is either completely 'full' or 'empty', 
    //  resulting in a simpler integral.
    for( int i = n-1; i >= 0; i-- )
    {
        delta = corner_max_phi( U, psis[i], xc );
        if( dabs( psis[i].eval( xc ) ) >= delta )
            if( ( ss[i]*psis[i].eval( xc ) > 0 ) ||
                ( ss[i]*psis[i].eval( xc ) >= 0 && S == false ) ) // Originally >=
            {
                // If this is true, the function psis[i] is located entirely within U,
                //  and does not need to be evaluated
                n -= 1;
                psis.erase(psis.begin() + i);
                ss.erase(ss.begin() + i);
            }
            else
                // If this is false, the function psis[i] has no portion interior to U,
                //  and the integral evaluates to zero.
                return 0;
    }
    if( n == 0 )
        // If every function psis[i] is "full", perform a tensor product gaussian quadrature
        return U.volume() * gauss_tensor_product( fp, U, q );

    // Evaluate the maximum of |phi_xj(x_c)| among each derivative direction j
    // The direction with the maximum, k, becomes the test height function direction
    for( int j = 0; j < d; j++ )
    {
        if( dabs( psis[0].eval_k( xc, j ) ) > 
            dabs( psis[0].eval_k( xc, k ) ) )
            k = j;
    }

    // Loop over each remaining level set function psi
    for( int i = 0; i < n; i++ )
    {
        g.clear();
        delta_k.clear();
        jacobian = 0;
        for( int j = 0; j < d; j++ )
        {
            g.push_back( psis[i].eval_k( xc, j ) );
            delta_k.push_back( corner_max_grad( U, psis[i], g[j], j) );
            jacobian += (g[j] + delta_k[j])*(g[j] + delta_k[j]);
        }
        jacobian /= (g[k] - delta_k[k])*(g[k] - delta_k[k]);

        // If these conditions are met, k is a good direction for the height function
        if( dabs(g[k]) > delta_k[k]  && jacobian < 20 )
        {
            // Fix the values of Phi along the max and min of U in the height direction,
            //  add them to the collection of level set functions
            psis_t.push_back( Psi( psis[i], U[k][0], k ) );
            ss_t.push_back( sgn( sign(g[k]), ss[i], S, -1 ) );

            psis_t.push_back( Psi( psis[i], U[k][1], k ) );
            ss_t.push_back( sgn( sign(g[k]), ss[i], S, 1 ) );
        }
        // If the conditions are not met, split the domain into 2 and seek a new height
        //  direction for each
        else
        {
            // If the recursive maximum depth is reached, 
            //  use a coarse approximation of the integral
            if( depth >= MAX_DEPTH )
                return U.volume() * fp->eval( xc );
            
            // split the domain, perform the algorithm on each half
            U.split( U1, U2 );
    
            double val1 = I( psis, ss, U1, fp, S, q, depth + 1 );
            double val2 = I( psis, ss, U2, fp, S, q, depth + 1 );
            return val1 + val2;
        }
    }

    // Define a new integrand with F_params, and evaluate it with
    //  the new collection of level sets.
    fpt = F_params( psis, ss, fp, U[k][0], U[k][1], k, q, S );
    
    U.erase( U3, k );
    
    return I( psis_t, ss_t, U3, &fpt, false, q, depth + 1 );
}
