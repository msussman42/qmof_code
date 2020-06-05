// File contains several utility classes that assist in the implementation
//  of Saye's algorithm
#include <stdio.h>
#include <set>
#include <cmath>
#include "saye_utils.h"
using namespace std;

// The 'Point' class is a maximally 3-long array of values.
//  Uninitialized values in a 'Point' are set to -1, and are marked
//  by the use of a 'd' variable that stores the current length
Point::Point()
{
	for( int i = 0; i < 3; i++ )
		vals[i] = -1;
	d = 0;
}

Point::Point( double xs[], int n )
{
	for( int i = 0; i < n; i++ )
		vals[i] = xs[i];
	for( int i = n; i < 3; i++ )
		vals[i] = -1;
	d = n;
}

Point::Point( const Point & x )
{
	for( int i = 0; i < 3; i++ )
		vals[i] = x.vals[i];
	d = x.d;
}

// Named for the similar procedure outlined in Saye's paper, 
//  where the notation x + ye_k denotes the point x with the value y 
//  inserted before the kth index
void Point::saye( double val, int k )
{
	for( int i = d; i > k; i-- )
		vals[i] = vals[i-1];
	vals[k] = val;
	d++;
}

double& Point::operator[]( int i )
{	return vals[i];	  }

Point& Point::operator=( const Point x )
{
	for( int i = 0; i < 3; i++ )
		vals[i] = x.vals[i];
	d = x.d;
    
    return *this;
}

int Point::getD()
{   return d;   }

///////////////////////////////////////

// The 'Box' class stores an array of ordered pairs denoting the max and min value in 
//  each of d directions. Contains a number of methods to facilitate the operations
//  of Saye's algorithm, such as splitting the box along its maximum length axis.

Box::Box()
{	d = 0;   }

Box::Box( const Box & U )
{
    d = U.d;
    for( int i = 0; i < d; i++ )
    {
        x[i][0] = U.x[i][0];
        x[i][1] = U.x[i][1];
    }
}

Box::Box( double nxL[], double nxU[], int nd )
{
	for( int i = 0; i < nd; i++ )
	{
		x[i][0] = nxL[i];
		x[i][1] = nxU[i];
	}

	d = nd;
}

const double* Box::operator[]( int dir ) const
{	return x[dir];   }

void Box::split( Box& U1, Box& U2 ) const
{
	int sd = 0;
	double x1L[3], x1U[3];
	double x2L[3], x2U[3];

	// Get the direction to split along
	for( int i = 0; i < d; i++ )
	{
		if( x[i][1] - x[i][0] > x[sd][1] - x[sd][0] )
			sd = i;
	}

	for( int i = 0; i < sd; i++ )
	{
		x1L[i] = x2L[i] = x[i][0];
		x1U[i] = x2U[i] = x[i][1];
	}
	
	x1L[sd] = x[sd][0];
	x1U[sd] = 0.5*( x[sd][0] + x[sd][1] );

	x2L[sd] = 0.5*( x[sd][0] + x[sd][1] );
	x2U[sd] = x[sd][1];

	for( int i = sd+1; i < d; i++ )
	{
		x1L[i] = x2L[i] = x[i][0];
		x1U[i] = x2U[i] = x[i][1];
	}
	
	U1 = Box( x1L, x1U, d );
	U2 = Box( x2L, x2U, d );
}

int Box::getD() const
{   return d;   }

double Box::volume() const
{
    double vol = 1;
    for( int i = 0; i < d; i++ )
        vol *= x[i][1] - x[i][0];
    return vol;
}

void Box::erase( Box& U, int k ) const
{
	double xL[3], xU[3];
    int nd = 0;

    for( int i = 0; i < k; i++ )
    {
        xL[nd] = x[i][0];
        xU[nd] = x[i][1];
        nd++;
    }

    for( int i = k + 1; i < d; i++ )
    {
        xL[nd] = x[i][0];
        xU[nd] = x[i][1];
        nd++;
    }

    U = Box( xL, xU, d - 1 );
}

//////////////////////////////////////
// Phi2D and Phi3D both inherit from the same abstract class Phi,
//  which contains methods for evaluation, root finding along an axis, and 
//  evaluation of partial derivatives, but with different definitions for the 
//  2 different sets of dimensions. Higher degree level sets can be accomodated 
//  by adding in the appropriate Phi4D class, etc.
// 'Phi2D's contain an array of 6 values as coefficients
Phi2D::Phi2D()
{ }

Phi2D::Phi2D( double nCoeff[] )
{
	for( int i = 0; i < 6; i++ )
		coeff[i] = nCoeff[i];
}

Phi2D::Phi2D( const Phi2D & phi )
{
    for( int i = 0; i < 6; i++ )
		coeff[i] = phi.coeff[i];

}

double Phi2D::eval( Point x )
{
    return coeff[0]*x[0]*x[0] + coeff[1]*x[1]*x[1] + 
	       coeff[2]*x[0]*x[1] + coeff[3]*x[0]      + 
           coeff[4]*x[1]      + coeff[5];
}

double Phi2D::eval_k( Point x, int k )
{
	if( k == 0 )
        return 2*coeff[0]*x[0] + coeff[2]*x[1] + 
				 coeff[3];
    if( k == 1 )
        return 2*coeff[1]*x[1] + coeff[2]*x[0] + 
				 coeff[4];
    return 99999;
}

set<double> Phi2D::roots( Point x, double x1, double x2, int k )
{
    double a, b, c;

	if(k == 0)
	{
		a = coeff[0];
        b = coeff[2]*x[1] + coeff[3];
        c = coeff[1]*x[1]*x[1] + coeff[4]*x[1] + 
		    coeff[5];
	}
	if(k == 1)
	{
		a = coeff[1];
        b = coeff[2]*x[0] + coeff[4];
        c = coeff[0]*x[0]*x[0] + coeff[3]*x[0] + 
		    coeff[5];
	}

	return quadratic( a, b, c, x1, x2 );
}


///////////////////////////////////////
// 'Phi3D's contain an array of 10 values as coefficients.
//  Otherwise, its definitions are analogous to the 2D case
Phi3D::Phi3D()
{ }

Phi3D::Phi3D( double nCoeff[] )
{
	for( int i = 0; i < 10; i++ )
		coeff[i] = nCoeff[i];
}

Phi3D::Phi3D( const Phi3D & phi )
{
    for( int i = 0; i < 10; i++ )
		coeff[i] = phi.coeff[i];

}

double Phi3D::eval( Point x )
{
    return coeff[0]*x[0]*x[0] + coeff[1]*x[1]*x[1] + 
	       coeff[2]*x[2]*x[2] + coeff[3]*x[0]*x[1] + 
		   coeff[4]*x[1]*x[2] + coeff[5]*x[0]*x[2] + 
           coeff[6]*x[0]      + coeff[7]*x[1]      + 
		   coeff[8]*x[2]      + coeff[9];
}

double Phi3D::eval_k( Point x, int k )
{
	if( k == 0 )
        return 2*coeff[0]*x[0] + coeff[3]*x[1] + 
				 coeff[5]*x[2] + coeff[6];
    if( k == 1 )
        return 2*coeff[1]*x[1] + coeff[3]*x[0] + 
				 coeff[4]*x[2] + coeff[7];
    if( k == 2 )
        return 2*coeff[2]*x[2] + coeff[4]*x[1] + 
				 coeff[5]*x[0] + coeff[8];
    return 99999;
}

set<double> Phi3D::roots( Point x, double x1, double x2, int k )
{
    double a, b, c;

	if(k == 0)
	{
		a = coeff[0];
        b = coeff[3]*x[1] + coeff[5]*x[2] + coeff[6];
        c = coeff[1]*x[1]*x[1] + coeff[2]*x[2]*x[2] + 
		    coeff[4]*x[1]*x[2] + coeff[7]*x[1] + 
			coeff[8]*x[2] + coeff[9];
	}
	if(k == 1)
	{
		a = coeff[1];
        b = coeff[3]*x[0] + coeff[4]*x[2] + coeff[7];
        c = coeff[0]*x[0]*x[0] + coeff[2]*x[2]*x[2] + 
		    coeff[5]*x[0]*x[2] + coeff[6]*x[0] + 
			coeff[8]*x[2] + coeff[9];
	}
	if(k == 2)
	{
		a = coeff[2];
        b = coeff[4]*x[1] + coeff[5]*x[0] + coeff[8];
        c = coeff[0]*x[0]*x[0] + coeff[1]*x[1]*x[1] + 
		    coeff[3]*x[0]*x[1] + coeff[6]*x[0] + 
			coeff[7]*x[1] + coeff[9];
	}

	return quadratic( a, b, c, x1, x2 );
}


////////////////////////////////////////////////
// Psi represents, mathematically, a Phi object with certain indexes fixed to certain values.
//  This is accomplished by holding a pointer to a Phi object (so that Phi2D or Phi3D objects
//  can be used interchangeably) that stores the proper coefficients. When such a Psi object
//  is evaluated, the full coordinates are reconstructed and passed to the Phi object
// Two arrays, 'dirs' and 'vals' are used to maintain these fixed values
Psi::Psi()
{ }

Psi::Psi( const Psi & psi )
{
    for( int i = 0; i < 2; i++ )
    {
        vals[i] = psi.vals[i];
        dirs[i] = psi.dirs[i];
    }
    
    phi0 = psi.phi0;
    d = psi.d;
}

// A dir value of -1 indicates a lack of fixing
Psi::Psi( Phi * phi )
{
    phi0 = phi;
	for( int i = 0; i < 2; i++ )
		vals[i] = dirs[i] = -1;
	d = 0;
}


Psi::Psi( Psi psi, double val, int dir )
{
	phi0 = psi.phi0;
    for( int i = 0; i < psi.getD(); i++ )
	{
		vals[i] = psi.getVal(i);
		dirs[i] = psi.getDir(i);
	}

	vals[psi.getD()] = val;
	dirs[psi.getD()] = dir;

	d = psi.getD() + 1;
}

// Original coordinates are reconstructed in reverse order
//  because directions can be repeated. For example, a Psi object with
//  dirs = [0, 0] and vals = [5, 7] evaluated at 3 will be reconstructed as 
//  evaluating Phi( 5, 7, 3 ).
double Psi::eval( Point x )
{
	for( int i = d - 1; i >= 0; i-- )
		x.saye( vals[i], dirs[i] );

	return phi0->eval( x );
}

double Psi::eval_k( Point x, int k )
{
	for( int i = d - 1; i >= 0; i-- )
	{
		if( dirs[i] <= k )
			k++;
		x.saye( vals[i], dirs[i] );
	}

	return phi0->eval_k( x, k );
    
}

set<double> Psi::roots( Point x, double x1, double x2, int k )
{
    // Need to insert a dummy value into the point to fix the indices
	x.saye( -999, k );
	for( int i = d - 1; i >= 0; i-- )
	{
		if( dirs[i] <= k )
			k++;
		x.saye( vals[i], dirs[i] );
	}

    return phi0->roots( x, x1, x2, k );
}


int Psi::getDir( int i )
{    return dirs[i];    }

double Psi::getVal( int i )
{   return vals[i];     }

int Psi::getD()
{   return d;   }

Phi* Psi::getPhi()
{   return phi0;    }

////////////////////////////////////////////////

// A stable quadrativ formula that returns the roots of a quadratic as a set.
set<double> quadratic( double a, double b, double c, double x1, double x2 )
{
    double meps = 1e-7; // "Machine epsilon"
	set<double> roots;
	double d = b*b - 4*a*c, r1, r2;

	if( d < 0 || ( dabs( a - b ) < meps && dabs( b ) < meps ) )
		return roots;
	else if( dabs( a ) < meps )
	{
		r1 = -c / b;
		r2 = r1;
	}
	else if( b > 0 )
	{
		r1 = ( -b - sqrt(d) )/2/a;
		r2 = 2 * c / ( -b - sqrt(d) );
	}
    else
    {
		r1 = ( -b + sqrt(d) )/2/a;
		r2 = 2 * c / ( -b + sqrt(d) );
    }

	if( x1 < r1 && r1 < x2 )
		roots.insert( r1 );
	if( x1 < r2 && r2 < x2 )
		roots.insert( r2 );
	
	return roots;
}

// The sign function defined by Saye
double sgn( int m, int s, bool S, int sigma )
{
    if( m == sigma*s || S == true )
        return sigma*m;
    return 0;
}

double sign( double d )
{
    if( d > 0 )
        return 1;
    else if( d < 0 )
        return -1;
    return 0;
}

// Helper print function
void print_arr( double a[], int N )
{
    for( int i = 0; i < N; i++ )
        printf("%.4f, ", a[i]);

    printf("\b\b\n");
}

// A double precision absolute value function
double dabs( double d )
{
    if( d < 0 )
        return -d;
    return d;
}



