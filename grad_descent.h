#ifndef GRADDESC
#define GRADDESC

#include "saye_utils.h"

double cost_PLIC( double a[], double ref[], Box U ); //
double cost_PQIC( double a[], double a0[], Box U ); //
double cost_PPIC( double p[], double ref[], Box U ); //

double alpha_forward_fdm( double cn[], double refs[], Box U, double h, double flood_tol );
double alpha_backward_fdm( double cn[], double refs[], Box U, double h, double flood_tol );
double alpha_centered_fdm( double cn[], double refs[], Box U, double h, double flood_tol );

void cost_fdm_PPIC( double cn[], double refs[], Box U, double grad[], double grad_tol, double flood_tol ); //

void cost_fdm_PLIC( double an[], double refs[], Box U, double grad[], double grad_tol, double flood_tol );

void cost_fdm_PQIC( double an[], double refs[], Box U, double grad[], double grad_tol, double flood_tol ); //

void flooding( double a[], double vof, Box U, double tol );
void para_flooding( double c[], double vof, Box U, double tol );
void init_bases();

void linear_MOF( double refs[], Box U, double error_tols[], double final_coeffs[], int max_iter );
void linear_MOF_bisection_unit( double unit_ref[], double error_tols[], double final_coeffs[], int max_iter );

void parabolic_MOF( double ref[], Box U, char guess_type, double error_tols[], double final_coeffs[], int max_iter, double init_coeffs[], int NN_data[] );
void quadratic_MOF( double ref[], Box U, char guess_type, double error_tols[], double final_coeffs[], int max_iter, double init_coeffs[] );

void get_refs( Box U, double coeffs[], double refs[] );
void transform_refs_to_unit( double refs[], double refs_t[], Box U );
void transform_coeffs_to_nonunit( double coeffs[], double coeffs_t[], Box U );

double quad_symdiff( double a[], double b[], Box U, int h );
double calc_curvature( double coeffs[], double x, double y );

void poly_coefficients( double a[], double c[] );
void line_to_param( double coeff[], double param[] );
void scaling( double a[] );
#endif

