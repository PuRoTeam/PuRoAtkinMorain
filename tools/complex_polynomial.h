#include "complessi.h"
#include <assert.h>

#ifndef COMPLEX_POLYNOMIAL_H
#define COMPLEX_POLYNOMIAL_H

//chiamare init_complex_polynomial per usarla correttamente
typedef struct complex_polynomial 
{
	complex_z* coefficients;
	int degree;
	int num_elements; //degree + 1
}c_polynomial;

typedef struct complex_fpolynomial 
{
	complex_f* coefficients;
	int degree;
	int num_elements; //degree + 1
}f_polynomial;

void init_complex_polynomial(complex_z* coef, int c_size, int* deg, int d_size, c_polynomial* poly);
void complex_polynomial_mul(c_polynomial* poly_mul, c_polynomial poly_A, c_polynomial poly_B);
void complex_polynomial_real_part(c_polynomial* real_poly, c_polynomial poly);

void init_complex_fpolynomial(complex_f* coef, int c_size, int* deg, int d_size, f_polynomial* poly);
void complex_fpolynomial_mul(f_polynomial* poly_mul, f_polynomial poly_A, f_polynomial poly_B);
void complex_fpolynomial_real_part(f_polynomial* real_poly, f_polynomial poly);
void complex_fpolynomial_copy(f_polynomial* T, f_polynomial temp_poly);
void complex_fpolynomial_clear(f_polynomial* poly);
void complex_fpolynomial_printf(f_polynomial* poly);

#endif