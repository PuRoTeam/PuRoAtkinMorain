#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>

#ifndef COMPLESSI_H
#define COMPLESSI_H

//#define PRECISION 10000
#define PRECISION 10000

typedef struct complex_z { 
    mpz_t real; 
	mpz_t img; 
}complex_z;

typedef struct complex_f { 
    mpf_t real; 
	mpf_t img; 
}complex_f;

//Init struct complex_z
void complex_zinit(struct complex_z* c);

//Init struct complex_f
void complex_finit(struct complex_f* c);

//Clear strcut complex_f
void complex_fclear(struct complex_f* c);

//Add two complex and return a struct complex
void add_zcomplex(struct complex_z* z3, struct complex_z z1, struct complex_z z2 );

//Add two complex and return a struct complex
void add_fcomplex(struct complex_f* result, struct complex_f z1, struct complex_f z2 );

//Mul two complex and return a struct complex
void mul_zcomplex(struct complex_z* z3, struct complex_z z1, struct complex_z z2);

void mul_fcomplex(struct complex_f* result, struct complex_f z1, struct complex_f z2);

void div_fcomplex(struct complex_f* result, struct complex_f z1, struct complex_f z2);

void pow_no_bin(struct complex_f* result, struct complex_f z1, int exp);

#endif