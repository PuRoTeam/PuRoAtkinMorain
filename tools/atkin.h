#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <gmp.h>
#include "utility.h"
#include "hilbert.h"
#include "elliptic_curve.h"
#include "factorize.h"

#ifndef ATKIN_H
#define ATKIN_H

#define NUM_PROBABILITY_TEST 10 //numero di test probabilistici di primalità eseguiti con mpz_probab_prime_p

typedef struct curve_orders
{
	mpz_t* orders; //da inizializzare!!!
	int num_elements;	
}curve_orders;

//i parametri a e b stanno in corrispondenza una a un: a[0] con b[0], a[1] con b[1] ecc...
typedef struct curve_parameters
{
	mpz_t* a;
	mpz_t* b;
	int num_elements;
	
}curve_parameters;


void init_curve_orders(mpz_t p, mpz_t D, struct resultset* r, curve_orders* c_orders);

void step_one(mpz_t p, mpz_t D, curve_orders* c_orders);

//restituisce -1 se non è stato possibile fattorizzare m (o è scaduto il timer, ANCORA NON IMPLEMENTATO)
//restituisce  1 se q è probabilmente primo
//restituisce  2 se q è sicuramente primo, e quindi anche p
int step_two(mpz_t m, mpz_t q, mpz_t p, mpz_t D, curve_orders* c_orders);

//Restituisce -1 se bisogna prendere il discriminante successivo
//Restituisce  0 se tutto è andato ok
int step_three(mpz_t p, mpz_t D, curve_parameters* paramSet);

//restituisce -2 se andato in loop infinito (programma termina)
//restituisce -1 se n è composto
//restituisce 0 se tutto ok
int step_four(curve_parameters* out_chosen_parameters, curve_point* P, mpz_t p, curve_parameters c_parameters, curve_orders c_orders, mpz_t m);

//restituisce 0 se n è primo
//restituisce 1 se bisogna andare allo step 4
//restituisce -1 se n è composto
//restituisce 2 se bisogna rieseguire Atkin-Morain
int step_five(mpz_t m, mpz_t q, curve_point P, curve_parameters out_chosen_parameters, mpz_t p);

//restituisce 0 se p è un numero primo
//restituisce -1 se p non è un numero primo
int Atkin_Morain(mpz_t p);

#endif
