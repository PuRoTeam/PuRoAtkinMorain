#include <gmp.h>
#include <math.h>
#include <complex.h>
#include "deltaprecision.h"
#include "complex_polynomial.h"
#include "polynomial.h"
#include "trigonometry.h"

#ifndef HILBERT_H
#define HILBERT_H
 
#define SET_SIZE 1000
#define CYCLES 1

struct abc{
	mpz_t a; 
	mpz_t b;
	mpz_t c;
};

typedef struct hilbert_result{
	mpz_t number_class;				//Numero di class h(D)
	r_big_polynomial Hpoly;				//Polinomio di classe Hilbert, con grado h(D)
	struct abc abc_set[SET_SIZE];	//Insieme di triple (a,b,c)
	int count;						//Numero di elementi nell'insieme
} Hresult;

/*Dato un discriminante D(negativo), questa funzione ritorna
 *una combinazione di:
 * -numero di classe h(D)
 * -il polinomio di classe Hilbert, con grado h(D)
 * -l'insieme delle forme ridotte (a,b,c) del discriminante D
*/
void hilbertclass(Hresult* res, mpz_t D, mpz_t n);

void exponential(complex_f* result, complex_f z);

void other_exponential(complex_f* result, complex_f z);

#endif
