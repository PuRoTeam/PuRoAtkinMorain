#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <signal.h>
#include <assert.h>
#include <unistd.h>

#ifndef FACTORIZE_H
#define FACTORIZE_H

#define FLAG_VERBOSE 0

#define WAITING_TIME 5 //numero di secondi di attesa prima che scada il timer

#define AUTO_TIMING 0 //se 1 il timer viene impostato in base alla grandezza del numero, se 0 viene preso il valore della variabile WAITING_TIME [CONSIGLIATO]

typedef struct prime_factors
{
	mpz_t* factors;
	int num_elements;	
}prime_factors;

void catch_alarm(int sig);

void factor_using_division (prime_factors* p_factors, mpz_t t, unsigned int limit);

void factor_using_division_2kp (prime_factors* p_factors, mpz_t t, unsigned int limit, unsigned long p);

void factor_using_pollard_rho (prime_factors* p_factors, mpz_t n, unsigned long a, unsigned long p);

//per fattorizzare richiamare factor (che richiama le altre funzioni):
//p_factor deve essere un array di almeno 30/50 elementi tanto non credo vengano generati più fattori (per evitare realloc) - num_elements deve essere posto a zero!
//t è il numero da fattorizzare
//p deve esser posto uguale a zero
//restituisce 0 se l'algoritmo è terminato e l'intero è stato totalmente fattorizzato
//restituisce -1 se il timer è scaduto e l'intero è stato fattorizzato in parte o per niente fattorizzato
int factor (prime_factors* p_factors, mpz_t t, unsigned long p);

#endif
