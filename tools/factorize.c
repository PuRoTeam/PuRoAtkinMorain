/* Factoring with Pollard's rho method.

Copyright 1995, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2005, 2009
Free Software Foundation, Inc.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see http://www.gnu.org/licenses/.  */

#include "factorize.h"

volatile sig_atomic_t keep_going = 1;

unsigned int TIMER;

static unsigned add[] = {4, 2, 4, 2, 4, 6, 2, 6};

void catch_alarm(int sig)
{
	signal(sig, SIG_IGN); //per portabilità
	keep_going = 0;
	signal(sig, catch_alarm); //per portabilità
}

void factor_using_division (prime_factors* p_factors, mpz_t t, unsigned int limit)
{
	mpz_t q, r;
	unsigned long int f;
	int ai;
	unsigned *addv = add;
	unsigned int failures;

	if (FLAG_VERBOSE > 0)
	{
		printf ("[trial division (%u)] ", limit);
		fflush (stdout);
	}

	mpz_init (q);
	mpz_init (r);

	f = mpz_scan1 (t, 0);
	mpz_div_2exp (t, t, f);
	while (f)
	{
		//printf ("2 ");     
		//fflush (stdout);

		mpz_set_ui(p_factors->factors[p_factors->num_elements], 2);
		p_factors->num_elements++;

		--f;
	}

	//for (;;)
	while(keep_going)	
	{
		mpz_tdiv_qr_ui (q, r, t, 3);
		if (mpz_cmp_ui (r, 0) != 0)
			break;
		mpz_set (t, q);
		//printf ("3 ");
		//fflush (stdout);

		mpz_set_ui(p_factors->factors[p_factors->num_elements], 3);
		p_factors->num_elements++;
	}

	if(!keep_going)
	{
		mpz_clear(q);
		mpz_clear(r);
		return;
	}

	//for (;;)
	while(keep_going)
	{
		mpz_tdiv_qr_ui (q, r, t, 5);
		if (mpz_cmp_ui (r, 0) != 0)
			break;
		mpz_set (t, q);
		//printf ("5 ");
		//fflush (stdout);

		mpz_set_ui(p_factors->factors[p_factors->num_elements], 5);
		p_factors->num_elements++;
	}

	if(!keep_going)
	{
		mpz_clear(q);
		mpz_clear(r);
		return;
	}

	failures = 0;
	f = 7;
	ai = 0;
	while (mpz_cmp_ui (t, 1) != 0)
	{
		mpz_tdiv_qr_ui (q, r, t, f);
		if (mpz_cmp_ui (r, 0) != 0)
		{
			f += addv[ai];
			if (mpz_cmp_ui (q, f) < 0)
				break;
			ai = (ai + 1) & 7;
			failures++;
			if (failures > limit)
				break;
		}
		else
		{
			mpz_swap (t, q);
			//printf ("%lu ", f);
			//fflush (stdout);

			mpz_set_ui(p_factors->factors[p_factors->num_elements], f);
			p_factors->num_elements++;

			failures = 0;
		}
	}

	mpz_clear(q);
	mpz_clear(r);
}

void factor_using_division_2kp (prime_factors* p_factors, mpz_t t, unsigned int limit, unsigned long p)
{
	mpz_t r;
	mpz_t f;
	unsigned int k;

	if (FLAG_VERBOSE > 0)
	{
		printf ("[trial division (%u)] ", limit);
		fflush (stdout);
	}

	mpz_init (r);
	mpz_init_set_ui (f, 2 * p);
	mpz_add_ui (f, f, 1);
	for (k = 1; k < limit; k++)
	{
		mpz_tdiv_r (r, t, f);
		while (mpz_cmp_ui (r, 0) == 0)
		{
			mpz_tdiv_q (t, t, f);
			mpz_tdiv_r (r, t, f);
			//mpz_out_str (stdout, 10, f);
			//fflush (stdout);
			//fputc (' ', stdout);

			mpz_set(p_factors->factors[p_factors->num_elements], f);
			p_factors->num_elements++;
		}
		mpz_add_ui (f, f, 2 * p);
	}

	mpz_clear(f);
	mpz_clear(r);
}

void factor_using_pollard_rho (prime_factors* p_factors, mpz_t n, unsigned long a, unsigned long p)
{
	mpz_t x, x1, y, P;
	mpz_t t1, t2;
	unsigned long long k, l, i;

	if(FLAG_VERBOSE > 0)
	{
		printf("[pollard-rho (%lu)] ", a);
		fflush(stdout);
	}

	mpz_init(t1);
	mpz_init(t2);

	mpz_init_set_si (y, 2);
	mpz_init_set_si (x, 2);
	mpz_init_set_si (x1, 2);
	mpz_init_set_ui (P, 1);
	k = 1;
	l = 1;

	while (mpz_cmp_ui (n, 1) != 0)
	{
		//for (;;)
		while(keep_going)
		{
			do
			{
				if (p != 0)
				{
					mpz_powm_ui (x, x, p, n);
					mpz_add_ui (x, x, a);
				}
				else
				{
					mpz_mul (t1, x, x);
					mpz_mod (x, t1, n);
					mpz_add_ui (x, x, a);
				}

				mpz_sub (t1, x1, x);
				mpz_mul (t2, P, t1);
				mpz_mod (P, t2, n);

				if (k % 32 == 1)
				{
					mpz_gcd (t1, P, n);
					if (mpz_cmp_ui (t1, 1) != 0)
						goto factor_found;
					mpz_set (y, x);
				}
			}
			while (--k != 0);

			mpz_gcd (t1, P, n);
			if (mpz_cmp_ui (t1, 1) != 0)
				goto factor_found;

			mpz_set (x1, x);
			k = l;
			l = 2 * l;
			for (i = 0; i < k; i++)
			{
				if (p != 0)
				{
					mpz_powm_ui (x, x, p, n);
					mpz_add_ui (x, x, a);
				}
				else
				{
					mpz_mul (t1, x, x);
					mpz_mod (x, t1, n);
					mpz_add_ui (x, x, a);
				}
			}
			mpz_set (y, x);
		}

		if(!keep_going)
		{
			mpz_clear(P);
			mpz_clear(t2);
			mpz_clear(t1);
			mpz_clear(x1);
			mpz_clear(x);
			mpz_clear(y);
			return;
		}		

		factor_found:
		do
		{
			if (p != 0)
			{
				mpz_powm_ui (y, y, p, n); mpz_add_ui (y, y, a);
			}
			else
			{
				mpz_mul (t1, y, y);
				mpz_mod (y, t1, n);
				mpz_add_ui (y, y, a);
			}
			mpz_sub (t1, x1, y);
			mpz_gcd (t1, t1, n);
		}
		while (mpz_cmp_ui (t1, 1) == 0);

		mpz_divexact (n, n, t1);	/* divide by t1, before t1 is overwritten */

		if (!mpz_probab_prime_p (t1, 25))
		{
			do
			{
				mp_limb_t a_limb;
				mpn_random (&a_limb, (mp_size_t) 1);
				a = a_limb;
			}
			while (a == 0);

			if (FLAG_VERBOSE > 0)
			{
				printf ("[composite factor--restarting pollard-rho] ");
				fflush (stdout);
			}
			factor_using_pollard_rho (p_factors, t1, a, p);
		}
		else
		{
			//mpz_out_str (stdout, 10, t1);
			//fflush (stdout);
			//fputc (' ', stdout);

			mpz_set(p_factors->factors[p_factors->num_elements], t1);
			p_factors->num_elements++;
		}
		mpz_mod (x, x, n);
		mpz_mod (x1, x1, n);
		mpz_mod (y, y, n);
		if (mpz_probab_prime_p (n, 25))
		{
			//mpz_out_str (stdout, 10, n);
			//fflush (stdout);
			//fputc (' ', stdout);

			mpz_set(p_factors->factors[p_factors->num_elements], n);
			p_factors->num_elements++;
			break;
		}
	}

  mpz_clear(P);
  mpz_clear(t2);
  mpz_clear(t1);
  mpz_clear(x1);
  mpz_clear(x);
  mpz_clear(y);
}

//restituisce 0 se l'algoritmo è terminato e l'intero è stato totalmente fattorizzato
//restituisce -1 se il timer è scaduto e l'intero è stato fattorizzato in parte o per niente fattorizzato
int factor(prime_factors* p_factors, mpz_t t_original, unsigned long p)
{
	assert(WAITING_TIME > 0);

	keep_going = 1; //FONDAMENTALE reimpostarlo ad ogni chiamata di factor!

	signal (SIGALRM, catch_alarm);
  	//alarm(WAITING_TIME);

	if(AUTO_TIMING)
	{
		//gmp_printf("t_original: %Zd\n", t_original);

		int num_digits = mpz_sizeinbase(t_original, 10);
		
		if(num_digits < 200)
			TIMER = 1;
		else if(num_digits >= 200 && num_digits < 250)
			TIMER = 5;
		else if(num_digits >= 250 && num_digits < 300)
			TIMER = 10;
		else //if(num_digits >= 300)
			TIMER = 25; 

		printf("TIMER: %d\n", TIMER);

		alarm(TIMER);
	}
	else
	{
		//printf("TIMER: %d\n", WAITING_TIME);
		alarm(WAITING_TIME);
	}


	mpz_t t;
	mpz_init_set(t, t_original); //faccio una copia della variabile originale

  	unsigned int division_limit;

	if(mpz_sgn (t) == 0)
    	{
		mpz_set_ui(p_factors->factors[p_factors->num_elements], 0);
		p_factors->num_elements++;
		return 0;
	}

  	/* Set the trial division limit according the size of t.  */
 	division_limit = mpz_sizeinbase (t, 2);
  	if(division_limit > 1000)
    		division_limit = 1000 * 1000;
  	else
    		division_limit = division_limit * division_limit;

  	if(p != 0)
    		factor_using_division_2kp (p_factors, t, division_limit / 10, p);
  	else
    		factor_using_division (p_factors, t, division_limit);

	if(!keep_going)
	{
		mpz_clear(t);
		return -1;
	}

  	if (mpz_cmp_ui (t, 1) != 0)
    	{
		if(FLAG_VERBOSE > 0)
		{
			printf ("[is number prime?] ");
			fflush (stdout);
		}
		if(mpz_probab_prime_p (t, 25))
		{	
			//mpz_out_str (stdout, 10, t);
			mpz_set(p_factors->factors[p_factors->num_elements], t);
			p_factors->num_elements++;
		}
		else
			factor_using_pollard_rho (p_factors, t, 1L, p);
    	}

	mpz_clear(t);

	if(!keep_going)		
		return -1;
  	
	return 0;
}

