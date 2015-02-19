#include "atkin.h"


void init_curve_orders(mpz_t p, mpz_t D, resultset* r, curve_orders* c_orders)
{
	mpz_t p_plus_1;
	mpz_init(p_plus_1);
	mpz_add_ui(p_plus_1, p, 1);

	if(mpz_cmp_si(D, -4) == 0)
	{
		mpz_t u;
		mpz_init_set(u, r->results[0]);
		mpz_t v;
		mpz_init_set(v, r->results[1]);

		mpz_t v_2;
		mpz_init(v_2);	
		mpz_mul_ui(v_2, v, 2);

		c_orders->num_elements = 4;
		c_orders->orders = (mpz_t*)malloc(sizeof(mpz_t)*c_orders->num_elements);		

		mpz_init(c_orders->orders[0]);
		mpz_add(c_orders->orders[0], p_plus_1, u);
		
		mpz_init(c_orders->orders[1]);
		mpz_sub(c_orders->orders[1], p_plus_1, u);

		mpz_init(c_orders->orders[2]);		
		mpz_add(c_orders->orders[2], p_plus_1, v_2);

		mpz_init(c_orders->orders[3]);
		mpz_sub(c_orders->orders[3], p_plus_1, v_2);

		mpz_clear(u);
		mpz_clear(v);
		mpz_clear(v_2);
	}
	else if(mpz_cmp_si(D, -3) == 0)
	{
		mpz_t u;
		mpz_init_set(u, r->results[0]);
		mpz_t v;
		mpz_init_set(v, r->results[1]);

		mpz_t v_3;
		mpz_init(v_3);	
		mpz_mul_ui(v_3, v, 3);

		mpz_t u_plus_3v_div2; //(u + 3v)/2
		mpz_init(u_plus_3v_div2);
		mpz_add(u_plus_3v_div2, u, v_3);
		mpz_tdiv_q_ui(u_plus_3v_div2, u_plus_3v_div2, 2);			

		mpz_t u_minus_3v_div2; //(u - 3v)/2
		mpz_init(u_minus_3v_div2);
		mpz_sub(u_minus_3v_div2, u, v_3);
		mpz_tdiv_q_ui(u_minus_3v_div2, u_minus_3v_div2, 2);

		c_orders->num_elements = 6;
		c_orders->orders = (mpz_t*)malloc(sizeof(mpz_t)*c_orders->num_elements);

		mpz_init(c_orders->orders[0]);
		mpz_add(c_orders->orders[0], p_plus_1, u);
		
		mpz_init(c_orders->orders[1]);
		mpz_sub(c_orders->orders[1], p_plus_1, u);

		mpz_init(c_orders->orders[2]);		
		mpz_add(c_orders->orders[2], p_plus_1, u_plus_3v_div2);

		mpz_init(c_orders->orders[3]);
		mpz_sub(c_orders->orders[3], p_plus_1, u_plus_3v_div2);

		mpz_init(c_orders->orders[4]);		
		mpz_add(c_orders->orders[4], p_plus_1, u_minus_3v_div2);		

		mpz_init(c_orders->orders[5]);		
		mpz_sub(c_orders->orders[5], p_plus_1, u_minus_3v_div2);	

		mpz_clear(u);
		mpz_clear(v);
		mpz_clear(v_3);
		mpz_clear(u_plus_3v_div2);
		mpz_clear(u_minus_3v_div2);
	}
	else if(mpz_cmp_si(D, -4) < 0)
	{
		mpz_t u;
		mpz_init_set(u, r->results[0]);
		mpz_t v;
		mpz_init_set(v, r->results[1]);

		c_orders->num_elements = 2;
		c_orders->orders = (mpz_t*)malloc(sizeof(mpz_t)*c_orders->num_elements);

		mpz_init(c_orders->orders[0]);
		mpz_add(c_orders->orders[0], p_plus_1, u);		

		mpz_init(c_orders->orders[1]);
		mpz_sub(c_orders->orders[1], p_plus_1, u);	
	}

	mpz_clear(p_plus_1);
}

void step_one(mpz_t p, mpz_t D, curve_orders* c_orders)
{	
	printf("STEP ONE\n\n");
	assert(mpz_cmp_ui(D, 0) < 0);

	mpz_t Dmod4;
	mpz_init(Dmod4);

	while(1)
	{
		gmp_printf("D: %Zd\n", D);

		mpz_mod_ui(Dmod4, D, 4);

		if(mpz_cmp_ui(Dmod4, 0) != 0 && mpz_cmp_ui(Dmod4, 1) != 0) //D congruo 0,1 mod 4
		{
			mpz_sub_ui(D, D, 1);
			continue;
		}
		
		if(mpz_jacobi(D, p) != 1)
		{
			mpz_sub_ui(D, D, 1);
			continue;
		}		

		resultset r;
		//printf("start\n");
		int ret = cornacchia_GMP(&r, p, D);
		//printf("ret: %d\n", ret);
		//printf("end\n");
		if(ret == 0)
		{
			mpz_sub_ui(D, D, 1);
			continue;
		}
		
		//printf("Cornacchia OK\n");

		init_curve_orders(p, D, &r, c_orders);

		//printf("-------------------------\n");
		int i;
		for(i = 0; i < r.num_elements; i++)
			gmp_printf("cornacchia[%d] %Zd\n", i, r.results[i]);
		printf("\n");
		//printf("-------------------------\n");		

		/*for(i = 0; i < c_orders->num_elements; i++)
			gmp_printf("orders[%d] %Zd\n", i, c_orders->orders[i]);
		printf("-------------------------\n");*/
		
		break;	
	}

	mpz_clear(Dmod4);
}

//Restituisce -1 se non è stato possibile fattorizzare m
//Restituisce  1 se q è probabilmente primo
//Restituisce  2 se q è sicuramente primo, e quindi anche p
int step_two(mpz_t m, mpz_t q, mpz_t p, mpz_t D, curve_orders* c_orders)
{
	printf("STEP TWO\n\n");
	int ret;
	int i;
	int j;
	int set_size = 100; //presumo che sia difficile che un numero abbia più di 100 fattori

	prime_factors p_factors;	
	p_factors.factors = (mpz_t*)malloc(sizeof(mpz_t)*set_size);	
	for(i = 0; i < set_size; i++)
		mpz_init(p_factors.factors[i]);

	mpz_t k;
	mpz_init(k);

	mpz_t lim_inf;
	mpz_init(lim_inf);
	mpz_root(lim_inf, p, 4); //p^(1/4)
	mpz_add_ui(lim_inf, lim_inf, 1);
	mpz_pow_ui(lim_inf, lim_inf, 2); //lim_inf = (n^1/4 + 1)^2

	//gmp_printf("lim_inf: %Zd\n", lim_inf);

	for(i = 0; i < c_orders->num_elements; i++) //ciclo fra tutti i possibili ordini della curva
	{
		p_factors.num_elements = 0; //devo riazzerarlo perchè uso lo stesso array di fattori per ogni ordine della curva e quindi voglio sovrascrivere i vecchi valori

		mpz_set(m, c_orders->orders[i]);		

		gmp_printf("Fattorizzazione ordine curva: %Zd\n", m);


		//dopo aver trovato gli n possibili ordini della curva, tenta per ogni ordine di effettuare una fattorizzazione. Ogni tentativo dura WAITING_TIME secondi, finiti i quali si prova con l'ordine successivo o si torna allo step 1
		int ret_factor = factor(&p_factors, m, 0); //fattorizzi m
		
		/*if(ret_factor == -1)
			printf("FATTORIZZAZIONE INCOMPLETA\n");
		else
			printf("FATTORIZZAZIONE COMPLETA\n");*/


		if(p_factors.num_elements == 0) //non sono riuscito a fattorizzare m
			continue; //provo con il prossimo ordine
		
		//arrivato a questo punto la fattorizzazione potrebbe non esser stata eseguita totalmente, ma comunque alcuni fattori sono stati trovati (controllo su num_elements != 0 sopra)

		for(j = 0; j < p_factors.num_elements; j++) //anche se non ho trovato tutti i fattori, potrei averne trovato uno che soddisfa le mie necessità
		{
			mpz_set(q, p_factors.factors[j]);

			//gmp_printf("q: %Zd\n", q);

			mpz_tdiv_q(k, m, q); //visto che q è un fattore di m, non calcolo il resto poichè è sicuramente nullo - FORSE NON SERVIREBBE PERCHE' CREDO FACTOR NON CONSIDERI IL FATTORE BANALE m, ma lasciamo per sicurezza
			
			if(mpz_cmp_ui(k, 1) == 0) //si verifica solo se q = m ... sarebbe forse meglio controllare solo se p_factors.num_elements = 1 prima del ciclo for interno
				continue; //prova prossimo fattore

			if(mpz_cmp(q, lim_inf) <= 0) //deve essere q > (n^1/4 + 1)^2 (strettamente maggiore)
				continue;

			ret = mpz_probab_prime_p(q, NUM_PROBABILITY_TEST);

			//gmp_printf("ret probab_prime(%Zd): %d\n", q, ret);

			if(ret == 0) //q è composto			
			{	
				continue;
			}
			else if(ret == 1) //q è probabilmente primo	
			{				
				gmp_printf("q: %Zd\n", q);
				gmp_printf("k: %Zd\n\n", k);

				for(i = 0; i < set_size; i++)
					mpz_clear(p_factors.factors[i]);
				mpz_clear(k);
				mpz_clear(lim_inf);
				return 1;
			}
			else if(ret == 2) //q è con certezza (!) primo
			{
				gmp_printf("q: %Zd\n", q);
				gmp_printf("k: %Zd\n\n", k);

				for(i = 0; i < set_size; i++)
					mpz_clear(p_factors.factors[i]);
				mpz_clear(k);
				mpz_clear(lim_inf);		
				return 2;
			}
		}
		
		if(ret_factor == -1) //se arrivo qui, almeno un fattore è stato trovato, ma non tutti quanti
		{
			//printf("PROVO CON L'ALTRO FATTORE\n");		

			mpz_set_ui(k, 1); //necessario sotto per la moltiplicazione dei fattori

			for(j = 0; j < p_factors.num_elements; j++)
			{
				//gmp_printf("fattore[%d]: %Zd\n", j, p_factors.factors[j]);
				mpz_mul(k, k, p_factors.factors[j]);
			}			

			if(mpz_cmp_ui(k, 1) == 0) //credo non si verifichi mai - q = m
				continue; //prova prossimo ordine

			if(mpz_cmp(k, m) == 0) //non si dovrebbe verificare mai, signfica che è scaduto il timer ma comunque ho fattorizzato tutto m - controllo per sicurezza - q = 1
				continue;

			mpz_tdiv_q(q, m, k);

			if(mpz_cmp(q, lim_inf) <= 0) //deve essere q > (n^1/4 + 1)^2 (strettamente maggiore)
				continue; //prova prossimo ordine

			ret = mpz_probab_prime_p(q, NUM_PROBABILITY_TEST);

			//gmp_printf("ret probab_prime(%Zd): %d\n", q, ret);

			

			if(ret == 0) //q è composto			
			{	
				continue;
			}
			else if(ret == 1) //q è probabilmente primo	
			{				
				gmp_printf("q: %Zd\n", q);
				gmp_printf("k: %Zd\n\n", k);

				for(i = 0; i < set_size; i++)
					mpz_clear(p_factors.factors[i]);
				mpz_clear(k);
				mpz_clear(lim_inf);
				return 1;
			}
			else if(ret == 2) //q è con certezza (!) primo
			{
				gmp_printf("q: %Zd\n", q);
				gmp_printf("k: %Zd\n\n", k);

				for(i = 0; i < set_size; i++)
					mpz_clear(p_factors.factors[i]);
				mpz_clear(k);
				mpz_clear(lim_inf);		
				return 2;
			}			
		}

	}

	for(i = 0; i < set_size; i++)
		mpz_clear(p_factors.factors[i]);
	mpz_clear(k);
	mpz_clear(lim_inf);

	//non è stato trovato l'ordine
	return -1;
}

/*
Algorithm 7.5.9 (CM method for generating curves and orders). We assume a list
 of fundamental discriminants {Dj < 0 : j = 1, 2, 3, . . .} ordered, say, 
 by increasing class number h(D), and within the same class number by increasing |D|. 
 We are given a prime p > 3. The algorithm reports (optionally)
possible curve orders or (also optionally) curve parameters for CM curves associ-
ated with the various Dj.*/

//Restituisce -1 se bisogna prendere il discriminante successivo
//Restituisce  0 se tutto è andato ok
int step_three(mpz_t p, mpz_t D, curve_parameters* paramSet){
	printf("STEP THREE\n\n");
	//gmp_printf("D: %Zd\n", D);
	//1. [Calculate nonresidue]
	mpz_t g;
	mpz_init(g);
	gmp_randstate_t state;
	gmp_randinit_default(state); //in utility.c ho usato gmp_randinit_mt, dovrebbe essere la stessa cosa (c'è scritto sono la stessa funzione)
	
	//mpz_t seed;
	//mpz_init_set_ui(seed, 5467799);
	//gmp_randseed(state, seed);	

	while(1){
		
		mpz_urandomm(g, state, p); //genera valore in [0, p-1]
		//gmp_printf("g: %Zd\n", g);
		
		//Find a random quadratic nonresidue g (mod p)
		if(mpz_jacobi(g, p) != -1)  //se (g/p) = -1 allora g non è un residuo quadratico
			continue; //ho trovato un residuo quadratico (a me serve NON quadratico) quindi estraggo il prossimo numero
		
		//dopo aver trovato un residuo non quadratico verifico l'altra condizione

		mpz_t res;
		mpz_init(res);
		mpz_mod_ui(res, p, 3);	//p mod 3

		// In case D = −3 is used, g must also be a noncube modulo - guarda wikipedia Cubic_reciprocity#Cubic_residue_character
		if(mpz_cmp_ui(res, 1) == 0) //p mod 3 = 1
		{
			mpz_t exp;
			mpz_init(exp);	
			mpz_sub_ui(exp, p, 1);			
			mpz_tdiv_q_ui(exp, exp, 3);			
			mpz_powm(exp, g, exp, p);	//g^((p-1)/3) mod p
			
			if(mpz_cmp_ui(exp, 1) == 0)
			{
				mpz_clear(res);
				mpz_clear(exp);	
				continue;
			}

			mpz_clear(exp);	
		}
	
		mpz_clear(res);
		break;
	}

	//2. [Discriminant loop]
	//in realtà questa istruzione nel caso nostro si può spostare prima dell'istruzione 1 per evitare calcoli inutili
	if(mpz_jacobi(D, p) != 1)
		return -1;
	
	//3. [Seek a quadratic form for 4p]
	resultset r_set; //Viene inizializzata dentro a cornacchia_GMP		
	if(!cornacchia_GMP(&r_set, p, D)) //return 1 se OK
		return -1;			
					
	//4. [Option: Curve orders]
	curve_orders* c_orders = (curve_orders*) malloc(sizeof(curve_orders));
	init_curve_orders(p, D, &r_set, c_orders);

	//5. [Option: Curve parameters]	
	if(mpz_cmp_si(D, -4) == 0){  //Four curves y^2 = x^3 − (g^k)*x    k = 0, 1, 2, 3
		int k;
		
		paramSet->num_elements = 4;
		paramSet->a = (mpz_t*) malloc(sizeof(mpz_t)*paramSet->num_elements);
		paramSet->b = (mpz_t*) malloc(sizeof(mpz_t)*paramSet->num_elements);
		
		mpz_t ginv;
		mpz_init(ginv);

		for(k=0; k<paramSet->num_elements; ++k){			
			
			mpz_set(ginv, g);
			mpz_pow_ui(ginv, ginv, k);
			mpz_mul_si(ginv, ginv, -1);			
			mpz_mod(ginv, ginv, p);
			
			mpz_init_set(paramSet->a[k], ginv);
			mpz_init_set_ui(paramSet->b[k], 0);					
		}
		
		printf("Possibili parametri curva:\n");
		for(k=0; k<paramSet->num_elements; ++k)
		{
			gmp_printf("Parametro A[%d]: %Zd\t", k, paramSet->a[k]);	
			gmp_printf("Parametro B[%d]: %Zd\n", k, paramSet->b[k]);
		}
		printf("\n");

		mpz_clear(g);
		gmp_randclear(state);
		mpz_clear(ginv);

		//return {(a, b)} = {(−g^k mod p, 0) : k = 0, 1, 2, 3};
		return 0;
	}
	if(mpz_cmp_si(D, -3) == 0){//Six curves y^2 = x^3 − (g^k)    k = 0, 1, 2, 3, 4, 5
		int k;
		paramSet->num_elements = 6; //quanto il numero di ordini possibili
		paramSet->a = (mpz_t*) malloc(sizeof(mpz_t)*paramSet->num_elements);
		paramSet->b = (mpz_t*) malloc(sizeof(mpz_t)*paramSet->num_elements);
		
		mpz_t ginv;
		mpz_init(ginv);		

		for(k=0; k<paramSet->num_elements; ++k){			
			mpz_set(ginv, g);
			mpz_pow_ui(ginv, ginv, k);
			mpz_mul_si(ginv, ginv, -1);			
			mpz_mod(ginv, ginv, p);
			
			mpz_init_set_ui(paramSet->a[k], 0);
			mpz_init_set(paramSet->b[k], ginv);					
		}

		printf("Possibili parametri curva:\n");
		for(k=0; k<paramSet->num_elements; ++k)
		{
			gmp_printf("Parametro A[%d]: %Zd\t", k, paramSet->a[k]);	
			gmp_printf("Parametro B[%d]: %Zd\n", k, paramSet->b[k]);
		}
		printf("\n");

		mpz_clear(g);
		gmp_randclear(state);
		mpz_clear(ginv);

		//return {(a, b)} = {(0, −g^k mod p) : k = 0, 1, 2, 3, 4, 5};
		return 0;
	} 	

	//6. [Continuation for D < −4]
	//Compute the Hilbert class polynomial T ∈ Z[X], via Algorithm 7.5.8;	
	Hresult* resHilbert = (Hresult*) malloc(sizeof(Hresult));
	hilbertclass(resHilbert, D, p); //Esegue anche la riduzione modulo p
	
	//printf("----------------------\n");
	printf("\nPolinomio di Hilbert ridotto modulo p:\n");
	r_big_polynomial_printf(&(resHilbert->Hpoly));
	printf("\n");
	//printf("----------------------\n");	

	//S = T mod p;	// Reduce to polynomial ∈ Fp [X]. Fatto dentro Hilbert
	
	//Obtain a root j ∈ Fp of S, via Algorithm 2.3.10;

	mpz_t one_root; //una radice del polinomio - inizializzata o con il primo metodo basato sul termine noto o con l'algoritmo pesante di fattorizzazione completa
	mpz_init(one_root);

	int success = 0; //sono riuscito a trovare una radice a partire dal termine noto (fallisce solo se termine noto = 0 - credo mai a parte D = -3 - e se non riesce a fattorizzare il termine noto)

	//provo a cercare una radice fra i fattori del termine noto, funziona sempre a meno che non fallisca nella fattorizzaione del termine noto //---////////////////////---//
	if((mpz_cmp_ui(resHilbert->Hpoly.degrees[0], 0) == 0 && mpz_cmp_ui(resHilbert->Hpoly.coefficients[resHilbert->Hpoly.num_elements - 1], 1) == 0) /*|| mpz_cmp_ui(resHilbert->Hpoly.degrees[resHilbert->Hpoly.num_elements - 1], 2) > 0*/) //se è presente il termine noto (CREDO SEMPRE A PARTE IL CASO D = -3 che però è trattato sopra) e il polinomio è monico
	{
		//se il coefficiente del termine di grado massimo fosse diverso da 1, allora bisognerebbe cercare le radici fra i divisori di p/q  (p termine noto, q coefficiente termine grado massimo) - non lo faccio perchè tanto è sempre monico il polinomio
		mpz_t low_coef; //termine noto
		mpz_init_set(low_coef, resHilbert->Hpoly.coefficients[0]);

		gmp_printf("Teorema Radici Razionali: cerco radice fra i divisori del termine noto: %Zd\n", low_coef);

		prime_factors prime;
		int prime_factors_size = 150; //per essere sicuri
		prime.factors = (mpz_t*)malloc(sizeof(mpz_t)*prime_factors_size);
		prime.num_elements = 0;		

		int i;
		int j;

		for(i = 0; i < prime_factors_size; i++)
			mpz_init(prime.factors[i]);

		//printf("Fattorizzo termine noto\n");
		factor(&prime, low_coef, 0);

		//for(i = 0; i < prime.num_elements; i++)
		//	gmp_printf("factors[%d]: %Zd\n", i, prime.factors[i]);

		mpz_t* dist_factors = (mpz_t*)malloc(sizeof(mpz_t)*prime.num_elements); //in realtà la dimensione finale è probabilmente molto minore - se tanti fattori uguali - poi rialloco

		int dist_factors_size = 0; //dimensione reale

		for(i = 0; i < prime.num_elements; i++) //scorro su tutti i fattori e prendo ogni fattore solo una volta
		{
			int found = 0;
			for(j = 0; j < dist_factors_size; j++)
			{
				if(mpz_cmp(prime.factors[i], dist_factors[j]) == 0) //ho già inserito questo fattore
				{
					found = 1;
					break;
				}
			}
			if(!found)
			{
				mpz_init_set(dist_factors[dist_factors_size], prime.factors[i]);
				dist_factors_size++;
			}
		}

		dist_factors = (mpz_t*)realloc(dist_factors, sizeof(mpz_t)*dist_factors_size);

		//for(i = 0; i < dist_factors_size; i++)
		//	gmp_printf("dist_factors[%d]: %Zd\n", i, dist_factors[i]);

		int* multiplicity = (int*)malloc(sizeof(int)*dist_factors_size);
		
		for(i = 0; i < dist_factors_size; i++)
		{
			int count = 0;
			for(j = 0; j < prime.num_elements; j++)
			{
				if(mpz_cmp(dist_factors[i], prime.factors[j]) == 0)
					count++;
			}
			multiplicity[i] = count;
		}

		//for(i = 0; i < dist_factors_size; i++)
		//	gmp_printf("multiplicity[%d]: %d\n", i, multiplicity[i]);

			
	

		mpz_t currentResult;
		mpz_init_set_ui(currentResult, 1);

		//time_t start, stop;
		//clock_t ticks;
		//time(&start);

		success = findOneRootInDivisors(one_root, resHilbert->Hpoly, p, dist_factors, multiplicity, dist_factors_size, 0, currentResult);
		//printf("success after: %d\n", success);
		//findDivisors(dist_factors, multiplicity, dist_factors_size, 0, currentResult);

		//ticks = clock();
		//time(&stop);
		//printf("Used %0.2f seconds of CPU time. \n", (double)ticks/CLOCKS_PER_SEC);
		//printf("Finished in about %.0f seconds. \n", difftime(stop, start));

		//clear
		for(i = 0; i < prime_factors_size; i++)
			mpz_clear(prime.factors[i]);		

		for(i = 0; i < dist_factors_size; i++) //ho riallocato a questa dimensione
			mpz_clear(dist_factors[i]);

		free(prime.factors);
		free(dist_factors);
		free(multiplicity);

		mpz_clear(low_coef);
	}  //---////////////////////---//
	
	//printf("success: %d\n", success);
	if(success)
	{
		gmp_printf("Radice polinomio: %Zd\n", one_root);
		printf("\n");
	}
	else
		printf("Teorema Radici Razionali fallito\n");	

	if(!success || mpz_cmp_ui(resHilbert->Hpoly.degrees[0], 0) != 0 || mpz_cmp_ui(resHilbert->Hpoly.coefficients[resHilbert->Hpoly.num_elements - 1], 1) != 0 /*|| mpz_cmp_ui(resHilbert->Hpoly.degrees[resHilbert->Hpoly.num_elements - 1], 2) <= 0*/) //se ho fallito prima nella ricerca della radice o se il termine noto è uguale a zero [in realtà nel secondo caso NON DOVREI USARE LA RADICE ZERO!! (non capita mai)] o se il polinomio non è monico [non capita mai]  ////////////////////////
	{
		/*if(mpz_cmp_ui(resHilbert->Hpoly.degrees[resHilbert->Hpoly.num_elements - 1], 2) > 0)			
		{
			printf("Ricerca radici polinomio di grado maggiore di 2\n");
			printf("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
			printf("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
			printf("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
			printf("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");

			//printf("ATTENZIONE: POTREBBE OCCUPARE TUTTA LA RAM - PREMERE INVIO PER PROCEDERE\n");
			//getchar();
		}*/


		/*if(mpz_cmp_ui(resHilbert->Hpoly.degrees[resHilbert->Hpoly.num_elements - 1], 2) > 0) //se il grado del polinomio è maggiore di 2, provo a cercare la radice testando tutti gli elementi del campo...lento ma meglio dell'algoritmo successivo
		{
			printf("Ricerca radice fra tutti gli elementi di Zp\n");
			find_one_root(one_root, resHilbert->Hpoly, p);

			if(mpz_cmp_si(one_root, -1) == 0)
			{
				printf("RADICE INESISTENTE - QUIT\n");
				exit(0);
			}
		}*/
		//else //grado minore uguale a 2
		{
			s_root_set* root_set = (s_root_set*) malloc(sizeof(s_root_set));
			roots(root_set, resHilbert->Hpoly, p);

			mpz_set(one_root, root_set->roots[0]);

			gmp_printf("Radice polinomio di Hilbert: %Zd\n", one_root);
			printf("\n");

			int i;
			for(i = 0; i < root_set->num_elements; i++) //dentro alla funzione roots, viene settata la dimensione di root_set_root a num_elements
				mpz_clear(root_set->roots[i]);

			free(root_set->roots);
			free(root_set);
		}
	}
	
	mpz_t c;
	mpz_init(c);
	mpz_t r;
	mpz_init(r);
	mpz_t s;
	mpz_init(s);
	
	//c = j*((j − 1728)^(−1)) mod p;	j è una radice
	mpz_t jtemp;
	mpz_init(jtemp);
	mpz_sub_ui(jtemp, one_root, 1728); //prendiamo la prima radice
	mpz_invert(jtemp, jtemp, p); //E SE INVERSIONE NON ANDATA A BUON FINE?
	mpz_mul(jtemp, one_root, jtemp); 
	mpz_mod(c, jtemp, p);
	//r = −3*c mod p;
	mpz_mul_si(r, c, -3);
	mpz_mod(r, r, p);
	//s = 2*c mod p;
	mpz_mul_ui(s, c, 2);
	mpz_mod(s, s, p);
	
	//7. [Return two curve-parameter pairs]
	paramSet->num_elements = 2;
	paramSet->a = (mpz_t*) malloc(sizeof(mpz_t)*paramSet->num_elements);
	paramSet->b = (mpz_t*) malloc(sizeof(mpz_t)*paramSet->num_elements);
	
	mpz_t g2;
	mpz_init(g2);
	mpz_t rg2;
	mpz_init(rg2);
	mpz_pow_ui(g2, g, 2);
	mpz_mul(rg2, r, g2);
	mpz_mod(rg2, rg2, p);
	
	mpz_t g3;
	mpz_init(g3);
	mpz_t sg3;
	mpz_init(sg3);
	mpz_pow_ui(g3, g, 3);
	mpz_mul(sg3, s, g3);
	mpz_mod(sg3, sg3, p);
	
	mpz_init(paramSet->a[0]);
	mpz_init(paramSet->b[0]);
	mpz_init(paramSet->a[1]);
	mpz_init(paramSet->b[1]);
	
	mpz_set(paramSet->a[0], r);
	mpz_set(paramSet->b[0], s);	
	
	mpz_set(paramSet->a[1], rg2);
	mpz_set(paramSet->b[1], sg3);

	printf("Possibili parametri curva:\n");
	int k;
	for(k=0; k<paramSet->num_elements; ++k)
	{
		gmp_printf("Parametro A[%d]: %Zd\t", k, paramSet->a[k]);	
		gmp_printf("Parametro B[%d]: %Zd\n", k, paramSet->b[k]);
	}
	printf("\n");

	mpz_clear(g);
	gmp_randclear(state);
	mpz_clear(c);
	mpz_clear(r);
	mpz_clear(s);
	mpz_clear(jtemp);
	mpz_clear(g2);
	mpz_clear(rg2);
	mpz_clear(g3);
	mpz_clear(sg3);

	//Two curves: y^2 = x^3 − 3c(g^(2k))x + 2c(g^(3k))    k = 0, 1
	//return {(a, b)} = {(r, s), (rg^2 mod p, s(g^3) mod p)};
	return 0;
}

//P e out_chosen_parameters sono parametri di output, il secondo indica i parametri di curva "scelti" fra le coppie in c_parameters fornito in input
//p, c_parameters, c_orders e m sono parametri di input
//restituisce -2 se andato in loop infinito (programma termina)
//restituisce -1 se n è composto
//restituisce 0 se tutto ok
int step_four(curve_parameters* out_chosen_parameters, curve_point* P, mpz_t p, curve_parameters c_parameters, curve_orders c_orders, mpz_t m)
{
	printf("STEP FOUR\n\n");

	////////////////////////////////////////////
	//sò che ho una curva di ordine m, ma non sò quali parametri (a, b) identificano tale curva. Cerco un punto sulla curva e verifico che appartanga a questa e solo questa
		
	curve_point points[2];
	init(&points[0]);
	init(&points[1]);

	curve_point P_mul;
	init(&P_mul);

	int i;
	int j;

	while(1) //magari ad una estrazione posso prendere un punto che appartiene su più curve, ma ci deve essere un punto che appartiene ad una e una sola curva, altrimenti le curve sarebbero uguali
	{
		int found = 0;
		for(i = 0; i < c_parameters.num_elements; i++)
		{
			find_two_points(points, c_parameters.a[i], c_parameters.b[i], p); //trovo due punti (x, y) e (x, -y) appartenenti alla curva y^2 = x^3 + Ax + B
		
			/*gmp_printf("MODULO: %Zd\n", p);
			gmp_printf("A: %Zd\n", c_parameters.a[i]);
			gmp_printf("B: %Zd\n", c_parameters.b[i]);
		
			gmp_printf("points[0].x: %Zd\n", points[0].x);
			gmp_printf("points[0].y: %Zd\n", points[0].y);
			gmp_printf("points[0].z: %Zd\n", points[0].z);

			gmp_printf("ORDINE m: %Zd\n", m);*/

			//testo solo il primo punto di points, tanto se (x,y) appartiene a E: y^2 = x^3 + Ax + B allora anche (x, -y) appartiene a E (c'è il termine y^2), quindi inutile continuare a testare per (x, -y)
			//prima testo se il punto (x, y) appartiene solo ad una curva, poi verifico che quella curva sia quella di ordine m
			int count = 0; //conta il numero di curve a cui appartiene il punto che stiamo testando	
			int index = -1;
			for(j = 0; j < c_orders.num_elements; j++)
			{
				int ret = mul(&P_mul, points[0], c_orders.orders[j], p, c_parameters.a[i]);
			
				if(ret == 1) //ho ottenuto il punto all'infinito, quindi il punto appartiene alla curva di ordine c_orders.orders[j]
				{
					//gmp_printf("ORDINE: %Zd\n", c_orders.orders[j]);		
					index = j; //ha importanza solo se count rimane uguale a 1
					count++;
				}
			}
			
			if(count == 1) //una volta che sò che il punto appartiene ad una sola curva, verifico che questa sia quella di ordine m
			{
				if(mpz_cmp(m, c_orders.orders[index]) == 0)
				{	
					//se tutto è andato a buon fine con questi parametri a e b, me li salvo per uso futuro (step 5 - somma di due punti su di una curva ellittica)
					out_chosen_parameters->a = (mpz_t*)malloc(sizeof(mpz_t));
					out_chosen_parameters->b = (mpz_t*)malloc(sizeof(mpz_t));
					out_chosen_parameters->num_elements = 1;

					mpz_init_set(out_chosen_parameters->a[0], c_parameters.a[i]);
					mpz_init_set(out_chosen_parameters->b[0], c_parameters.b[i]);

					printf("Parametri curva selezionata:\n");
					gmp_printf("Parametro A: %Zd\n", out_chosen_parameters->a[0]);
					gmp_printf("Parametro B: %Zd\n", out_chosen_parameters->b[0]);
					printf("\n");

					found = 1;
					break; //interrompo il for esterno
				}
			}
		}
		if(found == 1)
			break;
	}
	
	clear(&points[0]);
	clear(&points[1]);
	clear(&P_mul);
	////////////////////////////////////////////

	init(P); //faccio subito la init in modo tale che le clear lavorino sempre su coordinate di un punto inizializzato

	gmp_randstate_t state;
	gmp_randinit_mt(state);

	mpz_t x;
	mpz_init(x);
	mpz_t x_3;
	mpz_init(x_3);
	mpz_t ax;
	mpz_init(ax);


	mpz_t a;
	mpz_init_set(a, out_chosen_parameters->a[0]);	
	mpz_t b;
	mpz_init_set(b, out_chosen_parameters->b[0]);


	mpz_t Q;
	mpz_init(Q);
	mpz_t y;
	mpz_init(y);
	mpz_t y_2;
	mpz_init(y_2);
	mpz_t temp_mod;
	mpz_init(temp_mod);
	mpz_t p_2;
	mpz_init(p_2);
	mpz_mul_ui(p_2, p, 2);

	mpz_t counter;
	mpz_init_set_ui(counter, 0);

	int loop = 0;

	while(1)
	{
		//credo posso anche eliminare questo controllo, tanto a questo step si arriva per un p moooolto grande e quindi esiste per forza la x con le proprietà volute
		if(mpz_cmp(counter, p_2) == 0) //controllo in cui verifico se siano stati estratti tutti i numeri e in tal caso esco, per evitare cicli infiniti (metto p*2 per sicurezza = ogni numero estratto 2 volte)
		{
			loop = 1;
			break;
		}

		mpz_urandomm(x, state, p); //estrae numero casuale in [0, p-1], corretto mettere p come terzo parametro
			

		mpz_pow_ui(x_3, x, 3);
		mpz_mul(ax, x, a);
		mpz_add(Q, x_3, ax);
		mpz_add(Q, Q, b); //Q = x^3 + ax + b
		mpz_mod(Q, Q, p); //Q = (x^3 + ax + b) mod p

		//gmp_printf("LA SECONDA CONDIZIONE DELL'IF L'HO AGGIUNTA IO!!\n");

		if(mpz_jacobi(Q, p) == -1 /**/|| mpz_jacobi(Q, p) == 0/**/) //la seconda condizione l'ho aggiunta io, così mi assicuro che RADIX_GMP abbia soluzione (per p = 1004730010023465367117639 usciva Q = 0)
		{
			mpz_add_ui(counter, counter, 1);	
			continue;
		}

		//gmp_printf("Q: %Zd\n", Q);

		break;	
	}	

	if(loop) //sarei andato in loop infinito, senza trovare una x valida
	{
		gmp_randclear(state);
		mpz_clear(x);
		mpz_clear(x_3);
		mpz_clear(ax);
		mpz_clear(a);
		mpz_clear(b);
		mpz_clear(Q);
		mpz_clear(y);
		mpz_clear(y_2);
		mpz_clear(temp_mod);
		mpz_clear(p_2);
		mpz_clear(counter);

		return -2;
	}		

	//gmp_printf("p: %Zd\n", p);
	//gmp_printf("Q: %Zd\n", Q);

	int ret = radix_GMP(y, p, Q); //restituisce 0 in caso di errore

	if(ret == 0)
	{
		//printf("STEP FOUR - RADIX_GMP RET = 0\n");
		gmp_randclear(state);
		mpz_clear(x);
		mpz_clear(x_3);
		mpz_clear(ax);
		mpz_clear(a);
		mpz_clear(b);
		mpz_clear(Q);
		mpz_clear(y);
		mpz_clear(y_2);
		mpz_clear(temp_mod);
		mpz_clear(p_2);
		mpz_clear(counter);
		
		return -1; //n composite
	}		

	//gmp_printf("y: %Zd\n", y);
	//gmp_printf("Q: %Zd\n", Q);
	//gmp_printf("p: %Zd\n", p);
	
	mpz_pow_ui(y_2, y, 2); //y^2
	mpz_mod(temp_mod, y_2, p);		

	if(mpz_cmp(temp_mod, Q) != 0)
	{
		//printf("STEP FOUR - MOD != 0\n");
		gmp_randclear(state);
		mpz_clear(x);
		mpz_clear(x_3);
		mpz_clear(ax);
		mpz_clear(a);
		mpz_clear(b);
		mpz_clear(Q);
		mpz_clear(y);
		mpz_clear(y_2);
		mpz_clear(temp_mod);
		mpz_clear(p_2);
		mpz_clear(counter);
		
		return -1; //n composite
	}	

	//x e y per costruzione sono già ridotti modulo p
	mpz_set(P->x, x);
	mpz_set(P->y, y);
	mpz_set_ui(P->z, 1);
	//P->valid = 1;

	gmp_printf("P->x: %Zd\n", P->x);
	gmp_printf("P->y: %Zd\n", P->y);
	//gmp_printf("P->z: %Zd\n", P->z);
	printf("\n");

	gmp_randclear(state);
	mpz_clear(x);
	mpz_clear(x_3);
	mpz_clear(ax);
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(Q);
	mpz_clear(y);
	mpz_clear(y_2);
	mpz_clear(temp_mod);
	mpz_clear(p_2);
	mpz_clear(counter);

	return 0;
}

//restituisce 0 se n è primo
//restituisce 1 se bisogna andare allo step 4
//restituisce -1 se n è composto
//restituisce 2 se bisogna rieseguire Atkin-Morain
int step_five(mpz_t m, mpz_t q, curve_point P, curve_parameters out_chosen_parameters, mpz_t p_mod)
{
	printf("STEP FIVE\n\n");
	int ret;	

	mpz_t k;
	mpz_init(k);
	mpz_tdiv_q(k, m, q); //q divide m e k != 1 per costruzione

	//gmp_printf("P.x: %Zd\n", P.x);
	//gmp_printf("P.y: %Zd\n", P.y);
	//gmp_printf("P.z: %Zd\n", P.z);

	//gmp_printf("m/q = k: %Zd\n", k);

	//gmp_printf("A: %Zd\n", out_chosen_parameters.a[0]);

	curve_point U;
	init(&U);

	printf("Calcolo U = [m/q]P\n");

	ret = mul(&U, P, k, p_mod, out_chosen_parameters.a[0]); //U = [m/q]P

	//gmp_printf("U.x: %Zd\n", U.x);
	//gmp_printf("U.y: %Zd\n", U.y);
	//gmp_printf("U.z: %Zd\n", U.z);

	//gmp_printf("q: %Zd\n", q);	

	if(ret == -1) //ret = -1 errore inversione
	{
		printf("Errore inversione\n");
		mpz_clear(k);
		clear(&U);

		return -1; //n è composto
	}
	else if(ret == 1) //ret = 1 punto all'infinito
	{
		printf("U punto all'infinito\n");
		mpz_clear(k);
		clear(&U);

		return 1; //devo andare a step 4
	}
	//punto non all'infinito
	gmp_printf("U.x: %Zd\n", U.x);
	gmp_printf("U.y: %Zd\n", U.y);
	printf("\n");

	curve_point V;
	init(&V);

	printf("Calcolo V = [q]U\n");

	ret = mul(&V, U, q, p_mod, out_chosen_parameters.a[0]); //V = [q]U

	//gmp_printf("V.x: %Zd\n", V.x);
	//gmp_printf("V.y: %Zd\n", V.y);
	//gmp_printf("V.z: %Zd\n", V.z);


	//printf("RET MUL FIVE: %d\n", ret);
	if(ret == -1) //ret = -1 errore inversione
	{
		printf("Errore inversione\n");
		mpz_clear(k);
		clear(&U);
		clear(&V);

		return -1; //n è composto
	}
	else if(ret == 0) //V != O
	{
		printf("V non è punto all'infinito\n");
		gmp_printf("V.x: %Zd\n", V.x);
		gmp_printf("V.y: %Zd\n", V.y);

		mpz_clear(k);
		clear(&U);
		clear(&V);

		return -1; //n è composto
	}
	printf("V punto all'infinito\n");
	printf("\n");

	printf("IF Q IS PRIME, THEN N IS PRIME - Ricorsione\n\n");

	gmp_printf("Q: %Zd\n", q);
	gmp_printf("N: %Zd\n", p_mod);
	printf("\n");

	ret = mpz_probab_prime_p(p_mod, NUM_PROBABILITY_TEST);

	if(ret == 0) //n è composto
	{
		mpz_clear(k);
		clear(&U);
		clear(&V);

		return -1; //n è composto
	}
	else if(ret == 2) //n è sicuramente primo
	{
		mpz_clear(k);
		clear(&U);
		clear(&V);

		return 0; //n è sicuramente primo, posso terminare il programma
	}
	else if(ret == 1) //n è probabilmente primo
	{
		mpz_clear(k);
		clear(&U);
		clear(&V);

		return 2; //n è probabilmente primo, lancio un'altra volta Atkin-Morain
	}

	return -3; //istruzione mai raggiunta
}

//restituisce 0 se p è un numero primo
//restituisce -1 se p non è un numero primo
int Atkin_Morain(mpz_t p)
{
	gmp_printf("P=%Zd E' UN NUMERO PRIMO?\n", p);
	//manca controllo gcd(p, 6) = 1...serve!?	
	int ret;

	ret = mpz_probab_prime_p(p, NUM_PROBABILITY_TEST);

	if(ret == 0) //n è composto
	{
		gmp_printf("N = %Zd NON è un numero primo\n", p);
		return -1;
	}
	else if(ret == 2) //n è sicuramente primo
	{
		gmp_printf("N = %Zd è un numero primo\n", p);
		return 0;
	}

	//printf("Superato primo controllo\n");

	//se ret == 1, n è probabilmente primo, quindi andiamo avanti con il test

	mpz_t D;
	mpz_init_set_str(D, "-3", 10);
	mpz_t Dmod4;
	mpz_init(Dmod4);

	mpz_t m;
	mpz_init(m);
	mpz_t q;
	mpz_init(q);	

	curve_parameters c_parameters; //sono presenti più coppie a e b
	curve_orders c_orders; //più ordini, nello step 4 si associa all'ordine m la relativa coppia (a, b)

	//alla fine del ciclo vengono generate delle coppie (a, b)
	while(1)
	{
		gmp_printf("D: %Zd\n", D);
		mpz_mod_ui(Dmod4, D, 4);

		if(mpz_cmp_ui(Dmod4, 0) != 0 && mpz_cmp_ui(Dmod4, 1) != 0) //D congruo 0,1 mod 4
		{
			mpz_sub_ui(D, D, 1);
			continue;
		}		

		//1. [Choose discriminant]		
		step_one(p, D, &c_orders);
		//gmp_printf("D: %Zd\n", D);
		//2. [Factor orders]
		ret = step_two(m, q, p, D, &c_orders); //restituisce 0 se non è riuscito a trovare l'ordine m, 1 se ha avuto successo, 2 se q e quindi n sono con certezza primi
			
		//printf("RET STEP TWO: %d\n", ret);
		
		if(ret == -1) //torna al passo 1. [Choose discriminant]
		{
			mpz_sub_ui(D, D, 1);
			continue;
		}
		else if(ret == 2)
		{
			//gmp_printf("D: %Zd\n", D);
			gmp_printf("N = %Zd è un numero primo\n", p); //se q è primo (con certezza, non probabilisticamente!), allora anche n lo è (con certezza)

			mpz_clear(D);
			mpz_clear(Dmod4);
			mpz_clear(m);
			mpz_clear(q);

			//c_parameters NON è da "cancellare" in quanto viene inizializzato nello step three e quindi qui ancora non è inizializzato
			int i;
			for(i = 0; i < c_orders.num_elements; i++)
				mpz_clear(c_orders.orders[i]);
			

			return 0;
		}		
		//else if(ret == 1) //trovato un q probabilmente primo
		
		//3. [Obtain curve parameters]
		ret = step_three(p, D, &c_parameters);
		
		//printf("RET STEP THREE: %d\n", ret);
		
		if(ret == -1) //errore su jacobi o cornacchia
		{
			mpz_sub_ui(D, D, 1);
			continue;
		}
		else
			break; //ho trovato i parametri della curva
	
	}

	

	curve_point P;
	curve_parameters out_chosen_parameters;

	int ret_Atkin;

	while(1)
	{	
		//4. [Choose point on Ea,b(Zn)]
		ret = step_four(&out_chosen_parameters, &P, p, c_parameters, c_orders, m);

		if(ret == -1)
		{
			gmp_printf("N = %Zd NON è un numero primo\n", p);
			
			ret_Atkin = -1;
			break;
		}
		else if(ret == -2)
		{
			gmp_printf("LOOP INFINITO\n", p);
			exit(0);
		}

		

		//5. [Operate on point]
		ret = step_five(m, q, P, out_chosen_parameters, p); //gestione errori
		//restituisce 0 se n è primo, 1 se bisogna andare allo step 4, -1 se n è composto, 2 se bisogna rieseguire Atkin-Morain

		if(ret == -1)
		{
			gmp_printf("N = %Zd NON è un numero primo\n", p);
			ret_Atkin = -1;
			break;	
		}
		else if(ret == 0)
		{
			gmp_printf("N = %Zd è un numero primo\n", p);
			ret_Atkin = 0;
			break;
		}
		else if(ret == 1)
		{
			continue; //procedi allo step 4
		}
		else if(ret == 2) //rieseguo Atkin Morain (che potrà richiamare Atkin altre volte) e poi interrompo l'Atkin superiore
		{
			ret_Atkin = Atkin_Morain(q);

			if(ret_Atkin == 0)	
				gmp_printf("N = %Zd è un numero primo\n", p);
			else if(ret_Atkin == -1)
				gmp_printf("N = %Zd NON è un numero primo\n", p);

			break;
		}
	}
	
	//Clear
	int i;
	for(i = 0; i < c_parameters.num_elements; i++)
	{
		mpz_clear(c_parameters.a[i]);
		mpz_clear(c_parameters.b[i]);
	}
	
	for(i = 0; i < c_orders.num_elements; i++)
		mpz_clear(c_orders.orders[i]);

	for(i = 0; i < out_chosen_parameters.num_elements; i++)
	{
		mpz_clear(out_chosen_parameters.a[i]);
		mpz_clear(out_chosen_parameters.b[i]);
	}
	
	mpz_clear(P.x);
	mpz_clear(P.y);
	
	mpz_clear(D);
	mpz_clear(Dmod4);
	mpz_clear(m);
	mpz_clear(q);

	return ret_Atkin;
}

