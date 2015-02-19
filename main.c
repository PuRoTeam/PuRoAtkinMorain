#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "tools/polynomial.h"
#include "tools/complex_polynomial.h"
#include "tools/atkin.h"
#include "tools/hilbert.h"
#include "tools/factorize.h"
//#include <glib.h>

int test_cornacchia();
int test_gcd();
void test_roots();
void test_pow();
void test_fast_mod();
void test_hilbert();
void test_fast_mul();
void test_fast_add();
void test_atkin_morain();
void test_delta_precision();
void test_pow_bin();
void test_pow_no_bin();
void test_div();
void test_exponential();
void test_struct();
void test_modify();
void test_bit_operation();
void test_probabile_prime();
void test_factorization();
void test_factorize();
void test_elliptic();
void test_radix_gmp();
void test_quicksort();
void test_new_mul();
void test_new_add();
void test_step_three();
void test_one_root_finding();
void test_poly_iterative_pow();
void test_polynomial_radix();
void test_new_bit_operation();
void test_simple_mul();
void test_prova();

int main()
{
	mpf_set_default_prec(PRECISION);
	//printf("PRECISION: %d\n", mpf_get_default_prec());
	/* **********TIMING********** */
	time_t start, stop;
	clock_t ticks;
	time(&start);
	/* **********TIMING********** */

	test_atkin_morain();
	//test_fast_mul();
	//test_gcd();
	//test_fast_add();
	//test_fast_mul();
	//test_roots();	
	//test_hilbert();	
	//test_cornacchia();
	//test_delta_precision();
	//test_pow_no_bin();
	//test_pow_bin();
	//test_div();
	//test_exponential();
	//test_struct();
	//test_bit_operation();
	//test_probabile_prime();
	//test_factorization();
	//test_factorize();
	//test_elliptic();
	//test_radix_gmp();
	//test_quicksort();
	//test_new_mul();
	//test_new_add();
	//test_step_three();
	//test_one_root_finding();
	//test_poly_iterative_pow();
	//test_polynomial_radix();
	//test_new_bit_operation();
	//test_simple_mul();
	//test_prova();

	/* **********TIMING********** */
	ticks = clock();
	time(&stop);
	printf("Used %0.2f seconds of CPU time. \n", (double)ticks/CLOCKS_PER_SEC);
	printf("Finished in about %.0f seconds. \n", difftime(stop, start));
	/* **********TIMING********** */

	if(!AUTO_TIMING)
		printf("\nTIMER: %d\n\n", WAITING_TIME); //vedi factorize.h
	else
		printf("\nTIMER AUTOMATICO\n");

	return 0;
}

void test_atkin_morain()
{
	mpz_t p;
	mpz_init(p);	

	//mpz_init_set_str (p, "319705304701141539155720137200974664666792526059405792539680974929469783512821793995613718943171723765238853752439032835985158829038528214925658918372196742089464683960239919950882355844766055365179937610326127675178857306260955550407044463370239890187189 750909036833976197804646589380690779463976173", 10);

	mpz_init_set_str (p, "745213698191737003631319694753125429293968166002970537936165661845575001172678049743806549549977234670072449443569701103", 10);
	
	Atkin_Morain(p);	
}

void test_prova()
{
	mpz_t n;
	mpz_init(n);

	mpz_init_set_str (n, "71", 10); //coefficienti modulo n //prova 13, 7, 2
	gmp_printf("CALCOLI MODULO: %Zd\n", n);	

	r_big_polynomial poly_A;

	int num_el_A = 3;

	mpz_t* coef_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	//x^5 + 4x^4 + 2x^3 - 4x^2 - 3x = x(x-1)(x+1)(x+1)(x+3)

	mpz_init_set_str(coef_A[0], "1", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str(coef_A[1], "2", 10);
	mpz_init_set_str(coef_A[2], "6", 10);

	mpz_t* deg_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(deg_A[0], "1", 10);
	mpz_init_set_str(deg_A[1], "2", 10);
	mpz_init_set_str(deg_A[2], "3", 10);	

	init_real_big_polynomial(coef_A, num_el_A, deg_A, num_el_A, &poly_A, n);

	printf("POLINOMIO A\n");
	r_big_polynomial_printf(&poly_A);
	printf("POLINOMIO A\n");


	if(mpz_cmp_ui(poly_A.degrees[0], 0) == 0) 
	{
		mpz_sub_ui(poly_A.coefficients[0], poly_A.coefficients[0], 1);
		mpz_mod(poly_A.coefficients[0], poly_A.coefficients[0], n);
	}
	else
	{
		poly_A.num_elements = poly_A.num_elements + 1;
		poly_A.coefficients = (mpz_t*)realloc(poly_A.coefficients, sizeof(mpz_t)*poly_A.num_elements);
		poly_A.degrees = (mpz_t*)realloc(poly_A.degrees, sizeof(mpz_t)*poly_A.num_elements);

		mpz_init(poly_A.coefficients[poly_A.num_elements - 1]);
		mpz_init(poly_A.degrees[poly_A.num_elements - 1]);

		int i;
		for(i = poly_A.num_elements - 2; i >=0; i--)
		{
			mpz_set(poly_A.coefficients[i+1], poly_A.coefficients[i]);
			mpz_set(poly_A.degrees[i+1], poly_A.degrees[i]);
		}
		mpz_set_si(poly_A.coefficients[0], -1);
		mpz_mod(poly_A.coefficients[0], poly_A.coefficients[0], n);
		mpz_set_ui(poly_A.degrees[0], 0);
	}


	printf("POLINOMIO A\n");
	r_big_polynomial_printf(&poly_A);
	printf("POLINOMIO A\n");

}

void test_simple_mul()
{
	mpz_t n;
	mpz_init(n);

	mpz_init_set_str (n, "7000000000000000000000000000000000", 10); //coefficienti modulo n //prova 13, 7, 2
	gmp_printf("CALCOLI MODULO: %Zd\n", n);	

	r_big_polynomial poly_A;

	int num_el_A = 3;

	mpz_t* coef_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	//x^5 + 4x^4 + 2x^3 - 4x^2 - 3x = x(x-1)(x+1)(x+1)(x+3)

	mpz_init_set_str(coef_A[0], "3", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str(coef_A[1], "2", 10);
	mpz_init_set_str(coef_A[2], "6", 10);

	mpz_t* deg_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(deg_A[0], "0", 10);
	mpz_init_set_str(deg_A[1], "2", 10);
	mpz_init_set_str(deg_A[2], "3", 10);	

	init_real_big_polynomial(coef_A, num_el_A, deg_A, num_el_A, &poly_A, n);

	printf("POLINOMIO A\n");
	r_big_polynomial_printf(&poly_A);
	printf("POLINOMIO A\n");

	r_big_polynomial poly_B;

	int num_el_B = 3;

	mpz_t* coef_B = (mpz_t*)malloc(sizeof(mpz_t)*num_el_B);

	mpz_init_set_str (coef_B[0], "4", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str (coef_B[1], "3", 10);
	mpz_init_set_str (coef_B[2], "11", 10);

	mpz_t* deg_B = (mpz_t*)malloc(sizeof(mpz_t)*num_el_B);

	mpz_init_set_str (deg_B[0], "0", 10);
	mpz_init_set_str (deg_B[1], "2", 10);
	mpz_init_set_str (deg_B[2], "3", 10);

	init_real_big_polynomial(coef_B, num_el_B, deg_B, num_el_B, &poly_B, n);
	
	printf("POLINOMIO B\n");
	r_big_polynomial_printf(&poly_B);
	printf("POLINOMIO B\n");

	r_big_polynomial poly_mul;

	r_big_polynomial_simple_mul(&poly_mul, poly_A, poly_B, n);

	printf("POLINOMIO MUL\n");
	r_big_polynomial_printf(&poly_mul);
	printf("POLINOMIO MUL\n");

	r_big_polynomial poly_add;

	r_big_polynomial_simple_add(&poly_add, poly_A, poly_B, n);

	printf("POLINOMIO ADD\n");
	r_big_polynomial_printf(&poly_add);
	printf("POLINOMIO ADD\n");

}

void test_polynomial_radix()
{
	mpz_t n;
	mpz_init(n);
	mpz_init_set_str (n, "1349912222257818701234", 10); //coefficienti modulo n //prova 13, 7, 2
	gmp_printf("CALCOLI MODULO: %Zd\n", n);	

	r_big_polynomial poly_A;

	int num_el_A = 4;

	mpz_t* coef_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(coef_A[0], "5", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str(coef_A[1], "-10", 10);
	mpz_init_set_str(coef_A[2], "4", 10);
	mpz_init_set_str(coef_A[3], "2", 10);

	mpz_t* deg_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(deg_A[0], "0", 10);
	mpz_init_set_str(deg_A[1], "1", 10);
	mpz_init_set_str(deg_A[2], "2", 10);
	mpz_init_set_str(deg_A[3], "3", 10);

	init_real_big_polynomial(coef_A, num_el_A, deg_A, num_el_A, &poly_A, n);

	printf("POLINOMIO A\n");
	r_big_polynomial_printf(&poly_A);
	printf("POLINOMIO A\n");

	mpz_t root;
	mpz_init(root);

	find_one_root(root, poly_A, n);

	gmp_printf("root: %Zd\n", root);

	/*s_root_set root_set;	

	roots(&root_set, poly_A, n);

	if(root_set.num_elements == 0)
		printf("No radici intere\n");	
	else
	{
		printf("Radici:\n");
		int i;
		for(i = 0; i < root_set.num_elements; i++)
			gmp_printf("%Zd\n", root_set.roots[i]);		
	}*/

}

void test_poly_iterative_pow()
{
	mpz_t n;
	mpz_init(n);

	mpz_init_set_str (n, "7", 10); //coefficienti modulo n //prova 13, 7, 2
	gmp_printf("CALCOLI MODULO: %Zd\n", n);	

	r_big_polynomial poly_A;

	int num_el_A = 2;

	mpz_t* coef_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(coef_A[0], "6", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str(coef_A[1], "1", 10);

	mpz_t* deg_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(deg_A[0], "0", 10);
	mpz_init_set_str(deg_A[1], "1", 10);

	init_real_big_polynomial(coef_A, num_el_A, deg_A, num_el_A, &poly_A, n);

	printf("POLINOMIO A\n");
	r_big_polynomial_printf(&poly_A);
	printf("POLINOMIO A\n");

	r_big_polynomial poly_pow;
	r_big_polynomial_copy(&poly_pow, poly_A);

	mpz_t exp;
	mpz_init_set_str(exp, "377762837468273", 10);

	r_big_polynomial_iterative_pow(&poly_pow, exp, n);

	printf("POLINOMIO POW\n");
	r_big_polynomial_printf(&poly_pow);
	printf("POLINOMIO POW\n");
}

void test_one_root_finding()
{
	mpz_t n;
	mpz_init(n);

	mpz_init_set_str (n, "41", 10); //coefficienti modulo n //prova 13, 7, 2 PER 41 LA TROVA (in quanto sarebbe zero, e il termine noto vale 0 mod 41), MA NON E' DETTO LA TROVI PER ALTRI MODULI...
	gmp_printf("CALCOLI MODULO: %Zd\n", n);	

	r_big_polynomial poly_A;

	int num_el_A = 4;

	mpz_t* coef_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(coef_A[0], "4671133182399954782798673154437441310949376", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str(coef_A[1], "1363005552434666078211357566795785667232380388321042307465", 10);
	mpz_init_set_str(coef_A[2], "3005101108071026200706725969920", 10);
	mpz_init_set_str(coef_A[3], "1", 10);

	mpz_t* deg_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(deg_A[0], "0", 10);
	mpz_init_set_str(deg_A[1], "1", 10);
	mpz_init_set_str(deg_A[2], "2", 10);
	mpz_init_set_str(deg_A[3], "3", 10);

	init_real_big_polynomial(coef_A, num_el_A, deg_A, num_el_A, &poly_A, n);

	printf("POLINOMIO A\n");
	r_big_polynomial_printf(&poly_A);
	printf("POLINOMIO A\n");

	mpz_t one_root;
	mpz_init(one_root);

	int size = 7;

	mpz_t primeDivisors[7];
	mpz_init_set_ui(primeDivisors[0], 2);
	mpz_init_set_ui(primeDivisors[1], 3);
	mpz_init_set_ui(primeDivisors[2], 17);
	mpz_init_set_ui(primeDivisors[3], 23);
	mpz_init_set_ui(primeDivisors[4], 41);
	mpz_init_set_ui(primeDivisors[5], 71);
	mpz_init_set_ui(primeDivisors[6], 83);

	int multiplicity[7];
	multiplicity[0] = 48;
	multiplicity[1] = 9;
	multiplicity[2] = 3;
	multiplicity[3] = 3;
	multiplicity[4] = 3;
	multiplicity[5] = 3;
	multiplicity[6] = 3;

	mpz_t currentResult;
	mpz_init_set_ui(currentResult, 1);

	int success = findOneRootInDivisors(one_root, poly_A, n, primeDivisors, multiplicity, size, 0, currentResult);

	printf("Success: %d\n", success);
	gmp_printf("one_root: %Zd\n", one_root);


	/*mpz_t n;
	mpz_init(n);

	mpz_init_set_str (n, "41", 10); //coefficienti modulo n //prova 13, 7, 2
	gmp_printf("CALCOLI MODULO: %Zd\n", n);	

	r_big_polynomial poly_A;

	int num_el_A = 4;

	mpz_t* coef_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(coef_A[0], "21", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str(coef_A[1], "31", 10);
	mpz_init_set_str(coef_A[2], "11", 10);
	mpz_init_set_str(coef_A[3], "1", 10);

	mpz_t* deg_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(deg_A[0], "0", 10);
	mpz_init_set_str(deg_A[1], "1", 10);
	mpz_init_set_str(deg_A[2], "2", 10);
	mpz_init_set_str(deg_A[3], "3", 10);

	init_real_big_polynomial(coef_A, num_el_A, deg_A, num_el_A, &poly_A, n);

	printf("POLINOMIO A\n");
	r_big_polynomial_printf(&poly_A);
	printf("POLINOMIO A\n");

	mpz_t one_root;
	mpz_init(one_root);

	int size = 2;

	mpz_t primeDivisors[2];
	mpz_init_set_ui(primeDivisors[0], 2);
	mpz_init_set_ui(primeDivisors[1], 3);

	int multiplicity[2];
	multiplicity[0] = 2;
	multiplicity[1] = 1;

	mpz_t currentResult;
	mpz_init_set_ui(currentResult, 1);

	int success = findOneRootInDivisors(one_root, poly_A, n, primeDivisors, multiplicity, size, 0, currentResult);

	printf("Success: %d\n", success);
	gmp_printf("one_root: %Zd\n", one_root);*/
}

void test_step_three()
{
	mpz_t p;
	mpz_init(p);
	mpz_init_set_str (p, "1363005552434666078217421284621279933627102780881053358473", 10);

	mpz_t D;
	mpz_init(D);
	mpz_init_set_str(D, "-499", 10);

	curve_parameters c_parameters;
	
	int ret = step_three(p, D, &c_parameters);

	printf("ret: %d\n", ret);
}

void test_quicksort()
{
	int size = 9;
	mpz_t* deg = (mpz_t*)malloc(sizeof(mpz_t)*size);

	mpz_init_set_ui(deg[0], 2);
	mpz_init_set_ui(deg[1], 5);
	mpz_init_set_ui(deg[2], 2);
	mpz_init_set_ui(deg[3], 5);
	mpz_init_set_ui(deg[4], 8);
	mpz_init_set_ui(deg[5], 5);
	mpz_init_set_ui(deg[6], 7);
	mpz_init_set_ui(deg[7], 9);
	mpz_init_set_ui(deg[8], 7);

	quicksort_bottom_up(deg, 0, size - 1) ;

	int i;
	for(i = 0; i < size; i++)
		gmp_printf("deg[%d]: %Zd\n", i, deg[i]);

}

void test_radix_gmp()
{
	mpz_t p;
	mpz_init(p);
	//mpz_init_set_str (p, "11", 10);
	mpz_init_set_str (p, "393050634124102232869567034555427371542904833", 10); //1009
	gmp_printf("CALCOLI MODULO: %Zd\n", p);

	mpz_t Q;
	mpz_init(Q);
	//mpz_init_set_str (Q, "1", 10);
	mpz_init_set_str (Q, "-3", 10);
	gmp_printf("Q: %Zd\n", Q);

	mpz_t y;
	mpz_init(y);

	int ret = radix_GMP(y, p, Q);

	if(ret != 0)
		gmp_printf("y: %Zd\n", y);
	else
		printf("Oh nooooooo\n");

}

void test_elliptic()
{
	mpz_t p;
	mpz_init(p);

	mpz_init_set_str (p, "618970019642690137449562111", 10);
	gmp_printf("CALCOLI MODULO: %Zd\n", p);

	mpz_t A; 
	mpz_init_set_str(A, "0", 10); //X^3 + Ax + B  A = -1  B = 1

	
	curve_point P;
	init(&P);

	mpz_set_str(P.x, "2", 10);
	mpz_set_str(P.y, "2", 10);
	mpz_set_str(P.z, "1", 10);

	curve_point Q;
	init(&Q);

	int ret = doubleP(&Q, P, p, A);	
	gmp_printf("DOUBLE P(%Zd, %Zd, %Zd): \n", P.x, P.y, P.z);

	if(ret == 0)
	{
		gmp_printf("\tQ.x: %Zd\n", Q.x);
		gmp_printf("\tQ.y: %Zd\n", Q.y);
		gmp_printf("\tQ.z: %Zd\n", Q.z);
	}
	else if(ret == 1)
		printf("\tOttenuto punto all'infinito\n");
	else if(ret == -1)
		printf("\tInversione impossibile\n");

	curve_point P1;
	init(&P1);

	mpz_set_str(P1.x, "2", 10);
	mpz_set_str(P1.y, "2", 10);
	mpz_set_str(P1.z, "1", 10);

	curve_point P2;
	init(&P2);

	mpz_set_str(P2.x, "32", 10);
	mpz_set_str(P2.y, "3", 10);
	mpz_set_str(P2.z, "1", 10);

	ret = add(&Q, P1, P2, p, A);
	gmp_printf("P1(%Zd, %Zd, %Zd) + P2(%Zd, %Zd, %Zd): \n", P1.x, P1.y, P1.z, P2.x, P2.y, P2.z);

	if(ret == 0)
	{
		gmp_printf("\tQ.x: %Zd\n", Q.x);
		gmp_printf("\tQ.y: %Zd\n", Q.y);
		gmp_printf("\tQ.z: %Zd\n", Q.z);
	}
	else if(ret == 1)
		printf("\tOttenuto punto all'infinito\n");
	else if(ret == -1)
		printf("\tInversione impossibile\n");

	curve_point P3;
	init(&P3);

	mpz_set_str(P3.x, "174298697913514185939789940", 10);
	mpz_set_str(P3.y, "373006277778030021814199169", 10);
	mpz_set_str(P3.z, "1", 10);

	mpz_t k;
	mpz_init_set_str(k, "279057964", 10);

	ret = mul(&Q, P3, k, p, A);
	gmp_printf("[%Zd]*P(%Zd, %Zd, %Zd): \n", k, P3.x, P3.y, P3.z);

	if(ret == 0)
	{
		gmp_printf("\tQ.x: %Zd\n", Q.x);
		gmp_printf("\tQ.y: %Zd\n", Q.y);
		gmp_printf("\tQ.z: %Zd\n", Q.z);
	}
	else if(ret == 1)
		printf("\tOttenuto punto all'infinito\n");
	else if(ret == -1)
		printf("\tInversione impossibile\n");	


}

void test_hilbert()
{
	Hresult* res = malloc(sizeof(Hresult));	
	mpz_t n;
	mpz_init_set_str(n, "618970019642690137449562111", 10);
	gmp_printf("MODULO %Zd\n", n);
	mpz_t D;
	mpz_init_set_si(D, -512);
	gmp_printf("DISCRIMINANTE: %Zd\n", D);
	
	hilbertclass(res, D, n);
	
	printf("---------------\n");
	int i;
	for(i = 0; i < res->count; i++)
	{
		gmp_printf("A: %Zd\n", res->abc_set[i].a);
		gmp_printf("B: %Zd\n", res->abc_set[i].b);
		gmp_printf("C: %Zd\n\n", res->abc_set[i].c);
	}
}

void test_factorize()
{
	/*mpz_t p;
	mpz_init_set_str(p, "599444814781680415023062503435349875984323897837383", 10);

	gmp_printf("p: %Zd\n\n", p);

	prime_factors p_factors;

	int set_size = 30;
	p_factors.factors = (mpz_t*)malloc(sizeof(mpz_t)*set_size);	
	int i;
	for(i = 0; i < set_size; i++)
		mpz_init(p_factors.factors[i]);

	p_factors.num_elements = 0;

	factor(&p_factors, p, 20);


	printf("\n------------------------------------\n");
	for(i = 0; i < p_factors.num_elements; i++)
		gmp_printf("%Zd\n", p_factors.factors[i]);*/

}

void test_factorization()
{
	mpz_t n;
	mpz_init_set_str(n, "12", 10);
	mpz_t f;
	mpz_init_set_str(f, "3", 10);
	mpz_t result;
	mpz_init(result);

	int num_occur = mpz_remove(result, n, f);
	gmp_printf("N: %Zd\n", n);
	gmp_printf("Numero occorrenze di %Zd in %Zd: %d\n", f, n, num_occur);
	gmp_printf("N senza fattori: %Zd\n", result);

}

void test_probabile_prime()
{
	mpz_t n;
	mpz_init_set_str(n, "618970019642690137449562111", 10);

	int NUM_TIMES = 10;

	int ret = mpz_probab_prime_p(n, NUM_TIMES);

	printf("ret: %d\n", ret);
}

//funzionano anche con indice di bit superiore al numero di bit corrente!!!!!!!!!!!!
void test_bit_operation()
{
	mpz_t n;
	mpz_init_set_str(n, "8", 10);

	int num_bit = mpz_sizeinbase(n, 2);
	gmp_printf("num_bit(%Zd): %d\n", n, num_bit);

	mpz_setbit(n, 3); //clrbit tstbit combit
	gmp_printf("n: %Zd\n", n);

	int ret = mpz_tstbit(n, 4);
	printf("ret test: %d\n", ret);
}

void test_struct()
{
	mpz_t n;
	mpz_init_set_str(n, "4", 10);
	test_modify(n);
	gmp_printf("%Zd\n", n);
}

void test_modify(mpz_t n)
{
	mpz_set_str(n, "400", 10);
}

void test_exponential()
{
	complex_f result;
	complex_finit(&result);
	complex_f z;
	complex_finit(&z);

	mpf_set_str(z.real, "2", 10);
	mpf_set_str(z.img, "3", 10);


	
	double complex a = 2 + 3*I;
	exponential(&result, z);
	
	printf("COM:\n\treal->%.*f\n\timg->%.*f\n", 100, creal(cexp(a)), 100, cimag(cexp(a)));
	gmp_printf("MIO:\n\treal->%.*Ff\n\timg->%.*Ff\n", 100, result.real, 100, result.img);

}

void test_pow_no_bin()
{
	struct complex_f com;
	complex_finit(&com);	
	mpf_set_str(com.real, "2.3123456", 10);
	mpf_set_str(com.img, "3.3123456", 10);

	struct complex_f result; //uso result per evitare di sovrascrivere la variabile originaria
	complex_finit(&result);	
	mpf_set(result.real, com.real);
	mpf_set(result.img, com.img);

	int exp = 3;

	pow_no_bin(&result, com, exp);

	gmp_printf("result.real: %Ff\n", result.real);	
	gmp_printf("result.img: %Ff\n", result.img);

	gmp_printf("com.real: %Ff\n", com.real);	
	gmp_printf("com.img: %Ff\n", com.img);
}

void test_div()
{
	struct complex_f z1;
	complex_finit(&z1);	
	mpf_set_str(z1.real, "0", 10);
	mpf_set_str(z1.img, "4", 10);	

	struct complex_f z2;
	complex_finit(&z2);
	mpf_set_str(z2.real, "4", 10);
	mpf_set_str(z2.img, "0.0", 10);

	struct complex_f z3;
	complex_finit(&z3);	

	div_fcomplex(&z3, z1, z2);

	gmp_printf("z3.real: %Ff\n", z3.real);	
	gmp_printf("z3.img: %Ff\n", z3.img);
}

void test_pow_bin()
{
	struct complex_f res_pow_bin;
	complex_finit(&res_pow_bin);

	struct complex_f com;
	//mpf_init(com.real);
	//mpf_init(com.img);
	complex_finit(&com);	

	gmp_printf("com.real: %.*Ff\n", 100, com.real);
	gmp_printf("com.img: %.*Ff\n", 100, com.img);


	//mpf_set_d(com.real, 2.3123456);
	//mpf_set_d(com.img, 3.3123456);
	mpf_set_str(com.real, "2.3123456", 10);
	mpf_set_str(com.img, "3.3123456", 10);

	gmp_printf("com.real: %.*Ff\n", 100, com.real);
	gmp_printf("com.img: %.*Ff\n", 100, com.img);

	powbin(&res_pow_bin, com, 40);
	gmp_printf("POWBIN_REAL: %.*Ff\n", 300, res_pow_bin.real);
	gmp_printf("POWBIN_IMG: %.*Ff\n\n\n\n", 300, res_pow_bin.img);

	/*char* a = (char*)malloc(sizeof(char)*10);
	double b = 1.2345;
	sprintf(a, "%f\n", b);	
	printf("stringa: %s\n", a);*/
}

void test_delta_precision()
{	
	struct complex_f c, result;
	
	complex_finit(&c);
	complex_finit(&result);

	mpf_set_str(c.real, "0.5", 10);
	mpf_set_str(c.img, "0.1", 10);

	gmp_printf("real: %.*Ff\n\n\n\n", 50, c.real);
	gmp_printf("img: %.*Ff\n\n\n\n", 50, c.img);

	//mpf_set_d(c.real, 0.00000000002701383278246102271392992686439995337830);
	//mpf_set_d(c.img, 0.00000000000000000000000000661626225329939945883642);

	/////////////7//deltaprecision(cexp(c));
	deltaprecision(&result, c);
	gmp_printf("real: %.*Ff\n\n\n\n", 50, result.real);
	gmp_printf("img: %.*Ff\n\n\n\n", 50, result.img);
	printf("Buon compleanno Mr. Pup!!!!\n");
}

int test_cornacchia()
{
	mpz_t p;
	mpz_t D;
	mpz_t radix;

	mpz_init_set_str (p, "170141183460469231731687303715884105727", 10); //618970014642690137449562111 che stà sul sito non funzione per via di jacobi(D/p) --> credo sia toppato sul sito!
	mpz_init_set_str (p, "9", 10);
	mpz_init_set_si(D, -4);
	mpz_init(radix);

	struct resultset r;

	cornacchia_GMP(&r, p, D);

	int i;
	for(i = 0; i < r.num_elements; i++)
		gmp_printf("r[%d] %Zd\n", i, r.results[i]);
	return 0;

}


int test_gcd()
{
	mpz_t n;
	mpz_init(n);
	mpz_init_set_str (n, "2", 10); //coefficienti modulo n //prova 13, 7, 2
	gmp_printf("CALCOLI MODULO: %Zd\n", n);	

	r_big_polynomial poly_A;

	int num_el_A = 4;

	mpz_t* coef_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str (coef_A[0], "1", 10); //i coefficienti sono ridotti modulo n!!!!!!!!! //-1
	mpz_init_set_str (coef_A[1], "7", 10);
	mpz_init_set_str (coef_A[2], "1", 10);
	mpz_init_set_str (coef_A[3], "7", 10);
	//mpz_init_set_str (coef_A[4], "7", 10);

	mpz_t* deg_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str (deg_A[0], "0", 10);
	mpz_init_set_str (deg_A[1], "2", 10);
	mpz_init_set_str (deg_A[2], "9", 10);
	mpz_init_set_str (deg_A[3], "11", 10);
	//mpz_init_set_str (deg_A[4], "6", 10);

	int c_size_A = num_el_A;
	int d_size_A = c_size_A; //dimensione array gradi

	init_real_big_polynomial(coef_A, c_size_A, deg_A, d_size_A, &poly_A, n);

	printf("POLINOMIO A\n");
	r_big_polynomial_printf(&poly_A);
	printf("POLINOMIO A\n");

	r_big_polynomial poly_B;

	int num_el_B = 4;

	mpz_t* coef_B = (mpz_t*)malloc(sizeof(mpz_t)*num_el_B);

	mpz_init_set_str (coef_B[0], "1", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str (coef_B[1], "7", 10);
	mpz_init_set_str (coef_B[2], "-1", 10);
	mpz_init_set_str (coef_B[3], "-7", 10);

	mpz_t* deg_B = (mpz_t*)malloc(sizeof(mpz_t)*num_el_B);

	mpz_init_set_str (deg_B[0], "0", 10);
	mpz_init_set_str (deg_B[1], "2", 10);
	mpz_init_set_str (deg_B[2], "5", 10);
	mpz_init_set_str (deg_B[3], "7", 10);

	int c_size_B = num_el_B;
	int d_size_B = c_size_B; //dimensione array gradi

	init_real_big_polynomial(coef_B, c_size_B, deg_B, d_size_B, &poly_B, n);
	
	printf("POLINOMIO B\n");
	r_big_polynomial_printf(&poly_B);
	printf("POLINOMIO B\n");	

	r_big_polynomial poly_gcd;
	r_big_polynomial_gcd(&poly_gcd, poly_A, poly_B, n);

	printf("POLINOMIO GCD\n");
	r_big_polynomial_printf(&poly_gcd);
	printf("POLINOMIO GCD\n");
	
	
	/*r_big_polynomial poly_rem;
	r_big_polynomial poly_quoz;

	r_big_polynomial_fast_mod(&poly_rem, &poly_quoz, poly_A, poly_B, n);

	printf("POLINOMIO REMAINDER\n");
	r_big_polynomial_printf(&poly_rem);
	printf("POLINOMIO REMAINDER\n");	

	printf("POLINOMIO QUOZIENT\n");
	r_big_polynomial_printf(&poly_quoz);
	printf("POLINOMIO QUOZIENT\n");	

	gmp_printf("CALCOLI MODULO: %Zd\n", n);*/


	/*r_big_polynomial poly_mul;

	r_big_polynomial_mul(&poly_mul, poly_A, poly_B, n);

	printf("POLINOMIO MUL\n");
	r_big_polynomial_printf(&poly_mul);
	printf("POLINOMIO MUL\n");*/

	/*
	mpz_t degree;
	mpz_init_set_ui(degree, 8);

	r_big_polynomial poly_trunc_recip;

	r_big_polynomial_fast_inversion(&poly_trunc_recip, poly_A, degree, n);

	printf("POLINOMIO C\n");
	r_big_polynomial_printf(&poly_trunc_recip);
	printf("POLINOMIO C\n");*/


	return 0;
}

void test_roots()
{
	mpz_t n;
	mpz_init(n);

	//mpz_init_set_str (n, "7111122233", 10); //coefficienti modulo n, OCCHIO CHE DEVE ESSERE DISPARI!!
	mpz_init_set_str (n, "3567300520022142673540701160992979258890578006839077408987159726327696072974110783861838260391358445311165633846403085217066652747473563968420047909272241155503481353", 10);
	gmp_printf("CALCOLI MODULO: %Zd\n", n);	

	r_big_polynomial poly_A;

	int num_el_A = 4;

	mpz_t* coef_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(coef_A[0], "12771880859375", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str(coef_A[1], "-5151296875", 10);
	mpz_init_set_str(coef_A[2], "3491750", 10);
	mpz_init_set_str(coef_A[3], "1", 10);

	//x^5 + 4x^4 + 2x^3 - 4x^2 - 3x = x(x-1)(x+1)(x+1)(x+3)
	/*mpz_init_set_str(coef_A[0], "-3", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str(coef_A[1], "-4", 10);
	mpz_init_set_str(coef_A[2], "2", 10);
	mpz_init_set_str(coef_A[3], "4", 10);
	mpz_init_set_str(coef_A[4], "1", 10);*/

	mpz_t* deg_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(deg_A[0], "0", 10);
	mpz_init_set_str(deg_A[1], "1", 10);
	mpz_init_set_str(deg_A[2], "2", 10);
	mpz_init_set_str(deg_A[3], "3", 10);

	/*mpz_init_set_str(deg_A[0], "1", 10);
	mpz_init_set_str(deg_A[1], "2", 10);
	mpz_init_set_str(deg_A[2], "3", 10);
	mpz_init_set_str(deg_A[3], "4", 10);
	mpz_init_set_str(deg_A[4], "5", 10);*/

	init_real_big_polynomial(coef_A, num_el_A, deg_A, num_el_A, &poly_A, n);

	printf("POLINOMIO DI CUI TROVARE LE RADICI\n");
	r_big_polynomial_printf(&poly_A);
	printf("POLINOMIO DI CUI TROVARE LE RADICI\n");

	s_root_set root_set;	
	getchar();
	roots(&root_set, poly_A, n);
	
	if(root_set.num_elements == 0)
		printf("NESSUNA RADICE :( FORSE COMPLESSE O FLOAT ;)\n");
	else
	{
		printf("LA MOLTEPLICITA' DELLE RADICI NON E' RIPORTATA (PER INTERO) - NON RIMUOVO EVENTUALI DOPPIONI\n");
		int i;
		for(i = 0; i < root_set.num_elements; i++)
		{
			gmp_printf("radix_set[%d]: %Zd\n", i, root_set.roots[i]);
		}
	}
}

void test_pow()
{
	mpz_t n;
	mpz_init(n);
	mpz_init_set_str (n, "7", 10);
	gmp_printf("CALCOLI MODULO: %Zd\n", n);	

	r_big_polynomial poly_A;

	int num_el_A = 2;
	mpz_t* coef_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);
	mpz_init_set_str (coef_A[0], "4", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str (coef_A[1], "2", 10);

	mpz_t* deg_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);
	mpz_init_set_str (deg_A[0], "0", 10);
	mpz_init_set_str (deg_A[1], "1", 10);

	init_real_big_polynomial(coef_A, num_el_A, deg_A, num_el_A, &poly_A, n);

	printf("POLINOMIO A\n");
	r_big_polynomial_printf(&poly_A);
	printf("POLINOMIO A\n");

	r_big_polynomial poly_pow;
	r_big_polynomial_copy(&poly_pow, poly_A);

	mpz_t exp;
	mpz_init_set_ui(exp, 3);
	
	r_big_polynomial_recursive_pow(&poly_pow, exp, n);

	printf("POLINOMIO POW\n");
	r_big_polynomial_printf(&poly_pow);
	printf("POLINOMIO POW\n");
}

void test_fast_mod()
{
	mpz_t n;
	mpz_init(n);
	mpz_init_set_str (n, "7", 10); //coefficienti modulo n //prova 13, 7, 2
	gmp_printf("CALCOLI MODULO: %Zd\n", n);	

	r_big_polynomial poly_A;

	int num_el_A = 4;

	mpz_t* coef_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str (coef_A[0], "1", 10); //i coefficienti sono ridotti modulo n!!!!!!!!! //-1
	mpz_init_set_str (coef_A[1], "2", 10);
	mpz_init_set_str (coef_A[2], "15", 10);
	mpz_init_set_str (coef_A[3], "10", 10);
	//mpz_init_set_str (coef_A[4], "7", 10);

	mpz_t* deg_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str (deg_A[0], "0", 10);
	mpz_init_set_str (deg_A[1], "1", 10);
	mpz_init_set_str (deg_A[2], "3", 10);
	mpz_init_set_str (deg_A[3], "5", 10);
	//mpz_init_set_str (deg_A[4], "6", 10);

	init_real_big_polynomial(coef_A, num_el_A, deg_A, num_el_A, &poly_A, n);

	printf("POLINOMIO DIVIDENDO\n");
	r_big_polynomial_printf(&poly_A);
	printf("POLINOMIO DIVIDENDO\n");

	r_big_polynomial poly_B;

	int num_el_B = 2;

	mpz_t* coef_B = (mpz_t*)malloc(sizeof(mpz_t)*num_el_B);

	mpz_init_set_str (coef_B[0], "3", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str (coef_B[1], "2", 10);
	//mpz_init_set_str (coef_B[2], "-1", 10);
	//mpz_init_set_str (coef_B[3], "-7", 10);

	mpz_t* deg_B = (mpz_t*)malloc(sizeof(mpz_t)*num_el_B);

	mpz_init_set_str (deg_B[0], "0", 10);
	mpz_init_set_str (deg_B[1], "2", 10);
	//mpz_init_set_str (deg_B[2], "5", 10);
	//mpz_init_set_str (deg_B[3], "7", 10);

	init_real_big_polynomial(coef_B, num_el_B, deg_B, num_el_B, &poly_B, n);
	
	printf("POLINOMIO DIVISORE\n");
	r_big_polynomial_printf(&poly_B);
	printf("POLINOMIO DIVISORE\n");

	r_big_polynomial poly_rem;
	r_big_polynomial poly_quoz;

	r_big_polynomial_fast_mod(&poly_rem, &poly_quoz, poly_A, poly_B, n);
	
	printf("POLINOMIO QUOZIENTE\n");
	r_big_polynomial_printf(&poly_quoz);
	printf("POLINOMIO QUOZIENTE\n");

	printf("POLINOMIO RESTO\n");
	r_big_polynomial_printf(&poly_rem);
	printf("POLINOMIO RESTO\n");
	

}

void test_fast_mul()
{
	mpz_t n;
	mpz_init(n);

	mpz_init_set_str (n, "7000000000000000000000000000000000", 10); //coefficienti modulo n //prova 13, 7, 2
	gmp_printf("CALCOLI MODULO: %Zd\n", n);	

	r_big_polynomial poly_A;

	int num_el_A = 2;

	mpz_t* coef_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	//x^5 + 4x^4 + 2x^3 - 4x^2 - 3x = x(x-1)(x+1)(x+1)(x+3)

	mpz_init_set_str(coef_A[0], "3", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str(coef_A[1], "2", 10);

	mpz_t* deg_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(deg_A[0], "0", 10);
	mpz_init_set_str(deg_A[1], "6000000", 10);

	init_real_big_polynomial(coef_A, num_el_A, deg_A, num_el_A, &poly_A, n);

	printf("POLINOMIO A\n");
	r_big_polynomial_printf(&poly_A);
	printf("POLINOMIO A\n");

	r_big_polynomial poly_B;

	int num_el_B = 2;

	mpz_t* coef_B = (mpz_t*)malloc(sizeof(mpz_t)*num_el_B);

	mpz_init_set_str (coef_B[0], "4", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str (coef_B[1], "3", 10);

	mpz_t* deg_B = (mpz_t*)malloc(sizeof(mpz_t)*num_el_B);

	mpz_init_set_str (deg_B[0], "0", 10);
	mpz_init_set_str (deg_B[1], "1", 10);

	init_real_big_polynomial(coef_B, num_el_B, deg_B, num_el_B, &poly_B, n);
	
	printf("POLINOMIO B\n");
	r_big_polynomial_printf(&poly_B);
	printf("POLINOMIO B\n");

	r_big_polynomial poly_mul;

	r_big_polynomial_fast_mul(&poly_mul, poly_A, poly_B, n);

	printf("POLINOMIO PRODOTTO\n");
	r_big_polynomial_printf(&poly_mul);
	printf("POLINOMIO PRODOTTO\n");
}

void test_fast_add()
{
	mpz_t n;
	mpz_init(n);

	mpz_init_set_str (n, "70000000000000000000000000000000000", 10); //coefficienti modulo n //prova 13, 7, 2
	gmp_printf("CALCOLI MODULO: %Zd\n", n);	

	r_big_polynomial poly_A;

	int num_el_A = 2;

	mpz_t* coef_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	//x^5 + 4x^4 + 2x^3 - 4x^2 - 3x = x(x-1)(x+1)(x+1)(x+3)

	mpz_init_set_str(coef_A[0], "3", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str(coef_A[1], "2", 10);

	mpz_t* deg_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(deg_A[0], "0", 10);
	mpz_init_set_str(deg_A[1], "700", 10);

	init_real_big_polynomial(coef_A, num_el_A, deg_A, num_el_A, &poly_A, n);

	printf("POLINOMIO A\n");
	r_big_polynomial_printf(&poly_A);
	printf("POLINOMIO A\n");

	r_big_polynomial poly_B;

	int num_el_B = 2;

	mpz_t* coef_B = (mpz_t*)malloc(sizeof(mpz_t)*num_el_B);

	mpz_init_set_str (coef_B[0], "5", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str (coef_B[1], "3", 10);

	mpz_t* deg_B = (mpz_t*)malloc(sizeof(mpz_t)*num_el_B);

	mpz_init_set_str (deg_B[0], "0", 10);
	mpz_init_set_str (deg_B[1], "1", 10);

	init_real_big_polynomial(coef_B, num_el_B, deg_B, num_el_B, &poly_B, n);
	
	printf("POLINOMIO B\n");
	r_big_polynomial_printf(&poly_B);
	printf("POLINOMIO B\n");

	r_big_polynomial poly_add;

	r_big_polynomial_fast_add(&poly_add, poly_A, poly_B, n);

	printf("POLINOMIO SOMMA\n");
	r_big_polynomial_printf(&poly_add);
	printf("POLINOMIO SOMMA\n");


}
/*
void test_new_mul()
{
	mpz_t n;
	mpz_init(n);

	mpz_init_set_str (n, "700000000000", 10); //coefficienti modulo n //prova 13, 7, 2
	gmp_printf("CALCOLI MODULO: %Zd\n", n);	

	r_big_polynomial poly_A;

	int num_el_A = 4;

	mpz_t* coef_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(coef_A[0], "8", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str(coef_A[1], "10", 10);
	mpz_init_set_str(coef_A[2], "5", 10);
	mpz_init_set_str(coef_A[3], "3", 10);	

	mpz_t* deg_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(deg_A[0], "0", 10);
	mpz_init_set_str(deg_A[1], "2", 10);
	mpz_init_set_str(deg_A[2], "5", 10);
	mpz_init_set_str(deg_A[3], "20", 10);

	init_real_big_polynomial(coef_A, num_el_A, deg_A, num_el_A, &poly_A, n);

	printf("POLINOMIO A\n");
	r_big_polynomial_printf(&poly_A);
	printf("POLINOMIO A\n");

	r_big_polynomial poly_B;

	int num_el_B = 3;

	mpz_t* coef_B = (mpz_t*)malloc(sizeof(mpz_t)*num_el_B);

	mpz_init_set_str (coef_B[0], "134", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str (coef_B[1], "100", 10);
	mpz_init_set_str (coef_B[2], "1", 10);

	mpz_t* deg_B = (mpz_t*)malloc(sizeof(mpz_t)*num_el_B);

	mpz_init_set_str (deg_B[0], "1", 10);
	mpz_init_set_str (deg_B[1], "8", 10);
	mpz_init_set_str (deg_B[2], "10", 10);

	init_real_big_polynomial(coef_B, num_el_B, deg_B, num_el_B, &poly_B, n);
	
	printf("POLINOMIO B\n");
	r_big_polynomial_printf(&poly_B);
	printf("POLINOMIO B\n");

	r_big_polynomial poly_mul;

	r_big_polynomial_new_mul(&poly_mul, poly_A, poly_B, n);

	printf("POLINOMIO PRODOTTO\n");
	r_big_polynomial_printf(&poly_mul);
	printf("POLINOMIO PRODOTTO\n");

	printf("POLINOMIO A\n");
	r_big_polynomial_printf(&poly_A);
	printf("POLINOMIO A\n");

	printf("POLINOMIO B\n");
	r_big_polynomial_printf(&poly_B);
	printf("POLINOMIO B\n");
}
*/
/*
void test_new_add()
{
	mpz_t n;
	mpz_init(n);

	mpz_init_set_str (n, "7", 10); //coefficienti modulo n //prova 13, 7, 2
	gmp_printf("CALCOLI MODULO: %Zd\n", n);	

	r_big_polynomial poly_A;

	int num_el_A = 3;

	mpz_t* coef_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(coef_A[0], "1", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str(coef_A[1], "1", 10);
	mpz_init_set_str(coef_A[2], "2", 10);	

	mpz_t* deg_A = (mpz_t*)malloc(sizeof(mpz_t)*num_el_A);

	mpz_init_set_str(deg_A[0], "1", 10);
	mpz_init_set_str(deg_A[1], "2", 10);
	mpz_init_set_str(deg_A[2], "5", 10);

	init_real_big_polynomial(coef_A, num_el_A, deg_A, num_el_A, &poly_A, n);

	printf("POLINOMIO A\n");
	r_big_polynomial_printf(&poly_A);
	printf("POLINOMIO A\n");

	r_big_polynomial poly_B;

	int num_el_B = 3;

	mpz_t* coef_B = (mpz_t*)malloc(sizeof(mpz_t)*num_el_B);

	mpz_init_set_str (coef_B[0], "1", 10); //i coefficienti sono ridotti modulo n!!!!!!!!!
	mpz_init_set_str (coef_B[1], "2", 10);
	mpz_init_set_str (coef_B[2], "1", 10);

	mpz_t* deg_B = (mpz_t*)malloc(sizeof(mpz_t)*num_el_B);

	mpz_init_set_str (deg_B[0], "0", 10);
	mpz_init_set_str (deg_B[1], "1", 10);
	mpz_init_set_str (deg_B[2], "2", 10);

	init_real_big_polynomial(coef_B, num_el_B, deg_B, num_el_B, &poly_B, n);
	
	printf("POLINOMIO B\n");
	r_big_polynomial_printf(&poly_B);
	printf("POLINOMIO B\n");

	r_big_polynomial poly_add;

	r_big_polynomial_new_add(&poly_add, poly_A, poly_B, n);

	printf("POLINOMIO SOMMA\n");
	r_big_polynomial_printf(&poly_add);
	printf("POLINOMIO SOMMA\n");


}
*/




