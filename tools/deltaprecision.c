#include "deltaprecision.h"

//La funzione calcola il fattoriale di un numero con GMP
void mpf_factorial(mpf_t res, mpf_t n)
{
	//gmp_printf("%Ff fattoriale: ", n);	

	if(mpf_cmp_ui(n, 0) == 0)
	{
		mpf_set_ui(res, 1);
		return;
	}

	mpf_t fact;
	mpf_init2(fact, PRECISION);
	unsigned long i;
	unsigned long nn;
	nn = mpf_get_ui(n);
	mpf_set_ui(fact, 1); //Initialise fact to 1, default is 0
	
	for (i = 1; i <= nn; i++) {
		mpf_mul_ui(fact, fact, i);
	}
		
	mpf_set(res, fact);
	//gmp_printf("%Ff\n", res);	

	mpf_clear(fact);
	
	return;
}


//La funzione è un implementazione del Teorema Binomiale che 
//esprime lo sviluppo della potenza n-esima di un binomio 
//operando e risultato possono essere uguali
void powbin(struct complex_f* result, struct complex_f c, unsigned long n){

	mpf_t a, b, ab, countreal, countimg;
	mpf_init2(a, PRECISION);
	mpf_init2(b, PRECISION);
	mpf_init2(countreal, PRECISION);
	mpf_init2(countimg, PRECISION);
	mpf_init2(ab, PRECISION);
	
	mpf_t menouno;
	mpf_init(menouno);
	mpf_set_si(menouno, -1);
	
	mpf_t den;
	mpf_init2(den, PRECISION);
	mpf_t binomial;
	mpf_init2(binomial, PRECISION);
	
	mpf_t fact_k;
	mpf_init2(fact_k, PRECISION);
	mpf_t fact_nk;
	mpf_init2(fact_nk, PRECISION);
	mpf_t fact_n;
	mpf_init2(fact_n, PRECISION);
	
	mpf_t k_count;
	mpf_init2(k_count, PRECISION);
	mpf_t nk_count;
	mpf_init2(nk_count, PRECISION);
	mpf_t n_count;
	mpf_init2(n_count, PRECISION);
	
	unsigned long k;
	for(k=0; k<=n; ++k){
	
		//Inizializzo n, k e n-k
		mpf_set_ui(k_count, k);
		mpf_set_ui(n_count, n);
		mpf_sub(nk_count, n_count, k_count);

		//Calcolo il prodotto tra i fattoriali di k e n-k
		mpf_factorial(fact_k, k_count);
		mpf_factorial(fact_nk, nk_count);
		mpf_mul(den, fact_k, fact_nk);
		//gmp_printf("fact(%d): %Ff\n",k, fact_k);
		//gmp_printf("fact(%d): %Ff\n",n-k, fact_nk);	
		//gmp_printf("den: %Ff\n", den);

		//Calcolo il rapporto tra i fattoriali di n e k*n-k
		mpf_factorial(fact_n, n_count);
		mpf_div(binomial, fact_n, den);
		//gmp_printf("fact(%d): %Ff\n",n, fact_n);	
		//gmp_printf("den: %Ff\n", den);
		//gmp_printf("binomial: %Ff\n", binomial);	

		//gmp_printf("BINOMIALE con K=%.*Ff: %.*Ff\n", 0, k_count, 0, binomial);
		
		//Elevo a potenza la parte reale e quella immaginaria	
		mpf_pow_ui(a, c.real, n-k);
		mpf_pow_ui(b, c.img, k);		
		/*gmp_printf("\n\nc.real: %.*Ff\n", 50, c.real);		
		gmp_printf("n-k: %d\n", n-k);
		gmp_printf("a: %.*Ff\n", 50, a);
		gmp_printf("c.img: %.*Ff\n", 50, c.img);		
		gmp_printf("k: %d\n", k);
		gmp_printf("b: %.*Ff\n\n", 50, b);*/

		//Controllo sul segno della parte immaginaria
		//Se k%4 == 2 o 3, allora la i alla k produce segno negativo 
		if(k%4 == 2 || k%4 == 3)
			mpf_mul(b, b, menouno);
		
		//Calcolo il prodotto tra a e b
		mpf_mul(ab, a, b);

		/*gmp_printf("\n\na: %.*Ff\n", 50, a);	
		gmp_printf("b: %.*Ff\n", 50, b);	
		gmp_printf("ab: %.*Ff\n\n", 50, ab);*/

		//gmp_printf("\n\nab: %.*Ff\n", 50, ab);
		mpf_mul(ab, ab, binomial);

			
		/*gmp_printf("binomial: %.*Ff\n", 50, binomial);	
		gmp_printf("ab: %.*Ff\n\n", 50, ab);*/


		//Controllo su cosa sia il risultato "ab", reale o immaginaria?
		if(k%4 == 1 || k%4 == 3)
		{
			//gmp_printf("\n\ncountimg :%.*Ff\n", 50,countimg);
			//gmp_printf("ab: %.*Ff\n", 50, ab);
			mpf_add(countimg, countimg, ab);
			//gmp_printf("countimg :%.*Ff\n\n", 50,countimg);
		}		
		else
		{
			//gmp_printf("\n\ncountreal :%.*Ff\n", 50,countreal);
			//gmp_printf("ab: %.*Ff\n", 50, ab);
			mpf_add(countreal, countreal, ab);
			//gmp_printf("countreal :%.*Ff\n\n", 50,countreal);
		}
	}
	
	//Assegno i valori accumulati alla struct complex_f di ritorno
	mpf_set(result->real, countreal);
	mpf_set(result->img, countimg); 
	
	mpf_clear(a);
	mpf_clear(b);
	mpf_clear(ab);
	mpf_clear(countreal);
	mpf_clear(countimg);	
	mpf_clear(menouno);
	mpf_clear(den);
	mpf_clear(binomial);
	mpf_clear(fact_k);
	mpf_clear(fact_nk);
	mpf_clear(fact_n);
	mpf_clear(k_count);
	mpf_clear(nk_count);
	mpf_clear(n_count);
	
	return ;
}

void deltaprecision(complex_f* result, complex_f c){

	//Inizializzazione variabili e stuct
	double n;
	int a = -1;
	unsigned long exp1;
	unsigned long exp2;	
	
	complex_f aexp;
	complex_finit(&aexp);
	complex_f res1;
	complex_finit(&res1);
	complex_f res2;	
	complex_finit(&res2);
	complex_f sum;
	complex_finit(&sum);
	complex_f temp_part_sum;
	complex_finit(&temp_part_sum);
	complex_f part_sum;
	complex_finit(&part_sum);
	complex_f temp_result;
	complex_finit(&temp_result);		
	
	//Inizio calcolo sommatoria
	for(n=1; n<=LAPNUMBER; ++n){
		
//		gmp_printf("\nITERAZIONE %.*f:\n", 0, n);
		
//		gmp_printf("PARTSUM:\n\tre-->%.*Ff\n\timg-->%.*Ff\n", 0, part_sum.real, 0, part_sum.img);
		mpf_set_d(aexp.real, pow(a, n));
		//gmp_printf("aexp.real: %Ff\n", aexp.real);
		//gmp_printf("aexp.img: %Ff\n", aexp.img);

		//Gli esponenti risultano sempre interi
		exp1 = ((n * ((3 * n) - 1)) / 2);
		exp2 = ((n * ((3 * n) + 1)) / 2);
//		printf("EXP1: %lu\n", exp1);
//		printf("EXP2: %lu\n", exp2);
		
		powbin(&res1, c, exp1);//////////////////////
		powbin(&res2, c, exp2);/////////////////////
		
		//pow_no_bin(&res1, c, exp1);
		//pow_no_bin(&res2, c, exp2);

		/*printf("exp1: %d\n", exp1);

		gmp_printf("res1.real: %.*Ff\n", 1000, res1.real);		
		gmp_printf("res1.img:  %.*Ff\n", 1000, res2.img);*/

		/*gmp_printf("\n\nc.real: %.*Ff\n", 50, c.real);
		gmp_printf("c.img: %.*Ff\n", 50, c.img);
		gmp_printf("exp1: %d\n", exp1);
		gmp_printf("res1.real: %.*Ff\n", 500, res1.real);
		gmp_printf("res1.img: %.*Ff\n\n", 500, res1.img);*/

//		gmp_printf("Q1:\n\tre-->%.*Ff\n\timg-->%.*Ff\n", 0, res1.real, 0, res1.img);
//		gmp_printf("Q2:\n\tre-->%.*Ff\n\timg-->%.*Ff\n", 0, res2.real, 0, res2.img);
		
		//Calcolo di sum += aexp * (q1 + q2);
		add_fcomplex(&temp_part_sum, res1, res2);
//		gmp_printf("Q1+Q2:\n\tre-->%.*Ff\n\timg-->%.*Ff\n", 0, part_sum.real, 0, part_sum.img);

		//gmp_printf("\n\naexp.real: %Ff  aexp.img: %Ff\n", aexp.real, aexp.img);
		//gmp_printf("part_sum.real: %Ff  part_sum.img: %Ff\n", part_sum.real, part_sum.img);
		mul_fcomplex(&part_sum, aexp, temp_part_sum);
		//gmp_printf("part_sum.real: %Ff  part_sum.img: %Ff\n\n", part_sum.real, part_sum.img);

		//gmp_printf("\n\nsum.real: %Ff  sum.img: %Ff\n", sum.real, sum.img);		
		//gmp_printf("part_sum.real: %Ff  part_sum.img: %Ff\n", part_sum.real, part_sum.img);
		add_fcomplex(&sum , sum, part_sum);
		//gmp_printf("sum.real: %Ff  sum.img: %Ff\n\n", sum.real, sum.img);
		
//		gmp_printf("SUM:\n\tre-->%.*Ff\n\timg-->%.*Ff\n", 0, sum.real, 0, sum.img);
	}
	
	//gmp_printf("SOMMA.real: %Ff\n", sum.real);
	//gmp_printf("SOMMA.img: %Ff\n", sum.img);

//	gmp_printf("TOT-SUM:\n\tre-->%.*Ff\n\timg-->%.*Ff\n", 0, sum.real, 0, sum.img);

	mpf_add_ui(sum.real, sum.real, 1);
	
	//gmp_printf("SOMMA.real: %Ff\n", sum.real);
	//gmp_printf("SOMMA.img: %Ff\n", sum.img);


//	gmp_printf("TOT-SUM+1:\n\tre-->%.*Ff\n\timg-->%.*Ff\n", 0, sum.real, 0, sum.img);
		
//	//Inizializzo variabili per la potenza 24-esima
//	complex_f tot;
//	complex_finit(&tot);
	
	//Inizio elevazione a potenza 24-esima
	//gmp_printf("\n\nsum.real: %Ff\nsum.img: %Ff\n", sum.real, sum.img);


	powbin(&temp_result, sum, 24);//////////////////////
	//pow_no_bin(&temp_result, sum, 24);

	//gmp_printf("result->real: %Ff\n\n", result->real);
	//gmp_printf("result->img: %Ff\n\n", result->img);


	//gmp_printf("result->real: %Ff\nresult->img: %Ff\n\n", result->real, result->img);
//	gmp_printf("TOT-SUM^24:\n\tre-->%.*Ff\n\timg-->%.*Ff\n", 0, result->real, 0, result->img);

	//Calcolo il prodotto con il numero complesso c
	mul_fcomplex(result, c, temp_result); //SBAGLIATA! :D

	//gmp_printf("result->real: %Ff\n\n", result->real);
	//gmp_printf("result->img: %Ff\n\n", result->img);


//	gmp_printf("RISULTATO:\nre-->%.*Ff\nimg-->%.*Ff\n", 0, tot.real, 0, tot.img);
	
	complex_fclear(&aexp);	
	complex_fclear(&res1);
	complex_fclear(&res2);
	complex_fclear(&sum);
	complex_fclear(&part_sum);
	complex_fclear(&temp_part_sum);
	complex_fclear(&temp_result);

	return;
}
/*
int main(){

	struct complex_f c, result;
	
	complex_finit(&c);
	complex_finit(&result);
	
	mpf_set_d(c.real, 2);
	mpf_set_d(c.img, 3);
	
//	powbin(c, 150, &result);
//	gmp_printf("POWBIN:\n%.*Ff\n%.*Ff\n", 1, result.real, 1, result.img);
	deltaprecision(cexp(c));

	mpf_t n;
	mpf_init2(n, PRECISION);
	mpf_t res;
	mpf_init2(res, PRECISION);
	mpf_set_ui(n, 150);
	mpf_factorial(res, n);
	
	gmp_printf("Fattoriale di %.*Ff: %.*Ff\n", 0, n, 1, res);
	
	//printf("Buon compleanno Mr. Pup!!!!\n");
	
	return 0;
}
*/
