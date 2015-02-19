/* hilbert.c */

#include "hilbert.h"

/*Dato un discriminante D(negativo), questa funzione ritorna
 *una combinazione di:
 * -numero di classe h(D)
 * -il polinomio di classe Hilbert, con grado h(D)
 * -l'insieme delle forme ridotte (a,b,c) del discriminante D
*/
void hilbertclass(Hresult* res, mpz_t D, mpz_t n){

	assert(mpz_cmp_ui(D, 0) < 0);

	mpz_t Dmod4;
	mpz_init(Dmod4);
	mpz_mod_ui(Dmod4, D, 4);

	assert(mpz_cmp_ui(Dmod4, 0) == 0 || mpz_cmp_ui(Dmod4, 1) == 0);
	mpz_clear(Dmod4);

	//1.[Initialize]
	res->count = 0;
	//Inizializzo il polinomio T = 1
	int i;	
	complex_f* coef_T = (complex_f*)malloc(sizeof(complex_f));	
	complex_finit(&coef_T[0]);	
	mpf_init_set_str(coef_T[0].real, "1", 10);
	mpf_init_set_str(coef_T[0].img, "0", 10);
	
	int* deg_T = (int*)malloc(sizeof(int));
	deg_T[0] = 0;
	
	int num_el_T = 1;
	f_polynomial T;
	init_complex_fpolynomial(coef_T, num_el_T, deg_T, num_el_T, &T);
	
	//complex_fpolynomial_printf(&T);
	
	//Inizializzo altre variabili
	mpz_t b;
	mpz_init(b);
	mpz_mod_ui(b, D, 2);			//b = D mod 2
	
	mpz_t r;
	mpz_init(r);					//r
	mpz_t Dabs;
	mpz_init(Dabs);
	mpz_abs(Dabs, D);				//Dabs = |D|
	mpf_t Dabsfloat;					
	mpf_init(Dabsfloat);			//Dasbf
	mpf_set_z(Dabsfloat, Dabs);		//Dabsfloat = Dabs
	mpf_t Ddiv3;
	mpf_init(Ddiv3);				//Ddiv3; 
	mpf_div_ui(Ddiv3, Dabsfloat, 3);//Dvid3 = |D|/3	
	mpf_t rfloat;
	mpf_init(rfloat);				//rfloat
	mpf_sqrt(rfloat, Ddiv3);		//rfloat = sqrt(|D|/3)
	mpf_floor(rfloat, rfloat);		//rfloat = floor(sqrt(|D|/3)
	mpz_set_f(r, rfloat);			//r = floor(sqrtDz);
	
	//gmp_printf("r: %Zd\n", r);

	mpz_t h;
	mpz_init_set_ui(h, 0);			//h
	
	mpz_t red[50];	//non credo vengano superate le 50 triple <a,b,c>
	int k;
	for(k=0; k < 50; ++k)
		mpz_init(red[k]);
	
	//2.[Outer loop on b]
	while(mpz_cmp(b, r) <= 0) //in realtà viene eseguito solo una iterazione
	{

		//m = (b^2 - D)/4
		mpz_t m;						
		mpz_init(m);				//float m;
	
		mpz_t b2;
		mpz_init(b2);
		mpz_pow_ui(b2, b, 2);		//b2 = pow(b, 2);
		
		mpz_t B2subD;
		mpz_init(B2subD);
		mpz_sub(B2subD, b2, D);		//B2subD = b^2 - D;
				
		mpz_div_ui(m, B2subD, 4);	//m = B2subD / 4;

		mpf_t afloat;
		mpf_init(afloat);			//afloat = a;
		
		mpf_t c;
		mpf_init(c);				//float c;
		
		mpz_t cint;
		mpz_init(cint);

		mpz_t aint;
		mpz_init(aint);				//aint;
		
		unsigned int a = 1;			//a
		mpz_t mmoda;	
		mpz_init(mmoda);			//mmoda;
		double complex tau;			//tau <complex.h>

		mpf_t mfloat;
		mpf_init(mfloat);
		mpf_set_z(mfloat, m);
		
		//a = max(b, 1)
		unsigned int max_b_1;

		if(mpz_cmp_ui(b, 1) > 0)
			max_b_1 = mpz_get_ui(b);
		else
			max_b_1 = 1;

//		printf("WHILEWHILE\n");
		for(a = max_b_1; mpf_cmp_d(mfloat, pow(a, 2)) >= 0; ++a)
		//for(a=1; pow(a, 2) <= mpf_get_d(mfloat); ++a)
		{	
//			printf("CMP: %d\n", mpf_cmp_d(mfloat, pow(a, 2)));		

			mpz_set_ui(aint, a);				//aint = a;
			mpz_mod(mmoda, m, aint);
			
			if(mpz_cmp_ui(mmoda, 0) != 0)
				continue;
			
			mpf_set_z(afloat, aint);			//afloat = aint;
			mpf_div(c, mfloat, afloat);		
			
			mpz_set_ui(cint, mpf_get_ui(c));

			if(mpz_cmp(b, aint) > 0)
				continue;
		
			//printf("A: %d\n", a);
			//gmp_printf("B: %Zd\n", b);
			//gmp_printf("C: %Ff\n", c);
			//gmp_printf("MF: %.*Ff\n", 10, mfloat);
	
			mpz_t gcd;
			mpz_init(gcd);
			
			mpz_gcd(gcd, aint, b);

			mpz_gcd(gcd, gcd, cint);

			if(mpz_cmp_ui(gcd, 1) != 0)
				continue;

			//3.[Optional polynomial setup]
			double Ddouble = mpf_get_d(Dabsfloat);
			double bdouble = mpz_get_d(b);
			double adouble = mpf_get_d(afloat);
			
			double complex num;					//num
			double complex den;					//den
			
			/*La stampa tronca il numero che ho inserito facendo degli arrotondamenti
			double complex prova = 3.56670959665405490549540505354054945048504580485408540584058405805;
			printf("PROVA: %.*f\n", 100, creal(prova));
			*/
//#############################################################
			
			//Creazione del denominatore 2a
			mpf_t afloat2;
			mpf_init(afloat2);
			mpf_mul_ui(afloat2, afloat, 2);
			complex_f comp_afloat2;
			complex_finit(&comp_afloat2);
			mpf_set(comp_afloat2.real, afloat2);
			mpf_set_ui(comp_afloat2.img, 0);

			//faccio diventare b un mpf_t
			mpf_t bgmp;
			mpf_init(bgmp);
			mpf_set_z(bgmp, b);
			//bmeno per il risultato
			complex_f bmeno;
			complex_finit(&bmeno);
			//-1 gmp
			mpf_t menoone;
			mpf_init_set_si(menoone, -1);
			
			mpf_mul(bmeno.real, bgmp, menoone); //bmeno = -b
			//mpf_set_ui(bmeno.img, 0); 
			
			//sqrtDabs la radice del valore assoluto di D
			complex_f sqrtDabs;
			complex_finit(&sqrtDabs);
		
			mpf_sqrt(sqrtDabs.img, Dabsfloat);
			
			//I complesso gmp
			complex_f comp_i;
			complex_finit(&comp_i);
			mpf_set_ui(comp_i.real, 0);
			mpf_set_ui(comp_i.img, 1);
			
			//mul_fcomplex(&sqrtDabs, sqrtDabs, comp_i);	//sqrt(|D|)*I
			
			complex_f tempnum;
			complex_finit(&tempnum);
			add_fcomplex(&tempnum, bmeno, sqrtDabs);		//sommo i termini al numeratore
			
			/*
			gmp_printf("TEMPNUM-GMP:\n\treal->%.*Ff\n\timg->%.*Ff\n", 100, tempnum.real, 100, tempnum.img);
			gmp_printf("AFLOAT-GMP:\n\treal->%.*Ff\n\timg->%.*Ff\n", 100, comp_afloat2.real, 100, comp_afloat2.img);
			*/
			
			//tau = risultato div
			complex_f comp_tau;
			complex_finit(&comp_tau);
			div_fcomplex(&comp_tau, tempnum, comp_afloat2);
			
			//pigreco gmp
			complex_f pigreco;
			complex_finit(&pigreco);
			mpf_set_str(pigreco.real, "3.1415926535897932384626433832795028841971693993751058209749445923078164", 10);
			//gmp_printf("PIGRECO-GMP:\n\treal->%.*Ff\n", 100, pigreco.real);
			
			//numgmp e dengmp
			complex_f numgmp;
			complex_f dengmp;
			complex_finit(&numgmp);
			complex_finit(&dengmp);
			
			//costruzione di 4*PI*I*tau e 2*PI*I*tau 
			mpf_set(numgmp.real, comp_tau.real);
			mpf_set(numgmp.img, comp_tau.img);
			mpf_set(dengmp.real, comp_tau.real);
			mpf_set(dengmp.img, comp_tau.img);
			
			mpf_mul_ui(numgmp.real, numgmp.real, 4);
			mpf_mul_ui(numgmp.img, numgmp.img, 4);
			mpf_mul_ui(dengmp.real, dengmp.real, 2);
			mpf_mul_ui(dengmp.img, dengmp.img, 2);
			
			mul_fcomplex(&numgmp, numgmp, pigreco);
			mul_fcomplex(&dengmp, dengmp, pigreco);
	
			mul_fcomplex(&numgmp, numgmp, comp_i);
			mul_fcomplex(&dengmp, dengmp, comp_i);
			
			complex_f expnumgmp;
			complex_finit(&expnumgmp);
			complex_f expdengmp;
			complex_finit(&expdengmp);			
			
			other_exponential(&expdengmp, dengmp); 
			mul_fcomplex(&expnumgmp, expdengmp, expdengmp);	


			//gmp_printf("TAU-GMP:\n\treal->%.*Ff\n\timg->%.*Ff\n", 100, comp_tau.real, 100, comp_tau.img);

			//gmp_printf("NUM-GMP:\n\treal->%.*Ff\n\timg->%.*Ff\n", 100, numgmp.real, 100, numgmp.img);
			//gmp_printf("DEN-GMP:\n\treal->%.*Ff\n\timg->%.*Ff\n", 100, dengmp.real, 100, dengmp.img);

			//gmp_printf("NUMEXP-GMP:\n\treal->%.*Ff\n\timg->%.*Ff\n", 100, expnumgmp.real, 100, expnumgmp.img);
			//gmp_printf("EXPDEN-GMP:\n\treal->%.*Ff\n\timg->%.*Ff\n", 100, expdengmp.real, 100, expdengmp.img);
//#############################################################

			tau = (-bdouble+sqrt(Ddouble)*I)/(2*adouble);
			
			num = tau*4*M_PI*I;
			den = tau*2*M_PI*I;
			
			//printf("PIGRECO:\n\treal->%.*f\n", 100, M_PI);
			//printf("TAU:\n\treal->%.*f\n\timag->%.*f\n", 100, creal(tau), 100, cimag(tau));
			//printf("NUM: \n\treal->%.*f\n\timag->%.*f\n", 100, creal(num), 100, cimag(num));
			//printf("DEN: \n\treal->%.*f\n\timag->%.*f\n", 100, creal(den), 100, cimag(den));
			
			complex_f resnum;
			complex_finit(&resnum);
			complex_f resden;
			complex_finit(&resden);
			complex_f expnum;
			complex_finit(&expnum);
			complex_f expden;
			complex_finit(&expden);			

			mpf_set_d(expnum.real, creal(cexp(num)));
			mpf_set_d(expnum.img, cimag(cexp(num)));
			mpf_set_d(expden.real, creal(cexp(den)));
			mpf_set_d(expden.img, cimag(cexp(den)));
			
			//printf("NUMEXP:\n\treal->%.*f\n\timag->%.*f\n", 100, creal(cexp(num)), 100, cimag(cexp(num))); ////////////
			/*
			gmp_printf("mpf_creal: %.*Ff\n", 1000, expnum.real);			
			gmp_printf("mpf_img: %.*Ff\n", 1000, expnum.img);
			*/

			deltaprecision(&resnum, expnumgmp);/////////////////////
			deltaprecision(&resden, expdengmp);////////////////////
			
			//deltaprecision(&resnum, expnum);
			//deltaprecision(&resden, expden);
			
			
			//gmp_printf("expnum.real: %.*Ff\n", 50, expnum.real);			
			//gmp_printf("expnum.img: %.*Ff\n", 50, expnum.img);

			//gmp_printf("expden.real: %.*Ff\n", 50, expden.real);			
			//gmp_printf("expden.img: %.*Ff\n", 50, expden.img);

			//gmp_printf("num: %.*Ff +j%.*Ff\n", 50, resnum.real, 50, resnum.img);
			//gmp_printf("den: %.*Ff +j%.*Ff\n", 50, resden.real, 50, resden.img);
			
			
			//Calcolo rapporto tra i risultati delle deltaprecision()
			complex_f f;
			complex_finit(&f);
			div_fcomplex(&f, resnum, resden);
			//gmp_printf("\nRapporto deltaprecision:\nF: %.*Ff +j%.*Ff\n", 50, f.real, 50, f.img);
			
			complex_f temp_j;
			complex_finit(&temp_j);
			add_fcomplex(&temp_j, temp_j, f);
			mpf_mul_ui(temp_j.real, temp_j.real, 256);
			mpf_mul_ui(temp_j.img, temp_j.img, 256);
			mpf_add_ui(temp_j.real, temp_j.real, 1);

			powbin(&temp_j, temp_j, 3);
			
			complex_f j;
			complex_finit(&j);
			div_fcomplex(&j, temp_j, f);
				
			//Inizializzo il polinomio X - j
			complex_f menoj;
			complex_finit(&menoj);
			mpf_t menouno;
			mpf_init_set_si(menouno, -1);
			mpf_mul(menoj.real, j.real, menouno);
			mpf_mul(menoj.img, j.img, menouno);
			
			int num_el_Xmenoj = 2;
			complex_f* coef_Xmenoj = (complex_f*)malloc(sizeof(complex_f)*num_el_Xmenoj);
			int i;
			for(i = 0; i < num_el_Xmenoj; i++)	
				complex_finit(&coef_Xmenoj[i]);	
	
			mpf_init_set_str(coef_Xmenoj[0].real, "1", 10);
			mpf_init_set_str(coef_Xmenoj[0].img, "0", 10);
			mpf_init_set(coef_Xmenoj[1].real, menoj.real);
			mpf_init_set(coef_Xmenoj[1].img, menoj.img);
	
			int* deg_Xmenoj = (int*)malloc(sizeof(int)*num_el_Xmenoj);
			deg_Xmenoj[0] = 1;
			deg_Xmenoj[1] = 0;

			f_polynomial Xmenoj;

			init_complex_fpolynomial(coef_Xmenoj, num_el_Xmenoj, deg_Xmenoj, num_el_Xmenoj, &Xmenoj);

			//Clear
			for(i = 0; i < num_el_Xmenoj; i++)	
				complex_fclear(&coef_Xmenoj[i]);
			free(coef_Xmenoj);
			free(deg_Xmenoj);
			
			//Inizializzo il polinomio di 2 grado alternativo
			//Creo -2Re(j)
			complex_f menoj2;
			complex_finit(&menoj2);
			mpf_mul_ui(menoj2.real, menoj.real, 2);
			//Creo |j|
			mpf_t jabs2;
			mpf_init(jabs2);
			mpf_t powr;
			mpf_init(powr);
			mpf_t powi;
			mpf_init(powi);
			mpf_pow_ui(powr, j.real, 2);
			mpf_pow_ui(powi, j.img, 2);
			mpf_add(jabs2, powr, powi);
			
			int num_el_alt = 3;
			complex_f* coef_altPoly = (complex_f*)malloc(sizeof(complex_f)*num_el_alt);
			int l;
			for(l = 0; l < num_el_alt; l++)	
				complex_finit(&coef_altPoly[i]);	
	
			mpf_init_set_str(coef_altPoly[0].real, "1", 10);
			mpf_init_set_str(coef_altPoly[0].img, "0", 10);
			mpf_init_set(coef_altPoly[1].real, menoj2.real);
			mpf_init_set_str(coef_altPoly[1].img, "0", 10);
			mpf_init_set(coef_altPoly[2].real, jabs2);
			mpf_init_set_str(coef_altPoly[2].img, "0", 10);
	
			int* deg_altPoly = (int*)malloc( sizeof(int)*num_el_alt);
			deg_altPoly[0] = 2;
			deg_altPoly[1] = 1;
			deg_altPoly[2] = 0;
			
			f_polynomial altPoly;
			init_complex_fpolynomial(coef_altPoly, num_el_alt, deg_altPoly, num_el_alt, &altPoly);
			
			//Clear
			for(i = 0; i < num_el_alt; i++)	
				complex_fclear(&coef_altPoly[i]);
			
			free(coef_altPoly);	
			free(deg_altPoly);
			
			//Polinomio float
			f_polynomial temp_poly;
			
			//4.[Begin divisor test]
			struct abc triple;
			mpz_init_set(triple.a, aint);
			mpz_init_set(triple.b, b);
			mpz_init_set(triple.c, cint);

			if(mpz_cmp(b, aint) == 0 || mpf_cmp(c, afloat) == 0 || mpz_cmp_ui(b, 0) == 0)
			{
				//printf("------------------IF------------------\n");
				//Prodotto tra polinomi				
				complex_fpolynomial_mul(&temp_poly, T, Xmenoj);
				complex_fpolynomial_clear(&T);
				complex_fpolynomial_copy(&T, temp_poly);
				
				//Incremento class number
				mpz_add_ui(h, h, 1);						//h = h+1
				
				//Inizializzo la tripla con a, b,c
				mpz_init_set(res->abc_set[res->count].a, triple.a);
				mpz_init_set(res->abc_set[res->count].b, triple.b);
				mpz_init_set(res->abc_set[res->count].c, triple.c);//c'è un problema qui che non riesco a capire
				res->count = res->count + 1;
			}
			else
			{
				//printf("------------------ELSE------------------\n");	
				//Prodotto tra polinomi
				complex_fpolynomial_mul(&temp_poly, T, altPoly);
				complex_fpolynomial_clear(&T);
				complex_fpolynomial_copy(&T, temp_poly);
				
				//Incremento class number
				mpz_add_ui(h, h, 2);						//h = h+2
				
				//Creo -b
				mpz_t menob;
				mpz_init_set(menob, b);
				mpz_mul_si(menob, menob, -1);
				
				//Inizializzo la tripla con a,-b,c
				struct abc triple_minusb;
				mpz_init_set(triple_minusb.a, aint);
				mpz_init_set(triple_minusb.b, menob);
				mpz_init_set(triple_minusb.c, cint);
			
				mpz_init_set(res->abc_set[res->count].a, triple.a);
				mpz_init_set(res->abc_set[res->count].b, triple.b);
				mpz_init_set(res->abc_set[res->count].c, triple.c);
				res->count = res->count + 1; 
				mpz_init_set(res->abc_set[res->count].a, triple_minusb.a);
				mpz_init_set(res->abc_set[res->count].b, triple_minusb.b);
				mpz_init_set(res->abc_set[res->count].c, triple_minusb.c);
				res->count = res->count + 1; 		
				
				//Clear
				mpz_clear(menob);
				mpz_clear(triple_minusb.a);
				mpz_clear(triple_minusb.b);
				mpz_clear(triple_minusb.c);
			}
			/*printf("---------------------\n");
			printf("intermedio\n");			
			complex_fpolynomial_printf(&T);
			printf("---------------------\n");*/
			
			//Clear

			mpz_clear(gcd);
			mpf_clear(afloat2);
			complex_fclear(&comp_afloat2);
			mpf_clear(bgmp);
			complex_fclear(&bmeno);
			mpf_clear(menoone);
			complex_fclear(&sqrtDabs);
			complex_fclear(&comp_i);
			complex_fclear(&tempnum);
			complex_fclear(&comp_tau);
			complex_fclear(&pigreco);
			complex_fclear(&numgmp);
			complex_fclear(&dengmp);
			complex_fclear(&expnumgmp);
			complex_fclear(&expdengmp);	


			complex_fclear(&resnum);
			complex_fclear(&resden);
			complex_fclear(&expnum);
			complex_fclear(&expden);
			complex_fclear(&f);
			complex_fclear(&temp_j);
			complex_fclear(&j);
			complex_fclear(&menoj);
			mpf_clear(menouno);
			complex_fpolynomial_clear(&Xmenoj);
			complex_fclear(&menoj2);
			mpf_clear(jabs2);
			mpf_clear(powr);
			mpf_clear(powi);
			complex_fpolynomial_clear(&altPoly);
			complex_fpolynomial_clear(&temp_poly);
			mpz_clear(triple.a);
			mpz_clear(triple.b);
			mpz_clear(triple.c);
		}
		
		//Clear
		mpz_clear(m);				//float m;
		mpz_clear(b2);
		mpz_clear(B2subD);
		mpf_clear(afloat);			//afloat = a;
		mpf_clear(c);				//float c;
		mpz_clear(cint);
		mpz_clear(aint);			//aint;
		mpz_clear(mmoda);			//mmoda;
		mpf_clear(mfloat);
		
		mpz_add_ui(b, b, 2);
		//break;
	}
	
	//5.[Return values of interest]
	mpz_init_set(res->number_class, h);	//Variabile di ritorno del class number per D
	
	//Inizializzo il polinomio di tipo real_big
	mpz_t* coef = (mpz_t*)malloc(sizeof(mpz_t)*(T.num_elements));

	mpf_t roundcoef;
	mpf_init2(roundcoef, PRECISION); //arrotondo i coefficienti all'intero più vicino
	mpf_t one_half; //0.5
	mpf_init2(one_half, PRECISION);
	mpf_init_set_str(one_half, "0.5", 10);

	/*printf("\n---------------------\n");
	printf("senza arrotondare\n");
	complex_fpolynomial_printf(&T);
	printf("---------------------\n");*/
	
	for(i = 0; i < T.num_elements; i++)
	{
		mpf_set(roundcoef, T.coefficients[T.num_elements - 1 - i].real);
		//gmp_printf("roundcoef: %Ff\n", roundcoef);
		if(mpf_cmp_ui(roundcoef, 0) > 0) //se uguale a zero lo lascio invariato
		{
			mpf_add(roundcoef, roundcoef, one_half);
			mpf_floor(roundcoef, roundcoef);
		}
		else if(mpf_cmp_ui(roundcoef, 0) < 0) 
		{
			mpf_sub(roundcoef, roundcoef, one_half);
			mpf_ceil(roundcoef, roundcoef);
		}
		//gmp_printf("roundcoef: %Ff\n", roundcoef);
		mpz_init(coef[i]);
		mpz_set_f(coef[i], roundcoef);	
	}
	
	//Inizializzazione polinomio REAL-BIG
	mpz_t* deg = (mpz_t*)malloc(sizeof(mpz_t)*(T.num_elements));
	
	for(i = 0; i < T.num_elements; i++)	//I r_big_polynomial hanno il grado più basso alla prima posizione
		mpz_init_set_ui(deg[i], i);

	/////////////// da commentare 	
	init_real_big_polynomial_no_mod(coef, T.num_elements, deg, T.num_elements, &(res->Hpoly));
		
	
	printf("Polinomio di Hilbert Originale\n");
	r_big_polynomial_printf(&(res->Hpoly));
	
	r_big_polynomial_clear(&(res->Hpoly));
	/////////////// da commentare

	init_real_big_polynomial(coef, T.num_elements, deg, T.num_elements, &(res->Hpoly), n);
	//printf("\n---------------------\n");
	//gmp_printf("Polinomio ridotto modulo \"%Zd\"\n", n);
	//r_big_polynomial_printf(&(res->Hpoly));
	//printf("---------------------\n");


	//Clear
	complex_fclear(&coef_T[0]);

	free(coef_T);
	free(deg_T);

	complex_fpolynomial_clear(&T);
	
	for(i = 0; i < (res->Hpoly).num_elements; i++)	
		mpz_clear(coef[i]);

	free(coef);
	free(deg);
	
	mpz_clear(b);
	mpz_clear(r);
	mpz_clear(Dabs);
	mpf_clear(Dabsfloat);
	mpf_clear(Ddiv3);
	mpf_clear(rfloat);
	mpz_clear(h);
	for(k=0; k < 50; ++k)
		mpz_clear(red[k]);
	mpf_clear(roundcoef);
	mpf_clear(one_half);

	return;
}

void exponential(complex_f* result, complex_f z)
{
	int i;
	complex_f powres;
	complex_finit(&powres);
	complex_f temp;
	complex_finit(&temp);
	complex_f fact;
	complex_finit(&fact);
	complex_f sum;
	complex_finit(&sum);
	mpf_t igmp;
	mpf_init2(igmp, PRECISION);
	
	for(i=0; i<CYCLES; ++i)
	{	
		mpf_set_ui(igmp, i);
		mpf_factorial(fact.real, igmp);
		//gmp_printf("fact(%d)->%.*Ff\n", i, 5, fact.real);
		
		powbin(&powres, z, i);
		div_fcomplex(&temp, powres, fact);
		add_fcomplex(result, *result, temp);
	}

	complex_fclear(&powres);
	complex_fclear(&temp);
	complex_fclear(&fact);
	complex_fclear(&sum);
	mpf_clear(igmp);
}

void other_exponential(complex_f* result, complex_f z)
{
	mpf_t u[2];
	mpf_init(u[0]);
	mpf_init(u[1]);
	mpf_set(u[0], z.real);
	mpf_set(u[1], z.img);

	mpf_t w[2];
	mpf_init(w[0]);
	mpf_init(w[1]);

	complex_exp(u, w);

	mpf_set(result->real, w[0]);
	mpf_set(result->img, w[1]);

	mpf_clear(u[0]);
	mpf_clear(u[1]);
	mpf_clear(w[0]);
	mpf_clear(w[1]);
}



