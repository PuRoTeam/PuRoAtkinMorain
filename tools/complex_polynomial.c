#include <stdio.h>
#include <stdlib.h>
#include "complex_polynomial.h"

void init_complex_polynomial(complex_z* coef, int c_size, int* deg, int d_size, c_polynomial* poly)
{
	//CONTROLLI DI CONSISTENZA
	assert(c_size == d_size && d_size > 0);
	int i;	
	
	for(i = 0; i < d_size; i++)
	{
		assert(deg[i] >= 0); //no frazioni! 1/x
		if(i != 0)
			assert(deg[i] < deg[i-1]); //strettamente minore
	}
	//CONTROLLI DI CONSISTENZA
 	
	int index = -1; //indice primo coefficiente non nullo
	
	for(i = 0; i < c_size; i++)
	{		
		if(mpz_cmp_ui(coef[i].real, 0) != 0 || mpz_cmp_ui(coef[i].img, 0) != 0) //coef[i] != 0
		{			
			index = i;
			break;
		}
	}
	
	if(index == -1)//polinomio si è annullato completamente, crea polinomio nullo
	{
		poly->coefficients = (complex_z*)malloc(sizeof(complex_z));	//un solo coefficiente
		complex_zinit(&(poly->coefficients[0]));		
		poly->num_elements = 1;
		poly->degree = 0;
		mpz_set_ui(poly->coefficients[0].real, 0);	
		mpz_set_ui(poly->coefficients[0].img, 0);	
		return;
	}

	if(index > 0)//almeno il primo coefficiente è diventato nullo
	{		
		int new_size = c_size - index;
		complex_z* new_coef = (complex_z*)malloc(sizeof(complex_z)*new_size);		
	
		int k = 0;
		for(i = index; i < c_size; i++)
		{
			mpz_init_set(new_coef[k].real, coef[i].real);
			mpz_init_set(new_coef[k].img, coef[i].img);
			k++;
		}

		k = 0;
		int* new_deg = (int*)malloc(sizeof(int)*new_size);
		for(i = index; i < d_size; i++)
		{
			new_deg[k] = deg[i];
			k++;
		}
		init_complex_polynomial(new_coef, new_size, new_deg, new_size, poly);
		return;
	}

	int poly_size = deg[0] + 1; //deg(x+3)=1 -- > poly_size = 2, perchè contiene termine di grado 1 e termine di grado 0
	
	poly->coefficients = (complex_z*)malloc(sizeof(complex_z)*poly_size);

	for(i = 0; i < poly_size; i++) 
		complex_zinit(&(poly->coefficients[i])); //richiama sulla parte reale e immaginaria le funzioni mpz_t


	poly->num_elements = poly_size;
	poly->degree = deg[0]; //poly_size - 1

	int k = 0; //elemento di deg correntemente scansionato	

	for(i = 0; i < poly_size; i++)
	{
		if(k < d_size)	
		{
			if(i == deg[0] - deg[k])
			{
				//mpz_mod(coef[k], coef[k], n);
				mpz_set(poly->coefficients[i].real, coef[k].real);
				mpz_set(poly->coefficients[i].img, coef[k].img);
				k++;
			}
			else //manca termine intermedio
			{
				mpz_set_ui(poly->coefficients[i].real, 0);
				mpz_set_ui(poly->coefficients[i].img, 0);
			}
		}
		else //manca termine intermedio
		{
			mpz_set_ui(poly->coefficients[i].real, 0);
			mpz_set_ui(poly->coefficients[i].img, 0);		
		}
	}	
}

void complex_polynomial_mul(c_polynomial* poly_mul, c_polynomial poly_A, c_polynomial poly_B)
{
	int mul_deg = (poly_A.degree) + (poly_B.degree); //grado polinomio moltiplicazione
	int num_elements = mul_deg + 1;	

	int i;
	int j;
	
	complex_z* coefficients = (complex_z*)malloc(sizeof(complex_z)*num_elements);

	for(i = 0; i < num_elements; i++)
		complex_zinit(&(coefficients[i]));

	int* degrees = (int*)malloc(sizeof(int)*num_elements);	
	
	for(i = 0; i < num_elements; i++)	
		degrees[i] = mul_deg - i;
	
	//for(i = 0; i < num_elements; i++)
	//	complex_zinit(poly_mul->coefficients[i]);
	
	c_polynomial* max_poly;
	c_polynomial* min_poly;
	
	//int max_degree;
	int max_num_elements;
	int min_num_elements;

	if(poly_A.degree >= poly_B.degree)
	{
		max_poly = &poly_A;
		min_poly = &poly_B;
		//max_degree = poly_A->degree;
		max_num_elements = poly_A.num_elements;
		min_num_elements = poly_B.num_elements;
	}
	else
	{
		max_poly = &poly_B;
		min_poly = &poly_A;
		//max_degree = poly_B->degree;
		max_num_elements = poly_B.num_elements;
		min_num_elements = poly_A.num_elements;
	}


	for(i = 0; i < max_num_elements; i++)
	{
		for(j = 0; j < min_num_elements; j++)
		{
			int max_cur_deg = max_poly->degree - i;
			int min_cur_deg = min_poly->degree - j;			
			int deg = max_cur_deg + min_cur_deg;
			int pos = mul_deg - deg;
			//printf(/*"max_num: %d min_num: %d i: %d j: %d */"max: %d min: %d deg: %d\n",/* max_num_elements, min_num_elements, i, j,*/ max_cur_deg, min_cur_deg, deg);

			complex_z mul;
			complex_zinit(&mul);
			mul_zcomplex(&mul, max_poly->coefficients[i], min_poly->coefficients[j]);
			//gmp_printf("%Zd  %Zd\n", mul.real, mul.img);

			add_zcomplex(&coefficients[pos], coefficients[pos], mul);
		}	
	}

	//for(i = 0; i < num_elements; i++)
	//	gmp_printf("%Zd  %Zd\n", coefficients[i].real, coefficients[i].img);


	init_complex_polynomial(coefficients, num_elements, degrees, num_elements, poly_mul); //trovati i coefficienti, inizializzo il polinomio

}

void complex_polynomial_real_part(c_polynomial* real_poly, c_polynomial poly)
{
	int i;
	int num_elements = poly.num_elements;

	real_poly->coefficients = (complex_z*)malloc(sizeof(complex_z)*num_elements);		
	for(i = 0; i < num_elements; i++)
		complex_zinit(&(real_poly->coefficients[i]));

	real_poly->degree = num_elements - 1;
	real_poly->num_elements = num_elements;	

	
	for(i = 0; i < num_elements; i++)
	{
		mpz_set(real_poly->coefficients[i].real, poly.coefficients[i].real);
		mpz_set_ui(real_poly->coefficients[i].img, 0); 
	}
}

void init_complex_fpolynomial(complex_f* coef, int c_size, int* deg, int d_size, f_polynomial* poly)
{
	//CONTROLLI DI CONSISTENZA
	assert(c_size == d_size && d_size > 0);
	int i;	
	
	for(i = 0; i < d_size; i++)
	{
		assert(deg[i] >= 0); //no frazioni! 1/x
		if(i != 0)
			assert(deg[i] < deg[i-1]); //strettamente minore
	}
	//CONTROLLI DI CONSISTENZA
 	
	int index = -1; //indice primo coefficiente non nullo
	
	for(i = 0; i < c_size; i++)
	{		
		if(mpf_cmp_ui(coef[i].real, 0) != 0 || mpf_cmp_ui(coef[i].img, 0) != 0) //coef[i] != 0
		{			
			index = i;
			break;
		}
	}
	
	
	if(index == -1)//polinomio si è annullato completamente, crea polinomio nullo
	{
		poly->coefficients = (complex_f*)malloc(sizeof(complex_f));	//un solo coefficiente
		complex_finit(&(poly->coefficients[0]));		
		poly->num_elements = 1;
		poly->degree = 0;
		mpf_set_ui(poly->coefficients[0].real, 0);	
		mpf_set_ui(poly->coefficients[0].img, 0);	
		return;
	}

	if(index > 0)//almeno il primo coefficiente è diventato nullo
	{		
		//printf("blablablabla\n");
		int new_size = c_size - index;
		complex_f* new_coef = (complex_f*)malloc(sizeof(complex_f)*new_size);		
	
		int k = 0;
		for(i = index; i < c_size; i++)
		{
			mpf_init2(new_coef[k].real, PRECISION);
			mpf_set(new_coef[k].real, coef[i].real);
			mpf_init2(new_coef[k].img, PRECISION);
			mpf_set(new_coef[k].img, coef[i].img);
			k++;
		}

		k = 0;
		int* new_deg = (int*)malloc(sizeof(int)*new_size);
		for(i = index; i < d_size; i++)
		{
			new_deg[k] = deg[i];
			k++;
		}
		init_complex_fpolynomial(new_coef, new_size, new_deg, new_size, poly);
		
		//Clear
		for(i = 0; i < new_size; i++)
		{
			mpf_clear(new_coef[i].real);
			mpf_clear(new_coef[i].img);
		}
		free(new_coef);
		free(new_deg);
		return;
	}

	int poly_size = deg[0] + 1; //deg(x+3)=1 -- > poly_size = 2, perchè contiene termine di grado 1 e termine di grado 0
	
	poly->coefficients = (complex_f*)malloc(sizeof(complex_f)*poly_size);

	for(i = 0; i < poly_size; i++) 
		complex_finit(&(poly->coefficients[i])); //richiama sulla parte reale e immaginaria le funzioni mpz_t


	poly->num_elements = poly_size;
	poly->degree = deg[0]; //poly_size - 1

	int k = 0; //elemento di deg correntemente scansionato	

	for(i = 0; i < poly_size; i++)
	{
		if(k < d_size)	
		{
			if(i == deg[0] - deg[k])
			{
				//mpz_mod(coef[k], coef[k], n);
				mpf_set(poly->coefficients[i].real, coef[k].real);
				mpf_set(poly->coefficients[i].img, coef[k].img);
				k++;
			}
			else //manca termine intermedio
			{
				mpf_set_ui(poly->coefficients[i].real, 0);
				mpf_set_ui(poly->coefficients[i].img, 0);
			}
		}
		else //manca termine intermedio
		{
			mpf_set_ui(poly->coefficients[i].real, 0);
			mpf_set_ui(poly->coefficients[i].img, 0);		
		}
	}	
}

void complex_fpolynomial_mul(f_polynomial* poly_mul, f_polynomial poly_A, f_polynomial poly_B)
{
	int mul_deg = (poly_A.degree) + (poly_B.degree); //grado polinomio moltiplicazione
	int num_elements = mul_deg + 1;	

	int i;
	int j;
	
	complex_f* coefficients = (complex_f*)malloc(sizeof(complex_f)*num_elements);

	for(i = 0; i < num_elements; i++)
		complex_finit(&(coefficients[i]));

	int* degrees = (int*)malloc(sizeof(int)*num_elements);	
	
	for(i = 0; i < num_elements; i++)	
		degrees[i] = mul_deg - i;
	
	//for(i = 0; i < num_elements; i++)
	//	complex_zinit(poly_mul->coefficients[i]);
	
	f_polynomial* max_poly;
	f_polynomial* min_poly;
	
	//int max_degree;
	int max_num_elements;
	int min_num_elements;

	if(poly_A.degree >= poly_B.degree)
	{
		max_poly = &poly_A;
		min_poly = &poly_B;
		//max_degree = poly_A->degree;
		max_num_elements = poly_A.num_elements;
		min_num_elements = poly_B.num_elements;
	}
	else
	{
		max_poly = &poly_B;
		min_poly = &poly_A;
		//max_degree = poly_B->degree;
		max_num_elements = poly_B.num_elements;
		min_num_elements = poly_A.num_elements;
	}
	
	complex_f mul;
	complex_finit(&mul);

	for(i = 0; i < max_num_elements; i++)
	{
		for(j = 0; j < min_num_elements; j++)
		{
			int max_cur_deg = max_poly->degree - i;
			int min_cur_deg = min_poly->degree - j;			
			int deg = max_cur_deg + min_cur_deg;
			int pos = mul_deg - deg;
			//printf(/*"max_num: %d min_num: %d i: %d j: %d */"max: %d min: %d deg: %d\n",/* max_num_elements, min_num_elements, i, j,*/ max_cur_deg, min_cur_deg, deg);

			mul_fcomplex(&mul, max_poly->coefficients[i], min_poly->coefficients[j]);
			//complex_f mul = mul_fcomplex(max_poly->coefficients[i], min_poly->coefficients[j]);
			//gmp_printf("%Zd  %Zd\n", mul.real, mul.img);
			add_fcomplex(&coefficients[pos], coefficients[pos], mul);
			//coefficients[pos] = add_fcomplex(coefficients[pos], mul);
			
		}	
		
	}
	
	complex_fclear(&mul);
	//for(i = 0; i < num_elements; i++)
	//	gmp_printf("%Zd  %Zd\n", coefficients[i].real, coefficients[i].img);

	init_complex_fpolynomial(coefficients, num_elements, degrees, num_elements, poly_mul); //trovati i coefficienti, inizializzo il polinomio
	for(i = 0; i < num_elements; i++)
		complex_fclear(&(coefficients[i]));
	free(degrees);
	free(coefficients);
}

void complex_fpolynomial_real_part(f_polynomial* real_poly, f_polynomial poly)
{
	int i;
	int num_elements = poly.num_elements;

	real_poly->coefficients = (complex_f*)malloc(sizeof(complex_f)*num_elements);		
	for(i = 0; i < num_elements; i++)
		complex_finit(&(real_poly->coefficients[i]));

	real_poly->degree = num_elements - 1;
	real_poly->num_elements = num_elements;	
	
	for(i = 0; i < num_elements; i++)
	{
		mpf_set(real_poly->coefficients[i].real, poly.coefficients[i].real);
		mpf_set_ui(real_poly->coefficients[i].img, 0); 
	}
}

void complex_fpolynomial_copy(f_polynomial* T, f_polynomial temp_poly){
	
	T->num_elements = temp_poly.num_elements;
	
	T->coefficients = (complex_f*)malloc(sizeof(complex_f)*T->num_elements);
	T->degree = temp_poly.degree;

	int i;
	for(i = 0; i < temp_poly.num_elements; i++)
	{		
		mpf_init_set(T->coefficients[i].real, temp_poly.coefficients[i].real);
		mpf_init_set(T->coefficients[i].img, temp_poly.coefficients[i].img);
	}		
	
	return;
}

void complex_fpolynomial_clear(f_polynomial* poly)
{
	int i;
	for(i = 0; i < poly->num_elements; i++)
	{
		complex_fclear(&(poly->coefficients[i]));
	}
	free(poly->coefficients);
	poly->degree  = 0;
	poly->num_elements = 0;
}


void complex_fpolynomial_printf(f_polynomial* poly)
{
	int i = 0;
	for(i = 0; i < poly->num_elements; i++)
	{	
		gmp_printf("(%Ff + j%Ff)*X^%d \n", poly->coefficients[i].real, poly->coefficients[i].img, poly->degree - i);
	}
}