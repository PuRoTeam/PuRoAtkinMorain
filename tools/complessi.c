#include "complessi.h"

//Init struct complex_z
void complex_zinit(struct complex_z* c){
	
	mpz_init(c->real);
	mpz_set_ui(c->real, 0);
	mpz_init(c->img);
	mpz_set_ui(c->img, 0);
	return; 
}

//Clear struct complex_z
void complex_zclear(struct complex_z* c){
	mpz_clear(c->real);
	mpz_clear(c->img);
}

//Init struct complex_f
void complex_finit(struct complex_f* c){
	
	mpf_init2(c->real, PRECISION);
	mpf_set_ui(c->real, 0);
	mpf_init2(c->img, PRECISION);
	mpf_set_ui(c->img, 0);
	return; 
}

//Clear struct complex_f
void complex_fclear(struct complex_f* c){
	mpf_clear(c->real);
	mpf_clear(c->img);
}

//Add two complex and return a struct complex
void add_zcomplex(struct complex_z* z3, struct complex_z z1, struct complex_z z2 ){ 
	
	mpz_t zr;
	mpz_init(zr);
	mpz_t zi;
	mpz_init(zi);
	
	mpz_add(zr, z1.real, z2.real);
	mpz_set(z3->real, zr);
	mpz_add(zi, z1.img, z2.img);
	mpz_set(z3->img, zi);
	
	mpz_clear(zr);
	mpz_clear(zi);
	
	return; 
}

//Add two complex and return a struct complex
void add_fcomplex(struct complex_f* result, struct complex_f z1, struct complex_f z2 ){ 
	
	mpf_add(result->real, z1.real, z2.real);
	mpf_add(result->img, z1.img, z2.img);
	
	return; 
}

//Mul two complex and return a struct complex
//ora è possibile usare la stessa variabile per operando e risultato
void mul_zcomplex(struct complex_z* z3, struct complex_z z1, struct complex_z z2){ 

	struct complex_z temp_result;
	complex_zinit(&temp_result);
	
	mpz_t rr;
	mpz_init(rr);
	mpz_t ii;
	mpz_init(ii);
	mpz_t ri;
	mpz_init(ri);
	mpz_t ir;
	mpz_init(ir);

	mpz_mul(rr, z1.real, z2.real);
	mpz_mul(ii, z1.img, z2.img);
	mpz_sub(temp_result.real, rr, ii);
	
	mpz_mul(ri, z1.real, z2.img);
	mpz_mul(ir, z1.img, z2.real);
	mpz_add(temp_result.img, ri, ir);

	mpz_set(z3->real, temp_result.real);
	mpz_set(z3->img, temp_result.img);

	complex_zclear(&temp_result);
	
	mpz_clear(rr);
	mpz_clear(ii);
	mpz_clear(ri);
	mpz_clear(ir);
	
	return; 
}

//ora è possibile usare la stessa variabile per operando e risultato
void mul_fcomplex(struct complex_f* result, struct complex_f z1, struct complex_f z2){ 

	struct complex_f temp_result;
	complex_finit(&temp_result);	

	mpf_t rr;
	mpf_init2(rr, PRECISION);
	mpf_t ii;
	mpf_init2(ii, PRECISION);
	mpf_t ri;
	mpf_init2(ri, PRECISION);
	mpf_t ir;
	mpf_init2(ir, PRECISION);

	mpf_mul(rr, z1.real, z2.real);
	mpf_mul(ii, z1.img, z2.img);
	mpf_sub(temp_result.real, rr, ii);
	
	mpf_mul(ri, z1.real, z2.img);
	mpf_mul(ir, z1.img, z2.real);
	mpf_add(temp_result.img, ri, ir);
	
	mpf_set(result->real, temp_result.real);
	mpf_set(result->img, temp_result.img);

	complex_fclear(&temp_result);
	mpf_clear(rr);
	mpf_clear(ii);
	mpf_clear(ri);
	mpf_clear(ir);

	return; 
}

//ora è possibile usare la stessa variabile per operando e risultato
void div_fcomplex(struct complex_f* result, struct complex_f z1, struct complex_f z2){ 
	
	struct complex_f temp_result;
	complex_finit(&temp_result);	

	mpf_t rr;
	mpf_init2(rr, PRECISION);
	mpf_t ii;
	mpf_init2(ii, PRECISION);
	mpf_t ri;
	mpf_init2(ri, PRECISION);
	mpf_t ir;
	mpf_init2(ir, PRECISION);
	mpf_t powr;
	mpf_init2(powr, PRECISION);
	mpf_t powi;
	mpf_init2(powi, PRECISION);
	mpf_t z2abs;
	mpf_init2(z2abs, PRECISION);
	
	mpf_mul(rr, z1.real, z2.real);
	mpf_mul(ii, z1.img, z2.img);
	mpf_add(temp_result.real, rr, ii);
	
	mpf_mul(ri, z1.real, z2.img);
	mpf_mul(ir, z1.img, z2.real);
	mpf_sub(temp_result.img, ir, ri);
	
	mpf_pow_ui(powr, z2.real, 2);
	mpf_pow_ui(powi, z2.img, 2);
	mpf_add(z2abs, powr, powi);
	
	mpf_div(temp_result.real, temp_result.real, z2abs);
	mpf_div(temp_result.img, temp_result.img, z2abs);
	
	mpf_set(result->real, temp_result.real);
	mpf_set(result->img, temp_result.img);

	complex_fclear(&temp_result);
	mpf_clear(rr);
	mpf_clear(ii);
	mpf_clear(ri);
	mpf_clear(ir);
	mpf_clear(powr);
	mpf_clear(powi);
	mpf_clear(z2abs);
	
	return; 
}
//elevamento a potenza di un complesso, senza usare teorema binomiale
//ora è possibile usare la stessa variabile per operando e risultato
void pow_no_bin(struct complex_f* result, struct complex_f z1, int exp)
{
	mpf_set(result->real, z1.real);
	mpf_set(result->img, z1.img);

	int i;
	for(i = 0; i < exp - 1; i++) //a^N = N-1  moltiplicazioni
	{
		mul_fcomplex(result, *result, z1);
	}
}