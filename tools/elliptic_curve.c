/*F of characteristic not equal to 2 or 3.
in poche parole p > 3 
(p = 2: 1 + 1 = 0, quindi il campo ha characteristic = 2)
(p = 3: 1 + 1 + 1 = 0, quindi il campo ha characteristic = 3)*/

#include "elliptic_curve.h"

void init(curve_point* P)
{
	mpz_init(P->x);
	mpz_init(P->y);
	mpz_init(P->z);
}

void set(curve_point* Q, mpz_t c[3])
{
	mpz_set(Q->x, c[0]);
	mpz_set(Q->y, c[1]);
	mpz_set(Q->z, c[2]);
}

void copy(curve_point* Q, curve_point P)
{
	mpz_set(Q->x, P.x);
	mpz_set(Q->y, P.y);
	mpz_set(Q->z, P.z);
}

void clear(curve_point* P)
{
	mpz_clear(P->x);
	mpz_clear(P->y);
	mpz_clear(P->z);
}

//This algorithm returns two points (x, y) on E = x^3 + Ax + B
//points deve essere inizializzato come un array a due elementi, ogni elemento (un punto di una curva ellittica) a sua volta inizializzato con init(curve_point* P)
//se (x,y) appartiene a E: y^2 = x^3 + Ax + B allora anche (x, -y) appartiene a E (c'è il termine y^2)
void find_two_points(curve_point points[2], mpz_t A, mpz_t B, mpz_t p_mod)
{
	//1. [Loop]
	gmp_randstate_t state;
	gmp_randinit_mt(state);

	mpz_t x;
	mpz_init(x);
	mpz_t t;
	mpz_init(t);

	mpz_t sqrt_t;
	mpz_init(sqrt_t);

	while(1)
	{
		mpz_urandomm(x, state, p_mod); //Choose random x ∈ [0, p − 1]; - corretto mettere p come terzo parametro
		
		mpz_mul(t, x, x);
		mpz_add(t, t, A);
		mpz_mul(t, t, x);
		mpz_add(t, t, B);
		mpz_mod(t, t, p_mod); // Affine cubic form in x. - t = (x(x^2 + A) + B)) mod p

		if(mpz_jacobi(t, p_mod) == -1)
			continue;

		int ret = radix_GMP(sqrt_t, p_mod, t);

		if(ret == 0) //non dovrebbe capitare mai, visto che sopra c'è il controllo su Jacobi
			continue;

		mpz_set(points[0].x, x);
		mpz_set(points[0].y, sqrt_t);
		mpz_set_ui(points[0].z, 1);

		mpz_neg(sqrt_t, sqrt_t); //-sqrt(t)
		mpz_mod(sqrt_t, sqrt_t, p_mod); //-sqrt(t) mod p =  p - sqrt(t) (visto che sqrt(t) < p)

		mpz_set(points[1].x, x);
		mpz_set(points[1].y, sqrt_t);
		mpz_set_ui(points[1].z, 1);
		
		break;
	}

	gmp_randclear(state);
	mpz_clear(x);
	mpz_clear(t);
	mpz_clear(sqrt_t);

	return;
}

void negate(curve_point* Q, curve_point P) //Q = -P = (x, -y, z)
{
	mpz_set(Q->x, P.x);
	mpz_neg(Q->y, P.y); ////
	mpz_set(Q->z, P.z);
}

int doubleP(curve_point* Q, curve_point P, mpz_t p_mod, mpz_t param_a) //Q = P + P
{
	return add(Q, P, P, p_mod, param_a);
}

//FORSE DEVO FORZARE LA RIDUZIONE MODULO P
//restituisce 0 se è andato tutto ok, 1 se è stato ottenuto il punto all'infinito, -1 in caso di inversione impossibile
int add(curve_point* Q, curve_point P1, curve_point P2, mpz_t p_mod, mpz_t param_a) //Q = P1 + P2
{
	mpz_t m;
	mpz_init(m);

	curve_point Q_temp;
	init(&Q_temp);

	if(mpz_cmp_ui(P1.z, 0) == 0) // Q = O + P2 = P2
	{
		copy(Q, P2);

		mpz_clear(m);
		clear(&Q_temp);
		return 0;	
	}
	if(mpz_cmp_ui(P2.z, 0) == 0) // Q = P1 + O = P1
	{
		copy(Q, P1);

		mpz_clear(m);
		clear(&Q_temp);
		return 0;
	}
	if(mpz_cmp(P1.x, P2.x) == 0) //x1 = x2
	{
		mpz_t temp;
		mpz_init(temp);
		mpz_t inverse;
		mpz_init(inverse);	
		
		mpz_add(temp, P1.y, P2.y);
		mpz_mod(temp, temp, p_mod);

		if(mpz_cmp_ui(temp, 0) == 0) //Q = O
		{
			//mpz_set_ui(Q->x, 0);
			//mpz_set_ui(Q->y, 1);
			//mpz_set_ui(Q->z, 0);
			
			mpz_clear(m);
			clear(&Q_temp);
			mpz_clear(temp);
			mpz_clear(inverse);
			return 1; //PUNTO INFINITO
		}
		
		mpz_mul(temp, P1.x, P1.x);
		mpz_mul_ui(temp, temp, 3);
		mpz_add(temp, temp, param_a); //3*(x1^2) + a

		mpz_mul_ui(inverse, P1.y, 2);
		int ret = mpz_invert(inverse, inverse, p_mod); //(2*y1)^-1

		if(ret == 0) //INVERSO INESISTENTE
		{
			mpz_clear(m);
			clear(&Q_temp);
			mpz_clear(temp);
			mpz_clear(inverse);
			return -1;
		}

		mpz_mul(m, temp, inverse); //(3*(x1^2) + a)(2*y1)^-1
		mpz_mod(m, m, p_mod);

		mpz_clear(temp);
		mpz_clear(inverse);
	}
	else
	{
		mpz_t temp;
		mpz_init(temp);
		mpz_t inverse;
		mpz_init(inverse);

		mpz_sub(temp, P2.y, P1.y);

		mpz_sub(inverse, P2.x, P1.x);

		int ret = mpz_invert(inverse, inverse, p_mod);

		if(ret == 0) //INVERSO INESISTENTE
		{
			mpz_clear(m);
			clear(&Q_temp);
			mpz_clear(temp);
			mpz_clear(inverse);
			return -1;
		}
	
		mpz_mul(m, temp, inverse); //(y2 - y1)((x2 - x1)^(-1))
		mpz_mod(m, m, p_mod);

		mpz_clear(temp);
		mpz_clear(inverse);
	}	

	mpz_mul(Q_temp.x, m, m);
	mpz_sub(Q_temp.x, Q_temp.x, P1.x);
	mpz_sub(Q_temp.x, Q_temp.x, P2.x);
	mpz_mod(Q_temp.x, Q_temp.x, p_mod);

	mpz_sub(Q_temp.y, P1.x, Q_temp.x);
	mpz_mul(Q_temp.y, Q_temp.y, m);
	mpz_sub(Q_temp.y, Q_temp.y, P1.y);
	mpz_mod(Q_temp.y, Q_temp.y, p_mod);

	mpz_set_ui(Q_temp.z, 1);
	
	copy(Q, Q_temp);

	clear(&Q_temp);
	mpz_clear(m);
	return 0;
}

int sub(curve_point* Q, curve_point P1, curve_point P2, mpz_t p_mod, mpz_t param_a) //Q = P1 + (-P2)
{
	curve_point P2_neg;
	init(&P2_neg);
	negate(&P2_neg, P2);

	int ret = add(Q, P1, P2_neg, p_mod, param_a);

	clear(&P2_neg);
	return ret;
}

//restituisce 0 se è andato tutto ok, 1 se è stato ottenuto il punto all'infinito, -1 in caso di inversione impossibile
int mul(curve_point* Q, curve_point P_original, mpz_t n_original, mpz_t p_mod, mpz_t param_a) //Q = [n]P
{
	int ret;

	//1. [Initialize]
	if(mpz_cmp_ui(n_original, 0) == 0) //punto all'infinito
		return 1;
	
	curve_point P;
	init(&P);

	mpz_t n;
	mpz_init(n);	

	if(mpz_cmp_ui(n, 0) < 0) //Q = [-n]P = [n]P
	{	
		negate(&P, P_original); //-P
		mpz_neg(n, n_original); //n lo trasformo positivo
	}
	else
	{
		copy(&P, P_original);
		mpz_set(n, n_original);
	}

	mpz_t m;
	mpz_init(m);
	mpz_mul_ui(m, n, 3);

	int B = mpz_sizeinbase(m, 2); //B numero di bit di m = 3*n

	copy(Q, P);
	
	//2. [Compare bits of 3n, n]
	int j;
	for(j = B - 2; j >= 1; j--)
	{		
		ret = doubleP(Q, *Q, p_mod, param_a);
		
		if(ret != 0) //ret = 0 tutto ok, ret = 1 punto all'infinito, ret = -1 errore inversione
		{
			mpz_clear(m);
			return ret;
		}

		if(mpz_tstbit(m, j) == 1 && mpz_tstbit(n, j) == 0)
		{
			ret = add(Q, *Q, P, p_mod, param_a);

			if(ret != 0) //ret = 0 tutto ok,  ret = 1 punto all'infinito, ret = -1 errore inversione
			{
				mpz_clear(m);
				return ret;
			}
		}
		else if(mpz_tstbit(m, j) == 0 && mpz_tstbit(n, j) == 1)
		{
			ret = sub(Q, *Q, P, p_mod, param_a);

			if(ret != 0) //ret = 0 tutto ok,  ret = 1 punto all'infinito, ret = -1 errore inversione
			{
				mpz_clear(m);
				return ret;
			}
		}

	}

	clear(&P);
	mpz_clear(n);	
	mpz_clear(m);
	return 0;
}