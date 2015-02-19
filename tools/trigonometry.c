#include "trigonometry.h"

//w = exp(u) 
//u = a + i*b
//w = e^u = e^a * (cos b + i sin b)
void complex_exp(mpf_t * u, mpf_t * w)
{
	mpf_t x1;
	mpf_t e;

	//mpf_set_default_prec(num_of_digits);

	mpf_init(x1);
	mpf_init(e);

	myexp(&e, &u[0]);

	mycos(&x1, &u[1]);

	mpf_mul(w[0], e, x1);

	mysin(&x1, &u[1]);
	mpf_mul(w[1], e, x1);

	mpf_clear(x1);
	mpf_clear(e);

}

/* calculate w = exp(power) using the taylor series */
void myexp(mpf_t *w, mpf_t *power)
{
  mpf_t w1;
  mpf_t w2;
  mpf_t w3;
  mpf_t mye;
  mpf_t power1;
  mpf_t b;

  long i, flag = 0;

  //mpf_set_default_prec(num_of_digits);

  mpf_init(w1);
  mpf_init(w2);
  mpf_init(w3);
  mpf_init(mye);
  mpf_init(power1);
  mpf_init(b);

  create_e(&mye);


  mpf_set(power1, *power);

  if( mpf_sgn(power1) == -1) //se esponente negativo
  {
    mpf_neg(power1, power1); //cambio segno dell'esponente
    myexp(&w2, &power1); //calcolo con potenza positiva e poi faccio 1/risultato
    mpf_ui_div(*w, 1, w2); //1/risultato

    flag = 1;
  }

  if( (flag == 0) && mpf_cmp_ui(power1, 1) == -1) //non ho eseguito l'altro if e l'esponente è minore di 1
  {
    mpf_set_ui(w1, 1);
    mpf_set_ui(w3, 0);

   for(i=1; ; i++)
   {
     mpf_add(w2, w3, w1);

     if(mpf_cmp(w3, w2) == 0)
     	break;

     mpf_set(w3, w2);
     mpf_mul(w1, w1, power1);
     mpf_div_ui(w1, w1, i);
   }

   mpf_set(*w, w3);

   flag = 1;
  }

  if( (flag == 0) && mpf_cmp_ui(power1, 0) == 0) //non ho eseguito gli altri if e l'esponente è uguale a 0
  {
   mpf_set_ui(*w, 1); //w = e^0 = 1
   flag = 1;
  }

  if( (flag == 0) && mpf_cmp_ui(power1, 1) == 0) //non ho eseguito gli altri if e l'esponente è uguale a 1
  {
   mpf_set(*w, mye); //w = e^1 = e

   flag = 1;
  }

  if( (flag == 0) && mpf_cmp_ui(power1, 1) == 1 ) //non ho eseguito gli altri if e l'esponente è maggiore di 1
  {
    mpf_set_ui(w2, 1);

    mpf_trunc(w3, power1);

    mpf_set_ui(w1, 2);

   while( mpf_cmp_ui(w3, 0) == 1 )
   {
     mymod(&b, &w3, &w1);

     mpf_div_ui(w3, w3, 2);
     mpf_trunc(w3,w3);

     if( mpf_cmp_ui(b, 1) == 0 )
       mpf_mul(w2, mye, w2);

     mpf_mul(mye, mye, mye);
   }

    mpf_trunc(w3, power1);

    mpf_sub(power1, power1, w3);
    myexp(&w3, &power1);
    mpf_mul(*w, w2, w3);

  }

  mpf_clear(w1);
  mpf_clear(w2);
  mpf_clear(w3);
  mpf_clear(mye);
  mpf_clear(power1);
  mpf_clear(b);

}

/* calculate w = cos(x) using the taylor series */
void mycos(mpf_t *w, mpf_t *x)
{
   long i;

   mpf_t pi, t1, f, n;
   mpf_t help;
   mpf_t t2, s, s1, t;

   //mpf_set_default_prec(num_of_digits);


   mpf_init(pi);
   mpf_init(t1);
   mpf_init(f);
   mpf_init(n);
   mpf_init(help);
   mpf_init(t2);
   mpf_init(s);
   mpf_init(s1);
   mpf_init(t);


   if(mpf_sgn(*x) == 0) //x = 0 ---> w = cos(0) = 1
   {
   	mpf_set_ui(*w, 1);
	return;
   }


   else
   {
   	create_pi(&pi);

	for(;;)
	{
		mpf_div(t1, *x, pi); //t1 = x/pi
		mpf_floor(n, t1); //n = floor(t1)

		mpf_set_str(help, "0.5e0", 10); //0.5

		mpf_sub(f, t1, n); //f = t1 - n = x/pi - floor(x/pi)
		mpf_sub(f, f, help); //f = f - help

		break;

	}

	mpf_mul(f, f, pi);

	//if n is even, we negate f, which negates sin(x)
	mpf_div_ui(help, n, 2);
	mpf_floor(help, help);
	mpf_mul_ui(help, help, 2);

	if(mpf_cmp(help, n) == 0)
		mpf_neg(f, f);


   }

   mpf_set_ui(s, 0);

   mpf_set(t, f);

   for(i = 3; ; i = i+2)
   {
	mpf_add(s1, s, t);

	if(mpf_cmp(s, s1) == 0)
		break;

	mpf_set(s, s1);
	mpf_mul(t, t, f);
	mpf_mul(t, t, f);
	mpf_div_ui(t, t, i-1);
	mpf_div_ui(t, t, i);
	mpf_neg(t, t);

   }

   mpf_set(*w, s);


   mpf_clear(pi);
   mpf_clear(t1);
   mpf_clear(f);
   mpf_clear(n);
   mpf_clear(help);
   mpf_clear(t2);
   mpf_clear(s);
   mpf_clear(s1);
   mpf_clear(t);

}


/* calculate w = sin(x) using the taylor series */
void mysin(mpf_t *w, mpf_t *x)
{
   long i;

   mpf_t pi, t1, f, n;
   mpf_t help;
   mpf_t t2, s, s1, t;


   //mpf_set_default_prec(num_of_digits);

   mpf_init(pi);
   mpf_init(t1);
   mpf_init(f);
   mpf_init(n);
   mpf_init(help);
   mpf_init(t2);
   mpf_init(s);
   mpf_init(s1);
   mpf_init(t);


   if(mpf_sgn(*x) == 0)
   {
   	mpf_set_ui(*w, 0);
	return;
   }

   mpf_pow_ui(help, *x, 2);

   if(mpf_cmp_ui(help, 3) == -1)
   {
	mpf_set(f, *x);

   }

   else
   {
   	create_pi(&pi);
	mpf_set_str(help, "0.5e0", 10);

	for(;;)
	{
		mpf_div(t1, *x, pi);
		mpf_floor(n, t1);
		mpf_sub(f, t1, n);


		if(mpf_cmp(f, help) == 1)
		{
			mpf_add_ui(n, n, 1);
			mpf_sub(f, t1, n);
		}

		break;

	}

	mpf_mul(f, f, pi);

	//if n is odd, we negate f, which negates sin(x)
	mpf_div_ui(help, n, 2);
	mpf_floor(help, help);
	mpf_mul_ui(help, help, 2);

	if(mpf_cmp(help, n) != 0)
		mpf_neg(f, f);


   }

   mpf_set_ui(s, 0);

   mpf_set(t, f);

   for(i = 3; ; i = i+2)
   {
	mpf_add(s1, s, t);

	if(mpf_cmp(s, s1) == 0)
		break;

	mpf_set(s, s1);
	mpf_mul(t, t, f);
	mpf_mul(t, t, f);
	mpf_div_ui(t, t, i-1);
	mpf_div_ui(t, t, i);
	mpf_neg(t, t);

   }

   mpf_set(*w, s);


   mpf_clear(pi);
   mpf_clear(t1);
   mpf_clear(f);
   mpf_clear(n);
   mpf_clear(help);
   mpf_clear(t2);
   mpf_clear(s);
   mpf_clear(s1);
   mpf_clear(t);

}

/* creates e */
void create_e(mpf_t *e)
{
    //mpf_set_str(*e, "2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200305992181741359662904357290033429526059563073813232862794349076323382988075319525101901157383418793070215408914993488416750924476146066808226480016847741185374234544243710753907774499206955170276183860626133138458300075204493382656029760673711320070932870912744374704723069697720931014169283681902551510865746377211125238978442505695369677078544996996794686445490598793163688923009879312773617821542499922957635148220826989519366803318252886939849646510582093923982948879332036250944311730123819706841614039701983767932068328237646480429531180232878250981945581530175671736133206981125099618188159304169035159888851934580727386673858942287922849989208680582574927961048419844436346324496848756023362482704197862320900", 10);
    long i;
    mpf_t s, s1, t;

    mpf_init(s);
    mpf_init(s1);
    mpf_init(t);

    mpf_set_ui(s, 1);
    mpf_set_ui(t, 1);


    for(i = 2; ; i++)
    {
	mpf_add(s1, s, t);

	if(mpf_cmp(s1, s) == 0) break;

	mpf_set(s, s1);
	mpf_div_ui(t, t, i);
    }

    mpf_set(*e, s);

    mpf_clear(s);
    mpf_clear(s1);
    mpf_clear(t);

}

/* creates pi */
void create_pi(mpf_t *pi)
{
    //mpf_set_str(*pi, "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196442881097566593344612847564823378678316527120190914564856692346034861045432664821339360726024914127372458700660631558817488152092096282925409171536436789259036001133053054882046652138414695194151160943305727036575959195309218611738193261179310511854807446237996274956735188575272489122793818301194912983367336244065664308602139494639522473719070217986094370277053921717629317675238467481846766940513200056812714526356082778577134275778960917363717872146844090122495343014654958537105079227968925892354201995611212902196086403441815981362977477130996051870721134999999837297804995105973173281609631859502445945534690830264252230825", 10);
    long i;
    mpf_t sum1, s, s1, t, t1, t25;
    mpf_t g;

    mpf_init(sum1);
    mpf_init(s);
    mpf_init(s1);
    mpf_init(t);
    mpf_init(t1);
    mpf_init(t25);
    mpf_init(g);

    mpf_set_ui(s, 0);
    mpf_set_str(t, "0.5e0", 10);
    mpf_set_str(t1, "0.5e0", 10);
    mpf_set_str(t25, "0.25e0", 10);

    mpf_neg(t25, t25);

    for(i = 3; ; i+=2)
    {
	mpf_add(s1, s, t);

	if(mpf_cmp(s1, s) == 0) break;

	mpf_set(s, s1);

	mpf_mul(t1, t1, t25);
	mpf_div_ui(t, t1, i);
    }

    mpf_set(sum1, s);

    mpf_set_ui(t25, 1);
    mpf_div_ui(g, t25, 3);

    mpf_set_ui(s, 0);
    mpf_set(t, g);
    mpf_set(t1, g);

    mpf_pow_ui(g, g, 2);
    mpf_neg(g, g);

    for(i = 3; ; i+=2)
    {
	mpf_add(s1, s, t);

	if(mpf_cmp(s1, s) == 0) break;

	mpf_set(s, s1);

	mpf_mul(t1, t1, g);
	mpf_div_ui(t, t1, i);
    }

    mpf_add(s, s, sum1);
    mpf_mul_ui(s, s, 4);

    mpf_set(*pi, s);


    mpf_clear(sum1);
    mpf_clear(s);
    mpf_clear(s1);
    mpf_clear(t);
    mpf_clear(t1);
    mpf_clear(t25);
    mpf_clear(g);

}

/* calculate w = u modulo p, where u and w floating point numbers */
void mymod(mpf_t * w, mpf_t * u, mpf_t * p)
{

	mpf_t w1;

	//mpf_set_default_prec(num_of_digits);

	mpf_init(w1);

	mpf_div(w1, *u, *p);

	mpf_trunc(w1, w1);
	mpf_mul(w1, *p, w1);
	mpf_sub(*w, *u, w1);

	mpf_clear(w1);
}
