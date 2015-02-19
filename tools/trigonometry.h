#include <gmp.h>

#ifndef TROGONOMETRY_H
#define TROGONOMETRY_H

/* creates e */
void create_e(mpf_t *e);
/* creates pi */
void create_pi(mpf_t *pi);
/* calculate w = exp(power) using the taylor series */
void myexp(mpf_t *w, mpf_t *power);
/* calculate w = cos(x) using the taylor series */
void mysin(mpf_t *w, mpf_t *x);
/* calculate w = sin(x) using the taylor series */
void mycos(mpf_t *w, mpf_t *x);
/* w = exp(u) */
void complex_exp(mpf_t *u, mpf_t *w);
/* calculate w = u modulo p, where u and w floating point numbers */
void mymod(mpf_t *w, mpf_t *u, mpf_t *p);

#endif
