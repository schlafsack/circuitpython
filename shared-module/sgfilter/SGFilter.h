#ifndef MICROPY_INCLUDED_SGFILTER_SGFILTER_H
#define MICROPY_INCLUDED_SGFILTER_SGFILTER_H

#include <stdbool.h>
#include <stdint.h>

#include "py/obj.h"

typedef struct {
    mp_obj_base_t base;
    mp_int_t nl;
    mp_int_t nr;
    mp_int_t np;
    mp_int_t ld;
    mp_int_t m;
    mp_float_t * coeffs;
} sgfilter_sgfilter_obj_t;

mp_int_t* ivector(mp_int_t nl, mp_int_t nh);
void free_ivector(mp_int_t* v, mp_int_t nl, mp_int_t nh);
mp_float_t* fvector(mp_int_t nl, mp_int_t nh);
void free_fvector(mp_float_t* v, mp_int_t nl, mp_int_t nh);
mp_float_t **fmatrix(mp_int_t nrl, mp_int_t nrh, mp_int_t ncl, mp_int_t nch);
void free_fmatrix(mp_float_t** m, mp_int_t nrl, mp_int_t nrh, mp_int_t ncl, mp_int_t nch);
void lubksb(mp_float_t** a, mp_int_t n, mp_int_t* indx, mp_float_t b[]);
void ludcmp(mp_float_t** a, mp_int_t n, mp_int_t* indx, mp_float_t* d);
void four1(mp_float_t data[], mp_uint_t nn, mp_int_t isign);
void twofft(mp_float_t data1[], mp_float_t data2[], mp_float_t fft1[], mp_float_t fft2[], mp_uint_t n);
void realft(mp_float_t data[], mp_uint_t n, mp_int_t isign);
void convlv(mp_float_t data[], mp_uint_t n, mp_float_t respns[], mp_uint_t m, mp_int_t isign, mp_float_t ans[]);
void sgcoeff(mp_float_t c[], mp_int_t np, mp_int_t nl, mp_int_t nr, mp_int_t ld, mp_int_t m);
void sgfilter(mp_float_t c[], mp_float_t yr[], mp_float_t yf[], mp_int_t mm, mp_int_t nl, mp_int_t nr);

#endif // MICROPY_INCLUDED_SGFILTER_SGFILTER_H
