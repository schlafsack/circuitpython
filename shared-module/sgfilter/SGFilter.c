#include <math.h>
#include "py/runtime.h"
#include "SGFilter.h"

#define EPSILON ((float)(1.0e-20))
#define NCHMAX (256)
#define SWAP(a,b) tempr= (a);(a)= (b);(b)= tempr

void shared_module_sgfilter_sgfilter_construct(sgfilter_sgfilter_obj_t* self, mp_int_t nl, mp_int_t nr, mp_int_t ld, mp_int_t m) {
  self->nl = nl;
  self->nr = nr;
  self->np = nl + nr + 1;
  self->ld = ld;
  self->m = m;
  self->coeffs = fvector(1, self->np);
  sgcoeff(self->coeffs, self->np, self->nl, self->nr, self->ld, self->m);
}

bool shared_module_sgfilter_sgfilter_deinited(sgfilter_sgfilter_obj_t* self) {
  return self->coeffs == mp_const_none;
}

void shared_module_sgfilter_sgfilter_deinit(sgfilter_sgfilter_obj_t* self) {
  if (shared_module_sgfilter_sgfilter_deinited(self)) {
    return;
  }
  free_fvector(self->coeffs, 1, self->nl + self->nr + 1);
  self->coeffs = mp_const_none;
}

void shared_module_sgfilter_sgfilter_filter(sgfilter_sgfilter_obj_t* self, mp_float_t in[], mp_float_t out[], size_t len) {
    sgfilter(self->coeffs, in, out, len, self->nl, self->nr);
}

mp_int_t* ivector(mp_int_t nl, mp_int_t nh) {
    mp_int_t * v;
    mp_int_t k;
    mp_int_t size = (nh - nl + 2);
    v = m_new(mp_int_t, size);
    for (k = 0; k < size; k++) v[k] = 0;
    return v - nl + 1;
}

void free_ivector(mp_int_t * v, mp_int_t nl, mp_int_t nh) {
    m_del(mp_int_t, (char * )(v + nl - 1), (nh - nl + 2));
}

mp_float_t * fvector(mp_int_t nl, mp_int_t nh) {
    mp_float_t * v;
    mp_int_t k;
    mp_int_t size = (nh - nl + 2);
    v = m_new(mp_float_t, size);
    for (k = 0; k < size; k++) v[k] = 0.0;
    return v - nl + 1;
}

void free_fvector(mp_float_t * v, mp_int_t nl, mp_int_t nh) {
    m_del(mp_float_t, (char * )(v + nl - 1), (nh - nl + 2));
}

mp_float_t** fmatrix(mp_int_t nrl, mp_int_t nrh, mp_int_t ncl, mp_int_t nch) {
    mp_int_t i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    mp_float_t ** m;
    m = m_new(mp_float_t *, (nrow + 1));
    m += 1;
    m -= nrl;
    m[nrl] = m_new(mp_float_t, (nrow * ncol + 1));
    m[nrl] += 1;
    m[nrl] -= ncl;
    for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;
    return m;
}

void free_fmatrix(mp_float_t ** m, mp_int_t nrl, mp_int_t nrh, mp_int_t ncl, mp_int_t nch) {
    mp_int_t nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    m_del(mp_float_t, (char * )(m[nrl] + ncl - 1), (nrow * ncol + 1));
    m_del(mp_float_t *, (char * )(m + nrl - 1), (nrow + 1));
}

void lubksb(mp_float_t** a, mp_int_t n, mp_int_t* indx, mp_float_t b[]) {
    mp_int_t i, ii = 0, ip, j;
    mp_float_t sum;

    for (i = 1; i <= n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii)
            for (j = ii; j <= i - 1; j++) sum -= a[i][j] * b[j];
        else if (FP_ZERO != fpclassify(sum))
            ii = i;
        b[i] = sum;
    }
    for (i = n; i >= 1; i--) {
        sum = b[i];
        for (j = i + 1; j <= n; j++) sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}

void ludcmp(mp_float_t** a, mp_int_t n, mp_int_t* indx, mp_float_t* d) {
    mp_int_t i, imax = 0, j, k;
    mp_float_t big, dum, sum, temp;
    mp_float_t * vv;

    vv = fvector(1, n);
    *d = 1.0;
    for (i = 1; i <= n; i++) {
        big = 0.0;
        for (j = 1; j <= n; j++)
            if ((temp = fabsf(a[i][j])) > big)
                big = temp;
        if (FP_ZERO == fpclassify(big)) {
            mp_raise_RuntimeError(translate("Singular matrix found in routine ludcmp()."));
        }
        vv[i] = 1.0 / big;
    }
    for (j = 1; j <= n; j++) {
        for (i = 1; i < j; i++) {
            sum = a[i][j];
            for (k = 1; k < i; k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }
        big = 0.0;
        for (i = j; i <= n; i++) {
            sum = a[i][j];
            for (k = 1; k < j; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            if ((dum = vv[i] * fabsf(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        if (j != imax) {
            for (k = 1; k <= n; k++) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            * d = -( * d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (FP_ZERO == fpclassify(a[j][j]))
            a[j][j] = EPSILON;
        if (j != n) {
            dum = 1.0 / (a[j][j]);
            for (i = j + 1; i <= n; i++) a[i][j] *= dum;
        }
    }
    free_fvector(vv, 1, n);
}

void four1(mp_float_t data[], mp_uint_t nn, mp_int_t isign) {
    mp_uint_t n, mmax, m, j, istep, i;
    mp_float_t wtemp, wr, wpr, wpi, wi, theta;
    mp_float_t tempr, tempi;

    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {
        if (j > i) {
            SWAP(data[j], data[i]);
            SWAP(data[j + 1], data[i + 1]);
        }
        m = n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax = 2;
    while (n > mmax) {
        istep = mmax << 1;
        theta = isign * (6.28318530717959 / mmax);
        wtemp = sinf(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sinf(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                j = i + mmax;
                tempr = wr * data[j] - wi * data[j + 1];
                tempi = wr * data[j + 1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j + 1] = data[i + 1] - tempi;
                data[i] += tempr;
                data[i + 1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}

void twofft(mp_float_t data1[], mp_float_t data2[], mp_float_t fft1[], mp_float_t fft2[], mp_uint_t n) {
    mp_uint_t nn3, nn2, jj, j;
    mp_float_t rep, rem, aip, aim;
    nn3 = 1 + (nn2 = 2 + n + n);
    for (j = 1, jj = 2; j <= n; j++, jj += 2) {
        fft1[jj - 1] = data1[j];
        fft1[jj] = data2[j];
    }
    four1(fft1, n, 1);
    fft2[1] = fft1[2];
    fft1[2] = fft2[2] = 0.0;
    for (j = 3; j <= n + 1; j += 2) {
        rep = 0.5 * (fft1[j] + fft1[nn2 - j]);
        rem = 0.5 * (fft1[j] - fft1[nn2 - j]);
        aip = 0.5 * (fft1[j + 1] + fft1[nn3 - j]);
        aim = 0.5 * (fft1[j + 1] - fft1[nn3 - j]);
        fft1[j] = rep;
        fft1[j + 1] = aim;
        fft1[nn2 - j] = rep;
        fft1[nn3 - j] = -aim;
        fft2[j] = aip;
        fft2[j + 1] = -rem;
        fft2[nn2 - j] = aip;
        fft2[nn3 - j] = rem;
    }
}

void realft(mp_float_t data[], mp_uint_t n, mp_int_t isign) {
    mp_uint_t i, i1, i2, i3, i4, np3;
    mp_float_t c1 = 0.5, c2, h1r, h1i, h2r, h2i;
    mp_float_t wr, wi, wpr, wpi, wtemp, theta;
    theta = 3.141592653589793 / (mp_float_t)(n >> 1);
    if (isign == 1) {
        c2 = -0.5;
        four1(data, n >> 1, 1);
    } else {
        c2 = 0.5;
        theta = -theta;
    }
    wtemp = sinf(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sinf(theta);
    wr = 1.0 + wpr;
    wi = wpi;
    np3 = n + 3;
    for (i = 2; i <= (n >> 2); i++) {
        i4 = 1 + (i3 = np3 - (i2 = 1 + (i1 = i + i - 1)));
        h1r = c1 * (data[i1] + data[i3]);
        h1i = c1 * (data[i2] - data[i4]);
        h2r = -c2 * (data[i2] + data[i4]);
        h2i = c2 * (data[i1] - data[i3]);
        data[i1] = h1r + wr * h2r - wi * h2i;
        data[i2] = h1i + wr * h2i + wi * h2r;
        data[i3] = h1r - wr * h2r + wi * h2i;
        data[i4] = -h1i + wr * h2i + wi * h2r;
        wr = (wtemp = wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
    }
    if (isign == 1) {
        data[1] = (h1r = data[1]) + data[2];
        data[2] = h1r - data[2];
    } else {
        data[1] = c1 * ((h1r = data[1]) + data[2]);
        data[2] = c1 * (h1r - data[2]);
        four1(data, n >> 1, -1);
    }
}

void convlv(mp_float_t data[], mp_uint_t n, mp_float_t respns[], mp_uint_t m, mp_int_t isign, mp_float_t ans[]) {
    mp_uint_t i, no2;
    mp_float_t dum, mag2, * fft;

    fft = fvector(1, n << 1);
    for (i = 1; i <= (m - 1) / 2; i++)
        respns[n + 1 - i] = respns[m + 1 - i];
    for (i = (m + 3) / 2; i <= n - (m - 1) / 2; i++)
        respns[i] = 0.0;
    twofft(data, respns, fft, ans, n);
    no2 = n >> 1;
    for (i = 2; i <= n + 2; i += 2) {
        if (isign == 1) {
            ans[i - 1] = (fft[i - 1] * (dum = ans[i - 1]) - fft[i] * ans[i]) / no2;
            ans[i] = (fft[i] * dum + fft[i - 1] * ans[i]) / no2;
        } else if (isign == -1) {
            if (FP_ZERO == fpclassify(mag2 = ans[i - 1] * ans[i - 1] + ans[i] * ans[i])) {
                mp_raise_RuntimeError(translate("Attempt of deconvolving at zero response in convlv()."));
            }
            ans[i - 1] = (fft[i - 1] * (dum = ans[i - 1]) + fft[i] * ans[i]) / mag2 / no2;
            ans[i] = (fft[i] * dum - fft[i - 1] * ans[i]) / mag2 / no2;
        } else {
            mp_raise_RuntimeError(translate("No meaning for isign in convlv()."));
        }
    }
    ans[2] = ans[n + 1];
    realft(ans, n, -1);
    free_fvector(fft, 1, n << 1);
}

void sgcoeff(mp_float_t c[], mp_int_t np, mp_int_t nl, mp_int_t nr, mp_int_t ld, mp_int_t m) {
    mp_int_t imj, ipj, j, k, kk, mm, * indx;
    mp_float_t d, fac, sum, ** a, * b;
    if (np < nl + nr + 1 || nl < 0 || nr < 0 || ld > m || nl + nr < m) {
        mp_raise_RuntimeError(translate("Inconsistent arguments detected in routine sgcoeff."));
    }
    indx = ivector(1, m + 1);
    a = fmatrix(1, m + 1, 1, m + 1);
    b = fvector(1, m + 1);
    for (ipj = 0; ipj <= (m << 1); ipj++) {
        sum = (ipj ? 0.0 : 1.0);
        for (k = 1; k <= nr; k++) sum += powf((mp_float_t) k, (mp_float_t) ipj);
        for (k = 1; k <= nl; k++) sum += powf((mp_float_t) - k, (mp_float_t) ipj);
        mm = (ipj < 2 * m - ipj ? ipj : 2 * m - ipj);
        for (imj = -mm; imj <= mm; imj += 2) a[1 + (ipj + imj) / 2][1 + (ipj - imj) / 2] = sum;
    }
    ludcmp(a, m + 1, indx, & d);
    for (j = 1; j <= m + 1; j++) b[j] = 0.0;
    b[ld + 1] = 1.0;
    lubksb(a, m + 1, indx, b);
    for (kk = 1; kk <= np; kk++) c[kk] = 0.0;
    for (k = -nl; k <= nr; k++) {
        sum = b[1];
        fac = 1.0;
        for (mm = 1; mm <= m; mm++) sum += b[mm + 1] * (fac *= k);
        kk = ((np - k) % np) + 1;
        c[kk] = sum;
    }
    free_fvector(b, 1, m + 1);
    free_fmatrix(a, 1, m + 1, 1, m + 1);
    free_ivector(indx, 1, m + 1);
}

void sgfilter(mp_float_t c[], mp_float_t yr[], mp_float_t yf[], mp_int_t mm, mp_int_t nl, mp_int_t nr) {
    mp_int_t j, k;
    for (k = 1; k <= nl; k++) {
        for (yf[k] = 0.0, j = -nl; j <= nr; j++) {
            if (k + j >= 1) {
                yf[k] += c[(j >= 0 ? j + 1 : nr + nl + 2 + j)] * yr[k + j];
            }
        }
    }
    for (k = nl + 1; k <= mm - nr; k++) {
        for (yf[k] = 0.0, j = -nl; j <= nr; j++) {
            yf[k] += c[(j >= 0 ? j + 1 : nr + nl + 2 + j)] * yr[k + j];
        }
    }
    for (k = mm - nr + 1; k <= mm; k++) {
        for (yf[k] = 0.0, j = -nl; j <= nr; j++) {
            if (k + j <= mm) {
                yf[k] += c[(j >= 0 ? j + 1 : nr + nl + 2 + j)] * yr[k + j];
            }
        }
    }
}
