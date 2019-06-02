
#ifndef PCH_H
#define PCH_H

#define _USE_MATH_DEFINES
#include <complex>

extern void fft(double *x_in, std::complex<double> *x_out, int N);
void fft_rec(std::complex<double> *x, int N);

#endif //PCH_H
