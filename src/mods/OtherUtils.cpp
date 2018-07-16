/**
 *	@file		src/mods/OtherUtils.cpp
 *	@date		21/03/2018
 *	@author	ATroja
 */

#include "../include/OtherUtils.h"

#include <cmath>

#include "gsl/gsl_sf_hyperg.h"
#include "gsl/gsl_sf_gamma.h"

void OtherUtils::Swap (double* xp, double* yp) {

	double temp = *xp;
	*xp = *yp;
	*yp = temp;

}

void OtherUtils::BubbleSort_2vec (double* x, double* y, int size) {

	int i, j;
	for (i = 0; i < size-1; i++) {
		for (j = 0; j < size-i-1; j++) {
			if (x[j] > x[j+1]) {
				Swap(&x[j], &x[j+1]);
				Swap(&y[j], &y[j+1]);
			}
		}
	}
}

double OtherUtils::hyp2f1 (double a, double b, double c, double x_hyp) {

	double prefac1 = gsl_sf_gamma (b-a) * gsl_sf_gamma (c) * pow (-x_hyp, -a) / gsl_sf_gamma (b) / gsl_sf_gamma (c-a);
	double prefac2 = gsl_sf_gamma (a-b) * gsl_sf_gamma (c) * pow (-x_hyp, -b) / gsl_sf_gamma (a) / gsl_sf_gamma (c-b);

	double inv_x_hyp = 1. / x_hyp;

	return prefac1 * gsl_sf_hyperg_2F1 (a, a-c+1, a-b+1, inv_x_hyp) + prefac2 * gsl_sf_hyperg_2F1 (b, b-c+1, b-a+1, inv_x_hyp);
}

void OtherUtils::ComplexPow (double x, std::complex<double> y, std::complex<double>& z    ) {

	double x2y = pow (x, y.real ());
	double lnx = log (x);
	 
	z.real (x2y * cos (y.imag () * lnx));
	z.imag (x2y * sin (y.imag () * lnx));
	 
}

