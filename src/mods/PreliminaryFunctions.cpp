/**
 *	@file		src/mods/PreliminaryFunctions.cpp
 *	@date		22/03/2018
 *	@author	ATroja
 */

#include "../include/PreliminaryFunctions.h"

double ElementaryFunctions::GammaC (double z) {

	double result; 
	double q0 = 75122.6331530;
	double q1 = 80916.6278952;
	double q2 = 36308.2951477;
	double q3 = 8687.24529705;
	double q4 = 1168.92649479;
	double q5 = 83.8676043424;
	double q6 = 2.50662827511;

	double p1, p2;
	if (z >= 0.) {

		p1 = ( q0 
				 + q1 * z
				 + q2 * z * z
				 + q3 * z * z * z
				 + q4 * z * z * z * z
				 + q5 * z * z * z * z * z
				 + q6 * z * z * z * z * z * z
				 ) /
				 ( z
			   * (z + 1.)
			   * (z + 2.)
			   * (z + 3.)
			   * (z + 4.)
			   * (z + 5.)
			   * (z + 6.)
				 );

		result = p1 * pow ((z + 5.5), (z + 0.5)) * exp (-z-5.5);
	}
	else {

		p1 = ( q0 
				 + q1 * (1. - z)
				 + q2 * (1. - z) * (1. - z)
				 + q3 * (1. - z) * (1. - z) * (1. - z)
				 + q4 * (1. - z) * (1. - z) * (1. - z) * (1. - z)
				 + q5 * (1. - z) * (1. - z) * (1. - z) * (1. - z) * (1. - z)
				 + q6 * (1. - z) * (1. - z) * (1. - z) * (1. - z) * (1. - z) * (1. - z)
				 ) /
				 ( (1. - z)
			   * (2. - z)
			   * (3. - z)
			   * (4. - z)
			   * (5. - z)
			   * (6. - z)
			   * (7. - z)
				 );

		p2 = p1 * pow ((1. - z + 5.5), (1. - z + 0.5)) * exp (-1. + z - 5.5);
		
		result = M_PI / sin (M_PI * z) / p2;
	}

	return result;
}

double ElementaryFunctions::GammaRatioC (double z1, double z2) {

	int i;
	double result = 1.;
	double q0 = 0.0075122633153;
	double q1 = 0.0809166278952;
	double q2 = 0.363082951477;
	double q3 = 0.868724529705;
	double q4 = 1.16892649479;
	double q5 = 0.838676043424;
	double q6 = 0.250662827511;

#pragma omp parallel for private (i) reduction (*:result)
	for (i = 0; i < 7; i++) {
		result *= (z2 + i) / (z1 + i);
	}

	result *= exp (z2 - z1);

	double p1, p2;
	p1 = (z1 + 0.5) * log (z1 + 5.5) - (z2 + 0.5) * log (z2 + 5.5);
	result *= exp (p1);

	double z1t = z1 / 10.;
	double z2t = z2 / 10.;

	p1 = q0
		 + q1 * z1t
		 + q2 * z1t * z1t
		 + q3 * z1t * z1t * z1t
		 + q4 * z1t * z1t * z1t * z1t
		 + q5 * z1t * z1t * z1t * z1t * z1t
		 + q6 * z1t * z1t * z1t * z1t * z1t * z1t;
	p2 = q0
		 + q1 * z2t
		 + q2 * z2t * z2t
		 + q3 * z2t * z2t * z2t
		 + q4 * z2t * z2t * z2t * z2t
		 + q5 * z2t * z2t * z2t * z2t * z2t
		 + q6 * z2t * z2t * z2t * z2t * z2t * z2t;

	result *= p1 / p2;

	return result;

}

double ElementaryFunctions::Hyp2F1basic (double a, double b, double c, double z) {

	double s = 0.;
	double p = 1.;
	double eps = 1.;
	double n = 0.;
	double sold;

	while (eps > 0.0000000001) {

		sold = s;
		s = s + p;
		p *= (a+n) * (b+n) / (c+n) / (n+1.) * z;
		eps = fabs ((s-sold)/s);
		n++;
	
	}

	return s;
}
