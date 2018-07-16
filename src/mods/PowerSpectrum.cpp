/**
 *	@file		src/mods/PowerSpectrum.cpp
 *	@date		27/04/2018
 *	@author	ATroja
 */

#include "../include/PowerSpectrum.h"
#include "../include/OtherUtils.h"

#include <memory>

#include <cmath>

double deltaP (double k, double pk) {
	return 1./k * pk * exp (-k*k*0.01);
}

double deltaPFT (double k, int halfNmax, const std::complex<double>* cn, const std::complex<double>* nun) {

	std::complex<double> dpft;
	std::complex<double> sum (0.,0.);
	std::complex<double> mickey, goofy;

	std::shared_ptr <OtherUtils> oth;

	for (int n = 1; n <= halfNmax; n++) {
		oth -> ComplexPow (k, nun [halfNmax + n], mickey);
		sum += cn [halfNmax + n] * mickey;
	}

	oth -> ComplexPow (k, nun [halfNmax], goofy);
	dpft = cn [halfNmax] * goofy + 2. * sum;

	return dpft.real ();
}


void PowerSpectrum::FFTlogPS () {

	double kmin = 1e-8;
	double kmax = 52.;
	double bias = -1.9;

	double cst = deltaP (2.*kmin, Plin (2.*kmin)) * pow (2., -bias);

	int Nmax = 200;
	int halfNmax = Nmax / 2;

	std::complex<double>* cn = new std::complex<double> [Nmax+1];
	std::complex<double>* nun = new std::complex<double> [Nmax+1];
	CoeffTransfer (deltaP, bias, cst, Nmax, kmin, kmax, cn, nun);

	std::cout << "nu_max = " << fabs (nun [0].imag ()) << std::endl;

	double prefac, dp, dpft;
	for (double k = 1e-4; k <= 10; k+= 1e-4) {
		prefac = pow (k/kmin, -bias);
		dp = deltaP (k, Plin (k));
		dpft = deltaPFT (k, halfNmax, cn, nun); 

		std::cout << k << "\t" << prefac * dp + cst << "\t" << prefac * dpft << "\t" << dpft / dp - 1. << std::endl;
	}

	delete [] cn;
	delete [] nun;

}

double deltaP_Dirac (double k, double pk) {
	return pk * exp (-k*k*0.01);
}

void PowerSpectrum::FFTlogPS_Dirac () {

	double kmin = 1e-8;
	double kmax = 200.;
	double bias = -1.5;

	double cst = deltaP_Dirac (2.*kmin, Plin (2.*kmin)) * pow (2., -bias);

	int Nmax = 200;
	int halfNmax = Nmax / 2;

	std::complex<double>* cn = new std::complex<double> [Nmax+1];
	std::complex<double>* nun = new std::complex<double> [Nmax+1];
	CoeffTransfer (deltaP_Dirac, bias, cst, Nmax, kmin, kmax, cn, nun);

	std::cout << "nu_max = " << fabs (nun [0].imag ()) << std::endl;

	double prefac, dp, dpft;
	for (double k = 1e-4; k <= 10.; k+= 1e-4) {
		prefac = pow (k/kmin, -bias);
		dp = deltaP_Dirac (k, Plin (k));
		dpft = deltaPFT (k, halfNmax, cn, nun); 

		std::cout << k << "\t" << prefac * dp + cst << "\t" << prefac * dpft << "\t" << dpft / dp - 1. << std::endl;
	}

	delete [] cn;
	delete [] nun;

}

void PowerSpectrum::FFTlogPS (double kmin, double kmax, double bias, int Nmax, std::complex<double>* cn, std::complex<double>* nun) {

	double cst = deltaP (2.*kmin, Plin (2.*kmin)) * pow (2., -bias);

	CoeffTransfer (deltaP, bias, cst, Nmax, kmin, kmax, cn, nun);

}
