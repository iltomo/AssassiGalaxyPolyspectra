/**
 *	@file		src/include/PowerSpectrum.h
 *	@date		27/04/2018
 *	@author	ATroja
 */

#ifndef POWERSPECTRUM_H
#define POWERSPECTRUM_H

#include <iostream>

#include "PreliminaryFunctions.h"

class PowerSpectrum:public FourierTransform {

	public:

		PowerSpectrum (std::string pfile):FourierTransform (pfile) {};

		void FFTlogPS ();

		void FFTlogPS_Dirac ();

		void FFTlogPS (double kmin, double kmax, double bias, int halfNmax, std::complex<double>* cn, std::complex<double>* nun);

		~PowerSpectrum () {};
};


#endif
