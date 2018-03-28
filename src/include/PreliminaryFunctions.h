/**
 *	@file		src/include/PreliminaryFunctions.h
 *	@date		22/03/2018
 *	@author	ATroja
 */

#ifndef PRELIMINARYFUNCTIONS_H
#define	PRELIMINARYFUNCTIONS_H

#include <iostream>
#include <complex>

// ----- Elementary Functions ----- //

class ElementaryFunctions {

	public:

		/**
		 *	@brief	Default empty Constructor
		 */
		ElementaryFunctions () = default;

		double GammaC (double z);
		void GammaC (std::complex<double> z, std::complex<double>& result);

		double GammaRatioC (double z1, double z2);
		void GammaRatioC (std::complex<double> z1, std::complex<double> z2, std::complex<double>& result);

		double Hyp2F1basic (double a, double b, double c, double z);
		void Hyp2F1basic (std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> z, std::complex<double>& s);

		double MinC (double a, double b);

		double MaxC (double a, double b);


		/**
		 *	@brief	Empty destructor
		 */
		~ElementaryFunctions () {};

};

//  ----- Special Functions ----- //

class SpecialFunctions {

	public:

		SpecialFunctions () = default;

		double tminC (int l, std::complex<double> nu);

		void Il (int l, std::complex<double> nu, double t, std::complex <double>& Il);

		~SpecialFunctions () {};

};

#endif
