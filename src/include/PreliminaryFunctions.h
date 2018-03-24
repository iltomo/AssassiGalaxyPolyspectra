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

		double GammaRatioC (double z1, double z2);

		double Hyp2F1basic (double a, double b, double c, double z);


		/**
		 *	@brief	Empty destructor
		 */
		~ElementaryFunctions () {};

};

#endif
