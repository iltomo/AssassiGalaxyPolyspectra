/**
 *	@file		src/include/OtherUtils.h
 *	@date		21/03/2018
 *	@author	ATroja
 */

#ifndef OTHERUTILS_H
#define OTHERUTILS_H

#include <iostream>

class OtherUtils {

	public:

		/**
		 *	@brief	Swap two variables
		 *	@param	x,y the variables to swap
		 */
		void Swap (double* xp, double* yp);

		/**
		 *	@brief	Function for the sort of two linked array
		 *	@param	x,y, the array to sort
		 */
		void BubbleSort_2vec (double* x, double* y, int size);

};

#endif

