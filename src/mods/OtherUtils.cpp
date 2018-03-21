/**
 *	@file		src/mods/OtherUtils.cpp
 *	@date		21/03/2018
 *	@author	ATroja
 */

#include "../include/OtherUtils.h"

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
