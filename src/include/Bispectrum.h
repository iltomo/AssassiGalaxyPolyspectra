/**
 *	@file		src/include/Bispectrum.h
 *	@date		25/05/2018
 *	@author	ATroja
 */

#ifndef BISPECTRUM_H
#define BISPECTRUM_H

#include <iostream>

#include "PowerSpectrum.h"
#include "WindowFunction.h"

class Bispectrum:public PowerSpectrum {

	private:

		int GLt;
		std::string GQW;
		double* ttab;
		double* wttab;
		std::string fname;

	public:

		Bispectrum (std::string  pfile);

		void IntIlgalaxynu (std::complex<double> nu, double r, WindowParameters& WP, WindowFunction& WF, int n, std::complex<double>& result);

		double Ilgalaxy (WindowParameters& wp, double r, int n, WindowFunction& WF);

		double Bgalaxy (int l1, int l2, int l3, double Chi1av, double sigmaChi1, double Chi2av, double sigmaChi2, double Chi3av, double sigmaChi3, int n1, int n2, int n3);

		~Bispectrum ();

};

#endif
