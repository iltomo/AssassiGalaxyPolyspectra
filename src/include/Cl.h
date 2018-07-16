/**
 *	@file		src/include/Cl.h
 *	@date		02/05/2018
 *	@author	ATroja
 */

#ifndef CL_H
#define CL_H

#include <iostream>

#include "PowerSpectrum.h"
#include "WindowFunction.h"

class Cl_withoutRSD:public PowerSpectrum {

	private:

		int GLt;
		std::string GQW;
		double* ttab;
		double* wttab;
		std::string fname;

	public:

		Cl_withoutRSD (std::string pfile);

		void Clgnu (std::complex<double> nu, WindowParameters& wp1, WindowParameters& wp2, WindowFunction& WF, std::complex<double>& result);

		double Clg (int l, double Chi1av, double sigmaChi1, double Chi2av, double sigmaChi2, double d);

		~Cl_withoutRSD ();

};

class Cl_withRSD:public PowerSpectrum {

	private:

		int GLt;
		std::string GQW;
		double* ttab;
		double* wttab;
		std::string fname;

	public:

		Cl_withRSD (std::string pfile);

		void Clgnu (std::complex<double> nu, WindowParameters& wp1, WindowParameters& wp2, WindowFunction& WF, std::complex<double>& result);

		double Clg (int l, double Chi1av, double sigmaChi1, double Chi2av, double sigmaChi2, double d);

		~Cl_withRSD ();

};

#endif
