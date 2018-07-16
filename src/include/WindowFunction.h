/**
 *	@file		src/include/WindowFunction.h
 *	@date		03/05/2018
 *	@author	ATroja
 */

#ifndef WINDOWFUNCTION_H
#define	WINDOWFUNCTION_H

#include <iostream>
#include <complex>

#include "Cosmology.h"

// ----- Galaxy Clustering Window Function ----- //

struct WindowsParameters {

	int l;
	double Chiav;
	double sigmaChi;
	double DC;
	double fC;
	//double (GrowthFactor::*f) (double);
};

class WindowParameters:public GrowthFactor {

	public:

		WindowParameters (std::string pfile):GrowthFactor (pfile) {};
		int l;
		double Chiav;
		double sigmaChi;
		~WindowParameters () {};
};

class WindowFunction:public GrowthFactor {

	public:

		/*
		struct WindowsParameters WP;

		double Wg1;
		double DWg1;
		double DDWg1;
		double DDDWg1;

		double Wg2;
		double DWg2;
		double DDWg2;

		double WgRSD1;

		double GLChi;
		double GQW;
		double Chitab;
		double wChitab;
		*/

		int GLChi;
		double* ChiChi;
		double* wChi;

		//WindowFunction (std::string pfile, int ell, double sc, double cav);
		WindowFunction (std::string pfile, std::string GLQ);

		double Wgalaxy (int l, double Chi, double Chiav, double sigmaChi, int nD, int nG, int nRSD);

		double Wgalaxy_noRSD_1 (double Chi, WindowParameters& WP, int n);
		double Wgalaxy_noRSD_10 (double Chi, WindowParameters& WP);
		double Wgalaxy_noRSD_11 (double Chi, WindowParameters& WP);
		double Wgalaxy_noRSD_12 (double Chi, WindowParameters& WP);
		double Wgalaxy_noRSD_1above (double Chi, WindowParameters& WP);

		double Wgalaxy_noRSD_2 (double Chi, WindowParameters& WP, int n);
		double Wgalaxy_noRSD_20 (double Chi, WindowParameters& WP);
		double Wgalaxy_noRSD_21 (double Chi, WindowParameters& WP);
		double Wgalaxy_noRSD_22 (double Chi, WindowParameters& WP);


		void intWgalaxy (std::complex <double> nu, WindowParameters& WP1, WindowParameters& WP2, int nRSD1, int nRSD2, double t, std::string GLQ, std::complex <double>& result);


		~WindowFunction ();

};

#endif
