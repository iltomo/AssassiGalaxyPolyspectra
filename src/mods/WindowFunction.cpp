/**
 *	@file		src/mods/PreliminaryFunctions.cpp
 *	@date		22/03/2018
 *	@author	ATroja
 */

#include "../include/WindowFunction.h"
#include "../include/FileUtils.h"
#include "../include/OtherUtils.h"

#include "boost/multiprecision/cpp_dec_float.hpp"

#include "gsl/gsl_math.h"
#include "gsl/gsl_deriv.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"

#include <memory>

WindowFunction::WindowFunction (std::string pfile, std::string GLQ):GrowthFactor (pfile) {

	GLChi = 50;
	ChiChi = new double [GLChi];
	wChi = new double [GLChi];
	std::shared_ptr <FileUtils> file;
	file -> read_GLQ (GLQ, ChiChi, wChi, GLChi);

};

WindowFunction::~WindowFunction () {
	delete [] ChiChi;
	delete [] wChi;
}

// ----- Derive Function f NOT using GSL ----- //
double newDeriveFunc (double x, void* WP, double (*F) (double, void*)) {
	
	double h = 1e-2;
	double inv2h = 0.5e2;

	return (F(x+h, WP) - F(x-h, WP) )*inv2h;
}
double newDeriveFunc2 (double x, void* WP, double (*F) (double, void*)) {

	double h = 1e-1;
	double inv2h = 0.5e1;

	return  (F(x+2.*h, WP) + F(x-2.*h, WP) - 2.*F(x, WP)) * inv2h * inv2h;
}
void return_f_dfdx_d2fdx2 (double x, void* WP, double (*F) (double, void*), double& f, double& f1, double& f2) {

	/* // Cl
	double h = 1e-5;
	double inv2h = 0.5e5;

	double h2 = 1e-3;
	double inv2h2 = 0.5e3;
	*/

	/* //Quattro termini di bisp
	double h = 1e-3;
	double inv2h = 0.5e3;

	double h2 = 1e-1;
	double inv2h2 = 0.5e1;
	*/

	 //Tutti i termini di bisp //CONTROLLARE SE FUNZIONA ANCHE CON I CL, IN TAL CASO, PROBLEMA RISOLTO
	double h = 1e-2;
	double inv2h = 0.5e2;

	double h2 = 1e-0;
	double inv2h2 = 0.5e0;

	f = F (x, WP);
	f1 = (F(x+h, WP) - F(x-h, WP)) * inv2h;
	f2 = (F(x+2.*h2, WP) + F(x-2.*h2, WP) - 2.*f) * inv2h2 * inv2h2;

}
double newWindowFunctionPolynomial (double (*F) (double, void*), double x, void* WP, int l) {

	double inv_x = 1. / x;

	double f, dfdx, d2fdx2;
	return_f_dfdx_d2fdx2 (x, WP, F, f, dfdx, d2fdx2);
	double result = d2fdx2 - 2. * inv_x * dfdx - (l * (l + 1.) - 2.) * inv_x * inv_x * f;

	return -result;
}

double rs (double s, void* param) {

	return 1. / sqrt (2*M_PI) / sigmaz * exp (- (zz-s) * (zz-s) / 2. / sigmaz / sigmaz);
}
double newWg1_f (double x, void* param) {

	WindowParameters *par = (WindowParameters *) param;
	double DC =  (par -> DChi (x));

	/*
	double sigmaChi = (par -> sigmaChi);
	double Chiav = (par -> Chiav);
	*/

	Survey survey ("Assassi");

	double z = survey.zChi (x);
	double zmin = survey.zmin_bin;
	double zmax = survey.zmax_bin;
	double sigmaz = survey.sigmaz;

	double halfz = z/0.5;
	double dndz = halfz * halfz * exp (-pow (halfz, 1.5));

	gsl_function RS;
	RS.function = &rs;
	RS.params = &param;

	double res, err;
	gsl_integration_workspace *ww = gsl_integration_workspace_alloc (1000);
	gsl_integration_qags (&RS, zmin, zmax, 0, 1e-7, 1000, ww, &res, &err);
	gsl_integration_workspace_free (ww);

	return survey.inv_radFunc_norm_inv_factor * dndz * res * DC * survey.HH (z) * survey.c;
	
}

double newDWg1_f (double x, void* WP) {

	WindowParameters *par = (WindowParameters *) WP;
	int l = (par -> l);

	return newWindowFunctionPolynomial (newWg1_f, x, WP, l);

}

double newDDWg1_f (double x, void* WP) {

	WindowParameters *par = (WindowParameters *) WP;
	int l = (par -> l);

	return newWindowFunctionPolynomial (newDWg1_f, x, WP, l);

}

double newDDDWg1_f (double x, void* WP) {

	WindowParameters *par = (WindowParameters *) WP;
	int l = (par -> l);

	return newWindowFunctionPolynomial (newDDWg1_f, x, WP, l);

}

double newWg1FChi_f (double x, void* WP) {

	WindowParameters *par = (WindowParameters *) WP;
	double fC = (par -> fChi (x));

	return newWg1_f (x, WP) * fC;
}

double newWgRSD1_f (double x, void* WP) {

	return newDeriveFunc2 (x, WP, newWg1FChi_f);
}

double newWg2_f (double x, void* param) {

	WindowParameters *par = (WindowParameters *) param;
	double DC = (par -> DChi (x));

	double sigmaChi = (par -> sigmaChi);
	double Chiav = (par -> Chiav);
	
	return (1. / sqrt (2. * M_PI) / sigmaChi) * exp ( - (x - Chiav) * (x - Chiav) / (2. * sigmaChi * sigmaChi) ) * DC * DC;
}

double newDWg2_f (double x, void* WP) {

	WindowParameters *par = (WindowParameters *) WP;
	int l = (par -> l);

	return newWindowFunctionPolynomial (newWg2_f, x, WP, l);

}

double newDDWg2_f (double x, void* WP) {

	WindowParameters *par = (WindowParameters *) WP;
	int l = (par -> l);

	return newWindowFunctionPolynomial (newDWg2_f, x, WP, l);

}

double Wgalaxy_RSD_11 (double Chi, WindowParameters& WP) {

	return newDWg1_f (Chi, &WP) - newWgRSD1_f (Chi, &WP);
}

double WindowFunction::Wgalaxy_noRSD_1 (double Chi, WindowParameters& WP, int n) {

	double res;
	if (n == 0) {
		res = newWg1_f (Chi, &WP);
	}
	else {
		if (n == 1) {
			res = newDWg1_f (Chi, &WP);
		}
		else {
			if (n == 2) {
				res = newDDWg1_f (Chi, &WP);
			}
			else {
				res = newDDDWg1_f (Chi, &WP);
			}
		}
	}

	return res;

}

double WindowFunction::Wgalaxy_noRSD_2 (double Chi, WindowParameters& WP, int n) {

	double res;
	if (n == 0) {
		res = newWg2_f (Chi, &WP);
	}
	else {
		if (n == 1) {
			res = newDWg2_f (Chi, &WP);
		}
		else {
			res = newDDWg2_f (Chi, &WP);
		}
	}

	return res;

}

double WindowFunction::Wgalaxy_noRSD_20 (double Chi, WindowParameters& WP) {

	return newWg2_f (Chi, &WP);
}
double WindowFunction::Wgalaxy_noRSD_21 (double Chi, WindowParameters& WP) {

	return newDWg2_f (Chi, &WP);
}
double WindowFunction::Wgalaxy_noRSD_22 (double Chi, WindowParameters& WP) {

	return newDWg2_f (Chi, &WP);
}


double WindowFunction::Wgalaxy_noRSD_10 (double Chi, WindowParameters& WP) {

	return newWg1_f (Chi, &WP);
}
double WindowFunction::Wgalaxy_noRSD_11 (double Chi, WindowParameters& WP) {

	return newDWg1_f (Chi, &WP);
}
double WindowFunction::Wgalaxy_noRSD_12 (double Chi, WindowParameters& WP) {

	return newDDWg1_f (Chi, &WP);
}
double WindowFunction::Wgalaxy_noRSD_1above (double Chi, WindowParameters& WP) {

	return newDDDWg1_f (Chi, &WP);
}

// ----- Integral eq. (3.5) ----- //
void WindowFunction::intWgalaxy (std::complex <double> nu, WindowParameters& WP1, WindowParameters& WP2, int nRSD1, int nRSD2, double t, std::string GLQ, std::complex <double>& result) {

	double Chimin = WP2.Chiav - 5. * WP2.sigmaChi;
	double Chimax = WP2.Chiav + 5. * WP2.sigmaChi;

	double dChi = (Chimax - Chimin) / 2.;
	double Chiav = (Chimax + Chimin) / 2.;

	std::complex <double> mickey (0., 0.);
	std::complex <double> goofy (0., 0.);
	std::complex <double> sum (0., 0.);
	std::shared_ptr <OtherUtils> other;


	if (nRSD1 == 0 && nRSD2 == 0) {

		int i;
		for (i = 0; i < GLChi; i++) {

			other -> ComplexPow ((dChi * ChiChi [i] + Chiav), 1. - nu, mickey);
			other -> ComplexPow (t, nu - 2., goofy);
		
			sum += wChi [i] * mickey *
						 Wgalaxy_noRSD_11 (dChi * ChiChi [i] + Chiav, WP2) * ( 
						 Wgalaxy_noRSD_11 ((dChi * ChiChi [i] + Chiav) * t, WP1)  + goofy * 
						 Wgalaxy_noRSD_11 ((dChi * ChiChi [i] + Chiav) / t, WP1)   ); 

		}
	}

	if (nRSD1 == 1 && nRSD2 == 1) {

		int i;
		for (i = 0; i < GLChi; i++) {

			other -> ComplexPow ((dChi * ChiChi [i] + Chiav), 1. - nu, mickey);
			other -> ComplexPow (t, nu - 2., goofy);
		
			sum += wChi [i] * mickey *
						 Wgalaxy_RSD_11 (dChi * ChiChi [i] + Chiav, WP2) * ( 
						 Wgalaxy_RSD_11 ((dChi * ChiChi [i] + Chiav) * t, WP1)  + goofy * 
						 Wgalaxy_RSD_11 ((dChi * ChiChi [i] + Chiav) / t, WP1)   ); 

		}
	}

	result = dChi * sum;
}
