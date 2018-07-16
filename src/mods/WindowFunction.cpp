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

// DERIVATION USING BOOST
boost::multiprecision::cpp_dec_float_50 boostDeriveFunc (double x, void* WP, boost::multiprecision::cpp_dec_float_50 (*F) (double, void*)) {
	
	double h = 1e-2;
	double inv2h = 0.5e2;

	return (F(x+h, WP) - F(x-h, WP) )*inv2h;
}
boost::multiprecision::cpp_dec_float_50 boostDeriveFunc2 (double x, void* WP, boost::multiprecision::cpp_dec_float_50 (*F) (double, void*)) {

	double h = 1e-1;
	double inv2h = 0.5e1;

	return  (F(x+2.*h, WP) + F(x-2.*h, WP) - 2.*F(x, WP)) * inv2h * inv2h;
}
void boostreturn_f_dfdx_d2fdx2 (double x, void* WP, boost::multiprecision::cpp_dec_float_50 (*F) (double, void*), boost::multiprecision::cpp_dec_float_50& f, boost::multiprecision::cpp_dec_float_50& f1, boost::multiprecision::cpp_dec_float_50& f2) {

	double h = 1e-5;
	double inv2h = 0.5e5;

	double h2 = 1e-3;
	double inv2h2 = 0.5e3;

	f = F (x, WP);
	f1 = (F(x+h, WP) - F(x-h, WP)) * inv2h;
	f2 = (F(x+2.*h2, WP) + F(x-2.*h2, WP) - 2.*f) * inv2h2 * inv2h2;

}
boost::multiprecision::cpp_dec_float_50 boostWindowFunctionPolynomial (boost::multiprecision::cpp_dec_float_50 (*F) (double, void*), double x, void* WP, int l) {

	double inv_x = 1. / x;

	boost::multiprecision::cpp_dec_float_50 f, dfdx, d2fdx2;
	boostreturn_f_dfdx_d2fdx2 (x, WP, F, f, dfdx, d2fdx2);
	boost::multiprecision::cpp_dec_float_50 result = d2fdx2 - 2. * inv_x * dfdx - (l * (l + 1.) - 2.) * inv_x * inv_x * f;

	return -result;
}
boost::multiprecision::cpp_dec_float_50 boostWg1_f (double x, void* param) {

	WindowParameters *par = (WindowParameters *) param;
	double DC =  (par -> DChi (x));

	double sigmaChi = (par -> sigmaChi);
	double Chiav = (par -> Chiav);
	
	return (1. / sqrt (2. * M_PI) / sigmaChi) * exp ( - (x - Chiav) * (x - Chiav) / (2. * sigmaChi * sigmaChi) ) * DC;
}

boost::multiprecision::cpp_dec_float_50 boostDWg1_f (double x, void* WP) {

	WindowParameters *par = (WindowParameters *) WP;
	int l = (par -> l);

	return boostWindowFunctionPolynomial (boostWg1_f, x, WP, l);

}

boost::multiprecision::cpp_dec_float_50 boostDDWg1_f (double x, void* WP) {

	WindowParameters *par = (WindowParameters *) WP;
	int l = (par -> l);

	return boostWindowFunctionPolynomial (boostDWg1_f, x, WP, l);

}

boost::multiprecision::cpp_dec_float_50 boostDDDWg1_f (double x, void* WP) {

	WindowParameters *par = (WindowParameters *) WP;
	int l = (par -> l);

	return boostWindowFunctionPolynomial (boostDDWg1_f, x, WP, l);

}

boost::multiprecision::cpp_dec_float_50 boostWg1FChi_f (double x, void* WP) {

	WindowParameters *par = (WindowParameters *) WP;
	double fC = (par -> fChi (x));

	return boostWg1_f (x, WP) * fC;
}

boost::multiprecision::cpp_dec_float_50 boostWgRSD1_f (double x, void* WP) {

	return boostDeriveFunc2 (x, WP, boostWg1FChi_f);
}

boost::multiprecision::cpp_dec_float_50 boostWg2_f (double x, void* param) {

	WindowParameters *par = (WindowParameters *) param;
	double DC = (par -> DChi (x));

	double sigmaChi = (par -> sigmaChi);
	double Chiav = (par -> Chiav);
	
	return (1. / sqrt (2. * M_PI) / sigmaChi) * exp ( - (x - Chiav) * (x - Chiav) / (2. * sigmaChi * sigmaChi) ) * DC * DC;
}

boost::multiprecision::cpp_dec_float_50 boostDWg2_f (double x, void* WP) {

	WindowParameters *par = (WindowParameters *) WP;
	int l = (par -> l);

	return boostWindowFunctionPolynomial (boostWg2_f, x, WP, l);

}

boost::multiprecision::cpp_dec_float_50 boostDDWg2_f (double x, void* WP) {

	WindowParameters *par = (WindowParameters *) WP;
	int l = (par -> l);

	return boostWindowFunctionPolynomial (boostDWg2_f, x, WP, l);

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
double newWg1_f (double x, void* param) {

	WindowParameters *par = (WindowParameters *) param;
	double DC =  (par -> DChi (x));

	double sigmaChi = (par -> sigmaChi);
	double Chiav = (par -> Chiav);
	
	return (1. / sqrt (2. * M_PI) / sigmaChi) * exp ( - (x - Chiav) * (x - Chiav) / (2. * sigmaChi * sigmaChi) ) * DC;
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


//**********************************************
struct DeriveFuncparams {gsl_function func;};

double DeriveFunc (double x, void* param) {

	struct DeriveFuncparams *par = (struct DeriveFuncparams *) param;

	gsl_function f = (par -> func);
	double res, err;

	gsl_deriv_central (&f, x, 1e-5, &res, &err);

	return res;
}

double DeriveFunc2 (double x, void* param) {

	double res, err;

	gsl_function df;
	df.function = &DeriveFunc;
	df.params = param;

	gsl_deriv_central (&df, x, 1e-5, &res, &err);

	return res;
}

double WindowFunctionPolynomial (gsl_function f, double x, int l) {


	struct DeriveFuncparams DFpar {f};

	double inv_x = 1. / x;
	double result = DeriveFunc2 (x, &DFpar);
	result -= 2. * inv_x * DeriveFunc (x, &DFpar);
	result -= (l * (l + 1.) - 2.) * inv_x * inv_x * GSL_FN_EVAL (&f, x);

	return -result;

}

// ----- Gaussian Window function x D(Chi) ----- //
double Wg1_f (double x, void* param) {

	//struct WindowsParameters *par = (struct WindowsParameters *) param;
	//double DC = (par -> DC);
	WindowParameters *par = (WindowParameters *) param;
	double DC =  (par -> DChi (x));

	double sigmaChi = (par -> sigmaChi);
	double Chiav = (par -> Chiav);
	
	return (1. / sqrt (2. * M_PI) / sigmaChi) * exp ( - (x - Chiav) * (x - Chiav) / (2. * sigmaChi * sigmaChi) ) * DC;
}

double DWg1_f (double x, void* WP) {

	gsl_function F;
	F.function = &Wg1_f;
	F.params = WP;

	//struct WindowsParameters *par = (struct WindowsParameters *) WP;
	WindowParameters *par = (WindowParameters *) WP;
	int l = (par -> l);

	return WindowFunctionPolynomial (F, x, l);

}

double DDWg1_f (double x, void* WP) {

	gsl_function F;
	F.function = &DWg1_f;
	F.params = WP;

	//struct WindowsParameters *par = (struct WindowsParameters *) WP;
	WindowParameters *par = (WindowParameters *) WP;
	int l = (par -> l);

	return WindowFunctionPolynomial (F, x, l);

}

double DDDWg1_f (double x, void* WP) {

	gsl_function F;
	F.function = &DDWg1_f;
	F.params = WP;

	//struct WindowsParameters *par = (struct WindowsParameters *) WP;
	WindowParameters *par = (WindowParameters *) WP;
	int l = (par -> l);

	return WindowFunctionPolynomial (F, x, l);

}

double Wg2_f (double x, void* param) {

	//struct WindowsParameters *par = (struct WindowsParameters *) param;
	//double DC = (par -> DC);
	WindowParameters *par = (WindowParameters *) param;
	double DC = (par -> DChi (x));

	double sigmaChi = (par -> sigmaChi);
	double Chiav = (par -> Chiav);
	
	return (1. / sqrt (2. * M_PI) / sigmaChi) * exp ( - (x - Chiav) * (x - Chiav) / (2. * sigmaChi * sigmaChi) ) * DC * DC;
}

double DWg2_f (double x, void* WP) {

	gsl_function F;
	F.function = &Wg2_f;
	F.params = WP;

	//struct WindowsParameters *par = (struct WindowsParameters *) WP;
	WindowParameters *par = (WindowParameters *) WP;
	int l = (par -> l);

	return WindowFunctionPolynomial (F, x, l);

}

double DDWg2_f (double x, void* WP) {

	gsl_function F;
	F.function = &DWg2_f;
	F.params = WP;

	//struct WindowsParameters *par = (struct WindowsParameters *) WP;
	WindowParameters *par = (WindowParameters *) WP;
	int l = (par -> l);

	return WindowFunctionPolynomial (F, x, l);

}

double Wg1FChi_f (double x, void* param) {

	//struct WindowsParameters *par = (struct WindowsParameters *) param;
	//double DC = (par -> DC);
	//double fC = (par -> fC);
	WindowParameters *par = (WindowParameters *) param;
	double DC = (par -> DChi (x));
	double fC = (par -> fChi (x));

	double sigmaChi = (par -> sigmaChi);
	double Chiav = (par -> Chiav);
	
	return (1. / sqrt (2. * M_PI) / sigmaChi) * exp ( - (x - Chiav) * (x - Chiav) / (2. * sigmaChi * sigmaChi) ) * DC * fC;

}

double WgRSD1_f (double x, void* WP) {

	gsl_function F;
	F.function = &Wg1FChi_f;
	F.params = WP;

	struct DeriveFuncparams DFpar {F};

	double result = DeriveFunc2 (x, &DFpar);

	return result ;
}

// ----- GalaxyCluster Constructor ----- //
	/*
WindowFunction::WindowFunction (std::string pfile, int ell, double sc, double cav):GrowthFactor (pfile) {

	WP.l = ell;
	WP.sigmaChi = sc;
	WP.Chiav = cav;
	WP.DC = DChi (Chi);
	WP.fC = fChi (Chi);

	Wg1 = Wg1_f (Chi, &WP);
	DWg1 = DWg1_f (Chi, &WP);
	DDWg1 = DDWg1_f (Chi, &WP);
	DDDWg1 = DDDWg1_f (Chi, &WP);

	Wg2 = Wg2_f (Chi, &WP);
	DWg2 = DWg2_f (Chi, &WP);
	DDWg2 = DDWg2_f (Chi, &WP);

	WgRSD1 = WgRSD1_f (Chi, &WP);

}
	*/

	
// ----- Window Functions for galaxies ----- //
double WindowFunction::Wgalaxy (int l, double Chi, double Chiav, double sigmaChi, int nD, int nG, int nRSD) {

	double res = 0.;

	struct WindowsParameters WP;
	WP.l = l;
	WP.sigmaChi = sigmaChi;
	WP.Chiav = Chiav;
	WP.DC = DChi (Chi);
	WP.fC = fChi (Chi);

	if (nRSD == 1) {
		res = DWg1_f (Chi, &WP) - WgRSD1_f (Chi, &WP);
	}
	else {
		if (nG == 1) {
			if (nD == 0) {
				res = Wg1_f (Chi, &WP);
			}
			if (nD == 1) {
				res = DWg1_f (Chi, &WP);
			}
			if (nD == 2) {
				res = DDWg1_f (Chi, &WP);
			}
			else {
				res = DDDWg1_f (Chi, &WP);
			}
		}
		else {
			if (nD == 0) {
				res = Wg2_f (Chi, &WP);
			}
			if (nD == 1) {
				res = DWg2_f (Chi, &WP);
			}
			else {
				res = DDWg2_f (Chi, &WP);
			}
		}
	}
	return res;
}

double Wgalaxy_RSD_11 (double Chi, WindowParameters& WP) {

	return newDWg1_f (Chi, &WP) - newWgRSD1_f (Chi, &WP);
}

/*
double Wgalaxy_noRSD_11 (double Chi, WindowParameters& WP) {

	return newDWg1_f (Chi, &WP);
}
*/

double WindowFunction::Wgalaxy_noRSD_1 (double Chi, WindowParameters& WP, int n) {

	/*
	double res;

	if (n == 0) res = newWg1_f (Chi, &WP);
	if (n == 1) res = newDWg1_f (Chi, &WP);
	if (n == 2) res = newDDWg1_f (Chi, &WP);
	if (n > 2) res = newDDDWg1_f (Chi, &WP);

	return res;
	*/

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

	/*
	boost::multiprecision::cpp_dec_float_50 res;

	if (n == 0) res = boostWg1_f (Chi, &WP);
	if (n == 1) res = boostDWg1_f (Chi, &WP);
	if (n == 2) res = boostDDWg1_f (Chi, &WP);
	if (n > 2) res = boostDDDWg1_f (Chi, &WP);

	return res.convert_to<double> ();
	*/

}

double WindowFunction::Wgalaxy_noRSD_2 (double Chi, WindowParameters& WP, int n) {

	/*
	double res;

	if (n == 0) res = newWg2_f (Chi, &WP);
	if (n == 1) res = newDWg2_f (Chi, &WP);
	if (n > 1) res = newDDWg2_f (Chi, &WP);

	return res;
	*/

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

	/*
	boost::multiprecision::cpp_dec_float_50 res;

	if (n == 0) res = boostWg2_f (Chi, &WP);
	if (n == 1) res = boostDWg2_f (Chi, &WP);
	if (n > 1) res = boostDDWg2_f (Chi, &WP);

	return res.convert_to<double> ();
	*/
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
