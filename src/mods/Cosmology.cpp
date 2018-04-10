/**
 *	@file		src/mods/Cosmology.cpp
 *	@date		20/03/2018
 *	@author	ATroja
 */

#include "../include/Cosmology.h"
#include "../include/FileUtils.h"
#include "../include/OtherUtils.h"

#include <memory>
#include <cmath>

#include "gsl/gsl_integration.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_sf_hyperg.h"
#include "gsl/gsl_deriv.h"

#include "omp.h"

/// ----- Cosmology constructor ----- //
Cosmology::Cosmology (std::string pfile) {

	if (pfile == "Assassi") {

		c = 3.e8;
		zCMB = 1100.;
		zmax = 1500;

		h = 0.6711;
		H0 = 100.*h*1000.;

		Och2 = 0.12029;
		Oc = Och2 / h / h;
		Obh2 = 0.02207;
		Ob = Obh2 / h / h;
		Om = Ob + Oc;

		ns = 0.96;
		A0 = 2.21536e-9;

		Pfname = "input/PCAMBz0.txt";

		GLQ20 = "input/GaussianQuadratureWeightsAndAbscissae_lmax20.dat";
		GLQ50 = "input/GaussianQuadratureWeightsAndAbscissae_lmax50.dat";
	}

	inv_c = 1./c;
	inv_h = 1./h;
	inv_H0 = 1./H0;

	ChiCMB = Chiz (zCMB);
	Chimax = Chiz (zmax);
	inv_Chimax = 1. / Chimax;

	nChi = 0;
	for (double z = 0.; z <= 99.9; z+=0.1)
		nChi++;
	for (double z = 100; z <= zmax; z+=1.)
		nChi++;

	zarr = new double [nChi];
	Chiarr = new double [nChi];

	int count = 0;
	for (double z = 0.; z <= 99.9; z+=0.1){
		zarr [count] = z;
		Chiarr [count] = Chiz (z);
		count++;
	}
	for (double z = 100.; z <= zmax; z+=1.){
		zarr [count] = z;
		Chiarr [count] = Chiz (z);
		count++;
	}

}

double Cosmology::HH (double z) {

	return H0 * inv_c * inv_h * sqrt (Om * (1.+z) * (1.+z) * (1.+z) + 1. - Om);

}

double inv_HH (double z, void* par) {

	double Om = *(double *) par;
	return 1. / sqrt (Om * (1.+z) * (1.+z) * (1.+z) + 1. - Om);

}

double Cosmology::Chiz (double z) {

	double res, err;

	gsl_function Chiz;
	Chiz.function = &inv_HH;
	Chiz.params = &Om;

	gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
	gsl_integration_qags (&Chiz, 0., z, 0, 1e-7, 1000, w, &res, &err);
	gsl_integration_workspace_free (w);

	return c*h*inv_H0*res;
}

double Cosmology::zChi (double Chi) {

	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *zChi = gsl_spline_alloc (gsl_interp_cspline, nChi);
	gsl_spline_init (zChi, Chiarr, zarr, nChi);

	double res = gsl_spline_eval (zChi, Chi, acc);

	gsl_interp_accel_free (acc);
	gsl_spline_free (zChi);

	return res;
}

double Cosmology::Dz (double z) {

	double a = 1./3.;
	double b = 1.;
	double c = 11./6.;

	double x = (Om - 1.) / Om;
	double x1 = x / (1. + z) / (1. + z) / (1. + z);
	double num, den;

	std::shared_ptr <OtherUtils> other;

	if (fabs (x) < 1)
		den = gsl_sf_hyperg_2F1 (a, b, c, x);
	else
		den = other -> hyp2f1 (a, b, c, x);
	if (fabs (x1) < 1)
		num = gsl_sf_hyperg_2F1 (a, b, c, x1);
	else
		num = other -> hyp2f1 (a, b, c, x1);

	if (z == 0.)
		return 1.;
	else
		return num / den / (1 + z);
}

// ----- Cosmology destructor ----- //
Cosmology::~Cosmology () {

	delete [] zarr;
	delete [] Chiarr;

}

// ----- Geometry constructor ----- //
Geometry::Geometry (std::string pfile):Cosmology (pfile) {

	if (pfile == "Assassi") {

		La = 20;
		GLQ = GLQ20;

	}

	ua = new double [La];
	wa = new double [La];
	ca = new double [La+1];
	dataa = new double [La];

	std::shared_ptr <FileUtils> file;
	file -> read_GLQ (GLQ, ua, wa, La);
	std::shared_ptr <OtherUtils> other;
	other -> BubbleSort_2vec (ua, wa, La);

	int i;
#pragma omp parallel for private (i)
	for (i = 0; i < La; i++) {

		dataa [i] = 1. + zChi ( Chimax * 0.5 * (ua [i] + 1.) );
		dataa [i] = 1./dataa[i];

	}

	double suma;
	for (int l = 0; l < La + 1; l++) {
		suma = 0.;
#pragma omp parallel for private (i) reduction (+:suma)
		for (i = 0; i < La; i++) {
			suma += wa [i] * dataa[i] * gsl_sf_legendre_Pl (l, ua [i]);
		}
		ca [l] = (2.*l+1.) * 0.5 * suma;
	}
}

double Geometry::aChi (double Chi) {

	int i;
	double u = 2. * Chi * inv_Chimax - 1.;
	double aChi =0.;
#pragma omp parallel for private (i) reduction (+:aChi)
	for (i = 0; i < La; i++) {
		aChi += ca [i] * gsl_sf_legendre_Pl (i, u);
	}

	return aChi;
}

// ----- Cosmology destructor ----- //
Geometry::~Geometry () {

	delete [] ua;
	delete [] wa;
	delete [] ca;
	delete [] dataa;

}

// ----- GrowthFactor constructor ----- //
GrowthFactor::GrowthFactor (std::string pfile):Geometry (pfile) {

	if (pfile == "Assassi") {

		LD = 20;
		GLQ = GLQ20;

	}

	uD = new double [LD];
	wD = new double [LD];
	cD = new double [LD+1];
	dataChi = new double [LD];

	std::shared_ptr <FileUtils> file;
	file -> read_GLQ (GLQ, uD, wD, LD);
	std::shared_ptr <OtherUtils> other;
	other -> BubbleSort_2vec (uD, wD, LD);

	int i;
#pragma omp parallel for private (i)
	for (i = 0; i < LD; i++) {

		dataChi [i] = Dz (zChi ( Chimax * 0.5 * (uD [i] + 1.) ));

	}

	double sumD;
	for (int l = 0; l < LD + 1; l++) {
		sumD = 0.;
#pragma omp parallel for private (i) reduction (+:sumD)
		for (i = 0; i < LD; i++) {
			sumD += wD [i] * dataChi[i] * gsl_sf_legendre_Pl (l, uD [i]);
		}
		cD [l] = (2.*l+1.) * 0.5 * sumD;
	}

	pfname = pfile;
}

double GrowthFactor::DChi (double Chi) {

	int i;
	double u = 2. * Chi * inv_Chimax - 1.;
	double DChi =0.;
#pragma omp parallel for private (i) reduction (+:DChi)
	for (i = 0; i < LD; i++) {
		DChi += cD [i] * gsl_sf_legendre_Pl (i, u);
	}

	return DChi;
}

// ----- It follows from https://stackoverflow.com/questions/19450198/calling-gsl-function-inside-a-class-in-a-shared-library , required by fChi function below----- //
struct gsl_f_pars {
	GrowthFactor * pt_MyClass;
};

double gslClassWrapper1 (double x, void* pp) {
	gsl_f_pars* p = (gsl_f_pars *) pp;
	double dc = p->pt_MyClass->DChi (x);
	return log (dc);
}
double gslClassWrapper2 (double x, void* pp) {
	gsl_f_pars* p = (gsl_f_pars *) pp;
	double ac = p->pt_MyClass->aChi (x);
	return log (ac);
}

double GrowthFactor::fChi (double Chi) {

	GrowthFactor gf (pfname);
	struct gsl_f_pars ppp {&gf};

	gsl_function dc;
	dc.function = &gslClassWrapper1;
	dc.params = &ppp;

	double res1, err1;
	gsl_deriv_central (&dc, Chi, 10., &res1, &err1);

	gsl_function ac;
	ac.function = &gslClassWrapper2;
	ac.params = &ppp;
	
	double res2, err2;
	gsl_deriv_central (&ac, Chi, 10., &res2, &err2);

	return  res1/res2;
}



// ----- Cosmology destructor ----- //
GrowthFactor::~GrowthFactor () {

	delete [] uD;
	delete [] wD;
	delete [] cD;
	delete [] dataChi;

}

// ----- LinearMatterPowerSpectrum constructor ----- //
LinearMatterPowerSpectrum::LinearMatterPowerSpectrum (std::string pfile):GrowthFactor (pfile) {

	std::shared_ptr <FileUtils> file;
	nkinput = file -> read_number_of_lines (Pfname);
	kinput = new double [nkinput];
	Pinput = new double [nkinput];
	file -> read_power_spectrum_from_file (Pfname, kinput, Pinput, nkinput);

	kp = kinput [0];
	Pp = Pinput [0];

	nklow = 0;
	for (double k = -10; k < log(kp)/log(10) - 0.01; k += 0.01) {
		nklow ++;
	}

	klow = new double [nklow];
	Plow = new double [nklow];

	double ten2k;
	double inv_kp = 1. / kp;
	int count = 0;
	for (double k = -10; k < log(kp)/log(10) - 0.01; k += 0.01) {

		ten2k = pow (10, k);
		klow [count] = ten2k;
		Plow [count] = Pp * pow ( ten2k *  inv_kp, ns);
		count++;

	}

	nk = nkinput + nklow;
	karr = new double [nk];
	Parr = new double [nk];

	for (int k = 0; k < nklow; k++) {

		karr [k] = klow [k];
		Parr [k] = Plow [k];

	}
	for (int k = nklow; k < nk; k++) {

		karr [k] = kinput [k - nklow];
		Parr [k] = Pinput [k - nklow];

	}


}

double LinearMatterPowerSpectrum::Plin (double k) {

	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *Plin = gsl_spline_alloc (gsl_interp_cspline, nk);
	gsl_spline_init (Plin, karr, Parr, nk);

	double res = gsl_spline_eval (Plin, k, acc);

	gsl_interp_accel_free (acc);
	gsl_spline_free (Plin);

	return res;
}

// ----- Cosmology destructor ----- //
LinearMatterPowerSpectrum::~LinearMatterPowerSpectrum () {

	delete [] kinput;
	delete [] Pinput;
	delete [] klow;
	delete [] Plow;
	delete [] karr;
	delete [] Parr;

}
