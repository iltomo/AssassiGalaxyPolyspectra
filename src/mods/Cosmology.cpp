/**
 *	@file		src/mods/Cosmology.cpp
 *	@date		20/03/2018
 *	@author	ATroja
 */

#include "../include/Cosmology.h"

#include <cmath>

#include "gsl/gsl_integration.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_sf_legendre.h"

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

		GLQ20 = "input/GaussianQuadratureWeightsAndAbscissae_lmax20.dat";
		GLQ50 = "input/GaussianQuadratureWeightsAndAbscissae_lmax50.dat";
	}

	inv_c = 1./c;
	inv_h = 1./h;
	inv_H0 = 1./H0;

	ChiCMB = Chiz (zCMB);
	Chimax = Chiz (zmax);

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

	read_GLQ (GLQ, ua, wa, La);


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

// ----- Cosmology destructor ----- //
Geometry::~Geometry () {

	delete [] ua;
	delete [] wa;
	delete [] ca;
	delete [] dataa;

}
