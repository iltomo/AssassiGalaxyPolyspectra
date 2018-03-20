/**
 *	@file		src/mods/Cosmology.cpp
 *	@date		20/03/2018
 *	@author	ATroja
 */

#include "../include/Cosmology.h"

#include <cmath>

#include "gsl/gsl_integration.h"

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

	}

	inv_c = 1./c;
	inv_h = 1./h;
	inv_H0 = 1./H0;

	ChiCMB = Chiz (zCMB);
	Chimax = Chiz (zmax);
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


	return 1.;
}

/// ----- Cosmology destructor ----- //
Cosmology::~Cosmology () {

}
