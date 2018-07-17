#include "../include/radialselection.h"
#include "../include/distance.h"
#include "../include/growth.h"

#include <iostream>
#include <math.h>

#include <gsl/gsl_integration.h>
#include "gsl/gsl_sf_bessel.h"

using namespace std;

double zmin_rs, zmax_rs, sigma_z;
double zmin_whole, zmax_whole;
double radfunc_norm_inv_factor = 1.;

double inv_radfunc_norm_inv_factor = 1.;

double rs (double s, void* pparams)
{
	double zz = *(double *) pparams;
	return 1. / sqrt (2*M_PI) / sigma_z * exp (- (zz-s) * (zz-s) / 2. / sigma_z / sigma_z);
}

double RadialSelection (double z)
{
	/* RADIAL SELECTION PER IL CASO GENERALE DELLE MICE*/
	double halfz = z/0.5;
	double dndz = halfz * halfz * exp (-pow (halfz,1.5));

	gsl_function RS;
	RS.function = &rs;
	RS.params = &z;

	double res, err;

	gsl_integration_workspace *ww = gsl_integration_workspace_alloc (1000);
	gsl_integration_qags (&RS, zmin_rs, zmax_rs, 0, 1e-7, 1000, ww, &res, &err);
	gsl_integration_workspace_free (ww);

	return inv_radfunc_norm_inv_factor * dndz * res;

}

double radselec (double z, void *params)
{
	return RadialSelection (z);
}

double integrate_RadialSelection (double zmin_now, double zmax_now)
{
	double result, error, alpha = 1.;

	gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);

	gsl_function F;
	F.function = &radselec;
	F.params = &alpha;

	gsl_integration_qags (&F, zmin_now, zmax_now, 0, 1e-7, 1000, w, &result, &error);

	gsl_integration_workspace_free (w);

	return result;
}

void initialize_phi (double zmin, double zmax, double zm, double zw, double sigma_z_in)
{
	sigma_z = sigma_z_in;
	zmin_rs = zmin;
	zmax_rs = zmax;

	zmin_whole = zm-5.*sigma_z;
	zmax_whole = zm+5.*sigma_z;
	
	inv_radfunc_norm_inv_factor = 1./integrate_RadialSelection (zmin_whole, zmax_whole);

}
