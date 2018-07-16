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
	/*
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
	*/

	/* RADIAL SELECTION PER IL CASO PARTICOLARE SENZA PHOTO-Z E SENZA RSD */
	//return inv_radfunc_norm_inv_factor * z * z;
	return inv_radfunc_norm_inv_factor * R(z) * R(z) * 299792.458 / Hz(z);
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

/* SELECTION FUNCTION WITHOUT RSD */
struct my_func_params {double struct_l; double struct_k;};

double selfunc (double z_sf, void *p_sf)
{
	struct my_func_params * par_sf = (struct my_func_params *) p_sf;//it seems to be correct without the * before (struct my_func_params *)
	double l_self = (par_sf -> struct_l);
	double k_self = (par_sf -> struct_k);

	double R_z = R (z_sf);
	double selfunc_res = RadialSelection (z_sf) * GrowthFunc (z_sf) * gsl_sf_bessel_jl (l_self, k_self * R_z);
	
	if (isnan (selfunc_res))
		return 0.;
	else
		return selfunc_res;
}

double SelectionFunction (double l_sf, double k_sf)
{
	double res_sf, err_sf;
	struct my_func_params params_sf = {l_sf, k_sf};
	
	gsl_function F_sf;
	F_sf.function = &selfunc;
	F_sf.params = &params_sf;

	gsl_integration_workspace *w_sf = gsl_integration_workspace_alloc (1000);
	gsl_integration_qags (&F_sf, zmin_whole, zmax_whole, 0, 1e-7, 1000, w_sf, &res_sf, &err_sf);

	gsl_integration_workspace_free (w_sf);

	return res_sf;
}

/* SELECTION FUNCTION WITH RSD */
/*
struct my_func_params {double struct_l; double struct_k; int struct_type;};

double selfunc (double z_sf, void *p_sf)
{
	struct my_func_params * par_sf = (struct my_func_params *) p_sf;//it seems to be correct without the * before (struct my_func_params *)
	double l_self = (par_sf -> struct_l);
	double k_self = (par_sf -> struct_k);

	int type = (par_sf -> struct_type);
	double RSD_Factor;
	if (type == 0)
		RSD_Factor = 1.;
	if (type == 1)
		RSD_Factor = (2.*l_self*l_self+2.*l_self-1.)/((2.*l_self+3.)*(2.*l_self-1.));
	if (type == 2){
		RSD_Factor = (l_self*(l_self-1.))/((2.*l_self-1.)*(2.*l_self+1.));
		l_self -=2.;}
	if (type == 3){
		RSD_Factor = ((l_self+1.)*(l_self+2.))/((2.*l_self+3.)*(2.*l_self+1.));
		l_self +=2.;}

	double R_z = R (z_sf);

	double selfunc_res = RadialSelection (z_sf) * GrowthFunc (z_sf) * gsl_sf_bessel_jl (l_self, k_self * R_z) * RSD_Factor;
	
	if (isnan (selfunc_res))
		return 0.;
	else
		return selfunc_res;
}

double SelectionFunction (double l_sf, double k_sf)
{
	gsl_integration_workspace *w_sf = gsl_integration_workspace_alloc (1000);

	double err_sf;
	double res_sf0, res_sf1, res_sf2, res_sf3;
	
	gsl_function F_sf;
	F_sf.function = &selfunc;

	struct my_func_params params_sf = {l_sf, k_sf,0};
	F_sf.params = &params_sf;
	gsl_integration_qags (&F_sf, zmin_whole, zmax_whole, 0, 1e-7, 1000, w_sf, &res_sf0, &err_sf);

	struct my_func_params params_sf1 = {l_sf, k_sf,1};
	F_sf.params = &params_sf1;
	gsl_integration_qags (&F_sf, zmin_whole, zmax_whole, 0, 1e-7, 1000, w_sf, &res_sf1, &err_sf);

	struct my_func_params params_sf2 = {l_sf, k_sf,2};
	F_sf.params = &params_sf2;
	gsl_integration_qags (&F_sf, zmin_whole, zmax_whole, 0, 1e-7, 1000, w_sf, &res_sf2, &err_sf);

	struct my_func_params params_sf3 = {l_sf, k_sf,3};
	F_sf.params = &params_sf3;
	gsl_integration_qags (&F_sf, zmin_whole, zmax_whole, 0, 1e-7, 1000, w_sf, &res_sf3, &err_sf);

	gsl_integration_workspace_free (w_sf);

	double beta_factor = 0.7047;

	return res_sf0 + beta_factor * (res_sf1 - res_sf2 - res_sf3);
}
*/
