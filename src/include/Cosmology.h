/**
 *	@file		src/include/Cosmology.h
 *	@date		20/03/2018
 *	@author	ATroja
 */

#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include <iostream>

class Cosmology {

	public:

		/**
		 *	@brief	Initialized Constructor
		 *	@param	pfile, the name of the parameter file
		 *	@return	Initialized Cosmology object
		 */
		Cosmology (std::string pfile);

		/// Speed of light in m/s
		double c;
		/// Inverse of the speed of light in s/m
		double inv_c;
		/// CMB redshift
		double zCMB;
		/// Maximum redshift used in this analysis
		double zmax;

		/// Hubble parameter
		double h;
		/// Inverse of the Hubble parameter
		double inv_h;
		/// Hubble constant in 100 h m / s / Mpc
		double H0;
		/// Inverse of the Hubble constant in (100 h m / s / Mpc)^-1
		double inv_H0;

		/// Cold Dark Matter density in unity of h^2
		double Och2;
		/// Cold Dark Matter density
		double Oc;
		/// Baryon  density in unity of h^2
		double Obh2;
		/// Baryon density
		double Ob;
		/// Total Matter density
		double Om;

		/// Spectral index
		double ns;
		/// Spectral amplitude A0 = k^3/(2*M_PI^2) * P_xi (k)
		double A0;

		///	Filename of the linear matter power spectrum in input
		std::string Pfname;

		/**
		 *	@brief 	Hubble as a function of the redshift
		 *	@param	z, the redshift
		 *	@return	The Hubble parameter at redshift z
		 */
		double HH (double z);
		/**
		 *	@brief	Comoving distance as a function of the redshift
		 *	@param	z, the redshift
		 *	@return	The comoving distance at redshift z
		 */
		double Chiz (double z);

		/// Comoving distance to the CMB
		double ChiCMB;
		///	Comoving distance to the maximum redshift in this analysis
		double Chimax;
		///	Inverse of the comoving distance to the maximum redshift in this analysis
		double inv_Chimax;

		/// Size of the z and Chi(z) arrays
		int nChi;
		/// Array containing values of redshift
		double* zarr;
		/// Array containing values of Chi (z)
		double* Chiarr;

		/**
		 *	@brief	Interpolation of the comoving distance over the redshift
		 *	@param	Chi, the comoving distance
		 *	@return	The redshift as a function of the comoving distance
		 */
		double zChi (double Chi);

		/// Filename for the Gauss-Legendre quadrature with lmax = 20
		std::string GLQ20;
		/// Filename for the Gauss-Legendre quadrature with lmax = 50
		std::string GLQ50;

		/**
		 *	@brief	Function that evaluates the growth factor
		 *	@param	z, the redshift
		 *	@return	The growth factor as a function of redshift z
		 */
		double Dz (double z);

		/**
		 *	@brief	Destructor
		 */
		~Cosmology ();

};

class Geometry:public Cosmology {

	private:

		/// Maximum degree in the Legendre expansion
		int La;
		/// Filename for the Gauss-Legendre quadrature with lmax = 20
		std::string GLQ;
		/// Abscissae for the Gauss-Legendre quadrature
		double* ua;
		/// Weight for the Gauss-Legendre quadrature
		double* wa;
		/// Coefficients for the Legendre polynomials
		double* ca;
		/// List of scale factor for the abscissae ua
		double* dataa;

	public:

		/**
		 *	@brief	Initialized Constructor
		 *	@param	pfile, the name of the parameter file
		 *	@return	Initialized Geometry object
		 */
		Geometry (std::string pfile);

		/**
		 *	@brief	This function evaluates the scale factor at a given comoving distance
		 *	@param	Chi, the comoving distance
		 *	@return	The scale factor at Chi
		 */
		double aChi (double Chi);


		/**
		 *	@brief	Destructor
		 */
		~Geometry ();

};

class GrowthFactor:public Geometry {

	private:

		/// Maximum degree in the Legendre expansion
		int LD;
		/// Filename for the Gauss-Legendre quadrature with lmax = 20
		std::string GLQ;
		/// Abscissae for the Gauss-Legendre quadrature
		double* uD;
		/// Weight for the Gauss-Legendre quadrature
		double* wD;
		/// Coefficients for the Legendre polynomials
		double* cD;
		/// List of scale factor for the abscissae ua
		double* dataChi;

		/// The name of the parameter file is stored for being used in fChi function
		std::string pfname;

	public:

		/**
		 *	@brief	Initialized Constructor
		 *	@param	pfile, the name of the parameter file
		 *	@return	Initialized GrowthFactor object
		 */
		GrowthFactor (std::string pfile);

		/**
		 *	@brief	This function evaluates the growth factor at a given comoving distance
		 *	@param	Chi, the comoving distance
		 *	@return	The growth factor at Chi
		 */
		double DChi (double Chi);

		double fChi (double Chi);


		/**
		 *	@brief	Destructor
		 */
		~GrowthFactor ();

};

class LinearMatterPowerSpectrum:public GrowthFactor {

	private:

		///	Number of k in the Power Spectrum file given in input
		int nkinput;
		///	Array containing the k in the Power Spectrum file
		double* kinput;
		///	Array containing the P(k) in the Power Spectrum file
		double* Pinput;

		/// Minimum k in the input Power Spectrum file
		double kp;
		/// Minimum P(k) in the input Power Spectrum file
		double Pp;

		/// Number of k at larger scale than one in input
		int nklow;
		/// Array of k smaller than k in input
		double* klow;
		/// Array of P(k) for k smaller than k in input
		double* Plow;

		/// Total number of k used for interpolation
		int nk;
		/// Array of k for interpolation
		double* karr;
		/// Array of P(k) for interpolation
		double* Parr;

	public:

		/**
		 *	@brief	Initialized Constructor
		 *	@param	pfile, the name of the parameter file
		 *	@return	Initialized LinearMatterPowerSpectrum object
		 */
		LinearMatterPowerSpectrum (std::string pfile);

		/**
		 *	@brief	This function interpolate the Power Spectrum at a certain k
		 *	@param	k, the wave number
		 *	@return	The Power spectrum at k, P(k)
		 */
		double Plin (double k);


		/**
		 *	@brief	Destructor
		 */
		~LinearMatterPowerSpectrum ();

};

#endif
