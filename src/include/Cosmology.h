/**
 *	@file		src/include/Cosmology.h
 *	@date		20/03/2018
 *	$author	ATroja
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
		 *	@brief	Destructor
		 */
		~Geometry ();

};


#endif
