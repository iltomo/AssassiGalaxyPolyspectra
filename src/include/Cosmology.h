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

		/**
		 *	@brief	Interpolation of the comoving distance over the redshift
		 *	@param	Chi, the comoving distance
		 *	@return	The redshift as a function of the comoving distance
		 */
		double zChi (double Chi);


		/**
		 *	@brief	Destructor
		 */
		~Cosmology ();

};

#endif
