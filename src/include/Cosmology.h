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
		/// CMB redshift
		double zCMB;
		/// Maximum redshift used in this analysis
		double zmax;

		/// Hubble parameter
		double h;
		/// Hubble constant in 100 h m / s / Mpc
		double H0;

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
		 *	@brief	Destructor
		 */
		~Cosmology ();

};

#endif
