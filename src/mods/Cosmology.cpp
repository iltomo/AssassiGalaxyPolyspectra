/**
 *	@file		src/mods/Cosmology.cpp
 *	@date		20/03/2018
 *	@author	ATroja
 */

#include "../include/Cosmology.h"

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
		Obh2 = 0.00207;
		Ob = Obh2 / h / h;
		Om = Ob + Oc;

		ns = 0.96;
		A0 = 2.21536e-9;

	}
}

/// ----- Cosmology destructor ----- //
Cosmology::~Cosmology () {

}
