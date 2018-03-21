/**
 *	@file		src/include/TestAssassi.cpp
 *	@date		21/03/2018
 *	@author	ATroja
 */

#include "../include/TestAssassi.h"
#include "../include/Cosmology.h"

void TestAssassi::TestGeometry (std::string pfile) {
	
	Geometry cosmo (pfile);
	double a,z;
	for (double x = 0; x <= cosmo.Chimax; x+=1.) {
		a = cosmo.aChi (x);
		z = 1./(1.+cosmo.zChi (x));
		std::cout << a << "\t" << z << "\t" << 100. * (a-z) / z << std::endl;
	} 
}
