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

void TestAssassi::TestGrowthFactor (std::string pfile) {
	
	GrowthFactor cosmo (pfile);
	double a,z;
	for (double x = 0; x <= cosmo.Chimax; x+=1.) {
		a = cosmo.DChi (x);
		z = cosmo.Dz (cosmo.zChi(x));
		std::cout << a << "\t" << z << "\t" << 100. * (a-z) / z << std::endl;
	} 
}

void TestAssassi::TestLMPS (std::string pfile) {

	LinearMatterPowerSpectrum cosmo (pfile);
	
	for (double k = 1e-5; k < 2; k += 1e-4)
		std::cout << k << "\t" << cosmo.Plin (k) << std::endl;
}
