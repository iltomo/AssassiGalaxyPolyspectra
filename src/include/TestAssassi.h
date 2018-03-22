/**
 *	@file		src/include/TestAssassi.h
 *	@date		21/03/2018
 *	@author	ATroja
 */

#ifndef TESTASSASSI_H
#define TESTASSASSI_H

#include <iostream>
#include "../include/Cosmology.h"

class TestAssassi {

	public:

	TestAssassi () = default;

	void TestGeometry (std::string pfile);

	void TestGrowthFactor (std::string pfile);

	void TestLMPS (std::string pfile);

	~TestAssassi () {};

};

#endif
