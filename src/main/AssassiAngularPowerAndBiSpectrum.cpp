/**
 *	@file 	src/main/AssassiAngularPowerAndBiSpectrum.cpp
 *	@date		20/03/2018
 *	@author Atroja
 */

#include <iostream>
#include <chrono>
#include <memory>

#include "../include/Cosmology.h"
#include "../include/PreliminaryFunctions.h"
#include "../include/TestAssassi.h"

#include "mpi.h"

class Assassi {

	public:

		void mainMethod (int argc, char** argv) {

			// ----- start time -----//
			auto ts = std::chrono::high_resolution_clock::now ();

			if (argc <= 1) {
				std::cout << "ERROR, missing parameter file" << std::endl;
				exit (1);
			}

			MPI::Init (argc, argv);

			// ----- Testing functions ----- //
			//std::shared_ptr <TestAssassi> test;
			//test -> TestGeometry (argv [1]);
			//test -> TestGrowthFactor (argv [1]);
			//test -> TestLMPS (argv [1]);
			std::shared_ptr <SpecialFunctions> func;
			std::complex<double> nu (20., 2.);
			std::complex<double> res;
			func -> Il (10, nu, 1.2, res);
			std::cout << res << std::endl;
			/*
			std::shared_ptr <SpecialFunctions> spec;
			std::cout << spec -> tminC (25, nu) << std::endl;
			spec -> Il (10, nu, 10., nu);
			*/


			MPI::Finalize ();

			// ----- end time -----//
			auto te = std::chrono::high_resolution_clock::now ();
			std::cout << "The program took ";
			std::cout << std::chrono::duration <double> (te-ts).count ();
			std::cout << std::endl;
		}
};


int main (int argc, char** argv) {

	std::shared_ptr <Assassi> main;

	main -> mainMethod (argc, argv);

	return 0;
}
