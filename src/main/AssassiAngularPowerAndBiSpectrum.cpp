/**
 *	@file 	src/main/AssassiAngularPowerAndBiSpectrum.cpp
 *	@date		20/03/2018
 *	@author Atroja
 */

#include <iostream>
#include <chrono>
#include <memory>

#include "../include/Cosmology.h"
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

			std::shared_ptr <TestAssassi> test;
			test -> TestGeometry (argv [1]);


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
