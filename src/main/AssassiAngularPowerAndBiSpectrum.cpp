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
#include "../include/PowerSpectrum.h"
#include "../include/TestAssassi.h"
#include "../include/WindowFunction.h"
#include "../include/Cl.h"
#include "../include/Bispectrum.h"

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

			Bispectrum bisp (argv [1]);
			double a = 41./14. * bisp.Bgalaxy (55, 55, 55, 1500., 80., 1500., 80., 1500., 80., 0, 0, 0);
			double b = 3. / 4. * bisp.Bgalaxy (55, 55, 55, 1500., 80., 1500., 80., 1500., 80., 1, -1, 0);
			double c = 3. / 4. * bisp.Bgalaxy (55, 55, 55, 1500., 80., 1500., 80., 1500., 80., -1, 1, 0);
			double d = bisp.Bgalaxy (55, 55, 55, 1500., 80., 1500., 80., 1500., 80., -1, 0, 1);
			double e = bisp.Bgalaxy (55, 55, 55, 1500., 80., 1500., 80., 1500., 80., 0, -1, 1);
			double f = 1. / 4. * bisp.Bgalaxy (55, 55, 55, 1500., 80., 1500., 80., 1500., 80., 1, 1, 2);

			std::cout << a + b + c + d + e + f << std::endl;
			//std::cout << a + d + e + f << std::endl;
			//std::cout << b + c << std::endl;


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
