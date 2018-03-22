/**
 *	@file		src/mods/FileUtils.cpp
 *	@date		21/03/2018
 *	@author	ATroja
 */

#include "../include/FileUtils.h"

#include <fstream>

int FileUtils::read_number_of_lines (std::string fname) {

	int nlines = 0;
	std::string line;
	std::ifstream infile (fname.c_str ());
	if (infile.is_open ()) {

		for (; getline (infile, line);) {

			if (line.at(0) == '#')
				continue;

			nlines++;

		}
	}
	else {
		std::cout << "ERROR in opening " << fname << std::endl;
		exit (1);
	}

	return nlines;
}

void FileUtils::read_GLQ (std::string fname, double* x, double* w, int size){

	/// Dummy variable in order to store the variables I don't care
	double dummy;
	int ncomments = 0;
	std::string line;
	std::ifstream infile (fname.c_str ());
	if (infile.is_open ()) {

		for (; getline (infile, line);) {

			if (line.at(0) == '#') 
				ncomments++;
			else
				break;
		}

		infile.seekg(std::ios::beg);

		for (int i = 0; i < ncomments; i++)
			getline (infile, line);

		for (int i = 0; i < size; i++) {

			infile >> dummy;
			infile >> w [i];
			infile >> x [i];

		}
	}
	else {

		std::cout << "ERROR in opening " << fname << std::endl;
		exit (1);

	}
}

void FileUtils::read_power_spectrum_from_file (std::string Pfname, double* k, double* Pk, int size) {

	int ncomments = 0;
	std::string line;
	std::ifstream infile (Pfname.c_str ());
	if (infile.is_open ()) {

		for (; getline (infile, line);) {

			if (line.at(0) == '#') 
				ncomments++;
			else
				break;
		}

		infile.seekg(std::ios::beg);

		for (int i = 0; i < ncomments; i++)
			getline (infile, line);

		for (int i = 0; i < size; i++) {

			infile >> k [i];
			infile >> Pk [i];

		}
	}
	else {

		std::cout << "ERROR in opening " << Pfname << std::endl;
		exit (1);

	}
}
