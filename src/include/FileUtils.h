/**
 *	@file		src/include/FileUtils.h
 *	@date		21/03/2018
 *	@author	ATroja
 */

#ifndef FILEUTILS_H
#define FILEUTILS_H

#include <iostream>

class FileUtils {

	public:

		/**
		 *	@brief	Default empty  constructor
		 *	@return	FileUtils object
		 */
		FileUtils () = default;

		int read_number_of_lines (std::string fname);

		void read_GLQ (std::string fname, double* x, double* w, int size);

		/**
		 *	@brief	Empty destructor
		 */
		~FileUtils () {};

};

#endif
