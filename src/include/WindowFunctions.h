/**
 *	@file		src/include/WindowFunctions.h
 *	@date		29/03/2018
 *	@author	ATroja
 */

#ifndef WINDOWFUNCTIONS_H
#define WINDOWFUNCTIONS_H

#include <iostream>

class WindowFunctions {

	private:
		std::string pfile;
		int l;
		double sigmaChi;
		double Chi;
		double Chiav;
		double DC;

	public:

		WindowFunctions (std::string pfile, int lin, double sigmaChiin, double Chiin, double Chiavin);

		double Wg1 ();
		double DWg1 ();
		double DDWg1 ();
		double DDDWg1 ();

		double Wg2 ();
		double DWg2 ();
		double DDWg2 ();

		double WgRSD1 ();
		
		~WindowFunctions () {};

};

#endif
