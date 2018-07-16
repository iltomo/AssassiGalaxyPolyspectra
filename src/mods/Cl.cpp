/**
 *	@file		src/mods/Cl.cpp
 *	@date		02/05/2018
 *	@author	ATroja
 */

#include "../include/Cl.h"
//#include "../include/WindowFunction.h"
#include "../include/FileUtils.h"
#include "../include/OtherUtils.h"

#include <memory>

Cl_withoutRSD::Cl_withoutRSD (std::string pfile):PowerSpectrum (pfile) {

	if (pfile == "Assassi") {

		GLt = 50;
		GQW = GLQ50;

	}

	ttab = new double [GLt];
	wttab = new double [GLt];

	std::shared_ptr <FileUtils> file;
	file -> read_GLQ (GQW, ttab, wttab, GLt);
	std::shared_ptr <OtherUtils> other;
	other -> BubbleSort_2vec (ttab, wttab, GLt);

	fname = pfile;

}

void Cl_withoutRSD::Clgnu (std::complex<double> nu, WindowParameters& wp1, WindowParameters& wp2, WindowFunction& WF, std::complex<double>& result) {

	double Chi1min = wp1.Chiav - 5. * wp1.sigmaChi;
	double Chi1max = wp1.Chiav + 5. * wp1.sigmaChi;

	double Chi2min = wp2.Chiav - 5. * wp2.sigmaChi;
	double Chi2max = wp2.Chiav + 5. * wp2.sigmaChi;

	double tmin, tmax;
	if (wp1.l > 10) {
		std::shared_ptr <SpecialFunctions> sf;
		double zmin = sf -> tminC (wp1.l, nu);
		std::shared_ptr <ElementaryFunctions> ef;
		tmin = ef -> MaxC (Chi1min / Chi2max, zmin);
		tmax = ef -> MinC (Chi1max / Chi2min, 1.);
	}
	else {
		tmin = 0.;
		tmax = 1.;
	}

	double deltat = (tmax - tmin) * 0.5;
	double tav = (tmax + tmin) * 0.5;

	std::complex<double> res (0., 0.);

	std::shared_ptr <SpecialFunctions> sf;

	std::complex<double> il, iw;
	for (int it = 0; it < GLt; it++) {

		sf -> Il (wp1.l, nu + 1e-10, deltat * ttab[it] + tav, il);

		WF.intWgalaxy (nu, wp1, wp2, 0, 0, deltat * ttab [it] + tav, GQW, iw);
		
		res += wttab [it] * iw * il;

//std::cout << it << "\t" << iw << std::endl;
	}
//exit (1);
		
	result = deltat * res;
}

double Cl_withoutRSD::Clg (int l, double Chi1av, double sigmaChi1, double Chi2av, double sigmaChi2, double d) {

	int Nmax = 200;
	std::complex<double>* coef = new std::complex<double> [Nmax+1];
	std::complex<double>* nutab = new std::complex<double> [Nmax+1];
	FFTlogPS (1e-8, 52., -1.9, Nmax, coef, nutab);
	//std::cout << cn [0] << "\t" << nun [0] << std::endl;

	int NN2 = Nmax / 2;

	WindowFunction wf (fname, GQW);
	WindowParameters WP1 ("Assassi");
	WP1.l = l;
	WP1.Chiav = Chi1av;
	WP1.sigmaChi = sigmaChi1;
	WindowParameters WP2 ("Assassi");
	WP2.l = l;
	WP2.Chiav = Chi2av;
	WP2.sigmaChi = sigmaChi2;

	std::complex<double> res (0., 0.);
	std::complex<double> mickey;

	for (int nn = 0; nn < NN2; nn++) {

//std::cout << nn << std::endl;

		Clgnu (nutab [nn] + d, WP1, WP2, wf, mickey);
		res += 2. * coef [nn] * mickey;
		/*
std::cout << nn << "\t";
//std::cout << coef [nn] << "\t";
std::cout << mickey << "\t";
std::cout << std::endl;
//std::cout << res << std::endl;
*/
	}

	Clgnu (nutab [NN2] + d, WP1, WP2, wf, mickey);
	res += coef [NN2] * mickey;

	/*
std::cout << NN2 << "\t";
//std::cout << coef [NN2] << "\t";
std::cout << mickey << "\t";
std::cout << std::endl;
//std::cout << res << std::endl;
*/

	return res.real () * 0.5 * M_1_PI * M_1_PI; 
	
}

Cl_withoutRSD::~Cl_withoutRSD () {

	delete [] ttab;
	delete [] wttab;

}

//**********************************************************
//**********************************************************
//**********************************************************
//**********************************************************
//**********************************************************
//**********************************************************
//**********************************************************
//**********************************************************
//**********************************************************
//**********************************************************

Cl_withRSD::Cl_withRSD (std::string pfile):PowerSpectrum (pfile) {

	if (pfile == "Assassi") {

		GLt = 50;
		GQW = GLQ50;

	}

	ttab = new double [GLt];
	wttab = new double [GLt];

	std::shared_ptr <FileUtils> file;
	file -> read_GLQ (GQW, ttab, wttab, GLt);
	std::shared_ptr <OtherUtils> other;
	other -> BubbleSort_2vec (ttab, wttab, GLt);

	fname = pfile;

}

void Cl_withRSD::Clgnu (std::complex<double> nu, WindowParameters& wp1, WindowParameters& wp2, WindowFunction& WF, std::complex<double>& result) {

	double Chi1min = wp1.Chiav - 5. * wp1.sigmaChi;
	double Chi1max = wp1.Chiav + 5. * wp1.sigmaChi;

	double Chi2min = wp2.Chiav - 5. * wp2.sigmaChi;
	double Chi2max = wp2.Chiav + 5. * wp2.sigmaChi;

	double tmin, tmax;
	if (wp1.l > 10) {
		std::shared_ptr <SpecialFunctions> sf;
		double zmin = sf -> tminC (wp1.l, nu);
		std::shared_ptr <ElementaryFunctions> ef;
		tmin = ef -> MaxC (Chi1min / Chi2max, zmin);
		tmax = ef -> MinC (Chi1max / Chi2min, 1.);
	}
	else {
		tmin = 0.;
		tmax = 1.;
	}

	double deltat = (tmax - tmin) * 0.5;
	double tav = (tmax + tmin) * 0.5;

	std::complex<double> res (0., 0.);

	std::shared_ptr <SpecialFunctions> sf;

	std::complex<double> il, iw;
	for (int it = 0; it < GLt; it++) {

		sf -> Il (wp1.l, nu + 1e-10, deltat * ttab[it] + tav, il);

		WF.intWgalaxy (nu, wp1, wp2, 1, 1, deltat * ttab [it] + tav, GQW, iw);
		
		res += wttab [it] * iw * il;

	}
		
	result = deltat * res;
}

double Cl_withRSD::Clg (int l, double Chi1av, double sigmaChi1, double Chi2av, double sigmaChi2, double d) {

	int Nmax = 200;
	std::complex<double>* coef = new std::complex<double> [Nmax+1];
	std::complex<double>* nutab = new std::complex<double> [Nmax+1];
	FFTlogPS (1e-8, 52., -1.9, Nmax, coef, nutab);

	int NN2 = Nmax / 2;

	WindowFunction wf (fname, GQW);
	WindowParameters WP1 ("Assassi");
	WP1.l = l;
	WP1.Chiav = Chi1av;
	WP1.sigmaChi = sigmaChi1;
	WindowParameters WP2 ("Assassi");
	WP2.l = l;
	WP2.Chiav = Chi2av;
	WP2.sigmaChi = sigmaChi2;

	std::complex<double> res (0., 0.);
	std::complex<double> mickey;

	for (int nn = 0; nn < NN2; nn++) {

		Clgnu (nutab [nn] + d, WP1, WP2, wf, mickey);
		res += 2. * coef [nn] * mickey;

	}

	Clgnu (nutab [NN2] + d, WP1, WP2, wf, mickey);
	res += coef [NN2] * mickey;

	return res.real () * 0.5 * M_1_PI * M_1_PI; 
	
}

Cl_withRSD::~Cl_withRSD () {

	delete [] ttab;
	delete [] wttab;

}
