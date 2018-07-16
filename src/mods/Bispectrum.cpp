/**
 *	@file		src/mods/Bispectrum.cpp
 *	@date		25/05/2018
 *	@author	ATroja
 */

#include "../include/Bispectrum.h"
#include "../include/FileUtils.h"
#include "../include/OtherUtils.h"

#include <memory>

Bispectrum::Bispectrum (std::string pfile):PowerSpectrum (pfile) {


	if (pfile == "Assassi") {

		GLt = 50;
		GQW = GLQ50;
	}

	ttab =  new double [GLt];
	wttab = new double [GLt];

	std::shared_ptr <FileUtils> file;
	file -> read_GLQ (GQW, ttab, wttab, GLt);
	std::shared_ptr <OtherUtils> other;
	other -> BubbleSort_2vec (ttab, wttab, GLt);

	fname = pfile;

}
	
void newcomplexpow (double x, std::complex<double> y, std::complex<double>& z) {

	double x2y = pow (x, y.real ());
	double lnx = log (x);

	z.real (x2y * cos (y.imag () * lnx));
	z.imag (x2y * sin (y.imag () * lnx));
} 

void Bispectrum::IntIlgalaxynu (std::complex<double> nu, double r, WindowParameters& wp, WindowFunction& WF, int n, std::complex<double>& result) {

	double Chimin = wp.Chiav - 7. * wp.sigmaChi;
	double Chimax = wp.Chiav + 7. * wp.sigmaChi;

	std::shared_ptr <SpecialFunctions> sf;
	double zmin = sf -> tminC (wp.l, nu);
	std::shared_ptr <ElementaryFunctions> ef;
	double tmin = ef -> MaxC (Chimin / r, zmin);
	double tmax = ef -> MinC (Chimax / r, 1./zmin);

	double deltat = (tmax - tmin) * 0.5;
	double tav = (tmax + tmin) * 0.5;

	std::complex <double> res (0., 0.);

	std::complex<double> il, mickey;
	double iw;
	for (int it = 0; it < GLt; it++) {

		newcomplexpow (deltat * ttab [it] + tav, -nu, mickey);

		sf -> Il (wp.l, nu + 1e-7, deltat * ttab [it] + tav, il);

		iw = WF.Wgalaxy_noRSD_1 (r * (deltat * ttab [it] + tav), wp, n + 2);

		res += wttab [it] * mickey * iw  * il;
	}

	newcomplexpow (r, 1.-nu, mickey);

	result = deltat * mickey *res;

}

double Bispectrum::Ilgalaxy (WindowParameters& wp, double r, int n, WindowFunction& wf) {

	int Nmax = 200;
	std::complex<double>* coef = new std::complex<double> [Nmax+1];
	std::complex<double>* nutab = new std::complex<double> [Nmax+1];
	FFTlogPS (1e-6, 52., -1.9, Nmax, coef, nutab);

	int NN2 = Nmax / 2;

	std::complex<double> res (0., 0.);
	std::complex<double> mickey;
	for (int nn = 0; nn < NN2; nn++) {

		IntIlgalaxynu (nutab [nn], r, wp, wf, n, mickey);
		res += 2. * coef [nn] * mickey;
	}

	IntIlgalaxynu (nutab [NN2], r, wp, wf, n, mickey);
	res += coef [NN2] * mickey;

	return res.real ();
}

double Bispectrum::Bgalaxy (int l1, int l2, int l3, double Chi1av, double sigmaChi1, double Chi2av, double sigmaChi2, double Chi3av, double sigmaChi3, int n1, int n2, int n3) {

	double rmin = Chi3av - 5. * sigmaChi3;
	double rmax = Chi3av + 5. * sigmaChi3;

	double deltar = (rmax - rmin) * 0.5;
	double rav = (rmax + rmin) * 0.5;

	WindowFunction wf (fname, GQW);
	WindowParameters wp1 ("Assassi");
	wp1.l = l1;
	wp1.Chiav = Chi1av;
	wp1.sigmaChi = sigmaChi1;
	WindowParameters wp2 ("Assassi");
	wp2.l = l2;
	wp2.Chiav = Chi2av;
	wp2.sigmaChi = sigmaChi2;
	WindowParameters wp3 ("Assassi");
	wp3.l = l3;
	wp3.Chiav = Chi3av;
	wp3.sigmaChi = sigmaChi3;
	

	std::complex<double> res (0., 0.);

	for (int i = 0; i < GLt; i++) {
		
		res += wttab [i]  * Ilgalaxy (wp1, deltar * ttab [i] + rav, n1, wf)
										* Ilgalaxy (wp2, deltar * ttab [i] + rav, n2, wf)
										* wf.Wgalaxy_noRSD_2 (deltar * ttab [i] + rav, wp3, n3);
	}

	/*
	if (n3 == 0) {
		for (int i = 0; i < GLt; i++) {
		
			res += wttab [i]  * Ilgalaxy (wp1, deltar * ttab [i] + rav, n1, wf)
											* Ilgalaxy (wp2, deltar * ttab [i] + rav, n2, wf)
											* wf.Wgalaxy_noRSD_20 (deltar * ttab [i] + rav, wp3);
		}
	}
	if (n3 == 1) {
		for (int i = 0; i < GLt; i++) {
		
			res += wttab [i] * Ilgalaxy (wp1, deltar * ttab [i] + rav, n1, wf)
											* Ilgalaxy (wp2, deltar * ttab [i] + rav, n2, wf)
											* wf.Wgalaxy_noRSD_21 (deltar * ttab [i] + rav, wp3);
		}
	}
	if (n3 == 2) {
		for (int i = 0; i < GLt; i++) {
		
			res += wttab [i] * Ilgalaxy (wp1, deltar * ttab [i] + rav, n1, wf)
											* Ilgalaxy (wp2, deltar * ttab [i] + rav, n2, wf)
											* wf.Wgalaxy_noRSD_22 (deltar * ttab [i] + rav, wp3);
		}
	}
	*/

	return deltar * res.real() * 0.25 * M_1_PI * M_1_PI * M_1_PI * M_1_PI;
}

Bispectrum::~Bispectrum () {

	delete [] ttab;
	delete [] wttab;
}
