#include "stdafx.h"

//  main.cpp
//  Multivariate Gaussian via Gibbs Sampling ----------------------------------------
//
//  Function skeletons created by Prof.RS, filled in by yantongz.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
#include "include.h"
#include "newmat.h"
#include "newmatio.h"   // necessary for printing out data structures in newmat library

#define PI 3.141592654

using namespace std;

// cf http://www.cplusplus.com/reference/random/uniform_real_distribution/operator()/
// If you want to set a seed -- do it only after debug phase is completed
// otherwise errors will not be repeatable.

// initially just use 
default_random_engine generator;

// unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
// default_random_engine generator(seed);

double get_gaussian(double mean, double standard_deviation)
{
	std::normal_distribution<double> distribution(mean, standard_deviation);
	double number = distribution(generator);
	return (number);
}

// this is just a hint -- if you are confused with what I suggest, you are welcome to take a
// different approach... no compulsion to use this thought-process, per se!
ColumnVector Gibbs_sampler_for_Multivariate_Gaussian(int index, ColumnVector Previous_value, SymmetricMatrix C, ColumnVector mean)
{
	// we want the conditional distribution of the "index"-th rv, conditioned on the
	// values of the other indices (i.e. "index"-th rv is random; all others are fixed).
	// See desription and linked YouTube Video from Prof. Niemi.

	// Took the formulae from http://fourier.eng.hmc.edu/e161/lectures/gaussianprocess/node7.html
	// Keep in mind that the expression for the correlation-matrix for the conditional probability
	// has a minor-typo at this reference (he has transposed "i" and "j"; you could have guessed it
	// by checking for dimensional consistency of the resulting correlation-matrix.

	// In the following -- I am assuming the original correlation matrix (Sigma) is partioned as
	// [Sigma11 Sigma12; Sigma21 Sigma22] (in MATLAB notation).  Sigma11 is an (n-1) x (n-1) matrix
	// Sigma12 is an (n-1)-Columnvector; Sigma21 is an (n-1)-Rowvector; Sigma22 is a scalar.

	// The import being -- we are keeping all the "top" (n-1)-many variables constant, and picking
	// a value for the bottom variable (i.e conditioning the bottom-variable on the previous (n-1)
	// variables.

	// Since the Gibbs-Sampler is not always going to sample the bottom-value, you need to go through
	// some clerical-work to make things work.

	// Write code to construct Sigma11 here

	int n = C.nrows();
	Matrix Sigma11(n-1, n-1);

	for (int i = 1; i <= n - 1; i++) {
		for (int j = 1; j <= n - 1; j++) {
			Sigma11(i, j) = C(i, j);
		}
	}

	// Write code to construct Sigma12 here
	ColumnVector Sigma12(n - 1);
	
	for (int k = 1; k <= n - 1; k++){
		Sigma12(k) = C(k, n);
	}

	// Write code to construct Sigma21 here
	RowVector Sigma21(n - 1);

	for (int l = 1; l <= n - 1; l++) {
		Sigma21(l) = C(n, l);
	}

	// Write code to construct Sigma22 here
	Real Sigma22 = C(n, n);

	// Write some account-keeping code here... in my "teliology" https://www.merriam-webster.com/dictionary/teleology
	// xx is a ColumnVector, which starts off with xx = Previous_value; and then the "index"-th variable is an
	// appropriately generated Univariate-Normal RV (called "x")

	ColumnVector xx(mean.nrows());
	xx = Previous_value;

	double new_mean, new_sd;
	ColumnVector x1(n - 1);
	ColumnVector mean1(n - 1);

	for (int ii = 1; ii <= n - 1; ii++) {
		x1(ii) = Previous_value(ii);
		mean1(ii) = mean(ii);
	}

	Matrix sig1 = Sigma21 * Sigma11.i()*(x1 - mean1);
	Real sigma1 = sig1.as_scalar();
	new_mean = sigma1;

	Matrix sig2 = Sigma22 - Sigma12.t()*Sigma11.i()*Sigma12;
	new_sd = sig2.as_scalar();

	xx(index) = get_gaussian(new_mean, new_sd);

	return xx;
}

double Theoretical_PDF(ColumnVector x, SymmetricMatrix C, ColumnVector mean)
{
	// write C++ code that computes the expression of equation 1 of the assignment description
	double fx, coef;

	// coef: base
	coef = 1 / pow(pow(2 * PI, mean.nrows())*C.determinant(), 1 / 2);

	// exponent
	Matrix expo = (x - mean).t()*C.i()*(x - mean);
	Real expon = expo.as_scalar();

	fx = exp(-0.5*expon)*coef;

	return fx;
}

int main(int argc, char* argv[])
{
	ColumnVector y_prev, y;
	Matrix count(100, 100);
	int no_of_trials, dimension;

	// 2D case
	dimension = 2;

	sscanf_s(argv[1], "%d", &no_of_trials);
	ofstream pdf_data(argv[2]);
	ofstream pdf_theory(argv[3]);

	// The covariance matrix
	SymmetricMatrix C(2);
	C(1, 1) = 0.75;
	C(1, 2) = 0.25;
	C(2, 1) = 0.25;
	C(2, 2) = 0.5;

	// The mean vector
	ColumnVector mean(2);
	mean(1) = 1.0;
	mean(2) = 2.0;

	cout << "Multivariate Gaussian Generator using Gibbs Sampling" << endl;
	cout << "Dimension = " << mean.nrows() << endl;
	cout << endl << "Mean Vector = " << endl;
	cout << mean;
	cout << endl << "Covariance Matrix = " << endl << C;

	for (int i = 1; i <= 100; i++)
		for (int j = 1; j <= 100; j++)
			count(i, j) = 0.0;

	y_prev = mean;
	for (int i = 0; i < no_of_trials; i++)
	{
		y = Gibbs_sampler_for_Multivariate_Gaussian(i % (mean.nrows()) + 1, y_prev, C, mean);
		for (int j = 1; j <= 100; j++) {
			for (int k = 1; k <= 100; k++) {
				if ((y(1) >= ((double)(j - 52) / 10)) && (y(1) < ((float)(j - 51) / 10)) &&
					(y(2) >= ((double)(k - 52) / 10)) && (y(2) < ((float)(k - 51) / 10))) {
					// count(j, k)++;
					cout << "entered loop 2" << endl; count(j, k)++;
					cout << j << ", " << k << ", " << count(j, k) << endl;
				}
			}
		}
		y_prev = y;
	}

	cout << count << endl;

	for (int j = 1; j <= 100; j++) {
		for (int k = 1; k <= 100; k++) {
			if (k < 100)
				pdf_data << count(j, k) / ((double)no_of_trials) << ", ";
			if (k == 100)
				pdf_data << count(j, k) / ((double)no_of_trials) << endl;
		}
	}

	double x1, x2;
	for (int j = 1; j <= 100; j++) {
		x1 = ((double)(j - 51) / 10);
		for (int k = 1; k <= 100; k++) {
			x2 = ((double)(k - 51) / 10);
			ColumnVector x(2);
			x(1) = x1;
			x(2) = x2;
			if (k < 100)
				pdf_theory << Theoretical_PDF(x, C, mean)*0.01 << ", ";
			if (k == 100)
				pdf_theory << Theoretical_PDF(x, C, mean)*0.01 << endl;
		}
	}
}
