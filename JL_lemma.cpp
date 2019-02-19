#include "stdafx.h"

//  main.cpp
//  Johnson_Lindenstrauss Lemma
//
//  Function skeletons created by Prof.RS, filled in by yantongz.
//
#include <iostream>
#include <fstream>
#include <time.h>
#include <random>
#include <chrono>
#include <string.h>
#include <stdlib.h>
#include "include.h"
#include "newmat.h"
#include "newmatio.h"
#include "newmatap.h"

typedef std::pair<double, double> Box_Muller_Pair;

using namespace std;

// cf http://www.cplusplus.com/reference/random/uniform_real_distribution/operator()/
// If you want to set a seed -- do it only after debug phase is completed
// otherwise errors will not be repeatable.

unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator(seed);

// initially just use 
// default_random_engine generator;

// U.I.I.D. RV generator
double get_uniform()
{
	std::uniform_real_distribution <double> distribution(0.0, 1.0);
	double number = distribution(generator);
	return number;
}

// Using the Gaussian generator that is part of the C++ STL
double get_gaussian(double mean, double standard_deviation)
{
	std::normal_distribution<double> distribution(mean, standard_deviation);
	double number = distribution(generator);
	return number;
}

Matrix Generate_Random_Projection(int k, int d)
{
	// write code here that follows what we discussed in class

	Matrix jl(k, d); jl = 0.0;

	// d-dim gaussian RV's
	for (int ki = 1; ki <= k; ki++) {

		for (int di = 1; di <= d; di++) {

			double x = get_gaussian(0, 1 / sqrt(k));
			jl(ki, di) = x;
		}
	}

	return jl;
}


int main(int argc, const char * argv[])
{
	int original_dimension, reduced_dimension, no_of_cols, no_of_trials;
	double epsilon, delta, diff, diff2;
	clock_t time_before, time_after;
	clock_t time_before2, time_after2;


	sscanf_s(argv[1], "%d", &original_dimension);
	sscanf_s(argv[2], "%d", &reduced_dimension);
	sscanf_s(argv[3], "%d", &no_of_cols);
	sscanf_s(argv[4], "%lf", &epsilon);
	sscanf_s(argv[5], "%lf", &delta);
	ifstream input_file(argv[6]);
	sscanf_s(argv[7], "%d", &no_of_trials);

	cout << "Johnson-Lindenstrauss Lemma Demo" << endl;
	cout << "Reading a (" << original_dimension << " x " << no_of_cols << ") Matrix from file '";
	cout << argv[6] << "'" << endl;
	cout << "Reduced Dimension = " << reduced_dimension << endl;
	cout << "epsilon = " << epsilon << endl;
	cout << "delta = " << delta << endl;
	cout << "Reduced Dimension (i.e. " << reduced_dimension << ") should be >= " << ceil(1 / (epsilon*epsilon)*log(((double)no_of_cols) / delta));
	cout << " for the bound to hold with probability " << 1 - delta << endl;

	Matrix R = Generate_Random_Projection(reduced_dimension, original_dimension);

	Matrix A(original_dimension, no_of_cols);

	time_before = clock(); // recording time before we started to read data
	for (int j = 1; j <= original_dimension; j++)
	{
		for (int i = 1; i <= no_of_cols; i++)
		{
			double x;
			input_file >> x;
			A(j, i) = x;
		}
	}
	time_after = clock(); // recording time after testing is complete
	diff = ((double)time_after - (double)time_before);
	cout << "It took " << (diff / CLOCKS_PER_SEC) / 60.0 << " minutes to read data from file '" << argv[6] << "'" << endl;


	// testing Johnson-Lindenstrauss Lemma ------------------------------------------------

	int no_of_hits = 0;

	cout << "#Trails for the testing-phase = " << no_of_trials << endl;

	// this is the reduced-dimension representation of the x's (i.e the matrix of y's)
	// Matrix C(reduced_dimension, original_dimension);
	Matrix C;
	C = R * A;

	time_before2 = clock(); // recording time before the testing starts


// write code here for verification of JL-Lemma ------------------------------------------
	cout << "time_bf2 : " << time_before2 << endl;

	for (int ii = 1; ii <= no_of_trials; ii++) {
		
		cout << "trials : " << ii << endl;

		int rand1 = (rand() % no_of_cols) + 1;
		int rand2 = (rand() % no_of_cols) + 1;
		while (rand1 == rand2) {

			rand2 = (rand() % no_of_cols) + 1;
		}

		double diff_x = 0.000000;
		double xx = 0.000000;
    
		// pick a pair of vector x1 & x2 in dim d
		for (int jj = 1; jj <= original_dimension; jj++) {

			xx = A(jj, rand1) - A(jj, rand2);
			// cout << "xx : " << xx << endl;

			diff_x = diff_x + pow(xx, 2);
		}
		double diff_x_norm = sqrt(diff_x);
		cout << "diff_x_norm : " << diff_x_norm << endl;


		double diff_y = 0.000000;
		for (int kk = 1; kk <= reduced_dimension; kk++) {

			diff_y = diff_y + pow( C(kk, rand1)-C(kk, rand2), 2);
		}
		double diff_y_norm = sqrt(diff_y);
		cout << "diff_y_norm : " << diff_y_norm << endl;

		if ( ((1 - epsilon)*diff_x_norm <= diff_y_norm) && (diff_y_norm <= (1 + epsilon)*diff_x_norm)) {

			// cout << "in ifs" <<  endl;

			no_of_hits++;
		}
		

	} // end for

	time_after2 = clock(); // recording time after testing is complete

	diff2 = ((double)time_after2 - (double)time_before2);
	cout << "It took " << (diff2 / CLOCKS_PER_SEC) / 60.0 << " minutes for testing to be completed" << endl;

	cout << "Johnson-Lindenstrauss Lemma is satisfied " << no_of_hits << "-many times over ";
	cout << no_of_trials << " attempts" << endl;
	cout << "Empirical Probability = " << ((double)no_of_hits / no_of_trials);
	cout << "Theory says it should be at least:" << (1 - epsilon);

	return 0;
}
