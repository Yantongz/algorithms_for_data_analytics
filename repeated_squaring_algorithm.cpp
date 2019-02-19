#include "stdafx.h"

// Function skeletons created by Prof.RS, filled in by yantongz.

//  main.cpp
//  Repeated Squaring
//

#include <iostream>
#include <fstream>
#include <time.h>
#include <chrono>
#include <string.h>
#include <stdio.h> 
#include <cmath>
#include <stdlib.h> 
#include <random>

#include "include.h"
#include "newmat.h"
#include "newmatio.h"
#include "newmatap.h"

using namespace std;

/*
double get_uniform()
{
	return((double)rand() / RAND_MAX);
}

double rand_en()
{	// low + (high-low)*get_uniform()
	return( 10.0*get_uniform() - 5);
}
*/

double ran_en()
{
	//  generate random numbers from A to B
	int a = -5;
	int b = 5;

	double entry = (double)((b - a)*((double)rand() / RAND_MAX) + a);
	return entry;
}

// main function
Matrix repeated_squaring(Matrix B, int exponent, int no_rows) {

	IdentityMatrix I(no_rows);
	if (exponent == 0) return I;

	else if (exponent % 2 == 1) {     // odd 

		return (B*repeated_squaring(B*B, ((exponent - 1) / 2), no_rows));
	}
	else {   // even
		return (repeated_squaring(B*B, (exponent / 2), no_rows));
	}
}

// brute force multiplication using Newmat???
Matrix brute_force_multiply(Matrix B, int exponent, int no_rows) {

	Matrix C(no_rows, no_rows);
	C = B;
	for (int i = 1; i < exponent; i++)
		C = B * C;
	return C;
}

int main(int argc, const char * argv[])
{
	/* initialize random seed: */
	srand((unsigned)time(NULL));  // it only works in main() 

	// You can't run code outside main(); You can only define variables

	int num_row, exponent;
	clock_t time_before1, time_after1;
	clock_t time_before2, time_after2;

	sscanf_s(argv[1], "%d", &exponent);
	sscanf_s(argv[2], "%d", &num_row);

	// only for plotting -----------------------------------------
	// ofstream outfile_1(argv[3]);
	// ofstream outfile_2(argv[4]);

	cout << "The number of rows/columns in the square matrix is: " << num_row << endl;
	cout << "The exponent is: " << exponent << endl;

	// create and fill up matrix b first
	Matrix b(num_row, num_row);

	for (int i = 1; i <= num_row; i++) {

		for (int j = 1; j <= num_row; j++) {

			b(i, j) = ran_en() ;
		}
	}

	cout << endl;
	cout << "Generating " << num_row << " x " << num_row << " random matrix is :" << endl;
	cout << b << endl;

	// repeated sq timing
	cout << endl;
	cout << "Repeated Squaring Result: " << endl;

	time_before1 = clock();

	Matrix sq = repeated_squaring(b, exponent, num_row);

	time_after1 = clock();

	double diff1 = difftime(time_after1, time_before1);

	// 
	cout << sq << endl;
	cout << "It took " << diff1 / CLOCKS_PER_SEC << " seconds to complete" << endl;

	// brute force multiply timing
	cout << endl;
	cout << "Direct Multiplication Result:" << endl;

	time_before2 = clock();

	Matrix bf = brute_force_multiply(b, exponent, num_row);

	time_after2 = clock();

	double diff2 = difftime(time_after2, time_before2);
	
	// 
	cout << bf << endl;
	cout << "It took " << diff2 / CLOCKS_PER_SEC << " seconds to complete" << endl;
	

	// save results for plot -----------------------------------------------
	/*
	clock_t bf1, bf2, af1, af2;

	for (int expo = 1; expo <= exponent; expo++) {

		bf1 = clock();
		Matrix s = repeated_squaring(b, expo, num_row);
		af1 = clock();
		double diff_s = ((double)af1 - (double)bf1);
		outfile_1 << diff_s/CLOCKS_PER_SEC << endl;

	}

	for (int expo = 1; expo <= exponent; expo++) {

		bf2 = clock();
		Matrix s = brute_force_multiply(b, expo, num_row);
		af2 = clock();
		double diff_b = ((double)af2 - (double)bf2);
		outfile_2 << diff_b / CLOCKS_PER_SEC << endl;

	}
	*/
}

