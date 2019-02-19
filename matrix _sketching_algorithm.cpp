
#include "stdafx.h"

//
//  main.cpp
//  Matrix Sketch
//
//  Function skeletons created by Prof.RS, filled in by yantongz.
//
//  Edo Liberty's paper is about using a streaming algorithm for
//  Matrix Sketching.  But we are going to implement/test this in
//  in a non-streaming context.  "A" is a large (n x m) data-matrix
//  where m >>> n. We want to find a smaller (n x k) "B" matrix, where
//  k = 2/epsilon, such that the Frobenius-Norm of (A*A.t() - B*B.t())
//  is <= epsilon * the Frobenius-Norm of (A*A.t())

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include "include.h"
#include "newmat.h"
#include "newmatio.h"
#include "newmatap.h"

using namespace std;

double Frobenious_Norm(Matrix Data)
{
	double f_norm = Data.NormFrobenius();
	return f_norm;
}

Matrix Matrix_Sketch(Matrix Data, double epsilon)
{
	// Edo Liberty's Matrix Sketch will have the same number of rows
	// as the original data; the #cols is ceil(2.0/epsilon)
	int cols_of_sketch = ceil(2.0 / epsilon);

	if (cols_of_sketch < Data.nrows())    // deep matrix, fine for Newmat SVD
	{
		Matrix Result(Data.nrows(), cols_of_sketch);

		DiagonalMatrix D;
		Matrix U, V;

		Result = Data.submatrix(1, Data.nrows(), 1, cols_of_sketch);

		SVD(Result, D, U, V);

		double zero_out = D.Minimum();
		double delta = zero_out * zero_out;

		IdentityMatrix I(cols_of_sketch);
		D = D * D - delta * I;

		for (int diag = 1; diag <= D.nrows(); diag++) {

			D(diag, diag) = pow(D(diag, diag), 0.5);
		}
		Result = U * D;
		
		for (int l = cols_of_sketch + 1; l <= Data.ncols(); l++) {

			DiagonalMatrix D;
			Matrix U, V;

			// update the last col-------------------------------
			for (int i = 1; i <= Data.nrows(); i++)
			{
				Result(i, cols_of_sketch) = Data(i, l);
			}
			SVD(Result, D, U, V);

			double zero_out = D.Minimum();
			double delta = zero_out * zero_out;

			IdentityMatrix I(cols_of_sketch);
			D = D * D - delta * I;

			for (int diag = 1; diag <= D.nrows(); diag++) {

				D(diag, diag) = pow(D(diag, diag), 0.5);
			}
			Result = U * D;

		}
		cout << "Deep" << endl;
		return Result;
	}
	else  // long matrix, need to do transpose for svd in Newmat -----------------------------
	{
		Matrix Result1(cols_of_sketch, Data.nrows());
		DiagonalMatrix D;
		Matrix U, V;

		Result1 = Data.submatrix(1, Data.nrows(), 1, cols_of_sketch);

		SVD(Result1.t(), D, U, V);  // svd decomp format ---------------------------

		double zero_out = D.Minimum();
		double delta = zero_out * zero_out;

		IdentityMatrix I(Data.nrows());
		D = D * D - delta * I;
		
		for (int diag = 1; diag <= D.nrows(); diag++) {

			D(diag, diag) = pow(D(diag, diag), 0.5);
		}
		Result1 = D * V.t();

		for (int l = cols_of_sketch+1; l <= Data.ncols(); l++) {

			DiagonalMatrix D;
			Matrix U, V;

			for (int i = 1; i <= Data.nrows(); i++)
			{
				Result1(Data.nrows(), i) = Data(i, l );
			}
			SVD(Result1, D, U, V);


			double zero_out = D.Minimum();
			double delta = zero_out * zero_out;

			IdentityMatrix I(Data.nrows());
			D = D * D - delta * I;

			for (int diag = 1; diag <= D.nrows(); diag++) {

				D(diag, diag) = pow(D(diag, diag), 0.5);
			}
			Result1 = D * V.t();

		}
		return Result1;
	}
}

int main(int argc, char* argv[])
{
	int dimension, no_of_data_points;
	double epsilon;

	sscanf_s(argv[1], "%d", &dimension);
	sscanf_s(argv[2], "%d", &no_of_data_points);
	sscanf_s(argv[3], "%lf", &epsilon);
	ifstream input_file(argv[4]);
	ofstream output_file(argv[5]);

	Matrix Data(dimension, no_of_data_points);

	cout << "Edo Liberty's Matrix Sketching Algorithm" << endl;
	cout << "----------------------------------------" << endl;
	cout << "Original Data-Matrix has " << dimension << "-rows & " << no_of_data_points << "-cols" << endl;
	cout << "Epsilon = " << epsilon << " (i.e. max. of " << 100 * epsilon << "% reduction of  Frobenius-Norm of the Sketch Matrix)" << endl;
	cout << "Input File = " << argv[4] << endl;

	// Read the Data
	for (int i = 1; i <= dimension; i++)
		for (int j = 1; j <= no_of_data_points; j++)
		{
			double x;
			input_file >> x;
			Data(i, j) = x;
		}

	// Compute the Frobenius-Norm of the original Data-Matrix
	double Data_Forbenius_Norm = Frobenious_Norm(Data*Data.t());
	cout << "Frobenius Norm of the (" << Data.nrows() << " x " << Data.ncols() << ") Data Matrix = ";
	cout << Data_Forbenius_Norm << endl;

	Matrix Sketch(dimension, min(dimension, (int)ceil(2 / epsilon)));
	Sketch = Matrix_Sketch(Data, epsilon);
	double Sketch_Forbenius_Norm = Frobenious_Norm(Sketch*Sketch.t());
	cout << "Frobenius Norm of the (" << Sketch.nrows() << " x " << Sketch.ncols() << ") Sketch Matrix = ";
	cout << Sketch_Forbenius_Norm << endl;
	cout << "Change in Frobenius-Norm between Sketch & Original  = ";
	cout << setprecision(3) << 100 * (Sketch_Forbenius_Norm - Data_Forbenius_Norm) / Data_Forbenius_Norm << "%" << endl;

	output_file << Sketch;
	cout << "File `" << argv[5] << "' contains a (" << Sketch.nrows() << " x " << Sketch.ncols();
	cout << ") Matrix-Sketch" << endl;


}

