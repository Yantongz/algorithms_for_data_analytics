#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "C:\Users\yantongz\Desktop\531\lp_solve_5.5\lp_lib.h"
#include "stdafx.h"
using namespace std;

int main (int argc, char* argv[]) 
{

	lprec *lp;
	REAL solution[2];
	
	// setting the problem up: 2 real variables
	lp = make_lp(0,2);
	
	// first constraint x1 + 2x2 <= 6 
	{
		double row[] = {0, 1, 2};
		add_constraint(lp, row, LE, 6);
	}
	// second constraint 4x1 + 3x2 <= 12
	{
		double row[] = {0, 4, 3};
		add_constraint(lp, row, LE, 12);
	}
	// objective function 7x1+5x2... since the regular lp_solve 
	// call minimizes the cost, and I am interested in maximizing
	// I have to flip the signs of the coeeficients... I have to 
	// flip the sign of the optimal value back to get the maximum
	// optimal value. 
	{
		double row[] = {0, -7, -5};
		set_obj_fn(lp, row);
	}
	
	// solve the lp
	solve(lp);
	
	// print optimal value
	cout << "Optimal Value is: " << -1*get_objective(lp) << endl;
	
	// print the optimizing values of the variables
	get_variables(lp, solution);
	cout << "Optimal solution:" << endl;
	for (int i = 0; i < 2; i++)
		cout << "x[" << i+1 << "] = " << solution[i] << endl;
	
	// delete the lp to release memory
	delete_lp(lp);
	
	return(0);
}