// hw1code.cpp : Defines the entry point for the console application.
//

//  main.cpp
//  Generating Exponential and Normal RVs from an unknown discrete distribution
//
//  Function skeletons created by Prof.RS, filled in by yantongz.
//

#include "stdafx.h"

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <chrono>
#include <cmath>
#include <fstream>


using namespace std;

double z1, z2;
double prob1, prob2, prob3;

// cf http://www.cplusplus.com/reference/random/uniform_real_distribution/operator()/
// If you want to set a seed -- do it only after debug phase is completed
// otherwise errors will not be repeatable.

// initially just use 
// default_random_engine generator;

unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator(seed);

// I need this uniform i.i.d. r.v. generator to simulate the unfair-coin in
// simulate_coin_toss(double heads_probability)
double get_uniform()
{
	std::uniform_real_distribution <double> distribution(0.0, 1.0);
	double number = distribution(generator);
	return (number);
}

// to assign unpredictable-values for Prob(EVENT1), Prob(EVENT2) and Prob(EVENT3)
// with the proviso that Prob(EVENT1) + Prob(EVENT2) + Prob(EVENT3) = 1.
void assign_random_probabilities_to_three_events()
{
	// randomly assign values for prob1, prob2 and prob3 such that the 3-tuple
	// (prob1, prob2, prob3) is uniformly-distributed over the surface defined by
	// the equation prob1 + prob2 + prob3 = 1...
	// sec 1.4 of the Prob. & Stats review for ways of generating constrained RVs.

	// to generate uni dist random points on the hyperplane x_1+...+x_n = M,
	// make n-1 calls to a uniform random number generator in [0,M], sort them z_i
	// x_1 = z_1, x_2 = z_2 - z_1, x_n = M-z_(n-1).

	// define prob1-3 here

	double u1 = get_uniform();
	double u2 = get_uniform();
	z1 = min(u1, u2);
	z2 = max(u1, u2);

	prob1 = z1;
	prob2 = z2 - z1;
	prob3 = 1 - z2;

}

// Discrete Probability Distribution Generator
// Generates EVENT1 with probability "prob1", EVENT2 with prob "prob2", EVENT3 with prob "prob3"
int three_event_generator(double prob1, double prob2, double prob3)
{
	// returning an event (num) here?
	double prob = get_uniform();
	//cout << "prob "<<prob<<endl;
	if (prob <= z1) return 1;
	else if (prob > z1 && prob <= z2) return 2;
	else return 3;

}


// This procedure converts the three-discrete-event-generator (where
// Prob(EVENT1) = prob1, Prob(EVENT2) = prob2, Prob(EVENT3) = 1-prob1-prob2
// into a fair-coin where Prob(HEADS) = Prob(TAILS) = 1/2.
int simulate_fair_coin_toss_from_three_event_generator()
{
	// 3 tosses, 27 possibilities, rejects none-123
	// head : 123, 132, 321; tail: 213, 231, 312.

	int first_toss = three_event_generator(prob1, prob2, 1 - prob1 - prob2);
	int sec_toss = three_event_generator(prob1, prob2, 1 - prob1 - prob2);
	int thr_toss = three_event_generator(prob1, prob2, 1 - prob1 - prob2);

	if (first_toss == 1 && sec_toss == 2 && thr_toss == 3) {
		return 0;
	}

	else if (first_toss == 1 && sec_toss == 3 && thr_toss == 2) {
		return 0;
	}
	else if (first_toss == 2 && sec_toss == 1 && thr_toss == 3) {
		return 1;
	}

	else if (first_toss == 2 && sec_toss == 3 && thr_toss == 1) {
		return 1;
	}
	else if (first_toss == 3 && sec_toss == 1 && thr_toss == 2) {
		return 1;
	}

	else if (first_toss == 3 && sec_toss == 2 && thr_toss == 1) {
		return 0;
	}
	else {
		first_toss = 0; sec_toss = 0; thr_toss = 0;
		simulate_fair_coin_toss_from_three_event_generator();
	}

}

// This procedure gets 32 fair-tosses from (> 32 tosses) of the unfair-coin.
// For every "fair-heads" (resp. "fair-tails") we write down a "1" (resp. "0")
// and interpret the 32-bit unsigned integer as a random unsigned int between
// 0 and 2^32-1.
double get_uniform_from_three_event_generator(double prob1, double prob2, double prob3)
{
	// write code here

	int tosses;
	double sum = 0.0;
	for (int i = 0; i < 32; i++) {
		tosses = simulate_fair_coin_toss_from_three_event_generator();
		cout << "tosses " << tosses << endl;
		sum += (tosses * pow(2, i));
	}

	double uiid = sum / (pow(2, 32) - 1);
	return uiid;
}

// This generates i.i.d. rv.s that are exponentially distributed with rate lambda
// using an unfair-coin
double get_exp_from_unfair_coin(double lambda, double prob1, double prob2)
{
	cout << "uniform random variable value: " << get_uniform_from_three_event_generator(prob1, prob2, 1 - prob1 - prob2) << endl;
	return ((-1.0 / lambda)*log(get_uniform_from_three_event_generator(prob1, prob2, 1 - prob1 - prob2)));
}

// This generates i.i.d. rv.s that are unit-normal distributed using an unfair-coin
// (I am using the Box-Muller method here)
double get_gaussian_from_unfair_coin(double prob1, double prob2)
{
	return (sqrt(-2.0*log(get_uniform_from_three_event_generator(prob1, prob2, 1 - prob1 - prob2)))*cos(2 * 3.141592654*get_uniform_from_three_event_generator(prob1, prob2, 1 - prob1 - prob2)));
}


int main(int argc, char* argv[])
{

	double lambda;
	long no_of_trials;

	sscanf_s(argv[1], "%ld", &no_of_trials);
	sscanf_s(argv[2], "%lf", &lambda);
	ofstream outfile_1(argv[3]);
	ofstream outfile_2(argv[4]);

	assign_random_probabilities_to_three_events();

	cout << "Generating i.i.d. Unit-Normals and Exponentials from a 3-Event Generator" << endl;
	cout << "Probability of EVENT1                     = " << prob1 << endl;
	cout << "Probability of EVENT2                     = " << prob2 << endl;
	cout << "Probability of EVENT3                     = " << prob3 << endl;
	cout << "Number of trials                          = " << no_of_trials << endl;
	cout << "Output File Name for the Unit-Normal Data = " << argv[3] << endl;
	cout << "Output File Name for the Exponential Data = " << argv[4] << endl;
	cout << "Rate of Exponential RVs                   = " << lambda << endl;

	// i.i.d. exponential r.v. generation
	int count[100];
	for (int i = 0; i < 100; i++)
		count[i] = 0;

	for (int i = 0; i < no_of_trials; i++)
	{
		double y = get_exp_from_unfair_coin(lambda, prob1, prob2);

		// debug print all  ----------------------------------
		//cout << "y, " << y << endl;
		for (int j = 0; j < 100; j++)
			if (y < ((double)j / 10)) {
				count[j]++;
				//cout << j << "count[j], " << count[j] << endl;
			}
	}

	// Exponential CDF -- experimental vs. theory
	for (int j = 0; j < 100; j++)
		outfile_2 << ((double)j / 10) << ", " << ((double)count[j]) / ((double)no_of_trials) << ", " << (1.0 - exp(-1.0*lambda*((double)j / 10))) << endl;

	// i.i.d. unit-normal generation
	double cdf_gaussian[100];
	for (int i = 0; i < 100; i++)
	{
		count[i] = 0;
		cdf_gaussian[i] = 0.0;
	}

	for (int i = 0; i < no_of_trials; i++)
	{
		double y = get_gaussian_from_unfair_coin(prob1, prob2);
		for (int j = 0; j < 100; j++)
			if (y < ((float)(j - 50) / 10))
				count[j]++;
	}
	// Unit-Normal CDF -- experimental vs. theory
	for (int j = 1; j < 100; j++)
		cdf_gaussian[j] = cdf_gaussian[j - 1] +
		((1 / sqrt(2 * 3.141592654))*exp(-1.0*(((double)(j - 50) / 10)*
		((double)(j - 50) / 10)) / 2) / 10.0);

	for (int j = 0; j < 100; j++)
		outfile_1 << ((double)(j - 50) / 10) << ", " <<
		((double)count[j]) / ((double)no_of_trials) << ", " << cdf_gaussian[j]
		<< endl;

}

