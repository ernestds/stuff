#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <Eigen/Dense>
#include "newmatrix.h"
#include "myTimer.h"
double Exp2(double x);
double Exp(double x);

class Sonig
{
public:
	
	Sonig(MatrixXd hyp_lx, double hyp_ly, MatrixXd Xu, MatrixXd fu_mean, MatrixXd fu_cov, MatrixXd Kuu);
	~Sonig();
	distribution sonigStochasticPrediction(MatrixXd xDist_mean, MatrixXd xDist_cov);
	MatrixXd A; // number of inputs, number of inputs
	MatrixXd diffNormalized;
	MatrixXd fu_mean;
	MatrixXd fu_cov, Kuu;
	MatrixXd xDist_cov; //number of inputs , number of inputs
	MatrixXd hyp_lx; //number of inputs , 1
	MatrixXd C;
	MatrixXd Kuu_inverse;
	MatrixXd fCovPart1;
	double hyp_ly;
	int nu; //number of induction points
	MatrixXd q_; //for Q(::::)
	MatrixXd E;
	MatrixXd Xu;// number of inputs, number of induction points
	MatrixXd diff;
	
	MatrixXd xDist_mean;//number of inputs, number of inputs
	MatrixXd Q_, D;

};

