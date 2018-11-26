#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <stdexcept>
#include <cmath>
#include "myTimer.h"
#include <mutex>
#include <thread>
#include <atomic>
#include <future>
using namespace Eigen;
using namespace std;

typedef std::chrono::high_resolution_clock::time_point TimeVar;

#define duration(a) std::chrono::duration_cast<std::chrono::microseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()
class Matrix4D
{
public:
	Eigen::MatrixXd **m;
	int size[4], order[4];
	Matrix4D(int p1, int p2, int p3, int p4);
	//Matrix4D(int p1, int p2, int p3, int p4, int o1, int o2, int o3, int o4);
	Matrix4D();

	void deleteMatrix();

	int getSize(int i);
	Matrix4D operator+(const Matrix4D& b);
	double& operator() (unsigned p1, unsigned p2, unsigned p3, unsigned p4);
	MatrixXd& operator() (unsigned p1, unsigned p2);
	Matrix4D permute(int order[4]);
	Matrix4D repmat(int order[4]);
	MatrixXd return2DMatrix(int dimensions[2], int index[2]);
};

Matrix4D mmat(Matrix4D* A, Matrix4D* B);
MatrixXd repmat(Eigen::MatrixXd matrix, double order[2]);
Matrix4D repmat(Matrix4D matrix, int order[4]);
Matrix4D repmat(Eigen::MatrixXd matrix, int order[4]);

struct distribution
{
	MatrixXd mean;
	MatrixXd cov;
};
