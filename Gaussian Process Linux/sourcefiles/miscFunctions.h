#pragma once
#include <Eigen/Dense>
#include<fstream>

#define MAXBUFSIZE  ((int) 1e6)
using namespace Eigen;
using namespace std;
MatrixXd readMatrix(const char *filename);

class Transformation
{
public:
	MatrixXd RotationMatrix, RotationMatrix_Inv;
	Vector3d offset;
	double scaling;
	Transformation(MatrixXd rot, Vector3d oset, double scale);
	Transformation() {};
	Eigen::Vector3d apply(Eigen::Vector3d data);//applies transformation
	Eigen::Vector3d invert(Eigen::Vector3d data);//makes the inverse transformation
};
