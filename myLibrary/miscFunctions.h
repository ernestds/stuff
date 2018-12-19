#pragma once
#include <Eigen/Dense>
#include<fstream>

#define MAXBUFSIZE  ((int) 1e6)
//using namespace Eigen;
using namespace std;
Eigen::MatrixXd readMatrix(const char *filename);

class Transformation
{
public:
	Eigen::MatrixXd RotationMatrix, RotationMatrix_Inv;
	Eigen::Vector3d offset;
	double scaling;
	Transformation(Eigen::MatrixXd rot, Eigen::Vector3d oset, double scale);
	Transformation() {};
	Eigen::Vector3d apply(Eigen::Vector3d data);//applies transformation
	Eigen::Vector3d invert(Eigen::Vector3d data);//makes the inverse transformation
};

template<class dataType>
class recorder
{

	std::vector<dataType> data;
public:
	void push_back(dataType a);
	void save(string name);
};
template<>
class recorder <Vector3d>
{
	std::vector<Vector3d> data;
public:
	void push_back(Vector3d a);
	void save(string name);
};