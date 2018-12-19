#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <unsupported/Eigen/Polynomials>
#include <unsupported/Eigen/Splines>
using namespace std;
using namespace Eigen;
class Position {
public:
	double mPositionTime; // Time of the position
	Eigen::Vector3d mLocation; // 3D Location
	Position(double t, double x, double y, double z);
	Position();
};

typedef Eigen::Spline<double, 4> Spline4d;

const Position calculatePos(double positionTime, std::vector<Position> mPositionList);
MatrixXd generateTrajectory(MatrixXd intermediatePoints, int pointsNumber);