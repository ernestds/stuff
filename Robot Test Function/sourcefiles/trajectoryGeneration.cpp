#include "trajectoryGeneration.h"



Position::Position(double t, double x, double y, double z)
{
	mPositionTime = t;
	mLocation << x, y, z;
}
Position::Position() {}
const Position calculatePos(double positionTime, std::vector<Position> mPositionList) {
	Position calcPosition;
	calcPosition.mPositionTime = positionTime;
	// Build up the knots along the spline. The first row is the msecs offset from the
	// first stored position. The second, third and fourth rows are X, Y, and Z
	// respectively.
	MatrixXd points(4, mPositionList.size());
	for (int i = 0; i < mPositionList.size(); i++) {
		// The spline position will be milliseconds.
		points(0, i) = mPositionList.at(i).mPositionTime - mPositionList.front().mPositionTime; //time difference between i point and first point
		//https://doc.qt.io/qt-5/qdatetime.html#secsTo
		points(1, i) = mPositionList.at(i).mLocation(0); // X
		points(2, i) = mPositionList.at(i).mLocation(1); // Y
		points(3, i) = mPositionList.at(i).mLocation(2); // Z
	}
	// The degree of the interpolating spline needs to be one less than the number of points
	// that are fitted to the spline.
	const Spline4d spline = SplineFitting<Spline4d>::Interpolate(points, mPositionList.size() - 1);
	const Vector4d values = spline((positionTime - mPositionList.front().mPositionTime - points.row(0).minCoeff()) / (points.row(0).maxCoeff() - points.row(0).minCoeff()));
	//cout << ((positionTime - mPositionList.front().mPositionTime - points.row(0).minCoeff()) / (points.row(0).maxCoeff() - points.row(0).minCoeff())) << endl;

	// We already know the mPositionTime since that is how we
	// calculated the spline in the first place.
	calcPosition.mLocation(0) = values(1);
	calcPosition.mLocation(1) = values(2);
	calcPosition.mLocation(2) = values(3);

	return calcPosition;
}



MatrixXd generateTrajectory(MatrixXd intermediatePoints, int pointsNumber)
{
	//intermediatePoints include first and last points
	std::vector<Position> mPositionList(intermediatePoints.rows());
	for (int i = 0; i < intermediatePoints.rows(); i++)
	{
		mPositionList.at(i).mLocation(0) = intermediatePoints(i, 0);
		mPositionList.at(i).mLocation(1) = intermediatePoints(i, 1);
		mPositionList.at(i).mLocation(2) = intermediatePoints(i, 2);
		mPositionList.at(i).mPositionTime = (i) * 1 / (intermediatePoints.rows() - 1);
	}

	MatrixXd trajectoryPositions(pointsNumber, 3);

	for (int i = 0; i < pointsNumber; i++)
	{
		trajectoryPositions.row(i) = calculatePos((i + 1) * 1 / (double)pointsNumber, mPositionList).mLocation.transpose();
	}

	return trajectoryPositions;
};

