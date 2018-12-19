#include "miscFunctions.h"


Transformation::Transformation(MatrixXd rot, Vector3d oset, double scale)
{
	RotationMatrix = rot;
	offset = oset;
	scaling = scale;
	RotationMatrix_Inv = rot.inverse();
};
Eigen::Vector3d Transformation::apply(Eigen::Vector3d data) //data is (1,3) [x,y,z]
{
	//cout << RotationMatrix << endl << data << endl << scaling << endl << offset << endl;
	return RotationMatrix * data*scaling + offset;
	//xx = Transform.R*Yl'*Transform.s + Transform.t;
};
Eigen::Vector3d Transformation::invert(Eigen::Vector3d data) //data is (1,3) [x,y,z]
{
	return RotationMatrix_Inv * data / scaling - RotationMatrix_Inv * offset / scaling;
	//Transform.R^-1*Xl'*Transform.s^-1 - Transform.R^-1*Transform.t*Transform.s^-1
};


MatrixXd readMatrix(const char *filename)
{//
	ifstream fin(filename);

	int nrows, ncols;
	fin >> nrows;
	fin >> ncols;
	//cout << nrows << " " << ncols << endl;
	MatrixXd X = MatrixXd::Zero(nrows, ncols);


	if (fin.is_open())
	{
		for (int row = 0; row < nrows; row++)
			for (int col = 0; col < ncols; col++)
			{

				double item = 0.0;
				fin >> item;
				X(row, col) = (double)item;
			}
		fin.close();
	}
	return X;
};
