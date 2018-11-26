#include "newmatrix.h"


Matrix4D::Matrix4D(int p1, int p2, int p3, int p4)
{
	size[0] = p1;
	size[1] = p2;
	size[2] = p3;
	size[3] = p4;

	order[0] = 0;
	order[1] = 1;
	order[2] = 2;
	order[3] = 3;
	/*m = new Eigen::MatrixXd *[p1];
	for (int i = 0; i < p1; i++)
	{
		m[i] = new Eigen::MatrixXd[p2];
	}*/
	vector< MatrixXd > temp(p2,MatrixXd(p3,p4));
	m = vector< vector<MatrixXd> >(p1, temp);
}

Matrix4D::Matrix4D()
{}
Matrix4D::~Matrix4D()
{
	//this->deleteMatrix();
}
void Matrix4D::deleteMatrix()
{
	/*for (int i = 0; i < getSize(0); i++)
	{
		delete [] m[i];
	}
	delete [] m;*/
	//cout << "adsf";
	for (int i = 0; i<this->getSize(0); i++)
	m[i].clear();
}
int Matrix4D::getSize(int i)
{
	return size[order[i]];
}
Matrix4D Matrix4D::operator+(const Matrix4D& b) {
	Matrix4D result;
	int idx[4];

	for (idx[0] = 0; idx[0] < this->getSize(0); idx[0]++)
		for (idx[1] = 0; idx[1] < this->getSize(1); idx[1]++)
			result.m[idx[0]][idx[1]] = this->m[idx[0]][idx[1]] + b.m[idx[0]][idx[1]];


	return result;
}
double& Matrix4D::operator() (unsigned p1, unsigned p2, unsigned p3, unsigned p4)
{
	int temp[4] = { p1,p2,p3,p4 };
	//cout << "wants this coef " << temp[0] << " " << temp[1] << " " << temp[2] << " " << temp[3] << endl << "the stored sizes are " << getSize(0) << " " << getSize(1) << " " << getSize(2) << " " << getSize(3) << endl << "the order is " << order[0] << " " << order[1] << " " << order[2] << " " << order[3] << endl;
	//printf("that means we are doing m[%d][%d](%d, %d)\n", temp[order[0]], temp[order[1]], temp[order[2]], temp[order[3]]);
	return this->m[temp[order[0]]][temp[order[1]]](temp[order[2]], temp[order[3]]);
}
MatrixXd& Matrix4D::operator() (unsigned p1, unsigned p2)
{
	if (p1 >= getSize(0) || p2 >= getSize(1))
		throw std::out_of_range("tried access matrix out of range");
	return this->m[p1][p2];
}

Matrix4D Matrix4D::permute(int order[4]) {
	
	int idx[4];

	//makeTimer(decl);
	//Matrix4D *newmatrix = new Matrix4D(this->size[order[0]], this->size[order[1]], this->size[order[2]], this->size[order[3]]);
	Matrix4D newmatrix = Matrix4D(this->size[order[0]], this->size[order[1]], this->size[order[2]], this->size[order[3]]);
	//decl.stop();

	//cout << this->m[0][0] << endl;
	//transfers the data to the new matrix with the new indexes
	//makeTimer(set);
	//MatrixXd temp(this->size[order[2]], this->size[order[3]]);
	for (idx[0] = 0; idx[0] < this->getSize(0); idx[0]++)
	{
		for (idx[1] = 0; idx[1] < this->getSize(1); idx[1]++)
		{
			for (idx[2] = 0; idx[2] < this->getSize(2); idx[2]++)
			{
				for (idx[3] = 0; idx[3] < this->getSize(3); idx[3]++)
				{
					//cout << idx[order[2]] << ' ' << idx[order[3]] << ' ' << temp.rows() << ' ' << temp.cols() << endl;
					//temp(idx[order[2]], idx[order[3]]) = this->m[idx[0]][idx[1]](idx[2], idx[3]);
					//(*newmatrix)(idx[order[0]], idx[order[1]], idx[order[2]], idx[order[3]]) = this->m[idx[0]][idx[1]](idx[2], idx[3]);
					newmatrix(idx[order[0]], idx[order[1]], idx[order[2]], idx[order[3]]) = this->m[idx[0]][idx[1]](idx[2], idx[3]);

				}
			}
			//(*newmatrix).m[idx[order[0]]][idx[order[1]]] = temp;
		}

	}
	//	set.stop();
	return newmatrix;

}

Matrix4D Matrix4D::repmat(int order[4]) {
	int idx[4];
	int newsizes[4] = { order[0] * this->getSize(0), order[1] * this->getSize(1), order[2] * this->getSize(2), order[3] * this->getSize(3) };
	Matrix4D newmatrix = Matrix4D(newsizes[0], newsizes[1], newsizes[2], newsizes[3]);
	
	for (idx[0] = 0; idx[0] < newmatrix.getSize(0); idx[0]++)
	{
		for (idx[1] = 0; idx[1] < newmatrix.getSize(1); idx[1]++)
		{
			for (idx[2] = 0; idx[2] < newmatrix.getSize(2); idx[2]++)
				for (idx[3] = 0; idx[3] < newmatrix.getSize(3); idx[3]++)
				{

					newmatrix.m[idx[0]][idx[1]](idx[2], idx[3]) = this->m[idx[0] % this->getSize(0)][idx[1] % this->getSize(1)](idx[2] % this->getSize(2), idx[3] % this->getSize(3));

				}
		}
	}
	return newmatrix;
}

Matrix4D repmat(Eigen::MatrixXd matrix, int order[4]) {
	int idx[4];
	int newsizes[4] = { order[0] * matrix.rows(), order[1] * matrix.cols(), order[2], order[3] };
	Matrix4D newmatrix = Matrix4D(newsizes[0], newsizes[1], newsizes[2], newsizes[3]);

	//cout << "a = " << newmatrix->getSize(0) << " b = " << newmatrix->getSize(1) << " c = " << newmatrix->getSize(2) << " d = " << newmatrix->getSize(3) << endl;
	//cout << matrix << endl;
	for (idx[0] = 0; idx[0] < newmatrix.getSize(0); idx[0]++)
	{
		for (idx[1] = 0; idx[1] < newmatrix.getSize(1); idx[1]++)
		{
			for (idx[2] = 0; idx[2] < newmatrix.getSize(2); idx[2]++)
				for (idx[3] = 0; idx[3] < newmatrix.getSize(3); idx[3]++)
				{
					
					newmatrix.m[idx[0]][idx[1]](idx[2], idx[3]) = matrix(idx[0] % matrix.rows(), idx[1] % matrix.cols());
				}
		}
	}
	//cout << "newmatrix " << endl << newmatrix->m[0][0] << endl;
	return newmatrix;
}

Matrix4D repmat(Matrix4D matrix, int order[4]) {
	int idx[4];
	int newsizes[4] = { order[0] * matrix.getSize(0), order[1] * matrix.getSize(1), order[2] * matrix.getSize(2), order[3] * matrix.getSize(3) };
	Matrix4D newmatrix = Matrix4D(newsizes[0], newsizes[1], newsizes[2], newsizes[3]);

	//cout << "a = " << newmatrix->getSize(0) << " b = " << newmatrix->getSize(1) << " c = " << newmatrix->getSize(2) << " d = " << newmatrix->getSize(3) << endl;
	//newmatrix->m[0][0](1, 2) = (matrix)(remainder(idx[0], matrix.rows()), remainder(idx[1], matrix.cols()));

	for (idx[0] = 0; idx[0] < newmatrix.getSize(0); idx[0]++)
	{
		for (idx[1] = 0; idx[1] < newmatrix.getSize(1); idx[1]++)
		{
			for (idx[2] = 0; idx[2] < newmatrix.getSize(2); idx[2]++)
				for (idx[3] = 0; idx[3] < newmatrix.getSize(3); idx[3]++)
				{

					newmatrix.m[idx[0]][idx[1]](idx[2], idx[3]) = matrix.m[idx[0] % matrix.getSize(0)][idx[1] % matrix.getSize(1)](idx[2] % matrix.getSize(2), idx[3] % matrix.getSize(3));
					
				}
		}
	}
	return newmatrix;
}

MatrixXd repmat(Eigen::MatrixXd matrix, double order[2]) {
	int idx[2];
	int newsizes[2] = { order[0] * matrix.rows(), order[1] * matrix.cols() };
	MatrixXd newmatrix = MatrixXd::Zero(newsizes[0], newsizes[1]);
	//cout << *newmatrix << endl;
	for (idx[0] = 0; idx[0] < newmatrix.rows(); idx[0]++)
	{
		for (idx[1] = 0; idx[1] < newmatrix.cols(); idx[1]++)
		{
			(newmatrix)(idx[0], idx[1]) = matrix(idx[0] % matrix.rows(), idx[1] % matrix.cols());
		}
	}
	return newmatrix;
}

Matrix4D mmat(Matrix4D* A, Matrix4D* B) { //mmat(A,B,[1 2]) default
	Matrix4D m3, m4;
	int temp[4] = { 2, 3 ,0 ,1 };
	m3 = (*A).permute(temp);
	m4 = (*B).permute(temp);

	/*auto future = std::async(&Matrix4D::permute, A, new int[4]{ 2, 3 ,0 ,1 });
	auto future2 = std::async(&Matrix4D::permute, B, new int[4]{ 2, 3 ,0 ,1 });

	m3 = future.get();
	m4 = future2.get();
	*/
	//timePermute.stop();
	//makeTimer(timeMult);

	//cout << m3.getSize(0) << endl << m3.getSize(1) << endl;
	//MatrixXd AA,BB;
	
	//AA.resize(m3(0, 0).rows()*m3.getSize(0)*m3.getSize(1), m3(0, 0).cols());
	//BB.resize(m4(0, 0).rows(), m4(0, 0).cols()*m4.getSize(1)*m4.getSize(0));
	//#pragma omp parallel for collapse(2) num_threads(1)
	for (int i = 0; i < m3.getSize(0); i++)
	{
		for (int j = 0; j < m3.getSize(1); j++)
		{
			//AA.block(i*m3(0, 0).rows(), 0, m3(0, 0).rows(), m3(0, 0).cols()) = m3(i, j);
			//BB.block(0, j*m3(0, 0).rows(), m3(0, 0).rows(), m3(0, 0).cols()) = m4(i, j);
			//cout << m3(i, j) << endl << m4(i, j) << endl;	
			//makeTimer(mult);
			m3(i,j) = m3(i, j) * m4(i, j); //matrix multiplication defined by EIGEN
			//mult.stop();
			/*MatrixXd newMatrix = MatrixXd::Zero(m3(i, j).rows(), m4(i, j).cols());
			for (int m = 0; m < m4(i,j).cols(); m++)
				for (int n = 0; n < m3(i, j).rows(); n++)
				{
					for (int k = 0; k < m4(i, j).rows(); k++)
						newMatrix(n, m) += m3(i, j)(n,k)*m4(i, j)(k,m);
				}
			m3(i, j) = newMatrix;*/
			
			
		}
	}
	
	//timeMult.stop();
	//cout << "mmat time = " << duration(timeNow() - t1) << "microseconds" << endl;
	//cout << "mmat() end " << endl;
	m3.size[2] = m3.m[0][0].rows();
	m3.size[3] = m3.m[0][0].cols();
	return m3.permute(temp);
}

MatrixXd Matrix4D::return2DMatrix(int dimensions[2], int index[2])
{
	//Q(2,:,:,2) = return2DMatrix({0,3},{1,1})
	//returns 2D matrix as if it were used in matlab as above
	int permArg[4] = { 0,0,0,0 };
	permArg[0] = dimensions[0];
	permArg[1] = dimensions[1];
	int j = 1;
	for (int i = 0; i < 4; i++)
	{
		if (i != dimensions[0] && i != dimensions[1])
			permArg[++j] = i;
	}
	return this->permute(permArg)(index[0], index[1]);
}

