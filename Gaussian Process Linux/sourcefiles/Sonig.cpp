#include "Sonig.h"

double Exp2(double x)
{
	return std::exp(x / 2);
}
double Exp(double x)
{
	return std::exp(x);
}


Sonig::Sonig(MatrixXd hyp_lx, double hyp_ly, MatrixXd Xu, MatrixXd fu_mean, MatrixXd fu_cov, MatrixXd Kuu)
{
	int inputSize = hyp_lx.rows();
	this->nu = Xu.cols();
	cout << "nu "<< nu << endl;
	 // number of inputs, number of inputs
	diffNormalized = MatrixXd::Zero(3, 3);
	xDist_cov = MatrixXd::Zero(inputSize, inputSize); //number of inputs , number of inputs
	this->hyp_lx = hyp_lx;//number of inputs , 1
	this->hyp_ly = hyp_ly;
	this->fu_cov = fu_cov;
	this -> fu_mean = fu_mean;
	this->Kuu = Kuu;
	C = MatrixXd::Zero(inputSize, inputSize);
	q_ = MatrixXd::Zero(1, inputSize); //for Q(::::)
	E = MatrixXd::Zero(inputSize, inputSize);
	this->Xu = Xu;// number of inputs, number of induction points
	diff = MatrixXd::Zero(nu, nu);
	xDist_mean = MatrixXd::Zero(inputSize, 1);//number of inputs,1 
	Kuu_inverse = Kuu.inverse();
	fCovPart1 = Kuu.fullPivHouseholderQr().solve((Kuu - fu_cov)) * Kuu.inverse();
}


Sonig::~Sonig()
{
}

distribution Sonig::sonigStochasticPrediction(MatrixXd xDist_mean, MatrixXd xDist_cov)
{
	int debug = 0;
	double B;
	double *temp = new double[2]{ 1, (double)nu };
	//makeTimer(t1);
	diff = Xu - repmat(xDist_mean, temp);
	//t1.stop();
	//q(:,i) = sonig.hyp.ly(i)^2/sqrt(det(xDist.cov)*det(inv(xDist.cov) + diag(1./sonig.hyp.lx(:).^2)))*exp(-1/2*sum((diff'/(diag(sonig.hyp.lx(:).^2) + xDist.cov)).*diff', 2));
	
	//makeTimer(t2);
	A = MatrixXd::Zero(hyp_lx.rows(), hyp_lx.rows());
	C = MatrixXd::Zero(hyp_lx.rows(), hyp_lx.rows());

	A.diagonal() = hyp_lx.cwiseProduct(hyp_lx).cwiseInverse();


	A = A + xDist_cov.inverse();
	C.diagonal() = hyp_lx.cwiseProduct(hyp_lx);
	
	

	//  Insert the code that will be timed

	

	// Store the time difference between start and end
	
	q_ = 1 / sqrt(A.determinant()*xDist_cov.determinant()) * hyp_ly * hyp_ly * (((diff.transpose()*(xDist_cov + C).inverse()).cwiseProduct(diff.transpose())).rowwise().sum() * -0.5).unaryExpr(&Exp);
	//t2.stop();
	//makeTimer(t3);
	//xn
	//G[0] = xn = repmat(permute(diffNormalized(:, : , i), [1, 3, 2, 4]), [1, 1, 1, sonig.nu]) + repmat(permute(diffNormalized(:, : , j), [1, 3, 4, 2]), [1, 1, sonig.nu, 1]);
	diffNormalized = diff;
	for (int i = 0; i < hyp_lx.rows(); i++)
	{
		diffNormalized.row(i) = diff.row(i) / pow(hyp_lx(i), 2);
	}
	//t3.stop();
	//makeTimer(t4);

	G[0] = repmat(diffNormalized, new int[4]{ 1, 1,1,1 });
	G[1] = G[0];

	G[0] = G[0].permute(new int[4]{ 0,2,1,3 });
	G[1] = G[1].permute(new int[4]{ 0,2,3,1 });
	G[0] = G[0].repmat(new int[4]{ 1,1,1,(int)nu });
	G[1] = G[1].repmat(new int[4]{ 1,1,(int)nu,1 });
	//t4.stop();
	//makeTimer(t5);
	//G[0] = G[0] + G[1];
	int idx[4];
	//doesnt check for size
	for (idx[0] = 0; idx[0] < G[0].size[0]; idx[0]++)
		for (idx[1] = 0; idx[1] < G[0].size[1]; idx[1]++)
			G[0].m[idx[0]][idx[1]] = G[0].m[idx[0]][idx[1]] + G[1].m[idx[0]][idx[1]];
	//t5.stop();
	
	//Q(:,:,i,j)
	//Q(:,:,i,j) = sonig.hyp.ly(i)^2*sonig.hyp.ly(j)^2/sqrt(det(xDist.cov)*det(inv(xDist.cov) + diag(1./sonig.hyp.lx(:).^2) + diag(1./sonig.hyp.lx(:).^2)))*(exp(-1/2*sum(diffNormalized(:,:,i).*diff,1))'*exp(-1/2*sum(diffNormalized(:,:,j).*diff,1))).*permute(exp(1/2*mmat(mmat(permute(xn,[2,1,3,4]), repmat(inv(inv(xDist.cov) + diag(1./sonig.hyp.lx(:).^2) + diag(1./sonig.hyp.lx(:).^2)), [1,1,sonig.nu,sonig.nu])), xn)),[3,4,1,2]);
	//makeTimer(t6);
	A = MatrixXd::Zero(hyp_lx.rows(), hyp_lx.rows());
	A.diagonal() = hyp_lx.cwiseProduct(hyp_lx).cwiseInverse();
	
	A = 2 * A + xDist_cov.inverse();
	B = 1 / sqrt(A.determinant()*xDist_cov.determinant()) * pow(hyp_ly, 4); 
	C = diffNormalized.cwiseProduct(diff);
	
	D = -C.colwise().sum() * 1 / 2;
	D = D.unaryExpr(&Exp);
	E = D.transpose()*D; 
	F = repmat(A.inverse(), new int[4]{ 1,1,(int)nu,(int)nu }); 
	//t6.stop();
	//makeTimer(t7);
	Matrix4D H = G[0].permute(new int[4]{ 1, 0, 2, 3 }); //permute(xn,[2,1,3,4])

	
	H = mmat(&H, &F); 

	H = mmat(&H, &G[0]); 
	//t7.stop();
	
	//assuming this whole funciton is 100% optimized and only the time used for the previous mmat multiplications is significant (no permute and so on), this
	// would take about 22ms to run, compared to matlab where the whole function takes slightly less than 1ms to run
	// btw this is missing fCov(i,i) = max(fCov(i,i),1e-16); at the end of the function
	//makeTimer(t8);
	for (idx[0] = 0; idx[0] < H.size[0]; idx[0]++)
		for (idx[1] = 0; idx[1] < H.size[1]; idx[1]++)
			H.m[idx[0]][idx[1]] = H.m[idx[0]][idx[1]].unaryExpr(&Exp2); 
	
	Q_ = H.return2DMatrix(new int[2]{ 0,1 }, new int[2]{ 0,0 });

	Q_ = B * E.cwiseProduct(Q_); 
	//t8.stop();
	//makeTimer(t9);
	MatrixXd fMean = ((q_.transpose()*Kuu_inverse)*fu_mean);
	MatrixXd fCov = MatrixXd::Zero(1, 1);// ((hyp_ly * hyp_ly) - ((Kuu.fullPivHouseholderQr().solve((Kuu - fu_cov)) * Kuu_inverse * Q_)).trace());
	fCov << ((hyp_ly * hyp_ly) - ((fCovPart1 * Q_)).trace());
	//t9.stop();
	//makeTimer(t10);
	distribution Distribution;
	Distribution.mean = fMean;
	Distribution.cov = fCov + ((fu_mean.transpose()*Kuu_inverse*Q_*Kuu_inverse*fu_mean) - (fMean * fMean));
	//t10.stop();

	return Distribution;

	
}
