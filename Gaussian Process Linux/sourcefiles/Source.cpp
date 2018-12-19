#include <iostream>

#include "newmatrix.h"
#include <algorithm> 
#include <string> 
#include<fstream>
#include<stdio.h>
#include "Sonig.h"
#include "myTimer.h"
#include "miscFunctions.h"
#include <string>
// ...

#define PRED_STEPS 30
using namespace std;

void curve_plot(MatrixXd pred, MatrixXd trueTraj)


{
	string command_filename;
	ofstream command_unit;
	string data_filename;
	ofstream data_unit;
	int i;
	string plot_filename;

	//  Write the data file.
	//cout << "asfdi";
	data_unit.open("misc/pred");
	for (i = 0; i < pred.rows(); i++)
	{
		data_unit << pred(i,1) << "  "
			<< pred(i,0) << "  " << pred(i,0) + 2*sqrt(pred(i,2)) << "  " << pred(i, 0) - 2 * sqrt(pred(i, 2)) << "\n";
	}
	data_unit.close();

	data_unit.open("misc/trueTraj");
	for (i = 0; i < trueTraj.rows(); i++)
	{
		data_unit << trueTraj(i, 1) << "  "
			<< trueTraj(i, 0) << "\n";
	}
	data_unit.close();

	//  Write the command file.
	
	command_filename = "misc/gnu_commands.txt";
	command_unit.open(command_filename.c_str());
	//command_unit << "set term png\n";
	//plot_filename = name + ".png";
	//command_unit << "set output \"" << plot_filename << "\"\n";
	command_unit << "set grid\n";
	//command_unit << "set style data lines\n";
	//command_unit << "unset key\n";
	//command_unit << "set xlabel '<---X--->'\n";
	//command_unit << "set ylabel '<---Y--->'\n";
	//command_unit << "set timestamp\n";
	command_unit << "plot \"" << "misc/pred" << "\" using 1:2 title \"Prediction Mean\" with lines lw 3, \\\n \"misc/trueTraj\"  using 1:2 title \"True Trajectory\" with lines lw 3, \\\n \"misc/pred\"  using 1:3 title \"2 Std. Variations\" with lines lw 3, \\\n \"misc/pred\"  using 1:4  title \"2 Std. Variations\" with lines lw 3";

	//command_unit << "quit\n";
	command_unit.close();


	string cmd1 = "gnuplot -persist " + command_filename;

	system(cmd1.c_str());
	return;
}

#define MAXBUFSIZE  ((int) 1e6)
template<class dataType>
class recorder
{
	
	std::vector<dataType> data;
public:
	void push_back(dataType a)
	{
		data.push_back(a);
	}
	void save(string name)
	{
		std::ofstream file(name);
		for (const auto &temp : data)
		{
			file << temp;
			file << "\n";
		}
		file.close();
	}
};
template<>
class recorder <Vector3d>
{
	std::vector<Vector3d> data;
public:
	void push_back(Vector3d a)
	{
		data.push_back(a);
	}
	void save(string name)
	{
		std::ofstream file(name);
		for (const auto &temp : data)
		{
			file << temp.transpose();
			file << "\n";
		}
		file.close();
	}
};
int main()
{
	recorder<Vector3d> oi;
	oi.push_back(Vector3d::Zero()); oi.push_back(Vector3d::Zero()); oi.push_back(Vector3d::Zero());
	oi.save("teste");
	while (true);
	//read values for the sonigX
	MatrixXd Xhyp_lx = readMatrix("sonig/sonigX.hyp.lx");
	//cout << "hyp_lx " << hyp_lx << endl;
	double Xhyp_ly = readMatrix("sonig/sonigX.hyp.ly")(0, 0);
	//cout << "hyp_ly " << hyp_ly << endl;
	MatrixXd XXu = readMatrix("sonig/sonigX.Xu");
	//cout << "Xu " << Xu << endl;
	MatrixXd Xfu_mean = readMatrix("sonig/Xfu_mean");
	//cout << "fu_mean " << fu_mean << endl;
	MatrixXd Xfu_cov = readMatrix("sonig/Xfu_cov");
	//cout << "fu_cov " << fu_cov << endl;
	MatrixXd XKuu = readMatrix("sonig/XKuu");
	//cout << "Kuu " << Kuu << endl;

	//read values for the sonigV
	MatrixXd Vhyp_lx = readMatrix("sonig/sonigV.hyp.lx");
	//cout << "hyp_lx " << hyp_lx << endl;
	double Vhyp_ly = readMatrix("sonig/sonigV.hyp.ly")(0, 0);
	//cout << "hyp_ly " << hyp_ly << endl;
	MatrixXd VXu = readMatrix("sonig/sonigV.Xu");
	//cout << "Xu " << Xu << endl;
	MatrixXd Vfu_mean = readMatrix("sonig/Vfu_mean");
	//cout << "fu_mean " << fu_mean << endl;
	MatrixXd Vfu_cov = readMatrix("sonig/Vfu_cov");
	//cout << "fu_cov " << fu_cov << endl;
	MatrixXd VKuu = readMatrix("sonig/VKuu");
	//cout << "Kuu " << Kuu << endl;
	//read trajectories
	string s = "trajectory";
	//cout << "Type the filename for the trajectory (if none entered assumes trajectory)" << endl;
	//getline(cin, s);
	MatrixXd trajectoryV = readMatrix(s.c_str()).row(1);
	MatrixXd trajectoryX = readMatrix(s.c_str()).row(0);
	Sonig sonigX(Xhyp_lx,Xhyp_ly,XXu, Xfu_mean, Xfu_cov,XKuu);
	Sonig sonigV(Vhyp_lx, Vhyp_ly, VXu, Vfu_mean, Vfu_cov, VKuu);
	distribution xPDist, vPDist;
	MatrixXd vMean = MatrixXd::Zero(2,1), vCov = MatrixXd::Zero(2, 2), xMean = MatrixXd::Zero(3, 1), xCov = MatrixXd::Zero(3, 3);
	
	vPDist.mean = MatrixXd::Zero(1, 1);
	vPDist.cov = MatrixXd::Zero(1, 1);
	vPDist.mean << trajectoryV(0, 1);
	vPDist.cov << 0.0000000001;

	xPDist.mean = MatrixXd::Zero(1, 1);
	xPDist.cov = MatrixXd::Zero(1, 1);
	xPDist.mean << trajectoryX(0, 1);
	xPDist.cov << 0.0000000001;
	
	MatrixXd dataPred(PRED_STEPS,3), dataTrue;

	dataTrue = trajectoryX.transpose();
	dataTrue.conservativeResize(dataTrue.rows(), 2);
	
	for (int i = 0; i < dataTrue.rows(); i++)
		dataTrue(i, 1) = i;

	int temp = dataTrue.rows();
	

	for (int iter = 0; iter < dataTrue.rows(); iter = iter + 15)
	{
		xPDist.mean << trajectoryX(0, iter);
		xPDist.cov << 0.0000000001;
		vPDist.mean << trajectoryV(0, iter);
		vPDist.cov << 0.0000000001;
			makeTimer(timeXPred);
			for (int i = iter; i < min(temp, iter + PRED_STEPS); i++)
			{
				//cout << "     Iteration " << i << " -------------------------- " << endl << endl;
				vMean << i + 1, vPDist.mean;
				//cout << "vMean " << endl << vMean << endl << endl;
				vCov << 0.0000001, 0, 0, vPDist.cov;
				//cout << "vCov " << endl << vCov << endl << endl;
				vPDist = sonigV.sonigStochasticPrediction(vMean, vCov);
				xMean << vPDist.mean, i + 1, xPDist.mean;
				//cout << "xMean " << endl << xMean << endl << endl;
				xCov << 0.03, 0, 0, 0, 0.0000001, 0, 0, 0, xPDist.cov;
				//cout << "xCov " << endl << xCov << endl << endl;

				
				//makeTimer(predX);
				xPDist = sonigX.sonigStochasticPrediction(xMean, xCov);
				//predX.stop();

				dataPred(i - iter, 0) = xPDist.mean(0);
				dataPred(i - iter, 1) = i;
				dataPred(i - iter, 2) = xPDist.cov(0);
			}
			timeXPred.stop();
			cout << dataPred(PRED_STEPS-1, 0) << endl;
		curve_plot(dataPred, dataTrue);
		cin.get();
	}
	system("pause");
}
