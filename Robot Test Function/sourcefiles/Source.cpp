#include <iostream>

#include "newmatrix.h"
#include <algorithm> 
#include<fstream>
#include<stdio.h>
#include "Sonig.h"
#include "myTimer.h"
#include "miscFunctions.h"
#include "trajectoryGeneration.h"
#include <string>

#include <mutex>
#include <thread>
#include <atomic>

// ...

#define PRED_STEPS 30
using namespace std;
using namespace Eigen;


struct
{
	Eigen::Vector3d position = Vector3d::Zero(), velocity = Vector3d::Zero(), acceleration = Vector3d::Zero();
	bool newData;
	mutex mtx;
} dataNow{};

struct Predictions
{

	std::vector<Vector3d> mean, variance;
	bool newPred = false;
	mutex mtx;
	Predictions() : mean(PRED_STEPS, Vector3d::Zero()),
					variance(PRED_STEPS, Vector3d::Zero())
					{};
};

Predictions prediction;
Transformation sensorToDataset, sensorToRobot;
std::atomic_bool running{ true };


void make_prediction()
{
	
	MatrixXd Xhyp_lx = readMatrix("sonig/sonigX.hyp.lx");
	double Xhyp_ly = readMatrix("sonig/sonigX.hyp.ly")(0, 0);
	MatrixXd XXu = readMatrix("sonig/sonigX.Xu");
	MatrixXd Xfu_mean = readMatrix("sonig/Xfu_mean");
	MatrixXd Xfu_cov = readMatrix("sonig/Xfu_cov");
	MatrixXd XKuu = readMatrix("sonig/XKuu");
	//read values for the sonigV
	MatrixXd Vhyp_lx = readMatrix("sonig/sonigV.hyp.lx");
	double Vhyp_ly = readMatrix("sonig/sonigV.hyp.ly")(0, 0);
	MatrixXd VXu = readMatrix("sonig/sonigV.Xu");
	MatrixXd Vfu_mean = readMatrix("sonig/Vfu_mean");
	MatrixXd Vfu_cov = readMatrix("sonig/Vfu_cov");
	MatrixXd VKuu = readMatrix("sonig/VKuu");
	//read trajectories
	
	Sonig sonigX(Xhyp_lx, Xhyp_ly, XXu, Xfu_mean, Xfu_cov, XKuu);
	Sonig sonigV(Vhyp_lx, Vhyp_ly, VXu, Vfu_mean, Vfu_cov, VKuu);
	distribution xPDist, vPDist;
	MatrixXd vMean = MatrixXd::Zero(2, 1), vCov = MatrixXd::Zero(2, 2), xMean = MatrixXd::Zero(3, 1), xCov = MatrixXd::Zero(3, 3);

	vPDist.mean = MatrixXd::Zero(1, 1);
	vPDist.cov = MatrixXd::Zero(1, 1);

	xPDist.mean = MatrixXd::Zero(1, 1);
	xPDist.cov = MatrixXd::Zero(1, 1);

	MatrixXd predMean(PRED_STEPS, 3);
	MatrixXd predCov(PRED_STEPS, 3);
	
	while (true)
	{
		//
		//if (dataNow.newData)
		{
			makeTimer(predTime);

			dataNow.mtx.lock();
			xPDist.mean << dataNow.position(0);
			vPDist.mean << dataNow.velocity(0);
			dataNow.newData = false;
			dataNow.mtx.unlock();
			xPDist.cov << 0.0000000001;
			vPDist.cov << 0.0000000001;
			
			//later we have to convert actual time to timesteps here
			int calculatedTimestep = 0;
			for (int i = 0; i < (PRED_STEPS); i++)
			{
				
				vMean << calculatedTimestep + 1, vPDist.mean;
				vCov << 0.0000001, 0, 0, vPDist.cov;
				vPDist = sonigV.sonigStochasticPrediction(vMean, vCov);
				//vPDist = sonigV.sonigStochasticPrediction(vMean, vCov);
				//vPDist = sonigV.sonigStochasticPrediction(vMean, vCov);

				xMean << vPDist.mean, calculatedTimestep + 1, xPDist.mean;
				xCov << 0.03, 0, 0, 0, 0.0000001, 0, 0, 0, xPDist.cov;

				//xPDist = sonigX.sonigStochasticPrediction(xMean, xCov);
				//xPDist = sonigX.sonigStochasticPrediction(xMean, xCov);
				xPDist = sonigX.sonigStochasticPrediction(xMean, xCov);
				predMean(i, 0) = xPDist.mean(0);
				predMean(i, 1) = 0;
				predMean(i, 2) = 0;
				predCov(i, 0) = xPDist.cov(0);
				predCov(i, 1) = 0;
				predCov(i, 2) = 0; //can i apply the transformation to the variance just like that

				calculatedTimestep++;
			}
			
			// dataPred: each row is a time step, col are [x y z]
			prediction.mtx.lock();
			prediction.newPred = true;
			for (int i = 0; i < PRED_STEPS; i++) //applies the dataset -> sensor -> robot coordinate transformation to all predidictions
			{
				
				prediction.mean[i](0) = predMean(i, 0);
				prediction.mean[i](1) = predMean(i, 1);
				prediction.mean[i](2) = predMean(i, 2);
				prediction.variance[i](0) = predCov(i, 0);
				prediction.variance[i](1) = predCov(i, 1);
				prediction.variance[i](2) = predCov(i, 2);

			}
			prediction.mtx.unlock();
			predTime.stop();
			cout << "\n" << predMean(PRED_STEPS-1, 0) << "\n";
			//cout << prediction.mean[PRED_STEPS - 1](0) << endl; 
		}
	} //while(running)

}

void lala()
{
	while (true)
	{
		if (prediction.newPred);
		{
			//cout << "new prediction!" << endl;
			//prediction.newPred = false;
		}
	}
}
int main()
{

	sensorToDataset = Transformation(MatrixXd::Identity(3,3), Vector3d::Zero(),1);
	sensorToRobot = Transformation(MatrixXd::Identity(3, 3), Vector3d::Zero(), 1);
	string s = "trajectory";
	MatrixXd trajectoryV = readMatrix(s.c_str()).row(1);
	MatrixXd trajectoryX = readMatrix(s.c_str()).row(0);
	
	dataNow.mtx.lock();
	
	dataNow.position(0) = trajectoryX.col(0)(0);
	dataNow.velocity(0) = trajectoryV.col(0)(0);
	dataNow.newData = true;
	dataNow.mtx.unlock();

	running = true;
	cout << "before thread start " << endl;


	std::thread prediction_thread = thread(&make_prediction);
	std::thread t2= thread(lala);

	prediction_thread.join();
	running = false;
	prediction_thread.join();


	system("pause");
}
