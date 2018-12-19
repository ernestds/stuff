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
	int timeStep = 0;
	bool newData = false;
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

#define SENSOR_PERIOD (1.0 / 120.0)
#define FILTER_FREQ 20
void read_sensors() //this should read the sensors
{
	string s = "trajectory";
	bool periodOver = false;
	MatrixXd trajectoryVX = readMatrix(s.c_str()).row(1);
	MatrixXd trajectoryX = readMatrix(s.c_str()).row(0);
	string s1 = "trajectoryY";
	MatrixXd trajectoryVY = readMatrix(s1.c_str()).row(1);
	MatrixXd trajectoryY = readMatrix(s1.c_str()).row(0);
	string s2 = "trajectoryZ";
	MatrixXd trajectoryVZ = readMatrix(s2.c_str()).row(1);
	MatrixXd trajectoryZ = readMatrix(s2.c_str()).row(0);
	int i = 0;

	Vector3d lastPos, lastVel;
	lastPos << trajectoryX.col(i)(0), trajectoryY.col(i)(0), trajectoryZ.col(i)(0);
	lastVel << 0, 0, 0;

	dataNow.position(0) = trajectoryX.col(i)(0);
	dataNow.position(1) = trajectoryY.col(i)(0);
	dataNow.position(2) = trajectoryZ.col(i)(0);
	//cout << trajectoryZ << endl;
	while (running)
	{
		periodOver = false;
		auto t1 = std::chrono::high_resolution_clock::now();
		while (!periodOver)
		{
			auto t2 = std::chrono::high_resolution_clock::now();
			auto diff = t2 - t1;
			if (chrono::duration <double, milli>(diff).count() > SENSOR_PERIOD)
				periodOver = true;
		}


		dataNow.mtx.lock();
		if (dataNow.timeStep < 120)
		{
			dataNow.timeStep++;
			dataNow.newData = true;
			lastPos = dataNow.position;
			lastVel = dataNow.velocity;
		}

		dataNow.position(0) = trajectoryX.col(i)(0);
		dataNow.position(1) = trajectoryY.col(i)(0);
		dataNow.position(2) = trajectoryZ.col(i)(0);

		
		dataNow.velocity = (dataNow.position - lastPos) / SENSOR_PERIOD;
		Vector3d unfilt_acc = (dataNow.velocity - lastVel) / SENSOR_PERIOD;

		dataNow.acceleration = -dataNow.acceleration*(SENSOR_PERIOD* FILTER_FREQ - 1) + unfilt_acc * FILTER_FREQ * SENSOR_PERIOD;
		
		//cout << i << endl << dataNow.acceleration << endl << endl;
		//Yllfilt(i) = Yllonline(i)*w*Ts - Yllfilt(i-1)*(Ts*w-1);


		dataNow.mtx.unlock();
		if (i < (trajectoryVX.cols() - 1))
			i++;
	}
}
MatrixXd polynomialTraj(Vector3d iPos, Vector3d iVel, Vector3d iAcc, Vector3d fPos, double tf, int num_points)
{
	//makeTimer(aa);
	MatrixXd A(6, 6);
	VectorXd b(6);
	MatrixXd result(num_points, 3);
	double t = tf;
	A << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, t, t *t, t*t*t, t*t*t*t, t*t*t*t*t, 0, 1, 2 * t, 3 * t*t, 4 * t*t*t, 5 * t*t*t*t, 0, 0, 2, 6 * t, 12 * t*t, 20 * t*t*t;
	/*A <<
		1, 0, 0, 0, 0, 0,
		0, 1, 0, 0, 0, 0,
		0, 0, 2, 0, 0, 0,
		1, tf, tf *tf, tf*tf*tf, tf*tf*tf*tf, tf*tf*tf*tf*tf, 
		0, 1, 2 * tf, 3 * tf*tf, 4 * tf*tf*tf, 5 * tf*tf*tf*tf,
		0, 0, 2, 6 * tf, 12 * tf*tf, 20 * tf*tf*tf;*/

	MatrixXd temp(num_points, 6);
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < num_points; j++)
			temp(j, i) = pow(((double)(j + 1) * tf / (double)num_points), i);

	

	
	for (int n = 0; n < 3; n++) {
		//[trajectoryY(1,i); Ylonline(i); Yllonline(i); trajectoryY(1,50); 0; 0];
		b << iPos(n), iVel(n), iAcc(n), fPos(n), 0, 0;
		result.col(n) = temp * (A.inverse()*b);

	}
	//aa.stop();
	return result;
}
int main()
{
	sensorToDataset = Transformation((MatrixXd(3, 3) << 1, 0, 0, 0, 0, 1, 0, -1, 0).finished(), (Vector3d() << 0, 0, 0).finished(), 1);
	sensorToRobot = Transformation((MatrixXd(3,3) << 1, 0, 0, 0,0,1, 0,-1,0).finished(), (Vector3d() << 1, 1, 1).finished(), 1);
	std::thread t2 = thread(read_sensors);

	// A = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 2 0 0 0; 1 t t ^ 2 t ^ 3 t ^ 4 t ^ 5; 0 1 2 * t 3 * t ^ 2 4 * t ^ 3  5 * t ^ 4; 0 0 2 6 * t 12 * t ^ 2 20 * t ^ 3];
	// b = [trajectoryY(1, i); Ylonline(i); Yllonline(i); trajectoryY(1, 50); 0; 0];

	Vector3d pos, vel, acc, posf;
	pos << -0.2513, -0.1365, 0;
	vel << -0.1280, 0.1248, 0;
	acc << -2.9376, 1.1808, 0;
	posf << -0.4611, 0.0339, 0;
	polynomialTraj(pos, vel, acc, posf, 2, 21);
	polynomialTraj(pos, vel, acc, posf, 2, 22);
	cout << polynomialTraj(pos, vel, acc, posf, 2, 20) << endl;
	

	/*
	for i=1:6
    k=0;
    for j=2/100:2/100:2
        k=k+1;
        vec(k,i) = j^(i-1);
    end
end*/

	//cout << ans;
	while (true) {
		if (dataNow.newData)
		{
			//dataNow.newData = false;
			//cout << dataNow.velocity << endl;
		}
	};

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

	prediction_thread.join();
	running = false;
	prediction_thread.join();


	system("pause");
}
