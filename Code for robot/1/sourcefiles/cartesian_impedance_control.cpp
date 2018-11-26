// Copyright (c) 2017 Franka Emika GmbH
// Use of this source code is governed by the Apache-2.0 license, see LICENSE
#include <array>
#include <functional>
#include <franka/duration.h>
#include <franka/exception.h>
#include <franka/model.h>
#include <franka/robot.h>
#include <franka/rate_limiting.h>
#include "examples_common.h"

#include <iomanip>
#include <fstream>
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <algorithm> 
#include <iostream>
#include <time.h>
#include <mutex>
#include <thread>
#include <atomic>
#include <chrono>

#include <algorithm>
#include<stdio.h>

#include <string>

#include <Eigen/Dense>
#include <unsupported/Eigen/Polynomials>
#include <unsupported/Eigen/Splines>

#include "Sonig.h"
#include "myTimer.h"
#include "miscFunctions.h"
#include "trajectoryGeneration.h"

using namespace Eigen;
using namespace std;



#define MIDDLE_POINTS_NUMBER 5
#define PRED_STEPS 30



struct
{
	Eigen::Vector3d position = Vector3d::Zero(), velocity = Vector3d::Zero(), acceleration = Vector3d::Zero();
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
std::vector<Vector3d> record_xdot;
std::vector<Vector3d> record_x;

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
	
	while (running)
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
			//cout << prediction.mean[PRED_STEPS - 1](0) << endl; 
		}
	} //while(running)

}


void read_sensors() //this should read the sensors
{
	string s = "trajectory";
	bool periodOver = false;
	MatrixXd trajectoryV = readMatrix(s.c_str()).row(1);
	MatrixXd trajectoryX = readMatrix(s.c_str()).row(0);
	int i = 0;
	while (running)
	{
		periodOver = false;
		auto t1 = std::chrono::high_resolution_clock::now();
		while (!periodOver)
		{
			auto t2 = std::chrono::high_resolution_clock::now();
			auto diff = t2 - t1;
			if (chrono::duration <double, milli>(diff).count() < 1000.0)
				periodOver = true;
		}
		dataNow.mtx.lock();
		dataNow.position(0) = trajectoryX.col(i)(0);
		dataNow.velocity(0) = trajectoryV.col(i)(0);
		dataNow.newData = true;
		dataNow.mtx.unlock();
		if (i < (trajectoryV.cols()-1))
			i++;
	}
}

int main(int argc, char** argv) {
	// Check whether the required arguments were passed
cout << "teste";
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " <robot-hostname>" << std::endl;
		return -1;
	}
	// Compliance parameters
	const double translational_stiffness{ 150.0 };
	const double rotational_stiffness{ 10.0 };
	Eigen::MatrixXd stiffness(6, 6), damping(6, 6), stiffness_original;
	stiffness.setZero();
	//
	stiffness.topLeftCorner(3, 3) << 40, 0, 0, 0, 34, 0, 0, 0, 14;
	stiffness.bottomRightCorner(3, 3) << rotational_stiffness * Eigen::MatrixXd::Identity(3, 3);
	damping.setZero();
	damping.topLeftCorner(3, 3) << 4.177, 0, 0, 0, 3.799, 0, 0, 0.2373;
	damping.bottomRightCorner(3, 3) << 2.0 * sqrt(rotational_stiffness) * Eigen::MatrixXd::Identity(3, 3);

	stiffness_original = stiffness;
	
	sensorToDataset = Transformation(MatrixXd::Identity(3,3), Vector3d::Zero(),1);
	sensorToRobot = Transformation(MatrixXd::Identity(3, 3), Vector3d::Zero(), 1);
	
	try {
		// connect to robot
		franka::Robot robot(argv[1]);
		setDefaultBehavior(robot);
		// load the kinematics and dynamics model
		franka::Model model = robot.loadModel();
		franka::RobotState initial_state = robot.readOnce();
		// equilibrium point is the initial position
		Eigen::Affine3d initial_transform(Eigen::Matrix4d::Map(initial_state.O_T_EE.data()));
		Eigen::Vector3d position_d_initial(initial_transform.translation());
		Eigen::Vector3d position_d = position_d_initial;




		//creates a moving reference to simulate as if it were reading data
		Eigen::Vector3d endpoint[10];
		Eigen::MatrixXd points(2, 3);
		points.row(0)(0) = 0;
		points.row(0)(1) = -0.5;
		points.row(0)(2) = 0.3;

		points.row(1)(0) = 0.5;
		points.row(1)(1) = -0.5;
		points.row(1)(2) = 0.3;
		for (int i = 0; i < 10; i++)
		{
			endpoint[i](0) = generateTrajectory(points, 10)(i, 0);
			endpoint[i](1) = generateTrajectory(points, 10)(i, 1);
			endpoint[i](2) = generateTrajectory(points, 10)(i, 2);
		}


		Eigen::Quaterniond orientation_d(initial_transform.linear());
		// set collision behavior
		robot.setCollisionBehavior({ {100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0} },
			{ {100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0} },
			{ {100.0, 100.0, 100.0, 100.0, 100.0, 100.0} },
			{ {100.0, 100.0, 100.0, 100.0, 100.0, 100.0} });
		// define callback for the torque control loop

		
		double time = 0.0;
		double timeTraj = 0.0;
		double timeRef = 0.0;
		int timeStep = 0;
		int refIndex = 0;
		double uncertainty = 0; //maximum 100

		//generates initial trajectory
		points.row(0)(0) = position_d_initial(0);
		points.row(0)(1) = position_d_initial(1);
		points.row(0)(2) = position_d_initial(2);

		points.row(1)(0) = endpoint[refIndex](0);
		points.row(1)(1) = endpoint[refIndex](1);
		points.row(1)(2) = endpoint[refIndex](2);


		MatrixXd trajectory(MIDDLE_POINTS_NUMBER, 3); //trajectory are all the intermediate points that the robot goes through when going to the reference position

		trajectory = generateTrajectory(points, MIDDLE_POINTS_NUMBER);
		position_d = trajectory.row(0);
		position_d = position_d_initial;



		
		std::function<franka::Torques(const franka::RobotState&, franka::Duration)>
			impedance_control_callback = [&](const franka::RobotState& robot_state,
				franka::Duration period/*duration*/) -> franka::Torques {
			// get state variables
				cout << "Control loop " << endl;
			std::array<double, 7> coriolis_array = model.coriolis(robot_state);
			std::array<double, 42> jacobian_array =
				model.zeroJacobian(franka::Frame::kEndEffector, robot_state);
			// convert to Eigen
			Eigen::Map<const Eigen::Matrix<double, 7, 1> > coriolis(coriolis_array.data());
			Eigen::Map<const Eigen::Matrix<double, 6, 7> > jacobian(jacobian_array.data());
			Eigen::Map<const Eigen::Matrix<double, 7, 1> > q(robot_state.q.data());
			Eigen::Map<const Eigen::Matrix<double, 7, 1> > dq(robot_state.dq.data());
			Eigen::Affine3d transform(Eigen::Matrix4d::Map(robot_state.O_T_EE.data()));
			Eigen::Vector3d position(transform.translation());
			Eigen::Quaterniond orientation(transform.linear());
			// compute error to desired equilibrium pose
			// position error
			Eigen::Matrix<double, 6, 1> error;
			error.head(3) << position - position_d;


			//mine

			time += period.toSec();
			timeTraj += period.toSec();
			timeRef += period.toSec();

			//reduces stiffness according to uncertainty
			stiffness = stiffness_original * (100 - uncertainty) / 100;


			//updates reference and creates new trajectory
			
			//updates the "target value" used to generate the error used to generate torque
			/*
			if (timeTraj > 0.5 && timeStep < MIDDLE_POINTS_NUMBER)
			{
				//position_d = trajectory.row(timeStep);
				position_d = position_d_initial;
				timeStep++;
				timeTraj = 0;
			}*/

			position_d = position_d_initial;
			//calculates velocity
			Vector3d x_dot = jacobian.topRows(3)*dq;

			//saves data for recording
			
			record_xdot.push_back(x_dot);
			record_x.push_back(position);
			//end mine

			

			// orientation error
			// "difference" quaternion
			if (orientation_d.coeffs().dot(orientation.coeffs()) < 0.0) {
				orientation.coeffs() << -orientation.coeffs();
			}
			Eigen::Quaterniond error_quaternion(orientation * orientation_d.inverse());
			// convert to axis angle
			Eigen::AngleAxisd error_quaternion_angle_axis(error_quaternion);
			// compute "orientation error"
			error.tail(3) << error_quaternion_angle_axis.axis() * error_quaternion_angle_axis.angle();
			// compute control
			Eigen::VectorXd tau_task(7), tau_d(7);
			// Spring damper system with damping ratio=1

			//stiffness = stiffness * (100-uncertainty)/100;
			tau_task << jacobian.transpose() * (-stiffness * error - damping * (jacobian * dq));
			tau_d << tau_task + coriolis;
			std::array<double, 7> tau_d_array{};

			// print stuffff
			/*cout << "Position " << endl << position << endl;
			cout << "Velocity " << endl << x_dot << endl;
			cout << "Error " << endl << error.head(3) << endl;
			cout << "tau_task " << endl << tau_task << endl;*/

			Eigen::VectorXd::Map(&tau_d_array[0], 7) = tau_d;

			franka::Torques tau_final = franka::limitRate(franka::kMaxTorqueRate, tau_d_array, robot_state.tau_J_d);
			if (timeRef > 10)
			{
				running = false;
				std::cout << std::endl << "Finished motion, shutting down example" << std::endl;
        		return franka::MotionFinished(tau_final);
			}
			return tau_final;
		};
		// start real-time control loop
		std::cout << "WARNING: Collision thresholds are set to high values. "
			<< "Make sure you have the user stop at hand!" << std::endl
			<< "After starting try to push the robot and see how it reacts." << std::endl
			<< "Press Enter to continue..." << std::endl;
		std::cin.ignore();
		std::thread prediction_thread;
		std::thread sensor_thread;
		prediction_thread = thread(&make_prediction);
		sensor_thread = thread(&read_sensors);
		robot.control(impedance_control_callback);
		running = false;
		prediction_thread.join();
		sensor_thread.join();
      }
	catch (const franka::Exception& ex) {
		// print exception
		running = false;
		std::cout << ex.what() << std::endl;
	}

	cout << "Saving data..." << endl;
	//saves data
	std::ofstream output_file_xdot;
	output_file_xdot.open("data_x_dot_record.txt");
	for (unsigned long i = 0; i < record_xdot.size(); i++)
	{
		for (unsigned int ii = 0; ii < 3; ii++)
			output_file_xdot << "\t" << record_xdot[i][ii];
		output_file_xdot << "\r\n";
	}
	output_file_xdot.close();

	std::ofstream output_file_x;
	output_file_xdot.open("data_x_record.txt");
	for (unsigned long i = 0; i < record_x.size(); i++)
	{
		for (unsigned int ii = 0; ii < 3; ii++)
			output_file_x << "\t" << record_x[i][ii];
		output_file_x << "\r\n";
	}
	output_file_x.close();

	return 0;
}
