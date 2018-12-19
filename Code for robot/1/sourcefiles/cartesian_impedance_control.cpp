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
#define PRED_STEPS 10



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

std::vector<Vector3d> record_prediction(0,Eigen::Vector3d());
std::vector<Vector3d> record_gen_pos_ref(0,Eigen::Vector3d());
std::vector<Vector3d> record_covariance(0,Eigen::Vector3d());
std::vector<Vector3d> record_K_values(0,Eigen::Vector3d());
std::vector<Vector3d> record_position(0,Eigen::Vector3d());

double totalTime = 0;
int totalIter = 0;
void make_prediction()
{
	
	MatrixXd Xhyp_lx = readMatrix("sonig/sonigX.hyp.lx");
	double Xhyp_ly = readMatrix("sonig/sonigX.hyp.ly")(0, 0);
	MatrixXd XXu = readMatrix("sonig/sonigX.Xu");
	MatrixXd Xfu_mean = readMatrix("sonig/Xfu_mean");
	MatrixXd Xfu_cov = readMatrix("sonig/Xfu_cov");
	MatrixXd XKuu = readMatrix("sonig/XKuu");
	//read values for the sonigVX
	MatrixXd VXhyp_lx = readMatrix("sonig/sonigV.hyp.lx");
	double VXhyp_ly = readMatrix("sonig/sonigV.hyp.ly")(0, 0);
	MatrixXd VXXu = readMatrix("sonig/sonigV.Xu");
	MatrixXd VXfu_mean = readMatrix("sonig/Vfu_mean");
	MatrixXd VXfu_cov = readMatrix("sonig/Vfu_cov");
	MatrixXd VXKuu = readMatrix("sonig/VKuu");

	
	MatrixXd Yhyp_lx = readMatrix("sonig/Y/sonigX.hyp.lx");
	double Yhyp_ly = readMatrix("sonig/Y/sonigX.hyp.ly")(0, 0);
	MatrixXd YXu = readMatrix("sonig/Y/sonigX.Xu");
	MatrixXd Yfu_mean = readMatrix("sonig/Y/Xfu_mean");
	MatrixXd Yfu_cov = readMatrix("sonig/Y/Xfu_cov");
	MatrixXd YKuu = readMatrix("sonig/Y/XKuu");
	//read values for the sonigVX
	MatrixXd VYhyp_lx = readMatrix("sonig/Y/sonigV.hyp.lx");
	double VYhyp_ly = readMatrix("sonig/Y/sonigV.hyp.ly")(0, 0);
	MatrixXd VYXu = readMatrix("sonig/Y/sonigV.Xu");
	MatrixXd VYfu_mean = readMatrix("sonig/Y/Vfu_mean");
	MatrixXd VYfu_cov = readMatrix("sonig/Y/Vfu_cov");
	MatrixXd VYKuu = readMatrix("sonig/Y/VKuu");

	
	MatrixXd Zhyp_lx = readMatrix("sonig/Z/sonigX.hyp.lx");
	double Zhyp_ly = readMatrix("sonig/Z/sonigX.hyp.ly")(0, 0);
	MatrixXd ZXu = readMatrix("sonig/Z/sonigX.Xu");
	MatrixXd Zfu_mean = readMatrix("sonig/Z/Xfu_mean");
	MatrixXd Zfu_cov = readMatrix("sonig/Z/Xfu_cov");
	MatrixXd ZKuu = readMatrix("sonig/Z/XKuu");
	//read values for the sonigVX
	MatrixXd VZhyp_lx = readMatrix("sonig/Z/sonigV.hyp.lx");
	double VZhyp_ly = readMatrix("sonig/Z/sonigV.hyp.ly")(0, 0);
	MatrixXd VZXu = readMatrix("sonig/Z/sonigV.Xu");
	MatrixXd VZfu_mean = readMatrix("sonig/Z/Vfu_mean");
	MatrixXd VZfu_cov = readMatrix("sonig/Z/Vfu_cov");
	MatrixXd VZKuu = readMatrix("sonig/Z/VKuu");
	
	Sonig sonigX(Xhyp_lx, Xhyp_ly, XXu, Xfu_mean, Xfu_cov, XKuu);
	Sonig sonigVX(VXhyp_lx, VXhyp_ly, VXXu, VXfu_mean, VXfu_cov, VXKuu);
	
	Sonig sonigY(Yhyp_lx, Yhyp_ly, YXu, Yfu_mean, Yfu_cov, YKuu);
	Sonig sonigVY(VYhyp_lx, VYhyp_ly, VYXu, VYfu_mean, VYfu_cov, VYKuu);	
	
	Sonig sonigZ(Zhyp_lx, Zhyp_ly, ZXu, Zfu_mean, Zfu_cov, ZKuu);
	Sonig sonigVZ(VZhyp_lx, VZhyp_ly, VZXu, VZfu_mean, VZfu_cov, VZKuu);
	
	distribution xPDist, vxPDist, yPDist, vyPDist, zPDist, vzPDist;
	MatrixXd vMean = MatrixXd::Zero(2, 1), vCov = MatrixXd::Zero(2, 2), xMean = MatrixXd::Zero(3, 1), xCov = MatrixXd::Zero(3, 3);

	vxPDist.mean = MatrixXd::Zero(1, 1);
	vxPDist.cov = MatrixXd::Zero(1, 1);

	xPDist.mean = MatrixXd::Zero(1, 1);
	xPDist.cov = MatrixXd::Zero(1, 1);
	
	vyPDist.mean = MatrixXd::Zero(1, 1);
	vyPDist.cov = MatrixXd::Zero(1, 1);

	yPDist.mean = MatrixXd::Zero(1, 1);
	yPDist.cov = MatrixXd::Zero(1, 1);
	
	vzPDist.mean = MatrixXd::Zero(1, 1);
	vzPDist.cov = MatrixXd::Zero(1, 1);

	zPDist.mean = MatrixXd::Zero(1, 1);
	zPDist.cov = MatrixXd::Zero(1, 1);

	MatrixXd predMean(PRED_STEPS, 3);
	MatrixXd predCov(PRED_STEPS, 3);
	bool periodOver = false;
	while (running)
	{
		//
		//cout << "sfa" << endl;
		auto t3 = std::chrono::high_resolution_clock::now();
		if (dataNow.newData)
		{

			//cout << dataNow.newData << endl;
			auto t1 = std::chrono::high_resolution_clock::now();
			dataNow.mtx.lock();
			int calculatedTimestep = dataNow.timeStep;
			xPDist.mean << dataNow.position(0);
			vxPDist.mean << dataNow.velocity(0);
			yPDist.mean << dataNow.position(1);
			vyPDist.mean << dataNow.velocity(1);
			zPDist.mean << dataNow.position(2);
			vzPDist.mean << dataNow.velocity(2);
			dataNow.newData = false;
			dataNow.mtx.unlock();
			xPDist.cov << 0.0000000001;
			vxPDist.cov << 0.0000000001;

			yPDist.cov << 0.0000000001;
			vyPDist.cov << 0.0000000001;
						zPDist.cov << 0.0000000001;
			vzPDist.cov << 0.0000000001;
			
			//later we have to convert actual time to timesteps here
				
			for (int i = 0; i < (PRED_STEPS); i++)
			{
				
				vMean << calculatedTimestep + 1, vxPDist.mean;
				vCov << 0.0000001, 0, 0, vxPDist.cov;
				vxPDist = sonigVX.sonigStochasticPrediction(vMean, vCov);

				xMean << vxPDist.mean, calculatedTimestep + 1, xPDist.mean;
				xCov << 0.03, 0, 0, 0, 0.0000001, 0, 0, 0, xPDist.cov;
				xPDist = sonigX.sonigStochasticPrediction(xMean, xCov);
			

				vMean << calculatedTimestep + 1, vyPDist.mean;
				vCov << 0.0000001, 0, 0, vyPDist.cov;
				vyPDist = sonigVY.sonigStochasticPrediction(vMean, vCov);

				xMean << vyPDist.mean, calculatedTimestep + 1, yPDist.mean;
				xCov << 0.03, 0, 0, 0, 0.0000001, 0, 0, 0, yPDist.cov;
				yPDist = sonigY.sonigStochasticPrediction(xMean, xCov);
				
				vMean << calculatedTimestep + 1, vzPDist.mean;
				vCov << 0.0000001, 0, 0, vzPDist.cov;
				vzPDist = sonigVZ.sonigStochasticPrediction(vMean, vCov);

				xMean << vzPDist.mean, calculatedTimestep + 1, zPDist.mean;
				xCov << 0.03, 0, 0, 0, 0.0000001, 0, 0, 0, zPDist.cov;
				zPDist = sonigZ.sonigStochasticPrediction(xMean, xCov);
				
				predMean(i, 0) = xPDist.mean(0);
				predMean(i, 1) = yPDist.mean(0);
				predMean(i, 2) = zPDist.mean(0);
				predCov(i, 0) = xPDist.cov(0);
				predCov(i, 1) = yPDist.cov(0);
				predCov(i, 2) = zPDist.cov(0); //can i apply the transformation to the variance just like that
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
			
			auto t2 = std::chrono::high_resolution_clock::now();
			auto diff = t2 - t1;
			totalTime += chrono::duration <double, milli>(diff).count();
			totalIter++;
			//cout << endl << prediction.mean[0] << endl;
			//cout << prediction.mean[PRED_STEPS - 1](0) << endl;
			//cout << prediction.mean[PRED_STEPS-1] << endl;
		}
				periodOver = false;
		while (!periodOver)
		{
			auto t4 = std::chrono::high_resolution_clock::now();
			auto difff = t4 - t3;
			if (chrono::duration <double, milli>(difff).count() > 100.0)
				periodOver = true;
		}
	} //while(running)

}


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

	//cout << trajectoryZ << endl;
	while (running)
	{
		periodOver = false;
		auto t1 = std::chrono::high_resolution_clock::now();
		while (!periodOver)
		{
			auto t2 = std::chrono::high_resolution_clock::now();
			auto diff = t2 - t1;
			if (chrono::duration <double, milli>(diff).count() > 100.0)
				periodOver = true;
		}
		dataNow.mtx.lock();
		dataNow.position(0) = trajectoryX.col(i)(0);
		dataNow.velocity(0) = trajectoryVX.col(i)(0);
		dataNow.position(1) = trajectoryY.col(i)(0);
		dataNow.velocity(1) = trajectoryVY.col(i)(0);
		dataNow.position(2) = trajectoryZ.col(i)(0);
		dataNow.velocity(2) = trajectoryVZ.col(i)(0);
		dataNow.timeStep++;
		dataNow.newData = true;
		dataNow.mtx.unlock();
		if (i < (trajectoryVX.cols()-1))
			i++;
	}
}

int main(int argc, char** argv) {
	// Check whether the required arguments were passed
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " <robot-hostname>" << std::endl;
		return -1;
	}
	record_prediction.reserve(20000);
	record_gen_pos_ref.reserve(20000);
	record_covariance.reserve(20000);
	record_K_values.reserve(20000);
	record_position.reserve(20000);
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
	
	//NAO ESQUECER DAQUI AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAaa
	//AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
	//AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
	Eigen::Vector3d offset1 = Vector3d::Zero();
	offset1(0) = 1;
	offset1(2) = -0.4;
	sensorToDataset = Transformation(MatrixXd::Identity(3,3), Vector3d::Zero(),1);
	sensorToRobot = Transformation(MatrixXd::Identity(3, 3), offset1, 1);
	
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

		points.row(1)(0) = position_d_initial(0);
		points.row(1)(1) = position_d_initial(1);
		points.row(1)(2) = position_d_initial(2);


		MatrixXd trajectory(MIDDLE_POINTS_NUMBER, 3); //trajectory are all the intermediate points that the robot goes through when going to the reference position

		trajectory = generateTrajectory(points, MIDDLE_POINTS_NUMBER);
		position_d = trajectory.row(0);
		position_d = position_d_initial;




		
		std::function<franka::Torques(const franka::RobotState&, franka::Duration)>
			impedance_control_callback = [&](const franka::RobotState& robot_state,
				franka::Duration period/*duration*/) -> franka::Torques {
			// get state variables
				//cout << "Control loop " << endl;
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
			stiffness(0,0) = stiffness_original(0,0) * (1 - prediction.variance[PRED_STEPS-1](0)/0.0015);
			stiffness(1,1) = stiffness_original(1,1) * (1 - prediction.variance[PRED_STEPS-1](1)/0.0023);
			//stiffness(2,2) = stiffness_original(2,2) * (1 - prediction.variance[PRED_STEPS-1](2));

			//updates reference and creates new trajectory
			
			//updates the "target value" used to generate the error used to generate torque
			
			if (prediction.newPred)
			{
				Vector3d dist_to_pred = prediction.mean[PRED_STEPS-1] + offset1 - position;

				points.row(0)(0) = position(0);
				points.row(0)(1) = position(1);
				points.row(0)(2) = position(2);

				points.row(1)(0) = prediction.mean[PRED_STEPS-1](0);
				points.row(1)(1) = prediction.mean[PRED_STEPS-1](1);
				points.row(1)(2) = prediction.mean[PRED_STEPS-1](2);
				points.row(1) += offset1;
				//MatrixXd trajectory(MIDDLE_POINTS_NUMBER, 3); //trajectory are all the intermediate points that the robot goes through when going to the reference position
				trajectory = generateTrajectory(points, ceil(dist_to_pred.norm()*20));
				prediction.newPred = false;
				//cout << trajectory << endl << endl;
				//timeTraj = 0;
				//cout << endl << prediction.variance[PRED_STEPS-1] << endl;
				timeStep = 1;
			}
			
			if(timeTraj > 0.05)
			{
				timeStep++;
				if (timeStep > (MIDDLE_POINTS_NUMBER-1))
					timeStep = (MIDDLE_POINTS_NUMBER-1);
				
				Vector3d temp1 = trajectory.row(timeStep);
				//position_d(0) = trajectory(timeStep,0);
				//position_d(1) = trajectory(timeStep,1);
				//position_d(2) = trajectory(timeStep,2);	
				//position_d = sensorToRobot.apply(temp1);
				timeTraj = 0;
				//
				position_d = temp1;
				//cout << endl << temp1 << endl;
				
				
			}
			//cout << "\n" << position_d << "\n";
			//position_d = position_d_initial;
			
			//calculates velocity
			Vector3d x_dot = jacobian.topRows(3)*dq;

			//saves data for recording
			
			record_position.push_back(position);
			record_prediction.push_back(prediction.mean[PRED_STEPS-1]);
			record_covariance.push_back(prediction.variance[PRED_STEPS-1]);
			Vector3d tempK;
			tempK << stiffness(0,0),stiffness(1,1),stiffness(2,2);
			record_K_values.push_back(tempK);
			record_gen_pos_ref.push_back(position_d);
			
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
			if (timeRef > 7)
			{
				running = false;
				std::cout << std::endl << "Finished motion, shutting down example" << std::endl;
        		return franka::MotionFinished(tau_final);
			}
			return tau_final;
		};
		// start real-time control loop

		std::thread prediction_thread;
		std::thread sensor_thread;
						std::cout << "WARNING: Collision thresholds are set to high values. "
			<< "Make sure you have the user stop at hand!" << std::endl
			<< "After starting try to push the robot and see how it reacts." << std::endl
			<< "Press Enter to continue..." << std::endl;
		
			std::cin.ignore();
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
	std::ofstream output_file_position;

	output_file_position.open("record_position.txt");
	for (unsigned long i = 0; i < record_position.size(); i++)
	{
		for (unsigned int ii = 0; ii < 3; ii++)
			output_file_position << "\t" << record_position[i][ii];
		output_file_position << "\r\n";
	}
	output_file_position.close();

	std::ofstream output_file_covariance;
	output_file_covariance.open("record_covariance.txt");
	for (unsigned long i = 0; i < record_covariance.size(); i++)
	{
		for (unsigned int ii = 0; ii < 3; ii++)
			output_file_covariance << "\t" << record_covariance[i][ii];
		output_file_covariance << "\r\n";
	}
	output_file_covariance.close();

	std::ofstream output_file_prediction;
	output_file_prediction.open("record_prediction.txt");
	for (unsigned long i = 0; i < record_prediction.size(); i++)
	{
		for (unsigned int ii = 0; ii < 3; ii++)
			output_file_prediction << "\t" << record_prediction[i][ii];
		output_file_prediction << "\r\n";
	}
	output_file_prediction.close();

	std::ofstream output_file_K_values;
	output_file_K_values.open("record_K_values.txt");
	for (unsigned long i = 0; i < record_K_values.size(); i++)
	{
		for (unsigned int ii = 0; ii < 3; ii++)
			output_file_K_values << "\t" << record_K_values[i][ii];
		output_file_K_values << "\r\n";
	}
	output_file_K_values.close();

	std::ofstream output_file_reference;
	output_file_reference.open("record_reference.txt");
	for (unsigned long i = 0; i < record_gen_pos_ref.size(); i++)
	{
		for (unsigned int ii = 0; ii < 3; ii++)
			output_file_reference << "\t" << record_gen_pos_ref[i][ii];
		output_file_reference << "\r\n";
	}
	output_file_reference.close();

	return 0;
}
