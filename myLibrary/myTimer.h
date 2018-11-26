#pragma once
#include <string>
#include <chrono>
#include <iostream>
using namespace std;
typedef std::chrono::high_resolution_clock::time_point TimeVar;
#define makeTimer(Variable) myTimer Variable((string)#Variable)

class myTimer
{
public:
	TimeVar t1, t2;
	string name;
	myTimer(string s);
	void start();
	void stop();
};

