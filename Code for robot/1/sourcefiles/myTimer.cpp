#include "myTimer.h"





myTimer::myTimer(string s)
	{
		name = s;
		t1 = std::chrono::high_resolution_clock::now();
	};
void myTimer::start()
	{
		t1 = std::chrono::high_resolution_clock::now();
	};
void myTimer::stop()
	{
		t2 =std::chrono::high_resolution_clock::now();
		auto diff = t2 - t1;
		cout << chrono::duration <double, milli>(diff).count() << " ms time total for timer " << name << endl;
};
