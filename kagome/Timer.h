#pragma once

#include <chrono>
#include <iostream>

//simple wrapper to make std::chrono::stead_clock less tedious to use
class Timer
{	
	std::chrono::time_point<std::chrono::steady_clock> startTime, endTime;

public:
	typedef std::chrono::seconds s;
	typedef std::chrono::milliseconds ms;
	typedef std::chrono::microseconds us;

	void start();
	void end();

	template<class T>
	int read();
};

template<class T>
inline int Timer::read()
{
	return std::chrono::duration_cast<T>(endTime - startTime).count();
}
