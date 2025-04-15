#include "Timer.h"

void Timer::start() {
	startTime = std::chrono::steady_clock::now();
}

void Timer::end() {
	endTime = std::chrono::steady_clock::now();
}
