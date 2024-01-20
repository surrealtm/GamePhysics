#include "OpenProjectSimulator.h"

OpenProjectSimulator::OpenProjectSimulator() {

}

const char * OpenProjectSimulator::getTestCasesStr() {
	return "Open Project";
}

void OpenProjectSimulator::initUI(DrawingUtilitiesClass * DUC) {
	// @Incomplete: Maybe add a timestep variable.
}

void OpenProjectSimulator::reset() {
	// @Incomplete.
}

void OpenProjectSimulator::drawFrame(ID3D11DeviceContext * context) {
	// @Incomplete.
}

void OpenProjectSimulator::notifyCaseChanged(int testCase) {} // We don't have different test cases, so just ignore this.

void OpenProjectSimulator::simulateTimestep(float timeStep) {
	// @Incomplete
}

void OpenProjectSimulator::externalForcesCalculations(float timeStep) {
	// @Incomplete
}

void OpenProjectSimulator::onClick(int x, int y) {}

void OpenProjectSimulator::onMouse(int x, int y) {}