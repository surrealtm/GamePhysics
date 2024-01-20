#pragma once

#include "Simulator.h"

class OpenProjectSimulator : public Simulator {
public:
	//
	// Constructor
	//
	OpenProjectSimulator();

	//
	// Inherited functions.
	//
	const char *getTestCasesStr();
	void initUI(DrawingUtilitiesClass * DUC);
	void reset();
	void drawFrame(ID3D11DeviceContext * context);
	void notifyCaseChanged(int testCase);
	void simulateTimestep(float timestep);
	void externalForcesCalculations(float timeElapsed);
	void onClick(int x, int y);
	void onMouse(int x, int y);

private:
	
};