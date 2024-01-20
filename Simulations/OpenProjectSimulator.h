#pragma once

#include "Simulator.h"

#define SIMULATOR_UPDATES_PER_SECOND 100
#define FIXED_DT (1.0f / SIMULATOR_UPDATES_PER_SECOND)
#define USE_FIXED_DT true // :TimeStep

#define MAX_MASSPOINTS 32
#define MAX_SPRINGS    16

#define DRAW_SPHERE_RADII 0.1f // For drawing masspoints

#define PRINT_FIXED_FLOAT "%2.05f"
#define PRINT_FIXED_VEC3  "{ " PRINT_FIXED_FLOAT ", " PRINT_FIXED_FLOAT ", " PRINT_FIXED_FLOAT " }"

//
// Utility Functions.
//

Vec3 lerp(Vec3 from, Vec3 to, float t);
float clamp_float(float value, float min, float max);
float random_float(float min, float max);
int random_int(int min, int max);

double get_current_time_in_milliseconds();

//
// Mass Spring System
//

struct Spring {
	int a, b; // Indices into the Simulator's masspoint array.
	float stiffness;
	float initial_length, current_length;
	Vec3 current_force_from_a_to_b;
};

struct Masspoint {
	float inverse_mass; // 0 means this masspoint is fixed in space.
	Vec3 position;
	Vec3 velocity;
};

//
// Rigid Body
//

//
// Heat Diffusion
//

//
// Simulator
//

enum DrawRequest { // We may wish to draw certain things for debugging, but not for the final experience, e.g. the spring connections...
	DRAW_NOTHING      = 0x0,
	DRAW_SPRINGS      = 0x1,
	DRAW_RIGID_BODIES = 0x2,
	DRAW_HEAT_MAP     = 0x4,
	DRAW_EVERYTHING   = DRAW_SPRINGS | DRAW_RIGID_BODIES | DRAW_HEAT_MAP,
};

class OpenProjectSimulator : public Simulator {
public:
	//
	// Constructor.
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

	//
	// Simulator API.
	//

	int create_masspoint(Vec3 position, float mass);
	int create_spring(int a, int b, float initial_length, float stiffness);

	void apply_impulse_to_masspoint(int index, Vec3 force);

	void setup_game();
	void update_game(float dt);
	void debug_print();

private:
	//
	// :TimeStep
	// Since this is supposed to be a game, we definitely want the simulation to run
	// frame-rate-independent. Everything else is just really dumb, since a faster computer
	// might change the game's difficulty. Also, the simulation sometimes runs slower after
	// startup and then has a sudden burst of speed (maybe some D3D11 stuff?), which is just
	// very anti-awesome.
	// Therefore, this provides it's own timing to ensure that we run a specific amount of
	// updates per second, with the appropriate delta for these updates.
	//
	double time_of_previous_update;


	//
	// General stuff.
	//
	DrawingUtilitiesClass * duc;
	DrawRequest draw_requests;
	float gravity;

	//
	// Spring System
	//
	Masspoint masspoints[MAX_MASSPOINTS];
	Spring springs[MAX_SPRINGS];
	int masspoint_count;
	int spring_count;
	float spring_damping; // 0 means no damping, 1 means complete damping
};