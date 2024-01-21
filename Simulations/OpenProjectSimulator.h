#pragma once

#include "Simulator.h"

#define SIMULATOR_UPDATES_PER_SECOND 100
#define FIXED_DT (1.0f / SIMULATOR_UPDATES_PER_SECOND)
#define USE_FIXED_DT true // :TimeStep

#define MAX_MASSPOINTS   32
#define MAX_SPRINGS      16
#define MAX_RIGID_BODIES  8

#define PRINT_FIXED_FLOAT "%2.05f"
#define PRINT_FIXED_VEC3  "{ " PRINT_FIXED_FLOAT ", " PRINT_FIXED_FLOAT ", " PRINT_FIXED_FLOAT " }"

//
// Utility Functions.
//

struct Mat3 {
	Real _m[3][3];
};

Mat3 mat3_from_quaternion(Quat & q);
Mat3 mat3_mul_mat3(Mat3 & lhs, Mat3 & rhs);
Vec3 mat3_mul_vec3(Mat3 & lhs, Vec3 & rhs);
Mat3 mat3_tranpose(Mat3 & input);
Mat3 mat3_inverse(Mat3 & input);

Real turns_to_radians(Real value);
Quat quat_from_euler_angles(Real x, Real y, Real z);

Vec3 lerp(Vec3 from, Vec3 to, Real t);
float clamp_float(float value, float min, float max);
float random_float(float min, float max);
int random_int(int min, int max);

double get_current_time_in_milliseconds();


//
// Mass Spring System
//

struct Spring {
	int a, b; // Indices into the Simulator's masspoint array.
	Real stiffness;
	Real initial_length, current_length;
	Vec3 current_force_from_a_to_b;
};

struct Masspoint {
	Real inverse_mass; // 0 means this masspoint is fixed in space.
	Vec3 position;
	Vec3 velocity;
};

//
// Rigid Body
//

struct Rigid_Body {
	Vec3 center_of_mass;
	Quat orientation;
	Vec3 size;

	Vec3 force;
	Vec3 linear_impulse;
    Vec3 linear_velocity;
	
	Vec3 torque;
	Vec3 angular_impulse;
	Vec3 angular_velocity;
	Vec3 angular_momentum;
	
	Mat4 transformation;
	Mat3 inverse_inertia;

	Real inverse_mass;
	Mat3 inverse_I0;

	Vec3 albedo;

	void create(Vec3 size, Real mass);
	void warp(Vec3 center, Quat orientation);
	void build_transformation_matrix();
    void build_inertia_tensor();
    
	void apply_force(Vec3 world_space_position, Vec3 force);
	void apply_impulse(Vec3 world_space_position, Vec3 impulse);
    void apply_torque(Vec3 torque);
    
	Vec3 get_world_space_velocity_at(Vec3 world_space_position);
};

//
// Heat Diffusion
//

struct Heat_Grid {
	int width = 0, height = 0;
	float *values = NULL; // Array of size 'width * height'. @@Leak: Does not get freed at program step, but eh.

	void create(int width, int height);
	void destroy();
	void reset();
	void apply_boundary_condition();
	void randomize();
	void set(int x, int y, float value);
	float get(int x, int y);
};

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

	int create_masspoint(Vec3 position, Real mass);
	int create_spring(int a, int b, Real initial_length, Real stiffness);
	int create_rigid_body(Vec3 size, Real mass);

	void apply_impulse_to_masspoint(int index, Vec3 impulse);
	void apply_impulse_to_rigid_body(int index, Vec3 world_space_position, Vec3 impulse);
    void apply_torque_to_rigid_body(int index, Vec3 torque);
	void warp_rigid_body(int index, Vec3 position, Quat orientation);

	void setup_game();
	void update_game(float dt);
	void draw_game();

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

	//
	// Rigid Bodies.
	//
	Rigid_Body rigid_bodies[MAX_RIGID_BODIES];
	int rigid_body_count;

	//
	// Heat Diffusion.
	//
	Heat_Grid heat_grid;
	float heat_alpha; // How fast temperate is diffused. Higher value means faster diffusion.

    // @@Debug: For debugging the Rigid Body collision system.
    Vec3 contact_points_to_draw[4];
    int contact_points_to_draw_count;
};
