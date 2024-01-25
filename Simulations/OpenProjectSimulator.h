#pragma once

#include "Simulator.h"

#define SIMULATOR_UPDATES_PER_SECOND 100
#define FIXED_DT (1.0f / SIMULATOR_UPDATES_PER_SECOND)
#define USE_FIXED_DT true // :TimeStep

#define MAX_MASSPOINTS   32
#define MAX_SPRINGS      16
#define MAX_RIGID_BODIES 16

#define OFFSET_HEAT_GRID -1.f // Offset to heat grid for walls, ball, player rackets to the front
#define OFFSET_PLAYERACKETS 1.5f // Offset of player rackets to walls

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
	
	Mat4 transformation;
	Mat3 inverse_inertia;
	Mat3 inverse_I0;
	Real inverse_mass;
	
	Vec3 frame_force;
	Vec3 linear_velocity;
	Vec3 linear_factor; // Movement of this body will be multiplied by this factor, meaning we can restrict or increase movement along certain axis. @Cleanup: This does not seem to work apparently?

	Vec3 frame_torque;
	Vec3 angular_velocity; 
	Vec3 angular_momentum;
	Vec3 angular_factor; // See linear_factor
	
	Vec3 albedo;

    Real restitution;
	bool is_trigger;

	void create(Vec3 size, Real mass, Real restitution, bool is_trigger);
	void warp(Vec3 center, Quat orientation);
	void build_transformation_matrix();
    void build_inertia_tensor();
    
	void apply_force(Vec3 world_space_position, Vec3 force);
	void apply_impulse(Vec3 world_space_position, Vec3 impulse);
	void apply_angular_impulse(Vec3 impulse);
    void apply_torque(Vec3 torque);
    
	void set_linear_factor(Vec3 factor);
	void set_angular_factor(Vec3 factor);

	Vec3 get_world_space_velocity_at(Vec3 world_space_position);
};

// This records a collision between a trigger rigid body and another rigid body.
// It is only present the frame after the collision happened (since the collision
// detection happens after the game logic, but we want this info in the game logic).
struct Trigger_Collision {
	Rigid_Body *trigger;
	Rigid_Body *other;
};

struct Player_Racket
{
	Rigid_Body* platform;
	Spring* spring;
};

//
// Heat Diffusion
//

struct Heat_Grid {
	int width = 0, height = 0;
	int scale = 1; // scaling of each cell
	float *values = NULL; // Array of size 'width * height'. @@Leak: Does not get freed at program step, but eh.
	float heat_rise_by_ball = 1.0f;
	
	void create(int width, int height);
	void destroy();
	void reset();
	void apply_boundary_condition();
	void randomize();
	void set(int x, int y, float value);
	float get(int x, int y);
	float* get_cell_to_worldpos(float world_x, float world_y);
};

//
// Simulator
//

// We may wish to draw certain things for debugging, but not for the final experience, 
// e.g. the spring connections...
// Unfortunately the tweak bar UI works with indices, not enum values, so we cannot
// use cool bit-fiddeling here, and instead must do it like the poor peasants C++
// programmers are sometimes... Sadge.
enum DrawRequest {
	DRAW_SPRINGS      = 0x0,
	DRAW_HEAT_MAP     = 0x1,
	DRAW_RIGID_BODIES = 0x2,
	DRAW_EVERYTHING   = 0x3,
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
	int create_rigid_body(Vec3 size, Real mass, Real restitution, bool is_trigger);

	void apply_impulse_to_masspoint(int index, Vec3 impulse);
	void apply_impulse_to_rigid_body(int index, Vec3 world_space_position, Vec3 impulse);
    void apply_torque_to_rigid_body(int index, Vec3 torque);
	void warp_rigid_body(int index, Vec3 position, Quat orientation);

	void setup_demo_scene(); // Fuck you dennis.
	void setupHeatGrid();
	void setupWalls();
	void setupPlayerPlatforms();
	void setupBall();
	
	void set_default_camera_position();
	void setup_game();
	void game_logic(float dt);
	void update_game(float dt);
	void draw_game();

	void debug_print();

private:
    Rigid_Body * query_rigid_body(int index);
    Rigid_Body * create_and_query_rigid_body(Vec3 size, Real mass, Real restitution, bool is_trigger);

    Spring * query_spring(int index);
    Spring * create_and_query_spring(int a, int b, Real initial_length, Real stiffness);

    bool trigger_collision_occurred(Rigid_Body * trigger, Rigid_Body * other);
    
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
	double time_factor;
	bool running; // The UI can toggle this for debugging purposes. The game loop will not be executed when this is false.

	//
	// General stuff.
	//
	TwBar *tweak_bar;
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

	std::vector<Trigger_Collision> trigger_collisions;

	Rigid_Body* normal_walls[2];
	Rigid_Body* goals[2];
	Rigid_Body* ball;
	int score1;
	int score2;

	float heat_accelleration_for_ball;

	
	Player_Racket player_rackets[2];
    
	//
	// Heat Diffusion.
	//
	Heat_Grid heat_grid;
	float heat_alpha; // How fast temperate is diffused. Higher value means faster diffusion.
};
