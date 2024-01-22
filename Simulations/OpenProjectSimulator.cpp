#include "OpenProjectSimulator.h"
#include "pcgsolver.h"
#include "sat.h"

//
// Utility Functions
//

double win32_performance_frequency;


Mat3 mat3_from_quaternion(Quat & q) {
	Mat3 result;

	result._m[0][0] = 1.0f - 2.0f * q.y * q.y - 2.0f * q.z * q.z;
	result._m[0][1] = 2.0f * q.x * q.y - 2.0f * q.z * q.w;
	result._m[0][2] = 2.0f * q.x * q.z + 2.0f * q.y * q.w;

    result._m[1][0] = 2.0f * q.x * q.y + 2.0f * q.z * q.w;
	result._m[1][1] = 1.0f - 2.0f * q.x * q.x - 2.0f * q.z * q.z;
	result._m[1][2] = 2.0f * q.y * q.z - 2.0f * q.x * q.w;

    result._m[2][0] = 2.0f * q.x * q.z - 2.0f * q.y * q.w;
	result._m[2][1] = 2.0f * q.y * q.z + 2.0f * q.x * q.w;
	result._m[2][2] = 1.0f - 2.0f * q.x * q.x - 2.0f * q.y * q.y;

	return result;
}

Mat3 mat3_mul_mat3(Mat3 & lhs, Mat3 & rhs) {
	Mat3 result;

	result._m[0][0] = lhs._m[0][0] * rhs._m[0][0] + lhs._m[0][1] * rhs._m[1][0] + lhs._m[0][2] * rhs._m[2][0];
	result._m[0][1] = lhs._m[0][0] * rhs._m[0][1] + lhs._m[0][1] * rhs._m[1][1] + lhs._m[0][2] * rhs._m[2][1];
	result._m[0][2] = lhs._m[0][0] * rhs._m[0][2] + lhs._m[0][1] * rhs._m[1][2] + lhs._m[0][2] * rhs._m[2][2];
	
	result._m[1][0] = lhs._m[1][0] * rhs._m[0][0] + lhs._m[1][1] * rhs._m[1][0] + lhs._m[1][2] * rhs._m[2][0];
	result._m[1][1] = lhs._m[1][0] * rhs._m[0][1] + lhs._m[1][1] * rhs._m[1][1] + lhs._m[1][2] * rhs._m[2][1];
	result._m[1][2] = lhs._m[1][0] * rhs._m[0][2] + lhs._m[1][1] * rhs._m[1][2] + lhs._m[1][2] * rhs._m[2][2];
	
	result._m[2][0] = lhs._m[2][0] * rhs._m[0][0] + lhs._m[2][1] * rhs._m[1][0] + lhs._m[2][2] * rhs._m[2][0];
	result._m[2][1] = lhs._m[2][0] * rhs._m[0][1] + lhs._m[2][1] * rhs._m[1][1] + lhs._m[2][2] * rhs._m[2][1];
	result._m[2][2] = lhs._m[2][0] * rhs._m[0][2] + lhs._m[2][1] * rhs._m[1][2] + lhs._m[2][2] * rhs._m[2][2];

	return result;
}

Vec3 mat3_mul_vec3(Mat3 & lhs, Vec3 & rhs) {
	Vec3 result;

	result.x = lhs._m[0][0] * rhs.x + lhs._m[0][1] * rhs.y + lhs._m[0][2] * rhs.z;
	result.y = lhs._m[1][0] * rhs.x + lhs._m[1][1] * rhs.y + lhs._m[1][2] * rhs.z;
	result.z = lhs._m[2][0] * rhs.x + lhs._m[2][1] * rhs.y + lhs._m[2][2] * rhs.z;

	return result;
}

Mat3 mat3_tranpose(Mat3 & input) {
	Mat3 result;

	result._m[0][0] = input._m[0][0];
	result._m[0][1] = input._m[1][0];
	result._m[0][2] = input._m[2][0];

	result._m[1][0] = input._m[0][1];
	result._m[1][1] = input._m[1][1];
	result._m[1][2] = input._m[2][1];

	result._m[2][0] = input._m[0][2];
	result._m[2][1] = input._m[1][2];
	result._m[2][2] = input._m[2][2];

	return result;
}

Mat3 mat3_inverse(Mat3 & input) {        
    Real determinant = static_cast<Real>(1) /
        (input._m[0][0] * (input._m[1][1] * input._m[2][2] - input._m[2][1] * input._m[1][2]) -
         input._m[1][0] * (input._m[0][1] * input._m[2][2] - input._m[2][1] * input._m[0][2]) +
         input._m[2][0] * (input._m[0][1] * input._m[1][2] - input._m[1][1] * input._m[0][2]));
        
    Mat3 inverse;
    inverse._m[0][0] = +(input._m[1][1] * input._m[2][2] - input._m[2][1] * input._m[1][2]) * determinant;
    inverse._m[1][0] = -(input._m[1][0] * input._m[2][2] - input._m[2][0] * input._m[1][2]) * determinant;
    inverse._m[2][0] = +(input._m[1][0] * input._m[2][1] - input._m[2][0] * input._m[1][1]) * determinant;
    inverse._m[0][1] = -(input._m[0][1] * input._m[2][2] - input._m[2][1] * input._m[0][2]) * determinant;
    inverse._m[1][1] = +(input._m[0][0] * input._m[2][2] - input._m[2][0] * input._m[0][2]) * determinant;
    inverse._m[2][1] = -(input._m[0][0] * input._m[2][1] - input._m[2][0] * input._m[0][1]) * determinant;
    inverse._m[0][2] = +(input._m[0][1] * input._m[1][2] - input._m[1][1] * input._m[0][2]) * determinant;
    inverse._m[1][2] = -(input._m[0][0] * input._m[1][2] - input._m[1][0] * input._m[0][2]) * determinant;
    inverse._m[2][2] = +(input._m[0][0] * input._m[1][1] - input._m[1][0] * input._m[0][1]) * determinant;
    return inverse;  
}


Real turns_to_radians(Real value) {
    return value * static_cast<Real>(2 * 3.1415926535);
}

Quat quat_from_euler_angles(Real x, Real y, Real z) {
    Real rx = turns_to_radians(x / 2);
    Real ry = turns_to_radians(y / 2);
    Real rz = turns_to_radians(z / 2);

    Real cx = cos(rx);
    Real cy = cos(ry);
    Real cz = cos(rz);

    Real sx = sin(rx);
    Real sy = sin(ry);
    Real sz = sin(rz);

    Quat result;
    result.x = sx * cy * cz - cx * sy * sz;
    result.y = cx * sy * cz + sx * cy * sz;
    result.z = cx * cy * sz - sx * sy * cz;
    result.w = cx * cy * cz - sx * sy * sz;
    return result;
}


Vec3 lerp(Vec3 from, Vec3 to, Real t) {
    return Vec3(from.x * (1 - t) + to.x * t,
                from.y * (1 - t) + to.y * t,
                from.z * (1 - t) + to.z * t);
}

float clamp_float(float value, float min, float max) {
    if(value < min) return min;
    if(value > max) return max;
    return value;
}

float random_float(float min, float max) {
    float zero_to_one = std::rand() / (float) RAND_MAX;
    return zero_to_one * (max - min) + min;
}

int random_int(int min, int max) {
    int normalized = std::rand() % (max - min);
    return normalized + min;
}

void setup_timing() {
    LARGE_INTEGER frequency;
    QueryPerformanceFrequency(&frequency);
    win32_performance_frequency = (double) frequency.QuadPart;
}

double get_current_time_in_milliseconds() {
    LARGE_INTEGER counter;
    QueryPerformanceCounter(&counter);
    return (double) counter.QuadPart / win32_performance_frequency;
}

SAT_Input sat_input(Rigid_Body & body) {
    SAT_Input input;
    input.center      = sat_vec3(body.center_of_mass.x, body.center_of_mass.y, body.center_of_mass.z);
    input.orientation = sat_quat(body.orientation.x, body.orientation.y, body.orientation.z, body.orientation.w);
    input.size        = sat_vec3(body.size.x / 2, body.size.y / 2, body.size.z / 2); // The SAT expects the half_dimensions, whereas the transformation matrix expects the full_dimensions
    input.velocity    = sat_vec3(body.linear_velocity.x, body.linear_velocity.y, body.linear_velocity.z);
    return input;
}


//
// Rigid Bodies.
//

void Rigid_Body::create(Vec3 size, Real mass) {
    assert(mass >= 0.f && "Invalid Mass for Rigid Body.");

	this->size = size;
	this->inverse_mass = mass > 0.0 ? 1.0 / mass : 0.0;

	if(mass > 0) {
        Real lx = this->size.x, ly = this->size.y, lz = this->size.z; // The size vector is already in full extends.
        
		this->inverse_I0 = { 12.0 * this->inverse_mass / (ly * ly + lz * lz), 0, 0,
					         0, 12.0 * this->inverse_mass / (lx * lx + lz * lz), 0 ,
				             0, 0, 12.0 * this->inverse_mass / (lx * lx + ly * ly) };
	} else {
		this->inverse_I0 = { 1, 0, 0,
					         0, 1, 0 ,
				             0, 0, 1 };
	}

	this->warp(Vec3(0, 0, 0), Quat(0, 0, 0, 1));
}

void Rigid_Body::warp(Vec3 center, Quat orientation) {
    this->center_of_mass = center;
	this->orientation    = orientation.unit();

	this->frame_force          = Vec3(0, 0, 0);
	this->linear_velocity      = Vec3(0, 0, 0);
	
    this->frame_torque          = Vec3(0, 0, 0);
    this->angular_velocity      = Vec3(0, 0, 0);
	this->angular_momentum      = Vec3(0, 0, 0);

    this->build_inertia_tensor();
	this->build_transformation_matrix();
}

void Rigid_Body::build_transformation_matrix() {
	Mat4 translation_matrix, rotation_matrix, scale_matrix;
	translation_matrix.initTranslation(this->center_of_mass.x, this->center_of_mass.y, this->center_of_mass.z);
	rotation_matrix = this->orientation.getRotMat();
	scale_matrix.initScaling(this->size.x, this->size.y, this->size.z);

	this->transformation = (scale_matrix * rotation_matrix) * translation_matrix;
}

void Rigid_Body::build_inertia_tensor() {
    // Copied from https://github.com/bulletphysics/bullet3/blob/master/src/BulletDynamics/Dynamics/btRigidBody.cpp
    Mat3 rotation_matrix = mat3_from_quaternion(this->orientation);
    Mat3 transposed_rotation_matrix = mat3_tranpose(rotation_matrix);
    this->inverse_inertia = mat3_mul_mat3(mat3_mul_mat3(rotation_matrix, this->inverse_I0), transposed_rotation_matrix);
}

void Rigid_Body::apply_force(Vec3 world_space_position, Vec3 force) {
	Vec3 position_relative_to_center = world_space_position - this->center_of_mass;

	this->frame_force  += force;
	this->frame_torque += cross(position_relative_to_center, force);
}

void Rigid_Body::apply_impulse(Vec3 world_space_position, Vec3 impulse) {
	Vec3 position_relative_to_center = world_space_position - this->center_of_mass;

	this->linear_velocity  += impulse * this->inverse_mass;
	this->angular_momentum += cross(position_relative_to_center, impulse);
}

void Rigid_Body::apply_angular_impulse(Vec3 impulse) {
    this->angular_momentum += impulse;
}

void Rigid_Body::apply_torque(Vec3 torque) {
	this->frame_torque += torque;
}

Vec3 Rigid_Body::get_world_space_velocity_at(Vec3 world_space_position) {
	Vec3 position_relative_to_center = world_space_position - this->center_of_mass;

	return this->linear_velocity + cross(this->angular_velocity, position_relative_to_center);
}


//
// Heat Diffusion.
//

void Heat_Grid::create(int width, int height) {
    this->destroy();

    this->width  = width;
    this->height = height;
    this->values = (float *) malloc(this->width * this->height * sizeof(float));
    assert(this->values != NULL);
    this->reset();
}

void Heat_Grid::destroy() {
    if(this->values) free(this->values);
    this->values = NULL;
    this->width  = 0;
    this->height = 0;
}

void Heat_Grid::reset() {
    memset(this->values, 0, this->width * this->height * sizeof(float));
    this->apply_boundary_condition();
}

void Heat_Grid::apply_boundary_condition() {
    const float boundary_value = 0.f;

    for(int x = 0; x < this->width; ++x) {
        this->set(x, 0, boundary_value);
        this->set(x, this->height - 1, boundary_value);
    }

    for(int y = 1; y < this->height - 1; ++y) {
        this->set(0, y, boundary_value);
        this->set(this->width - 1, y, boundary_value);
    }
}

void Heat_Grid::randomize() {
    for(int y = 0; y < this->height; ++y) {
        for(int x = 0; x < this->width; ++x) {
            this->set(x, y, random_float(-1, 1));
        }
    }

    this->apply_boundary_condition();
}

void Heat_Grid::set(int x, int y, float value) {
    assert(x >= 0 && x < this->width);
    assert(y >= 0 && y < this->height);
    this->values[y * this->width + x] = value;
}

float Heat_Grid::get(int x, int y) {
    assert(x >= 0 && x < this->width);
    assert(y >= 0 && y < this->height);
    return this->values[y * this->width + x];
}


//
// Simulator
//

OpenProjectSimulator::OpenProjectSimulator() {
    setup_timing();
}

const char * OpenProjectSimulator::getTestCasesStr() {
    return "Open Project";
}

void OpenProjectSimulator::initUI(DrawingUtilitiesClass * DUC) {
    // @Incomplete: Maybe add a gravity, spring_damping, heat_alpha UI.
    this->DUC = DUC;
}

void OpenProjectSimulator::reset() {
    // @Incomplete.
    this->draw_requests = DRAW_EVERYTHING;
    this->gravity = -10.0f;

    this->masspoint_count = 0;
    this->spring_count    = 0;
    this->spring_damping  = 0.0f; // No damping

    this->rigid_body_count = 0;
    
    this->heat_alpha = 0.5f; // Half decay.
    
    this->setup_game();    
}

void OpenProjectSimulator::drawFrame(ID3D11DeviceContext * context) {
    this->draw_game();
}

void OpenProjectSimulator::notifyCaseChanged(int testCase) {
    // This is called when pressing the "Reset Scene" button in the UI...
    this->reset();
} 

void OpenProjectSimulator::simulateTimestep(float timestep) {
    //
    // :TimeStep
    //
#if USE_FIXED_DT
    double now = get_current_time_in_milliseconds();

    while(now - this->time_of_previous_update > FIXED_DT) {
        this->update_game(FIXED_DT);
        this->time_of_previous_update += FIXED_DT;
    }
#else
    this->update_game(timestep);
#endif
}

void OpenProjectSimulator::externalForcesCalculations(float timeStep) {
    // @Incomplete
}

void OpenProjectSimulator::onClick(int x, int y) {}

void OpenProjectSimulator::onMouse(int x, int y) {}


int OpenProjectSimulator::create_masspoint(Vec3 position, Real mass) {
    assert(this->masspoint_count < MAX_MASSPOINTS);
    int index = this->masspoint_count;

    Masspoint & m  = this->masspoints[index];
    m.position     = position;
    m.velocity     = Vec3(0, 0, 0);
    m.inverse_mass = mass > 0.0f ? 1.0f / mass : 0.0f;

    ++this->masspoint_count;
    return index;
}

int OpenProjectSimulator::create_spring(int a, int b, Real initial_length, Real stiffness) {
    assert(this->spring_count < MAX_SPRINGS);
    int index = this->spring_count;

    Spring & s = this->springs[index];
    s.a = a;
    s.b = b;
    s.stiffness = stiffness;
    s.initial_length = initial_length;
    s.current_length = 0;
    s.current_force_from_a_to_b = Vec3(0, 0, 0);

    ++this->spring_count;
    return index;
}

int OpenProjectSimulator::create_rigid_body(Vec3 size, Real mass) {
    assert(this->rigid_body_count < MAX_RIGID_BODIES);
    int index = this->rigid_body_count;

    Rigid_Body & body = this->rigid_bodies[index];
    body.create(size, mass);
    body.albedo = Vec3(random_float(0, 1), random_float(0, 1), random_float(0, 1));
    
    ++this->rigid_body_count;
    return index;
}

Rigid_Body & OpenProjectSimulator::query_rigid_body(int index) {
    assert(index >= 0 && index < this->rigid_body_count);
    return this->rigid_bodies[index];
}

Rigid_Body & OpenProjectSimulator::create_and_query_rigid_body(Vec3 size, Real mass) {
    return this->query_rigid_body(this->create_rigid_body(size, mass));
}

void OpenProjectSimulator::apply_impulse_to_masspoint(int index, Vec3 impulse) {
    assert(index >= 0 && index < this->masspoint_count);
    this->masspoints[index].velocity += impulse;
}

void OpenProjectSimulator::apply_impulse_to_rigid_body(int index, Vec3 world_space_position, Vec3 impulse) {
    assert(index >= 0 && index < this->rigid_body_count);
    this->rigid_bodies[index].apply_impulse(world_space_position, impulse);
}

void OpenProjectSimulator::apply_torque_to_rigid_body(int index, Vec3 torque) {
    assert(index >= 0 && index < this->rigid_body_count);
    this->rigid_bodies[index].apply_torque(torque);
}

void OpenProjectSimulator::warp_rigid_body(int index, Vec3 position, Quat orientation) {
    assert(index >= 0 && index < this->rigid_body_count);
    this->rigid_bodies[index].warp(position, orientation);
}

void OpenProjectSimulator::setup_game() {
    //
    // Set up the springs.
    //
    if(false) {
        int a = this->create_masspoint(Vec3(0, 0, 0), 10);
        int b = this->create_masspoint(Vec3(0, 2, 0), 10);
        this->apply_impulse_to_masspoint(a, Vec3(-1, 0, 0));
        this->apply_impulse_to_masspoint(b, Vec3( 1, 0, 0));
        this->create_spring(a, b, 1.f, 40.f);
    }
    
    //
    // Set up a rigid body.
    //
    if(true) {        
        Rigid_Body & floor      = this->create_and_query_rigid_body(Vec3(4, 1, 4), 0);
        Rigid_Body & wall_north = this->create_and_query_rigid_body(Vec3(4, 4, .2), 0);
        Rigid_Body & wall_south = this->create_and_query_rigid_body(Vec3(4, 4, .2), 0);
        Rigid_Body & wall_east  = this->create_and_query_rigid_body(Vec3(.2, 4, 4), 0);
        Rigid_Body & wall_west  = this->create_and_query_rigid_body(Vec3(.2, 4, 4), 0);
        
        wall_north.warp(Vec3(0, 2, -2), Quat(0, 0, 0, 1));
        wall_south.warp(Vec3(0, 2,  2), Quat(0, 0, 0, 1));
        wall_east .warp(Vec3(-2, 2, 0), Quat(0, 0, 0, 1));
        wall_west .warp(Vec3( 2, 2, 0), Quat(0, 0, 0, 1));
        
        
        Real min = -2 + 0.1 + 0.5;
        Real max =  2 - 0.1 - 0.5;

        for(int i = 0; i < 10; ++i) {
            Rigid_Body & ball = this->create_and_query_rigid_body(Vec3(.5, .5, .5), 1);
            ball.warp(Vec3(random_float(min, max), 3, random_float(min, max)), quat_from_euler_angles(random_float(0, 1), random_float(0, 1), random_float(0, 1)));
            ball.apply_angular_impulse(Vec3(1, 1, 1));
        }
       
        
        this->gravity = -10;
    }
    
    //
    // Set up the heat grid.
    //
    if(false) {
        this->heat_grid.create(16, 16);
        this->heat_grid.randomize();
    }

    //
    // Set up the timing info.
    //
    this->time_of_previous_update = get_current_time_in_milliseconds();
}

void OpenProjectSimulator::update_game(float dt) {
    // @Incomplete
    //
    // Update the mass-spring-system using the Midpoint method.
    // @Incomplete: Add external forces to each masspoint.
    //
    {
        Vec3 temp_positions[MAX_MASSPOINTS];  // x(t + h/2)
        Vec3 temp_velocities[MAX_MASSPOINTS]; // v(t + h/2)
        float dt_2 = 0.5f * dt;

        //
        // Calculate x(t + h/2) based on v(t) for all mass points. Apply gravity to v(t + h/2).
        //

        for(int i = 0; i < this->masspoint_count; ++i) {
            Masspoint & masspoint = this->masspoints[i];
            if(masspoint.inverse_mass == 0.0f) {
                temp_velocities[i] = Vec3(0, 0, 0);
                temp_positions[i]  = masspoint.position;
                continue;
            }

            temp_positions[i]  = masspoint.position + dt_2 * masspoint.velocity;
            temp_velocities[i] = masspoint.velocity + Vec3(0, -dt_2 * this->gravity, 0);
        }

        //
        // Calculate a(t) for all springs
        // Calculate v(t + h/2) based on a(t) of all springs for all masspoints
        //

        for(int i = 0; i < this->spring_count; ++i) {
            Spring & spring = this->springs[i];
            Masspoint & a = this->masspoints[spring.a], & b = this->masspoints[spring.b];

            // Calculate a(t) for all springs, based on x(t) and v(t) of a, b

            Vec3 direction = a.position - b.position;
            spring.current_length = norm(direction);
            spring.current_force_from_a_to_b = -spring.stiffness * (spring.current_length - spring.initial_length) * (direction / spring.current_length);

            // Calculate v(t + h/2) based on the spring forces for all masspoints

            if(a.inverse_mass != 0.f) temp_velocities[spring.a] += dt_2 * a.inverse_mass * (spring.current_force_from_a_to_b);
            if(b.inverse_mass != 0.f) temp_velocities[spring.b] -= dt_2 * b.inverse_mass * (spring.current_force_from_a_to_b);
        }

        //
        // Apply damping to v(t + h/2)
        //

        for(int i = 0; i < this->masspoint_count; ++i) {
            Masspoint & masspoint = this->masspoints[i];
            temp_velocities[i] -= dt * masspoint.inverse_mass * temp_velocities[i] * this->spring_damping;
        }

        //
        // Calculate x(t + h) using v(t + h/2). Apply gravity to v(t + h).
        //

        for(int i = 0; i < this->masspoint_count; ++i) {
            Masspoint & masspoint = this->masspoints[i];
            masspoint.position += dt * temp_velocities[i];
            masspoint.velocity += Vec3(0, -dt * this->gravity, 0);
        }

        //
        // Calculate a(t + h/2) based on x(t + h/2) and v(t + h/2) for all springs
        // Calculate v(t + h) based on the spring forces for all masspoints
        //

        for(int i = 0; i < this->spring_count; ++i) {
            // Calculate a(t + h/2) based on x(t + h/2) and v(t + h/2) for all springs

            Spring & spring = this->springs[i];
            Vec3 direction = (temp_positions[spring.a] + dt_2 * temp_velocities[spring.a]) - (temp_positions[spring.b] + dt_2 * temp_velocities[spring.b]);
            Real current_length = norm(direction);
            spring.current_force_from_a_to_b = -spring.stiffness * (current_length - spring.initial_length) * (direction / current_length);

            // Calculate v(t + h) based on the spring forces for all masspoints

            Masspoint & a = this->masspoints[spring.a], & b = this->masspoints[spring.b];

            if(a.inverse_mass != 0.f) a.velocity += dt * a.inverse_mass * spring.current_force_from_a_to_b;
            if(b.inverse_mass != 0.f) b.velocity -= dt * b.inverse_mass * spring.current_force_from_a_to_b;
        }

        //
        // Apply damping to v(t + h)
        //

        for(int i = 0; i < this->masspoint_count; ++i) {
            Masspoint & masspoint = this->masspoints[i];
            masspoint.velocity -= dt * masspoint.inverse_mass * masspoint.velocity * this->spring_damping;
        }
    }

    //
    // Update the heat grid using the implicit method.
    //
    if(this->heat_grid.width > 0 && this->heat_grid.height > 0) {
        //
        // The "world space" distance between each point, for now just one unit.
        //
        float dx = 1.f;
        float dy = 1.f;

        //
        // The parameters for the system of equations.
        //
        int N = this->heat_grid.height * this->heat_grid.width;
        SparseMatrix<float> A(N);
        std::vector<float> b(N, 0.0f);

        //
        // Prepare the matrix A and vector b, according to:
        // https://hplgit.github.io/fdm-book/doc/pub/diffu/html/._diffu-solarized001.html
        // a, b, c are constant for each entry.
        //

        float i_dt = 1.f / dt;
        float theta = 0.5f;
        float Fx = this->heat_alpha * dt / (dx * dx);
        float Fy = this->heat_alpha * dt / (dy * dy);


#define index(x, y) ((y) * this->heat_grid.width + (x))

        // Boundary condition on the top row.
        for(int x = 0; x < this->heat_grid.width; ++x) {
            int i = index(x, 0);
            A.set_element(i, i, 1.0f);
            b[i] = 0;
        }

        // Fill in all the inner rows.
        for(int y = 1; y < this->heat_grid.height - 1; ++y) {
            int x = 0;
            int i = index(x, y);
            A.set_element(i, i, 1.0f); // Boundary condition on the left column
            b[i] = 0.f;

            for(x = 1; x < this->heat_grid.width - 1; ++x) {
                int i = index(x, y);

                A.set_element(i, index(x, y - 1), -theta * Fx);
                A.set_element(i, index(x - 1, y), -theta * Fy);

                A.set_element(i, i, 1 + 2 * theta * (Fx + Fy));

                A.set_element(i, index(x + 1, y), -theta * Fx);
                A.set_element(i, index(x, y + 1), -theta * Fy);

                b[i] = this->heat_grid.get(x, y) +
                    (1 - this->heat_alpha) * (
                        Fx * (this->heat_grid.get(x + 1, y) - 2 * this->heat_grid.get(x, y) + this->heat_grid.get(x - 1, y)) +
                        Fy * (this->heat_grid.get(x, y + 1) - 2 * this->heat_grid.get(x, y) + this->heat_grid.get(x, y - 1))
                    );
            }

            x = this->heat_grid.width - 1;
            i = index(x, y);
            A.set_element(i, i, 1.0f); // Boundary condition on the right column
            b[i] = 0.f;
        }

        // Boundary condition on the bottom row.
        for(int x = 0; x < this->heat_grid.width; ++x) {
            int i = index(x, this->heat_grid.height - 1);
            A.set_element(i, i, 1.0f);
            b[i] = 0;
        }

#undef index

        //
        // Set up the solution vector.
        //
        std::vector<float> T(this->heat_grid.width * this->heat_grid.height, 0.f);

        //
        // Set up the matrix solver.
        //
        float pcg_target_residual = 1e-05f;
        float pcg_max_iterations  = 1000.f;
        float ret_pcg_residual    = 1e10f;
        int ret_pcg_iterations  = -1;

        SparsePCGSolver<float> solver;
        solver.set_solver_parameters(pcg_target_residual, pcg_max_iterations, 0.97f, 0.25f);

        //
        // Actually solve the system of equations.
        //
        solver.solve(A, b, T, ret_pcg_residual, ret_pcg_iterations, 1);

        //
        // Extract the new temperatur values from the solution vector.
        //
        for(int y = 0; y < this->heat_grid.height; ++y) {
            for(int x = 0; x < this->heat_grid.width; ++x) {
                int index = y * this->heat_grid.width + x;
                this->heat_grid.set(x, y, T[index]);
            }
        }

        this->heat_grid.apply_boundary_condition();
    }

    //
    // Update the rigid body system.
    //
    {
        //
        // Add velocity and external forces to all rigid bodies.
        //
        for(int i = 0; i < this->rigid_body_count; ++i) {
            Rigid_Body & body = this->rigid_bodies[i];
    	    if(body.inverse_mass == 0) continue;
            
            // @Incomplete: Add external forces.
            body.linear_velocity.y += dt * this->gravity; // Add gravity
        }

        //
        // Resolve collisions between the rigid bodies.
        //
        for(int i = 0; i < this->rigid_body_count; ++i) {
            Rigid_Body & lhs = this->rigid_bodies[i];
            SAT_Input sat_lhs = sat_input(lhs);

            for(int j = i + 1; j < this->rigid_body_count; ++j) {
                Rigid_Body & rhs = this->rigid_bodies[j];
                if(lhs.inverse_mass == 0 && rhs.inverse_mass == 0) continue; // Don't do collisions between two static objects.

                //
                // Check for a collision between the two bodies.
                //
                SAT_Input sat_rhs = sat_input(rhs);
                SAT_Result result = sat(sat_lhs, sat_rhs);
                if(!result.found_collision) continue;
                
                Real point_factor = 1.0 / (Real) result.world_space_position_count;
                Vec3 contact_normal = Vec3(result.normal.x, result.normal.y, result.normal.z);
                
                //
                // Handle each contact point.
                //
                for(int k = 0; k < result.world_space_position_count; ++k) {
                    Vec3 contact_point = Vec3(result.world_space_positions[k].x, result.world_space_positions[k].y, result.world_space_positions[k].z);
                    Vec3 contact_point_velocity_lhs = lhs.get_world_space_velocity_at(contact_point);
                    Vec3 contact_point_velocity_rhs = rhs.get_world_space_velocity_at(contact_point);

                    Vec3 relative_velocity = contact_point_velocity_lhs - contact_point_velocity_rhs;
                    Real rv_dot_normal = dot(relative_velocity, contact_normal);
                    if(rv_dot_normal > 0.0) continue; // Bodies are already separating, nothing to do.

                    //
                    // Move the two objects apart along the penetration normal by a tiny bit, to combat
                    // numerical errors over time when stacked.
                    //

                    Real total_mass = lhs.inverse_mass + rhs.inverse_mass;
                    Real correction_factor = 0.1;
                    Vec3 correction_vector = contact_normal * correction_factor;
                    lhs.center_of_mass += correction_vector * lhs.inverse_mass / total_mass;
                    rhs.center_of_mass -= correction_vector * rhs.inverse_mass / total_mass;

                    //
                    // Do a proper collision response, by calculating the impulse between the two bodies,
                    // and applying the respective part to each of them.
                    //

                    const Real restitution = 1;

                    Vec3 collision_point_on_lhs = contact_point - lhs.center_of_mass; // xa
                    Vec3 collision_point_on_rhs = contact_point - rhs.center_of_mass; // xb

                    Vec3 collision_point_on_lhs_cross_normal = cross(collision_point_on_lhs, contact_normal);
                    Vec3 collision_point_on_rhs_cross_normal = cross(collision_point_on_rhs, contact_normal);

                    Vec3 applied_inerta_on_lhs = cross(mat3_mul_vec3(lhs.inverse_inertia, collision_point_on_lhs_cross_normal), collision_point_on_lhs);
                    Vec3 applied_inerta_on_rhs = cross(mat3_mul_vec3(rhs.inverse_inertia, collision_point_on_rhs_cross_normal), collision_point_on_rhs);

                    Real impulse_magnitude_nominator   = -(1 + restitution) * rv_dot_normal;
                    Real impulse_magnitude_denominator = lhs.inverse_mass + rhs.inverse_mass + dot(applied_inerta_on_lhs + applied_inerta_on_rhs, contact_normal);

                    Real impulse_magnitude = impulse_magnitude_nominator / impulse_magnitude_denominator; // AKA Big 'J'.
                    
                    Vec3 impulse = contact_normal * impulse_magnitude * point_factor;

                    //
                    // Apply the impulse to both bodies.
                    //

                    lhs.apply_impulse(contact_point,  impulse);
                    rhs.apply_impulse(contact_point, -impulse);
                }
            }
        }

        //
        // Integrate all rigid bodies.
        //
        for(int i = 0; i < this->rigid_body_count; ++i) {
            Rigid_Body & body = this->rigid_bodies[i];
    	    if(body.inverse_mass == 0) continue;

            float dt_2 = dt * 0.5;

            //
            // Integrate the linear part using Euler.
            //

            body.center_of_mass   = body.center_of_mass  + body.linear_velocity * dt; // Integrate the position
            body.linear_velocity  = body.linear_velocity + body.frame_force * body.inverse_mass * dt; // Integrate the velocity.

            // Reset the accumulators every frame.
            body.frame_force = { 0 };
            
            //
            // Integrate the angular part.
            //

            // Integrate the angular momentum
            body.angular_momentum = body.angular_momentum + dt * body.frame_torque;
            
            // Calculate the new inertia tensor
            body.build_inertia_tensor();
                        
            // Calculate the angular velocity for this frame.
            body.angular_velocity = mat3_mul_vec3(body.inverse_inertia, body.angular_momentum);

            // Integrate the rotation
            body.orientation = body.orientation + dt_2 * Quat(body.angular_velocity.x, body.angular_velocity.y, body.angular_velocity.z, 0) * body.orientation;
            body.orientation = body.orientation.unit();

            // Reset the accumulators every frame.
            body.frame_torque = { 0 };
            
            //
            // Finally build the new transformation matrices.
            //
            
            body.build_transformation_matrix();
        }
    }

    //this->debug_print();
}

void OpenProjectSimulator::draw_game() {
    // @Incomplete.

    //
    // Draw all masspoints and springs if requested.
    //
    if(this->draw_requests & DRAW_SPRINGS) {
        const float sphere_radius = 0.1f;

        for(int i = 0; i < this->masspoint_count; ++i) {
            Masspoint & masspoint = this->masspoints[i];

            float r = (i % this->masspoint_count)  / (float) this->masspoint_count;
            float g = ((i + this->masspoint_count / 2) % this->masspoint_count) / (float) this->masspoint_count;
            float b = ((this->masspoint_count + 1 - i) % this->masspoint_count) / (float) this->masspoint_count;

            this->DUC->setUpLighting(Vec3(0, 0, 0), Vec3(r, g, b), 1, Vec3(r, g, b));
            this->DUC->drawSphere(masspoint.position, Vec3(sphere_radius, sphere_radius, sphere_radius));
        }

        this->DUC->beginLine();

        for(int i = 0; i < this->spring_count; ++i) {
            Spring & spring = this->springs[i];

            float tension = spring.current_length / spring.initial_length;

            Vec3 color = lerp(Vec3(0, 1, 0), Vec3(1, 0, 00), clamp_float(tension / 3, 0, 1));

            this->DUC->drawLine(this->masspoints[spring.a].position, color, this->masspoints[spring.b].position, color);
        }

        this->DUC->endLine();
    }

    //
    // Draw the heat grid if requested.
    //
    if(this->draw_requests & DRAW_HEAT_MAP) {
        float dominating_grid_size = (this->heat_grid.width > this->heat_grid.height) ? this->heat_grid.width : this->heat_grid.height;

        float visual_grid_width = 1, visual_grid_height = 1; // In world space.
        float visual_sphere_radius = 2 * (visual_grid_width / dominating_grid_size); // Calculated so that the spheres are always twice the "perfect" radius. "Perfect" meaning the spheres would just barely touch each other.
        float visual_grid_offset_x = -visual_grid_width / 2 + visual_sphere_radius / 2, visual_grid_offset_y = -visual_grid_height / 2 + visual_sphere_radius / 2;

        Vec3 medium_color = Vec3(0, 0, 0);
        Vec3 cold_color   = Vec3(.8, .2, .2);
        Vec3 hot_color    = Vec3(1, 1, 1);

        for(int y = 0; y < this->heat_grid.height; ++y) {
            for(int x = 0; x < this->heat_grid.width; ++x) {
                float visual_x = (x / (float) this->heat_grid.width)  * visual_grid_width  + visual_grid_offset_x;
                float visual_y = (y / (float) this->heat_grid.height) * visual_grid_height + visual_grid_offset_y;
                float value = this->heat_grid.get(x, y);

                Vec3 position = Vec3(visual_x, visual_y, 0);
                Vec3 scale    = Vec3(visual_sphere_radius, visual_sphere_radius, visual_sphere_radius);

                Vec3 color;
                if(value >= 0.0)
                    color = lerp(medium_color, hot_color, value);
                else
                    color = lerp(medium_color, cold_color, -value);

                this->DUC->setUpLighting(Vec3(0, 0, 0), color, 5, color);
                this->DUC->drawSphere(position, scale);
            }
        }
    }

    //
    // Draw all rigid bodies if requested.
    //
    if(this->draw_requests & DRAW_RIGID_BODIES) {
        for(int i = 0; i < this->rigid_body_count; ++i) {
            Rigid_Body & body = this->rigid_bodies[i];
            this->DUC->setUpLighting(Vec3(0, 0, 0), body.albedo, 0.2, body.albedo);
            this->DUC->drawRigidBody(body.transformation);
        }
    }
}

void OpenProjectSimulator::debug_print() {
    printf("================================================================================\n");

    if(this->draw_requests & DRAW_SPRINGS) {
        printf(" > Masspoints:\n");
        for(int i = 0; i < this->masspoint_count; ++i) {
            Masspoint & masspoint = this->masspoints[i];
            printf("    %u > Position = " PRINT_FIXED_VEC3 ", Velocity = " PRINT_FIXED_VEC3 ", Inverse Mass: " PRINT_FIXED_FLOAT "\n",
                   i, masspoint.position.x, masspoint.position.y, masspoint.position.z, masspoint.velocity.x, masspoint.velocity.y, masspoint.velocity.z, masspoint.inverse_mass);
        }

        printf(" > Springs:\n");
        for(int i = 0; i < this->spring_count; ++i) {
            Spring & spring = this->springs[i];
            printf("    %u > a = %u, b = %u, Initial Length = " PRINT_FIXED_FLOAT ", Current Length = " PRINT_FIXED_FLOAT ", Force: " PRINT_FIXED_VEC3 "\n",
                   i, spring.a, spring.b, spring.initial_length, spring.current_length, spring.current_force_from_a_to_b.x, spring.current_force_from_a_to_b.y, spring.current_force_from_a_to_b.z);
        }
    }

    if(this->draw_requests & DRAW_HEAT_MAP) {
        printf(" > Heat Grid: width = %u, height = %u.\n", this->heat_grid.width, this->heat_grid.height);

        for(int y = 0; y < this->heat_grid.height; ++y) {
            printf("    ");
            for(int x = 0; x < this->heat_grid.width; ++x) {
                printf(PRINT_FIXED_FLOAT, this->heat_grid.get(x, y));
                printf(" ");
            }

            printf("\n");
        }

        if(this->heat_grid.height) printf("\n");
    }

    if(this->draw_requests & DRAW_RIGID_BODIES) {
        printf(" > Rigid Bodies:\n");
        for(int i = 0; i < this->rigid_body_count; ++i) {
            Rigid_Body & body = this->rigid_bodies[i];
            printf("    %u > Position = " PRINT_FIXED_VEC3 ", Linear Velocity = " PRINT_FIXED_VEC3 ", Angular Momentum = " PRINT_FIXED_VEC3 "\n",
                i, body.center_of_mass.x, body.center_of_mass.y, body.center_of_mass.z, body.linear_velocity.x, body.linear_velocity.y, body.linear_velocity.z,
                body.angular_momentum.x, body.angular_momentum.y, body.angular_momentum.z);
            printf("         Orientation = { %f, %f, %f, %f }, Angular Velocity = " PRINT_FIXED_VEC3 "\n", 
                body.orientation.x, body.orientation.y, body.orientation.z, body.orientation.w, body.angular_velocity.x, body.angular_velocity.y, body.angular_velocity.z);
        }
    }

    printf("================================================================================\n");
}
