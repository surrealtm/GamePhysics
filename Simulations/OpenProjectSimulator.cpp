#include "OpenProjectSimulator.h"

//
// Utility Functions
//

Vec3 lerp(Vec3 from, Vec3 to, float t) {
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

//
// Simulator
//

OpenProjectSimulator::OpenProjectSimulator() {
	this->reset();
}

const char * OpenProjectSimulator::getTestCasesStr() {
	return "Open Project";
}

void OpenProjectSimulator::initUI(DrawingUtilitiesClass * DUC) {
	// @Incomplete: Maybe add a timestep, gravity UI.
	this->duc = DUC;
}

void OpenProjectSimulator::reset() {
	// @Incomplete.
	this->draw_requests = DRAW_EVERYTHING;
	this->gravity = 0.0f;

	this->masspoint_count = 0;
	this->spring_count    = 0;
	this->spring_damping  = 0.0f; // No damping

	this->setup_game();
}

void OpenProjectSimulator::drawFrame(ID3D11DeviceContext * context) {
	// @Incomplete.

	//
	// Draw all masspoints and springs if requested.
	//
	if(this->draw_requests & DRAW_SPRINGS) {
		for(int i = 0; i < this->masspoint_count; ++i) {
			Masspoint & masspoint = this->masspoints[i];

			float r = (i % this->masspoint_count)  / (float) this->masspoint_count;
			float g = ((i + this->masspoint_count / 2) % this->masspoint_count) / (float) this->masspoint_count;
			float b = ((this->masspoint_count + 1 - i) % this->masspoint_count) / (float) this->masspoint_count;
        
			this->duc->setUpLighting(Vec3(0, 0, 0), Vec3(r, g, b), 1, Vec3(r, g, b));
			this->duc->drawSphere(masspoint.position, Vec3(DRAW_SPHERE_RADII, DRAW_SPHERE_RADII, DRAW_SPHERE_RADII));
		}

		this->duc->beginLine();

		for(int i = 0; i < this->spring_count; ++i) {
			Spring & spring = this->springs[i];

			float tension = spring.current_length / spring.initial_length;

			Vec3 color = lerp(Vec3(0, 1, 0), Vec3(1, 0, 00), clamp_float(tension / 3, 0, 1));

			this->duc->drawLine(this->masspoints[spring.a].position, color, this->masspoints[spring.b].position, color);
		}

		this->duc->endLine();
	}
}

void OpenProjectSimulator::notifyCaseChanged(int testCase) {} // We don't have different test cases, so just ignore this.

void OpenProjectSimulator::simulateTimestep(float dt) {
	// @Incomplete

	//
	// Update the mass-spring-system using the Midpoint method.
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
            float current_length = norm(direction);
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

	//this->debug_print();
}

void OpenProjectSimulator::externalForcesCalculations(float timeStep) {
	// @Incomplete
}

void OpenProjectSimulator::onClick(int x, int y) {}

void OpenProjectSimulator::onMouse(int x, int y) {}


int OpenProjectSimulator::create_masspoint(Vec3 position, float mass) {
	assert(this->masspoint_count < MAX_MASSPOINTS);
	int index = this->masspoint_count;
	
	Masspoint & m  = this->masspoints[index];
	m.position     = position;
	m.velocity     = Vec3(0, 0, 0);
	m.inverse_mass = mass > 0.0f ? 1.0f / mass : 0.0f;

	++this->masspoint_count;
	return index;
}

int OpenProjectSimulator::create_spring(int a, int b, float initial_length, float stiffness) {
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

void OpenProjectSimulator::setup_game() {
	int a = this->create_masspoint(Vec3(0, 0, 0), 2);
	int b = this->create_masspoint(Vec3(1, 1, 0), 2);
	this->create_spring(a, b, 0.5f, 40.f);
}

void OpenProjectSimulator::debug_print() {
	printf("================================================================================\n");

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

	printf("================================================================================\n");
}