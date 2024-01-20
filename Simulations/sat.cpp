#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <limits> // For the max double value...

#include "sat.h"

/* =============================== SAT_Vec3 =============================== */

SAT_Scalar sat_abs(SAT_Scalar value) {
    return (value < 0) ? -value : value;
}

SAT_Scalar sat_turns_to_radians(SAT_Scalar value) {
    return value * static_cast<SAT_Scalar>(2 * 3.1415926535);
}

SAT_Vec3 sat_vec3(SAT_Scalar x, SAT_Scalar y, SAT_Scalar z) {
    return { x, y, z };
}

SAT_Vec3 sat_scaled_vec3(SAT_Vec3 vec, SAT_Scalar scale) {
    return { vec.x * scale, vec.y * scale, vec.z * scale };
}

SAT_Quat sat_quat_from_euler_angles(SAT_Scalar x, SAT_Scalar y, SAT_Scalar z) {
    SAT_Scalar rx = sat_turns_to_radians(x / 2);
    SAT_Scalar ry = sat_turns_to_radians(y / 2);
    SAT_Scalar rz = sat_turns_to_radians(z / 2);

    SAT_Scalar cx = cos(rx);
    SAT_Scalar cy = cos(ry);
    SAT_Scalar cz = cos(rz);

    SAT_Scalar sx = sin(rx);
    SAT_Scalar sy = sin(ry);
    SAT_Scalar sz = sin(rz);

    SAT_Quat result;
    result.x = sx * cy * cz - cx * sy * sz;
    result.y = cx * sy * cz + sx * cy * sz;
    result.z = cx * cy * sz - sx * sy * cz;
    result.w = cx * cy * cz - sx * sy * sz;
    return result;
}

SAT_Scalar sat_dot(SAT_Vec3 lhs, SAT_Vec3 rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

SAT_Vec3 sat_cross(SAT_Vec3 lhs, SAT_Vec3 rhs) {
    SAT_Vec3 result;
    result.x = lhs.y * rhs.z - lhs.z - rhs.y;
    result.y = lhs.z * rhs.x - lhs.x * rhs.z;
    result.z = lhs.x * rhs.y - lhs.y * rhs.x;
    return result;
}

SAT_Vec3 sat_rotate(SAT_Quat quat, SAT_Vec3 vec) {
    // The formal math to rotate a vector by a quaternion would be
    // v' = q * v * q'
    // This however is a pretty slow approach, as multiplying two quaternions is slow.
    // This approad instead uses two dot products and one cross product:
    // u is the axis of the quaternion
    // v is the incoming vector
    // result = 2 * dot(u, v) * u
    //          + (q.w * q.w - dot(u, u)) * v
    //          + 2 * q.w * cross(u, v);

    SAT_Scalar w_squared = quat.w * quat.w;
    SAT_Scalar u_dot_u = quat.x * quat.x + quat.y * quat.y + quat.z * quat.z;
    SAT_Scalar u_dot_v = quat.x * vec.x  + quat.y * vec.y  + quat.z * vec.z;

    SAT_Scalar u_cross_v_x = (quat.y * vec.z - quat.z * vec.y);
    SAT_Scalar u_cross_v_y = (quat.z * vec.x - quat.x * vec.z);
    SAT_Scalar u_cross_v_z = (quat.x * vec.y - quat.y * vec.x);
    
    SAT_Scalar rotated_x = 2 * u_dot_v * quat.x + (w_squared - u_dot_u) * vec.x + 2 * quat.w * u_cross_v_x;
    SAT_Scalar rotated_y = 2 * u_dot_v * quat.y + (w_squared - u_dot_u) * vec.y + 2 * quat.w * u_cross_v_y;
    SAT_Scalar rotated_z = 2 * u_dot_v * quat.z + (w_squared - u_dot_u) * vec.z + 2 * quat.w * u_cross_v_z;

    return { rotated_x, rotated_y, rotated_z };
}


/* =============================== Internal Logic =============================== */

bool __sat_is_separating_axis(SAT_State *state, SAT_Vec3 axis) {
    if(sat_dot(axis, axis) <= 0.00001f) return false; // The axis was spawned from a cross product, but the two vectors were parallel, so no valid axis was created.

    SAT_Scalar lhs = sat_abs(sat_dot(state->a_to_b, axis));
    SAT_Scalar rhs = 0;

    for(int i = 0; i < SAT_AXIS_COUNT; ++i) {
        rhs += sat_abs(sat_dot(state->scaled_axis_a[i], axis));
        rhs += sat_abs(sat_dot(state->scaled_axis_b[i], axis));
    }

    SAT_Scalar overlap = rhs - lhs; // A negative value means there is a gap, no overlap.

    if(overlap >= 0 && overlap < state->penetration_depth) {
        // We have found our new collision axis.
        state->penetration_depth = overlap;
        state->collision_normal  = axis;
    }
    
    return overlap < 0;
}


/* =============================== API implementation =============================== */

SAT_Result sat(SAT_Input lhs, SAT_Input rhs) {
    //
    // Prepare the result structure.
    //
    SAT_Result result = { 0 };

    //
    // Prepare the internal SAT state to save some precomputed values.
    //
    SAT_State state;
    state.center_a  = lhs.center;
    state.unit_axis_a[0]   = sat_rotate(lhs.orientation, { 1, 0, 0 });
    state.unit_axis_a[1]   = sat_rotate(lhs.orientation, { 0, 1, 0 });
    state.unit_axis_a[2]   = sat_rotate(lhs.orientation, { 0, 0, 1 });
    state.scaled_axis_a[0] = sat_scaled_vec3(state.unit_axis_a[0], lhs.size[0]);
    state.scaled_axis_a[1] = sat_scaled_vec3(state.unit_axis_a[1], lhs.size[1]);
    state.scaled_axis_a[2] = sat_scaled_vec3(state.unit_axis_a[2], lhs.size[2]);
    
    state.center_b  = rhs.center;
    state.unit_axis_b[0]   = sat_rotate(rhs.orientation, { 1, 0, 0 });
    state.unit_axis_b[1]   = sat_rotate(rhs.orientation, { 0, 1, 0 });
    state.unit_axis_b[2]   = sat_rotate(rhs.orientation, { 0, 0, 1 });
    state.scaled_axis_b[0] = sat_scaled_vec3(state.unit_axis_b[0], rhs.size[0]);
    state.scaled_axis_b[1] = sat_scaled_vec3(state.unit_axis_b[1], rhs.size[1]);
    state.scaled_axis_b[2] = sat_scaled_vec3(state.unit_axis_b[2], rhs.size[2]);
    
    state.a_to_b = sat_vec3(state.center_b.x - state.center_a.x, state.center_b.y - state.center_a.y, state.center_b.z - state.center_a.z);

    state.found_separating_axis = false;
    state.penetration_depth     = std::numeric_limits<double>::max(); // This is so fucking stupid...
    
    //
    // Check each axis to find a possible separating axis.
    //
    for(int i = 0; i < SAT_AXIS_COUNT; ++i) {
        // Check each of the boxes axis.
        state.found_separating_axis |= __sat_is_separating_axis(&state, state.unit_axis_a[i]);
        state.found_separating_axis |= __sat_is_separating_axis(&state, state.unit_axis_b[i]);

        // Check some edge-to-edge axis.
        state.found_separating_axis |= __sat_is_separating_axis(&state, sat_cross(state.unit_axis_a[i], state.unit_axis_b[i]));
    }
        
    // Check the remaining edge-to-edge axis.
    state.found_separating_axis |= __sat_is_separating_axis(&state, sat_cross(state.unit_axis_a[SAT_AXIS_X], state.unit_axis_b[SAT_AXIS_Y]));
    state.found_separating_axis |= __sat_is_separating_axis(&state, sat_cross(state.unit_axis_a[SAT_AXIS_X], state.unit_axis_b[SAT_AXIS_Z]));
    state.found_separating_axis |= __sat_is_separating_axis(&state, sat_cross(state.unit_axis_a[SAT_AXIS_Y], state.unit_axis_b[SAT_AXIS_Z]));

    state.found_separating_axis |= __sat_is_separating_axis(&state, sat_cross(state.unit_axis_a[SAT_AXIS_Y], state.unit_axis_b[SAT_AXIS_X]));
    state.found_separating_axis |= __sat_is_separating_axis(&state, sat_cross(state.unit_axis_a[SAT_AXIS_Z], state.unit_axis_b[SAT_AXIS_X]));
    state.found_separating_axis |= __sat_is_separating_axis(&state, sat_cross(state.unit_axis_a[SAT_AXIS_Z], state.unit_axis_b[SAT_AXIS_Y]));
    
    //
    // If we have not found a separating axis, it means that there is a collision.
    // Fill in the collision data.
    //
    if(!state.found_separating_axis) {
        result.found_collision = true;
        result.depth  = state.penetration_depth;
        result.normal = state.collision_normal;
    }
    
    //
    // Done.
    //
    
    return result;
}


/* =============================== Test Cases =============================== */

void __sat_tests() {
#define sat_test(name, exp) if(exp) { ++passed; printf("[SAT]: Passed test '%s'.\n", name); } else { ++failed; printf("[SAT]: Failed test '%s'.\n", name); }
    int passed = 0, failed = 0;
    
    {
        SAT_Input lhs = { { 0, 0, 0 }, { 0, 0, 0, 1 }, { 1, 1, 1 } };
        SAT_Input rhs = { { 0, 2, 0 }, { 0, 0, 0, 1 }, { 1, 1, 1 } };
        SAT_Result result = sat(lhs, rhs);
        sat_test("Simple True", result.found_collision);
    }

    {
        SAT_Input lhs = { { 0, 0, 0 },   { 0, 0, 0, 1 }, { 1, 1, 1 } };
        SAT_Input rhs = { { 0, 2.1, 0 }, { 0, 0, 0, 1 }, { 1, 1, 1 } };
        SAT_Result result = sat(lhs, rhs);
        sat_test("Simple False", !result.found_collision);
    }

    {
        SAT_Input lhs = { { 0, 0, 0 },   { 0, 0, 0, 1 }, { 1, 1, 1 } };
        SAT_Input rhs = { { 0, 2.1, 0 }, sat_quat_from_euler_angles(0, 0, 0.25), { 2, 1, 1 } };
        SAT_Result result = sat(lhs, rhs);
        sat_test("Rotated True", result.found_collision);
    }

    {
        SAT_Input lhs = { { 0, 0, 0 },   { 0, 0, 0, 1 }, { 1, 1, 1 } };
        SAT_Input rhs = { { 0, 2.1, 0 }, sat_quat_from_euler_angles(0, 0.25, 0), { 2, 1, 1 } };
        SAT_Result result = sat(lhs, rhs);
        sat_test("Rotated False", !result.found_collision);
    }

    printf("Passed tests: %u, Failed tests: %u\n", passed, failed);

#undef sat_assert
}
