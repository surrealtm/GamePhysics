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

SAT_Vec3 sat_negate(SAT_Vec3 vec) {
    return { -vec.x, -vec.y, -vec.z };
} 

SAT_Scalar sat_dot(SAT_Vec3 lhs, SAT_Vec3 rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

SAT_Vec3 sat_normalized_cross(SAT_Vec3 lhs, SAT_Vec3 rhs) {
    SAT_Vec3 result;
    result.x = lhs.y * rhs.z - lhs.z * rhs.y;
    result.y = lhs.z * rhs.x - lhs.x * rhs.z;
    result.z = lhs.x * rhs.y - lhs.y * rhs.x;


    SAT_Scalar magnitude = sqrt(result.x * result.x + result.y * result.y + result.z * result.z);

    if(magnitude >= 0.0001) {
        result.x /= magnitude;
        result.y /= magnitude;
        result.z /= magnitude;
    }
    
    return result;
}

bool sat_ray_plane_intersection(SAT_Vec3 plane_center, SAT_Vec3 plane_normal, SAT_Vec3 ray_origin, SAT_Vec3 ray_direction, SAT_Vec3 *result) {
    SAT_Scalar denom = sat_dot(plane_normal, ray_direction);
    if(sat_abs(denom) <= 0.0001) return false;

    SAT_Vec3 delta = sat_vec3(plane_center.x - ray_origin.x, plane_center.y - ray_origin.y, plane_center.z - ray_origin.z);
    SAT_Scalar t = sat_dot(delta, plane_normal) / denom;

    if(sat_abs(t) <= 0.0001) return false;

    result->x = ray_origin.x + t * ray_direction.x;
    result->y = ray_origin.y + t * ray_direction.y;
    result->z = ray_origin.z + t * ray_direction.z;
    return true;
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

SAT_Quat sat_quat(SAT_Scalar x, SAT_Scalar y, SAT_Scalar z, SAT_Scalar w) {
    return { x, y, z, w };
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


/* =============================== Internal Logic =============================== */

// The collision normal is always supposed to be the face normal of LHS. If we are checking a normal
// of RHS, then we need to flip that normal around for the result.
bool sat_is_separating_axis(SAT_State *state, SAT_Vec3 axis) {
    if(sat_dot(axis, axis) <= 0.00001) return false; // The axis was spawned from a cross product, but the two vectors were parallel, so no valid axis was created.

    SAT_Scalar lhs = sat_abs(sat_dot(state->lhs_to_rhs, axis));
    SAT_Scalar rhs = 0;

    for(int i = 0; i < SAT_AXIS_COUNT; ++i) {
        rhs += sat_abs(sat_dot(state->scaled_axis[SAT_A][i], axis));
        rhs += sat_abs(sat_dot(state->scaled_axis[SAT_B][i], axis));
    }

    SAT_Scalar overlap = rhs - lhs; // A negative value means there is a gap, no overlap.
    SAT_Scalar axis_dot_relative_velocity = sat_dot(axis, state->relative_velocity);

    // There needs to be actual overlap between the boxes on this axis, and the
    // axis must also not be in the same direction as the relative velocity. This
    // removes collisions with "backfaces", improving the tunneling situation a bit.
    bool is_collision_axis = overlap > 0;

    bool is_better_collision_axis = is_collision_axis;

#if SAT_USE_MINIMAL_TRANSLATION_VECTOR
    is_better_collision_axis &= overlap < state->penetration_depth;
#elif SAT_USE_PROJECTED_VELOCITY
    is_better_collision_axis &= axis_dot_relative_velocity > state->collision_normal_dot_relative_velocity;
#endif

    if(is_better_collision_axis) {
        // We have found our new collision axis.
        state->collision_normal_dot_relative_velocity = axis_dot_relative_velocity;
        state->penetration_depth = overlap;
        state->collision_normal  = axis;

        // Ensure the collision normal always goes from RHS to LHS
        if(sat_dot(state->lhs_to_rhs, axis) > 0.0) state->collision_normal = sat_negate(state->collision_normal);
    }
    
    return !is_collision_axis;
}

void sat_find_significant_face(SAT_State *state, SAT_Significant_Face face_index, SAT_Vec3 face_normal, SAT_Vec3 scaled_normal, SAT_Vec3 scaled_u, SAT_Vec3 scaled_v, SAT_Vec3 center) {
    SAT_Vec3 signed_normal = (face_index == SAT_SIGNIFICANT_FACE_incident) ? sat_negate(state->collision_normal) : state->collision_normal; // The reference face should point the exact opposite way.

    SAT_Scalar dot = sat_dot(face_normal, signed_normal);
    if(dot > state->significant_face[face_index].face_normal_dot_collision_normal) {
        //
        // We have found a face more parallel to the normal than the previous one.
        //

        state->significant_face[face_index].face_normal = face_normal;
        state->significant_face[face_index].face_center = { center.x + scaled_normal.x, center.y + scaled_normal.y, center.z + scaled_normal.z };
        state->significant_face[face_index].face_normal_dot_collision_normal = dot;

        //
        // Calculate the corner positions of this face, based on the given normal, u and v directions.
        //
        
        state->significant_face[face_index].corners[0] = {
            center.x + scaled_normal.x + scaled_u.x + scaled_v.x,
            center.y + scaled_normal.y + scaled_u.y + scaled_v.y,
            center.z + scaled_normal.z + scaled_u.z + scaled_v.z,
        };

        state->significant_face[face_index].corners[1] = {
            center.x + scaled_normal.x + scaled_u.x - scaled_v.x,
            center.y + scaled_normal.y + scaled_u.y - scaled_v.y,
            center.z + scaled_normal.z + scaled_u.z - scaled_v.z,
        };

        state->significant_face[face_index].corners[2] = {
            center.x + scaled_normal.x - scaled_u.x - scaled_v.x,
            center.y + scaled_normal.y - scaled_u.y - scaled_v.y,
            center.z + scaled_normal.z - scaled_u.z - scaled_v.z,
        };

        state->significant_face[face_index].corners[3] = {
            center.x + scaled_normal.x - scaled_u.x + scaled_v.x,
            center.y + scaled_normal.y - scaled_u.y + scaled_v.y,
            center.z + scaled_normal.z - scaled_u.z + scaled_v.z,
        };
    }
}

SAT_Vec3 sat_clip_point_against_plane(SAT_Vec3 point, SAT_Vec3 edge_direction, SAT_Vec3 plane_center, SAT_Vec3 plane_normal) {
    //
    // Check if the point is "inside" the plane like this:
    //   P *      | -> Face
    //
    SAT_Vec3 center_to_point = { point.x - plane_center.x, point.y - plane_center.y, point.z - plane_center.z };
    SAT_Scalar dot = sat_dot(center_to_point, plane_normal);
    if(dot <= 0.0) return point; // No clipping required.

    //
    // The point is "outside" the plane, clip it onto the plane.
    //
    SAT_Vec3 projected_point;
    bool valid_intersection = sat_ray_plane_intersection(plane_center, plane_normal, point, edge_direction, &projected_point);
    if(!valid_intersection) return point; // No clipping possible.
    
    return projected_point;
}

void sat_clip_incident_corners_against_plane(SAT_State *state, SAT_Vec3 plane_center, SAT_Vec3 plane_normal) {
    SAT_Vec3 *corners = state->significant_face[SAT_SIGNIFICANT_FACE_incident].corners;
    SAT_Vec3 u = sat_vec3(corners[1].x - corners[0].x, corners[1].y - corners[0].y, corners[1].z - corners[0].z);
    SAT_Vec3 v = sat_vec3(corners[3].x - corners[0].x, corners[3].y - corners[0].y, corners[3].z - corners[0].z);
    
    for(int i = 0; i < 4; ++i) {
        corners[i] = sat_clip_point_against_plane(corners[i], u, plane_center, plane_normal);
        corners[i] = sat_clip_point_against_plane(corners[i], v, plane_center, plane_normal);
    }
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
    state.center[SAT_A]         = lhs.center;
    state.unit_axis[SAT_A][0]   = sat_rotate(lhs.orientation, { 1, 0, 0 });
    state.unit_axis[SAT_A][1]   = sat_rotate(lhs.orientation, { 0, 1, 0 });
    state.unit_axis[SAT_A][2]   = sat_rotate(lhs.orientation, { 0, 0, 1 });
    state.scaled_axis[SAT_A][0] = sat_scaled_vec3(state.unit_axis[SAT_A][0], lhs.size[0]);
    state.scaled_axis[SAT_A][1] = sat_scaled_vec3(state.unit_axis[SAT_A][1], lhs.size[1]);
    state.scaled_axis[SAT_A][2] = sat_scaled_vec3(state.unit_axis[SAT_A][2], lhs.size[2]);
    
    state.center[SAT_B]         = rhs.center;
    state.unit_axis[SAT_B][0]   = sat_rotate(rhs.orientation, { 1, 0, 0 });
    state.unit_axis[SAT_B][1]   = sat_rotate(rhs.orientation, { 0, 1, 0 });
    state.unit_axis[SAT_B][2]   = sat_rotate(rhs.orientation, { 0, 0, 1 });
    state.scaled_axis[SAT_B][0] = sat_scaled_vec3(state.unit_axis[SAT_B][0], rhs.size[0]);
    state.scaled_axis[SAT_B][1] = sat_scaled_vec3(state.unit_axis[SAT_B][1], rhs.size[1]);
    state.scaled_axis[SAT_B][2] = sat_scaled_vec3(state.unit_axis[SAT_B][2], rhs.size[2]);
    
    state.lhs_to_rhs = sat_vec3(state.center[SAT_B].x - state.center[SAT_A].x,
                                state.center[SAT_B].y - state.center[SAT_A].y,
                                state.center[SAT_B].z - state.center[SAT_A].z);

    state.relative_velocity = { rhs.velocity.x - lhs.velocity.x, rhs.velocity.y - lhs.velocity.y, rhs.velocity.z - lhs.velocity.z };

    state.found_separating_axis = false;
    state.penetration_depth     = std::numeric_limits<SAT_Scalar>::max(); // This is so fucking stupid...
    state.collision_normal_dot_relative_velocity = std::numeric_limits<SAT_Scalar>::lowest();

    //
    // Check each axis to find a possible separating axis.
    // @@Speed: We could stop checking as soon as we found one separating axis, since a collision cannot
    // occur at that point.
    //
    for(int i = 0; i < SAT_AXIS_COUNT; ++i) {
        // Check each of the boxes axis.
        state.found_separating_axis |= sat_is_separating_axis(&state, state.unit_axis[SAT_A][i]);
        state.found_separating_axis |= sat_is_separating_axis(&state, state.unit_axis[SAT_B][i]);

        // Check some edge-to-edge axis.
        state.found_separating_axis |= sat_is_separating_axis(&state, sat_normalized_cross(state.unit_axis[SAT_A][i], state.unit_axis[SAT_B][i]));
    }
        
    // Check the remaining edge-to-edge axis.
    state.found_separating_axis |= sat_is_separating_axis(&state, sat_normalized_cross(state.unit_axis[SAT_B][SAT_AXIS_X], state.unit_axis[SAT_A][SAT_AXIS_Y]));
    state.found_separating_axis |= sat_is_separating_axis(&state, sat_normalized_cross(state.unit_axis[SAT_B][SAT_AXIS_X], state.unit_axis[SAT_A][SAT_AXIS_Z]));
    state.found_separating_axis |= sat_is_separating_axis(&state, sat_normalized_cross(state.unit_axis[SAT_B][SAT_AXIS_Y], state.unit_axis[SAT_A][SAT_AXIS_Z]));

    state.found_separating_axis |= sat_is_separating_axis(&state, sat_normalized_cross(state.unit_axis[SAT_B][SAT_AXIS_Y], state.unit_axis[SAT_A][SAT_AXIS_X]));
    state.found_separating_axis |= sat_is_separating_axis(&state, sat_normalized_cross(state.unit_axis[SAT_B][SAT_AXIS_Z], state.unit_axis[SAT_A][SAT_AXIS_X]));
    state.found_separating_axis |= sat_is_separating_axis(&state, sat_normalized_cross(state.unit_axis[SAT_B][SAT_AXIS_Z], state.unit_axis[SAT_A][SAT_AXIS_Y]));

    //
    // If we have not found a separating axis, it means that there is a collision.
    // Fill in the collision data.
    //
    if(!state.found_separating_axis) {
        //
        // Find the significant faces of both boxes. Each Box consists of 6 faces, the normals being the
        // axis both in positive and in negative direction.
        //

        state.significant_face[SAT_SIGNIFICANT_FACE_reference].face_normal_dot_collision_normal = 0;
        state.significant_face[SAT_SIGNIFICANT_FACE_incident].face_normal_dot_collision_normal  = 0;

        for(int i = 0; i < SAT_BOX_COUNT; ++i) {
            for(int j = 0; j < SAT_AXIS_COUNT; ++j) {
                int j_plus_1 = (j + 1) % SAT_AXIS_COUNT;
                int j_plus_2 = (j + 2) % SAT_AXIS_COUNT;

                SAT_Vec3 normal        = state.unit_axis[i][j];
                SAT_Vec3 scaled_normal = state.scaled_axis[i][j];
                SAT_Vec3 scaled_u      = state.scaled_axis[i][j_plus_1];
                SAT_Vec3 scaled_v      = state.scaled_axis[i][j_plus_2];
                
                sat_find_significant_face(&state, (SAT_Significant_Face) i, normal, scaled_normal, scaled_u, scaled_v, state.center[i]);
                sat_find_significant_face(&state, (SAT_Significant_Face) i, sat_negate(normal), sat_negate(scaled_normal), scaled_u, scaled_v, state.center[i]);
            }
        }
        
        SAT_Box_Index clipping_box = SAT_B;

        if(state.significant_face[SAT_SIGNIFICANT_FACE_incident].face_normal_dot_collision_normal > state.significant_face[SAT_SIGNIFICANT_FACE_reference].face_normal_dot_collision_normal) {
            // If the current incident face is more parallel to the collision normal than
            // the reference face, swap these two for better collision resolution.
            SAT_Face tmp = state.significant_face[SAT_SIGNIFICANT_FACE_incident];
            state.significant_face[SAT_SIGNIFICANT_FACE_incident]  = state.significant_face[SAT_SIGNIFICANT_FACE_reference];
            state.significant_face[SAT_SIGNIFICANT_FACE_reference] = tmp;
            clipping_box = SAT_A;
        }

        //
        // Do the Sutherland-Hodgman Clipping to find the polygon which actually overlaps the two boxes.
        // This polygon then determines the contact points.
        // 
        
        // 
        // We look at all the adjacent faces of the reference face, and check if any of the incident corners
        // lay outside of the planes built from the adjacent faces. If so, that corner point gets moved to
        // be on the edge of the face plane.
        //      |    |-----------
        //      |****|****|
        //      ^    |
        //      |    |-----------
        //      This corner needs to be clipped towards the right.


        for(int i = 0; i < SAT_AXIS_COUNT; ++i) {
            SAT_Vec3 positive_plane_normal = state.unit_axis[clipping_box][i];
            SAT_Vec3 positive_plane_center = { state.center[clipping_box].x + state.scaled_axis[clipping_box][i].x,
                                               state.center[clipping_box].y + state.scaled_axis[clipping_box][i].y,
                                               state.center[clipping_box].z + state.scaled_axis[clipping_box][i].z };

            SAT_Vec3 negative_plane_normal = sat_negate(state.unit_axis[clipping_box][i]);
            SAT_Vec3 negative_plane_center = { state.center[clipping_box].x - state.scaled_axis[clipping_box][i].x,
                                               state.center[clipping_box].y - state.scaled_axis[clipping_box][i].y,
                                               state.center[clipping_box].z - state.scaled_axis[clipping_box][i].z };

            sat_clip_incident_corners_against_plane(&state, positive_plane_center, positive_plane_normal);
            sat_clip_incident_corners_against_plane(&state, negative_plane_center, negative_plane_normal);
        }
        
        //
        // We finally clip all points against the reference face, so that corners of the incident face that
        // extend over the reference face are clipped as well:
        //     /*********/  < These corners need to be clipped down.
        // ---/---------/
        //   /*********/
        //

        {
            SAT_Vec3 plane_center = state.significant_face[SAT_SIGNIFICANT_FACE_reference].face_center;
            SAT_Vec3 plane_normal = state.significant_face[SAT_SIGNIFICANT_FACE_reference].face_normal;
            sat_clip_incident_corners_against_plane(&state, plane_center, plane_normal);
        }

        //
        // Fill in the result state.
        //
        
        result.found_collision = true;
        result.depth  = state.penetration_depth;
        result.normal = state.significant_face[SAT_SIGNIFICANT_FACE_reference].face_normal;
        
        // Ensure the collision normal always goes from RHS to LHS
        if(sat_dot(state.lhs_to_rhs, result.normal) > 0.0) result.normal = sat_negate(result.normal);
    
#if false
        result.world_space_position_count = 4; // @Incomplete: When an edge collides with a face, we only have two points. When a corner collides with a face, we only have one point.
        for(int i = 0; i < 4; ++i) result.world_space_positions[i] = state.significant_face[SAT_SIGNIFICANT_FACE_incident].corners[i];
#else
        //
        // Find the one average position of the contact area for a simpler collision resolution.
        //
        result.world_space_position_count = 1;
        result.world_space_positions[0] = sat_vec3(0, 0, 0);
        for(int i = 0; i < 4; ++i) {
            result.world_space_positions[0].x += state.significant_face[SAT_SIGNIFICANT_FACE_incident].corners[i].x / 4;
            result.world_space_positions[0].y += state.significant_face[SAT_SIGNIFICANT_FACE_incident].corners[i].y / 4;
            result.world_space_positions[0].z += state.significant_face[SAT_SIGNIFICANT_FACE_incident].corners[i].z / 4;
        }
#endif
    }
    
    //
    // Done.
    //
    
    return result;
}
