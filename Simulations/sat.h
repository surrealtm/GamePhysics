//
// This file implements a collision check using SAT.
// It correctly handles the (very common case) where the  collision area is not only a single point (e.g. when
// two faces collide head-on), by generating a collection  of collision points, which then all cause a collision
// resolution.
// The following papers were taken as a guide:
// https://research.ncl.ac.uk/game/mastersdegree/gametechnologies/previousinformation/physics5collisionmanifolds/2017%20Tutorial%205%20-%20Collision%20Manifolds.pdf
// https://jkh.me/files/tutorials/Separating%20Axis%20Theorem%20for%20Oriented%20Bounding%20Boxes.pdf
//

#define SAT_USE_MINIMAL_TRANSLATION_VECTOR true // Resolve along the axis with the least penetration
//#define SAT_USE_PROJECTED_VELOCITY true // Resolve along the axis most parallel to the relative velocity

/* =============================== SAT API =============================== */

typedef double SAT_Scalar;

union SAT_Vec3 {
    struct {
        SAT_Scalar x;
        SAT_Scalar y;
        SAT_Scalar z;
    };

    SAT_Scalar _[3];

    SAT_Scalar operator[](int index) { return this->_[index]; }
};

struct SAT_Quat {
    SAT_Scalar x;
    SAT_Scalar y;
    SAT_Scalar z;
    SAT_Scalar w;
};

struct SAT_Input {
    SAT_Vec3 center;
    SAT_Quat orientation;
    SAT_Vec3 size;
    SAT_Vec3 velocity;
};

struct SAT_Result {
    bool found_collision;
    SAT_Scalar depth;
    SAT_Vec3 normal; // The normal from RHS towards LHS.
    SAT_Vec3 world_space_positions[4]; // A collision will always occur on an area in 3D space, and these positions are the corners of that area. All of them should cause a collision resolution, so that the whole area is reflected and not just a single point. Sometimes however, that area is so small that we represent it as either 2 points (e.g. edge-vs-face), or 1 point (e.g. corner-vs-face).
    int world_space_position_count;
};


SAT_Vec3 sat_vec3(SAT_Scalar x, SAT_Scalar y, SAT_Scalar z);
SAT_Quat sat_quat(SAT_Scalar x, SAT_Scalar y, SAT_Scalar z, SAT_Scalar w);
SAT_Result sat(SAT_Input lhs, SAT_Input rhs);



/* =============================== Internal State =============================== */

enum SAT_Box_Index {
    SAT_A = 0x0,
    SAT_B = 0x1,
    SAT_BOX_COUNT = 0x2,
};

enum SAT_Axis {
    SAT_AXIS_X = 0x0,
    SAT_AXIS_Y = 0x1,
    SAT_AXIS_Z = 0x2,
    SAT_AXIS_COUNT = 0x3,
};

enum SAT_Significant_Face {
    SAT_SIGNIFICANT_FACE_incident  = 0x0,
    SAT_SIGNIFICANT_FACE_reference = 0x1,
    SAT_SIGNIFICANT_FACE_COUNT     = 0x2,
};

struct SAT_Face {
    SAT_Vec3 face_normal;
    SAT_Vec3 face_center; // For clipping corners against this face.
    SAT_Scalar face_normal_dot_collision_normal;
    SAT_Vec3 corners[4];
};

struct SAT_State {
    SAT_Vec3 center[SAT_BOX_COUNT];
    SAT_Vec3 unit_axis[SAT_BOX_COUNT][SAT_AXIS_COUNT];
    SAT_Vec3 scaled_axis[SAT_BOX_COUNT][SAT_AXIS_COUNT];
    SAT_Vec3 lhs_to_rhs; // T := rhs.center - lhs.center
    SAT_Vec3 relative_velocity;

    SAT_Face significant_face[SAT_SIGNIFICANT_FACE_COUNT];
    
    bool found_separating_axis; // If this is set to true, then there is no collision.
    SAT_Scalar penetration_depth; // The smallest overlap value on any given axis.
    
    // We want the collision normal to be as close to the "impact velocity" as possible,
    // to improve the response in cases like this:
    // ┌──────┐              ┌──────┐
    // │      │ ---------→   │      │
    // │      │       ┌──────┼──────┼───────┐ ↑
    // └──────┘       │      └──────┘       │ ↓
    //                └─────────────────────┘
    // We want the collision normal to be the axis to the left, not the axis to the top,
    // since the object came from the left. The MTV (smallest penetration vector) would
    // move the object upwards.
    SAT_Scalar collision_normal_dot_relative_velocity;
    SAT_Vec3 collision_normal; // The axis which had a collision with the smallest overlap.
};

