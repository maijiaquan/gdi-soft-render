#pragma once

#include <math.h>
#include <string.h>
typedef unsigned short USHORT;

#define RENDERLIST4DV1_MAX_POLYS 32768 // 16384
#define POLY4DV1_STATE_ACTIVE 0x0001
#define CAM_MODEL_EULER 0x0008
#define CAM_MODEL_UVN 0x0010
#define PI ((float)3.141592654f)
#define _RGB16BIT565(r, g, b) ((b & 31) + ((g & 63) << 5) + ((r & 31) << 11))
#define DEG_TO_RAD(ang) ((ang)*PI / 180.0)

// defines for small numbers
#define EPSILON_E3 (float)(1E-3)
#define EPSILON_E4 (float)(1E-4)
#define EPSILON_E5 (float)(1E-5)
#define EPSILON_E6 (float)(1E-6)

#define WINDOW_WIDTH 400 // size of window
#define WINDOW_HEIGHT 400

// transformation control flags
#define TRANSFORM_LOCAL_ONLY 0 // perform the transformation in place on the \
                               // local/world vertex list
#define TRANSFORM_TRANS_ONLY 1 // perfrom the transformation in place on the \
                               // "transformed" vertex list

#define TRANSFORM_LOCAL_TO_TRANS 2

// defines for camera rotation sequences
#define CAM_ROT_SEQ_XYZ 0
#define CAM_ROT_SEQ_YXZ 1
#define CAM_ROT_SEQ_XZY 2
#define CAM_ROT_SEQ_YZX 3
#define CAM_ROT_SEQ_ZYX 4
#define CAM_ROT_SEQ_ZXY 5

// states of polygons and faces
#define POLY4DV1_STATE_ACTIVE 0x0001
#define POLY4DV1_STATE_CLIPPED 0x0002
#define POLY4DV1_STATE_BACKFACE 0x0004

// storage for our lookup tables
extern float cos_look[361]; // 1 extra so we can store 0-360 inclusive
extern float sin_look[361]; // 1 extra so we can store 0-360 inclusive

// 3D vector, point without the w ////////////////////////
typedef struct VECTOR3D_TYP
{
    union {
        float M[3];
        // explicit names
        struct
        {
            float x, y, z;
        };
    };

} VECTOR3D, POINT3D, *VECTOR3D_PTR, *POINT3D_PTR;

//顶点的结构
typedef struct VECTOR4D_TYP
{
    union {
        float M[4];

        struct
        {
            float x, y, z, w;
        };
    };
} VECTOR4D, POINT4D, *VECTOR4D_PTR, *POINT4D_PTR;

//Polygon多边形，其实就是一个三角形
typedef struct POLYF4DV1_TYP
{
    int state;
    int attr;
    int color;

    POINT4D vlist[3];  // 该三角形的三个顶点
    POINT4D tvlist[3]; // 该三角形变换后的顶点

    POLYF4DV1_TYP *next;
    POLYF4DV1_TYP *prev;
} POLYF4DV1, *POLYF4DV1_PTR;

//渲染列表，其实就是三角形数组
typedef struct RENDERLIST4DV1_TYP
{
    int state; 
    int attr;  

    POLYF4DV1_PTR poly_ptrs[RENDERLIST4DV1_MAX_POLYS];  //索引列表
    POLYF4DV1 poly_data[RENDERLIST4DV1_MAX_POLYS];

    int num_polys; // 三角形的数量
} RENDERLIST4DV1, *RENDERLIST4DV1_PTR;

// 3D 平面 ///////////////////////////////////////////////////
typedef struct PLANE3D_TYP
{
    POINT3D p0; // point on the plane
    VECTOR3D n; // normal to the plane (not necessarily a unit vector)
} PLANE3D, *PLANE3D_PTR;

//矩阵
// 4x4 matrix /////////////////////////////////////////////
typedef struct MATRIX4X4
{
    union {
        float M[4][4]; // array indexed data storage

        // storage in row major form with explicit names
        struct
        {
            float M00, M01, M02, M03;
            float M10, M11, M12, M13;
            float M20, M21, M22, M23;
            float M30, M31, M32, M33;
        }; // end explicit names

    }; // end union

} MATRIX4X4, *MATRIX4X4_PTR;

//矩阵

// 4x4 identity matrix
const MATRIX4X4 IMAT_4X4 = {1, 0, 0, 0,
                            0, 1, 0, 0,
                            0, 0, 1, 0,
                            0, 0, 0, 1};

// camera version 1
typedef struct CAM4DV1_TYP
{
    int state; // state of camera
    int attr;  // camera attributes

    POINT4D pos; // world position of camera used by both camera models

    VECTOR4D dir; // angles or look at direction of camera for simple
                  // euler camera models, elevation and heading for
                  // uvn model

    VECTOR4D u; // extra vectors to track the camera orientation
    VECTOR4D v; // for more complex UVN camera model
    VECTOR4D n;

    VECTOR4D target; // look at target

    float view_dist; // focal length

    float fov; // field of view for both horizontal and vertical axes

    // 3d clipping planes
    // if view volume is NOT 90 degree then general 3d clipping
    // must be employed
    float near_clip_z; // near z=constant clipping plane
    float far_clip_z;  // far z=constant clipping plane

    PLANE3D rt_clip_plane; // the right clipping plane
    PLANE3D lt_clip_plane; // the left clipping plane
    PLANE3D tp_clip_plane; // the top clipping plane
    PLANE3D bt_clip_plane; // the bottom clipping plane

    float viewplane_width;  // width and height of view plane to project onto
    float viewplane_height; // usually 2x2 for normalized projection or
                            // the exact same size as the viewport or screen window

    // remember screen and viewport are synonomous
    float viewport_width; // size of screen/viewport
    float viewport_height;
    float viewport_center_x; // center of view port (final image destination)
    float viewport_center_y;

    // aspect ratio
    float aspect_ratio;

    // these matrices are not necessarily needed based on the method of
    // transformation, for example, a manual perspective or screen transform
    // and or a concatenated perspective/screen, however, having these
    // matrices give us more flexibility

    MATRIX4X4 mcam; // storage for the world to camera transform matrix
    MATRIX4X4 mper; // storage for the camera to perspective transform matrix
    MATRIX4X4 mscr; // storage for the perspective to screen transform matrix

} CAM4DV1, *CAM4DV1_PTR;

//内联函数
inline void VECTOR3D_INITXYZ(VECTOR3D_PTR v, float x, float y, float z)
{
    (v)->x = (x);
    (v)->y = (y);
    (v)->z = (z);
}
inline void VECTOR4D_COPY(VECTOR4D_PTR vdst, VECTOR4D_PTR vsrc)
{
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
    (vdst)->z = (vsrc)->z;
    (vdst)->w = (vsrc)->w;
}

inline void VECTOR4D_INITXYZ(VECTOR4D_PTR v, float x, float y, float z)
{
    (v)->x = (x);
    (v)->y = (y);
    (v)->z = (z);
    (v)->w = 1.0;
}

inline void POINT3D_COPY(POINT3D_PTR vdst, POINT3D_PTR vsrc)
{
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
    (vdst)->z = (vsrc)->z;
}

inline void VECTOR4D_ZERO(VECTOR4D_PTR v)
{
    (v)->x = (v)->y = (v)->z = 0.0;
    (v)->w = 1.0;
}

inline void VECTOR3D_COPY(VECTOR3D_PTR vdst, VECTOR3D_PTR vsrc)
{
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
    (vdst)->z = (vsrc)->z;
}

inline void VECTOR3D_ZERO(VECTOR3D_PTR v)
{
    (v)->x = (v)->y = (v)->z = 0.0;
}

USHORT RGB16Bit565(int r, int g, int b);

void Reset_RENDERLIST4DV1(RENDERLIST4DV1_PTR renderList);
void Init_CAM4DV1(CAM4DV1_PTR cam, int attr, POINT4D_PTR cam_pos,
                  VECTOR4D_PTR cam_dir, VECTOR4D_PTR cam_target,
                  float near_clip_z, float far_clip_z, float fov,
                  float viewport_width, float viewport_height);

void MAT_IDENTITY_4X4(MATRIX4X4_PTR m);
void PLANE3D_Init(PLANE3D_PTR plane, POINT3D_PTR p0,
                  VECTOR3D_PTR normal, int normalize);

void VECTOR3D_Normalize(VECTOR3D_PTR va);
void VECTOR3D_Normalize(VECTOR3D_PTR va, VECTOR3D_PTR vn);
float VECTOR3D_Length(VECTOR3D_PTR va);

int Insert_POLYF4DV1_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, POLYF4DV1_PTR poly);

void Build_XYZ_Rotation_MATRIX4X4(float theta_x, // euler angles
                                  float theta_y,
                                  float theta_z,
                                  MATRIX4X4_PTR mrot);

float Fast_Sin(float theta);
float Fast_Cos(float theta);
void Mat_Mul_4X4(MATRIX4X4_PTR ma, MATRIX4X4_PTR mb, MATRIX4X4_PTR mprod);

void Mat_Init_4X4(MATRIX4X4_PTR ma,
                  float m00, float m01, float m02, float m03,
                  float m10, float m11, float m12, float m13,
                  float m20, float m21, float m22, float m23,
                  float m30, float m31, float m32, float m33);

#define MAT_COPY_4X4(src_mat, dest_mat)                                   \
    {                                                                     \
        memcpy((void *)(dest_mat), (void *)(src_mat), sizeof(MATRIX4X4)); \
    }

void Transform_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, // render list to transform
                              MATRIX4X4_PTR mt,             // transformation matrix
                              int coord_select);

void Model_To_World_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, POINT4D_PTR world_pos,
                                   int coord_select = TRANSFORM_LOCAL_TO_TRANS);

void Build_CAM4DV1_Matrix_Euler(CAM4DV1_PTR cam, int cam_rot_seq);

void World_To_Camera_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, CAM4DV1_PTR cam);

void Camera_To_Perspective_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, CAM4DV1_PTR cam);

void Perspective_To_Screen_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, CAM4DV1_PTR cam);

void Mat_Mul_VECTOR4D_4X4(VECTOR4D_PTR va, MATRIX4X4_PTR mb, VECTOR4D_PTR vprod);

void VECTOR4D_Add(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vsum);
VECTOR4D VECTOR4D_Add(VECTOR4D_PTR va, VECTOR4D_PTR vb);