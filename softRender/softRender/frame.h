#pragma once
#define RENDERLIST4DV1_MAX_POLYS 32768 // 16384
#define POLY4DV1_STATE_ACTIVE 0x0001
#define CAM_MODEL_EULER 0x0008
#define CAM_MODEL_UVN 0x0010
#define PI         ((float)3.141592654f)

typedef unsigned short USHORT;
#define _RGB16BIT565(r, g, b) ((b & 31) + ((g & 63) << 5) + ((r & 31) << 11))

#include <math.h>


// storage for our lookup tables
float cos_look[361]; // 1 extra element so we can store 0-360 inclusive
float sin_look[361];
#define WINDOW_WIDTH 400 // size of window
#define WINDOW_HEIGHT 400


// 3D vector, point without the w ////////////////////////
typedef struct VECTOR3D_TYP
{
    union {
        float M[3]; // array indexed storage

        // explicit names
        struct
        {
            float x, y, z;
        }; // end struct

    }; // end union

} VECTOR3D, POINT3D, *VECTOR3D_PTR, *POINT3D_PTR;

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
    int state; // state information
    int attr;  // physical attributes of polygon
    int color; // color of polygon

    POINT4D vlist[3];  // the vertices of this triangle
    POINT4D tvlist[3]; // the vertices after transformation if needed

    POLYF4DV1_TYP *next; // pointer to next polygon in list??
    POLYF4DV1_TYP *prev; // pointer to previous polygon in list??

} POLYF4DV1, *POLYF4DV1_PTR;

//渲染列表，其实就是三角形数组
typedef struct RENDERLIST4DV1_TYP
{
    int state; // state of renderlist ???
    int attr;  // attributes of renderlist ???

    POLYF4DV1_PTR poly_ptrs[RENDERLIST4DV1_MAX_POLYS];
    POLYF4DV1 poly_data[RENDERLIST4DV1_MAX_POLYS];

    int num_polys; // number of polys in render list

} RENDERLIST4DV1, *RENDERLIST4DV1_PTR;

// 3D 平面 ///////////////////////////////////////////////////
typedef struct PLANE3D_TYP
{
    POINT3D p0; // point on the plane
    VECTOR3D n; // normal to the plane (not necessarily a unit vector)
} PLANE3D, *PLANE3D_PTR;

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



inline void VECTOR4D_COPY(VECTOR4D_PTR vdst, VECTOR4D_PTR vsrc) 
{(vdst)->x = (vsrc)->x; (vdst)->y = (vsrc)->y;  
(vdst)->z = (vsrc)->z; (vdst)->w = (vsrc)->w;  }

inline void VECTOR4D_INITXYZ(VECTOR4D_PTR v, float x,float y,float z) 
{(v)->x = (x); (v)->y = (y); (v)->z = (z); (v)->w = 1.0;}
USHORT RGB16Bit565(int r, int g, int b);

// function ptr to RGB16 builder
USHORT (*RGB16Bit)
(int r, int g, int b) = nullptr ;

void Reset_RENDERLIST4DV1(RENDERLIST4DV1_PTR renderList);
void Build_Sin_Cos_Tables(void);
void Init_CAM4DV1(CAM4DV1_PTR cam, int attr, POINT4D_PTR cam_pos,
                  VECTOR4D_PTR cam_dir, VECTOR4D_PTR cam_target,
                  float near_clip_z, float far_clip_z, float fov,
                  float viewport_width, float viewport_height);