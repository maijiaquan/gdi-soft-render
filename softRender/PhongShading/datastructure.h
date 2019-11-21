#pragma once

//#include "frame.h"
#include <math.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <windows.h> // include important windows stuff
#include <windowsx.h>
#include <mmsystem.h>
#include <objbase.h>
#include <iostream> // include important C/C++ stuff
#include <conio.h>
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <io.h>
#include <fcntl.h>
#include <direct.h>
#include <wchar.h>

#include "device.h"
// defines for objects version 1

// used to compute the min and max of two expresions
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define SWAP(a, b, t) \
	{                 \
		t = a;        \
		a = b;        \
		b = t;        \
	}

// transformation control flags
#define TRANSFORM_LOCAL_ONLY 0
#define TRANSFORM_TRANS_ONLY 1
#define TRANSFORM_LOCAL_TO_TRANS 2

#define OBJECT4DV1_MAX_VERTICES 1024   // 64
#define OBJECT4DV1_MAX_POLYS 1024	  // 128
#define RENDERLIST4DV1_MAX_POLYS 32768 // 16384
#define PI ((float)3.141592654f)
#define DEG_TO_RAD(ang) ((ang)*PI / 180.0)

// states of polygons and faces
#define POLY4DV1_STATE_ACTIVE 0x0001
#define POLY4DV1_STATE_CLIPPED 0x0002
#define POLY4DV1_STATE_BACKFACE 0x0004

// states for objects
#define OBJECT4DV1_STATE_ACTIVE 0x0001
#define OBJECT4DV1_STATE_VISIBLE 0x0002
#define OBJECT4DV1_STATE_CULLED 0x0004

// shading mode of polygon
#define PLX_SHADE_MODE_PURE_FLAG 0x0000		 // this poly is a constant color
#define PLX_SHADE_MODE_CONSTANT_FLAG 0x0000  // alias
#define PLX_SHADE_MODE_FLAT_FLAG 0x2000		 // this poly uses flat shading
#define PLX_SHADE_MODE_GOURAUD_FLAG 0x4000   // this poly used gouraud shading
#define PLX_SHADE_MODE_PHONG_FLAG 0x6000	 // this poly uses phong shading
#define PLX_SHADE_MODE_FASTPHONG_FLAG 0x6000 // this poly uses phong shading (alias)

// attributes of polygons and polygon faces
#define POLY4DV1_ATTR_2SIDED 0x0001
#define POLY4DV1_ATTR_TRANSPARENT 0x0002
#define POLY4DV1_ATTR_8BITCOLOR 0x0004
#define POLY4DV1_ATTR_RGB16 0x0008
#define POLY4DV1_ATTR_RGB24 0x0010

#define POLY4DV1_ATTR_SHADE_MODE_PURE 0x0020
#define POLY4DV1_ATTR_SHADE_MODE_CONSTANT 0x0020 // (alias)
#define POLY4DV1_ATTR_SHADE_MODE_FLAT 0x0040
#define POLY4DV1_ATTR_SHADE_MODE_GOURAUD 0x0080
#define POLY4DV1_ATTR_SHADE_MODE_PHONG 0x0100
#define POLY4DV1_ATTR_SHADE_MODE_FASTPHONG 0x0100 // (alias)
#define POLY4DV1_ATTR_SHADE_MODE_TEXTURE 0x0200

#define PLX_2SIDED_FLAG 0x1000 // this poly is double sided
#define PLX_1SIDED_FLAG 0x0000 // this poly is single sided

#define PLX_COLOR_MODE_RGB_FLAG 0x8000	 // this poly uses RGB color
#define PLX_COLOR_MODE_INDEXED_FLAG 0x0000 // this poly uses an indexed 8-bit color

#define PLX_RGB_MASK 0x8000		   // mask to extract RGB or indexed color
#define PLX_SHADE_MODE_MASK 0x6000 // mask to extract shading mode
#define PLX_2SIDED_MASK 0x1000	 // mask for double sided
#define PLX_COLOR_MASK 0x0fff	  // xxxxrrrrggggbbbb, 4-bits per channel RGB \
                                   // xxxxxxxxiiiiiiii, indexed mode 8-bit index

#define RESET_BIT(word, bit_flag) ((word) = ((word) & (~bit_flag)))

// 相机旋转顺序
#define CAM_ROT_SEQ_XYZ 0
#define CAM_ROT_SEQ_YXZ 1
#define CAM_ROT_SEQ_XZY 2
#define CAM_ROT_SEQ_YZX 3
#define CAM_ROT_SEQ_ZYX 4
#define CAM_ROT_SEQ_ZXY 5

USHORT RGB16Bit565(int r, int g, int b);
#define _RGB16BIT565(r, g, b) ((b & 31) + ((g & 63) << 5) + ((r & 31) << 11))

// storage for our lookup tables
// float cos_look[361]; // 1 extra element so we can store 0-360 inclusive
// float sin_look[361];

// extern float cos_look[361]; // 1 extra element so we can store 0-360 inclusive
// extern float sin_look[361];

// extern float *cos_look; // 1 extra element so we can store 0-360 inclusive
// extern float *sin_look;

// float *cos_look;
// float *sin_look;

typedef unsigned short USHORT;

#define POLY4DV1_STATE_ACTIVE 0x0001
#define CAM_MODEL_EULER 0x0008
#define CAM_MODEL_UVN 0x0010

// general culling flags
#define CULL_OBJECT_X_PLANE 0x0001 // cull on the x clipping planes
#define CULL_OBJECT_Y_PLANE 0x0002 // cull on the y clipping planes
#define CULL_OBJECT_Z_PLANE 0x0004 // cull on the z clipping planes
#define CULL_OBJECT_XYZ_PLANES (CULL_OBJECT_X_PLANE | CULL_OBJECT_Y_PLANE | CULL_OBJECT_Z_PLANE)

#define RAND_RANGE(x, y) ((x) + (rand() % ((y) - (x) + 1)))
float Fast_Sin(float theta);
float Fast_Cos(float theta);

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

// 基于外部顶点列表的多边形
typedef struct POLY4DV1_TYP
{
	int state; // 状态
	int attr;  // 属性
	int color; // 颜色

	POINT4D_PTR vlist; // 顶点指针
	int vert[3];	   // 三个顶点的索引
} POLY4DV1, *POLY4DV1_PTR;

//物体
typedef struct OBJECT4DV1_TYP
{
	int id;
	char name[64];
	int state;
	int attr;
	float avg_radius; // average radius of object used for collision detection
	float max_radius; // maximum radius of object

	POINT4D world_pos;

	VECTOR4D dir;

	VECTOR4D ux, uy, uz;

	int num_vertices;

	POINT4D vlist_local[OBJECT4DV1_MAX_VERTICES]; //局部顶点列表
	POINT4D vlist_trans[OBJECT4DV1_MAX_VERTICES]; //变换后的顶点列表

	int num_polys;
	POLY4DV1 plist[OBJECT4DV1_MAX_POLYS]; //多边形数组

} OBJECT4DV1, *OBJECT4DV1_PTR;

//渲染列表，其实就是三角形数组
typedef struct RENDERLIST4DV1_TYP
{
	int state;
	int attr;

	POLYF4DV1_PTR poly_ptrs[RENDERLIST4DV1_MAX_POLYS]; //索引列表
	POLYF4DV1 poly_data[RENDERLIST4DV1_MAX_POLYS];

	int num_polys; // 三角形的数量
} RENDERLIST4DV1, *RENDERLIST4DV1_PTR;

// 3D 平面 ///////////////////////////////////////////////////
typedef struct PLANE3D_TYP
{
	POINT3D p0; // point on the plane
	VECTOR3D n; // normal to the plane (not necessarily a unit vector)
} PLANE3D, *PLANE3D_PTR;

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

//矩阵
// 4x4 matrix /////////////////////////////////////////////
typedef struct MATRIX4X4
{
	union {
		float M[4][4];
		struct
		{
			float M00, M01, M02, M03;
			float M10, M11, M12, M13;
			float M20, M21, M22, M23;
			float M30, M31, M32, M33;
		};
	};

} MATRIX4X4, *MATRIX4X4_PTR;

void Build_XYZ_Rotation_MATRIX4X4(float theta_x, // euler angles
								  float theta_y,
								  float theta_z,
								  MATRIX4X4_PTR mrot);

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

// 4x4 identity matrix
const MATRIX4X4 IMAT_4X4 = {1, 0, 0, 0,
							0, 1, 0, 0,
							0, 0, 1, 0,
							0, 0, 0, 1};

void MAT_IDENTITY_4X4(MATRIX4X4_PTR m);

void PLANE3D_Init(PLANE3D_PTR plane, POINT3D_PTR p0,
				  VECTOR3D_PTR normal, int normalize);

void VECTOR3D_Normalize(VECTOR3D_PTR va);
void VECTOR3D_Normalize(VECTOR3D_PTR va, VECTOR3D_PTR vn);
float VECTOR3D_Length(VECTOR3D_PTR va);
void VECTOR4D_Build(VECTOR4D_PTR init, VECTOR4D_PTR term, VECTOR4D_PTR result);

void Mat_Mul_VECTOR4D_4X4(VECTOR4D_PTR va, MATRIX4X4_PTR mb, VECTOR4D_PTR vprod);

void VECTOR4D_Add(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vsum);
VECTOR4D VECTOR4D_Add(VECTOR4D_PTR va, VECTOR4D_PTR vb);

float VECTOR4D_Dot(VECTOR4D_PTR va, VECTOR4D_PTR vb);
void VECTOR4D_Cross(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vn);
VECTOR4D VECTOR4D_Cross(VECTOR4D_PTR va, VECTOR4D_PTR vb);

//矩阵
//相机
typedef struct CAM4DV1_TYP
{
	int state;
	int attr;

	POINT4D pos;  // 相机在世界坐标中的位置
	VECTOR4D dir; // 相机的注视方向

	VECTOR4D u;
	VECTOR4D v;
	VECTOR4D n;
	VECTOR4D target; // 目标位置

	float view_dist; // 视距

	float fov; // 水平方向和垂直方向的视野

	float near_clip_z; // 近剪裁面的z
	float far_clip_z;  // 远剪裁面的z

	PLANE3D rt_clip_plane; // 右剪裁面
	PLANE3D lt_clip_plane; // 左剪裁面
	PLANE3D tp_clip_plane; // 上剪裁面
	PLANE3D bt_clip_plane; // 下剪裁面

	float viewplane_width;  // 视平面的宽
	float viewplane_height; // 视平面的高

	float viewport_width;	// 屏幕的宽
	float viewport_height;   // 屏幕的高
	float viewport_center_x; // 屏幕中心的x
	float viewport_center_y; // 屏幕中心的y
	float aspect_ratio;		 // 屏幕宽高比

	MATRIX4X4 mcam; // 世界坐标到相机坐标的变换矩阵
	MATRIX4X4 mper; // 相机坐标到透视坐标的变换矩阵
	MATRIX4X4 mscr; // 透视坐标到屏幕坐标的变换矩阵

} CAM4DV1, *CAM4DV1_PTR;

// RGB+alpha color
typedef struct RGBAV1_TYP
{
	union {
		int rgba;		 // compressed format
		UCHAR rgba_M[4]; // array format
		struct
		{
			UCHAR a, b, g, r;
		}; // explict name format
	};	 // end union

} RGBAV1, *RGBAV1_PTR;

// first light structure
typedef struct LIGHTV1_TYP
{
	int state; // state of light
	int id;	// id of light
	int attr;  // type of light, and extra qualifiers

	RGBAV1 c_ambient;  // ambient light intensity
	RGBAV1 c_diffuse;  // diffuse light intensity
	RGBAV1 c_specular; // specular light intensity

	POINT4D pos;	  // position of light
	VECTOR4D dir;	 // direction of light
	float kc, kl, kq; // attenuation factors
	float spot_inner; // inner angle for spot light
	float spot_outer; // outer angle for spot light
	float pf;		  // power factor/falloff for spot lights

	int iaux1, iaux2; // auxiliary vars for future expansion
	float faux1, faux2;
	void *ptr;

} LIGHTV1, *LIGHTV1_PTR;

// void Reset_RENDERLIST4DV1(RENDERLIST4DV1_PTR renderList);
void Init_CAM4DV1(CAM4DV1_PTR cam, int attr, POINT4D_PTR cam_pos, VECTOR4D_PTR cam_dir, VECTOR4D_PTR cam_target, float near_clip_z, float far_clip_z, float fov, float viewport_width, float viewport_height);

int Insert_POLYF4DV1_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, POLYF4DV1_PTR poly);
void Build_CAM4DV1_Matrix_Euler(CAM4DV1_PTR cam, int cam_rot_seq);

//加载plg格式文件
int Load_OBJECT4DV1_PLG(OBJECT4DV1_PTR obj, char *filename, VECTOR4D_PTR scale, VECTOR4D_PTR pos, VECTOR4D_PTR rot);
char *Get_Line_PLG(char *buffer, int maxlength, FILE *fp);
float Compute_OBJECT4DV1_Radius(OBJECT4DV1_PTR obj);

#define SET_BIT(word, bit_flag) ((word) = ((word) | (bit_flag)))

void Build_Sin_Cos_Tables(void);

//光照宏定义
// create some constants for ease of access
#define AMBIENT_LIGHT_INDEX 0  // ambient light index
#define INFINITE_LIGHT_INDEX 1 // infinite light index
#define POINT_LIGHT_INDEX 2	// point light index
#define SPOT_LIGHT_INDEX 3	 // spot light index
#define SPOT_LIGHT1_INDEX 4	// point light index
#define SPOT_LIGHT2_INDEX 3	// spot light index

#define LIGHTV1_STATE_ON 1  // light on
#define LIGHTV1_STATE_OFF 0 // light off

// defines for light types
#define LIGHTV1_ATTR_AMBIENT 0x0001		// basic ambient light
#define LIGHTV1_ATTR_INFINITE 0x0002	// infinite light source
#define LIGHTV1_ATTR_DIRECTIONAL 0x0002 // infinite light source (alias)
#define LIGHTV1_ATTR_POINT 0x0004		// point light source
#define LIGHTV1_ATTR_SPOTLIGHT1 0x0008  // spotlight type 1 (simple)
#define LIGHTV1_ATTR_SPOTLIGHT2 0x0010  // spotlight type 2 (complex)
//光照计算
int Light_RENDERLIST4DV1_World16(RENDERLIST4DV1_PTR rend_list, // list to process
								 CAM4DV1_PTR cam,			   // camera position
								 LIGHTV1_PTR lights,		   // light list (might have more than one)
								 int max_lights);			   // maximum lights in list

#define MAX_LIGHTS 8 // good luck with 1!

int Reset_Lights_LIGHTV1(void);

// lighting system
int Init_Light_LIGHTV1(int index,		   // index of light to create (0..MAX_LIGHTS-1)
					   int _state,		   // state of light
					   int _attr,		   // type of light, and extra qualifiers
					   RGBAV1 _c_ambient,  // ambient light intensity
					   RGBAV1 _c_diffuse,  // diffuse light intensity
					   RGBAV1 _c_specular, // specular light intensity
					   POINT4D_PTR _pos,   // position of light
					   VECTOR4D_PTR _dir,  // direction of light
					   float _kc,		   // attenuation factors
					   float _kl,
					   float _kq,
					   float _spot_inner, // inner angle for spot light
					   float _spot_outer, // outer angle for spot light
					   float _pf);		  // power factor/falloff for spot lights

#define _RGBA32BIT(r, g, b, a) ((a) + ((b) << 8) + ((g) << 16) + ((r) << 24))

void VECTOR4D_Normalize(VECTOR4D_PTR va);
float VECTOR4D_Length_Fast(VECTOR4D_PTR va);
float Fast_Distance_3D(float x, float y, float z);

extern LIGHTV1 lights[MAX_LIGHTS]; // lights in system
extern int num_lights;			   // current number of lights

// flags for sorting algorithm
#define SORT_POLYLIST_AVGZ 0  // sorts on average of all vertices
#define SORT_POLYLIST_NEARZ 1 // sorts on closest z vertex of each poly
#define SORT_POLYLIST_FARZ 2  // sorts on farthest z vertex of each poly
void Sort_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, int sort_method = SORT_POLYLIST_AVGZ);

// avg z-compare
int Compare_AvgZ_POLYF4DV1(const void *arg1, const void *arg2);
// near z-compare
int Compare_NearZ_POLYF4DV1(const void *arg1, const void *arg2);
// far z-compare
int Compare_FarZ_POLYF4DV1(const void *arg1, const void *arg2);

typedef struct VECTOR2D_TYP
{
	union {
		float M[2]; // array indexed storage

		// explicit names
		struct
		{
			float x, y;
		}; // end struct

	}; // end union

} VECTOR2D, POINT2D, *VECTOR2D_PTR, *POINT2D_PTR;

typedef struct BITMAP_IMAGE_TYP
{
	int state;		   // state of bitmap
	int attr;		   // attributes of bitmap
	int x, y;		   // position of bitmap
	int width, height; // size of bitmap
	int num_bytes;	 // total bytes of bitmap
	int bpp;		   // bits per pixel
	UCHAR *buffer;	 // pixels of bitmap

} BITMAP_IMAGE, *BITMAP_IMAGE_PTR;
//num_lights
// a first version of a "material"
typedef struct MATV1_TYP
{
	int state;	 // state of material
	int id;		   // id of this material, index into material array
	char name[64]; // name of material
	int attr;	  // attributes, the modes for shading, constant, flat,
				   // gouraud, fast phong, environment, textured etc.
				   // and other special flags...

	RGBAV1 color;			 // color of material
	float ka, kd, ks, power; // ambient, diffuse, specular,
							 // coefficients, note they are
							 // separate and scalars since many
							 // modelers use this format
							 // along with specular power

	RGBAV1 ra, rd, rs; // the reflectivities/colors pre-
					   // multiplied, to more match our
					   // definitions, each is basically
					   // computed by multiplying the
					   // color by the k's, eg:
					   // rd = color*kd etc.

	char texture_file[80]; // file location of texture
	BITMAP_IMAGE texture;  // actual texture map (if any)

	int iaux1, iaux2; // auxiliary vars for future expansion
	float faux1, faux2;
	void *ptr;

} MATV1, *MATV1_PTR;

typedef struct VERTEX4DTV1_TYP
{
	union {
		float M[12]; // array indexed storage

		// explicit names
		struct
		{
			float x, y, z, w;	 // point
			float nx, ny, nz, nw; // normal (vector or point)
			float u0, v0;		  // texture coordinates

			float i;  // final vertex intensity after lighting
			int attr; // attributes/ extra texture coordinates
		};			  // end struct

		// high level types
		struct
		{
			POINT4D v;  // the vertex
			VECTOR4D n; // the normal
			POINT2D t;  // texture coordinates
		};

	}; // end union

} VERTEX4DTV1, *VERTEX4DTV1_PTR;

// a polygon ver 2.0 based on an external vertex list  //////////////////////////////////
typedef struct POLY4DV2_TYP
{
	int state;		  // state information
	int attr;		  // physical attributes of polygon
	int color;		  // color of polygon
	int lit_color[3]; // holds colors after lighting, 0 for flat shading
					  // 0,1,2 for vertex colors after vertex lighting

	BITMAP_IMAGE_PTR texture; // pointer to the texture information for simple texture mapping

	int mati; // material index (-1) no material (new)

	VERTEX4DTV1_PTR vlist; // the vertex list itself
	POINT2D_PTR tlist;	 // the texture list itself (new)
	int vert[3];		   // the indices into the vertex list
	int text[3];		   // the indices into the texture coordinate list (new)
	float nlength;		   // length of normal (new)

} POLY4DV2, *POLY4DV2_PTR;

typedef struct OBJECT4DV2_TYP
{
	int id;			   // numeric id of this object
	char name[64];	 // ASCII name of object just for kicks
	int state;		   // state of object
	int attr;		   // attributes of object
	int mati;		   // material index overide (-1) - no material (new)
	float *avg_radius; // [OBJECT4DV2_MAX_FRAMES];   // average radius of object used for collision detection
	float *max_radius; // [OBJECT4DV2_MAX_FRAMES];   // maximum radius of object

	POINT4D world_pos; // position of object in world

	VECTOR4D dir; // rotation angles of object in local
				  // cords or unit direction vector user defined???

	VECTOR4D ux, uy, uz; // local axes to track full orientation
						 // this is updated automatically during
						 // rotation calls

	int num_vertices;   // number of vertices per frame of this object
	int num_frames;		// number of frames
	int total_vertices; // total vertices, redudant, but it saves a multiply in a lot of places
	int curr_frame;		// current animation frame (0) if single frame

	VERTEX4DTV1_PTR vlist_local; // [OBJECT4DV1_MAX_VERTICES]; // array of local vertices
	VERTEX4DTV1_PTR vlist_trans; // [OBJECT4DV1_MAX_VERTICES]; // array of transformed vertices

	// these are needed to track the "head" of the vertex list for mult-frame objects
	VERTEX4DTV1_PTR head_vlist_local;
	VERTEX4DTV1_PTR head_vlist_trans;

	// texture coordinates list (new)
	POINT2D_PTR tlist; // 3*num polys at max

	BITMAP_IMAGE_PTR texture; // pointer to the texture information for simple texture mapping (new)

	int num_polys;		// number of polygons in object mesh
	POLY4DV2_PTR plist; // ptr to polygons (new)

	int ivar1, ivar2;   // auxiliary vars
	float fvar1, fvar2; // auxiliary vars

	// METHODS //////////////////////////////////////////////////

	// setting the frame is so important that it should be a member function
	// calling functions without doing this can wreak havok!
	int Set_Frame(int frame);

} OBJECT4DV2, *OBJECT4DV2_PTR;

// a 2D vertex
typedef struct VERTEX2DF_TYP
{
	float x, y; // the vertex
} VERTEX2DF, *VERTEX2DF_PTR;

#define PARSER_DEBUG_OFF // enables/disables conditional compilation

#define PARSER_STRIP_EMPTY_LINES 1 // strips all blank lines
#define PARSER_LEAVE_EMPTY_LINES 2 // leaves empty lines
#define PARSER_STRIP_WS_ENDS 4	 // strips ws space at ends of line
#define PARSER_LEAVE_WS_ENDS 8	 // leaves it
#define PARSER_STRIP_COMMENTS 16   // strips comments out
#define PARSER_LEAVE_COMMENTS 32   // leaves comments in

#define PARSER_BUFFER_SIZE 256 // size of parser line buffer
#define PARSER_MAX_COMMENT 16  // maximum size of comment delimeter string

#define PARSER_DEFAULT_COMMENT "#" // default comment string for parser

// pattern language
#define PATTERN_TOKEN_FLOAT 'f'
#define PATTERN_TOKEN_INT 'i'
#define PATTERN_TOKEN_STRING 's'
#define PATTERN_TOKEN_LITERAL '\''

// state machine defines for pattern matching
#define PATTERN_STATE_INIT 0

#define PATTERN_STATE_RESTART 1
#define PATTERN_STATE_FLOAT 2
#define PATTERN_STATE_INT 3
#define PATTERN_STATE_LITERAL 4
#define PATTERN_STATE_STRING 5
#define PATTERN_STATE_NEXT 6

#define PATTERN_STATE_MATCH 7
#define PATTERN_STATE_END 8

#define PATTERN_MAX_ARGS 16
#define PATTERN_BUFFER_SIZE 80

// parser class ///////////////////////////////////////////////
class CPARSERV1
{
public:
	// constructor /////////////////////////////////////////////////
	CPARSERV1();

	// destructor ///////////////////////////////////////////////////
	~CPARSERV1();

	// reset file system ////////////////////////////////////////////
	int Reset();

	// open file /////////////////////////////////////////////////////
	int Open(char *filename);

	// close file ////////////////////////////////////////////////////
	int Close();

	// get line //////////////////////////////////////////////////////
	char *Getline(int mode);

	// sets the comment string ///////////////////////////////////////
	int SetComment(char *string);

	// find pattern in line //////////////////////////////////////////
	int Pattern_Match(char *string, char *pattern, ...);

	// VARIABLE DECLARATIONS /////////////////////////////////////////

public:
	FILE *fstream;					  // file pointer
	char buffer[PARSER_BUFFER_SIZE];  // line buffer
	int length;						  // length of current line
	int num_lines;					  // number of lines processed
	char comment[PARSER_MAX_COMMENT]; // single line comment string

	// pattern matching parameter storage, easier that variable arguments
	// anything matched will be stored here on exit from the call to pattern()
	char pstrings[PATTERN_MAX_ARGS][PATTERN_BUFFER_SIZE]; // any strings
	int num_pstrings;

	float pfloats[PATTERN_MAX_ARGS]; // any floats
	int num_pfloats;

	int pints[PATTERN_MAX_ARGS]; // any ints
	int num_pints;

}; // end CLASS CPARSERV1 //////////////////////////////////////////////

int Load_OBJECT4DV2_COB(OBJECT4DV2_PTR obj,	// pointer to object
						char *filename,		   // filename of Caligari COB file
						VECTOR4D_PTR scale,	// initial scaling factors
						VECTOR4D_PTR pos,	  // initial position
						VECTOR4D_PTR rot,	  // initial rotations
						int vertex_flags = 0); // flags to re-order vertices

// defines for objects version 2
// objects use dynamic allocation now, but keep as max values
#define OBJECT4DV2_MAX_VERTICES 4096 // 64
#define OBJECT4DV2_MAX_POLYS 8192	// 128

// states for objects
#define OBJECT4DV2_STATE_NULL 0x0000
#define OBJECT4DV2_STATE_ACTIVE 0x0001
#define OBJECT4DV2_STATE_VISIBLE 0x0002
#define OBJECT4DV2_STATE_CULLED 0x0004

// new
#define OBJECT4DV2_ATTR_SINGLE_FRAME 0x0001 // single frame object (emulates ver 1.0)
#define OBJECT4DV2_ATTR_MULTI_FRAME 0x0002  // multi frame object for .md2 support etc.
#define OBJECT4DV2_ATTR_TEXTURES 0x0004		// flags if object contains textured polys?

// render list defines ver 2.0
#define RENDERLIST4DV2_MAX_POLYS 32768

// defines for vertices, these are "hints" to the transform and
// lighting systems to help determine if a particular vertex has
// a valid normal that must be rotated, or a texture coordinate
// that must be clipped etc., this helps us minmize load during lighting
// and rendering since we can determine exactly what kind of vertex we
// are dealing with, something like a (direct3d) flexible vertex format in
// as much as it can hold:
// point
// point + normal
// point + normal + texture coordinates
#define VERTEX4DTV1_ATTR_NULL 0x0000 // this vertex is empty
#define VERTEX4DTV1_ATTR_POINT 0x0001
#define VERTEX4DTV1_ATTR_NORMAL 0x0002
#define VERTEX4DTV1_ATTR_TEXTURE 0x0004

// these are some defines for conditional compilation of the new rasterizers
// I don't want 80 million different functions, so I have decided to
// use some conditionals to change some of the logic in each
// these names aren't necessarily the most accurate, but 3 should be enough
#define RASTERIZER_ACCURATE 0 // sub-pixel accurate with fill convention
#define RASTERIZER_FAST 1	 //
#define RASTERIZER_FASTEST 2

// set this to the mode you want the engine to use
#define RASTERIZER_MODE RASTERIZER_ACCURATE

#define VERTEX_FLAGS_INVERT_X 0x0001 // inverts the Z-coordinates
#define VERTEX_FLAGS_INVERT_Y 0x0002 // inverts the Z-coordinates
#define VERTEX_FLAGS_INVERT_Z 0x0004 // inverts the Z-coordinates
#define VERTEX_FLAGS_SWAP_YZ 0x0008  // transforms a RHS model to a LHS model
#define VERTEX_FLAGS_SWAP_XZ 0x0010
#define VERTEX_FLAGS_SWAP_XY 0x0020
#define VERTEX_FLAGS_INVERT_WINDING_ORDER 0x0040 // invert winding order from cw to ccw or ccw to cc

#define VERTEX_FLAGS_TRANSFORM_LOCAL 0x0200		  // if file format has local transform then do it!
#define VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD 0x0400 // if file format has local to world then do it!

void Print_Mat_4X4(MATRIX4X4_PTR ma, char *name);

int Init_OBJECT4DV2(OBJECT4DV2_PTR obj, // object to allocate
					int _num_vertices,
					int _num_polys,
					int _num_frames,
					int destroy = 0);

float Compute_OBJECT4DV2_Radius(OBJECT4DV2_PTR obj);

int ReplaceChars(char *string_in, char *string_out, char *replace_chars, char rep_char, int case_on = 1);

// states of polygons and faces
#define POLY4DV2_STATE_NULL 0x0000
#define POLY4DV2_STATE_ACTIVE 0x0001
#define POLY4DV2_STATE_CLIPPED 0x0002
#define POLY4DV2_STATE_BACKFACE 0x0004
#define POLY4DV2_STATE_LIT 0x0008

#define POLY4DV2_ATTR_2SIDED 0x0001
#define POLY4DV2_ATTR_TRANSPARENT 0x0002
#define POLY4DV2_ATTR_8BITCOLOR 0x0004
#define POLY4DV2_ATTR_RGB16 0x0008
#define POLY4DV2_ATTR_RGB24 0x0010

#define POLY4DV2_ATTR_SHADE_MODE_PURE 0x0020
#define POLY4DV2_ATTR_SHADE_MODE_CONSTANT 0x0020 // (alias)
#define POLY4DV2_ATTR_SHADE_MODE_EMISSIVE 0x0020 // (alias)

#define POLY4DV2_ATTR_SHADE_MODE_FLAT 0x0040
#define POLY4DV2_ATTR_SHADE_MODE_GOURAUD 0x0080
#define POLY4DV2_ATTR_SHADE_MODE_PHONG 0x0100
#define POLY4DV2_ATTR_SHADE_MODE_FASTPHONG 0x0100 // (alias)
#define POLY4DV2_ATTR_SHADE_MODE_TEXTURE 0x0200

// new
#define POLY4DV2_ATTR_ENABLE_MATERIAL 0x0800  // use a real material for lighting
#define POLY4DV2_ATTR_DISABLE_MATERIAL 0x1000 // use basic color only for lighting (emulate version 1.0)

#define MAX_MATERIALS 256
extern MATV1 materials[MAX_MATERIALS]; // materials in system
extern int num_materials;			   // current number of materials

// defines for materials, follow our polygon attributes as much as possible
#define MATV1_ATTR_2SIDED 0x0001
#define MATV1_ATTR_TRANSPARENT 0x0002
#define MATV1_ATTR_8BITCOLOR 0x0004
#define MATV1_ATTR_RGB16 0x0008
#define MATV1_ATTR_RGB24 0x0010

#define MATV1_ATTR_SHADE_MODE_CONSTANT 0x0020
#define MATV1_ATTR_SHADE_MODE_EMMISIVE 0x0020 // alias
#define MATV1_ATTR_SHADE_MODE_FLAT 0x0040
#define MATV1_ATTR_SHADE_MODE_GOURAUD 0x0080
#define MATV1_ATTR_SHADE_MODE_FASTPHONG 0x0100
#define MATV1_ATTR_SHADE_MODE_TEXTURE 0x0200

char *Extract_Filename_From_Path(char *filepath, char *filename);

extern char texture_path[80]; // root path to ALL textures, make current directory for now

// container structure for bitmaps .BMP file
typedef struct BITMAP_FILE_TAG
{
	BITMAPFILEHEADER bitmapfileheader; // this contains the bitmapfile header
	BITMAPINFOHEADER bitmapinfoheader; // this is all the info including the palette
	PALETTEENTRY palette[256];		   // we will store the palette here
	UCHAR *buffer;					   // this is a pointer to the data

} BITMAP_FILE, *BITMAP_FILE_PTR;

extern BITMAP_FILE bitmap16bit; // a 16 bit bitmap file
int Load_Bitmap_File(BITMAP_FILE_PTR bitmap, char *filename);
int Create_Bitmap(BITMAP_IMAGE_PTR image, int x, int y, int width, int height, int bpp = 8);
int Load_Image_Bitmap16(BITMAP_IMAGE_PTR image, BITMAP_FILE_PTR bitmap, int cx, int cy, int mode);

int Load_Image_Bitmap(BITMAP_IMAGE_PTR image, BITMAP_FILE_PTR bitmap, int cx, int cy, int mode);
#define BITMAP_EXTRACT_MODE_CELL 0
#define BITMAP_EXTRACT_MODE_ABS 1

int Unload_Bitmap_File(BITMAP_FILE_PTR bitmap);

extern int screen_width, // width of screen
	screen_height,		 // height of screen
	screen_bpp,			 // bits per pixel
	screen_windowed;	 // is this a windowed app?

#define SCREEN_WIDTH 600 // size of screen
#define SCREEN_HEIGHT 600
#define SCREEN_BPP 8 // bits per pixel
int RGBto8BitIndex(UCHAR r, UCHAR g, UCHAR b, LPPALETTEENTRY palette, int flush_cache);
extern PALETTEENTRY palette[256]; // color palette
#define MAX_COLORS_PALETTE 256
#define VERTEX_FLAGS_INVERT_TEXTURE_U 0x0080 // invert u texture coordinate
#define VERTEX_FLAGS_INVERT_TEXTURE_V 0x0100 // invert v texture coordinate
#define VERTEX_FLAGS_INVERT_SWAP_UV 0x0800   // swap u and v texture coordinates
int Compute_OBJECT4DV2_Poly_Normals(OBJECT4DV2_PTR obj);
int Compute_OBJECT4DV2_Vertex_Normals(OBJECT4DV2_PTR obj);
char *StringLtrim(char *string);
char *StringRtrim(char *string);

float IsFloat(char *fstring);

int IsInt(char *istring);
int Destroy_OBJECT4DV2(OBJECT4DV2_PTR obj);

// bitmap defines
#define BITMAP_ID 0x4D42 // universal id for a bitmap
#define BITMAP_STATE_DEAD 0
#define BITMAP_STATE_ALIVE 1
#define BITMAP_STATE_DYING 2
#define BITMAP_ATTR_LOADED 128

#define BITMAP_EXTRACT_MODE_CELL 0
#define BITMAP_EXTRACT_MODE_ABS 1
int Flip_Bitmap(UCHAR *image, int bytes_per_line, int height);
float VECTOR4D_Length(VECTOR4D_PTR va);

typedef struct POLYF4DV2_TYP
{
	int state;				  // state information
	int attr;				  // physical attributes of polygon
	int color;				  // color of polygon
	int lit_color[3];		  // holds colors after lighting, 0 for flat shading
							  // 0,1,2 for vertex colors after vertex lighting
	BITMAP_IMAGE_PTR texture; // pointer to the texture information for simple texture mapping

	int mati; // material index (-1) for no material  (new)

	float nlength;   // length of the polygon normal if not normalized (new)
	VECTOR4D normal; // the general polygon normal (new)

	float avg_z; // average z of vertices, used for simple sorting (new)

	VERTEX4DTV1 vlist[3];  // the vertices of this triangle
	VERTEX4DTV1 tvlist[3]; // the vertices after transformation if needed

	POLYF4DV2_TYP *next; // pointer to next polygon in list??
	POLYF4DV2_TYP *prev; // pointer to previous polygon in list??

} POLYF4DV2, *POLYF4DV2_PTR;

typedef struct RENDERLIST4DV2_TYP
{
	int state; // state of renderlist ???
	int attr;  // attributes of renderlist ???

	// the render list is an array of pointers each pointing to
	// a self contained "renderable" polygon face POLYF4DV2
	POLYF4DV2_PTR poly_ptrs[RENDERLIST4DV2_MAX_POLYS];

	// additionally to cut down on allocatation, de-allocation
	// of polygons each frame, here's where the actual polygon
	// faces will be stored
	POLYF4DV2 poly_data[RENDERLIST4DV2_MAX_POLYS];

	int num_polys; // number of polys in render list

} RENDERLIST4DV2, *RENDERLIST4DV2_PTR;
void Reset_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list);
void Reset_OBJECT4DV2(OBJECT4DV2_PTR obj);

void Transform_OBJECT4DV2(OBJECT4DV2_PTR obj, MATRIX4X4_PTR mt,
						  int coord_select, int transform_basis, int all_frames = 0);

void Model_To_World_OBJECT4DV2(OBJECT4DV2_PTR obj, int coord_select = TRANSFORM_LOCAL_TO_TRANS, int all_frames = 0);

int Insert_POLY4DV2_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
								   POLY4DV2_PTR poly);
int Insert_OBJECT4DV2_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
									 OBJECT4DV2_PTR obj,
									 int insert_local);

inline void VERTEX4DTV1_COPY(VERTEX4DTV1_PTR vdst, VERTEX4DTV1_PTR vsrc)
{
	*vdst = *vsrc;
}

void Remove_Backfaces_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list, CAM4DV1_PTR cam);

int Light_RENDERLIST4DV2_World16(RENDERLIST4DV2_PTR rend_list, // list to process
								 CAM4DV1_PTR cam,			   // camera position
								 LIGHTV1_PTR lights,		   // light list (might have more than one)
								 int max_lights);			   // maximum lights in list
void World_To_Camera_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
									CAM4DV1_PTR cam);
void Sort_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list, int sort_method);

void Camera_To_Perspective_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
										  CAM4DV1_PTR cam);

void Perspective_To_Screen_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
										  CAM4DV1_PTR cam);
// avg z-compare
int Compare_AvgZ_POLYF4DV2(const void *arg1, const void *arg2);

// near z-compare
int Compare_NearZ_POLYF4DV2(const void *arg1, const void *arg2);

// far z-compare
int Compare_FarZ_POLYF4DV2(const void *arg1, const void *arg2);

inline float VECTOR4D_Length_Fast2(VECTOR4D_PTR va)
{
	// this function computes the distance from the origin to x,y,z

	int temp;	// used for swaping
	int x, y, z; // used for algorithm

	// make sure values are all positive
	x = fabs(va->x) * 1024;
	y = fabs(va->y) * 1024;
	z = fabs(va->z) * 1024;

	// sort values
	if (y < x)
		SWAP(x, y, temp)
	if (z < y)
		SWAP(y, z, temp)
	if (y < x)
		SWAP(x, y, temp)

	int dist = (z + 11 * (y >> 5) + (x >> 2));

	// compute distance with 8% error
	return ((float)(dist >> 10));

} // end VECTOR4D_Length_Fast2

#define TRI_TYPE_NONE 0
#define TRI_TYPE_FLAT_TOP 1
#define TRI_TYPE_FLAT_BOTTOM 2
#define TRI_TYPE_FLAT_MASK 3
#define TRI_TYPE_GENERAL 4
#define INTERP_LHS 0
#define INTERP_RHS 1
#define MAX_VERTICES_PER_POLY 6

#define FIXP16_SHIFT 16
#define FIXP16_MAG 65536
#define FIXP16_DP_MASK 0x0000ffff
#define FIXP16_WP_MASK 0xffff0000
#define FIXP16_ROUND_UP 0x00008000
void Draw_Gouraud_Triangle16(device_t *device, POLYF4DV2_PTR face);
// void DrawPhongTriangle(device_t *device, POLYF4DV2_PTR face, LIGHTV1_PTR lights);
void DrawPhongTriangle(device_t *device, CAM4DV1_PTR ptrCam, POLYF4DV2_PTR face, POLYF4DV2_PTR faceInWorld, LIGHTV1_PTR lights);

void ComputePhongShadingPixelColor(int r_base, int g_base, int b_base, LIGHTV1_PTR prtlights, CAM4DV1_PTR ptrCam, VECTOR4D_PTR prtFragPos, VECTOR4D_PTR ptrFragNormal, IUINT32 &color);
