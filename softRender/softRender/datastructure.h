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

// defines for objects version 1

// transformation control flags
#define TRANSFORM_LOCAL_ONLY 0
#define TRANSFORM_TRANS_ONLY 1
#define TRANSFORM_LOCAL_TO_TRANS 2

#define OBJECT4DV1_MAX_VERTICES 1024   // 64
#define OBJECT4DV1_MAX_POLYS 1024	  // 128
#define RENDERLIST4DV1_MAX_POLYS 32768 // 16384
#define PI ((float)3.141592654f)
#define DEG_TO_RAD(ang) ((ang)*PI / 180.0)

// defines for small numbers
#define EPSILON_E3 (float)(1E-3)
#define EPSILON_E4 (float)(1E-4)
#define EPSILON_E5 (float)(1E-5)
#define EPSILON_E6 (float)(1E-6)

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

// �����ת˳��
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

#define WINDOW_WIDTH 400 // size of window
#define WINDOW_HEIGHT 400

// general culling flags
#define CULL_OBJECT_X_PLANE 0x0001 // cull on the x clipping planes
#define CULL_OBJECT_Y_PLANE 0x0002 // cull on the y clipping planes
#define CULL_OBJECT_Z_PLANE 0x0004 // cull on the z clipping planes
#define CULL_OBJECT_XYZ_PLANES (CULL_OBJECT_X_PLANE | CULL_OBJECT_Y_PLANE | CULL_OBJECT_Z_PLANE)


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

//����Ľṹ
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

//Polygon����Σ���ʵ����һ��������
typedef struct POLYF4DV1_TYP
{
	int state;
	int attr;
	int color;

	POINT4D vlist[3];  // �������ε���������
	POINT4D tvlist[3]; // �������α任��Ķ���

	POLYF4DV1_TYP *next;
	POLYF4DV1_TYP *prev;
} POLYF4DV1, *POLYF4DV1_PTR;

// �����ⲿ�����б�Ķ����
typedef struct POLY4DV1_TYP
{
	int state; // ״̬
	int attr;  // ����
	int color; // ��ɫ

	POINT4D_PTR vlist; // ����ָ��
	int vert[3];	   // �������������
} POLY4DV1, *POLY4DV1_PTR;

//����
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

	POINT4D vlist_local[OBJECT4DV1_MAX_VERTICES];	//�ֲ������б�
	POINT4D vlist_trans[OBJECT4DV1_MAX_VERTICES];	//�任��Ķ����б�

	int num_polys;
	POLY4DV1 plist[OBJECT4DV1_MAX_POLYS];	//���������

} OBJECT4DV1, *OBJECT4DV1_PTR;

//��Ⱦ�б���ʵ��������������
typedef struct RENDERLIST4DV1_TYP
{
	int state;
	int attr;

	POLYF4DV1_PTR poly_ptrs[RENDERLIST4DV1_MAX_POLYS]; //�����б�
	POLYF4DV1 poly_data[RENDERLIST4DV1_MAX_POLYS];

	int num_polys; // �����ε�����
} RENDERLIST4DV1, *RENDERLIST4DV1_PTR;

// 3D ƽ�� ///////////////////////////////////////////////////
typedef struct PLANE3D_TYP
{
	POINT3D p0; // point on the plane
	VECTOR3D n; // normal to the plane (not necessarily a unit vector)
} PLANE3D, *PLANE3D_PTR;

//��������
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

//����
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
void VECTOR4D_Cross(VECTOR4D_PTR va,VECTOR4D_PTR vb,VECTOR4D_PTR vn);
VECTOR4D VECTOR4D_Cross(VECTOR4D_PTR va, VECTOR4D_PTR vb);


//����
//���
typedef struct CAM4DV1_TYP
{
	int state;
	int attr;

	POINT4D pos;  // ��������������е�λ��
	VECTOR4D dir; // �����ע�ӷ���

	VECTOR4D u;
	VECTOR4D v;
	VECTOR4D n;
	VECTOR4D target; // Ŀ��λ��

	float view_dist; // �Ӿ�

	float fov; // ˮƽ����ʹ�ֱ�������Ұ

	float near_clip_z; // ���������z
	float far_clip_z;  // Զ�������z

	PLANE3D rt_clip_plane; // �Ҽ�����
	PLANE3D lt_clip_plane; // �������
	PLANE3D tp_clip_plane; // �ϼ�����
	PLANE3D bt_clip_plane; // �¼�����

	float viewplane_width;  // ��ƽ��Ŀ�
	float viewplane_height; // ��ƽ��ĸ�

	float viewport_width;	// ��Ļ�Ŀ�
	float viewport_height;   // ��Ļ�ĸ�
	float viewport_center_x; // ��Ļ���ĵ�x
	float viewport_center_y; // ��Ļ���ĵ�y
	float aspect_ratio;		 // ��Ļ��߱�

	MATRIX4X4 mcam; // �������굽�������ı任����
	MATRIX4X4 mper; // ������굽͸������ı任����
	MATRIX4X4 mscr; // ͸�����굽��Ļ����ı任����

} CAM4DV1, *CAM4DV1_PTR;

// void Reset_RENDERLIST4DV1(RENDERLIST4DV1_PTR renderList);
void Init_CAM4DV1(CAM4DV1_PTR cam, int attr, POINT4D_PTR cam_pos, VECTOR4D_PTR cam_dir, VECTOR4D_PTR cam_target, float near_clip_z, float far_clip_z, float fov, float viewport_width, float viewport_height);

int Insert_POLYF4DV1_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, POLYF4DV1_PTR poly);
void Build_CAM4DV1_Matrix_Euler(CAM4DV1_PTR cam, int cam_rot_seq);

//����plg��ʽ�ļ�
int Load_OBJECT4DV1_PLG(OBJECT4DV1_PTR obj, char *filename, VECTOR4D_PTR scale, VECTOR4D_PTR pos, VECTOR4D_PTR rot);
char *Get_Line_PLG(char *buffer, int maxlength, FILE *fp);
float Compute_OBJECT4DV1_Radius(OBJECT4DV1_PTR obj);

#define SET_BIT(word, bit_flag) ((word) = ((word) | (bit_flag)))

void Build_Sin_Cos_Tables(void);