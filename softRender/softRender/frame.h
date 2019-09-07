#pragma once
#define RENDERLIST4DV1_MAX_POLYS 32768 // 16384

typedef struct Vector4Type
{
	union {
		float M[4];

		struct
		{
			float x, y, z, w;
		};
	};	
} Vector4, Point4, *Vector4Ptr, *Point4Ptr;

//Polygon多边形，其实就是一个三角形
typedef struct Poly4Type
{
    int state; // state information
    int attr;  // physical attributes of polygon
    int color; // color of polygon

    Point4 vlist[3];  // the vertices of this triangle
    Point4 tvlist[3]; // the vertices after transformation if needed

    Poly4Type *next; // pointer to next polygon in list??
    Poly4Type *prev; // pointer to previous polygon in list??

} Poly4, *Poly4Ptr;

//渲染列表，其实就是三角形数组
typedef struct RenderList4Type
{
    int state; // state of renderlist ???
    int attr;  // attributes of renderlist ???

    Poly4Ptr polyPtrs[RENDERLIST4DV1_MAX_POLYS];

    Poly4 polyData[RENDERLIST4DV1_MAX_POLYS];

    int numPolys; // number of polys in render list

} RenderList4, *RenderList4Ptr;

void ResetRenderList(RenderList4Ptr renderList);
