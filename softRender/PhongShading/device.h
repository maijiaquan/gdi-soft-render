
#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <windows.h>
#include <tchar.h>

#include <iostream>

#include <windows.h>
#include "conio.h"

#define _RGB565FROM16BIT(RGB, r,g,b) { *r = ( ((RGB) >> 11) & 0x1f); *g = (((RGB) >> 5) & 0x3f); *b = ((RGB) & 0x1f); }

#define RENDER_STATE_WIREFRAME 1 // ????
#define RENDER_STATE_TEXTURE 2   // ????
#define RENDER_STATE_COLOR 4	 // ????

typedef unsigned int IUINT32;
//=====================================================================
// ????
//=====================================================================
typedef struct
{
	int width;			   // ????
	int height;			   // ????
	IUINT32 **framebuffer; // ?????framebuffer[y] ??? y?
	float **zbuffer;	   // ?????zbuffer[y] ?? y???
	IUINT32 **texture;	 // ??????????
	int tex_width;		   // ????
	int tex_height;		   // ????
	float max_u;		   // ???????tex_width - 1
	float max_v;		   // ???????tex_height - 1
	int render_state;	  // ????
	IUINT32 background;	// ????
	IUINT32 foreground;	// ????
} device_t;

void device_pixel(device_t *device, int x, int y, IUINT32 color);
void device_draw_line(device_t *device, int x1, int y1, int x2, int y2, IUINT32 c);

void device_init(device_t *device, int width, int height, void *fb);
void device_destroy(device_t *device);
void device_clear(device_t *device, int mode);

#define WINDOW_WIDTH 600 // size of window
#define WINDOW_HEIGHT 600

#define min_clip_x  0 // clipping rectangle
#define max_clip_x  (WINDOW_WIDTH - 1)
#define min_clip_y  0
#define max_clip_y  (WINDOW_HEIGHT - 1)

void DrawTopTriangle(device_t *device, int x1, int y1, int x2, int y2, int x3, int y3,  int color);
void DrawDownTriangle(device_t *device,int x1, int y1, int x2, int y2, int x3, int y3,  int color);
void DrawTrianglePureColor(device_t *device, int x1, int y1, int x2, int y2, int x3, int y3, int color);


void DrawTopTriangle2(device_t *device, float x1, float y1, float x2, float y2, float x3, float y3,  int color);
void DrawDownTriangle2(device_t *device,float x1, float y1, float x2, float y2, float x3, float y3,  int color);
void DrawTrianglePureColor2(device_t *device, float x1, float y1, float x2, float y2, float x3, float y3, int color);

// defines for small numbers
#define EPSILON_E3 (float)(1E-3)
#define EPSILON_E4 (float)(1E-4)
#define EPSILON_E5 (float)(1E-5)
#define EPSILON_E6 (float)(1E-6)
#define FCMP(a,b) ( (fabs(a-b) < EPSILON_E3) ? 1 : 0)
#define SWAP(a,b,t) {t=a; a=b; b=t;}

void RGBFrom565(int color, IUINT32 &c);