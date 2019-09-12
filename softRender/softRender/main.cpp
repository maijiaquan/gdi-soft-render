//=====================================================================
//
// mini3d.c - Mini Software Render All-In-One
//
// build:
//   mingw: gcc -O3 mini3d.c -o mini3d.exe -lgdi32
//   msvc:  cl -O2 -nologo mini3d.c
//
// history:
//   2007.7.01  skywind  create this file as a tutorial
//   2007.7.02  skywind  implementate texture and color render
//   2008.3.15  skywind  fixed a trapezoid issue
//   2015.8.09  skywind  rewrite with more comment
//   2015.8.12  skywind  adjust interfaces for clearity
//
//=====================================================================
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <windows.h>
#include <tchar.h>

#include <iostream>

#include "frame.h"
#include <windows.h>

typedef unsigned int IUINT32;

int CMID(int x, int min, int max) { return (x < min) ? min : ((x > max) ? max : x); }

//=====================================================================
// 渲染设备
//=====================================================================
typedef struct
{
	int width;			   // 窗口宽度
	int height;			   // 窗口高度
	IUINT32 **framebuffer; // 像素缓存：framebuffer[y] 代表第 y行
	float **zbuffer;	   // 深度缓存：zbuffer[y] 为第 y行指针
	IUINT32 **texture;	 // 纹理：同样是每行索引
	int tex_width;		   // 纹理宽度
	int tex_height;		   // 纹理高度
	float max_u;		   // 纹理最大宽度：tex_width - 1
	float max_v;		   // 纹理最大高度：tex_height - 1
	int render_state;	  // 渲染状态
	IUINT32 background;	// 背景颜色
	IUINT32 foreground;	// 线框颜色
} device_t;

#define RENDER_STATE_WIREFRAME 1 // 渲染线框
#define RENDER_STATE_TEXTURE 2   // 渲染纹理
#define RENDER_STATE_COLOR 4	 // 渲染颜色

// 设备初始化，fb为外部帧缓存，非 NULL 将引用外部帧缓存（每行 4字节对齐）
void device_init(device_t *device, int width, int height, void *fb)
{
	int need = sizeof(void *) * (height * 2 + 1024) + width * height * 8;
	char *ptr = (char *)malloc(need + 64);
	char *framebuf, *zbuf;
	int j;
	assert(ptr);
	device->framebuffer = (IUINT32 **)ptr;
	device->zbuffer = (float **)(ptr + sizeof(void *) * height);
	ptr += sizeof(void *) * height * 2;
	device->texture = (IUINT32 **)ptr;
	ptr += sizeof(void *) * 1024;
	framebuf = (char *)ptr;
	zbuf = (char *)ptr + width * height * 4;
	ptr += width * height * 8;
	if (fb != NULL)
		framebuf = (char *)fb;
	for (j = 0; j < height; j++)
	{
		device->framebuffer[j] = (IUINT32 *)(framebuf + width * 4 * j);
		device->zbuffer[j] = (float *)(zbuf + width * 4 * j);
	}
	device->texture[0] = (IUINT32 *)ptr;
	device->texture[1] = (IUINT32 *)(ptr + 16);
	memset(device->texture[0], 0, 64);
	device->tex_width = 2;
	device->tex_height = 2;
	device->max_u = 1.0f;
	device->max_v = 1.0f;
	device->width = width;
	device->height = height;
	device->background = 0xc0c0c0;
	device->foreground = 0;
	device->render_state = RENDER_STATE_WIREFRAME;
}

// 删除设备
void device_destroy(device_t *device)
{
	if (device->framebuffer)
		free(device->framebuffer);
	device->framebuffer = NULL;
	device->zbuffer = NULL;
	device->texture = NULL;
}

// 清空 framebuffer 和 zbuffer
void device_clear(device_t *device, int mode)
{
	int y, x, height = device->height;
	for (y = 0; y < device->height; y++)
	{
		IUINT32 *dst = device->framebuffer[y];
		// IUINT32 cc = (height - 1 - y) * 230 / (height - 1);
		IUINT32 cc = (height - 1 - 0) * 230 / (height - 1);
		cc = (cc << 16) | (cc << 8) | cc;
		// std::cout<<"cc = "<<cc<<std::endl;
		if (mode == 0)
			cc = device->background;
		for (x = device->width; x > 0; dst++, x--)
			dst[0] = cc;
	}
	for (y = 0; y < device->height; y++)
	{
		float *dst = device->zbuffer[y];
		for (x = device->width; x > 0; dst++, x--)
			dst[0] = 0.0f;
	}
}

// 画点
void device_pixel(device_t *device, int x, int y, IUINT32 color)
{
	if (((IUINT32)x) < (IUINT32)device->width && ((IUINT32)y) < (IUINT32)device->height)
	{
		device->framebuffer[y][x] = color;
	}
}

// 绘制线段
void device_draw_line(device_t *device, int x1, int y1, int x2, int y2, IUINT32 c)
{
	// std::cout<<"x1 = "<<x1<<"y1 = "<<y1<<std::endl;
	// std::cout<<"x2 = "<<x2<<"y2 = "<<y2<<std::endl;

	int x, y, rem = 0;
	if (x1 == x2 && y1 == y2)
	{
		device_pixel(device, x1, y1, c);
	}
	else if (x1 == x2)
	{
		int inc = (y1 <= y2) ? 1 : -1;
		for (y = y1; y != y2; y += inc)
			device_pixel(device, x1, y, c);
		device_pixel(device, x2, y2, c);
	}
	else if (y1 == y2)
	{
		int inc = (x1 <= x2) ? 1 : -1;
		for (x = x1; x != x2; x += inc)
			device_pixel(device, x, y1, c);
		device_pixel(device, x2, y2, c);
	}
	else
	{
		int dx = (x1 < x2) ? x2 - x1 : x1 - x2;
		int dy = (y1 < y2) ? y2 - y1 : y1 - y2;
		if (dx >= dy)
		{
			if (x2 < x1)
				x = x1, y = y1, x1 = x2, y1 = y2, x2 = x, y2 = y;
			for (x = x1, y = y1; x <= x2; x++)
			{
				device_pixel(device, x, y, c);
				rem += dy;
				if (rem >= dx)
				{
					rem -= dx;
					y += (y2 >= y1) ? 1 : -1;
					device_pixel(device, x, y, c);
				}
			}
			device_pixel(device, x2, y2, c);
		}
		else
		{
			if (y2 < y1)
				x = x1, y = y1, x1 = x2, y1 = y2, x2 = x, y2 = y;
			for (x = x1, y = y1; y <= y2; y++)
			{
				device_pixel(device, x, y, c);
				rem += dx;
				if (rem >= dy)
				{
					rem -= dy;
					x += (x2 >= x1) ? 1 : -1;
					device_pixel(device, x, y, c);
				}
			}
			device_pixel(device, x2, y2, c);
		}
	}
}

//=====================================================================
// Win32 窗口及图形绘制：为 device 提供一个 DibSection 的 FB
//=====================================================================
int screen_w, screen_h, screen_exit = 0;
int screen_mx = 0, screen_my = 0, screen_mb = 0;
int screen_keys[512];			  // 当前键盘按下状态
static HWND screen_handle = NULL; // 主窗口 HWND
static HDC screen_dc = NULL;	  // 配套的 HDC
static HBITMAP screen_hb = NULL;  // DIB
static HBITMAP screen_ob = NULL;  // 老的 BITMAP
unsigned char *screen_fb = NULL;  // frame buffer
long screen_pitch = 0;

int screen_init(int w, int h, const TCHAR *title); // 屏幕初始化
int screen_close(void);							   // 关闭屏幕
void screen_dispatch(void);						   // 处理消息
void screen_update(void);						   // 显示 FrameBuffer

// win32 event handler
static LRESULT screen_events(HWND, UINT, WPARAM, LPARAM);

#ifdef _MSC_VER
#pragma comment(lib, "gdi32.lib")
#pragma comment(lib, "user32.lib")
#endif

// 初始化窗口并设置标题
int screen_init(int w, int h, const TCHAR *title)
{
	WNDCLASS wc = {CS_BYTEALIGNCLIENT, (WNDPROC)screen_events, 0, 0, 0,
				   NULL, NULL, NULL, NULL, _T("SCREEN3.1415926")};
	BITMAPINFO bi = {{sizeof(BITMAPINFOHEADER), w, -h, 1, 32, BI_RGB,
					  w * h * 4, 0, 0, 0, 0}};
	RECT rect = {0, 0, w, h};
	int wx, wy, sx, sy;
	LPVOID ptr;
	HDC hDC;

	screen_close();

	wc.hbrBackground = (HBRUSH)GetStockObject(BLACK_BRUSH);
	wc.hInstance = GetModuleHandle(NULL);
	wc.hCursor = LoadCursor(NULL, IDC_ARROW);
	if (!RegisterClass(&wc))
		return -1;

	screen_handle = CreateWindow(_T("SCREEN3.1415926"), title,
								 WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU | WS_MINIMIZEBOX,
								 0, 0, 0, 0, NULL, NULL, wc.hInstance, NULL);
	if (screen_handle == NULL)
		return -2;

	screen_exit = 0;
	hDC = GetDC(screen_handle);
	screen_dc = CreateCompatibleDC(hDC);
	ReleaseDC(screen_handle, hDC);

	screen_hb = CreateDIBSection(screen_dc, &bi, DIB_RGB_COLORS, &ptr, 0, 0);
	if (screen_hb == NULL)
		return -3;

	screen_ob = (HBITMAP)SelectObject(screen_dc, screen_hb);
	screen_fb = (unsigned char *)ptr;
	screen_w = w;
	screen_h = h;
	screen_pitch = w * 4;

	AdjustWindowRect(&rect, GetWindowLong(screen_handle, GWL_STYLE), 0);
	wx = rect.right - rect.left;
	wy = rect.bottom - rect.top;
	sx = (GetSystemMetrics(SM_CXSCREEN) - wx) / 2;
	sy = (GetSystemMetrics(SM_CYSCREEN) - wy) / 2;
	if (sy < 0)
		sy = 0;
	SetWindowPos(screen_handle, NULL, sx, sy, wx, wy, (SWP_NOCOPYBITS | SWP_NOZORDER | SWP_SHOWWINDOW));
	SetForegroundWindow(screen_handle);

	ShowWindow(screen_handle, SW_NORMAL);
	screen_dispatch();

	memset(screen_keys, 0, sizeof(int) * 512);
	memset(screen_fb, 0, w * h * 4);

	return 0;
}

int screen_close(void)
{
	if (screen_dc)
	{
		if (screen_ob)
		{
			SelectObject(screen_dc, screen_ob);
			screen_ob = NULL;
		}
		DeleteDC(screen_dc);
		screen_dc = NULL;
	}
	if (screen_hb)
	{
		DeleteObject(screen_hb);
		screen_hb = NULL;
	}
	if (screen_handle)
	{
		CloseWindow(screen_handle);
		screen_handle = NULL;
	}
	return 0;
}

static LRESULT screen_events(HWND hWnd, UINT msg,
							 WPARAM wParam, LPARAM lParam)
{
	switch (msg)
	{
	case WM_CLOSE:
		screen_exit = 1;
		break;
	case WM_KEYDOWN:
		screen_keys[wParam & 511] = 1;
		break;
	case WM_KEYUP:
		screen_keys[wParam & 511] = 0;
		break;
	default:
		return DefWindowProc(hWnd, msg, wParam, lParam);
	}
	return 0;
}

void screen_dispatch(void)
{
	MSG msg;
	while (1)
	{
		if (!PeekMessage(&msg, NULL, 0, 0, PM_NOREMOVE))
			break;
		if (!GetMessage(&msg, NULL, 0, 0))
			break;
		DispatchMessage(&msg);
	}
}

void screen_update(void)
{
	HDC hDC = GetDC(screen_handle);
	BitBlt(hDC, 0, 0, screen_w, screen_h, screen_dc, 0, 0, SRCCOPY);
	ReleaseDC(screen_handle, hDC);
	screen_dispatch();
}

//=====================================================================
// 主程序
//=====================================================================

//--------------------------------------------------------------
// initialize camera position and direction
POINT4D cam_pos = {0, 0, -100, 1};
VECTOR4D cam_dir = {0, 0, 0, 1};

// all your initialization code goes here...
VECTOR4D vscale = {.5, .5, .5, 1},
		 vpos = {0, 0, 0, 1},
		 vrot = {0, 0, 0, 1};

RENDERLIST4DV1 rend_list;			// the single renderlist
POLYF4DV1 poly1;					// our lonely polygon
CAM4DV1 cam;						// the single camera
POINT4D poly1_pos = {0, 0, 100, 1}; // world position of polygon

device_t device;

void GameInit();
void GameMain();

USHORT(*RGB16Bit)
(int r, int g, int b) = nullptr;

void Build_Sin_Cos_Tables(void);

void Build_Sin_Cos_Tables(void)
{
	for (int ang = 0; ang <= 360; ang++)
	{
		float theta = (float)ang * PI / (float)180;
		cos_look[ang] = cos(theta);
		sin_look[ang] = sin(theta);
	}
}

void GameInit()
{
	RGB16Bit = RGB16Bit565;

	Build_Sin_Cos_Tables();
	poly1.state = POLY4DV1_STATE_ACTIVE;
	poly1.attr = 0;
	poly1.color = RGB16Bit(0, 255, 0);

	poly1.vlist[0].x = 0;
	poly1.vlist[0].y = 50;
	poly1.vlist[0].z = 0;
	poly1.vlist[0].w = 1;

	poly1.vlist[1].x = 50;
	poly1.vlist[1].y = -50;
	poly1.vlist[1].z = 0;
	poly1.vlist[1].w = 1;

	poly1.vlist[2].x = -50;
	poly1.vlist[2].y = -50;
	poly1.vlist[2].z = 0;
	poly1.vlist[2].w = 1;

	poly1.next = poly1.prev = NULL;

	// initialize the camera with 90 FOV, normalized coordinates
	Init_CAM4DV1(&cam,			  // the camera object
				 CAM_MODEL_EULER, // euler camera model
				 &cam_pos,		  // initial camera position
				 &cam_dir,		  // initial camera angles
				 NULL,			  // no initial target
				 50.0,			  // near and far clipping planes
				 500.0,
				 90.0,		   // field of view in degrees
				 WINDOW_WIDTH, // size of final screen viewport
				 WINDOW_HEIGHT);
}

void GameMain()
{
	static MATRIX4X4 gMatrixRotate; // general rotation matrix

	static float ang_y = 0; // rotation angle

	Reset_RENDERLIST4DV1(&rend_list);
	Insert_POLYF4DV1_RENDERLIST4DV1(&rend_list, &poly1);
	Build_XYZ_Rotation_MATRIX4X4(0, ang_y, 0, &gMatrixRotate);

	ang_y += 10;
	if (ang_y >= 360.0)
		ang_y = 0;

	Transform_RENDERLIST4DV1(&rend_list, &gMatrixRotate, TRANSFORM_LOCAL_ONLY);
	Model_To_World_RENDERLIST4DV1(&rend_list, &poly1_pos);
	Build_CAM4DV1_Matrix_Euler(&cam, CAM_ROT_SEQ_ZYX);
	World_To_Camera_RENDERLIST4DV1(&rend_list, &cam);
	Camera_To_Perspective_RENDERLIST4DV1(&rend_list, &cam);
	Perspective_To_Screen_RENDERLIST4DV1(&rend_list, &cam);

	RENDERLIST4DV1_PTR rend_list_ptr = &rend_list;

	for (int idx_poly = 0; idx_poly < rend_list_ptr->num_polys; idx_poly++)
	{
		// std::cout << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[0].x << std::endl;
		// std::cout << "x1 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[0].x << std::endl;
		// std::cout << "y1 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[0].y << std::endl;
		// std::cout << "x2 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[1].x << std::endl;
		// std::cout << "y2 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[1].y << std::endl;
		// std::cout << "x3 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[2].x << std::endl;
		// std::cout << "y3 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[2].y << std::endl;

		float x1 = rend_list_ptr->poly_ptrs[idx_poly]->tvlist[0].x;
		float y1 = rend_list_ptr->poly_ptrs[idx_poly]->tvlist[0].y;
		float x2 = rend_list_ptr->poly_ptrs[idx_poly]->tvlist[1].x;
		float y2 = rend_list_ptr->poly_ptrs[idx_poly]->tvlist[1].y;
		float x3 = rend_list_ptr->poly_ptrs[idx_poly]->tvlist[2].x;
		float y3 = rend_list_ptr->poly_ptrs[idx_poly]->tvlist[2].y;

		device_draw_line(&device, x1, y1, x2, y2, device.foreground); //3 1
		device_draw_line(&device, x1, y1, x3, y3, device.foreground); //3 1
		device_draw_line(&device, x2, y2, x3, y3, device.foreground); //3 1

	} // end for poly
}

int main(void)
{
	DWORD time_start, time_end;

	bool isOnlyBox = false;
	int states[] = {RENDER_STATE_TEXTURE, RENDER_STATE_COLOR, RENDER_STATE_WIREFRAME};
	int indicator = 0;
	int kbhit = 0;
	float alpha = 1;
	float pos = 3.5;

	TCHAR *title = _T("Mini3d (software render tutorial) - ")
				   _T("Left/Right: rotation, Up/Down: forward/backward, Space: switch state");

	if (screen_init(800, 600, title))
		return -1;

	device_init(&device, 800, 600, screen_fb);

	//init_texture(&device);
	device.render_state = RENDER_STATE_TEXTURE;

	int tmp = 1;

	GameInit();

	int deltaTime = 0;

	while (screen_exit == 0 && screen_keys[VK_ESCAPE] == 0)
	{
		/* 获取开始时间 */
		time_start = GetTickCount(); //从操作系统启动经过的毫秒数

		screen_dispatch();
		device_clear(&device, 1);

		GameMain();

		screen_update();
		// Sleep(50);

		time_end = GetTickCount();
		deltaTime = (time_end - time_start);
		// std::cout<<"delta time = "<<deltaTime<<std::endl;
		std::cout << "delta time = " << 1000 / float(deltaTime) << std::endl;
	}
	return 0;
}
