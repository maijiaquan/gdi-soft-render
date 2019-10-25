

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <windows.h>
#include <tchar.h>

#include <iostream>

// #include "frame.h"
#include "transform.h"
#include "datastructure.h"

#include <windows.h>
#include "conio.h"

HANDLE hStdout;
//   光标位置
COORD cursorPos;

static DWORD time_start, time_end;

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
	device->foreground = 0xc0c0c0;
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
		// IUINT32 cc = (height - 1 - 0) * 255 / (height - 1);
		IUINT32 cc = (height - 1 - 0) * 0 / (height - 1); //黑色背景
		int R = 0;
		int G = 0;
		int B = 0;
		//cc = (cc << 16) | (cc << 8) | cc;
		cc = (R << 16) | (G << 8) | B;

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

static float gRotationAngle = 0; // rotation angle

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

void DrawText(char *text)
{
	HDC hDC = GetDC(screen_handle);

	SetTextColor(hDC, RGB(0, 255, 0));
	SetBkMode(hDC, TRANSPARENT);
	// char* text = "abc爱抚";
	// char* text = "";
	// _itoa(count++, text, 5);

	TextOutA(hDC, 10, 10, text, strlen(text));

	// ReleaseDC(screen_handle, hDC);
	// ReleaseDC(screen_handle, hDC);
}

void DrawTextOnScreen(char *text, int x, int y)
{
	HDC hDC = GetDC(screen_handle);

	SetTextColor(hDC, RGB(0, 255, 0));
	SetBkMode(hDC, TRANSPARENT);
	// char* text = "abc爱抚";
	// char* text = "";
	// _itoa(count++, text, 5);

	TextOutA(hDC, x, y, text, strlen(text));
}
void screen_update(void)
{
	HDC hDC = GetDC(screen_handle);
	BitBlt(hDC, 0, 0, screen_w, screen_h, screen_dc, 0, 0, SRCCOPY);

	ReleaseDC(screen_handle, hDC);

	// char text[100] = "Rotation Angle: ";
	// int tmp = gRotationAngle;
	// char textInt[10];
	// _itoa(tmp, textInt, 10);
	// strcat(text, textInt);
	// DrawText(text);
	screen_dispatch();
}

//=====================================================================
// 主程序
//=====================================================================

//--------------------------------------------------------------
// initialize camera position and direction

RENDERLIST4DV1 rend_list;			// the single renderlist
POLYF4DV1 poly1;					// our lonely polygon
CAM4DV1 cam;						// the single camera
POINT4D poly1_pos = {0, 0, 100, 1}; // world position of polygon
OBJECT4DV1 obj;						// used to hold our cube mesh
// all your initialization code goes here...

device_t device;

void GameInit();
void GameMain();

// USHORT(*RGB16Bit)(int r, int g, int b);
USHORT(*RGB16Bit)
(int r, int g, int b);

void InitDemo7_1();
void DrawDemo7_1();
void InitDemo7_2();
void DrawDemo7_2();

void Initmo7_4();
void DrawDemo7_4();

void Initmo7_6();
void DrawDemo7_6();

void InitDemo7_1()
{
	POINT4D cam_pos = {0, 0, -100, 1};
	VECTOR4D cam_dir = {0, 0, 0, 1};

	VECTOR4D vscale = {5.0, 5.0, 5.0, 1}, // scale of object
		vpos = {0, 0, 0, 1},			  // position of object
		vrot = {0, 0, 0, 1};			  // initial orientation of object

	RGB16Bit = RGB16Bit565;

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
	Init_CAM4DV1(&cam, CAM_MODEL_EULER, &cam_pos, &cam_dir, NULL, 50.0, 500.0, 90.0, WINDOW_WIDTH, WINDOW_HEIGHT);
}

void InitDemo7_2()
{
	POINT4D cam_pos = {0, 0, -100, 1};
	VECTOR4D cam_dir = {0, 0, 0, 1};

	VECTOR4D vscale = {5.0, 5.0, 5.0, 1}, // scale of object
		vpos = {0, 0, 0, 1},			  // position of object
		vrot = {0, 0, 0, 1};			  // initial orientation of object

	RGB16Bit = RGB16Bit565;

	Init_CAM4DV1(&cam, CAM_MODEL_EULER, &cam_pos, &cam_dir, NULL, 50.0, 500.0, 90.0, WINDOW_WIDTH, WINDOW_HEIGHT);

	// Load_OBJECT4DV1_PLG(&obj, "./plg/tank1.plg", &vscale, &vpos, &vrot);
	Load_OBJECT4DV1_PLG(&obj, "./plg/cube2.plg", &vscale, &vpos, &vrot);
	// Load_OBJECT4DV1_PLG(&obj, "cube1.plg", &vscale, &vpos, &vrot);

	obj.world_pos.x = 0;
	obj.world_pos.y = 0;
	obj.world_pos.z = 100;
}

void DrawDemo7_1()
{
	char text[100] = "Rotation Angle: ";
	int tmp = gRotationAngle;
	char textInt[10];
	_itoa(tmp, textInt, 10);
	strcat(text, textInt);
	DrawText(text);

	Sleep(20);

	static MATRIX4X4 mrot; // general rotation matrix
	gRotationAngle += 1;
	gRotationAngle = gRotationAngle >= 360.0 ? 0 : gRotationAngle;

	Reset_RENDERLIST4DV1(&rend_list);
	Insert_POLYF4DV1_RENDERLIST4DV1(&rend_list, &poly1); //每一次都从重新赋值，相当于每一次都重置了该三角形的顶点坐标

	Build_XYZ_Rotation_MATRIX4X4(0, gRotationAngle, 0, &mrot); //构造旋转矩阵，绕y轴旋转

	Transform_RENDERLIST4DV1(&rend_list, &mrot, TRANSFORM_LOCAL_ONLY); //在局部坐标中实现旋转矩阵的变换

	Model_To_World_RENDERLIST4DV1(&rend_list, &poly1_pos, TRANSFORM_LOCAL_TO_TRANS); //平移到世界坐标

	Build_CAM4DV1_Matrix_Euler(&cam, CAM_ROT_SEQ_ZYX); //构造欧拉相机矩阵
	World_To_Camera_RENDERLIST4DV1(&rend_list, &cam);  //世界坐标到相机坐标的变换

	Camera_To_Perspective_RENDERLIST4DV1(&rend_list, &cam); //相机坐标到透视坐标

	Perspective_To_Screen_RENDERLIST4DV1(&rend_list, &cam); //透视坐标到屏幕坐标

	RENDERLIST4DV1_PTR rend_list_ptr = &rend_list;

	for (int poly = 0; poly < rend_list_ptr->num_polys; poly++)
	{
		float x1 = rend_list_ptr->poly_ptrs[poly]->tvlist[0].x;
		float y1 = rend_list_ptr->poly_ptrs[poly]->tvlist[0].y;
		float x2 = rend_list_ptr->poly_ptrs[poly]->tvlist[1].x;
		float y2 = rend_list_ptr->poly_ptrs[poly]->tvlist[1].y;
		float x3 = rend_list_ptr->poly_ptrs[poly]->tvlist[2].x;
		float y3 = rend_list_ptr->poly_ptrs[poly]->tvlist[2].y;

		device_draw_line(&device, x1, y1, x2, y2, device.foreground); //3 1
		device_draw_line(&device, x1, y1, x3, y3, device.foreground); //3 1
		device_draw_line(&device, x2, y2, x3, y3, device.foreground); //3 1
	}
}

void DrawDemo7_2()
{
	char text[100] = "Rotation Angle: ";
	int tmp = gRotationAngle;
	char textInt[10];
	_itoa(tmp, textInt, 10);
	strcat(text, textInt);
	DrawText(text);

	Sleep(10);

	static MATRIX4X4 mrot; // general rotation matrix

	gRotationAngle = 1;

	Reset_OBJECT4DV1(&obj);

	Build_XYZ_Rotation_MATRIX4X4(0, gRotationAngle, 0, &mrot); //构造旋转矩阵，绕y轴旋转

	Transform_OBJECT4DV1(&obj, &mrot, TRANSFORM_LOCAL_ONLY, 1);

	Model_To_World_OBJECT4DV1(&obj, TRANSFORM_LOCAL_TO_TRANS);

	Remove_Backfaces_OBJECT4DV1(&obj, &cam);

	Build_CAM4DV1_Matrix_Euler(&cam, CAM_ROT_SEQ_ZYX);
	World_To_Camera_OBJECT4DV1(&obj, &cam);
	Camera_To_Perspective_OBJECT4DV1(&obj, &cam);
	Perspective_To_Screen_OBJECT4DV1(&obj, &cam);

	OBJECT4DV1_PTR obj_ptr = &obj;
	for (int poly = 0; poly < obj_ptr->num_polys; poly++)
	{
		if (!(obj_ptr->plist[poly].state & POLY4DV1_STATE_ACTIVE) || (obj_ptr->plist[poly].state & POLY4DV1_STATE_CLIPPED) || (obj_ptr->plist[poly].state & POLY4DV1_STATE_BACKFACE))
			continue;

		int vindex_0 = obj_ptr->plist[poly].vert[0];
		int vindex_1 = obj_ptr->plist[poly].vert[1];
		int vindex_2 = obj_ptr->plist[poly].vert[2];

		device_draw_line(&device, obj_ptr->vlist_trans[vindex_0].x, obj_ptr->vlist_trans[vindex_0].y, obj_ptr->vlist_trans[vindex_1].x, obj_ptr->vlist_trans[vindex_1].y, device.foreground); //3 1
		device_draw_line(&device, obj_ptr->vlist_trans[vindex_1].x, obj_ptr->vlist_trans[vindex_1].y, obj_ptr->vlist_trans[vindex_2].x, obj_ptr->vlist_trans[vindex_2].y, device.foreground); //3 1
		device_draw_line(&device, obj_ptr->vlist_trans[vindex_2].x, obj_ptr->vlist_trans[vindex_2].y, obj_ptr->vlist_trans[vindex_0].x, obj_ptr->vlist_trans[vindex_0].y, device.foreground); //3 1
	}
}

char *work_string;
char *textBuffer; // used to print text

void InitDemo7_4()
{
	POINT4D cam_pos = {0, 200, 0, 1};
	VECTOR4D cam_dir = {0, 0, 0, 1};

	// VECTOR4D vscale = {1.0, 1.0, 1.0, 1}, // scale of object
	VECTOR4D vscale = {10.0, 10.0, 10.0, 1}, // scale of object
		vpos = {0, 0, 0, 1},			  // position of object
		vrot = {0, 0, 0, 1};			  // initial orientation of object

	RGB16Bit = RGB16Bit565;
	Init_CAM4DV1(&cam, CAM_MODEL_EULER, &cam_pos, &cam_dir, NULL, 50.0, 1000.0, 90.0, WINDOW_WIDTH, WINDOW_HEIGHT);

	// load the object
    //Load_OBJECT4DV1_PLG(&obj, "./plg/cube1.plg", &vscale, &vpos, &vrot);
	Load_OBJECT4DV1_PLG(&obj, "./plg/cube2.plg", &vscale, &vpos, &vrot);
	// Load_OBJECT4DV1_PLG(&obj, "./plg/tank1.plg", &vscale, &vpos, &vrot);

	// set the default position of the object in the world
	obj.world_pos.x = 0;
	obj.world_pos.y = 0;
	obj.world_pos.z = 400;
}

void DrawDemo7_4()
{
	Sleep(20);

#define KEY_DOWN(vk_code) ((GetAsyncKeyState(vk_code) & 0x8000) ? 1 : 0)

	const int NUM_OBJECTS = 2;		// number of objects on a row
	const int OBJECT_SPACING = 250; // spacing between objects

	work_string = new char[256];
	textBuffer = new char[1024]; // used to print text

	strcpy(textBuffer, "Objects Culled: ");

	static MATRIX4X4 mrot;
	gRotationAngle = 1;
	Reset_RENDERLIST4DV1(&rend_list);

	// is user trying to rotate camera
	if (KEY_DOWN(VK_DOWN))
		cam.dir.x += 1;
	else if (KEY_DOWN(VK_UP))
		cam.dir.x -= 1;

	// is user trying to rotate camera
	if (KEY_DOWN(VK_RIGHT))
		cam.dir.y -= 1;
	else if (KEY_DOWN(VK_LEFT))
		cam.dir.y += 1;

	
	Build_XYZ_Rotation_MATRIX4X4(0, gRotationAngle, 0, &mrot);
	Transform_OBJECT4DV1(&obj, &mrot, TRANSFORM_LOCAL_ONLY, 1);
	int idxLine = 4;	//控制台打印行号
	cursorPos.X = 0;
	cursorPos.Y = idxLine;

	SetConsoleCursorPosition(hStdout, cursorPos);
	std::cout << "cam.far z = " << cam.far_clip_z << ", cam.near z =" << cam.near_clip_z  << std::endl;
	idxLine++;

	for (int x = -NUM_OBJECTS / 2; x < NUM_OBJECTS / 2; x++)
	{
		for (int z = -NUM_OBJECTS / 2; z < NUM_OBJECTS / 2; z++)
		{
			POINT4D sphere_pos;
			Mat_Mul_VECTOR4D_4X4(&obj.world_pos, &cam.mcam, &sphere_pos); //将包围球转换到相机空间

			char textCharArray[1024];

			strcpy(textCharArray, "Sphere pos: ");
			char intCharArray[10];

			_itoa(sphere_pos.x, intCharArray, 10);
			strcat(textCharArray, intCharArray);

			_itoa(sphere_pos.y, intCharArray, 10);
			strcat(textCharArray, intCharArray);

			_itoa(sphere_pos.z, intCharArray, 10);
			strcat(textCharArray, intCharArray);

			//DrawTextOnScreen(textCharArray, (x+1)*10+20, (z+1)*10+20);

			cursorPos.X = 0;
			cursorPos.Y = 4 + idxLine++;
			SetConsoleCursorPosition(hStdout, cursorPos);

			float z_test = (0.5) * cam.viewplane_width * sphere_pos.z / cam.view_dist;

			std::cout << x<<", "<< z ;
			std::cout<<" pos.z > " << ((sphere_pos.z - obj.max_radius) > cam.far_clip_z);
			std::cout <<" pos.z < " << ((sphere_pos.z + obj.max_radius) < cam.near_clip_z);
			std::cout <<" pos.x > " << ((sphere_pos.x - obj.max_radius) > z_test);
			std::cout <<" pos.x < " << ((sphere_pos.x + obj.max_radius) < -z_test);
			std::cout <<" pos.y > " << ((sphere_pos.y - obj.max_radius) > z_test);
			std::cout <<" pos.y < " << ((sphere_pos.y + obj.max_radius) < -z_test);

			std::cout<<std::endl;
			
			idxLine++;
			idxLine++;


			Reset_OBJECT4DV1(&obj);

			obj.world_pos.x = x * OBJECT_SPACING + OBJECT_SPACING / 2;
			obj.world_pos.y = 0;
			obj.world_pos.z = 500 + z * OBJECT_SPACING + OBJECT_SPACING / 2;

			if (Cull_OBJECT4DV1(&obj, &cam, CULL_OBJECT_XYZ_PLANES))
			{
				sprintf(work_string, "[%d, %d] ", x, z);
				strcat(textBuffer, work_string);
			}
			else
			{
				Model_To_World_OBJECT4DV1(&obj);
				Insert_OBJECT4DV1_RENDERLIST4DV1(&rend_list, &obj); //将物体插入到渲染列表
			}
		}
	}

	SetConsoleCursorPosition(hStdout, cursorPos);
	std::cout<<rend_list.num_polys<<std::endl;

	Build_CAM4DV1_Matrix_Euler(&cam, CAM_ROT_SEQ_ZYX);



	// remove backfaces
	Remove_Backfaces_RENDERLIST4DV1(&rend_list, &cam);

	// apply world to camera transform
	World_To_Camera_RENDERLIST4DV1(&rend_list, &cam);

	// apply camera to perspective transformation
	Camera_To_Perspective_RENDERLIST4DV1(&rend_list, &cam);

	// apply screen transform
	Perspective_To_Screen_RENDERLIST4DV1(&rend_list, &cam);

	RENDERLIST4DV1_PTR rend_list_ptr = &rend_list;

	for (int poly = 0; poly < rend_list_ptr->num_polys; poly++)
	{

		if (!(rend_list_ptr->poly_ptrs[poly]->state & POLY4DV1_STATE_ACTIVE) ||(rend_list_ptr->poly_ptrs[poly]->state & POLY4DV1_STATE_CLIPPED) ||(rend_list_ptr->poly_ptrs[poly]->state & POLY4DV1_STATE_BACKFACE))
			continue; 

		float x1 = rend_list_ptr->poly_ptrs[poly]->tvlist[0].x;
		float y1 = rend_list_ptr->poly_ptrs[poly]->tvlist[0].y;
		float x2 = rend_list_ptr->poly_ptrs[poly]->tvlist[1].x;
		float y2 = rend_list_ptr->poly_ptrs[poly]->tvlist[1].y;
		float x3 = rend_list_ptr->poly_ptrs[poly]->tvlist[2].x;
		float y3 = rend_list_ptr->poly_ptrs[poly]->tvlist[2].y;

		IUINT32 c = (255 << 16) | (255 << 8) | 255;

		device_draw_line(&device, x1, y1, x2, y2, c);
		device_draw_line(&device, x1, y1, x3, y3, c);
		device_draw_line(&device, x2, y2, x3, y3, c);
	}

	// char text[100] = "Rotation Angle: ";
	// int tmp = gRotationAngle;
	// char textInt[10];
	// _itoa(tmp, textInt, 10);
	// strcat(text, textInt);
}


// defines for objects
#define NUM_TOWERS        96
#define NUM_TANKS         24
#define TANK_SPEED        15

#define UNIVERSE_RADIUS   4000


#define POINT_SIZE        200
#define NUM_POINTS_X      (2*UNIVERSE_RADIUS/POINT_SIZE)
#define NUM_POINTS_Z      (2*UNIVERSE_RADIUS/POINT_SIZE)
#define NUM_POINTS        (NUM_POINTS_X*NUM_POINTS_Z)

OBJECT4DV1 obj_tower, // used to hold the master tower
	obj_tank,		  // used to hold the master tank
	obj_marker,		  // the ground marker
	obj_player;		  // the player object

POINT4D towers[NUM_TOWERS],
	tanks[NUM_TANKS];

void InitDemo7_6()
{

	// all your initialization code goes here...
	VECTOR4D vscale = {1.0, 1.0, 1.0, 1},
			 vpos = {0, 0, 0, 1},
			 vrot = {0, 0, 0, 1};

	// initialize camera position and direction
	POINT4D cam_pos = {0, 40, 0, 1};
	POINT4D cam_target = {0, 0, 0, 1};
	VECTOR4D cam_dir = {0, 0, 0, 1};

	int index; // looping var

	srand(13);
	Init_CAM4DV1(&cam, CAM_MODEL_EULER, &cam_pos, &cam_dir, &cam_target, 200.0, 12000.0, 120.0, WINDOW_WIDTH, WINDOW_HEIGHT);

	// load the master tank object
	VECTOR4D_INITXYZ(&vscale, 0.75, 0.75, 0.75);
	Load_OBJECT4DV1_PLG(&obj_tank, "./plg/tank2.plg", &vscale, &vpos, &vrot);

	// load player object for 3rd person view
	VECTOR4D_INITXYZ(&vscale, 0.75, 0.75, 0.75);
	Load_OBJECT4DV1_PLG(&obj_player, "./plg/tank3.plg", &vscale, &vpos, &vrot);

	// load the master tower object
	VECTOR4D_INITXYZ(&vscale, 1.0, 2.0, 1.0);
	Load_OBJECT4DV1_PLG(&obj_tower, "./plg/tower1.plg", &vscale, &vpos, &vrot);

	// load the master ground marker
	VECTOR4D_INITXYZ(&vscale, 3.0, 3.0, 3.0);
	Load_OBJECT4DV1_PLG(&obj_marker, "./plg/marker1.plg", &vscale, &vpos, &vrot);

	// position the tanks
	for (index = 0; index < NUM_TANKS; index++)
	{
		// randomly position the tanks
		tanks[index].x = RAND_RANGE(-UNIVERSE_RADIUS, UNIVERSE_RADIUS);
		tanks[index].y = 0; // obj_tank.max_radius;
		tanks[index].z = RAND_RANGE(-UNIVERSE_RADIUS, UNIVERSE_RADIUS);
		tanks[index].w = RAND_RANGE(0, 360);
	} // end for

	// position the towers
	for (index = 0; index < NUM_TOWERS; index++)
	{
		// randomly position the tower
		towers[index].x = RAND_RANGE(-UNIVERSE_RADIUS, UNIVERSE_RADIUS);
		towers[index].y = 0; // obj_tower.max_radius;
		towers[index].z = RAND_RANGE(-UNIVERSE_RADIUS, UNIVERSE_RADIUS);
	} // end for
}

void DrawDemo7_6()
{
	Sleep(20);
	static MATRIX4X4 mrot; // general rotation matrix

	static float view_angle = 0;
	static float camera_distance = 6000;
	static VECTOR4D pos = {0, 0, 0, 0};
	static float tank_speed;
	static float turning = 0;

	char work_string[256]; // temp string

	int index; // looping var

	//Draw_Rectangle(0, 0, WINDOW_WIDTH - 1, WINDOW_HEIGHT / 2, RGB16Bit(0, 140, 192), lpddsback);

	// draw the ground
	//Draw_Rectangle(0, WINDOW_HEIGHT / 2, WINDOW_WIDTH - 1, WINDOW_HEIGHT - 1, RGB16Bit(103, 62, 3), lpddsback);

	Reset_RENDERLIST4DV1(&rend_list);

	

	// turbo
	if (KEY_DOWN(VK_SPACE))
		tank_speed = 5 * TANK_SPEED;
	else
		tank_speed = TANK_SPEED;

	// forward/backward
	if (KEY_DOWN(VK_UP))
	{
		// move forward
		cam.pos.x += tank_speed * Fast_Sin(cam.dir.y);

		cam.pos.z += tank_speed * Fast_Cos(cam.dir.y);
	} // end if

	if (KEY_DOWN(VK_DOWN))
	{
		// move backward
		cam.pos.x -= tank_speed * Fast_Sin(cam.dir.y);
		cam.pos.z -= tank_speed * Fast_Cos(cam.dir.y);
	} // end if

	// rotate
	if (KEY_DOWN(VK_RIGHT))
	{
		cam.dir.y += 3;

		// add a little turn to object
		if ((turning += 2) > 15)
			turning = 15;

	} // end if

	if (KEY_DOWN(VK_LEFT))
	{
		cam.dir.y -= 3;

		// add a little turn to object
		if ((turning -= 2) < -15)
			turning = -15;

	}	// end if
	else // center heading again
	{
		if (turning > 0)
			turning -= 1;
		else if (turning < 0)
			turning += 1;

	} // end else

	// generate camera matrix
	Build_CAM4DV1_Matrix_Euler(&cam, CAM_ROT_SEQ_ZYX);

	// insert the tanks in the world
	for (index = 0; index < NUM_TANKS; index++)
	{
		// reset the object (this only matters for backface and object removal)
		Reset_OBJECT4DV1(&obj_tank);

		// generate rotation matrix around y axis
		Build_XYZ_Rotation_MATRIX4X4(0, tanks[index].w, 0, &mrot);

		// rotate the local coords of the object
		Transform_OBJECT4DV1(&obj_tank, &mrot, TRANSFORM_LOCAL_TO_TRANS, 1);

		// set position of tank
		obj_tank.world_pos.x = tanks[index].x;
		obj_tank.world_pos.y = tanks[index].y;
		obj_tank.world_pos.z = tanks[index].z;

		// attempt to cull object
		if (!Cull_OBJECT4DV1(&obj_tank, &cam, CULL_OBJECT_XYZ_PLANES))
		{
			// if we get here then the object is visible at this world position
			// so we can insert it into the rendering list
			// perform local/model to world transform
			Model_To_World_OBJECT4DV1(&obj_tank, TRANSFORM_TRANS_ONLY);

			// insert the object into render list
			Insert_OBJECT4DV1_RENDERLIST4DV1(&rend_list, &obj_tank);
		} // end if

	} // end for

	// insert the player into the world
	// reset the object (this only matters for backface and object removal)
	Reset_OBJECT4DV1(&obj_player);

	// set position of tank
	obj_player.world_pos.x = cam.pos.x + 300 * Fast_Sin(cam.dir.y);
	obj_player.world_pos.y = cam.pos.y - 70;
	obj_player.world_pos.z = cam.pos.z + 300 * Fast_Cos(cam.dir.y);

	// generate rotation matrix around y axis
	Build_XYZ_Rotation_MATRIX4X4(0, cam.dir.y + turning, 0, &mrot);

	// rotate the local coords of the object
	Transform_OBJECT4DV1(&obj_player, &mrot, TRANSFORM_LOCAL_TO_TRANS, 1);

	// perform world transform
	Model_To_World_OBJECT4DV1(&obj_player, TRANSFORM_TRANS_ONLY);

	// insert the object into render list
	Insert_OBJECT4DV1_RENDERLIST4DV1(&rend_list, &obj_player);

	// insert the towers in the world
	for (index = 0; index < NUM_TOWERS; index++)
	{
		// reset the object (this only matters for backface and object removal)
		Reset_OBJECT4DV1(&obj_tower);

		// set position of tower
		obj_tower.world_pos.x = towers[index].x;
		obj_tower.world_pos.y = towers[index].y;
		obj_tower.world_pos.z = towers[index].z;

		// attempt to cull object
		if (!Cull_OBJECT4DV1(&obj_tower, &cam, CULL_OBJECT_XYZ_PLANES))
		{
			// if we get here then the object is visible at this world position
			// so we can insert it into the rendering list
			// perform local/model to world transform
			Model_To_World_OBJECT4DV1(&obj_tower);

			// insert the object into render list
			Insert_OBJECT4DV1_RENDERLIST4DV1(&rend_list, &obj_tower);
		} // end if

	} // end for

	// seed number generator so that modulation of markers is always the same
	srand(13);

	// insert the ground markers into the world
	for (int index_x = 0; index_x < NUM_POINTS_X; index_x++)
		for (int index_z = 0; index_z < NUM_POINTS_Z; index_z++)
		{
			// reset the object (this only matters for backface and object removal)
			Reset_OBJECT4DV1(&obj_marker);

			// set position of tower
			obj_marker.world_pos.x = RAND_RANGE(-100, 100) - UNIVERSE_RADIUS + index_x * POINT_SIZE;
			obj_marker.world_pos.y = obj_marker.max_radius;
			obj_marker.world_pos.z = RAND_RANGE(-100, 100) - UNIVERSE_RADIUS + index_z * POINT_SIZE;

			// attempt to cull object
			if (!Cull_OBJECT4DV1(&obj_marker, &cam, CULL_OBJECT_XYZ_PLANES))
			{
				// if we get here then the object is visible at this world position
				// so we can insert it into the rendering list
				// perform local/model to world transform
				Model_To_World_OBJECT4DV1(&obj_marker);

				// insert the object into render list
				Insert_OBJECT4DV1_RENDERLIST4DV1(&rend_list, &obj_marker);
			} // end if

		} // end for

	// remove backfaces
	Remove_Backfaces_RENDERLIST4DV1(&rend_list, &cam);

	// apply world to camera transform
	World_To_Camera_RENDERLIST4DV1(&rend_list, &cam);

	// apply camera to perspective transformation
	Camera_To_Perspective_RENDERLIST4DV1(&rend_list, &cam);

	// apply screen transform
	Perspective_To_Screen_RENDERLIST4DV1(&rend_list, &cam);


	// render the object
	//Draw_RENDERLIST4DV1_Wire16(&rend_list, back_buffer, back_lpitch);


	RENDERLIST4DV1_PTR rend_list_ptr = &rend_list;
	std::cout<<rend_list_ptr->num_polys<<std::endl;

	for (int poly = 0; poly < rend_list_ptr->num_polys; poly++)
	{
		if (!(rend_list_ptr->poly_ptrs[poly]->state & POLY4DV1_STATE_ACTIVE) ||(rend_list_ptr->poly_ptrs[poly]->state & POLY4DV1_STATE_CLIPPED) ||(rend_list_ptr->poly_ptrs[poly]->state & POLY4DV1_STATE_BACKFACE))
			continue; 
		float x1 = rend_list_ptr->poly_ptrs[poly]->tvlist[0].x;
		float y1 = rend_list_ptr->poly_ptrs[poly]->tvlist[0].y;
		float x2 = rend_list_ptr->poly_ptrs[poly]->tvlist[1].x;
		float y2 = rend_list_ptr->poly_ptrs[poly]->tvlist[1].y;
		float x3 = rend_list_ptr->poly_ptrs[poly]->tvlist[2].x;
		float y3 = rend_list_ptr->poly_ptrs[poly]->tvlist[2].y;

		IUINT32 c = (255 << 16) | (255 << 8) | 255;

		device_draw_line(&device, x1, y1, x2, y2, c);
		device_draw_line(&device, x1, y1, x3, y3, c);
		device_draw_line(&device, x2, y2, x3, y3, c);
	}
}

int gDemoIndex = 6;

void GameInit()
{
	Build_Sin_Cos_Tables();

	switch (gDemoIndex)
	{
	case 1:
		InitDemo7_1();
		break;
	case 2:
		InitDemo7_2();
		break;
	case 4:
		InitDemo7_4();
		break;

	case 6:
		InitDemo7_6();
		break;
	default:
		break;
	}
}

void GameMain()
{
	switch (gDemoIndex)
	{
	case 1:
		DrawDemo7_1();
		break;
	case 2:
		DrawDemo7_2();
		break;
	case 4:
		DrawDemo7_4();
		break;
	case 6:
		DrawDemo7_6();
		break;
	default:
		break;
	}
}

int count = 10000000;
int main(void)
{
	hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
	bool isOnlyBox = false;
	int states[] = {RENDER_STATE_TEXTURE, RENDER_STATE_COLOR, RENDER_STATE_WIREFRAME};
	int indicator = 0;
	int kbhit = 0;
	float alpha = 1;
	float pos = 3.5;

	// TCHAR *title = _T("Mini3d (software render tutorial) - ") _T("Left/Right: rotation, Up/Down: forward/backward, Space: switch state");
	TCHAR *title = _T("Wireframe");

	if (screen_init(400, 400, title))
		return -1;

	device_init(&device, 400, 400, screen_fb);

	//init_texture(&device);
	device.render_state = RENDER_STATE_TEXTURE;

	int tmp = 1;

	GameInit();

	float deltaTime = 0;

	int lastCharLen = 0;
	while (screen_exit == 0 && screen_keys[VK_ESCAPE] == 0)
	{
		time_start = GetTickCount();

		screen_dispatch();
		device_clear(&device, 1);

		GameMain();

		// screen_update();
		HDC hDC = GetDC(screen_handle);
		BitBlt(hDC, 0, 0, screen_w, screen_h, screen_dc, 0, 0, SRCCOPY);

		//DrawText(textBuffer);

		ReleaseDC(screen_handle, hDC);

		screen_dispatch();

		time_end = GetTickCount();
		deltaTime = ((float)time_end - (float)time_start);

		//   标准输出句柄
		cursorPos.X = 0;
		cursorPos.Y = 1;
		SetConsoleCursorPosition(hStdout, cursorPos);
		if (deltaTime > 0)
			std::cout << "FPS : " << 1000 / deltaTime << std::endl;

		cursorPos.X = 0;
		cursorPos.Y = 2;
		SetConsoleCursorPosition(hStdout, cursorPos);
		//if (strlen(textBuffer) < lastCharLen)
			//std::cout << "Objects Culled:                                                                                    " << std::endl;

		// std::cout << "Delta Tiem : " << 1000/deltaTime <<"ms"<< std::endl;
		// 清空这一行

		cursorPos.X = 0;
		cursorPos.Y = 2;
		SetConsoleCursorPosition(hStdout, cursorPos);
		//std::cout << textBuffer << std::endl;
		//lastCharLen = strlen(textBuffer);
	}
	CloseHandle(hStdout);
	return 0;
}
