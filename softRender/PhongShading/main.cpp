
#include "device.h"

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

// typedef unsigned int IUINT32;

int CMID(int x, int min, int max) { return (x < min) ? min : ((x > max) ? max : x); }


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

void jjjDrawTextOnScreen(char *text, int x, int y)
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
RENDERLIST4DV2 rend_list2; // the render list
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

void Initmo7_4();
void DrawDemo7_4();

void Initmo7_6();

void InitDemo7_1()
{
}
void DrawDemo7_1()
{

}

char *work_string;
char *textBuffer; // used to print text

#define KEY_DOWN(vk_code) ((GetAsyncKeyState(vk_code) & 0x8000) ? 1 : 0)

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


void InitDemo8_4();
void DrawDemo8_4();


void InitDemo8_4()
{

}

void DrawDemo8_4()
{
}

void InitDemo8_6();
void DrawDemo8_6();

void InitDemo8_6()
{
}

void DrawDemo8_6()
{
}	

void InitDemo9_2();
void DrawDemo9_2();

OBJECT4DV2 obj_constant_water,
    obj_flat_water,
    obj_gouraud_water,
	obj_constant_light;

RGBAV1 white, gray, black, red, green, blue;
void InitDemo9_2()
{
	// POINT4D cam_pos = {0, 0, 0, 1};
	POINT4D cam_pos = {105, 0, 13, 1};
	POINT4D cam_target = {0, 0, 0, 1};
	VECTOR4D cam_dir = {0, -42, 0, 1};

	// all your initialization code goes here...
	VECTOR4D vscale = {1.0, 1.0, 1.0, 1},
			 vpos = {0, 0, 0, 1},
			 vrot = {0, 0, 0, 1};

   Init_CAM4DV1(&cam,            // the camera object
                CAM_MODEL_EULER, // the euler model
                &cam_pos,        // initial camera position
                &cam_dir,        // initial camera angles
                &cam_target,     // no target
                200.0,           // near and far clipping planes
                12000.0,
                120.0,        // field of view in degrees
                WINDOW_WIDTH, // size of final screen viewport
                WINDOW_HEIGHT);

	VECTOR4D_INITXYZ(&vscale, 10.00, 10.00, 10.00);

	// Load_OBJECT4DV2_COB(&obj_constant_water, "./cob/water_constant_01.cob",
	// Load_OBJECT4DV2_COB(&obj_constant_water, "./cob/cube_constant_01.cob",						&vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD);
   Load_OBJECT4DV2_COB(&obj_constant_water, "./cob/cube_flat_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD);


//    Load_OBJECT4DV2_COB(&obj_flat_water, "./cob/cube_flat_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD);
   // load flat shaded water
   VECTOR4D_INITXYZ(&vscale, 10.00, 10.00, 10.00);
//    Load_OBJECT4DV2_COB(&obj_flat_water, "./cob/water_flat_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD);
   Load_OBJECT4DV2_COB(&obj_flat_water, "./cob/cube_gouraud_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD); //修改颜色后的水分子

   // load gouraud shaded water
   VECTOR4D_INITXYZ(&vscale, 10.00, 10.00, 10.00);
//    Load_OBJECT4DV2_COB(&obj_gouraud_water, "./cob/water_flat_01_gouraud.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD); //修改颜色后的水分子
   Load_OBJECT4DV2_COB(&obj_gouraud_water, "./cob/cube_gouraud_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD); //修改颜色后的水分子

   VECTOR4D_INITXYZ(&vscale, 5.00, 5.00, 5.00);

   Load_OBJECT4DV2_COB(&obj_constant_light, "./cob/cube_flat_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD);
   Reset_Lights_LIGHTV1();

   // create some working colors
   white.rgba = _RGBA32BIT(255, 255, 255, 0);
   gray.rgba = _RGBA32BIT(100, 100, 100, 0);
   black.rgba = _RGBA32BIT(0, 0, 0, 0);
   red.rgba = _RGBA32BIT(255, 0, 0, 0);
   green.rgba = _RGBA32BIT(0, 255, 0, 0);
   blue.rgba = _RGBA32BIT(0, 0, 255, 0);

   // ambient light
   Init_Light_LIGHTV1(AMBIENT_LIGHT_INDEX,
					  LIGHTV1_STATE_ON,		// turn the light on
					  LIGHTV1_ATTR_AMBIENT, // ambient light type
					//   gray, black, black,   // color for ambient term only
					  black, black, black,   // color for ambient term only
					  NULL, NULL,			// no need for pos or dir
					  0, 0, 0,				// no need for attenuation
					  0, 0, 0);				// spotlight info NA

   VECTOR4D dlight_dir = {-1, 0, -1, 0};

   // directional light
   Init_Light_LIGHTV1(INFINITE_LIGHT_INDEX,
					  LIGHTV1_STATE_ON,		 // turn the light on
					  LIGHTV1_ATTR_INFINITE, // infinite light type
					//   black, gray, black,	// color for diffuse term only
					  gray, gray, gray,	// color for diffuse term only
					  NULL, &dlight_dir,	 // need direction only
					  0, 0, 0,				 // no need for attenuation
					  0, 0, 0);				 // spotlight info NA

   VECTOR4D plight_pos = {0, 200, 0, 0};

   // point light
   Init_Light_LIGHTV1(POINT_LIGHT_INDEX,
					  LIGHTV1_STATE_ON,	// turn the light on
					  LIGHTV1_ATTR_POINT,  // pointlight type
					//   black, green, black, // color for diffuse term only
					  black, white, white, // color for diffuse term only
					  &plight_pos, NULL,   // need pos only
					  0, .001, 0,		   // linear attenuation only
					  0, 0, 1);			   // spotlight info NA

   VECTOR4D slight2_pos = {0, 200, 0, 0};
   VECTOR4D slight2_dir = {-1, 0, -1, 0};

   // spot light2
//    Init_Light_LIGHTV1(SPOT_LIGHT2_INDEX,
// 					  LIGHTV1_STATE_ON,			  // turn the light on
// 					  LIGHTV1_ATTR_SPOTLIGHT2,	// spot light type 2
// 					//   black, red, black,		  // color for diffuse term only
// 					  black, white, black,		  // color for diffuse term only
// 					  &slight2_pos, &slight2_dir, // need pos only
// 					  0, .001, 0,				  // linear attenuation only
// 					  0, 0, 1);
}
void DrawDemo9_2()
{
	// Sleep(20);
	std::cout<<"sdfkjl"<<std::endl;
	SetConsoleCursorPosition(hStdout, cursorPos);
	// std::cout<<rend_list.num_polys<<std::endl;
	cursorPos.X = 0;
	cursorPos.Y = 2;
	SetConsoleCursorPosition(hStdout, cursorPos);
	std::cout << "                                                                                                                     " << std::endl;
	SetConsoleCursorPosition(hStdout, cursorPos);
	std::cout << "cam pos/dir :" 
	<< cam.pos.x << ", " << cam.pos.y << ", " << cam.pos.z <<"/"
	<< cam.dir.x << ", " << cam.dir.y << ", " << cam.dir.z 
	<< std::endl;

	cursorPos.X = 0;
	cursorPos.Y = 3;
	SetConsoleCursorPosition(hStdout, cursorPos);
	std::cout << "                                                                                                                     " << std::endl;
	SetConsoleCursorPosition(hStdout, cursorPos);
	std::cout << "plight_pos :" 
	<< lights[POINT_LIGHT_INDEX].pos.x << ", " << lights[POINT_LIGHT_INDEX].pos.y << ", " <<  lights[POINT_LIGHT_INDEX].pos.z<<" /"
	<< lights[POINT_LIGHT_INDEX].dir.x << ", " << lights[POINT_LIGHT_INDEX].dir.y << ", " <<  lights[POINT_LIGHT_INDEX].dir.z
	<< std::endl;

	cursorPos.X = 0;
	cursorPos.Y = 4;
	SetConsoleCursorPosition(hStdout, cursorPos);
	std::cout << "                                                                                                                     " << std::endl;
	SetConsoleCursorPosition(hStdout, cursorPos);
	std::cout << "spotlight_pos : "
	 << lights[SPOT_LIGHT2_INDEX].pos.x << ", " << lights[SPOT_LIGHT2_INDEX].pos.y << ", " <<  lights[SPOT_LIGHT2_INDEX].pos.z
	 <<"/ "<< lights[SPOT_LIGHT2_INDEX].dir.x << ", " << lights[SPOT_LIGHT2_INDEX].dir.y << ", " <<  lights[SPOT_LIGHT2_INDEX].dir.z
	 << std::endl;

	cursorPos.X = 0;
	cursorPos.Y = 5;
	SetConsoleCursorPosition(hStdout, cursorPos);
	std::cout << "                                                                                                                     " << std::endl;
	SetConsoleCursorPosition(hStdout, cursorPos);
	std::cout <<"obj_constant_water.world_pos :"<< obj_constant_water.world_pos.x << ", " << obj_constant_water.world_pos.y << ", " <<  obj_constant_water.world_pos.z<< std::endl;


	cursorPos.X = 0;
	cursorPos.Y = 6;
	SetConsoleCursorPosition(hStdout, cursorPos);
	std::cout << "                                                                                                                     " << std::endl;
	SetConsoleCursorPosition(hStdout, cursorPos);
	std::cout <<"obj_flat_water.world_pos :"<< obj_flat_water.world_pos.x << ", " << obj_flat_water.world_pos.y << ", " <<  obj_flat_water.world_pos.z<< std::endl;


	cursorPos.X = 0;
	cursorPos.Y = 7;
	SetConsoleCursorPosition(hStdout, cursorPos);
	std::cout << "                                                                                                                     " << std::endl;
	SetConsoleCursorPosition(hStdout, cursorPos);
	std::cout <<"obj_gouraud_water.world_pos :"<< obj_gouraud_water.world_pos.x << ", " << obj_gouraud_water.world_pos.y << ", " <<  obj_gouraud_water.world_pos.z<< std::endl;

	static MATRIX4X4 mrot; // general rotation matrix

	static float view_angle = 0;
	static float camera_distance = 6000;
	static VECTOR4D pos = {0, 0, 0, 0};
	static float tank_speed = 10;
	static float turning = 0;

	char work_string[256]; // temp string

	int index; // looping var

	static float plight_ang = 0, slight_ang = 0; // angles for light motion

	// move point light source in ellipse around game world
	lights[POINT_LIGHT_INDEX].pos.x = 1000 * Fast_Cos(plight_ang);
	lights[POINT_LIGHT_INDEX].pos.y = 100;
	// lights[POINT_LIGHT_INDEX].pos.z = 1000 * Fast_Sin(plight_ang);

	// if ((plight_ang += 3) > 360)
	// 	plight_ang = 0;

	// move spot light source in ellipse around game world
	lights[SPOT_LIGHT2_INDEX].pos.x = 1000 * Fast_Cos(slight_ang);
	lights[SPOT_LIGHT2_INDEX].pos.y = 200;
	lights[SPOT_LIGHT2_INDEX].pos.z = 1000 * Fast_Sin(slight_ang);

	// if ((slight_ang -= 5) < 0)
	// 	slight_ang = 360;
	//Reset_RENDERLIST4DV1(&rend_list);
	Reset_RENDERLIST4DV2(&rend_list2);
	if (KEY_DOWN(VK_NUMPAD8))
	{
		// move forward
		cam.pos.y += 1;
	} // end if

	if (KEY_DOWN(VK_NUMPAD2))
	{
		cam.pos.y -= 1;
	} // end if

	if (KEY_DOWN(VK_NUMPAD9))
	{
		lights[POINT_LIGHT_INDEX].pos.z += 10;
		// move forward
	} // end if

	if (KEY_DOWN(VK_NUMPAD3))
	{
		lights[POINT_LIGHT_INDEX].pos.z -= 10;
	} // end if
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
	Build_CAM4DV1_Matrix_Euler(&cam, CAM_ROT_SEQ_ZYX);
	static float x_ang = 0, y_ang = 0, z_ang = 0;

	Reset_OBJECT4DV2(&obj_constant_water);
	// set position of constant shaded water
	obj_constant_water.world_pos.x = -50;
	obj_constant_water.world_pos.y = 0;
	obj_constant_water.world_pos.z = 120;
	// generate rotation matrix around y axis
	Build_XYZ_Rotation_MATRIX4X4(x_ang, y_ang, z_ang, &mrot);
	// rotate the local coords of the object
	Transform_OBJECT4DV2(&obj_constant_water, &mrot, TRANSFORM_LOCAL_TO_TRANS, 1);
	// perform world transform
	Model_To_World_OBJECT4DV2(&obj_constant_water, TRANSFORM_TRANS_ONLY);
	// insert the object into render list
	Insert_OBJECT4DV2_RENDERLIST4DV2(&rend_list2, &obj_constant_water, 0);


	Reset_OBJECT4DV2(&obj_constant_light);
	// set position of constant shaded water
	obj_constant_light.world_pos.x = lights[POINT_LIGHT_INDEX].pos.x;
	obj_constant_light.world_pos.y = lights[POINT_LIGHT_INDEX].pos.y;
	obj_constant_light.world_pos.z = lights[POINT_LIGHT_INDEX].pos.z;
	// generate rotation matrix around y axis
	Build_XYZ_Rotation_MATRIX4X4(x_ang, y_ang, z_ang, &mrot);
	// rotate the local coords of the object
	Transform_OBJECT4DV2(&obj_constant_light, &mrot, TRANSFORM_LOCAL_TO_TRANS, 1);
	// perform world transform
	Model_To_World_OBJECT4DV2(&obj_constant_light, TRANSFORM_TRANS_ONLY);
	// insert the object into render list
	Insert_OBJECT4DV2_RENDERLIST4DV2(&rend_list2, &obj_constant_light, 0);


	Reset_OBJECT4DV2(&obj_flat_water);
	// set position of constant shaded water
	obj_flat_water.world_pos.x = 0;
	obj_flat_water.world_pos.y = 0;
	obj_flat_water.world_pos.z = 120;
	// generate rotation matrix around y axis
	Build_XYZ_Rotation_MATRIX4X4(x_ang, y_ang, z_ang, &mrot);
	// rotate the local coords of the object
	Transform_OBJECT4DV2(&obj_flat_water, &mrot, TRANSFORM_LOCAL_TO_TRANS, 1);
	// perform world transform
	Model_To_World_OBJECT4DV2(&obj_flat_water, TRANSFORM_TRANS_ONLY);
	// insert the object into render list
	Insert_OBJECT4DV2_RENDERLIST4DV2(&rend_list2, &obj_flat_water, 0);

	//////////////////////////////////////////////////////////////////////////
	// gouraud shaded water

	// reset the object (this only matters for backface and object removal)
	Reset_OBJECT4DV2(&obj_gouraud_water);

	// set position of constant shaded water
	obj_gouraud_water.world_pos.x = 50;
	obj_gouraud_water.world_pos.y = 0;
	obj_gouraud_water.world_pos.z = 120;

	// generate rotation matrix around y axis
	Build_XYZ_Rotation_MATRIX4X4(x_ang, y_ang, z_ang, &mrot);

	// rotate the local coords of the object
	Transform_OBJECT4DV2(&obj_gouraud_water, &mrot, TRANSFORM_LOCAL_TO_TRANS, 1);

	// perform world transform
	Model_To_World_OBJECT4DV2(&obj_gouraud_water, TRANSFORM_TRANS_ONLY);
	// 丢到了 obj->vlist_trans[vertex].v

	// insert the object into render list
	Insert_OBJECT4DV2_RENDERLIST4DV2(&rend_list2, &obj_gouraud_water, 0);

	RENDERLIST4DV2_PTR rend_list_ptr = &rend_list2;
	std::cout<<rend_list_ptr->num_polys<<std::endl;
	// update rotation angles
	// if ((x_ang += 1) > 360)
	// 	x_ang = 0;
	// if ((y_ang += 2) > 360)
		// y_ang = 0;
	// if ((z_ang += 3) > 360)
	// 	z_ang = 0;

	// remove backfaces
		Remove_Backfaces_RENDERLIST4DV2(&rend_list2, &cam);

	// light scene all at once
		Light_RENDERLIST4DV2_World16(&rend_list2, &cam, lights, 4);

	// apply world to camera transform
	World_To_Camera_RENDERLIST4DV2(&rend_list2, &cam);

	// sort the polygon list (hurry up!)
		Sort_RENDERLIST4DV2(&rend_list2, SORT_POLYLIST_AVGZ);

	// apply camera to perspective transformation
	Camera_To_Perspective_RENDERLIST4DV2(&rend_list2, &cam);

	// apply screen transform
	Perspective_To_Screen_RENDERLIST4DV2(&rend_list2, &cam);


	
	//  RENDERLIST4DV1_PTR rend_list_ptr = &rend_list;

	POLYF4DV2 face; // temp face used to render polygon

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

		int color = rend_list_ptr->poly_ptrs[poly]->lit_color[0];
		if ((rend_list_ptr->poly_ptrs[poly]->attr & POLY4DV2_ATTR_SHADE_MODE_FLAT) ||
			(rend_list_ptr->poly_ptrs[poly]->attr & POLY4DV2_ATTR_SHADE_MODE_CONSTANT))
		{

			IUINT32 c = (255 << 16) | (255 << 8) | 255;

		//  std::cout<<"color = "<<color<<std::endl;
		//  int color = rend_list_ptr->poly_ptrs[poly]->color;
		//  rend_list->poly_ptrs[poly]->lit_color[0]
			 // DrawTrianglePureColor(&device, x1, y1, x2, y2, x3, y3, color);
		 DrawTrianglePureColor2(&device, x1, y1, x2, y2, x3, y3, color);
		//  device_draw_line(&device, x1, y1, x2, y2, c);
		//  device_draw_line(&device, x1, y1, x3, y3, c);
		//  device_draw_line(&device, x2, y2, x3, y3, c);
		}
		else if (rend_list_ptr->poly_ptrs[poly]->attr & POLY4DV2_ATTR_SHADE_MODE_GOURAUD)
		{
			// {andre take advantage of the data structures later..}
			// set the vertices
			face.tvlist[0].x = (int)rend_list_ptr->poly_ptrs[poly]->tvlist[0].x;
			face.tvlist[0].y = (int)rend_list_ptr->poly_ptrs[poly]->tvlist[0].y;
			face.lit_color[0] = rend_list_ptr->poly_ptrs[poly]->lit_color[0];

			face.tvlist[1].x = (int)rend_list_ptr->poly_ptrs[poly]->tvlist[1].x;
			face.tvlist[1].y = (int)rend_list_ptr->poly_ptrs[poly]->tvlist[1].y;
			face.lit_color[1] = rend_list_ptr->poly_ptrs[poly]->lit_color[1];

			face.tvlist[2].x = (int)rend_list_ptr->poly_ptrs[poly]->tvlist[2].x;
			face.tvlist[2].y = (int)rend_list_ptr->poly_ptrs[poly]->tvlist[2].y;
			face.lit_color[2] = rend_list_ptr->poly_ptrs[poly]->lit_color[2];

			// draw the gouraud shaded triangle
			// Draw_Gouraud_Triangle16(&face, video_buffer, lpitch);

			// Draw_Gouraud_Triangle16(&device, &face);

			DrawPhongTriangle(&device, &cam,&face, rend_list_ptr->poly_ptrs[poly], lights);

			//     std::cout<<"...tu0 = "<< tu0    <<"tv0 = "<< tv0     <<"tw0 = "<< tw0
            //  <<"tu1 = "<< tu1     <<"tv1 = "<< tv1     <<"tw1 = "<< tw1
            //  <<"tu2 = "<< tu2     <<"tv2 = "<< tv2     <<"tw2 = "<< tw2<< std::endl;

			//          std::cout<<"...tpu0 = "<< tpu0    <<"tpv0 = "<< tpv0     <<"tpw0 = "<< tpw0
            //  <<"tpu1 = "<< tpu1     <<"tpv1 = "<< tpv1     <<"tpw1 = "<< tpw1
            //  <<"tpu2 = "<< tpu2     <<"tpv2 = "<< tpv2     <<"tpw2 = "<< tpw2<< std::endl;

			// DrawTrianglePureColor2(&device, x1, y1, x2, y2, x3, y3, color);
		} // end if gouraud
	}
}


int gDemoIndex = 92;

void GameInit()
{
	Build_Sin_Cos_Tables();

	switch (gDemoIndex)
	{
	case 1:
		InitDemo7_1();
		break;
	case 84:
		InitDemo8_4();
		break;
	case 86:
		InitDemo8_6();
		break;
	case 92:
		InitDemo9_2();
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
	case 84:
		DrawDemo8_4();
		break;
	case 86:
		Sleep(20);
		DrawDemo8_6();
		break;
	case 92:
		Sleep(20);
		DrawDemo9_2();
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

	if (screen_init(600, 600, title))
		return -1;

	device_init(&device, 600, 600, screen_fb);

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
