#define WIN32_LEAN_AND_MEAN  // just say no to MFC
 
#include "header.h"
#include "Vertex.h"

// DEFINES ////////////////////////////////////////////////
 
// defines for windows
#define WINDOW_CLASS_NAME _T("WINCLASS1")
 
#define KEYDOWN(vk_code) ((GetAsyncKeyState(vk_code) & 0x8000) ? 1 : 0)
#define KEYUP(vk_code)   ((GetAsyncKeyState(vk_code) & 0x8000) ? 0 : 1)
 

// GLOBALS ////////////////////////////////////////////////
HWND      main_window_handle = NULL; // globally track main window
HINSTANCE hinstance_app      = NULL; // globally track hinstance
 
char buffer[80];                     // general printing buffer
 
#define FRAME_PER_SECOND (20)
#define TIME_IN_FRAME (1000/FRAME_PER_SECOND)
#define WIDTH 800
#define HEIGHT 600
 

//函数的声明
void drawPixel(int x, int y, COLORREF color);
void drawLine(int x1,int y1,int x2,int y2,COLORREF  color);
void sortTriangleVector2( Vertex &v1, Vertex &v2, Vertex &v3);
void drawScanLine( Vertex &v1, Vertex &v2);
void drawTriangleBottomFlat( Vertex &v1, Vertex &v2, Vertex &v3);
void drawTriangleTopFlat(Vertex &v1, Vertex &v2, Vertex &v3);
void drawMesh();
// FUNCTIONS //////////////////////////////////////////////
LRESULT CALLBACK WindowProc(HWND hwnd,
                    UINT msg,
                    WPARAM wparam,
                    LPARAM lparam)
{
   // this is the main message handler of the system
   PAINTSTRUCT   ps;     // used in WM_PAINT
   HDC           hdc;  // handle to a device context
   char buffer[80] = {0};        // used to print strings
 
   // what is the message
   switch(msg)
   { 
   case WM_CREATE:
      {
        // do initialization stuff here
        // return success
        return(0);
      } break;
 
   case WM_PAINT:
      {
        // simply validate the window
        hdc = BeginPaint(hwnd,&ps);
 
        // end painting
        EndPaint(hwnd,&ps);
 
        // return success
        return(0);
      } break;
 
   case WM_DESTROY:
      {
 
        // kill the application, this sends a WM_QUIT message
        PostQuitMessage(0);
 
        // return success
        return(0);
      } break;
 
   default:break;
 
   } // end switch
 
   // process any messages that we didn't take care of
   return (DefWindowProc(hwnd, msg, wparam, lparam));
 
} // end WinProc
 
///////////////////////////////////////////////////////////
 


void Swap(int &a, int &b){
	int temp = a;
	a = b;
	b = temp;
}

int Game_Main(void *parms = NULL, int num_parms = 0){
   DWORD dwStartTime = GetTickCount();

 
   // for now test if user is hitting ESC and send WM_CLOSE
   if (KEYDOWN(VK_ESCAPE)){
      SendMessage(main_window_handle,WM_CLOSE,0,0);
   }
  
   //HDC hdc = GetDC(main_window_handle);
 
   // draw 1000 pixels
   for (int index = 0; index < 1000; index++){
      // get random position
      int x = rand()%WIDTH;
      int y = rand()%HEIGHT;
 
      COLORREF color = RGB(rand()%255,rand()%255,rand()%255);
      //SetPixel(hdc, x,y, color);
	   drawPixel(x, y, color);
   } 
 
   // release the dc
   //ReleaseDC(main_window_handle, hdc);
 
   while(GetTickCount() - dwStartTime < TIME_IN_FRAME){
      Sleep(1);
   }
 
   return(1);
 
} 


//画点
void drawPixel(int x, int y, COLORREF color){
	HDC hdc = GetDC(main_window_handle);
	SetPixel(hdc, x,y, color);
	ReleaseDC(main_window_handle, hdc);
}

//画线
void drawLine(int x1,int y1,int x2,int y2,COLORREF  color){
	HDC hdc = GetDC(main_window_handle);

	int dx = abs(x2 - x1);
	int dy = abs(y2 - y1);
	if (dx >= dy)  {   //斜率 <= 1
		if (x1 > x2){    //令x1在左边
			
			Swap(x1,x2);
			Swap(y1,y2);
		}
		float A = y2 - y1;
		float B = x1 - x2;
		float C = x2 * y1 - x1 * y2;
		int incrementY = (y2 > y1) ? 1 : -1;
		for (int x = x1,y = y1; x <= x2; ++x){   //从左到右画
			drawPixel(x, y, color);          //画点
			//SetPixel(hdc, x,y, color);

			float k = A*(x+1) + B*(y+incrementY) + C;
			if (k * incrementY >= 0){   //点在下面面 或者 点在直线上
				y += incrementY;
			}
		}
	}

	else{                //斜率 > 1
		if (y1>y2){
			Swap(x1,x2);
			Swap(y1,y2);
		}
		float A = y2-y1;
		float B = x1-x2;
		float C = x2*y1-x1*y2;
		int incrementX = (x2>x1)?1:-1;
		for (int x=x1,y=y1;y<=y2;++y)
		{
			//SetPixel(hdc, x,y, color);

			drawPixel(x,y,color);
			float k = A*(x+incrementX)+B*(y+1)+C;
			if (k*incrementX<=0)
			{
				x += incrementX;
			}
		}
	}	
	ReleaseDC(main_window_handle, hdc);

}
 
//画网格
void drawMesh(){
	for(int x = 0; x <= WIDTH; x++){
		if(x % 2 == 0){
			drawLine(x,0,x,HEIGHT ,RGB(255,255,255));
		}
	}

	for(int y = 0; y <= HEIGHT; y++){
		if(y % 2 == 0){
			drawLine(0, y, WIDTH,y,RGB(255,255,255));
		}
	}

}
//画三角形
void drawTriangle( Vertex &v1, Vertex &v2, Vertex &v3){
	sortTriangleVector2(v1,v2,v3);
	if (v1.isxy_same(v2) && v2.isxy_same(v3))	{
		//drawPixel(v1.position_.x_, v1.position_.y_, v1.color_);
	}
	else if (v1.isxy_same(v2))	{
		//drawLine(v1.position_.x_,v1.position_.y_,v3.position_.x_,v3.position_.y_,v1.color_);
	}
	else if(v1.isxy_same(v3))	{
		//drawLine(v1.position_.x_,v1.position_.y_,v2.position_.x_,v2.position_.y_,v1.color_);
	}
	else if (v2.isxy_same(v3))	{
		//drawLine(v1.position_.x_,v1.position_.y_,v3.position_.x_,v3.position_.y_,v1.color_);
	}
	else{
		if (v1.position_.y_ == v2.position_.y_)		{
			drawTriangleTopFlat(v1,v2,v3);
		}
		else if (v2.position_.y_==v3.position_.y_)		{
			drawTriangleBottomFlat(v1,v2,v3);
		}
		else{
			float factor = (v2.position_.y_-v1.position_.y_)/(v3.position_.y_-v1.position_.y_);
			Vertex v4 = v1.interp(v3,factor);
			drawTriangleBottomFlat(v1,v2,v4);
			//Sleep(500);
			drawTriangleTopFlat(v2,v4,v3);
			cout<<"draw!!!"<<endl;
		}
	}
}
 
//对三角形顶点以y从小到大排序  min = v1, max = v3 
void sortTriangleVector2( Vertex &v1, Vertex &v2, Vertex &v3){
	if (v1.position_.y_>v2.position_.y_)	{    //max = y_v2, mid = y_v1
		swap(v1,v2);
	}
	if(v3.position_.y_<v2.position_.y_)	{      //max = y_v3, mid = y_v2
		swap(v2,v3);
	}
	if (v1.position_.y_>v2.position_.y_)	{    //mid = y_v2, min = y_v1
		swap(v1,v2);
	}
}

//绘制插值（颜色变化）直线
void drawScanLine( Vertex &v1, Vertex &v2)
{
	if (v1.position_.x_>v2.position_.x_)
	{
		swap(v1,v2);
	}

	//换成加法后速度快了好多
	int x_start = v1.position_.x_+0.5;
	int x_end = v2.position_.x_+0.5;
	bool isZero = x_end-x_start?false:true;
	float df = 1.0f/(x_end-x_start);
	//depth
	float dd = (v2.position_.w_-v1.position_.w_)*df;
	float one_depth = v1.position_.w_;
	//u,v
	float w = v1.position_.w_;
	float dw = (v2.position_.w_-v1.position_.w_)*df;

	float uw = v1.u_*v1.position_.w_ ;
	float vw = v1.v_*v1.position_.w_;
	float duw = (v2.u_*v2.position_.w_-v1.u_*v1.position_.w_)*df;
	float dvw = (v2.v_*v2.position_.w_-v1.v_*v1.position_.w_)*df;

	//light
	float r = v1.light_.x_;
	float dr = (v2.light_.x_-v1.light_.x_)*df;
	float g = v1.light_.y_;
	float dg = (v2.light_.y_-v1.light_.y_)*df;
	float b = v1.light_.z_;
	float db = (v2.light_.z_-v1.light_.z_)*df;
	for (int i = x_start; i < x_end; ++i){		
		float factor = min(1.0f,df*(i-x_start));
		AColor color = v1.color_.interp(v2.color_,factor);
		drawPixel(i,v1.position_.y_,RGB(color.r_,color.g_,color.b_));


		one_depth += dd;
		uw += duw;
		vw += dvw;
		w += dw;
		r += dr;
		g += dg;
		b += db;

	}

}


// 绘制底平三角形	 v1为上顶点	(v2左下，v3右下)		
void drawTriangleBottomFlat( Vertex &v1, Vertex &v2, Vertex &v3)
{
	if (v2.position_.x_ > v3.position_.x_){
		swap(v2,v3);
	}
	int startY = v1.position_.y_+0.5;   //四舍五入，>=0.5进1，<0.5归0
	int endY = v2.position_.y_+0.5;
	for (int y = startY; y < endY; y++){
		float factor = 0;
		if (startY - endY != 0)		{
			factor = (float(float(y) + 0.5 - v1.position_.y_))/
				(v2.position_.y_ - v1.position_.y_);
		}
		Vertex vl = v1.interp(v2,factor);
		Vertex vr = v1.interp(v3,factor);

		drawScanLine(vl,vr);
		//Sleep(100);
		//cout<<"draw line"<<endl;
	}
}

//绘制顶平三角形		v3为底顶点	(v1左上，v2右上)		
void drawTriangleTopFlat(Vertex &v1, Vertex &v2, Vertex &v3){
	if (v1.position_.x_ > v2.position_.x_){
		swap(v1,v2);
	}
	int startY = v2.position_.y_+0.5;
	int endY = v3.position_.y_+0.5;
	for (int y=startY;y<endY;y++)
	{
		float factor =0;
		if (startY-endY!=0)
		{
			factor = (float(float(y)+0.5-v2.position_.y_))/(v3.position_.y_-v2.position_.y_);
		}
		Vertex vl = v1.interp(v3,factor);
		Vertex vr = v2.interp(v3,factor);
		drawScanLine(vl,vr);
		//Sleep(100);

		//cout<<"draw line"<<endl;
	}
}

int Game_Init(void *parms = NULL, int num_parms = 0){
   return(1);
}
 
 
int Game_Shutdown(void *parms = NULL, int num_parms = 0){
   return(1);
}
 

int WINAPI WinMain( HINSTANCE hinstance,
              HINSTANCE hprevinstance,
              LPSTR lpcmdline,
              int ncmdshow)
{
 
   WNDCLASSEX winclass; // this will hold the class we create
   HWND     hwnd; // generic window handle
   MSG        msg;    // generic message
   HDC        hdc;      // graphics device context
 
   // first fill in the window class stucture
   winclass.cbSize         = sizeof(WNDCLASSEX);
   winclass.style        = CS_DBLCLKS | CS_OWNDC |
      CS_HREDRAW | CS_VREDRAW;
   winclass.lpfnWndProc  = WindowProc;
   winclass.cbClsExtra   = 0;
   winclass.cbWndExtra   = 0;
   winclass.hInstance    = hinstance;
   winclass.hIcon        = LoadIcon(NULL, IDI_APPLICATION);
   winclass.hCursor    = LoadCursor(NULL, IDC_ARROW);
   winclass.hbrBackground = (HBRUSH)GetStockObject(BLACK_BRUSH);
   winclass.lpszMenuName = NULL;
   winclass.lpszClassName = WINDOW_CLASS_NAME;
   winclass.hIconSm        = LoadIcon(NULL, IDI_APPLICATION);
 
   // save hinstance in global
   hinstance_app = hinstance;
 
   // register the window class
   if (!RegisterClassEx(&winclass))
      return(0);
 
   // create the window
   if (!(hwnd = CreateWindowEx(NULL,                  // extended style
      WINDOW_CLASS_NAME,     // class
      _T("Show Point 0.1"), // title
      WS_OVERLAPPEDWINDOW | WS_VISIBLE,
      0,0,    // initial x,y
      WIDTH,HEIGHT,  // initial width, height
      NULL,   // handle to parent
      NULL,   // handle to menu
      hinstance,// instance of this application
      NULL))) // extra creation parms
      return(0);
 
   // save main window handle
   main_window_handle = hwnd;
 
   // initialize game here
   Game_Init();

   // enter main event loop
   while(TRUE){
      // test if there is a message in queue, if so get it
      if (PeekMessage(&msg,NULL,0,0,PM_REMOVE)){
        // test if this is a quit
        if (msg.message == WM_QUIT)
           break;
 
        // translate any accelerator keys
        TranslateMessage(&msg);
        // send the message to the window proc
        DispatchMessage(&msg);
      } 
	  //Game_Main();


	  //画东西
	  //drawMesh();

	  Vector3 v1 = Vector3(600,400,0);
	  Vector3 v2 = Vector3(800,500,0);
	  Vector3 v3 = Vector3(700,300,0);
	  drawTriangle(
		  Vertex(v1, AColor(0,0,255,0),0,0),
		  Vertex(v2,AColor(0,255,0,0),0,0),
		  Vertex(v3,AColor(0,0,0,255),0,0));   
	  drawLine(0,0,199,135,RGB(255,255,255));
}

   Game_Shutdown();
   return(msg.wParam);
} 