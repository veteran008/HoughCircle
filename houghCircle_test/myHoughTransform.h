#pragma once
#include <stdio.h>
#include <iostream>
#include <stdarg.h>
#include "ViType.h"
#include <ipp.h>
#include <vector>

#include "Vis_assist.h"

//#include <opencv2/opencv.hpp>

//#define DEBUG

//using namespace cv;
using namespace std;

/***			memory leak					****/
#ifdef _DEBUG  
#define new  new(_NORMAL_BLOCK, __FILE__, __LINE__)  
#endif 
#include <crtdbg.h>  
//在入口函数中包含 _CrtDumpMemoryLeaks();    
//即可检测到内存泄露
/**********************************************/

#define PI 3.14159265386
typedef struct
{
	IMG_COORD xyInteger; //像素点
	IMG_RCOORD xyDecimal;//亚像素点
	int gradient;
	float angle;
}edgeInformation;//边缘点

typedef struct taghoughCircle3f	//圆心半径
{
	float centerX;
	float centerY;
	float radius;

	taghoughCircle3f()
	{
		centerX = 0;
		centerY = 0;
		radius = 0;
	}
}houghCircle3f;


int my_vsprintf(char *format, ...);



#define PYRAMID_DOWN_LEVEL 2
class CVisHoughTransform
{
public:
	CVisHoughTransform();
	~CVisHoughTransform();

	int detectCircle(IMG_UBBUF srcRoi, vector<houghCircle3f> &bestCircles);

	int newDetectCircle(IMG_UBBUF ubbSrc);

	vector<houghCircle3f> getBestCircles() const { return bestCircles; }
	void setBestCircles(vector<houghCircle3f> val) { bestCircles = val; }

private:
	int houghCircle(edgeInformation *edgeArray, IMG_INT eNum, /*IMG_WBUF mag_wBuf,*/IMG_RBUF angle_rBuf, int voteScore, int minRadius, int maxRadius, int center_Dis, vector<houghCircle3f> &bestCircles);
	int inc_if_inside(int *** H, int x, int y, int height, int width, int r);

	IppStatus pyramid(IMG_UBBUF src,
		unsigned char* pDst,
		int &pyramid_width,
		int &pyramid_height,
		int level);

private:
	vector<houghCircle3f> bestCircles;
};




