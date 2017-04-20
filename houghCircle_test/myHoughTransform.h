﻿#pragma once
#include <stdio.h>
#include <iostream>
#include <stdarg.h>
#include "ViType.h"
#include <ipp.h>
#include <vector>
#include<algorithm>	//sort

#include "Vis_assist.h"


#define DEBUG


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

typedef struct taghoughCircle3i	//圆心半径
{
	IMG_COORD center;
	int radius;

	taghoughCircle3i()
	{
		center = { 0,0 };
		radius = 0;
	}
}houghCircle3i;


int my_vsprintf(char *format, ...);



class CVisHoughTransform
{
public:
	CVisHoughTransform();
	~CVisHoughTransform();
	
	///////////////		old		////////////////////////////////////////////	
	int detectCircle(IMG_UBBUF srcRoi, vector<houghCircle3f> &bestCircles);
	///////////////////////////////////////////////////////////////////////

	int newDetectCircle(IMG_UBBUF ubbSrc);
	vector<houghCircle3i> getBestCircles();

	void setParams(int downLevel,
		int sectors,
		float selectedRatio,
		int selectMin,
		int selectMax,
		float Ttheta,
		float Tshift,
		int localThreshMin,
		int radiusMin,
		int radiusMax,
		int voteScoreMin
	);
	
private:
	////////		old		//////////////////////////////////////////////
	int houghCircle(edgeInformation *edgeArray, IMG_INT eNum, /*IMG_WBUF mag_wBuf,*/IMG_RBUF angle_rBuf, int voteScore, int minRadius, int maxRadius, int center_Dis, vector<houghCircle3f> &bestCircles);
	int inc_if_inside(int *** H, int x, int y, int height, int width, int r);
	/////////////////////////////////////////////////////////////////////////

	IppStatus pyramid(IMG_UBBUF src,
		unsigned char* pDst,
		int &pyramid_width,
		int &pyramid_height,
		int level);

	void getGaussianKernel_dim2(IMG_LREAL ** gaus, const int size, const double sigma);
	int gaussfilter(IMG_UBBUF src, IMG_UBYTE * pDst, int kernalSize, double sigma);

	int gaussfilter_UWORD(IMG_UWBUF src, IMG_UWORD * pDst, int kernalSize, double sigma);

	int findLocalmaximum(IMG_UWBUF uwbSrc);
	void regionDFS(IMG_UWORD *pic, IMG_UWORD *label, int r, int c, int height, int width, int id, IMG_COORD * storeMax, int threshold);

private:
	//1
	int m_downLevel;

	//3
	int m_sectors;		//angle seperate
	float m_selectedRatio;
	int m_nSelectMin;
	int m_nSelectMax;		//每个扇区少量较大的梯度值点

	//4	
	float m_Ttheta;	//theta 范围内匹配
	float m_Tshift;		//判断点对的偏移

	//5 
	int m_localThreshMin;		//投票最小阈值

	//7
	int m_radiusMin;
	int m_radiusMax;

	//9
	int m_voteScoreMin;		//最终综合得分阈值(和金字塔层数有关，图像小梯度点少得分少)

	vector<vector<edgeInformation> > angleTable;
	vector<IMG_COORD> circleCenter;
	vector<houghCircle3i> bestCircles;
};




