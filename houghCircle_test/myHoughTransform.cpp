#include "stdafx.h"
#include "myHoughTransform.h"


const int printBufSize = 256;
int my_vsprintf(char *format, ...)
{
#ifdef DEBUG

	char buf[printBufSize];
	va_list st;
	va_start(st, format);
	//vsprintf(buf, format, st);
	vsnprintf(buf, printBufSize, format, st);
	/***************************************************************************/
	/*       函数名: vsprintf                                       　　　　　　
	/*       功 能: 送格式化输出到串中                                         　　
	/*       返回值: 正常情况下返回生成字串的长度(除去\0),错误情况返回负值
	/*       用 法: int vsprintf(char *string, char *format, va_list param);
	/*                将param 按格式format写入字符串string中
	/*       注: 该函数会出现内存溢出情况,建议使用vsnprintf                　　 　  　                                                     */
	/***************************************************************************/
	va_end(st);
	cout << buf;// << endl;
#endif

	return 0;
}

//2、边缘检测

///////////////////////////////////////////////////////////////////////////
//VisEdge_detection功能说明
//Input
//srcRoi   输入图像
//roiSize  输入图像的尺寸
//threshold  梯度强度的阈值（理论上是大于0小于1250的整数）。如果用户不知道阈值设为多少合适，可以输入0（算法自动获取合适的阈值）。
//
//output
//dstRoi  梯度强度
//edgeInformation *&edgeArray  边缘点信息，包括像素坐标、亚像素坐标、梯度强度、角度
//
//函数返回
//正常情况下返回1；
//如果用户输入阈值小于0，函数返回-1；
//如果输入参数不正确，包括图像尺寸不正确、srcRoi大小与尺寸不符合，函数返回-1。
//Author：姜贺/20170227
////////////////////////////////////////////////

////////////////////////	namespace ommTool	/////////////////////////
namespace ommTool 
{
	IMG_VVOID SobelFilter_8u16s_C1_5x5(IMG_UBYTE * pSrc, IppiSize roiSize, Ipp16s *pDst, Ipp32f *pAngle)
	{
		IppiMaskSize mask = ippMskSize5x5;
		IppiBorderType bordertype = ippBorderRepl; //Border is replicated from the edge pixels
		Ipp16s *pHoriz, *pVert;

		int srcStep = roiSize.width * sizeof(Ipp8u);

		int dstStep = roiSize.width * sizeof(Ipp16s);
		int angleStep = roiSize.width * sizeof(Ipp32f);
		int bufferSize;
		int bufLen = roiSize.width * roiSize.height;
		IppStatus statusVert, statusHoriz, status;
		Ipp8u *pBuffer;
		IppNormType normType = ippNormL2;//input gradient magnitude

		pVert = (Ipp16s *)malloc(sizeof(Ipp16s)*bufLen);
		pHoriz = (Ipp16s *)malloc(sizeof(Ipp16s)*bufLen);

		ippiGradientVectorGetBufferSize(roiSize, mask, ipp16s, 1, &bufferSize);
		pBuffer = (Ipp8u *)malloc(bufferSize);
		ippiGradientVectorSobel_8u16s_C1R(pSrc, srcStep, pVert, dstStep, pHoriz, dstStep, pDst, dstStep, pAngle, angleStep, roiSize, mask, normType, bordertype, NULL, pBuffer);

		free(pVert);
		free(pHoriz);
		free(pBuffer);
	}

	//给定一个值threshold，计算最大类间方差
	IMG_REAL getIntraClassVariance(Ipp16s* src, int srcRows, int srcCols, int &varTh)//int &varian)
	{
		//intra-class variance
		float varian = 0;

		int sumPixel = srcRows*srcCols;
		int sumGrayValue = 0;
		int average = 0;

		int sumApixel = 0;
		double PA = 0;
		int sumAgray = 0;
		int averageA = 0;

		int sumBpixel = 0;
		double PB = 0;
		int sumBgray = 0;
		int averageB = 0;

		for (int i = 0; i < sumPixel; i++)
		{
			sumGrayValue = sumGrayValue + src[i];
			if (src[i] < varTh)
			{
				sumApixel++;
				sumAgray = sumAgray + src[i];
			}
		}

		average = sumGrayValue / sumPixel;
		PA = (double)sumApixel / sumPixel;
		if (sumApixel > 0)
		{
			averageA = sumAgray / sumApixel;
		}
		else
		{
			averageA = 0;
		}

		sumBpixel = sumPixel - sumApixel;
		PB = 1.0 - PA;
		sumBgray = sumGrayValue - sumAgray;
		if (sumBpixel > 0)
		{
			averageB = sumBgray / sumBpixel;
		}
		else
		{
			averageB = 0;
		}

		//ICV = PA?(MA?M)2 + PB?(MB?M)2
		varian = PA * (averageA - average) * (averageA - average) + PB * (averageB - average) * (averageB - average);

		return varian;
	}

	IMG_INT VisEdge_detection(IMG_UBYTE *srcRoi, IMG_SIZE roiSize, int threshold, IMG_WORD *dstRoi, IMG_UBYTE *dstRoiE, Ipp32f *angAll,edgeInformation *&edgeArray, IMG_INT &eNum)
	{
		//如果阈值小于0，函数直接返回-1
		if (threshold < 0)
		{
			return -1;
		}

		int roiRows = roiSize.height;
		int roiCols = roiSize.width;

		//角度信息
//		Ipp32f *angAll;
		//angAll = (Ipp32f*)malloc(roiRows*roiCols * sizeof(Ipp32f));

		std::vector<edgeInformation> edgeInfor;
		edgeInformation edInf;

		int k = 0;//记录边缘点的个数
		Ipp16u k1;//抛物线拟合的三个已知点
		Ipp16u k2;
		Ipp16u k3;
		float deci;//抛物线拟合顶点的小数部分，即对应的亚像素
		float sumx = 0;//边缘点的x坐标之和
		float sumy = 0;
		int numberChannels = 1; //the source image is single channel

		IppiSize dstRoiSize = { roiCols,roiRows };

		SobelFilter_8u16s_C1_5x5(srcRoi, dstRoiSize, dstRoi, angAll);


		//把角度由[-PI，PI]变为[0，360]
		/*for (int i = 0; i < roiRows; i++)
		{
		for (int j = 0; j < roiCols; j++)
		{
		angAll[j + i * roiCols] = (float)(180 - angAll[j + i * roiCols] / PI * 180);
		}
		}*/

		/*
		FILE *sx;
		sx = fopen("E:\\ProjectMy\\a1sobel.txt", "w");
		FILE *ang;
		ang = fopen("E:\\ProjectMy\\a1ang.txt", "w");
		for (int i = 0; i<roiRows; i++)
		{
		for (int j = 0; j < roiCols; j++)
		{
		fprintf(sx, "%d   ", dstRoi[j+i*roiCols]);
		fprintf(ang, "%f   ", angAll[j+i*roiCols]);
		}
		fprintf(sx,"\n");
		fprintf(ang, "\n");
		}
		fclose(sx);
		fclose(ang);
		*/
		//自动获取梯度强度的阈值
		//什么情况下算是没有输入阈值呢？
		if (threshold == 0)
		{
			//Otsu法，遍历所有的灰度值，从1到255，使intra-class invariance最大的那个值，即为要求的阈值
			int varian = 0;
			int temp = 0;
			for (int p = 1; p < 800; p++)
			{
				temp = getIntraClassVariance(dstRoi, roiRows, roiCols, p);
				if (varian < temp)
				{
					varian = temp;
					threshold = p;
				}
			}
		}
		//printf("%d\n",threshold);

		//到亚像素
		for (int i = 1; i<roiRows - 1; i++)
		{
			for (int j = 1; j<roiCols - 1; j++)
			{
				if (dstRoi[j + i*roiCols] > threshold)
				{
					angAll[j + i * roiCols] = (float)180 - angAll[j + i * roiCols] / PI * 180;
					if ((angAll[j + i*roiCols]>22.5) && (angAll[j + i*roiCols] < 67.5))
					{
						if ((dstRoi[j + i*roiCols] > dstRoi[j - 1 + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j + 1) + (i + 1)*roiCols]))
						{
							k1 = dstRoi[j - 1 + (i - 1)*roiCols];
							k2 = dstRoi[j + i*roiCols];
							k3 = dstRoi[(j + 1) + (i + 1)*roiCols];
							deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));


							edInf.xyInteger.x = j;
							edInf.xyInteger.y = i;
							edInf.xyDecimal.x = j + deci;
							edInf.xyDecimal.y = i + deci;
							edInf.gradient = dstRoi[j + i*roiCols];
							edInf.angle = angAll[j + i*roiCols];
							edgeInfor.push_back(edInf);
							k++;
						}
					}
					else
					{
						if ((angAll[j + i*roiCols] > 202.5) && (angAll[j + i*roiCols] < 247.5))
						{
							if ((dstRoi[j + i*roiCols] > dstRoi[j - 1 + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j + 1) + (i + 1)*roiCols]))
							{
								k3 = dstRoi[j - 1 + (i - 1)*roiCols];
								k2 = dstRoi[j + i*roiCols];
								k1 = dstRoi[(j + 1) + (i + 1)*roiCols];
								deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

								edInf.xyInteger.x = j;
								edInf.xyInteger.y = i;
								edInf.xyDecimal.x = j - deci;
								edInf.xyDecimal.y = i - deci;
								edInf.gradient = dstRoi[j + i*roiCols];
								edInf.angle = angAll[j + i*roiCols];
								edgeInfor.push_back(edInf);
								k++;
							}
						}
						else
						{
							if ((angAll[j + i*roiCols] > 112.5) && (angAll[j + i*roiCols] < 157.5))
							{

								if ((dstRoi[j + i*roiCols] > dstRoi[(j + 1) + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j - 1) + (i + 1)*roiCols]))
								{
									k1 = dstRoi[(j + 1) + (i - 1)*roiCols];
									k2 = dstRoi[j + i*roiCols];
									k3 = dstRoi[(j - 1) + (i + 1)*roiCols];
									deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

									edInf.xyInteger.x = j;
									edInf.xyInteger.y = i;
									edInf.xyDecimal.x = j - deci;
									edInf.xyDecimal.y = i + deci;
									edInf.gradient = dstRoi[j + i*roiCols];
									edInf.angle = angAll[j + i*roiCols];
									edgeInfor.push_back(edInf);
									k++;
								}
							}
							else
							{
								if ((angAll[j + i*roiCols] > 292.5) && (angAll[j + i*roiCols] < 337.5))
								{
									if ((dstRoi[j + i*roiCols] > dstRoi[(j + 1) + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j - 1) + (i + 1)*roiCols]))
									{
										k3 = dstRoi[(j + 1) + (i - 1)*roiCols];
										k2 = dstRoi[j + i*roiCols];
										k1 = dstRoi[(j - 1) + (i + 1)*roiCols];
										deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

										edInf.xyInteger.x = j;
										edInf.xyInteger.y = i;
										edInf.xyDecimal.x = j + deci;
										edInf.xyDecimal.y = i - deci;
										edInf.gradient = dstRoi[j + i*roiCols];
										edInf.angle = angAll[j + i*roiCols];
										edgeInfor.push_back(edInf);
										k++;
									}
								}
								else
								{
									if (((angAll[j + i*roiCols] >= -1) && (angAll[j + i*roiCols] <= 22.5)) || ((angAll[j + i*roiCols] >= 337.5) && (angAll[j + i*roiCols] <= 361)))
									{
										if ((dstRoi[j + i*roiCols] > dstRoi[(j - 1) + i*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j + 1) + i*roiCols]))
										{
											k1 = dstRoi[(j - 1) + i*roiCols];
											k2 = dstRoi[j + i*roiCols];
											k3 = dstRoi[(j + 1) + i*roiCols];
											deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

											edInf.xyInteger.x = j;
											edInf.xyInteger.y = i;
											edInf.xyDecimal.x = j + deci;
											edInf.xyDecimal.y = i;
											edInf.gradient = dstRoi[j + i*roiCols];
											edInf.angle = angAll[j + i*roiCols];
											edgeInfor.push_back(edInf);
											k++;
										}
									}
									else
									{
										if ((angAll[j + i*roiCols] <= 202.5) && (angAll[j + i*roiCols] >= 157.5))
										{
											if ((dstRoi[j + i*roiCols] > dstRoi[(j - 1) + i*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j + 1) + i*roiCols]))
											{
												k3 = dstRoi[(j - 1) + i*roiCols];
												k2 = dstRoi[j + i*roiCols];
												k1 = dstRoi[(j + 1) + i*roiCols];
												deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

												edInf.xyInteger.x = j;
												edInf.xyInteger.y = i;
												edInf.xyDecimal.x = j - deci;
												edInf.xyDecimal.y = i;
												edInf.gradient = dstRoi[j + i*roiCols];
												edInf.angle = angAll[j + i*roiCols];
												edgeInfor.push_back(edInf);
												k++;
											}
										}
										else
										{
											if ((angAll[j + i*roiCols] >= 67.5) && (angAll[j + i*roiCols] <= 112.5))
											{

												if ((dstRoi[j + i*roiCols] > dstRoi[j + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[j + (i + 1)*roiCols]))
												{
													k1 = dstRoi[j + (i - 1)*roiCols];
													k2 = dstRoi[j + i*roiCols];
													k3 = dstRoi[j + (i + 1)*roiCols];
													deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

													edInf.xyInteger.x = j;
													edInf.xyInteger.y = i;
													edInf.xyDecimal.x = j;
													edInf.xyDecimal.y = i + deci;
													edInf.gradient = dstRoi[j + i*roiCols];
													edInf.angle = angAll[j + i*roiCols];
													edgeInfor.push_back(edInf);
													k++;
												}
											}
											else
											{
												if ((angAll[j + i*roiCols] >= 247.5) && (angAll[j + i*roiCols] <= 292.5))
												{
													if ((dstRoi[j + i*roiCols] > dstRoi[j + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[j + (i + 1)*roiCols]))
													{
														k3 = dstRoi[j + (i - 1)*roiCols];
														k2 = dstRoi[j + i*roiCols];
														k1 = dstRoi[j + (i + 1)*roiCols];
														deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

														edInf.xyInteger.x = j;
														edInf.xyInteger.y = i;
														edInf.xyDecimal.x = j;
														edInf.xyDecimal.y = i - deci;
														edInf.gradient = dstRoi[j + i*roiCols];
														edInf.angle = angAll[j + i*roiCols];
														edgeInfor.push_back(edInf);
														k++;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}


		eNum = k;
		//二值图
		for (int t = 0; t < roiCols*roiRows; t++)//二值图像，所有像素先都赋值为0，边缘点赋值255
		{
			dstRoiE[t] = 0;
		}
		for (int q = 0; q < k; q++)
		{
			dstRoiE[edgeInfor[q].xyInteger.x + edgeInfor[q].xyInteger.y * roiCols] = 255;
		}
		/*
		FILE *Binary;
		Binary = fopen("E:\\ProjectMy\\Binary.txt", "w");
		for (int i = 0; i<roiRows; i++)
		{
		for (int j = 0; j < roiCols; j++)
		{
		fprintf(Binary, "%d   ", dstRoiE[j + i*roiCols]);
		}
		fprintf(Binary, "\n");
		}
		fclose(Binary);
		*/
		//以数组的方式传出边缘信息
		edgeArray = (edgeInformation*)malloc(k * sizeof(edgeInformation));
		for (int i = 0; i < k; i++)
		{
			edgeArray[i] = edgeInfor[i];
		}


		/*FILE *e;
		e = fopen("E:\\project03\\e.txt", "w");
		for (int i = 0; i < k; i++)
		{
		fprintf(e,"%d   %d   \n", edgeInfor[i].xInteger, edgeInfor[i].yInteger);
		}
		*/

		//free(angAll);
		return 1;
	}


}



CVisHoughTransform::CVisHoughTransform()
{
	downLevel = 2;
}

CVisHoughTransform::~CVisHoughTransform()
{
}

int CVisHoughTransform::detectCircle(IMG_UBBUF srcBuf, vector<houghCircle3f> &bestCircles)
{
	int status = 0;
	VisBuf setVisbuf;

	//get edge points
	int threshold = 0;
	IMG_WORD *dstRoi = new IMG_WORD[srcBuf.size.height*srcBuf.size.width];
	IMG_UBYTE *dstRoiE = new IMG_UBYTE[srcBuf.size.height*srcBuf.size.width];
	Ipp32f *angAll = new Ipp32f[srcBuf.size.height*srcBuf.size.width];
	edgeInformation *edgeArray = NULL;
	IMG_INT eNum = 0;

	status = ommTool::VisEdge_detection(srcBuf.ptr, srcBuf.size, threshold, dstRoi, dstRoiE, angAll,edgeArray, eNum);

	//show edgeDetect result
	IMG_WBUF mag_wBuf;
	setVisbuf.set_IMG_WBUF(mag_wBuf, dstRoi, srcBuf.size, srcBuf.linestep * sizeof(IMG_WORD));
	IMG_RBUF angle_rBuf;
	setVisbuf.set_IMG_RBUF(angle_rBuf, angAll, srcBuf.size, srcBuf.linestep * sizeof(IMG_REAL));

	//detect circle
	int voteScore = 5;
	int minRadius = 40;
	int maxRadius = 100;
	int cen_Distance = 140;		//shoule be minR + maxR
	
	houghCircle(edgeArray, eNum,angle_rBuf, voteScore, minRadius, maxRadius, cen_Distance, bestCircles);

	//free
	delete[] dstRoi;
	delete[] dstRoiE;
	delete[] angAll;
	if (edgeArray != NULL)
	{
		free(edgeArray);
	}


	return 0;
}

int regionDFS(int *pic,int *label,int r, int c, int &cnt)
{

	return 0;
}

//////////////////////////////////////////////////
//houghCircle功能说明：霍夫圆检测
//Input
//	edgeInformation *edgeArray  边缘点信息
//	IMG_INT eNum  边缘点个数
//  IMG_RBUF angle_rBuf 梯度方向图
//  int voteScore 投票阈值
//  int minRadius 最小半径
//  int maxRadius 最大半径 
//  int center_Dis 圆心距离

//Output
//	vector<houghCircle3f> &bestCircles  输出所有最佳圆

//////////////////////////////////////////////////
int CVisHoughTransform::houghCircle(edgeInformation *edgeArray, IMG_INT eNum, /*IMG_WBUF mag_wBuf,*/IMG_RBUF angle_rBuf, int voteScore, int minRadius, int maxRadius, int center_Dis, vector<houghCircle3f> &bestCircles)
{
	VisBuf setVisbuf;

	int HEIGHT = angle_rBuf.size.height;
	int WIDTH = angle_rBuf.size.width;
	int DEPTH = maxRadius;

	int ***H;

	// Allocate memory
	H = new int**[HEIGHT + 1];
	for (int i = 0; i < HEIGHT; i++) 
	{
		H[i] = new int*[WIDTH + 1];
		for (int j = 0; j < WIDTH; j++)
		{
			H[i][j] = new int[DEPTH + 1];
			memset(H[i][j], 0, sizeof(int)* (DEPTH + 1) );
		}
	}

	//accumulate1
	for (int k = 0; k < eNum; k++)
	{
		for (int r = minRadius; r <= maxRadius; r++)
		{
			int temp_x = edgeArray[k].xyInteger.x;
			int temp_y = edgeArray[k].xyInteger.y;

			int x0 = (int)round(edgeArray[k].xyDecimal.x - r * cos(angle_rBuf.ptr[temp_y * WIDTH + temp_x] / 180 * PI));
			int y0 = (int)round(edgeArray[k].xyDecimal.y - r * sin(angle_rBuf.ptr[temp_y * WIDTH + temp_x] / 180 * PI));
			inc_if_inside(H, x0, y0, HEIGHT, WIDTH, r);
		}
	}

	//accumulate2
	//for (int y = 0; y < HEIGHT; y++)
	//{
	//	for (int x = 0; x < WIDTH; x++)
	//	{
	//		// printf("data point : %f\n", img_data.at<float>(y,x));
	//		if ((float)mag_wBuf.ptr[y * WIDTH + x] > 10000.0)  //threshold image  
	//		{
	//			for (int r = minRadius; r <= maxRadius; r++)
	//			{
	//				int x0 = round(x - r * cos(angle_rBuf.ptr[y * WIDTH + x] / 180 * PI) );
	//				int y0 = round(y - r * sin(angle_rBuf.ptr[y * WIDTH + x] / 180 * PI) );
	//				inc_if_inside(H, x0, y0, HEIGHT, WIDTH, r);
	//			}
	//		}
	//	}
	//}

	//show hAcc;
	IMG_WBUF hAcc_wBuf;
	IMG_WORD *pHAccData = new IMG_WORD[WIDTH * HEIGHT];
	memset(pHAccData, 0, sizeof(IMG_WORD) * WIDTH * HEIGHT);
	int *radius_ref = new int[HEIGHT * WIDTH];
	memset(radius_ref, 0, sizeof(int) * WIDTH * HEIGHT);
	int *haccLabel = new int[HEIGHT * WIDTH];
	memset(haccLabel, 0, sizeof(int) * HEIGHT * WIDTH);

	//int pos = 0;
	//for (int y0 = 0; y0 < HEIGHT; y0++) 
	//{
	//	for (int x0 = 0; x0 < WIDTH; x0++) 
	//	{
	//		for (int r = minRadius; r <= maxRadius; r++)
	//		{
	//			if (H[y0][x0][r] > pHAccData[pos + x0] )
	//			{
	//				pHAccData[pos + x0] = H[y0][x0][r];
	//				radius_ref[pos + x0] = r;
	//			}
	//		}
	//	}
	//	pos += WIDTH;
	//}

	//label region


	//setVisbuf.set_IMG_WBUF(hAcc_wBuf, pHAccData, { (IMG_UWORD)WIDTH,(IMG_UWORD)HEIGHT }, WIDTH * sizeof(IMG_WORD));


	//compute optimal circles
	for (int y0 = 0; y0 < HEIGHT; y0++) 
	{
		for (int x0 = 0; x0 < WIDTH; x0++) 
		{
			for (int r = minRadius; r <= maxRadius; r++) 
			{
				if (H[y0][x0][r] > voteScore)
				{
					houghCircle3f circle;
					circle.centerX = (float)x0;
					circle.centerY = (float)y0;
					circle.radius = (float)r;
					int i;
					for (i = 0; i < (int)bestCircles.size(); i++) {
						int xCoord = (int)bestCircles[i].centerX;
						int yCoord = (int)bestCircles[i].centerY;
						int radius = (int)bestCircles[i].radius;
						if (abs(xCoord - x0) < center_Dis && abs(yCoord - y0) < center_Dis)		//圆心距在一定范围（同圆）
						{
							if (H[y0][x0][r] > H[yCoord][xCoord][radius])	//票数高，替换
							{
								bestCircles.erase(bestCircles.begin() + i);
								bestCircles.insert(bestCircles.begin(), circle);
							}
							break;
						}
					}
					if (i == bestCircles.size()) {
						bestCircles.insert(bestCircles.begin(), circle);
					}
				}
			}
		}
	}

	//free
	for (int i = 0; i < HEIGHT; i++)
	{
		for (int j = 0; j < WIDTH; j++)
		{
			if (H[i][j] != NULL)
			{
				delete[] H[i][j];
			}

		}
	}
	for (int i = 0; i < HEIGHT; i++)
	{
		if (H[i] != NULL)
			delete[] H[i];
	}
	if(H != NULL)
		delete[] H;

	delete[] pHAccData;
	delete[] radius_ref;
	delete[] haccLabel;
	return 0;
}

int CVisHoughTransform::inc_if_inside(int *** H, int x, int y, int height, int width, int r)
{
	if (x >= 0 && x < width && y >= 0 && y < height)
	{
		H[y][x][r]++;
	}

	return 0;
}

IppStatus CVisHoughTransform::pyramid(IMG_UBBUF src, unsigned char * pDst, int & pyramid_width, int & pyramid_height, int level)
{
	IppStatus   status = ippStsNoErr;
	//ofstream outfile("pyramidData.txt");

	IppiPyramid*  pPyrStruct = NULL;
	unsigned char **  pPyrImage = NULL;

	//init var
	//double sigma = 3;
	unsigned char* pSrc = src.ptr;
	IMG_SIZE roiSize = src.size;

	if (!pSrc) { status = ippStsNoMemErr; goto exit; }

	int srcStep = roiSize.width * sizeof(unsigned char);
	float      rate = 2.f;                  // Neighbour levels ratio
	signed short kernel[3] = { 1,1,1 };		// Separable symmetric kernel of odd length

											//signed short *kernel = (signed short *)malloc(3 * sizeof(signed short));
											//__GetGaussianKernel_dim1(kernel, 3, sigma);	// preserved

	int pyrBufferSize = 0;
	int pyrStructSize = 0;
	unsigned char   *pPyrBuffer = NULL;
	unsigned char   *pPyrStrBuffer = NULL;

	int      pyrLStateSize = 0;
	int      pyrLBufferSize = 0;
	unsigned char   *pPyrLStateBuf = NULL;
	unsigned char   *pPyrLBuffer = NULL;

	// Computes the temporary work buffer size
	status = ippiPyramidGetSize(&pyrStructSize, &pyrBufferSize, level, { roiSize.width,roiSize.height }, rate);
	if (status < 0) goto exit;

	//pPyrBuffer = ippsMalloc_8u(pyrBufferSize);
	//pPyrStrBuffer = ippsMalloc_8u(pyrStructSize);
	pPyrBuffer = (unsigned char*)malloc(pyrBufferSize * sizeof(unsigned char));
	pPyrStrBuffer = (unsigned char*)malloc(pyrStructSize * sizeof(unsigned char));	//not pop
	if ((pyrBufferSize && !pPyrBuffer) || (pyrStructSize && !pPyrStrBuffer)) {
		status = ippStsNoMemErr; goto exit_proc;
	}

	// Initializes Gaussian structure for pyramids
	//pPyrStruct = (IppiPyramid*)malloc(level * sizeof(IppiPyramid));	
	status = ippiPyramidInit(&pPyrStruct, level, { roiSize.width,roiSize.height }, rate, pPyrStrBuffer, pPyrBuffer);
	if (status < 0) goto exit_proc;

	// ????????????????Correct maximum scale level 
	level = pPyrStruct->level;

	// Allocate structures to calculate pyramid layers
	status = ippiPyramidLayerDownGetSize_8u_C1R({ roiSize.width,roiSize.height }, rate, 3, &pyrLStateSize, &pyrLBufferSize);
	if (status < 0) goto exit_proc;

	//pPyrLStateBuf = ippsMalloc_8u(pyrLStateSize);
	//pPyrLBuffer = ippsMalloc_8u(pyrLBufferSize);
	pPyrLStateBuf = (unsigned char*)malloc(pyrLStateSize * sizeof(unsigned char));
	pPyrLBuffer = (unsigned char*)malloc(pyrLBufferSize * sizeof(unsigned char));
	if ((pyrLStateSize && !pPyrLStateBuf) || (pyrLBufferSize && !pPyrLBuffer)) { status = ippStsNoMemErr; goto exit; }

	// Initialize the structure for creating a lower pyramid layer
	status = ippiPyramidLayerDownInit_8u_C1R((IppiPyramidDownState_8u_C1R**)&pPyrStruct->pState, { roiSize.width,roiSize.height }, rate, kernel, 3, IPPI_INTER_LINEAR, pPyrLStateBuf, pPyrLBuffer);
	if (status < 0) goto exit_proc;

	// Allocate pyramid layers
	pPyrImage = pPyrStruct->pImage;
	pPyrImage[0] = pSrc;
	pPyrStruct->pStep[0] = srcStep;

	for (int i = 1; i <= level; i++)
	{
		//pPyrImage[i] = ippiMalloc_8u_C1(pPyrStruct->pRoi[i].width, pPyrStruct->pRoi[i].height, &pPyrStruct->pStep[i]);
		pPyrImage[i] = (unsigned char*)malloc((pPyrStruct->pRoi[i].width) * (pPyrStruct->pRoi[i].height) * sizeof(unsigned char));
		pPyrStruct->pStep[i] = (pPyrStruct->pRoi[i].width) * sizeof(unsigned char);
		if (!pPyrImage[i]) { status = ippStsNoMemErr; goto exit_proc; }
	}

	// Perform downsampling of the image with 5x5 Gaussian kernel
	for (int i = 1; i <= level; i++)
	{
		status = ippiPyramidLayerDown_8u_C1R(pPyrImage[i - 1], pPyrStruct->pStep[i - 1], pPyrStruct->pRoi[i - 1], pPyrImage[i], pPyrStruct->pStep[i], pPyrStruct->pRoi[i], (IppiPyramidDownState_8u_C1R*)pPyrStruct->pState);
		if (status < 0) goto exit_proc;

	}

	for (int i = 0; i < pPyrStruct->pRoi[level].height; i++)
	{
		for (int j = 0; j < pPyrStruct->pRoi[level].width; j++)
		{
			pDst[i * pPyrStruct->pRoi[level].width + j] = pPyrImage[level][i * pPyrStruct->pRoi[level].width + j];
		}
	}
	pyramid_width = pPyrStruct->pRoi[level].width;
	pyramid_height = pPyrStruct->pRoi[level].height;

	//	test //
	//for (int i = 0; i < pPyrStruct->pRoi[level].height; i++)
	//{
	//	for (int j = 0; j < pPyrStruct->pRoi[level].width; j++)
	//	{
	//		outfile << (float)pPyrImage[level][i * pPyrStruct->pRoi[level].width + j] << " ";
	//	}
	//	outfile << endl;
	//}
	//outfile.close();



exit_proc:
	for (int i = 1; i <= level; i++)
		free(pPyrImage[i]);
	free(pPyrStrBuffer);
	free(pPyrBuffer);
	free(pPyrLBuffer);
	free(pPyrLStateBuf);
	//for (int i = 1; i <= level; i++)
	//	ippiFree(pPyrImage[i]);
	//ippiFree(pPyrLStateBuf);
	//ippiFree(pPyrLBuffer);
	//ippiFree(pPyrStrBuffer);
	//ippiFree(pPyrBuffer);
	goto exit;

exit:
	//printf("pyramidStatus: %s\n", ippGetStatusString(status));
	if (!status)
	{
	}
	return status;
}




int CVisHoughTransform::newDetectCircle(IMG_UBBUF ubbSrc)
{
	int status = 0;
	VisBuf setVisbuf;
	int srcHeight = ubbSrc.size.height;
	int srcWidth = ubbSrc.size.width;

	//////////////////	pyramid layerDown	///////////////
	int downWidth = 0;
	int downHeight = 0;
	IMG_UBYTE *pPyramidData = new IMG_UBYTE[srcHeight * srcWidth];
	status = (int)pyramid(ubbSrc, pPyramidData, downWidth, downHeight, downLevel);
	//show pyramid result
	IMG_UBBUF ubbPyramid;
	setVisbuf.set_IMG_UBBUF(ubbPyramid, pPyramidData, { (IMG_UWORD)downWidth,(IMG_UWORD)downHeight }, downWidth);

	/////////////	sobel 5*5	///////////////////////
	int threshold = 0;
	IMG_WORD *dstRoi = new IMG_WORD[downWidth * downHeight];
	IMG_UBYTE *dstRoiE = new IMG_UBYTE[downWidth * downHeight];
	Ipp32f *angAll = new Ipp32f[downWidth * downHeight];
	edgeInformation *edgeArray = NULL;
	IMG_INT eNum = 0;
	status = ommTool::VisEdge_detection(ubbPyramid.ptr, ubbPyramid.size, threshold, dstRoi, dstRoiE, angAll, edgeArray, eNum);
	//show edgeDetect result
	IMG_WBUF wbGradMag;
	setVisbuf.set_IMG_WBUF(wbGradMag, dstRoi, ubbPyramid.size, ubbPyramid.linestep * sizeof(IMG_WORD));
	IMG_RBUF rbGradAngle;
	setVisbuf.set_IMG_RBUF(rbGradAngle, angAll, ubbPyramid.size, ubbPyramid.linestep * sizeof(IMG_REAL));

	///////////////////		seperate gradAngle		//////////////////


	//free
	delete[] pPyramidData;
	delete[] dstRoi;
	delete[] dstRoiE;
	delete[] angAll;
	if (edgeArray != NULL)
	{
		free(edgeArray);
	}


	return 0;
}

//int CVisHoughTransform::houghShape()