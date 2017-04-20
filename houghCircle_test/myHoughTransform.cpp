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
	/*       ������: vsprintf                                       ������������
	/*       �� ��: �͸�ʽ�����������                                         ����
	/*       ����ֵ: ��������·��������ִ��ĳ���(��ȥ\0),����������ظ�ֵ
	/*       �� ��: int vsprintf(char *string, char *format, va_list param);
	/*                ��param ����ʽformatд���ַ���string��
	/*       ע: �ú���������ڴ�������,����ʹ��vsnprintf                ���� ��  ��                                                     */
	/***************************************************************************/
	va_end(st);
	cout << buf;// << endl;
#endif

	return 0;
}

bool points_comp(const edgeInformation &a, const edgeInformation &b)
{
	return a.gradient > b.gradient;
}

//2����Ե���

///////////////////////////////////////////////////////////////////////////
//VisEdge_detection����˵��
//Input
//srcRoi   ����ͼ��
//roiSize  ����ͼ��ĳߴ�
//threshold  �ݶ�ǿ�ȵ���ֵ���������Ǵ���0С��1250��������������û���֪����ֵ��Ϊ���ٺ��ʣ���������0���㷨�Զ���ȡ���ʵ���ֵ����
//
//output
//dstRoi  �ݶ�ǿ��
//edgeInformation *&edgeArray  ��Ե����Ϣ�������������ꡢ���������ꡢ�ݶ�ǿ�ȡ��Ƕ�
//
//��������
//��������·���1��
//����û�������ֵС��0����������-1��
//��������������ȷ������ͼ��ߴ粻��ȷ��srcRoi��С��ߴ粻���ϣ���������-1��
//Author������/20170227
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

	//����һ��ֵthreshold�����������䷽��
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
		//�����ֵС��0������ֱ�ӷ���-1
		if (threshold < 0)
		{
			return -1;
		}

		int roiRows = roiSize.height;
		int roiCols = roiSize.width;

		//�Ƕ���Ϣ
//		Ipp32f *angAll;
		//angAll = (Ipp32f*)malloc(roiRows*roiCols * sizeof(Ipp32f));

		std::vector<edgeInformation> edgeInfor;
		edgeInformation edInf;

		int k = 0;//��¼��Ե��ĸ���
		Ipp16u k1;//��������ϵ�������֪��
		Ipp16u k2;
		Ipp16u k3;
		float deci;//��������϶����С�����֣�����Ӧ��������
		float sumx = 0;//��Ե���x����֮��
		float sumy = 0;
		int numberChannels = 1; //the source image is single channel

		IppiSize dstRoiSize = { roiCols,roiRows };

		SobelFilter_8u16s_C1_5x5(srcRoi, dstRoiSize, dstRoi, angAll);


		//�ѽǶ���[-PI��PI]��Ϊ[0��360]
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
		//�Զ���ȡ�ݶ�ǿ�ȵ���ֵ
		//ʲô���������û��������ֵ�أ�
		if (threshold == 0)
		{
			//Otsu�����������еĻҶ�ֵ����1��255��ʹintra-class invariance�����Ǹ�ֵ����ΪҪ�����ֵ
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

		//��������
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
		//��ֵͼ
		for (int t = 0; t < roiCols*roiRows; t++)//��ֵͼ�����������ȶ���ֵΪ0����Ե�㸳ֵ255
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
		//������ķ�ʽ������Ե��Ϣ
		edgeArray = (edgeInformation*)malloc(k * sizeof(edgeInformation));
		for (int i = 0; i < k; i++)
		{
			edgeArray[i] = edgeInfor[i];
		}

		//FILE *e;
		//e = fopen("e.txt", "w");
		//for (int i = 0; i < k; i++)
		//{
		//	fprintf(e,"%d   %d   \n", edgeArray[i].xyInteger.x, edgeArray[i].xyInteger.y);
		//}
		//fclose(e);

		//free(angAll);
		return 1;
	}


}



CVisHoughTransform::CVisHoughTransform()
{
	m_downLevel = 2;
	
	m_sectors = 36;	//10 du
	m_selectedRatio = (float)0.3;
	m_nSelectMin = 200 / pow(2, m_downLevel);		//��ʵӦ�ø��ܵı�Ե��������
	m_nSelectMax = 400 / pow(2, m_downLevel);

	m_Ttheta = 5;
	m_Tshift = 10;

	m_localThreshMin = 20;	//��ͨ������ֵ��̫С-��ͨ��̫�̫࣬��-����һ����ͨ�򱻷ָ��ˡ�

	m_radiusMax = 120 / pow(2,m_downLevel);
	m_radiusMin = 40 / pow(2, m_downLevel);

	m_voteScoreMin = m_nSelectMin;
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

void CVisHoughTransform::regionDFS(IMG_UWORD *pic, IMG_UWORD *label,int r, int c, int height,int width,int id, IMG_COORD * storeMax,int threshold)
{
	if (r < 0 || r >= height || c < 0 || c >= width)		return;
	if (pic[r * width + c] <= threshold || label[r * width + c] != 0)	return;

	label[r * width + c] = id;
	if (pic[r * width + c] > pic[storeMax[id].y * width + (storeMax[id].x)])
	{
		storeMax[id].y = r;
		storeMax[id].x = c;
	}
	for (int dr = -1; dr <= 1; dr++)
	{
		for (int dc = -1; dc <= 1; dc++)
		{
			if (dr != 0 || dc != 0)
			{
				regionDFS(pic, label, r + dr, c + dc, height, width, id, storeMax,threshold);
			}
		}
	}

}

//////////////////////////////////////////////////
//houghCircle����˵��������Բ���
//Input
//	edgeInformation *edgeArray  ��Ե����Ϣ
//	IMG_INT eNum  ��Ե�����
//  IMG_RBUF angle_rBuf �ݶȷ���ͼ
//  int voteScore ͶƱ��ֵ
//  int minRadius ��С�뾶
//  int maxRadius ���뾶 
//  int center_Dis Բ�ľ���

//Output
//	vector<houghCircle3f> &bestCircles  ����������Բ

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
						if (abs(xCoord - x0) < center_Dis && abs(yCoord - y0) < center_Dis)		//Բ�ľ���һ����Χ��ͬԲ��
						{
							if (H[y0][x0][r] > H[yCoord][xCoord][radius])	//Ʊ���ߣ��滻
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

void CVisHoughTransform::getGaussianKernel_dim2(IMG_LREAL **gaus, const int size, const double sigma)
{
	const double PIII = 4.0*atan(1.0); //Բ���ʦи�ֵ  
	IMG_INT center = size / 2;
	IMG_LREAL sum = 0;
	for (IMG_INT i = 0; i < size; i++)
	{
		for (IMG_INT j = 0; j < size; j++)
		{
			gaus[i][j] = (1 / (2 * PIII*sigma*sigma))*exp(-((i - center)*(i - center) + (j - center)*(j - center)) / (2 * sigma*sigma));
			sum += gaus[i][j];
		}
	}
	for (IMG_INT i = 0; i < size; i++)
	{
		for (IMG_INT j = 0; j < size; j++)
		{
			gaus[i][j] /= sum;
		}
	}

}

int CVisHoughTransform::gaussfilter(IMG_UBBUF src, IMG_UBYTE *pDst,int kernalSize,double sigma)
{
	int status = ippStsNoErr;
	//init var
	IMG_UBYTE *pSrc = src.ptr;
	IppiSize roiSize = {src.size.width,src.size.height };
	IMG_INT size = kernalSize;

	if (pSrc == NULL)
	{
		status = ippStsNoMemErr;
		goto exit;
	}

	//begin
	{
		IMG_INT srcStep = sizeof(IMG_UBYTE) * roiSize.width;
		IMG_INT dstStep = srcStep;

		IMG_LREAL **gaus = new IMG_LREAL *[size];
		for (int i = 0; i < size; i++)
		{
			gaus[i] = new IMG_LREAL[size];
		}

		if (gaus == NULL)
		{
			status = ippStsNoMemErr;
			goto exit;
		}

		getGaussianKernel_dim2(gaus, size, sigma);
		IppiSize  kernelSize = { size,size };
		IppiSize gaussize = { size,size };
		//IMG_WORD kernel[9];
		IMG_WORD *kernel = new IMG_WORD[size * size];
		for (IMG_INT i = 0; i < size; i++)
		{
			for (IMG_INT j = 0; j < size; j++)
			{
				kernel[i*size + j] = gaus[i][j] * 32767;		//-32767 ~ 32767
			}
		}
		/*for (IMG_INT i = 0; i<size; i++)
		{
		for (IMG_INT j = 0; j<size; j++)
		{
		printf("%d,", kernel[i*size + j]);
		}
		}*/
		IMG_INT divisor = 32767;//��һ��
		IMG_UBYTE *pBuffer = NULL;                /* Pointer to the work buffer */
		IppiFilterBorderSpec* pSpec = NULL;   /* context structure */
		IMG_INT iTmpBufSize = 0, iSpecSize = 0;   /* Common work buffer size */
		IppiBorderType borderType = ippBorderConst;
		IMG_UBYTE borderValue = 0;
		IMG_INT numChannels = 1;

		status = ippiFilterBorderGetSize(kernelSize, roiSize, ipp8u, ipp16s, numChannels, &iSpecSize, &iTmpBufSize);
		if (status < 0)	goto exit_proc;

		pSpec = (IppiFilterBorderSpec *)malloc(iSpecSize);
		pBuffer = new IMG_UBYTE[iTmpBufSize];//(IMG_UBYTE *)pool.PMalloc(iTmpBufSize);
		if (!pSpec || !pBuffer)
		{
			status = ippStsNoMemErr;
			goto exit_proc;
		}

		status = ippiFilterBorderInit_16s(kernel, kernelSize, divisor, ipp8u, numChannels, ippRndNear, pSpec);
		if (status < 0)	goto exit_proc;

		status = ippiFilterBorder_8u_C1R(pSrc, srcStep, pDst, dstStep, roiSize, borderType, &borderValue, pSpec, pBuffer);
		if (status < 0)	goto exit_proc;

	exit_proc:
		delete[] pBuffer;
		free(pSpec);
		delete[] kernel;
		for (int i = 0; i < size; i++)
		{
			delete[] gaus[i];
		}
		delete[] gaus;
		goto exit;
	}
exit:
	//printf("gaussianFilterStatus: %s\n", ippGetStatusString(status));
	return status;
}

int CVisHoughTransform::gaussfilter_UWORD(IMG_UWBUF src, IMG_UWORD *pDst, int kernalSize, double sigma)
{
	int status = ippStsNoErr;
	//init var
	IMG_UWORD *pSrc = src.ptr;
	IppiSize roiSize = { src.size.width,src.size.height };
	IMG_INT size = kernalSize;

	if (pSrc == NULL)
	{
		status = ippStsNoMemErr;
		goto exit;
	}

	//begin
	{
		IMG_INT srcStep = sizeof(IMG_UWORD) * roiSize.width;
		IMG_INT dstStep = srcStep;

		IMG_LREAL **gaus = new IMG_LREAL *[size];
		for (int i = 0; i < size; i++)
		{
			gaus[i] = new IMG_LREAL[size];
		}

		if (gaus == NULL)
		{
			status = ippStsNoMemErr;
			goto exit;
		}

		getGaussianKernel_dim2(gaus, size, sigma);
		IppiSize  kernelSize = { size,size };
		IppiSize gaussize = { size,size };
		//IMG_WORD kernel[9];
		IMG_WORD *kernel = new IMG_WORD[size * size];
		for (IMG_INT i = 0; i < size; i++)
		{
			for (IMG_INT j = 0; j < size; j++)
			{
				kernel[i*size + j] = gaus[i][j] * 32767;		//-32767 ~ 32767
			}
		}
		/*for (IMG_INT i = 0; i<size; i++)
		{
		for (IMG_INT j = 0; j<size; j++)
		{
		printf("%d,", kernel[i*size + j]);
		}
		}*/
		IMG_INT divisor = 32767;//��һ��
		IMG_UBYTE *pBuffer = NULL;                /* Pointer to the work buffer */
		IppiFilterBorderSpec* pSpec = NULL;   /* context structure */
		IMG_INT iTmpBufSize = 0, iSpecSize = 0;   /* Common work buffer size */
		IppiBorderType borderType = ippBorderConst;
		IMG_UWORD borderValue = 0;
		IMG_INT numChannels = 1;

		status = ippiFilterBorderGetSize(kernelSize, roiSize, ipp16u, ipp16s, numChannels, &iSpecSize, &iTmpBufSize);
		if (status < 0)	goto exit_proc;

		pSpec = (IppiFilterBorderSpec *)malloc(iSpecSize);
		pBuffer = new IMG_UBYTE[iTmpBufSize];//(IMG_UBYTE *)pool.PMalloc(iTmpBufSize);
		if (!pSpec || !pBuffer)
		{
			status = ippStsNoMemErr;
			goto exit_proc;
		}

		status = ippiFilterBorderInit_16s(kernel, kernelSize, divisor, ipp16u, numChannels, ippRndNear, pSpec);
		if (status < 0)	goto exit_proc;

		status = ippiFilterBorder_16u_C1R(pSrc, srcStep, pDst, dstStep, roiSize, borderType, &borderValue, pSpec, pBuffer);
		if (status < 0)	goto exit_proc;

	exit_proc:
		delete[] pBuffer;
		free(pSpec);
		delete[] kernel;
		for (int i = 0; i < size; i++)
		{
			delete[] gaus[i];
		}
		delete[] gaus;
		goto exit;
	}
exit:
	//printf("gaussianFilterStatus: %s\n", ippGetStatusString(status));
	return status;
}

int CVisHoughTransform::findLocalmaximum(IMG_UWBUF uwbSrc)
{
	int status = ippStsNoErr;
	VisBuf setVisbuf;

	int width = uwbSrc.size.width;
	int height = uwbSrc.size.height;
	IMG_UWORD *pFilter = new IMG_UWORD[width * height];

	//do gaussianFilter to find accuracy center region
	gaussfilter_UWORD(uwbSrc, pFilter, 3, 3);
	IMG_UWBUF uwbFilter;
	setVisbuf.set_IMG_UWBUF(uwbFilter, pFilter, uwbSrc.size, uwbSrc.linestep);
	

	//select
	//for (int k = 0; k < width * height; k++)
	//{
	//	if (pMarker[k] < m_localThreshMin)
	//	{
	//		pMarker[k] = 0;
	//	}
	//}
	//int markerStep = width * sizeof(Ipp16u);
	//int minLabel = 1;                      /* Minimal label value */
	//int maxLabel = 65534;                  /* Maximal label value */
	//int markersNum;                        /* Pointer to number of markers */
	//Ipp8u* pBuffer = NULL;                 /* Pointer to the work buffer */
	//int bufferSize = 0;
	///* Calculate size of temporary buffer */
	//status = ippiLabelMarkersGetBufferSize_16u_C1R({ width,height }, &bufferSize);
	//if (status < 0) goto exit;
	//pBuffer = ippsMalloc_8u(bufferSize);
	//if (pBuffer == NULL) { status = ippStsNoMemErr; goto exit; }
	///* Label connected non-zero components with different label values */
	//status = ippiLabelMarkers_16u_C1IR(pMarker, markerStep, { width,height }, minLabel, maxLabel, ippiNormInf, &markersNum, pBuffer);		//should be ippiNormInf
	//if (status < 0) goto exit;
	//IMG_UWBUF uwbMarker;
	//setVisbuf.set_IMG_UWBUF(uwbMarker, pMarker, uwbSrc.size, uwbSrc.linestep);

	//find local maximum
	IMG_COORD * storeMax = new IMG_COORD[width * height];
	for (int k = 0; k < width * height; k++)
	{
		storeMax[k].x = 0;
		storeMax[k].y = 0;
	}
	IMG_UWORD *pLabel = new IMG_UWORD[width * height];
	memset(pLabel, 0, width * height * sizeof(IMG_UWORD));
	
	int pos = 0;
	int cnt = 0;		//marker nums
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			if (pLabel[pos + j] == 0 && pFilter[pos + j] > m_localThreshMin)
			{
				regionDFS(pFilter, pLabel, i, j, height, width, ++cnt, storeMax,m_localThreshMin);
			}
		}
		pos += width;
	}
	IMG_UWBUF uwbLabel;
	setVisbuf.set_IMG_UWBUF(uwbLabel, pLabel, uwbSrc.size, uwbSrc.linestep);

	//show center points
	IMG_UWBUF uwbAllCenter;
	IMG_UWORD *pAllCenter = new IMG_UWORD[width * height];
	memset(pAllCenter, 0, width * height * sizeof(IMG_UWORD));
	circleCenter.clear();
	for (int k = 1; k <= cnt; k++)
	{
		IMG_WORD y = storeMax[k].y;
		IMG_WORD x = storeMax[k].x;
		pAllCenter[y * width + x] = 100;
		circleCenter.push_back({ x,y });
	}
	setVisbuf.set_IMG_UWBUF(uwbAllCenter, pAllCenter, uwbSrc.size, uwbSrc.linestep);

exit:
	//if (pBuffer)
	//{
	//	ippsFree(pBuffer);
	//}
	delete[] pFilter;
	delete[] pLabel;
	delete[] storeMax;
	delete[] pAllCenter;
	return status;
}

int CVisHoughTransform::newDetectCircle(IMG_UBBUF ubbSrc)
{
	int status = 0;
	VisBuf setVisbuf;
	int srcHeight = ubbSrc.size.height;
	int srcWidth = ubbSrc.size.width;

	//////////////////	pyramid layerDown (8ms)	///////////////
	int downWidth = 0;
	int downHeight = 0;
	IMG_UBYTE *pPyramidData = new IMG_UBYTE[srcHeight * srcWidth];
	if (m_downLevel)
	{
		status = (int)pyramid(ubbSrc, pPyramidData, downWidth, downHeight, m_downLevel);
	}
	else
	{
		memcpy(pPyramidData, ubbSrc.ptr, srcHeight * srcWidth);
		downHeight = srcHeight;
		downWidth = srcWidth;
	}
	//show pyramid result
	IMG_UBBUF ubbPyramid;
	setVisbuf.set_IMG_UBBUF(ubbPyramid, pPyramidData, { (IMG_UWORD)downWidth,(IMG_UWORD)downHeight }, downWidth);

	///////////////	sobel 5*5 (40+ ms)	///////////////////////
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

	///////////////////		seperate gradAngle	(ms)	//////////////////
	float range = (float)360.0 / m_sectors;
	float inv_range = (float)1.0 / range;
	angleTable.clear();
	angleTable.resize(m_sectors);
	for (int k = 0; k < eNum; k++)
	{
		int _angle = int(edgeArray[k].angle + range / 2) % 360;
		int angleIndex = (int)(_angle * inv_range);
		if (angleIndex == m_sectors)
		{
			angleIndex--;
		}
		angleTable[angleIndex].push_back(edgeArray[k]);
	}
	//sort
	for (int k = 0; k < m_sectors; k++)
	{
		sort(angleTable[k].begin(), angleTable[k].end(), points_comp);
		int newSize = (int)(angleTable[k].size() * m_selectedRatio);
		if (newSize < m_nSelectMin )
		{
			newSize = m_nSelectMin;
		}
		else if(newSize > m_nSelectMax)
		{
			newSize = m_nSelectMax;
		}
		angleTable[k].resize(newSize);
	}
//#ifdef DEBUG
//	//FILE *fTable;
//	//fTable = fopen("fTable.txt", "w");
//	//for (int i = 0; i < angleTable.size(); i++)
//	//{
//	//	for (int j = 0; j < angleTable[i].size(); j++)
//	//	{
//	//		fprintf(fTable, "%d   %d   \n", angleTable[i][j].xyInteger.x, angleTable[i][j].xyInteger.y);
//	//	}
//	//}
//	//fclose(fTable);
//#endif // DEBUG
//
//	
	//////////////////		find opposite points pair and accumulate center (50 ms) ///////////////////////////////
	IMG_UWBUF hAcc2_wBuf;
	IMG_UWORD *pHacc2 = new IMG_UWORD[downHeight * downWidth];
	memset(pHacc2, 0, sizeof(IMG_UWORD) * downWidth * downHeight);
	//int t = 0;
	for (int i = 0; i < angleTable.size();i++)
	{
		for (int j = 0; j < angleTable[i].size(); j++)
		{
			int opposite_sec = (i + m_sectors / 2) % m_sectors;
			for (int k = 0; k < angleTable[opposite_sec].size(); k++)
			{
				float lineTheta = atan2((angleTable[i][j].xyDecimal.y - angleTable[opposite_sec][k].xyDecimal.y), (angleTable[i][j].xyDecimal.x - angleTable[opposite_sec][k].xyDecimal.x)) / PI * 180;
				if (lineTheta < 0)
				{
					lineTheta += 360;
				}
				if (fabs(lineTheta - angleTable[i][j].angle) > m_Tshift && fabs(lineTheta - angleTable[opposite_sec][k].angle) > m_Tshift)
				{
					//cout << ++t << endl;
					continue;
				}
				if ( fabs( (int)(angleTable[i][j].angle + 180) % 360 - (int)angleTable[opposite_sec][k].angle )  > m_Ttheta)
				{
					continue;
				}
				//else , calculate center of points pair
				IMG_WORD centerX = round( (angleTable[i][j].xyDecimal.x + angleTable[opposite_sec][k].xyDecimal.x) / 2);
				IMG_WORD centerY = round( (angleTable[i][j].xyDecimal.y + angleTable[opposite_sec][k].xyDecimal.y) / 2);
				//accumulate
				pHacc2[centerY * downWidth + centerX] += 1;
			}
		}
	}
	//show pHacc2
	setVisbuf.set_IMG_UWBUF(hAcc2_wBuf, pHacc2, ubbPyramid.size, downWidth * sizeof(IMG_UWORD));
	
	///////////////////////////		find local max		//////////////////////////////////
	findLocalmaximum(hAcc2_wBuf);

	//////////////////////////		accumulate r		//////////////////////////////////
	IMG_UWORD **pHacc3 = new IMG_UWORD*[downHeight * downWidth];
	for (int i = 0; i <= downHeight * downWidth - 1; i++)
	{
		pHacc3[i] = new IMG_UWORD[m_radiusMax + 1];
		memset(pHacc3[i], 0, sizeof(IMG_UWORD) * (m_radiusMax + 1));
	}
	//cycle all center points candidate
	for (int k = 0; k < circleCenter.size(); k++)
	{
		for (int i = 0; i < angleTable.size(); i++)
		{
			for (int j = 0; j < angleTable[i].size(); j++)
			{
				int _distance = sqrt(pow((circleCenter[k].x - angleTable[i][j].xyDecimal.x), 2) + pow((circleCenter[k].y - angleTable[i][j].xyDecimal.y), 2));
				if (_distance >= m_radiusMin && _distance <= m_radiusMax)
				{
					pHacc3[circleCenter[k].y * downWidth + circleCenter[k].x][_distance]++;
				}
			}
		}
	}
	//show pHacc3 max r 
	IMG_UWBUF hAcc3_wBuf;
	IMG_UWORD *pHacc3_show = new IMG_UWORD[downHeight * downWidth];
	memset(pHacc3_show, 0, sizeof(IMG_UWORD) * downHeight * downWidth);

	bestCircles.clear();
	int pos = 0;
	for (int i = 0; i < downHeight; i++)
	{
		for (int j = 0; j < downWidth; j++)
		{
			int temp_r = 0;
			for (int r = m_radiusMin; r <= m_radiusMax; r++)
			{
				if (pHacc3[pos + j][r] > pHacc3_show[pos + j])
				{
					pHacc3_show[pos + j] = pHacc3[pos + j][r];
					temp_r = r;
				}
			}
			//push final center points using m_voteScoreMin
			if (pHacc3_show[pos + j] > m_voteScoreMin)
			{
				houghCircle3i _cicle;
				_cicle.center = { (IMG_WORD)j,(IMG_WORD)(pos / downWidth) };
				_cicle.radius = temp_r;
				bestCircles.push_back(_cicle);
			}
		}
		pos += downWidth;
	}
	setVisbuf.set_IMG_UWBUF(hAcc3_wBuf, pHacc3_show, ubbPyramid.size, downWidth * sizeof(IMG_UWORD));


	//free
	for (int i = 0; i <= downHeight * downWidth - 1; i++)
	{
		//cout << i << endl;
 		delete[] pHacc3[i];
	}

	delete[] pHacc3;
	delete[] pHacc3_show;
	delete[] pHacc2;
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

vector<houghCircle3i> CVisHoughTransform::getBestCircles() 
{ 
	for (int i = 0; i < bestCircles.size(); i++)
	{
		bestCircles[i].center.x *= pow(2, m_downLevel);
		bestCircles[i].center.y *= pow(2, m_downLevel);
		bestCircles[i].radius *= pow(2, m_downLevel);
	}
	return bestCircles; 
}

void CVisHoughTransform::setParams(int downLevel, 
								int sectors, 
								float selectedRatio, 
								int selectMin, 
								int selectMax, 
								float Ttheta, 
								float Tshift, 
								int localThreshMin,		
								int radiusMin, 
								int radiusMax, 
								int voteScoreMin)
{
	m_downLevel = downLevel;

	m_sectors = sectors;
	m_selectedRatio = selectedRatio;
	m_nSelectMin = selectMin;
	m_nSelectMax = selectMax;

	m_Ttheta = Ttheta;
	m_Tshift = Tshift;

	m_localThreshMin = localThreshMin;

	m_radiusMin = radiusMin / pow(2, m_downLevel);
	m_radiusMax = radiusMax / pow(2, m_downLevel);

	m_voteScoreMin = voteScoreMin;
}
