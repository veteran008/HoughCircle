﻿// houghCircle_test.cpp : 定义控制台应用程序的入口点。
//
#include "stdafx.h"
#include "houghCircle_test.h"



void sobel(Mat img, Mat &dx, Mat &dy, Mat &mag, Mat &dist)
{
	float acc_dx = 0, acc_dy = 0;         //accumulators
	float k1[] = { -1,-2,-1,0,0,0,1,2,1 }; //{-2,-4,-2,0,0,0,2,4,2};//{-1,-2,-1,0,0,0,1,2,1};    //sobel kernal dx
	float k2[] = { -1,0,1,-2,0,2,-1,0,1 };//{-2,0,2,-4,0,4,-2,0,2};//{-1,0,1,-2,0,2,-1,0,1};    //sobel kernal dy

	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			acc_dx = acc_dy = 0;

			//apply kernel/mask
			for (int nn = -1; nn < 2; nn++) {
				for (int mm = -1; mm < 2; mm++) {
					if (i + nn > 0 && i + nn < img.rows && j + mm > 0 && j + mm < img.cols) {
						acc_dx += (float)img.at<uchar>(i + nn, j + mm) * k1[((mm + 1) * 3) + nn + 1];
						acc_dy += (float)img.at<uchar>(i + nn, j + mm) * k2[((mm + 1) * 3) + nn + 1];
					}
				}
			}
			//write final values
			dx.at<float>(i, j) = acc_dx;
			dy.at<float>(i, j) = acc_dy;
			mag.at<float>(i, j) = (sqrtf(acc_dy*acc_dy + acc_dx*acc_dx)) > 100 ? 255 : 0;
			dist.at<float>(i, j) = atan2f(acc_dy, acc_dx);
			// printf("dist : %f \n", dist.at<float>(i,j) / 3.14159265f * 180 );
		}
	}
}



void inc_if_inside(double *** H, int x, int y, int height, int width, int r)
{
	if (x > 0 && x < width && y> 0 && y < height)
		H[y][x][r]++;
}


void hough(Mat &img_data, Mat &dist, double threshold, int minRadius, int maxRadius, double distance, Mat &h_acc, Mat &coins) 
{
	//int radiusRange = maxRadius - minRadius;
	int HEIGHT = img_data.rows;
	int WIDTH = img_data.cols;
	int DEPTH = maxRadius;

	double ***H;

	// Allocate memory
	H = new double**[HEIGHT];
	for (int i = 0; i < HEIGHT; ++i) {
		H[i] = new double*[WIDTH];

		for (int j = 0; j < WIDTH; ++j)
		{
			H[i][j] = new double[DEPTH];
			memset(H[i][j], 0, sizeof(double)*DEPTH);
		}
	}

	for (int y = 0; y < img_data.rows; y++)
	{
		for (int x = 0; x < img_data.cols; x++)
		{
			// printf("data point : %f\n", img_data.at<float>(y,x));
			if ((float)img_data.at<float>(y, x) > 250.0)  //threshold image  
			{
				for (int r = minRadius; r <= maxRadius; r++)
				{

					int x0 = round(x + r * cos(dist.at<float>(y, x)));		//-180 ~ 180
					int x1 = round(x - r * cos(dist.at<float>(y, x)));
					int y0 = round(y + r * sin(dist.at<float>(y, x)));
					int y1 = round(y - r * sin(dist.at<float>(y, x)));


					inc_if_inside(H, x0, y0, HEIGHT, WIDTH, r);		
					// inc_if_inside(H,x0,y1,HEIGHT, WIDTH, r);
					// inc_if_inside(H,x1,y0,HEIGHT, WIDTH, r);
					inc_if_inside(H, x1, y1, HEIGHT, WIDTH, r);			//?????????
				}
			}
		}
	}

	//create 2D image by summing values of the radius dimension
	for (int y0 = 0; y0 < HEIGHT; y0++) {
		for (int x0 = 0; x0 < WIDTH; x0++) {
			for (int r = minRadius; r <= maxRadius; r++) 
			{
				//zhan
				if (H[y0][x0][r] > h_acc.at<float>(y0, x0))
				{
					h_acc.at<float>(y0, x0) = H[y0][x0][r];
				}
				//h_acc.at<float>(y0, x0) += H[y0][x0][r];// > 1 ? 255 : 0;
										// printf("h : %d", H[y0][x0][r]);
			}
		}
	}

	std::vector<Point3f> bestCircles;

	//compute optimal circles
	for (int y0 = 0; y0 < HEIGHT; y0++) {
		for (int x0 = 0; x0 < WIDTH; x0++) {
			for (int r = minRadius; r <= maxRadius; r++) {
				if (H[y0][x0][r] > threshold) 
				{
					Point3f circle(x0, y0, r);
					int i;
					for (i = 0; i < bestCircles.size(); i++) {
						int xCoord = bestCircles[i].x;
						int yCoord = bestCircles[i].y;
						int radius = bestCircles[i].z;
						if (abs(xCoord - x0) < distance && abs(yCoord - y0) < distance)		//圆心距在一定范围（同圆）
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

	//	draw	//
	for (int i = 0; i < bestCircles.size(); i++) {
		int lineThickness = 1;
		int lineType = 10;
		int shift = 0;
		int xCoord = bestCircles[i].x;
		int yCoord = bestCircles[i].y;
		int radius = bestCircles[i].z;
		Point2f center(xCoord, yCoord);
		for (int dy = -2; dy <= 2; dy++)
		{
			for (int dx = 0; dx <= 0; dx++)
			{
				if ((xCoord + dx) >= 0 && (xCoord + dx) < WIDTH && (yCoord + dy) >= 0 && (yCoord + dy) < HEIGHT)
				{
					coins.at<Vec3b>(yCoord + dy, xCoord + dx)[0] = 0;
					coins.at<Vec3b>(yCoord + dy, xCoord + dx)[1] = 0;
					coins.at<Vec3b>(yCoord + dy, xCoord + dx)[2] = 255;
				}
			}
		}
		for (int dy = 0; dy <= 0; dy++)
		{
			for (int dx = -2; dx <= 2; dx++)
			{
				if ((xCoord + dx) >= 0 && (xCoord + dx) < WIDTH && (yCoord + dy) >= 0 && (yCoord + dy) < HEIGHT)
				{
					coins.at<Vec3b>(yCoord + dy, xCoord + dx)[0] = 0;
					coins.at<Vec3b>(yCoord + dy, xCoord + dx)[1] = 0;
					coins.at<Vec3b>(yCoord + dy, xCoord + dx)[2] = 255;
				}
			}
		}
		circle(coins, center, radius - 1, Scalar(255, 0, 0), lineThickness, lineType, shift);
	}

	//for (int i = 0; i < HEIGHT; i++)
	//{
	//	for (int j = 0; j < WIDTH; j++)
	//	{
	//		if (H[i][j] != NULL)
	//		{
	//			delete[] H[i][j];
	//			//cout << i << "," << j << endl;
	//		}

	//	}
	//}
}



