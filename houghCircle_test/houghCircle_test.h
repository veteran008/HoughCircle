#pragma once

#include "stdafx.h"
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <math.h>

//#include <highgui.h>
//#include <opencv2/opencv.hpp>
//#include <opencv2/core/core.hpp> 

using namespace cv;

void sobel(Mat img, Mat &dx, Mat &dy, Mat &mag, Mat &dist);

void inc_if_inside(double *** H, int x, int y, int height, int width, int r);

void hough(Mat &img_data, Mat &dist, double threshold, int minRadius, int maxRadius, double distance, Mat &h_acc, Mat &coins);

