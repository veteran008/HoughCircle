#include "stdafx.h"
#include "houghCircle_test.h"
#include "myHoughTransform.h"

using namespace cv;

int test1()
{

	char* imageName = "img/coins1.png";

	Mat image, img_grey, img_grey1;     //input mat
	Mat dx, dy, mag, dist;
	Mat dx_out, dy_out, dis_out;  //final output mat
	Mat h_acc, h_out;       //hough space matricies

	image = imread(imageName, 1);
	cvtColor(image, img_grey, COLOR_BGR2GRAY);

	dx.create(img_grey.rows, img_grey.cols, CV_32FC1);
	dy.create(img_grey.rows, img_grey.cols, CV_32FC1);
	mag.create(img_grey.rows, img_grey.cols, CV_32FC1);
	dist.create(img_grey.rows, img_grey.cols, CV_32FC1);

	sobel(img_grey, dx, dy, mag, dist);

	//normalize arrays with max and min values of 255 and 0
	normalize(dx, dx_out, 0, 255, NORM_MINMAX, -1, Mat());
	normalize(dy, dy_out, 0, 255, NORM_MINMAX, -1, Mat());
	normalize(dist, dis_out, 0, 255, NORM_MINMAX, -1, Mat());

	// double H_h = sqrt(2.0) * (double) (mag.rows>mag.cols ? mag.rows : mag.cols); //-r -> +r
	// double H_w = 180;  

	h_acc.create(mag.rows, mag.cols, CV_32FC1);
	memset(h_acc.data, 0, sizeof(float) * mag.rows * mag.cols);
	// threshold(h_acc,h_out,0,255,THRESH_TOZERO);
	// threshold(h_out,h_acc,0,255,THRESH_TOZERO);

	// medianBlur(mag,mag,1);
	//mag,dist,thresh,minRad,maxRad,Dist-circles,output_Hspace, final_result
	hough(mag, dist, 5, 15, 100, 45, h_acc, image);

	// normalize(h_acc, h_out, 0, 255, NORM_MINMAX, -1, Mat());
	// threshold(h_acc,h_out, 200,255,THRESH_TOZERO);

	//save images
	imwrite("dx.jpg", dx_out);
	imwrite("dy.jpg", dy_out);
	imwrite("mag.jpg", mag);
	imwrite("dist.jpg", dis_out);
	imwrite("h_space.jpg", h_acc);

	imwrite("result.png", image);
	return 0;
}



int test_detectCircle()
{
	VisBuf setVisBf;
	CVisHoughTransform houghObject;
	char* imageName = "img/4model.png";

	Mat image, img_grey;     //input mat

	image = imread(imageName, 1);
	cvtColor(image, img_grey, COLOR_BGR2GRAY);
	//img_grey.convertTo(img_grey, CV_16SC1);

	IMG_UBBUF srcRoi;
	setVisBf.set_IMG_UBBUF(srcRoi, img_grey.data, { (IMG_UWORD)img_grey.cols,(IMG_UWORD)img_grey.rows }, (IMG_UWORD)img_grey.cols);
	
	//detect
	vector<houghCircle3f> bestCircles;
	houghObject.detectCircle(srcRoi, bestCircles);

	//draw
#ifdef DEBUG
	int WIDTH = image.cols;
	int HEIGHT = image.rows;
	for (int i = 0; i < (int)bestCircles.size(); i++) 
	{
		int lineThickness = 1;
		int lineType = 10;
		int shift = 0;
		int xCoord = (int)bestCircles[i].centerX;
		int yCoord = (int)bestCircles[i].centerY;
		int radius = (int)bestCircles[i].radius;
		Point2i center(xCoord, yCoord);
		for (int dy = -2; dy <= 2; dy++)
		{
			for (int dx = 0; dx <= 0; dx++)
			{
				if ((xCoord + dx) >= 0 && (xCoord + dx) < WIDTH && (yCoord + dy) >= 0 && (yCoord + dy) < HEIGHT)
				{
					image.at<Vec3b>(yCoord + dy, xCoord + dx)[0] = 0;
					image.at<Vec3b>(yCoord + dy, xCoord + dx)[1] = 0;
					image.at<Vec3b>(yCoord + dy, xCoord + dx)[2] = 255;
				}
			}
		}
		for (int dy = 0; dy <= 0; dy++)
		{
			for (int dx = -2; dx <= 2; dx++)
			{
				if ((xCoord + dx) >= 0 && (xCoord + dx) < WIDTH && (yCoord + dy) >= 0 && (yCoord + dy) < HEIGHT)
				{
					image.at<Vec3b>(yCoord + dy, xCoord + dx)[0] = 0;
					image.at<Vec3b>(yCoord + dy, xCoord + dx)[1] = 0;
					image.at<Vec3b>(yCoord + dy, xCoord + dx)[2] = 255;
				}
			}
		}
		circle(image, center, radius - 1, Scalar(255, 0, 0), lineThickness, lineType, shift);
	}
#endif
	return 0;
}

int test_newDetectCircle()
{
	VisBuf setVisBf;
	CVisHoughTransform houghObject;
	char* imageName = "img/2017-2-28-22-52-31.bmp";

	Mat image;     //input mat

	image = imread(imageName, 0);

	IMG_UBBUF ubbSrc;
	setVisBf.set_IMG_UBBUF(ubbSrc, image.data, { (IMG_UWORD)image.cols,(IMG_UWORD)image.rows }, (IMG_UWORD)image.cols);

	houghObject.newDetectCircle(ubbSrc);

	return 0;
}


int main()
{
	double t = (double)cvGetTickCount();
	//test1();
	//test_detectCircle();
	test_newDetectCircle();
	t = (double)cvGetTickCount() - t;
	printf("exec time = %g ms\n", t / (cvGetTickFrequency() * 1000));
	_CrtDumpMemoryLeaks();
	return 0;
}