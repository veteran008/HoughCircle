#include "stdafx.h"
#include "houghCircle_test.h"
#include "myHoughTransform.h"


using namespace cv;

//int test1()
//{
//
//	char* imageName = "img/coins1.png";
//
//	Mat image, img_grey, img_grey1;     //input mat
//	Mat dx, dy, mag, dist;
//	Mat dx_out, dy_out, dis_out;  //final output mat
//	Mat h_acc, h_out;       //hough space matricies
//
//	image = imread(imageName, 1);
//	cvtColor(image, img_grey, COLOR_BGR2GRAY);
//
//	dx.create(img_grey.rows, img_grey.cols, CV_32FC1);
//	dy.create(img_grey.rows, img_grey.cols, CV_32FC1);
//	mag.create(img_grey.rows, img_grey.cols, CV_32FC1);
//	dist.create(img_grey.rows, img_grey.cols, CV_32FC1);
//
//	sobel(img_grey, dx, dy, mag, dist);
//
//	//normalize arrays with max and min values of 255 and 0
//	normalize(dx, dx_out, 0, 255, NORM_MINMAX, -1, Mat());
//	normalize(dy, dy_out, 0, 255, NORM_MINMAX, -1, Mat());
//	normalize(dist, dis_out, 0, 255, NORM_MINMAX, -1, Mat());
//
//	// double H_h = sqrt(2.0) * (double) (mag.rows>mag.cols ? mag.rows : mag.cols); //-r -> +r
//	// double H_w = 180;  
//
//	h_acc.create(mag.rows, mag.cols, CV_32FC1);
//	memset(h_acc.data, 0, sizeof(float) * mag.rows * mag.cols);
//	// threshold(h_acc,h_out,0,255,THRESH_TOZERO);
//	// threshold(h_out,h_acc,0,255,THRESH_TOZERO);
//
//	// medianBlur(mag,mag,1);
//	//mag,dist,thresh,minRad,maxRad,Dist-circles,output_Hspace, final_result
//	hough(mag, dist, 5, 15, 100, 45, h_acc, image);
//
//	// normalize(h_acc, h_out, 0, 255, NORM_MINMAX, -1, Mat());
//	// threshold(h_acc,h_out, 200,255,THRESH_TOZERO);
//
//	//save images
//	imwrite("dx.jpg", dx_out);
//	imwrite("dy.jpg", dy_out);
//	imwrite("mag.jpg", mag);
//	imwrite("dist.jpg", dis_out);
//	imwrite("h_space.jpg", h_acc);
//
//	imwrite("result.png", image);
//	return 0;
//}

int OpenCV_HoughCircle(int pic_num)
{
	printf("\n------kase:%d ------------\n", pic_num);
	char readBmpPath[100];
	char saveBmpPath[100];
	sprintf(readBmpPath, "D:/PICs/Hough_Circle/testCase/t (%u).bmp", pic_num);
	sprintf(saveBmpPath, "D:/PICs/Hough_Circle/testCase/Result/CV_res(%u).bmp", pic_num);

	Mat image, img_grey;     //input mat
	image = imread(readBmpPath, 1);
	if (image.empty())
	{
		printf("no image\n");
		return -1;
	}
	cvtColor(image, img_grey, COLOR_BGR2GRAY);
	GaussianBlur(img_grey, img_grey, Size(9, 9), 2, 2);
	int circleNum_Max = 20;
	vector<Vec3f> circles(circleNum_Max);		//init!!!!!!!!!!!!!!!!!
	////////////////////////////////////////////////////////////////////////////////////////
	double t = (double)cvGetTickCount();
	HoughCircles(img_grey, circles, CV_HOUGH_GRADIENT,1, img_grey.rows / 8, 200, 30,60,110);
	//HoughCircles(img_grey, circles, CV_HOUGH_GRADIENT, 1, 80, 200, 30, 60, 110);
	t = (double)cvGetTickCount() - t;
	printf("OpenCV houghCircle exec time = %g ms\n", t / (cvGetTickFrequency() * 1000));
	///////////////////////////////////////////////////////////////////////////////////////

	for (size_t i = 0; i < circles.size(); i++)
	{
		Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
		int radius = cvRound(circles[i][2]);
		// draw the circle center
		circle(image, center, 3, Scalar(0, 255, 0), -1, 8, 0);
		// draw the circle outline
		circle(image, center, radius, Scalar(0, 0, 255), 5, 8, 0);
	}

	imwrite(saveBmpPath, image);
	
	return 0;
}

int test_detectCircle()
{
	VisBuf setVisBf;
	CVisHoughCircle houghObject;
	char* imageName = "img/2017-2-28-22-52-31.bmp";

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

int test_newDetectCircle(int pic_num)
{
	printf("\n------kase:%d ------------\n", pic_num);
	char readBmpPath[100];
	char saveBmpPath[100];
	sprintf(readBmpPath, "D:/PICs/Hough_Circle/testCase/t (%u).bmp", pic_num);
	sprintf(saveBmpPath, "D:/PICs/Hough_Circle/testCase/Result/BM_res(%u).bmp", pic_num);

	VisBuf setVisBf;
	CVisHoughCircle houghObject;
	//char* imageName = "img/2017-2-28-23-54-10.bmp";
	Mat image, img_grey;     //input mat
	image = imread(readBmpPath, 1);
	if (image.data == NULL)
	{
		printf("no image\n");
		return -1;
	}
	cvtColor(image, img_grey, COLOR_BGR2GRAY);

	IMG_UBBUF ubbSrc;
	setVisBf.set_IMG_UBBUF(ubbSrc, img_grey.data, { (IMG_UWORD)img_grey.cols,(IMG_UWORD)img_grey.rows }, (IMG_UWORD)img_grey.cols);

	//////////////////////////////////////////////////////////////////////////////////////////////
	double t = (double)cvGetTickCount();
	houghObject.newDetectCircle(ubbSrc);			
	t = (double)cvGetTickCount() - t;
	printf("my houghCircle exec time = %g ms\n", t / (cvGetTickFrequency() * 1000));
	//////////////////////////////////////////////////////////////////////////////////////////////
	vector<houghCircle3i> bestCircles = houghObject.getBestCircles();

	//draw
#ifdef DEBUG
	int WIDTH = image.cols;
	int HEIGHT = image.rows;
	for (int i = 0; i < (int)bestCircles.size(); i++)
	{
		int lineThickness = 5;
		int lineType = 10;
		int shift = 0;
		int xCoord = (int)bestCircles[i].center.x;
		int yCoord = (int)bestCircles[i].center.y;
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
		circle(image, center, radius - 2, Scalar(255, 0, 0), lineThickness, lineType, shift);
	}
#endif

	imwrite(saveBmpPath, image);

	return 0;
}

int test_gaussianFilter()
{
	VisBuf setVisBf;
	CVisHoughCircle houghObject;
	char* imageName = "img/280.000-280.000-188.800.bmp";
	Mat image, img_grey;     //input mat
	image = imread(imageName, 1);
	cvtColor(image, img_grey, COLOR_BGR2GRAY);

	IMG_UBBUF ubbSrc;
	setVisBf.set_IMG_UBBUF(ubbSrc, img_grey.data, { (IMG_UWORD)img_grey.cols,(IMG_UWORD)img_grey.rows }, (IMG_UWORD)img_grey.cols);

	IMG_UBYTE *pDst = new IMG_UBYTE[img_grey.cols * img_grey.rows];
	//houghObject.gaussfilter(ubbSrc, pDst, 3, 3);

	img_grey.data = pDst;
	imwrite("img/gaussianTest.bmp",img_grey);

	delete[] pDst;
	return 0;
}




int main()
{
	
	//freopen("timeLog.txt", "w", stdout);

	//test1();
	//test_detectCircle();
	for (int i = 1; i <= 15; i++)
	{
		//test_newDetectCircle(i);
	}
	for (int i = 1; i <= 15; i++)
	{
		OpenCV_HoughCircle(i);
	}
	//test_gaussianFilter();

	system("pause");
	_CrtDumpMemoryLeaks();
	return 0;
}