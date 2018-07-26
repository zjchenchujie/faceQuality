#include <random>
#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <algorithm>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include "laplacian.hpp"
// #include "hello.hpp"

// using namespace std;
// int Laplacian(const unsigned char* IN src_, unsigned char* OUT dst_, int IN width_, int IN height_, int IN ksize_, float &var)
using namespace cv;

int main()
{
	float result{0};
	float result2{0};
	// sayhello();
	// float* result_ptr = &result;
	// unsigned char* dst = NULL;
	cv::Mat src = cv::imread("/home/chujie/Documents/x1.jpg", 0);
	if (!src.data || src.channels() != 1) 
	{
		fprintf(stderr, "read image fail\n");
		return -1;
	}
	
	int width = src.cols;
	int height = src.rows;
	std::unique_ptr<unsigned char[]> dst(new unsigned char[width * height]);
	std::unique_ptr<unsigned char[]> dst2(new unsigned char[width * height]);

	// Mat imgSrc_gray, dst;

	// int kernel_size = 3;
	// int scale = 1;
	// int delta = 0;
	// int ddepth = CV_16S;


	// // cvtColor(src, imgSrc_gray, CV_BGR2GRAY);
	//  // Mat abs_dst;
	// Laplacian(src.data, dst, ddepth, kernel_size, scale, delta, BORDER_DEFAULT);

	// Mat mean, stddev;
	// meanStdDev(dst, mean, stddev);

	// double result = stddev.at<double>(0,0);

	Laplacian(src.data, dst2.get(), width, height, 1);
	Laplacian_(src.data, dst.get(), width, height, 1);

	// Mat mean, stddev;
	// meanStdDev(dst, mean, stddev);

	// result = stddev.at<double>(0,0);
	result = variance(dst.get(), width, height);
	result2 = variance(dst2.get(), width, height);
	// display(dst.get(), width, height);

	std::cout<<width<<std::endl;
	std::cout<<height<<std::endl;
	std::cout<<"Lap: "<<result<<std::endl;
	std::cout<<"Lap2: "<<result2<<std::endl;
	// std::cout<<(int)(dst[10])<<std::endl;

	return 0;
}
