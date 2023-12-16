#include <iostream>

#include <opencv2/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
//#include <aruco/src/aruco.h>
#include <opencv2/aruco.hpp>
#include <stdio.h>

using namespace cv;
using namespace std;

//弧度制转角度制
double radian_to_angle(double rad)
{
	double PI = 3.1415926;
	return (rad*180)/PI;
	
}

//3D坐标：这个坐标指的是ARUCO的世界坐标系。
vector<cv::Point3f> Generate3DPoints()
{
	vector<cv::Point3f> points;

	float x,y,z;
	x = .5; y= 0 ;  z = 0;
	points.push_back(cv::Point3f(x,y,z));
	x = 0;  y = .5; z = 0;
	points.push_back(cv::Point3f(x,y,z));
	x = 0; y = 0;   z = .5;
	points.push_back(cv::Point3f(x,y,z));
	return points;
}


int main(int argc, char *argv[])

{
	int Marker_x,Marker_y,err_x,err_y; //Marker是ARUCO的中心坐标，err是ARUCO与图片中心的偏差值。
	float k,angle; //k是ARUCO与图片的斜率，angle是ARUCO与图片的角度
	cv::Mat image, image_copy;
	float K[3][3] = {3.06272049496904e+03, 0., 1.00674972631e+03, 0.,3.05644872162e+03, 1.71771441719e+03, 0., 0., 1.}; //标定参数：内参
	float J[1][5] = { 1.3110582e-01, -2.5083437e-01,7.5872e-04, -8.127e-04,4.4984408e-01}; //标定参数：外参
	cv::Mat camera_matrix = cv::Mat(3,3,CV_32FC1,K);
	cv::Mat dist_coeffs = cv::Mat(1,5,CV_32FC1,J);
	vector<cv::Point3f> objectPoints = Generate3DPoints();
	cout << "a = " << endl << camera_matrix << endl << endl;
	cout << "b = " << endl << dist_coeffs << endl << endl;

	//设置识别ARUCO的类别
	cv::Ptr<cv::aruco::Dictionary> dictionary =
	cv::aruco::getPredefinedDictionary( \
		cv::aruco::PREDEFINED_DICTIONARY_NAME(10));
	
	//图片的长宽
	int cols = 640;
	int rows = 480;

	cv::VideoCapture inputVideo;
	inputVideo.open(1);
	while (inputVideo.grab())
	{
		inputVideo.retrieve(image);//抓取视频中的一张照片
		image.copyTo(image_copy);

        std::vector<int> ids; //识别到ARUCO的数量
        std::vector<std::vector<cv::Point2f> > corners; //识别到的ARUCO在图片中的坐标
		cv::aruco::detectMarkers(image, dictionary, corners, ids); //识别图片中的ARUCO
		if (ids.size() > 0)
		{
			std::vector<cv::Vec3d> rvecs, tvecs;
			cv::aruco::estimatePoseSingleMarkers(corners, 1,camera_matrix, dist_coeffs, rvecs, tvecs);//用ARUCO在图片中的坐标对应世界坐标系，来计算旋转矩阵和平移向量。
			
			//计算ARUCO的图片坐标系的中心
			Marker_x = ((corners[0][0].x+corners[0][1].x)/2+(corners[0][2].x+corners[0][3].x)/2)/2;
			Marker_y = ((corners[0][0].y+corners[0][1].y)/2+(corners[0][2].y+corners[0][3].y)/2)/2;

			//计算ARUCOO的图片坐标系中心与相机中心的差值，单位是像素
			err_x = Marker_x - cols/2;
			err_y = Marker_y - rows/2;

			//再用旋转向量和平移向量求出ARUCO坐标系下的点在图像中的坐标
			std::vector<cv::Point2f> projectedPoints; 
			cv::projectPoints(objectPoints,rvecs,tvecs,camera_matrix,dist_coeffs,projectedPoints);
			cv::line(image,cv::Point(projectedPoints[1].x,projectedPoints[1].y),cv::Point(Marker_x,Marker_y),(0,0,255),5);
			k = (-(projectedPoints[1].y - Marker_y)) / (projectedPoints[1].x - Marker_x);
			angle = radian_to_angle(atan(k));
			if((projectedPoints[1].y > Marker_y) && (projectedPoints[1].x > Marker_x))
			{
				angle = 360 + angle;		
			}
			if((projectedPoints[1].y > Marker_y) && (projectedPoints[1].x < Marker_x))
			{
				angle = 180 + angle;
			}
			if((projectedPoints[1].y < Marker_y) && (projectedPoints[1].x < Marker_x))
			{
				angle = 180 + angle;
			}
			cout << "------------------>angle:    "<< angle << endl;
		}
		cv::line(image,cv::Point(330,240),cv::Point(310,240),(0,255,0),2);
		cv::line(image,cv::Point(320,250),cv::Point(320,230),(0,255,0),2);
		cv::imshow("out", image);
		cv::waitKey(50);
	}

	return 0;
}


