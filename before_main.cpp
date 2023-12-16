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
	String calib_images_path = "/home/casia/Downloads/aruco/calib_images";
	printf(CV_VERSION);
	cv::Mat image, image_copy;
	float K[3][3] = {3.06272049496904e+03, 0., 1.00674972631e+03, 0.,3.05644872162e+03, 1.71771441719e+03, 0., 0., 1.};
	float J[1][5] = { 1.3110582e-01, -2.5083437e-01,7.5872e-04, -8.127e-04,4.4984408e-01};
	cv::Mat camera_matrix = cv::Mat(3,3,CV_32FC1,K);
	cv::Mat dist_coeffs = cv::Mat(1,5,CV_32FC1,J);
	vector<cv::Point3f> objectPoints = Generate3DPoints();
	cout << "a = " << endl << camera_matrix << endl << endl;
	cout << "b = " << endl << dist_coeffs << endl << endl;
	std::ostringstream vector_to_marker;
	//camera_matrix[0][0] = 3062.72049496;
	//dist_coeffs = [ 1.3110582e-01, -2.5083437e-01,7.5872e-04, -8.127e-04,4.4984408e-01 ];
	cv::Ptr<cv::aruco::Dictionary> dictionary =
	cv::aruco::getPredefinedDictionary( \
		cv::aruco::PREDEFINED_DICTIONARY_NAME(10));
	int cols = 640;
	int rows = 480;

	//cv::FileStorage fs("calibration_params.yml", cv::FileStorage::READ);

	//fs["camera_matrix"] >> camera_matrix;
	//fs["distortion_coefficients"] >> dist_coeffs;

	cv::VideoCapture inputVideo;
	inputVideo.open(1);
	while (inputVideo.grab())
	{
		inputVideo.retrieve(image);//抓取视频中的一张照片
		image.copyTo(image_copy);

        	std::vector<int> ids;
        	std::vector<std::vector<cv::Point2f> > corners;
		cv::aruco::detectMarkers(image, dictionary, corners, ids);
		if (ids.size() > 0)
		{
			cv::aruco::drawDetectedMarkers(image_copy, corners, ids);
			std::vector<cv::Vec3d> rvecs, tvecs;
            		cv::aruco::estimatePoseSingleMarkers(corners, 1,
                    		camera_matrix, dist_coeffs, rvecs, tvecs);
			cv::line(image_copy,cv::Point(330,240),cv::Point(310,240),(33,33,133),2);
			//for(int i=0; i<ids.size(); i++)
				//cv::aruco::drawAxis(image_copy, camera_matrix, dist_coeffs, rvecs[i], tvecs[i], 0.1);
			//cout << "width:" << image_copy.cols << endl;
			//cout << "height:" << image_copy.rows << endl;
			//cout << "corners[0]: " << corners[0] << endl;
			//cout << "corners[0][0]: " << corners[0][0] << endl;
			//cout << "corners[0][2]x: " << corners[0][2].x << endl;
			//cout << "corners[0][2]y: " << corners[0][2].y << endl;
			//cout << "corners[0][1]: " << corners[0][1].x << endl;
			//cout << "shape " << corners[0][0].x << endl;
			//cout << "shape corners " << sizeof(corners[8]) << endl;
			Marker_x = ((corners[0][0].x+corners[0][1].x)/2+(corners[0][2].x+corners[0][3].x)/2)/2;
			err_x = Marker_x - cols/2;
			//cout << "------------>centre_x:" << rows/2 << endl;
			//cout << "------->Marker_x:" << Marker_x << endl;
			//cout << "err_x:" << err_x << endl;	


			Marker_y = ((corners[0][0].y+corners[0][1].y)/2+(corners[0][2].y+corners[0][3].y)/2)/2;
			err_y = Marker_y - rows/2;
			//cout << "------------------------------------->centre_y:" << cols/2 << endl;
			//cout << "--------------------------->err_y:" << err_y << endl;	
			//cout << "---------------->Marker_y:" << Marker_y << endl;
			
			std::vector<cv::Point2f> projectedPoints;
			cv::projectPoints(objectPoints,rvecs,tvecs,camera_matrix,dist_coeffs,projectedPoints);
			cv::line(image_copy,cv::Point(projectedPoints[1].x,projectedPoints[1].y),cv::Point(Marker_x,Marker_y),(0,0,255),5);
			k = (-(projectedPoints[1].y - Marker_y)) / (projectedPoints[1].x - Marker_x);
			angle = radian_to_angle(atan(k));
			cout << k << endl;
			
			cout << "------->p.y:     " << projectedPoints[1].y << endl;
			cout << "------------->marker_y:     " << Marker_y << endl;
			if((projectedPoints[1].y > Marker_y) && (projectedPoints[1].x > Marker_x))
			{
				cout << "HH" << endl;
				angle = 360 + angle;		
			}
			if((projectedPoints[1].y > Marker_y) && (projectedPoints[1].x < Marker_x))
			{
				cout << "HEHE" << endl;
				angle = 180 + angle;
			}
			if((projectedPoints[1].y < Marker_y) && (projectedPoints[1].x < Marker_x))
			{
				cout << "HAHA" << endl;
				angle = 180 + angle;
			}
			cout << "------------------>angle:    "<< angle << endl;
			//Marker_y = ((corners[0][0][0].y+corners[0][0][1][1])/2+(corners[0][0][3][1]+corners[0][0][2][1])/2)/2;
			//cout << "x:"<<Marker_x << endl;
			//cout << "y:"<<Marker_y <<endl;
		}
		cv::line(image_copy,cv::Point(330,240),cv::Point(310,240),(0,255,0),2);
		cv::line(image_copy,cv::Point(320,250),cv::Point(320,230),(0,255,0),2);
		cv::imshow("out", image_copy);
		cv::waitKey(50);
	}

	return 0;
/*
//内参与畸变矩阵，笔者在前面的博客已经给出求解方法，有需要的可以找找看看
    double fx,fy,cx,cy,k1,k2,k3,p1,p2;
    fx=955.8925;
    fy=955.4439;
    cx=296.9006;
    cy=215.9074;
    k1=-0.1523;
    k2=0.7722;
    k3=0;
    p1=0;
    p2=0;

  Mat cameraMatrix = (cv::Mat_<float>(3, 3) <<
        fx, 0.0, cx,
        0.0, fy, cy,
        0.0, 0.0, 1.0);
   Mat distCoeffs = (cv::Mat_<float>(5, 1) << k1, k2, p1, p2, k3);
   cv::VideoCapture inputVideo;
   inputVideo.open(0);
   //cv::Ptr<cv::aruco::Dictionary> dictionary = cv::aruco::getPredefinedDictionary(cv::aruco::DICT_6X6_250);
   cv::Ptr<aruco::Dictionary> dictionary = aruco::getPredefinedDictionary(aruco::DICT_6X6_250);
//= aruco::getPredefinedDictionary(aruco::DICT_6X6_250);

   while (inputVideo.grab()) {
       cv::Mat image, imageCopy;
       inputVideo.retrieve(image);//抓取视频中的一张照片
       image.copyTo(imageCopy);
       std::vector<int> ids;
       std::vector<std::vector<cv::Point2f>> corners;
       cv::aruco::detectMarkers(image, dictionary, corners, ids);//检测靶标
       // if at least one marker detected
       if (ids.size() > 0) {
           cv::aruco::drawDetectedMarkers(imageCopy, corners, ids);//绘制检测到的靶标的框
           std::vector<cv::Vec3d> rvecs, tvecs;
           cv::aruco::estimatePoseSingleMarkers(corners, 0.055, cameraMatrix, distCoeffs, rvecs, tvecs);//求解旋转矩阵rvecs和平移矩阵tvecs
           //cout<<"R :"<<rvecs[0]<<endl;
           cout<<"T :"<<tvecs[0]<<endl;
           // draw axis for each marker
           for(int i=0; i<ids.size(); i++)
               cv::aruco::drawAxis(imageCopy, cameraMatrix, distCoeffs, rvecs[i], tvecs[i], 0.1);
       }
       cv::imshow("out", imageCopy);
       cv::waitKey(50);
       //if (key == 27)1
       // break;
   }
*/
//return 0;
}


