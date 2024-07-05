#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>
typedef Eigen::Matrix<float, 3, 1> Vec3;

// Define the position of a light source 
static Vec3 lightSource(0.0f, 5.0f, 0.0f);
static Image world(cv::imread("./R3.png"));
// Utility class
class Util
{
public:
	// Static function to clamp an integer value to the range [0, 255]
	static unsigned char Clamp(int pixelCol)
	{
		if (pixelCol < 0) return 0;
		if (pixelCol >= 255) return 255;
		return pixelCol;
	}

};