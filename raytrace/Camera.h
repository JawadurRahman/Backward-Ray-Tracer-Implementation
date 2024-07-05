#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>

typedef Eigen::Matrix<float, 3, 1> Vec3;

// Camera class
class Camera
{
public:
	// Constructor that initializes a camera with its position
	Camera(Vec3 cameraPos) : mPosition_(cameraPos) {};

	// Method to generate a ray from the camera position to a given point
	inline Vec3 generateRay(Vec3 const& pt) {
		Vec3 direction = pt - mPosition_;
		return direction;
	}

	// Method to get the camera's current position
	inline Vec3 getPosition() const {
		return mPosition_;
	}

private:
	// Camera's position
	Vec3 mPosition_;
};