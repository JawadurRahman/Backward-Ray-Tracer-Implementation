#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "Image.h"
#include "Sphere.h"


using Colour = cv::Vec3b;

// Define two static RGB color constants
static Colour neonPurple(130, 0, 61);
static Colour lighterNeonPurple(255, 151, 255);

// Data types
typedef float Scalar;
typedef Eigen::Matrix<Scalar, 3, 1> Vec3;
typedef Eigen::Matrix<Scalar, 2, 1> Vec2;

// Triangle Class
class Triangle
{
public:
	// Constructor that initializes a triangle with three vertices
	Triangle(Vec3 v1, Vec3 v2, Vec3 v3) :vertex1_(v1), vertex2_(v2), vertex3_(v3) {	
		// Compute and store the normal vector of the triangle
		normal_ = ComputeNormal(v1, v2, v3);
	};

	// Calculate the barycentric coordinates of a point
	Vec3 BaryCentric(const Vec3&);		

	// Compute the normal vector of the triangle given its three vertices
	Vec3 ComputeNormal(const Vec3&, const Vec3&, const Vec3&);

	// Check if a ray intersects the triangle
	bool intersectRay(const Vec3& rayOrigin, const Vec3& rayDirection, float epsilon, Sphere sphere, Vec3& Intersection, Colour& floorColour);

public:
	// Vertice of the triangle
	Vec3 vertex1_, vertex2_, vertex3_;
	// Normal vector of the triangle
	Vec3 normal_;
	
};

