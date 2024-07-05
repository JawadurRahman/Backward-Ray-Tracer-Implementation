#include "Triangle.h"
#include <iostream>

Vec3 Triangle::BaryCentric(const Vec3& point)
{
	// Compute the normals from the point to the triangle's vertices
    Vec3 na = ComputeNormal(vertex2_, vertex3_, point);
    Vec3 nb = ComputeNormal(vertex3_, vertex1_, point);
    Vec3 nc = ComputeNormal(vertex1_, vertex2_, point);

	// Calculate and return the barycentric coordinates
    Scalar norm = normal_.squaredNorm();
    Scalar alpha = normal_.dot(na) / norm;
    Scalar beta = normal_.dot(nb) / norm;
    Scalar gamma = normal_.dot(nc) / norm;
    return Vec3(alpha, beta, gamma);
}


Vec3 Triangle::ComputeNormal(const Vec3& p1, const Vec3& p2, const Vec3& p3)
{
	// Compute and return the cross product of the two edges as the normal vector
    Vec3 v12 = p2 - p1;
    Vec3 v13 = p3 - p1;
    return v12.cross(v13);
}

bool Triangle::intersectRay(const Vec3& rayOrigin, const Vec3& rayDirection, float epsilon, Sphere sphere, Vec3& Intersection, Colour& floorColour)
{
	// Compute the normalized normal of the triangle
	Vec3 n = normal_;
	n.normalize();

	// Compute the intersection point of the ray-triangle intersection
	float t = (vertex1_ - rayOrigin).dot(n) / rayDirection.dot(n);
	Intersection = rayDirection * t + rayOrigin;
	
    // Compute the Barycentric coordinates of the Intersection Point
	Vec3 baryCoords = BaryCentric(Intersection);

	// Checking if the intersection point lies in the triangle using barycentric coordinates.
	// If it does then compute the colour of that point and return true.
	if (fabs(baryCoords.sum() - 1) < epsilon
		&& (baryCoords[0] <= 1 && baryCoords[0] >= 0)
		&& (baryCoords[1] <= 1 && baryCoords[1] >= 0)
		&& (baryCoords[2] <= 1 && baryCoords[2] >= 0)) {

		// Set triangle color to a grid pattern
		floorColour = neonPurple;		
		if (Intersection.x() - floor(Intersection.x()) <= 0.1)
		{
			floorColour = lighterNeonPurple;
		}
		if (Intersection.z() - floor(Intersection.z()) <= 0.1)
		{
			floorColour = lighterNeonPurple;
		}

		// Compute difffuse term
		Vec3 LightVector = lightSource - Intersection;
		LightVector.normalize();		
		float diffuseTerm = LightVector.dot(n);		
		if (diffuseTerm < 0) diffuseTerm = 0;

		// Calculate the specular term
		Vec3 ViewVector = -rayDirection;
		Vec3 HalfVec = (LightVector + ViewVector).normalized();
		float specularTerm = HalfVec.dot(n);
		if (specularTerm < 0) specularTerm = 0;
		
		// Specifying the ambient, diffuse and specular light colors 
		Vec3 ambLight(1.0f, 1.0f, 1.0f); // white
		Vec3 diffLight(1.0f, 1.0f, 1.0f); // white
		Vec3 specLight(1.0f, 1.0f, 1.0f); // white

		// The reflection coefficients of floor/triangle
		Vec3 k_a = Vec3(0.2f, 0.2f, 0.2f);
		Vec3 k_d = Vec3(0.3f, 0.8f, 1.0f);
		Vec3 k_s = Vec3(1.0f, 1.0f, 1.0f);

		// Check for shadow raw intersection with a sphere. Set specularTerm and diffuseTerm to 0 if so.
		float t = sphere.intersectRay(Intersection, LightVector);
		if (t > 0) {
			specularTerm = 0;
			diffuseTerm = 0;
		}

		// Calculate the final color using Phong shading model
		Colour colour(0, 0, 0);
		int shininess = 32;
		colour[0] = Util::Clamp(((k_a[0] * ambLight[0] + k_d[0] * diffLight[0] * diffuseTerm + k_s[0] * specLight[0] * pow(specularTerm, shininess)) * floorColour[0]) );
		colour[1] = Util::Clamp(((k_a[1] * ambLight[1] + k_d[1] * diffLight[1] * diffuseTerm + k_s[1] * specLight[1] * pow(specularTerm, shininess)) * floorColour[1]) );
		colour[2] = Util::Clamp(((k_a[2] * ambLight[2] + k_d[2] * diffLight[2] * diffuseTerm + k_s[2] * specLight[2] * pow(specularTerm, shininess)) * floorColour[2]) );
		floorColour = colour;
        return true;
	}
    
    return false;
}




